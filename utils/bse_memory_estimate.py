#!/usr/bin/env python3
"""Estimate the memory required by WanTiBEXOS BSE calculations.

The model follows the explicit allocatable and automatic arrays in the BSE
solvers.  It deliberately does not claim to include MPI, OpenMP, BLAS, ELPA,
or operating-system runtime allocations that are outside the source tree.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path


REAL_BYTES = 4
INTEGER_BYTES = 4
COMPLEX_BYTES = 8
MIB = 1024**2
GIB = 1024**3


@dataclass
class InputData:
    path: Path
    values: dict[str, str]

    def integer(self, key: str, default: int | None = None) -> int:
        value = self.values.get(key)
        if value is None:
            if default is None:
                raise ValueError(f"{key} is not set in {self.path}")
            return default
        try:
            return int(float(value))
        except ValueError as error:
            raise ValueError(f"{key}={value!r} is not an integer") from error

    def text(self, key: str, default: str | None = None) -> str:
        value = self.values.get(key, default)
        if value is None:
            raise ValueError(f"{key} is not set in {self.path}")
        return value.strip().strip("\"'")

    def logical(self, key: str, default: bool = False) -> bool:
        value = self.values.get(key)
        if value is None:
            return default
        normalized = value.strip().strip("\"'").upper().strip(".")
        if normalized in {"T", "TRUE", "1", "YES"}:
            return True
        if normalized in {"F", "FALSE", "0", "NO"}:
            return False
        raise ValueError(f"{key}={value!r} is not a logical value")


@dataclass
class HamiltonianInfo:
    path: Path
    basis: int
    nvec: int
    layout: str


def stripped_record(line: str) -> str:
    """Remove input-file comments while preserving Fortran values."""
    return line.split("!", 1)[0].split("#", 1)[0].strip()


def read_input(path: Path) -> InputData:
    values: dict[str, str] = {}
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$")
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(stripped_record(raw_line))
        if match:
            values[match.group(1).upper()] = match.group(2).strip()
    return InputData(path=path, values=values)


def resolve_from_input(input_path: Path, name: str) -> Path:
    candidate = Path(name)
    if candidate.is_absolute():
        return candidate
    beside_input = input_path.parent / candidate
    return beside_input if beside_input.exists() else candidate


def first_integer(record: str) -> int | None:
    try:
        return int(record.split()[0])
    except (IndexError, ValueError):
        return None


def read_hamiltonian(path: Path, dft: str) -> HamiltonianInfo:
    records = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(records) < 2:
        raise ValueError(f"Hamiltonian file {path} is too short")

    # Hückel conversion output begins directly with basis size and vector count.
    compact_basis = first_integer(records[0])
    compact_nvec = first_integer(records[1])
    if compact_basis is not None and compact_nvec is not None:
        return HamiltonianInfo(path, compact_basis, compact_nvec, "compact Hückel")

    if dft == "S":
        # hamiltonian_nort_input_read: systype, scs, efermi, 3 lattice rows,
        # w90basis, nvec, then a header record.
        if len(records) < 8:
            raise ValueError(f"Non-orthogonal Hamiltonian file {path} is too short")
        basis = first_integer(records[6])
        nvec = first_integer(records[7])
        if basis is None or nvec is None:
            raise ValueError(f"Could not read basis size and vector count from {path}")
        return HamiltonianInfo(path, basis, nvec, "non-orthogonal")

    # hamiltonian_input_read: systype, scs, efermi, 3 lattice rows, a label,
    # w90basis and nvec.  The SP form has two independent blocks and needs a
    # record-level skip through the spin-up matrix entries.
    if len(records) < 9:
        raise ValueError(f"Wannier Hamiltonian file {path} is too short")
    systype = records[0].split()[0].upper()
    if systype != "SP":
        basis = first_integer(records[7])
        nvec = first_integer(records[8])
        if basis is None or nvec is None:
            raise ValueError(f"Could not read basis size and vector count from {path}")
        return HamiltonianInfo(path, basis, nvec, "Wannier")

    basis_up = first_integer(records[7])
    nvec_up = first_integer(records[8])
    if basis_up is None or nvec_up is None:
        raise ValueError(f"Could not read spin-up dimensions from {path}")

    # The source consumes nvec_up degeneracy factors list-directed before one
    # record per matrix element.  Generated files normally use one factor row;
    # consume enough numeric tokens to permit wrapped factor lists as well.
    index = 9
    factors_seen = 0
    while index < len(records) and factors_seen < nvec_up:
        try:
            factors_seen += len([int(token) for token in records[index].split()])
        except ValueError as error:
            raise ValueError(f"Invalid spin-up degeneracy factors in {path}") from error
        index += 1
    if factors_seen < nvec_up:
        raise ValueError(f"Truncated spin-up degeneracy factors in {path}")
    index += nvec_up * basis_up * basis_up
    if index + 2 >= len(records):
        raise ValueError(f"Truncated spin-down block in {path}")
    basis_down = first_integer(records[index + 1])
    nvec_down = first_integer(records[index + 2])
    if basis_down is None or nvec_down is None:
        raise ValueError(f"Could not read spin-down dimensions from {path}")
    return HamiltonianInfo(path, basis_up + basis_down, nvec_up + nvec_down, "spin-polarized Wannier")


def qpoint_count(path: Path) -> int:
    records = [stripped_record(line) for line in path.read_text(encoding="utf-8").splitlines()]
    records = [record for record in records if record]
    if not records:
        raise ValueError(f"Q-point file {path} is empty")
    mode = records[0].split()[0].upper()
    if mode == "GRID":
        if len(records) < 2:
            raise ValueError(f"Q-grid file {path} has no dimensions")
        dimensions = [int(value) for value in records[1].split()[:3]]
        if len(dimensions) != 3 or any(value <= 0 for value in dimensions):
            raise ValueError(f"Invalid Q-grid dimensions in {path}")
        return math.prod(dimensions)
    if mode == "PATH":
        if len(records) < 3:
            raise ValueError(f"Q-path file {path} is incomplete")
        endpoints = int(records[1].split()[0])
        points_per_segment = int(records[2].split()[0])
    else:
        if len(records) < 2:
            raise ValueError(f"Legacy Q-path file {path} is incomplete")
        endpoints = int(records[0].split()[0])
        points_per_segment = int(records[1].split()[0])
    if endpoints <= 0 or endpoints % 2 or points_per_segment < 2:
        raise ValueError(f"Invalid Q-path dimensions in {path}")
    return (endpoints // 2) * points_per_segment


def human_bytes(value: int) -> str:
    if value >= GIB:
        return f"{value / GIB:,.3f} GiB"
    if value >= MIB:
        return f"{value / MIB:,.3f} MiB"
    if value >= 1024:
        return f"{value / 1024:,.3f} KiB"
    return f"{value:,} B"


def hamiltonian_resident_bytes(basis: int, nvec: int) -> int:
    # rvec, ffactor, hopmatrices, ihopmatrices and ovp.
    return nvec * (3 * REAL_BYTES + INTEGER_BYTES + 3 * basis * basis * REAL_BYTES)


def eigsys_thread_scratch_bytes(basis: int, nvec: int) -> int:
    """Automatic arrays plus eaux/vaux used by the OpenMP eigsys loops."""
    matrices = 5 * basis * basis * COMPLEX_BYTES
    lapack = (
        (1 + 5 * basis + 2 * basis * basis) * REAL_BYTES
        + (2 * basis + basis * basis) * COMPLEX_BYTES
        + (3 + 5 * basis) * INTEGER_BYTES
    )
    eigen_io = basis * REAL_BYTES + basis * basis * COMPLEX_BYTES
    kphase = nvec * COMPLEX_BYTES
    return matrices + lapack + eigen_io + kphase


def serial_workspace_bytes(algorithm: str, dimension: int) -> int:
    algorithm = algorithm.lower()
    if algorithm == "cheev":
        return max(0, 3 * dimension - 2) * REAL_BYTES + max(0, 2 * dimension - 1) * COMPLEX_BYTES
    if algorithm == "cheevr":
        return (
            dimension * dimension * COMPLEX_BYTES
            + 2 * dimension * COMPLEX_BYTES
            + 2 * dimension * INTEGER_BYTES
            + 10 * dimension * INTEGER_BYTES
            + 24 * dimension * REAL_BYTES
        )
    if algorithm == "cheevx":
        return (
            dimension * dimension * COMPLEX_BYTES
            + 2 * dimension * COMPLEX_BYTES
            + 5 * dimension * INTEGER_BYTES
            + 7 * dimension * REAL_BYTES
            + dimension * INTEGER_BYTES
        )
    if algorithm == "cheevd":
        return (
            (2 * dimension + dimension * dimension) * COMPLEX_BYTES
            + (3 + 5 * dimension) * INTEGER_BYTES
            + (1 + 5 * dimension + 2 * dimension * dimension) * REAL_BYTES
        )
    raise ValueError(f"Unsupported serial BSE_ALGO={algorithm!r}; use cheev, cheevr, cheevx, or cheevd")


def numroc(length: int, block: int, process: int, processes: int) -> int:
    full_blocks, remainder = divmod(length, block)
    owned = (full_blocks // processes) * block
    extra_blocks = full_blocks % processes
    if process < extra_blocks:
        owned += block
    elif process == extra_blocks:
        owned += remainder
    return owned


def process_grid(ranks: int) -> tuple[int, int]:
    rows = int(math.sqrt(ranks))
    while ranks % rows:
        rows -= 1
    return rows, ranks // rows


def local_matrix_bytes(dimension: int, ranks: int) -> list[int]:
    rows, columns = process_grid(ranks)
    sizes: list[int] = []
    for rank in range(ranks):
        row = rank // columns
        column = rank % columns
        local_rows = max(1, numroc(dimension, 64, row, rows))
        local_columns = max(1, numroc(dimension, 64, column, columns))
        sizes.append(local_rows * local_columns * COMPLEX_BYTES)
    return sizes


def base_optical_bytes(nk: int, dimension: int, basis: int, nvec: int, nc: int, nv: int, dft: str) -> int:
    value = hamiltonian_resident_bytes(basis, nvec)
    value += 3 * nk * REAL_BYTES  # kpt
    value += nk * (nc + nv) * REAL_BYTES  # eigv
    value += basis * (nc + nv) * nk * COMPLEX_BYTES  # vector
    value += nk * INTEGER_BYTES  # nocpk
    value += 4 * dimension * INTEGER_BYTES  # stt
    value += 5 * dimension * COMPLEX_BYTES  # hrx/hry/hrz/hrsp/hrsm
    value += 8 * dimension * REAL_BYTES  # optical tensor arrays
    if dft == "S":
        value += basis * basis * nk * COMPLEX_BYTES  # sk
    value += 4 * dimension * INTEGER_BYTES  # stt_bse
    value += 3 * nk * REAL_BYTES  # kpt_bse
    value += dimension * REAL_BYTES  # W
    return value


def base_qpath_bytes(nk: int, nq: int, dimension: int, basis: int, nvec: int, nc: int, nv: int, dft: str) -> int:
    value = hamiltonian_resident_bytes(basis, nvec)
    value += nq * 4 * REAL_BYTES  # qauxv
    value += 9 * nk * REAL_BYTES  # kpt, kpt_bse and qpt
    value += 2 * nk * (nc + nv) * REAL_BYTES  # energy and energyq
    value += 2 * basis * (nc + nv) * nk * COMPLEX_BYTES  # vector and vectorq
    value += 2 * nk * INTEGER_BYTES  # nocpk and nocpq
    value += dimension * nq * REAL_BYTES  # exk, replicated on every rank
    value += 8 * dimension * INTEGER_BYTES  # stt and stt_bse
    if dft == "S":
        value += 2 * basis * basis * nk * COMPLEX_BYTES  # sk and skq
    return value


def print_optical_estimate(
    nk: int,
    dimension: int,
    basis: int,
    nvec: int,
    nc: int,
    nv: int,
    dft: str,
    algorithm: str,
    ranks: int,
    threads: int,
    thread_reserve: int,
    safety_factor: float,
) -> None:
    dense = dimension * dimension * COMPLEX_BYTES
    base = base_optical_bytes(nk, dimension, basis, nvec, nc, nv, dft)
    thread_scratch = threads * (eigsys_thread_scratch_bytes(basis, nvec) + thread_reserve)
    print("\nOptical BSE (BSE=T)")
    print(f"  Dense Hamiltonian H({dimension},{dimension}): {human_bytes(dense)}")
    if ranks == 1:
        workspace = serial_workspace_bytes(algorithm, dimension)
        source_peak = base + dense + workspace
        planning_peak = math.ceil((source_peak + thread_scratch) * safety_factor)
        print(f"  BSE_ALGO={algorithm} source workspace: {human_bytes(workspace)}")
        print(f"  Explicit source peak, one rank: {human_bytes(source_peak)}")
        print(f"  Conservative plan with {threads} OpenMP threads: {human_bytes(planning_peak)}")
        return

    locals_ = local_matrix_bytes(dimension, ranks)
    source_peak = base + 2 * max(locals_)
    planning_peak = math.ceil((source_peak + thread_scratch) * safety_factor)
    rows, columns = process_grid(ranks)
    print(f"  MPI process grid: {rows} x {columns}; 64 x 64 ScaLAPACK blocks")
    print(f"  Local H storage: {human_bytes(min(locals_))} to {human_bytes(max(locals_))} per rank")
    print("  Peak includes H and one distributed eigenvector matrix; ELPA/ScaLAPACK internal workspace is excluded.")
    print(f"  Explicit source peak, most-loaded rank: {human_bytes(source_peak)}")
    print(f"  Conservative plan with {threads} OpenMP threads: {human_bytes(planning_peak)} per rank")
    print(f"  Conservative job total (all {ranks} ranks): {human_bytes(ranks * planning_peak)}")


def print_qpath_estimate(
    nk: int,
    nq: int,
    dimension: int,
    basis: int,
    nvec: int,
    nc: int,
    nv: int,
    dft: str,
    algorithm: str,
    ranks: int,
    threads: int,
    thread_reserve: int,
    safety_factor: float,
) -> None:
    dense = dimension * dimension * COMPLEX_BYTES
    base = base_qpath_bytes(nk, nq, dimension, basis, nvec, nc, nv, dft)
    workspace = serial_workspace_bytes(algorithm, dimension)
    thread_scratch = threads * (eigsys_thread_scratch_bytes(basis, nvec) + thread_reserve)
    source_peak = base + dense + workspace
    planning_peak = math.ceil((source_peak + thread_scratch) * safety_factor)
    active_ranks = min(ranks, nq)
    print("\nFinite-Q BSE path/grid (BSE_BND=T)")
    print(f"  Q points: {nq}; active ranks at once: {active_ranks} of {ranks}")
    print(f"  Dense Hamiltonian per active rank H({dimension},{dimension}): {human_bytes(dense)}")
    print(f"  BSE_ALGO={algorithm} source workspace: {human_bytes(workspace)}")
    print("  The full H is replicated for each active Q rank; MPI-over-Q does not distribute one H.")
    print(f"  Explicit source peak, active rank: {human_bytes(source_peak)}")
    print(f"  Conservative plan with {threads} OpenMP threads: {human_bytes(planning_peak)} per active rank")
    print(f"  Conservative job upper bound ({ranks} ranks): {human_bytes(ranks * planning_peak)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="WanTiBEXOS input file")
    parser.add_argument("--nodes", "--ranks", dest="ranks", type=int, default=1,
                        help="MPI ranks launched for the job (default: 1)")
    parser.add_argument("--threads", type=int, help="OpenMP threads per rank (default: NTHREADS in input)")
    parser.add_argument("--mode", choices=("auto", "optical", "qpath", "all"), default="auto",
                        help="calculation(s) to estimate (default: infer from BSE/BSE_BND)")
    parser.add_argument("--external-thread-mib", type=float, default=0.0,
                        help="extra unmodelled stack/BLAS reserve per OpenMP thread")
    parser.add_argument("--safety-factor", type=float, default=1.25,
                        help="multiply the source estimate for planning (default: 1.25)")
    args = parser.parse_args()

    if args.ranks <= 0 or args.safety_factor < 1.0 or args.external_thread_mib < 0.0:
        parser.error("ranks must be positive; safety factor must be at least 1; thread reserve cannot be negative")
    try:
        data = read_input(args.input)
        dft = data.text("DFT", "W").upper()[0]
        threads = args.threads if args.threads is not None else data.integer("NTHREADS", 1)
        if threads <= 0:
            raise ValueError("thread count must be positive")
        nx, ny, nz = (data.integer("NGX"), data.integer("NGY"), data.integer("NGZ"))
        nc, nv = data.integer("NBANDSC"), data.integer("NBANDSV")
        if min(nx, ny, nz, nc, nv) <= 0:
            raise ValueError("NGX, NGY, NGZ, NBANDSC and NBANDSV must be positive")
        params = resolve_from_input(args.input, data.text("PARAMS_FILE"))
        hamiltonian = read_hamiltonian(params, dft)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    nk = nx * ny * nz
    dimension = nk * nc * nv
    algorithm = data.text("BSE_ALGO", "cheev").lower()
    thread_reserve = round(args.external_thread_mib * MIB)
    optical_selected = data.logical("BSE")
    qpath_selected = data.logical("BSE_BND")
    if args.mode == "optical":
        optical_selected, qpath_selected = True, False
    elif args.mode == "qpath":
        optical_selected, qpath_selected = False, True
    elif args.mode == "all":
        optical_selected, qpath_selected = True, True
    if not optical_selected and not qpath_selected:
        print("error: the input enables neither BSE nor BSE_BND; use --mode to select an estimate", file=sys.stderr)
        return 2

    print(f"Input: {args.input}")
    print(f"Hamiltonian: {hamiltonian.path} ({hamiltonian.layout}; basis={hamiltonian.basis}, nvec={hamiltonian.nvec})")
    print(f"BSE basis: Nk={nk} ({nx}x{ny}x{nz}), Nv={nv}, Nc={nc}, dimension={dimension}")
    print(f"Execution model: {args.ranks} MPI rank(s), {threads} OpenMP thread(s) per rank")
    print("Assumptions: default REAL/INTEGER=4 B and COMPLEX=8 B; source arrays only.")
    print(f"Planning safety factor: {args.safety_factor:.2f}; external reserve: {args.external_thread_mib:.3f} MiB/thread")

    try:
        if optical_selected:
            print_optical_estimate(nk, dimension, hamiltonian.basis, hamiltonian.nvec, nc, nv, dft,
                                   algorithm, args.ranks, threads, thread_reserve, args.safety_factor)
        if qpath_selected:
            qfile = resolve_from_input(args.input, data.text("KPATH_BSE"))
            nq = qpoint_count(qfile)
            print_qpath_estimate(nk, nq, dimension, hamiltonian.basis, hamiltonian.nvec, nc, nv, dft,
                                 algorithm, args.ranks, threads, thread_reserve, args.safety_factor)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
