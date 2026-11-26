#!/usr/bin/env python3
"""Regression of the symmetry analysis (bin/wtb_symmetry.x, the same code prints it in log_bse_optics.dat) and of the
tolerant reading of the input file.

Usage: python3 tests/test_symmetry_app.py /absolute/path/to/wtb_symmetry.x

Generated one-orbital-per-site honeycomb lattices (DFT=S text files with a basis_set-NP file) are analysed:
  * graphene (equal sites)                      D6h, order 24, primitive cell
  * hBN (two species, two on-site energies)      D3h, order 12
  * a 2 x 1 supercell of graphene               two pure translations, D2h
  * graphene with one stronger bond (H only)     the structure is D6h, the Hamiltonian is not: a warning, smaller usable group
  * the k mesh: 12 x 12 keeps every operation, 12 x 10 does not keep the threefold axis, a shift of 1/36 breaks k -> -k
The input file is also written with tabs, without blanks around the = sign, with lower case keys and comments: the
result must not change.  No third-party Python packages are needed.
"""
import argparse
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

A = 2.46
BOND = A/math.sqrt(3)


def honeycomb(path, n1=1, n2=1, species=('C', 'C'), onsite=(0., 0.), bond_t=None, t=-2.7):
    """tb.dat and basis_set-NP of an n1 x n2 supercell of the honeycomb lattice, one orbital per atom"""
    a1 = (A, 0.)
    a2 = (-A/2, A*math.sqrt(3)/2)
    sa1 = (n1*a1[0], n1*a1[1])
    sa2 = (n2*a2[0], n2*a2[1])
    sites = []
    for i in range(n1):
        for j in range(n2):
            o = (i*a1[0]+j*a2[0], i*a1[1]+j*a2[1])
            sites.append((o[0], o[1], 0))
            sites.append((o[0]+(2*a1[0]+a2[0])/3, o[1]+(2*a1[1]+a2[1])/3, 1))
    norb = len(sites)
    rs = [(p, q) for p in range(-2, 3) for q in range(-2, 3)]
    entries = {}
    for s, (xs, ys, ks) in enumerate(sites):
        for u, (xu, yu, ku) in enumerate(sites):
            for (p, q) in rs:
                rx = p*sa1[0]+q*sa2[0]
                ry = p*sa1[1]+q*sa2[1]
                dx, dy = xu+rx-xs, yu+ry-ys
                if abs(math.hypot(dx, dy)-BOND) < 1e-3 and ks != ku:
                    tt = t
                    if bond_t is not None:
                        # the bond along the first A-B direction (30 degrees, 210 for B-A) is different
                        ang = math.degrees(math.atan2(dy, dx))
                        if abs(ang-30.) < 1. or abs(ang+150.) < 1.:
                            tt = bond_t
                    entries[(s, u, p, q)] = tt
    vecs = sorted({(p, q) for (_, _, p, q) in entries} | {(0, 0)})
    lines = ['NP', '0.0', '0.0', f'{sa1[0]} {sa1[1]} 0', f'{sa2[0]} {sa2[1]} 0', '0 0 20', str(norb), str(len(vecs)),
             '# x y z a b h ih s # #']
    for (p, q) in vecs:
        rx = p*sa1[0]+q*sa2[0]
        ry = p*sa1[1]+q*sa2[1]
        for u in range(norb):
            for s in range(norb):
                h = entries.get((s, u, p, q), 0.0)
                ov = 1.0 if (s == u and (p, q) == (0, 0)) else (0.05 if (s, u, p, q) in entries else 0.0)
                if s == u and (p, q) == (0, 0):
                    h = onsite[sites[s][2]]
                lines.append(f'{rx} {ry} 0 {s+1} {u+1} {h} 0 {ov}')
    (path/'tb.dat').write_text('\n'.join(lines)+'\n')
    basis = ['bindex aspecie ax ay az l m spin']
    for s, (x, y, k) in enumerate(sites):
        basis.append(f'{s+1} {species[k]} {x} {y} 10.0 1 0 0')
    (path/'basis_set-NP').write_text('\n'.join(basis)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    root = Path(tempfile.mkdtemp(prefix='wtb-symmetry-app-'))
    env = dict(os.environ, OMP_NUM_THREADS='1')
    runs = 0

    def analyse(name, mesh=(12, 12), shift=(0., 0.), style='plain', **model):
        nonlocal runs
        path = root/name
        path.mkdir()
        honeycomb(path, **model)
        items = [('DFT', '"S"'), ('PARAMS_FILE', '"tb.dat"'), ('PARAMS_FORMAT', '"text"'), ('NGX', mesh[0]), ('NGY', mesh[1]),
                 ('NGZ', 1), ('SHIFT_1', shift[0]), ('SHIFT_2', shift[1]), ('SHIFT_3', 0.0), ('NOCPF', 1), ('NBANDSC', 1),
                 ('NBANDSV', 1)]
        if style == 'plain':
            text = ''.join(f'{k}= {v}\n' for k, v in items)
        else:
            lines = ['# a comment line', '']
            for i, (k, v) in enumerate(items):
                k = k.lower() if i % 2 else k
                lines.append([f'{k}={v}', f'{k} = {v}', f'\t{k}\t=\t{v}\t# comment', f'  {k}= {v} ! comment', f'{k} ={v}'][i % 5])
            text = '\r\n'.join(lines)+'\r\n'
        (path/'input.dat').write_text(text)
        proc = subprocess.run([str(executable)], cwd=path, env=env, capture_output=True, text=True, input=text, timeout=120)
        (path/'run.log').write_text(proc.stdout+proc.stderr)
        assert proc.returncode == 0, (name, path/'run.log')
        runs += 1
        out = proc.stdout
        group = re.search(r'Point group: (\S+)\s+\(Hermann-Mauguin (\S+)\),\s+order (\d+)', out)
        assert group, (name, 'no point group', path/'run.log')
        result = {'group': group.group(1), 'order': int(group.group(3)), 'text': out}
        m = re.search(r'Pure translations of the cell: (\d+)', out)
        result['ntrans'] = int(m.group(1)) if m else None
        m = re.search(r'map the mesh onto itself: (\d+) of (\d+)', out)
        result['mesh'] = (int(m.group(1)), int(m.group(2)))
        result['tr'] = 'the mesh is symmetric under k -> -k' in out
        result['warning'] = 'WARNING: the Hamiltonian does not have the full symmetry' in out
        m = re.search(r'Group usable with this mesh and this Hamiltonian: (\S+)\s+\(\S+\), order (\d+)', out)
        result['usable'] = (m.group(1), int(m.group(2))) if m else None
        m = re.search(r'Hamiltonian check: largest deviation of the site-pair norms of H\(R\)\s+(\S+)', out)
        result['hdev'] = float(m.group(1)) if m else None
        return result

    r = analyse('graphene')
    assert (r['group'], r['order'], r['ntrans']) == ('D6h', 24, 1), r['text']
    assert r['mesh'] == (24, 24) and r['tr'] and not r['warning'] and r['usable'] == ('D6h', 24), r['text']
    assert r['hdev'] is not None and r['hdev'] < 1e-4, r['text']
    # the same input in a loose format
    r2 = analyse('graphene-loose', style='loose')
    assert (r2['group'], r2['order'], r2['mesh']) == (r['group'], r['order'], r['mesh']), r2['text']

    r = analyse('hbn', species=('B', 'N'), onsite=(-1.5, 1.5))
    assert (r['group'], r['order']) == ('D3h', 12), r['text']
    assert r['mesh'] == (12, 12) and r['usable'] == ('D3h', 12), r['text']

    r = analyse('supercell', n1=2, n2=1, mesh=(6, 12))
    assert r['ntrans'] == 2 and r['group'] == 'D2h' and r['order'] == 8, r['text']

    r = analyse('one-bond', bond_t=-3.2)
    assert r['group'] == 'D6h' and r['warning'], r['text']
    assert r['usable'] is not None and r['usable'][1] < 24, r['text']

    # the k mesh
    r = analyse('mesh-12x10', mesh=(12, 10))
    assert r['mesh'][0] < r['mesh'][1] and r['mesh'][0] >= 4 and r['tr'], r['text']
    r = analyse('mesh-shifted', shift=(1/36, 1/36))
    assert not r['tr'], r['text']
    r = analyse('mesh-half-step', shift=(1/24, 1/24))
    assert r['tr'], r['text']

    print(f'symmetry_app: PASS ({runs} runs); artifacts: {root}')


if __name__ == '__main__':
    main()
