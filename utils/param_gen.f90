!script to generate the tight binding parameters file
!from wannier90 output in the format that our code reads it
!wannier90 output files: wannier90_hr.dat, wannier90.win

!gfortran param_gen.f90 -o param_gen.x

program main

	implicit none

	integer :: i,tipo,erro
	double precision :: efermi,ediel,lc
	double precision :: scs
	
	OPEN(UNIT=100, FILE="tb_hr.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening output file"

	write(*,*) "Write type of DFT calculation"
	write(*,*) "1 - Non-Polarized"
	write(*,*) "2 - Spin-Polarized"
	write(*,*) "3 - Spin Orbit Coupling"
	read(*,*)  tipo

	select case (tipo)

		case(1)
			write(100,*) "NP"
		
		case(2)
			write(100,*) "SP"

		case(3)
			write(100,*) "SOC"

		case default
			write(*,*) "wrong value of DFT calculation type"
			STOP

	end select

	write(*,*) "Write Fermi level Energy (eV)"
	read(*,*) efermi

	!write(*,*) "Write Effective dielectric constant"
	!read(*,*) ediel

	!write(*,*) "Write LC (angstrom), for 3D systems put 0.0"
	!read(*,*) lc

	write(*,*) "Write Scissors Operator"
	read(*,*) scs

	!write(100,*) ediel, " ","!effective dielectric"
	write(100,*) scs, " ", "!scissors operator"
	write(100,*) efermi,"  ", "!efermi"
	!write(100,*) lc, " ","!LC"



	if (tipo .eq. 2) then

	call SYSTEM(" awk '/''begin unit_cell_cart''/{f=1;next} /''end unit_cell_cart''/{f=0} f' wannier90.up.win >> tb_hr.dat ")
	call system(" cat wannier90.up_hr.dat >> tb_hr.dat ")
	call system(" cat wannier90.dn_hr.dat >> tb_hr.dat ")

	else

	call SYSTEM(" awk '/''begin unit_cell_cart''/{f=1;next} /''end unit_cell_cart''/{f=0} f' wannier90.win >> tb_hr.dat ")
	call system(" cat wannier90_hr.dat >> tb_hr.dat ")

	end if

	close(100)



end program main

!awk 'NR==3{ print; }' tb_hr.dat | awk '{print $2;}'

