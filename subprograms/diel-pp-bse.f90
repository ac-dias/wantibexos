!ifort diel-pp-bse.f90 -o diel-pp-bse.x -qopenmp -mkl

subroutine bseoptprop(numbse,outputfolder)
 
	use omp_lib

	implicit none

	integer :: i,j,erro,dimpp

	double precision,allocatable,dimension(:) :: rxx,rxy,rxz
	double precision,allocatable,dimension(:) :: ryy,ryz,rzz

	double precision,allocatable,dimension(:) :: ixx,ixy,ixz
	double precision,allocatable,dimension(:) :: iyy,iyz,izz

	double precision,allocatable,dimension(:) :: energy

	double precision :: flag

	character(len=2) :: aread



	double precision :: refr_xx,refr_xy,refr_xz,refr_yy,refr_yz,refr_zz !indice refracao
	double precision :: ext_xx,ext_xy,ext_xz,ext_yy,ext_yz,ext_zz !indice extincao
	double precision :: refl_xx,refl_xy,refl_xz,refl_yy,refl_yz,refl_zz !reflectibilidade
	double precision :: abs_xx,abs_xy,abs_xz,abs_yy,abs_yz,abs_zz !coeficiente de absorcao
	double precision :: els_xx,els_xy,els_xz,els_yy,els_yz,els_zz !energy loss function	

	character(len=70) :: outputfolder
	double precision ::  numbse

	!call input_read

	dimpp=int(numbse)

	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_diel_xx.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xx input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"bse_diel_xy.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xy input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"bse_diel_xz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xz input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"bse_diel_yy.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_yy input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"bse_diel_yz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_yz input file"
	OPEN(UNIT=105, FILE=trim(outputfolder)//"bse_diel_zz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_zz input file"

	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"bse_refractive_index.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_refractive_index output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"bse_extinction_coef.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_extinction_coef output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"bse_reflectibility.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_reflectibility output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"bse_absorption_coef.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_absorption_coef output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"bse_en_loss_func.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_en_loss_func output file"	


	read(100,*) aread
	read(101,*) aread
	read(102,*) aread
	read(103,*) aread
	read(104,*) aread
	read(105,*) aread

	allocate(energy(dimpp))

	allocate(rxx(dimpp),ryy(dimpp),rzz(dimpp))
	allocate(rxy(dimpp),rxz(dimpp),ryz(dimpp))

	allocate(ixx(dimpp),iyy(dimpp),izz(dimpp))
	allocate(ixy(dimpp),ixz(dimpp),iyz(dimpp))

	do i=1,dimpp

		read(100,*) energy(i),rxx(i),ixx(i)
		read(101,*) flag,rxy(i),ixy(i)
		read(102,*) flag,rxz(i),ixz(i)
		read(103,*) flag,ryy(i),iyy(i)
		read(104,*) flag,ryz(i),iyz(i)
		read(105,*) flag,rzz(i),izz(i)

	end do

	write(200,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(300,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(400,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(500,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(600,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"	


	do i=1,dimpp

		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rxy(i),ixy(i),refr_xy)
		call refracao(rxz(i),ixz(i),refr_xz)
		call refracao(ryz(i),iyz(i),refr_yz)

		write(200,"(7F15.4)") energy(i),refr_xx,refr_yy,refr_zz,refr_xy,refr_xz,refr_yz

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rxy(i),ixy(i),ext_xy)
		call extincao(rxz(i),ixz(i),ext_xz)
		call extincao(ryz(i),iyz(i),ext_yz)

		write(300,"(7F15.4)") energy(i),ext_xx,ext_yy,ext_zz,ext_xy,ext_xz,ext_yz

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rxy(i),ixy(i),refl_xy)
		call reflectibilidade(rxz(i),ixz(i),refl_xz)
		call reflectibilidade(ryz(i),iyz(i),refl_yz)

		write(400,"(7F15.4)") energy(i),refl_xx,refl_yy,refl_zz,refl_xy,refl_xz,refl_yz

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rxy(i),ixy(i),energy(i),abs_xy)
		call abscoef(rxz(i),ixz(i),energy(i),abs_xz)
		call abscoef(ryz(i),iyz(i),energy(i),abs_yz)

		write(500,"(7F15.4)") energy(i),abs_xx,abs_yy,abs_zz,abs_xy,abs_xz,abs_yz
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rxy(i),ixy(i),els_xy)
		call enloss(rxz(i),ixz(i),els_xz)
		call enloss(ryz(i),iyz(i),els_yz)

		write(600,"(7F15.4)") energy(i),els_xx,els_yy,els_zz,els_xy,els_xz,els_yz		


	end do


	deallocate(energy)

	deallocate(rxx,ryy,rzz)
	deallocate(rxy,rxz,ryz)

	deallocate(ixx,iyy,izz)
	deallocate(ixy,ixz,iyz)
	
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(105)

	close(200)
	close(300)
	close(400)
	close(500)
	close(600)



end subroutine bseoptprop
