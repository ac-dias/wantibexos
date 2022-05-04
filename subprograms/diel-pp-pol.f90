
!ifort diel-pp.f90 -o diel-pp.x -qopenmp -mkl
!bseoptproppol
subroutine spoptproppol(numbse,outputfolder)

	use omp_lib

	implicit none

	integer :: i,j,erro,dimpp

	double precision,allocatable,dimension(:) :: rxx,ryy,rzz
	double precision,allocatable,dimension(:) :: rsp,rsm

	double precision,allocatable,dimension(:) :: ixx,iyy,izz
	double precision,allocatable,dimension(:) :: isp,ism

	double precision,allocatable,dimension(:) :: energy

	double precision :: flag

	character(len=2) :: aread

	double precision :: refr_xx,refr_sp,refr_sm,refr_yy,refr_zz !indice refracao
	double precision :: ext_xx,ext_sp,ext_sm,ext_yy,ext_zz !indice extincao
	double precision :: refl_xx,refl_sp,refl_sm,refl_yy,refl_zz !reflectibilidade
	double precision :: abs_xx,abs_sp,abs_sm,abs_yy,abs_zz !coeficiente de absorcao
	double precision :: els_xx,els_sp,els_sm,els_yy,els_zz !energy loss function	

	character(len=70) :: outputfolder
	double precision ::  numbse

	!call input_read

	dimpp=int(numbse)

	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"sp_diel-pol_x.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_x input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"sp_diel-pol_y.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_y input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"sp_diel-pol_z.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_z input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"sp_diel-pol_sp.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_sp input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"sp_diel-pol_sm.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_sm input file"


	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"sp_refractive_index-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_refractive_index-pol output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"sp_extinction_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_extinction_coef-pol output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"sp_reflectibility-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_reflectibility-pol output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"sp_absorption_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_absorption_coef-pol output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"sp_en_loss_func-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_en_loss_func-pol output file"		

	read(100,*) aread
	read(101,*) aread
	read(102,*) aread
	read(103,*) aread
	read(104,*) aread


	allocate(energy(dimpp))

	allocate(rxx(dimpp),ryy(dimpp),rzz(dimpp))
	allocate(rsp(dimpp),rsm(dimpp))

	allocate(ixx(dimpp),iyy(dimpp),izz(dimpp))
	allocate(isp(dimpp),ism(dimpp))

	do i=1,dimpp

		read(100,*) energy(i),rxx(i),ixx(i)
		read(101,*) flag,ryy(i),iyy(i)
		read(102,*) flag,rzz(i),izz(i)
		read(103,*) flag,rsp(i),isp(i)
		read(104,*) flag,rsm(i),ism(i)


	end do

	write(200,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(300,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(400,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(500,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(600,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"	

	do i=1,dimpp

		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rsp(i),isp(i),refr_sp)
		call refracao(rsm(i),ism(i),refr_sm)


		write(200,"(6F15.4)") energy(i),refr_xx,refr_yy,refr_zz,refr_sp,refr_sm

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rsp(i),isp(i),ext_sp)
		call extincao(rsm(i),ism(i),ext_sm)


		write(300,"(6F15.4)") energy(i),ext_xx,ext_yy,ext_zz,ext_sp,ext_sm

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rsp(i),isp(i),refl_sp)
		call reflectibilidade(rsm(i),ism(i),refl_sm)


		write(400,"(6F15.4)") energy(i),refl_xx,refl_yy,refl_zz,refl_sp,refl_sm

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rsp(i),isp(i),energy(i),abs_sp)
		call abscoef(rsm(i),ism(i),energy(i),abs_sm)


		write(500,"(6F15.4)") energy(i),abs_xx,abs_yy,abs_zz,abs_sp,abs_sm
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rsp(i),isp(i),els_sp)
		call enloss(rsm(i),ism(i),els_sm)


		write(600,"(6F15.4)") energy(i),els_xx,els_yy,els_zz,els_sp,els_sm		




	end do


	deallocate(energy)

	deallocate(rxx,ryy,rzz)
	deallocate(rsp,rsm)

	deallocate(ixx,iyy,izz)
	deallocate(isp,ism)
	
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)


	close(200)
	close(300)
	close(400)
	close(500)
	close(600)



end subroutine spoptproppol

