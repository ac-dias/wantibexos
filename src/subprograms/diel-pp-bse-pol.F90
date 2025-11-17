
!ifort diel-pp.f90 -o diel-pp.x -qopenmp -mkl
!bseoptproppol
subroutine bseoptproppol(numbse,outputfolder,tmax,ni,ns)

	use omp_lib

	implicit none

	integer :: i,j,erro,dimpp
	
	real :: tmax,ni,ns	

	real,allocatable,dimension(:) :: rxx,ryy,rzz
	real,allocatable,dimension(:) :: rsp,rsm

	real,allocatable,dimension(:) :: ixx,iyy,izz
	real,allocatable,dimension(:) :: isp,ism

	real,allocatable,dimension(:) :: energy

	real :: flag

	character(len=2) :: aread

	real :: refr_xx,refr_sp,refr_sm,refr_yy,refr_zz !indice refracao
	real :: ext_xx,ext_sp,ext_sm,ext_yy,ext_zz !indice extincao
	real :: refl_xx,refl_sp,refl_sm,refl_yy,refl_zz !reflectibilidade
	real :: abs_xx,abs_sp,abs_sm,abs_yy,abs_zz !coeficiente de absorcao
	real :: els_xx,els_sp,els_sm,els_yy,els_zz !energy loss function
	
	complex :: optc_xx,optc_sp,optc_sm,optc_yy,optc_zz
	
	real :: tr_xx,tr_yy,tr_zz,tr_sp,tr_sm
	real :: ref_xx,ref_yy,ref_zz,ref_sp,ref_sm
	real :: abt_xx,abt_yy,abt_zz,abt_sp,abt_sm		

	character(len=70) :: outputfolder
	real ::  numbse

	!call input_read

	dimpp=int(numbse)

	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_diel-pol_x.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_x input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"bse_diel-pol_y.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_y input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"bse_diel-pol_z.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_z input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"bse_diel-pol_sp.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sp input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"bse_diel-pol_sm.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sm input file"


	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"bse_refractive_index-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_refractive_index-pol output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"bse_extinction_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_extinction_coef-pol output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"bse_reflectibility-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_reflectibility-pol output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"bse_absorption_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_absorption_coef-pol output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"bse_en_loss_func-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_en_loss_func-pol output file"	
	OPEN(UNIT=700, FILE=trim(outputfolder)//"bse_opt_cond_real-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_opt_cond_real-pol output file"
	OPEN(UNIT=800, FILE=trim(outputfolder)//"bse_opt_cond_imag-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_opt_cond_imag-pol output file"
	OPEN(UNIT=900, FILE=trim(outputfolder)//"bse_transmittance-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_transmittance-pol.dat output file"
	OPEN(UNIT=901, FILE=trim(outputfolder)//"bse_reflectance-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_reflectance-pol output file"
	OPEN(UNIT=902, FILE=trim(outputfolder)//"bse_absorptance-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_absorptance-pol output file"

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
	write(700,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(800,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"	
	write(900,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(901,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(902,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"

	do i=1,dimpp

		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rsp(i),isp(i),refr_sp)
		call refracao(rsm(i),ism(i),refr_sm)


		write(200,"(6F15.6)") energy(i),refr_xx,refr_yy,refr_zz,refr_sp,refr_sm

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rsp(i),isp(i),ext_sp)
		call extincao(rsm(i),ism(i),ext_sm)


		write(300,"(6F15.6)") energy(i),ext_xx,ext_yy,ext_zz,ext_sp,ext_sm

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rsp(i),isp(i),refl_sp)
		call reflectibilidade(rsm(i),ism(i),refl_sm)


		write(400,"(6F15.6)") energy(i),refl_xx,refl_yy,refl_zz,refl_sp,refl_sm

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rsp(i),isp(i),energy(i),abs_sp)
		call abscoef(rsm(i),ism(i),energy(i),abs_sm)


		write(500,"(6F15.6)") energy(i),abs_xx,abs_yy,abs_zz,abs_sp,abs_sm
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rsp(i),isp(i),els_sp)
		call enloss(rsm(i),ism(i),els_sm)


		write(600,"(6F15.6)") energy(i),els_xx,els_yy,els_zz,els_sp,els_sm		


		call optcondcalc(rxx(i),ixx(i),energy(i),optc_xx)
		call optcondcalc(ryy(i),iyy(i),energy(i),optc_yy)
		call optcondcalc(rzz(i),izz(i),energy(i),optc_zz)
		call optcondcalc(rsp(i),isp(i),energy(i),optc_sp)
		call optcondcalc(rsm(i),ism(i),energy(i),optc_sm)
										

		write(700,"(1F15.6,5E15.6)") energy(i),real(optc_xx),real(optc_yy),real(optc_zz),real(optc_sp),real(optc_sm)
		write(800,"(1F15.6,5E15.6)") energy(i),aimag(optc_xx),aimag(optc_yy),aimag(optc_zz),aimag(optc_sp),aimag(optc_sm)

		call transmittance(ni,ns,tmax,optc_xx,tr_xx)
		call transmittance(ni,ns,tmax,optc_yy,tr_yy)
		call transmittance(ni,ns,tmax,optc_zz,tr_zz)
		call transmittance(ni,ns,tmax,optc_sp,tr_sp)
		call transmittance(ni,ns,tmax,optc_sm,tr_sm)
	
		
		write(900,"(1F15.6,5E15.6)")  energy(i),tr_xx,tr_yy,tr_zz,tr_sp,tr_sm							


		call reflectance(ni,ns,tmax,optc_xx,ref_xx)
		call reflectance(ni,ns,tmax,optc_yy,ref_yy)
		call reflectance(ni,ns,tmax,optc_zz,ref_zz)
		call reflectance(ni,ns,tmax,optc_sp,ref_sp)
		call reflectance(ni,ns,tmax,optc_sm,ref_sm)
	
		
		write(901,"(1F15.6,5E15.6)")  energy(i),ref_xx,ref_yy,ref_zz,ref_sp,ref_sm		

		abt_xx = 1.0 - ref_xx - tr_xx
		abt_yy = 1.0 - ref_yy - tr_yy
		abt_zz = 1.0 - ref_zz - tr_zz
		abt_sp = 1.0 - ref_sp - tr_sp
		abt_sm = 1.0 - ref_sm - tr_sm
										

		write(902,"(1F15.6,5E15.6)")  energy(i),abt_xx,abt_yy,abt_zz,abt_sp,abt_sm

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
	close(700)
	close(800)
	close(900)
	close(901)
	close(902)

end subroutine bseoptproppol

