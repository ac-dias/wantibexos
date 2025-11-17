!ifort diel-pp-bse.f90 -o diel-pp-bse.x -qopenmp -mkl

subroutine bseoptprop(numbse,outputfolder,tmax,ni,ns)
 
	use omp_lib

	implicit none

	integer :: i,j,erro,dimpp
	
	real :: tmax,ni,ns

	real,allocatable,dimension(:) :: rxx,rxy,rxz
	real,allocatable,dimension(:) :: ryy,ryz,rzz

	real,allocatable,dimension(:) :: ixx,ixy,ixz
	real,allocatable,dimension(:) :: iyy,iyz,izz

	real,allocatable,dimension(:) :: energy

	real :: flag

	character(len=2) :: aread



	real :: refr_xx,refr_xy,refr_xz,refr_yy,refr_yz,refr_zz !indice refracao
	real :: ext_xx,ext_xy,ext_xz,ext_yy,ext_yz,ext_zz !indice extincao
	real :: refl_xx,refl_xy,refl_xz,refl_yy,refl_yz,refl_zz !reflectibilidade
	real :: abs_xx,abs_xy,abs_xz,abs_yy,abs_yz,abs_zz !coeficiente de absorcao
	real :: els_xx,els_xy,els_xz,els_yy,els_yz,els_zz !energy loss function
	
	complex :: optc_xx,optc_xy,optc_xz,optc_yy,optc_yz,optc_zz		

	real :: tr_xx,tr_yy,tr_zz,tr_xy,tr_xz,tr_yz
	real :: ref_xx,ref_yy,ref_zz,ref_xy,ref_xz,ref_yz
	real :: abt_xx,abt_yy,abt_zz,abt_xy,abt_xz,abt_yz

	character(len=70) :: outputfolder
	real ::  numbse

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
	OPEN(UNIT=700, FILE=trim(outputfolder)//"bse_opt_cond_real.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_opt_cond_real output file"
	OPEN(UNIT=800, FILE=trim(outputfolder)//"bse_opt_cond_imag.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_opt_cond_imag output file"
	OPEN(UNIT=900, FILE=trim(outputfolder)//"bse_transmittance.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_transmittance.dat output file"
	OPEN(UNIT=901, FILE=trim(outputfolder)//"bse_reflectance.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_reflectance output file"
	OPEN(UNIT=902, FILE=trim(outputfolder)//"bse_absorptance.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_absorptance output file"	

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
	write(700,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(800,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"	
	write(900,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(901,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(902,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"

	do i=1,dimpp

		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rxy(i),ixy(i),refr_xy)
		call refracao(rxz(i),ixz(i),refr_xz)
		call refracao(ryz(i),iyz(i),refr_yz)

		write(200,"(7F15.6)") energy(i),refr_xx,refr_yy,refr_zz,refr_xy,refr_xz,refr_yz

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rxy(i),ixy(i),ext_xy)
		call extincao(rxz(i),ixz(i),ext_xz)
		call extincao(ryz(i),iyz(i),ext_yz)

		write(300,"(7F15.6)") energy(i),ext_xx,ext_yy,ext_zz,ext_xy,ext_xz,ext_yz

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rxy(i),ixy(i),refl_xy)
		call reflectibilidade(rxz(i),ixz(i),refl_xz)
		call reflectibilidade(ryz(i),iyz(i),refl_yz)

		write(400,"(7F15.6)") energy(i),refl_xx,refl_yy,refl_zz,refl_xy,refl_xz,refl_yz

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rxy(i),ixy(i),energy(i),abs_xy)
		call abscoef(rxz(i),ixz(i),energy(i),abs_xz)
		call abscoef(ryz(i),iyz(i),energy(i),abs_yz)

		write(500,"(7F15.6)") energy(i),abs_xx,abs_yy,abs_zz,abs_xy,abs_xz,abs_yz
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rxy(i),ixy(i),els_xy)
		call enloss(rxz(i),ixz(i),els_xz)
		call enloss(ryz(i),iyz(i),els_yz)

		write(600,"(7F15.6)") energy(i),els_xx,els_yy,els_zz,els_xy,els_xz,els_yz		

		call optcondcalc(rxx(i),ixx(i),energy(i),optc_xx)
		call optcondcalc(ryy(i),iyy(i),energy(i),optc_yy)
		call optcondcalc(rzz(i),izz(i),energy(i),optc_zz)
		call optcondcalc(rxy(i),ixy(i),energy(i),optc_xy)
		call optcondcalc(rxz(i),ixz(i),energy(i),optc_xz)
		call optcondcalc(ryz(i),iyz(i),energy(i),optc_yz)										

		write(700,"(1F15.6,6E15.6)") energy(i),real(optc_xx),real(optc_yy),real(optc_zz),real(optc_xy),real(optc_xz),real(optc_yz)
		write(800,"(1F15.6,6E15.6)") energy(i),aimag(optc_xx),aimag(optc_yy),aimag(optc_zz),aimag(optc_xy),aimag(optc_xz),aimag(optc_yz)

		call transmittance(ni,ns,tmax,optc_xx,tr_xx)
		call transmittance(ni,ns,tmax,optc_yy,tr_yy)
		call transmittance(ni,ns,tmax,optc_zz,tr_zz)
		call transmittance(ni,ns,tmax,optc_xy,tr_xy)
		call transmittance(ni,ns,tmax,optc_xz,tr_xz)
		call transmittance(ni,ns,tmax,optc_yz,tr_yz)	
		
		write(900,"(1F15.6,6E15.6)")  energy(i),tr_xx,tr_yy,tr_zz,tr_xy,tr_xz,tr_yz							


		call reflectance(ni,ns,tmax,optc_xx,ref_xx)
		call reflectance(ni,ns,tmax,optc_yy,ref_yy)
		call reflectance(ni,ns,tmax,optc_zz,ref_zz)
		call reflectance(ni,ns,tmax,optc_xy,ref_xy)
		call reflectance(ni,ns,tmax,optc_xz,ref_xz)
		call reflectance(ni,ns,tmax,optc_yz,ref_yz)	
		
		write(901,"(1F15.6,6E15.6)")  energy(i),ref_xx,ref_yy,ref_zz,ref_xy,ref_xz,ref_yz		

		abt_xx = 1.0 - ref_xx - tr_xx
		abt_yy = 1.0 - ref_yy - tr_yy
		abt_zz = 1.0 - ref_zz - tr_zz
		abt_xy = 1.0 - ref_xy - tr_xy
		abt_xz = 1.0 - ref_xz - tr_xz
		abt_yz = 1.0 - ref_yz - tr_yz										

		write(902,"(1F15.6,6E15.6)")  energy(i),abt_xx,abt_yy,abt_zz,abt_xy,abt_xz,abt_yz

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
	close(700)
	close(800)
	close(900)
	close(901)
	close(902)	

end subroutine bseoptprop
