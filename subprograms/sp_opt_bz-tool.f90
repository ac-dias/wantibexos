
!funciona apenas para semicondutores

!gfortran -mcmodel=large sp_opt_bz-tool.f90 -o sp_opt_bz-tool.x -llapack95 -lopenblas -fopenmp

subroutine spoptpolbz(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	double precision,parameter:: pi=acos(-1.)

	integer :: dimrpa !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector

	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:,:) :: stto

	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores


	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	double complex :: hxsp,hysp,hzsp
 

	double precision :: auxx,auxy,auxz,auxsp,auxsm


	double complex,parameter :: imag=cmplx(0.0,1.0)

	integer :: nk

	double precision,allocatable,dimension(:,:) :: output

	double precision,parameter :: norm=(1.0/sqrt(2.))

	integer :: ngkpt

	integer,allocatable,dimension(:) :: nocpk

	!modificacoes versao 2.1

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	!character(len=70) :: diein    !ambiente dieletrico
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	double precision,dimension(3) :: ediel

	!fim modificacoes versao 2.1
	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"


	!OUTPUT
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bz_act_x.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_x output file"
	OPEN(UNIT=302, FILE= trim(outputfolder)//"bz_act_y.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_y output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"bz_act_z.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_z output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"bz_act_sp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_sp output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"bz_act_sm.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_sm output file"
	OPEN(UNIT=306, FILE=trim(outputfolder)//"dichroism.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening dichroism output file"




	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimrpa = nc*nv

	allocate(kpt(ngkpt,3))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)


	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,w90basis),vector(ngkpt,w90basis,w90basis))
	allocate(nocpk(ngkpt))


	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


			do j=1,w90basis
	
				eigv(i,j)= eaux(j)

			end do
			


			do l=1,w90basis


				do h=1,w90basis

					vector(i,l,h)=vaux(l,h)


				end do
			

			end do

	
	end do
	! $omp end parallel do

	deallocate(eaux,vaux)



	!definindo os numeros quanticos dos estados




	allocate (stto(ngkpt,nc*nv,3))

do i=1,ngkpt

	counter=1
	
	ncaux=(nocpk(i)+1)+nc-1

	nvaux=nocpk(i)-nv+1

	do c=(nocpk(i)+1),ncaux


		do v=nvaux,nocpk(1)


			stto(i,counter,1) = counter !numero do estado

			stto(i,counter,2) = v   !numero da banda da valencia

			stto(i,counter,3) = c   !numero da banda da condução

			counter=counter+1
 			

		end do



	end do

end do

	counter=counter-1 !numero total de estados para equação bse

	!allocate(hoptxf(nk,dimrpa),hoptyf(nk,dimrpa),hoptspf(nk,dimrpa),hoptsmf(nk,dimrpa))

	!allocate(output(nk,7))

	!$omp do ordered
	do j=1,ngkpt

	auxx = 0.0
	auxy = 0.0
	auxz = 0.0
	auxsp = 0.0
	auxsm = 0.0

	! $omp critical
	do i=1,dimrpa

		call optspbz(vector(j,stto(j,i,2),:),vector(j,stto(j,i,3),:),&
		     kpt(j,1),kpt(j,2),kpt(j,3),&
		     ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)


		auxx = auxx+real(hxsp*conjg(hxsp))/(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)))	
		auxy = auxy+real(hysp*conjg(hysp))/(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)))
		auxz = auxz+real(hzsp*conjg(hzsp))/(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)))
		auxsp = auxsp&
		+real(norm*(hxsp+cmplx(0.,1.)*hysp)*conjg(norm*(hxsp+cmplx(0.,1.)*hysp)))/(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)))
	auxsm = auxsm&
	+real(norm*(hxsp-cmplx(0.,1.)*hysp)*conjg(norm*(hxsp-cmplx(0.,1.)*hysp)))/(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)))		



	end do
	! $omp end critical
		!output(j,1) = kpt(j,1)
		!output(j,2) = kpt(j,2)
		!output(j,3) = auxx
		!output(j,4) = auxy
		!output(j,5) = auxsp
		!output(j,6) = auxsm
		!output(j,7) = (auxsp-auxsm)/(auxsp+auxsm)
		!$omp ordered
		write(301,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxx
		write(302,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxy
		write(303,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxz
		write(304,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxsp
		write(305,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxsm
		write(306,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),(auxsp-auxsm)/(auxsp+auxsm)
		!$omp end ordered
	end do
	!$omp end do


	!do j=1,nk
	!	write(301,*) output(j,1),output(j,2),output(j,3)
	!	write(302,*) output(j,1),output(j,2),output(j,4)
	!	write(303,*) output(j,1),output(j,2),output(j,5)
	!	write(304,*) output(j,1),output(j,2),output(j,6)
	!	write(305,*) output(j,1),output(j,2),output(j,7)

	!end do



	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	!deallocate(hoptxf,hoptyf,hoptspf,hoptsmf)
	deallocate(stto)
	deallocate(kpt)
	!deallocate(output)





	close(200)

	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
	close(306)

	


end subroutine spoptpolbz
