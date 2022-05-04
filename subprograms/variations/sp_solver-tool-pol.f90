

!gfortran -mcmodel=large rpa_optics-tool.f90 -o rpa_opt-tool.x -llapack95 -lopenblas -fopenmp

subroutine spopticspol(nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype,tmcoef)


	use omp_lib
	use hamiltonian_input_variables

	implicit none

	double precision,parameter:: pi=acos(-1.)

	integer :: dimsp !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector

	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:) :: stt

	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	!double precision,allocatable,dimension(:) :: ensp

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	double complex:: matrizelbse

	double complex :: hxsp,hysp,hzsp

	double complex,allocatable,dimension(:) :: hoptxf,hoptyf,hoptzf,hoptspf,hoptsmf 


	double complex,parameter :: imag=cmplx(0.0,1.0)


	integer :: ngkpt
	!real,allocatable,dimension(:) :: auxx,auxy,auxz,auxsp,auxsm
	real,allocatable,dimension(:,:) :: vecres

	integer,allocatable,dimension(:) :: nocpk

	!definicoes diagonalizacao LA_HEEVR
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
	COMPLEX*16,allocatable,dimension(:,:) :: hbse

        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
        double complex,allocatable,dimension (:) :: WORK

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
	double precision :: ec,ev
	logical :: tmcoef

	!fim modificacoes versao 2.1

	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"


	!OUTPUT
	OPEN(UNIT=301, FILE=trim(outputfolder)//"sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida sp pol optics"
    	
    	if (tmcoef) then
   	OPEN(UNIT=302, FILE=trim(outputfolder)//"tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida sp tmcoef pol" 
    	else
    	 continue
    	end if    	



	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)
	ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimsp = ngkpt*nc*nv


	allocate(kpt(ngkpt,3))
	allocate(nocpk(ngkpt))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,w90basis),vector(ngkpt,w90basis,w90basis))


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



	allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)


	allocate(hoptxf(dimsp),hoptyf(dimsp),hoptspf(dimsp),hoptsmf(dimsp))
	allocate(hoptzf(dimsp))
	!allocate(ensp(dimsp))
	!allocate(auxx(dimsp),auxy(dimsp),auxz(dimsp))
	!allocate(auxsp(dimsp),auxsm(dimsp))

	allocate(vecres(dimsp,12))

	write(301,*) "#","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	
	
	if (tmcoef) then
	write(302,*) "number of kpoints:",ngkpt
	write(302,*) "number of conduction states",nc
	write(302,*) "number of valence states",nv
	write(302,*) "#"," ", "kx", " ", "ky"," ", "kz"," ", "nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"	
	else
	
	 continue
	end if
	!$omp do 

	do i=1,dimsp
	
		ec = eigv(stt(i,4),stt(i,3))
		ev = eigv(stt(i,4),stt(i,2))

		     
	       call optsp(eigv(stt(i,4),stt(i,2)),vector(stt(i,4),stt(i,2),:),&
		     eigv(stt(i,4),stt(i,3)),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)     


		hoptxf(i)= hxsp
		hoptyf(i)= hysp
		hoptzf(i)= hzsp
		hoptspf(i)= (hxsp+cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		hoptsmf(i)= (hxsp-cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		
		
		
		!ensp(i) = eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2))

		!auxx(i) = real(hoptxf(i)*conjg(hoptxf(i)))
		!auxy(i) = real(hoptyf(i)*conjg(hoptyf(i)))
		!auxz(i) = real(hoptzf(i)*conjg(hoptzf(i)))
		!auxsp(i) = real(hoptspf(i)*conjg(hoptspf(i)))
		!auxsm(i) = real(hoptsmf(i)*conjg(hoptsmf(i)))
		
		
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(stt(i,3))
		vecres(i,6) = real(stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hoptxf(i)*conjg(hoptxf(i)))		
		vecres(i,9) = real(hoptyf(i)*conjg(hoptyf(i)))		
		vecres(i,10) = real(hoptzf(i)*conjg(hoptzf(i)))
		vecres(i,11) = real(hoptspf(i)*conjg(hoptspf(i)))
		vecres(i,12) = real(hoptsmf(i)*conjg(hoptsmf(i)))
								
						
		! $omp ordered
		!write(301,"(6F15.4)") ensp(i),auxx(i),auxy(i),auxz(i),auxsp(i),auxsm(i)
		
		!if (tmcoef) then
	        !write(302,"(3F10.4,3I10.0,6F10.4)") kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),nocpk(stt(i,4)),stt(i,3),&
		!				  stt(i,2),ensp(i),auxx(i),auxy(i),auxz(i),auxsp(i),auxsm(i)
		! $omp end ordered
		!else
		! continue
		!end if 

	end do

	!$omp end do

	call Bubblem(7,12,vecres, dimsp)
	
	do i=1,dimsp
	
		write(301,"(6F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12)
		
		if (tmcoef) then
	        write(302,"(3F10.4,3I10.0,6F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12)
		end if
	
	end do


	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(hoptxf,hoptyf,hoptzf,hoptspf,hoptsmf)
	!deallocate(auxx,auxy,auxz,auxsp,auxsm)
	deallocate(nocpk)
	deallocate(stt)
	!deallocate(ensp)
	deallocate(kpt)
	deallocate(vecres)




	close(200)



	close(300)
	close(301)
	
	if (tmcoef) then
	close(302)
	else
	 continue
	end if
	


end subroutine spopticspol
