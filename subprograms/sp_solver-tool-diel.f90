!gfortran -mcmodel=large rpa_optics-tool-diel.f90 -o rpa_opt-tool-diel.x -llapack95 -lopenblas -fopenmp

subroutine spoptics(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
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

	!double precision,allocatable,dimension(:) :: esp

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	double complex:: matrizelbse
	
	double complex :: hxsp,hysp,hzsp,hrsp,hrsm

	!double complex :: exx,exy,exz,eyy,eyz,ezz


	double precision :: ec,ev

	double complex,parameter :: imag=cmplx(0.0,1.0)


	integer :: ngkpt
	!real,allocatable,dimension(:) :: auxx,auxy,auxz,auyy,auyz,auzz
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
	logical :: tmcoef

	!fim modificacoes versao 2.1
	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"


	!OUTPUT
	OPEN(UNIT=301, FILE=trim(outputfolder)//"sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=304, FILE=trim(outputfolder)//"tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=303, FILE=trim(outputfolder)//"sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics-pol output file"	
    	
  	



	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimsp = ngkpt*nc*nv


	allocate(kpt(ngkpt,3))
	allocate(nocpk(ngkpt))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))


	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


			do j=1,nc+nv
	
				eigv(i,j)= eaux(nocpk(i)-nv+j)

			end do
			


			do l=1,nc+nv


				do h=1,w90basis

					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)

				end do
			

			end do

	
	end do
	! $omp end parallel do
	deallocate(eaux,vaux)



	!definindo os numeros quanticos dos estados



	allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)

	!allocate(esp(dimsp))
	

	!allocate(auxx(dimsp),auxy(dimsp),auxz(dimsp))
	!allocate(auyy(dimsp),auyz(dimsp),auzz(dimsp))
	
	allocate(vecres(dimsp,15))

	write(301,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(303,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	
	if (tmcoef) then	
	write(302,*) "number of kpoints:",ngkpt
	write(302,*) "number of conduction states",nc
	write(302,*) "number of valence states",nv
	write(302,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(304,*) "number of kpoints:",ngkpt
	write(304,*) "number of conduction states",nc
	write(304,*) "number of valence states",nv
	write(304,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
		      
	else
	
	 continue
	end if
	! $omp do ordered
	!$omp do 

	do i=1,dimsp

		ec = eigv(stt(i,4),stt(i,3))
		ev = eigv(stt(i,4),stt(i,2))
		     
	       call optsp(eigv(stt(i,4),stt(i,2)),vector(stt(i,4),stt(i,2),:),&
		     eigv(stt(i,4),stt(i,3)),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)  
		
		hrsp = (hxsp+cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		hrsm = (hxsp-cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hxsp*conjg(hxsp))		
		vecres(i,9) = real(hysp*conjg(hysp))		
		vecres(i,10) = real(hzsp*conjg(hzsp))
		vecres(i,11) = real(hxsp*conjg(hysp))
		vecres(i,12) = real(hxsp*conjg(hzsp))
		vecres(i,13) = real(hysp*conjg(hzsp))		 
		vecres(i,14) = real(hrsp*conjg(hrsp))
		vecres(i,15) = real(hrsm*conjg(hrsm))		

		! $omp ordered
		!write(301,"(7F15.4)") esp(i),auxx(i),auyy(i),auzz(i),auxy(i),auxz(i),auyz(i)
		
		!if (tmcoef) then
		!write(302,"(3F10.4,3I10.0,7F10.4)") kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),nocpk(stt(i,4)),stt(i,3),&
		!				  stt(i,2),esp(i),auxx(i),auyy(i),auzz(i),auxy(i),auxz(i),auyz(i)
						  
		!else
		 !continue
		!end if 						  
		! $omp end ordered

		!write(2077,*) i,"/",dimrpa

	end do

	!$omp end do

	call Bubblem(7,15,vecres, dimsp)
	
	do i=1,dimsp
	
		write(301,"(7F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(303,"(6F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(302,"(3F10.4,3I10.0,7F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(304,"(3F10.4,3I10.0,6F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do

	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)

	!deallocate(auxx,auxy,auxz,auyy,auyz,auzz)

	deallocate(nocpk)
	deallocate(stt)
	!deallocate(esp)
	deallocate(kpt)
	deallocate(vecres)





	close(200)



	close(300)
	close(301)
	close(303)	
	
	if (tmcoef) then
	close(302)
	close(304)
	else
	 continue
	end if
	

	


end subroutine spoptics
