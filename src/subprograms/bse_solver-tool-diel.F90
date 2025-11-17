
subroutine bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef,nocpf,fermishift,bsealgo,dft,mag)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	real,parameter:: pi=acos(-1.)

	integer :: dimbse !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	real,allocatable,dimension(:,:) :: eigv
	complex,allocatable,dimension(:,:,:) :: vector

	real,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:) :: stt

	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	!real,allocatable,dimension(:,:) :: lcount !pertencimento a camada de um dado estado

	!real,allocatable,dimension(:) :: orbweight

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	complex:: matrizelbse


	complex,allocatable,dimension(:) :: hrx,hry,hrz,hrsp,hrsm 

	real,allocatable,dimension(:) :: actxx,actyy,actzz,actxy,actxz,actyz,actsp,actsm

	complex,parameter :: imag=cmplx(0.0,1.0)

	real :: a,r0,ed,ec,ev,egap

	integer :: ngkpt
	real,allocatable,dimension(:,:) :: vecres
		
	!real,dimension(3) :: shift !shift no mhkpack

	!variaveis relacionadas a marcacao do tempo

	real:: t0,tf
	integer,dimension(8) :: values,values2

	integer,allocatable,dimension(:) :: nocpk
	
	integer :: nocpf
	real :: fermishift
	character(len=12) :: bsealgo
	character(len=1) :: dft
	real,dimension(3) :: mag	
	
	complex,allocatable,dimension(:,:,:) :: sk

	!definicoes diagonalizacao 
	INTEGER   ::       ifail
	integer,allocatable,dimension(:) :: ifailx
	real,parameter :: ABSTOL=1.0e-6
	INTEGER  ::        INFO
	real,allocatable,dimension(:) :: W,RWORK
	COMPLEX,allocatable,dimension(:,:) :: hbse

        INTEGER :: LWMAX
   	INTEGER :: LWORK
	INTEGER :: LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        complex,allocatable,dimension (:) :: WORK

        !zheevr definitions
        real :: VL,VU
        integer :: IL,IU,M
        COMPLEX,allocatable,dimension(:,:) :: Z
        integer :: LDZ
        integer,allocatable,dimension(:) :: ISUPPZ

	!modificacoes versao 2.1

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	real :: edos0,edosf,numdos
	real :: ebse0,ebsef,numbse
	real :: sme,exc,rk
	real,dimension(3) :: mshift
	real :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	!character(len=70) :: diein    !ambiente dieletrico
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	real,dimension(3) :: ediel
	logical :: bsewf
	integer :: excwf0,excwff
	real :: ez,w1,lc
	
	logical :: cpol,dtfull
	logical :: tmcoef	

	!fim modificacoes versao 2.1

	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada orb weight"

	!OUTPUT

	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse-diel output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bse_oscf.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_oscf.dat output file"
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"bse_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_oscf-pol.dat output file"	


	OPEN(UNIT=401, FILE=trim(outputfolder)//"ipa_oscf.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ipa_oscf.dat output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=402, FILE=trim(outputfolder)//"tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=404, FILE=trim(outputfolder)//"tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=403, FILE=trim(outputfolder)//"ipa_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ipa_oscf-pol.dat output file"	


	!OPEN(UNIT=500, FILE=trim(outputfolder)//"log_bse-matrix.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Error opening bse hamiltonian matrix progress output file"

	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)
!#ifdef gnu
!	call OPENBLAS_SET_NUM_THREADS(nthreads)
!#endif



	!inicio leitura parametros

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select

	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!dtfull = .false.
	!cpol = .true.


	write(301,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","xy","  ","xz","  ","yz"
	write(302,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","sp","  ","sm"	


	!termino parametros calculo 

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if


	!Informações para o arquivo de log do calculo
	write(300,*)
	write(300,*)
	write(300,*)
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,"(A14,1E15.4)") 'ktol-coulomb:', ktol
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	!write(300,*) 'angulo theta:',beta,'rad'

	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'Coulomb Potential:',coultype
	write(300,*)
	write(300,*) 'BSE_ALGO:',bsealgo

	!write(300,*)
	!write(300,*) 'intensidade termo exchange:',exc
	!write(300,*)
	!write(300,*) 'intensidade termo exchange l:',excl
	!write(300,*)
	!write(300,*) 'intensidade campo eletrico:',ez
	!write(300,*)
	!write(300,*) 'direção da magnetização do termo de exchange:',mag(1),'x',mag(2),'y',mag(3),'z'
	!write(300,*)

	!write(300,*) 'arquivo de saida 1:',logbse
	!write(300,*) 'arquivo de saida 2:',bseoptx
	!write(300,*) 'arquivo de saida 3:',bseopty
	!write(300,*) 'arquivo de saida 4:',bseoptsp
	!write(300,*) 'arquivo de saida 5:',bseoptsm

	
	write(300,*)
	write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 

	call flush(300)

	allocate(kpt(ngkpt,3))

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)


	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))
	allocate(nocpk(ngkpt))

	allocate (stt(ngkpt*nc*nv,4))
	allocate(hrx(dimbse),hry(dimbse),hrz(dimbse),hrsp(dimbse),hrsm(dimbse))

	allocate(hbse(dimbse,dimbse))
	


	allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	allocate(actxz(dimbse),actyz(dimbse))
	allocate(actsp(dimbse),actsm(dimbse))	
	
	
	
	
	select case (bsealgo)
		
		case ("cheev")
		
		LWMAX = 2*dimbse-1
				
		LWORK = -1

		allocate (RWORK(3*dimbse-2))
		allocate(W(dimbse),WORK(LWMAX))

		case ("cheevr")
		
		!LWMAX = (1 + 5*dimbse + 2*dimbse**2)
				
		LWORK = -1
      		LIWORK = -1
      		LRWORK = -1
		
		allocate(W(dimbse))
		allocate(WORK(2*dimbse))
		allocate(ISUPPZ(2*dimbse))
		allocate(IWORK(10*dimbse))
		allocate(RWORK(24*dimbse))
		allocate(Z(dimbse,dimbse)) 
		
		case ("cheevx")
		
		!LWMAX = (1 + 5*dimbse + 2*dimbse**2)
				
		LWORK = -1
 
		
		allocate(W(dimbse))
		allocate(WORK(2*dimbse))
		allocate(IWORK(5*dimbse))
		allocate(RWORK(7*dimbse))
		allocate(Z(dimbse,dimbse)) 
		allocate(ifailx(dimbse))
		
		case ("cheevd")
		
		LWMAX = (2*dimbse-1)**2
		
		LWORK = -1
      		LIWORK = -1
      		LRWORK = -1
      		!LWORK = 2*dimbse + dimbse**2
      		!LIWORK = 3 + 5*dimbse
      		!LRWORK = 1 + 5*dimbse + 2*dimbse**2
		allocate(W(dimbse),WORK(2*dimbse + dimbse**2))
		allocate (IWORK(3 + 5*dimbse))
		allocate (RWORK(1 + 5*dimbse + 2*dimbse**2))
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
	end select

	!allocate(orbweight(w90basis))
	
	!do i=1,w90basis
	!	read(203,*) orbweight(i)
	!end do

	!allocate(lcount(ngkpt,w90basis))

	egap = 50.0

	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux)
	do i=1,ngkpt

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,ovp,efermi,eaux,vaux,nocpf,fermishift,mag)
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	

			do j=1,nc+nv
	
				eigv(i,j)= eaux(nocpk(i)-nv+j)

			end do
			
			!if (ntype .eq. 2) then

			!do kl = 1,w90basis

			!	call layercont(w90basis,vaux(kl,:),orbweight,lcount(i,kl))


			!end do

			!else 

			!	continue

			!end if

			if ((eigv(i,nv+1)-eigv(i,nv)) .le. egap) then

				egap = eigv(i,nv+1)-eigv(i,nv)

			else

				continue

			end if

			


			do l=1,nc+nv


				do h=1,w90basis

					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)


				end do
			

			end do

	
	end do
	!$omp end parallel do

	deallocate(eaux,vaux)
	write(300,*) 'direct gap:', egap
	write(300,*) 'eigenvalues and eigenvectors calculated'
	call flush(300)
	
	!write(*,*) "autovetores e autovalores"


	if (dft .eq. "S") then
	
		allocate(sk(ngkpt,w90basis,w90basis))
		
		do i=1,ngkpt
		
			call overlap(w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),sk(i,:,:))
		
		end do
		
		write(300,*) 'overlap matrices calculated'
		call flush(300)	
	else
	
	end if


	!definindo os numeros quanticos dos estados


	!allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)



	!allocate(hrx(dimbse),hry(dimbse),hrz(dimbse))

	write(300,*) 'quantum numbers for exciton basis set finished'
	call flush(300)


	
	allocate(vecres(dimbse,15))	
	
	write(401,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(403,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	
	if (tmcoef) then	
	write(402,*) "number of kpoints:",ngkpt
	write(402,*) "number of conduction states",nc
	write(402,*) "number of valence states",nv
	write(402,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ","ec "," ","ev"," ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(404,*) "number of kpoints:",ngkpt
	write(404,*) "number of conduction states",nc
	write(404,*) "number of valence states",nv
	write(404,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ","ec "," ","ev","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
		      
	else
	
	 continue
	end if	


	 ! $omp parallel default(shared) private(i,w90basis,rlat,rvec,hopmatrices,ihopmatrices)
	
	 !$omp parallel do 
	do i=1,dimbse

		!ec = eigv(stt(i,4),stt(i,3))
		!ev = eigv(stt(i,4),stt(i,2))

		call optsp(eigv(stt(i,4),stt(i,2)),vector(stt(i,4),stt(i,2),:),&
		     eigv(stt(i,4),stt(i,3)),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hrx(i),hry(i),hrz(i))
		     
		     hrsp(i) = (hrx(i)+cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     hrsm(i) = (hrx(i)-cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hrx(i)*conjg(hrx(i)))		
		vecres(i,9) = real(hry(i)*conjg(hry(i)))		
		vecres(i,10) = real(hrz(i)*conjg(hrz(i)))
		vecres(i,11) = real(hrx(i)*conjg(hry(i)))
		vecres(i,12) = real(hrx(i)*conjg(hrz(i)))
		vecres(i,13) = real(hry(i)*conjg(hrz(i)))		 
		vecres(i,14) = real(hrsp(i)*conjg(hrsp(i)))
		vecres(i,15) = real(hrsm(i)*conjg(hrsm(i)))				     


	end do
	 !$omp end parallel do


	 ! $omp end parallel
	
	call Bubblem(7,15,vecres, dimbse)
	
	do i=1,dimbse
	
		write(401,"(7F15.6)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(403,"(6F15.6)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(402,"(3F10.6,3I10.0,9F10.6)")  vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),eigv(stt(i,4),stt(i,3)),eigv(stt(i,4),stt(i,2)),vecres(i,7),&
						  vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(404,"(3F10.6,3I10.0,8F10.6)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),eigv(stt(i,4),stt(i,3)),eigv(stt(i,4),stt(i,2)),vecres(i,7),vecres(i,8),&
						  vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do	

	write(300,*) 'IPA single particle optics finished'
	call flush(300)

	deallocate(vecres)
	
	!go to 789

	!allocate(hbse(dimbse,dimbse),W(dimbse))

	hbse=0.0
	!counter = 0
	!$omp parallel do 
!collapse(2)

	do i=1,dimbse



		do j=i,dimbse



  hbse(i,j)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt(i,:),eigv(stt(i,4)&
  	    ,stt(i,3)),eigv(stt(i,4),stt(i,2)),vector(stt(i,4)&
            ,stt(i,3),:) ,vector(stt(i,4),stt(i,2),:),kpt(stt(i,4),:),stt(j,:),eigv(stt(j,4),stt(j,3))&
  	    ,eigv(stt(j,4),stt(j,2)) &
            ,vector(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:),dft,nvec,rvec,&
            sk(stt(i,4),:,:),sk(stt(j,4),:,:))

	
		!write(500,*) "i",i,"/",dimbse,"        ","j",j,"/",dimbse
		!call flush(500)

		end do



	end do

	!$omp end parallel do

	write(300,*) 'exciton Hamiltonian matrix finished'
	call flush(300)	

	!call OMP_SET_NUM_THREADS(nthreads)

	select case (bsealgo)
		
		case ("cheev")

     	
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		!LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      		IF(INFO .GT. 0 ) THEN
        	WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	STOP
      		END IF   

		case ("cheevr")
		
		CALL CHEEVR( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, ISUPPZ, WORK, LWORK,RWORK, LRWORK, IWORK, LIWORK,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 24*dimbse, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 10*dimbse, IWORK( 1 ) )	
		CALL CHEEVR( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, ISUPPZ, WORK, LWORK,RWORK, LRWORK, IWORK, LIWORK,INFO )
       	        IF( INFO.GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF  		
		
		IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF  
         	
         	hbse = Z
         	
         	deallocate(Z)
         	
         	case ("cheevx")
		
		CALL CHEEVX( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, WORK, LWORK,RWORK, IWORK, ifailx,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )

		CALL CHEEVX( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, WORK, LWORK,RWORK, IWORK, ifailx,INFO )
       	        IF( INFO.GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF  		
		
		IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF  
         	
         	hbse = Z
         	
         	deallocate(Z)
		
		case ("cheevd")
		
		CALL CHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK, INFO )
      		LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 1 + 5*dimbse + 2*dimbse**2, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 3 + 5*dimbse, IWORK( 1 ) )
      		!LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		!LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
      		!LIWORK = MIN( LWMAX, IWORK( 1 ) )
      		
		CALL CHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK, INFO )
       	        IF(INFO .GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF  
		
         	IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF  		
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
	end select			
	!allocate(pinter(dimbse),pintra(dimbse))

	!if (ntype .eq. 2) then

	!do i=1,dimbse

	!	call excitonil(i,w90basis,nc,nv,ngrid*ngrid,stt,dimbse,hbse(:,i),lcount,pinter(i),pintra(i))

	!	write(306,*) W(i),pinter(i),pintra(i)
	!	call flush(306)

	!end do 

	!else 

	!	continue

	!end if 
	
	if (bsewf) then
	
	      	do i=excwf0,excwff
      	
      			call excwf(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,hbse(:,i))
      	
      		end do
	
	else
	
	 continue
	
	end if

	write(300,*) 'exciton Hamiltonian diagonalized'
	call flush(300)
	write(300,*) "exciton ground state",W(1)

	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	
	
	select case (bsealgo)
	
		case ("cheev")
		
		deallocate(RWORK,WORK)
		
		case ("cheevr")
		
		deallocate(WORK,ISUPPZ,IWORK,RWORK)
		
		case ("cheevx")
		
		deallocate(WORK,IWORK,RWORK,ifailx)
		
		case ("cheevd")		
		
		deallocate(WORK,IWORK,RWORK)

		case default 
		
		 write(*,*) "please select a subroutine for diagonalization"
		 go to 789	
		 	
	end select	
	!deallocate (IWORK)

	!allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	!allocate(actxz(dimbse),actyz(dimbse))
	

	
	call dielbsev(nthreads,dimbse,hbse,hrx,actxx)
	write(300,*) 'xx tensor component'
	call flush(300)
		
	call dielbsev(nthreads,dimbse,hbse,hry,actyy)
	write(300,*) 'yy tensor component'
	call flush(300)
		
	call dielbsev(nthreads,dimbse,hbse,hrz,actzz)
	write(300,*) 'zz tensor component'
	call flush(300)
	
	if (dtfull) then
	
	call dielbsep(nthreads,dimbse,hbse,hrx,hry,actxy)
	write(300,*) 'xy tensor component'
	call flush(300)
	
	call dielbsep(nthreads,dimbse,hbse,hrx,hrz,actxz)
	write(300,*) 'xz tensor component'
	call flush(300)
		
	call dielbsep(nthreads,dimbse,hbse,hry,hrz,actyz)
	write(300,*) 'yz tensor component'
	call flush(300)
	
	else
	
	actxy = 0.0000
	write(300,*) 'xy tensor component set to 0'
	call flush(300)
	
	actxz = 0.0000
	write(300,*) 'xz tensor component set to 0'
	call flush(300)
		
	actyz = 0.0000
	write(300,*) 'yz tensor component set to 0'
	call flush(300)
	
	end if
	
	if (cpol) then
	
	call dielbsev(nthreads,dimbse,hbse,hrsp,actsp)
	write(300,*) 'sp polarization'
	call flush(300)
		
	call dielbsev(nthreads,dimbse,hbse,hrsm,actsm)
	write(300,*) 'sm polarization'
	call flush(300)
	
	else 
	
	actsp = 0.0000
	write(300,*) 'sp polarization set to 0'
	call flush(300)
		
	actsm = 0.0000
	write(300,*) 'sm polarization set to 0'
	call flush(300)
	
	end if	


	write(300,*) 'optics finished'
	call flush(300)

	! $omp do ordered
	do i=1,dimbse
		! $omp ordered
		write(301,"(7F15.6)") W(i),actxx(i),actyy(i),actzz(i),actxy(i),actxz(i),actyz(i)
		call flush(301)
		
		if (cpol) then
		write(302,"(6F15.6)") W(i),actxx(i),actyy(i),actzz(i),actsp(i),actsm(i)
		call flush(302)		
		end if

		! $omp end ordered
	end do
	! $omp end do
	
	deallocate(hrx,hry,hrz,hrsp,hrsm)
	deallocate(actxx,actxy,actxz,actyy,actyz,actzz)
	deallocate(actsp,actsm)
	deallocate(hbse,W,stt,nocpk)
	deallocate(ovp)

	deallocate(kpt)
	
	if (dft .eq. "S") then
	
	deallocate(sk)
	
	else
	
	continue
	
	end if	


789     continue

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)




	close(200)
	!close(203)



	close(300)
	close(301)
	close(302)
	
	close(401)
	close(402)
	close(403)
	close(404)
	
	!close(500)			


	


end subroutine bsesolver
