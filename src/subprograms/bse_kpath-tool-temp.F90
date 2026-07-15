!gfortran -mcmodel=large bse_kpath-tool.f90 -o bse_kpath-tool.x -llapack95 -lopenblas -fopenmp 

subroutine bsebndstemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     st,phavg,ta,temp,nocpf,fermishift,bsealgo,dft,mag)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	real,parameter:: pi=acos(-1.)
	integer :: dimbse
	real,allocatable,dimension(:,:) :: exk !(dimbse,nqpts*npath)


	integer :: ncaux,nvaux

	integer,allocatable,dimension(:,:) :: stt !(ngrid*ngrid*nc*nv,4)
	real,allocatable,dimension(:,:) :: kpt,qpt !pontos k do grid (ngrid*ngrid,2)

	real,dimension(4) :: q

	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	real,allocatable,dimension(:,:) :: energy !variavel que guarda os autovalores para cada ponto k
	complex,allocatable,dimension(:,:,:) :: vector !variavel que guarda os autovetores para cada ponto k

	real,allocatable,dimension(:,:) :: energyq !variavel que guarda os autovalores para cada ponto k
	complex,allocatable,dimension(:,:,:) :: vectorq !variavel que guarda os autovetores para cada ponto k

	real,allocatable,dimension(:,:) :: qauxv !(nqpts*npath,4) 

	integer:: counter,c,v,i,j,f,l,h,erro,i2,k

	complex:: matrizelbsekqtemp

	real :: a,r0,ed

	integer :: ngkpt

	!variaveis relacionadas a marcacao do tempo

	real:: t0,tf
	integer,dimension(8) :: values,values2

	!variaveis kpoints

	integer :: nks
	integer :: nkpts
	real,allocatable,dimension(:,:) :: ks

	integer,allocatable,dimension(:) :: nocpk,nocpq
	
	integer :: nocpf
	real :: fermishift
	character(len=12) :: bsealgo
	character(len=1) :: dft
	
	real,dimension(3) :: mag
	complex,allocatable,dimension(:,:,:) :: sk,skq		


	!definicoes diagonalizacao 
	INTEGER   ::       ifail
	integer,allocatable,dimension(:) :: ifailx
	real,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	real,allocatable,dimension(:) :: W,RWORK
	COMPLEX,allocatable,dimension(:,:) :: hbse

        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        complex,allocatable,dimension (:) :: WORK
        
        !zheevr definitions
        real :: VL,VU
        integer :: IL,IU,M
        complex,allocatable,dimension(:,:) :: Z
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
	character(len=2) :: ta
	real,dimension(3) :: ediel
	real :: ez,w1,lc
	logical :: bsewf
	integer :: excwf0,excwff
	
	real :: st,phavg,temp
	real :: gapcortemp,gapcortemp2,tcor
		
	real,allocatable,dimension(:) :: gwcor
	integer :: gwin
	real :: gwaux
	character(len=1) :: gwchar	

	!call input_read

	! INPUT 
	OPEN(UNIT=500, FILE= kpathsbse,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse-kpath input file"
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"gw_qp_energy_cor_avg.dat",STATUS='old', IOSTAT=gwin) 

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse_kpath output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"bands_bse.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands_bse output file"


	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select
	!ediel(2) = edielh

	read(500,*) nks
	read(500,*) nkpts

	allocate(ks(nks,3))

	do i=1,nks

		read(500,*) ks(i,1),ks(i,2),ks(i,3)
	
	end do



	!termino leitura parametros


	!parametros do calculo


	!termino parametros calculo 




	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)
	!lc = rlat3(3)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

	!definindo kpath
	!allocate(qauxv(nkpts*(nks-1),4))
	allocate(qauxv((nks/2)*nkpts,4))	

	call kpathbse(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,qauxv)

	!Informações para o arquivo de log do calculo

	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'temperature', temp
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,*) 'ktol-coulomb:', ktol
	write(300,*)
	write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'Coulomb Potential:',coultype
	write(300,*)
	write(300,*) 'BSE_ALGO:',bsealgo	
	write(300,*)
	write(300,*) 'number of kpoints in the path:','   ',(nks/2)*nkpts
	write(300,*)
	write(300,*)
	write(300,*) 'begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 

	call flush(300)

	allocate(kpt(ngkpt,3))
	
	!shift= 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)
	!call gridgenmhp(ngrid,rlat,kpt)

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(energy(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))

	allocate(nocpk(ngkpt))
	allocate(nocpq(ngkpt))

	select case (ta)
	
	case("FA")
	
		tcor = 0.0
	
	case("VE")
	
		tcor = gapcortemp(st,phavg,temp)
	
	case("BE")
	
		tcor = gapcortemp2(st,phavg,temp)
	
	case default
	
		tcor= 0.00
	
	end select
	
    	if (gwin .eq. 0) then
    	
    		allocate(gwcor(w90basis))
    		
    			read(304,*) gwchar
    		   		
    		do i=1,w90basis
    		
    			read(304,*) gwcor(i),gwaux,gwaux,gwaux,gwaux
    			
    		
    		end do
    		
    		write(2077,*) "G0W0 correction applied in BSE kpath"
    	
    	else
    	
    		continue
    	
    	end if  	

	!$omp parallel do default(shared) private(i,eaux,vaux,j,l,h)
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

		call eigsys(nthreads,dft,systype,scs+tcor,exc,nocpk(i),&
			    ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
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
	
				if (gwin .eq. 0) then
				
				energy(i,j)= eaux(nocpk(i)-nv+j)+gwcor(nocpk(i)-nv+j)
				
				else
	
				energy(i,j)= eaux(nocpk(i)-nv+j)

				end if

			end do


			do l=1,nc+nv


				do h=1,w90basis

					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)


				end do
			

			end do

	
	end do
	!$omp end parallel do

	deallocate(eaux,vaux)

	if (dft .eq. "S") then
	
		allocate(sk(ngkpt,w90basis,w90basis))
		
		do i=1,ngkpt
		
			call overlap(w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),sk(i,:,:))
		
		end do
		!write(300,*) 'overlap matrices calculated'
		!call flush(300)	
	
	end if	

	!definindo os numeros quanticos dos estados




	write(300,*) "quantum numbers for exciton basis set finished"
	call flush(300)	

	counter=counter-1 !numero total de estados para equação bse

	allocate(qpt(ngkpt,3))
	allocate(exk(dimbse,nkpts*(nks-1)))



	do i=1,(nks/2)*nkpts

		
		!definindo os pontos q

		q(1)= qauxv(i,1)
		q(2)= qauxv(i,2) 
		q(3)= qauxv(i,3)
		q(4)= qauxv(i,4)

		!gerando o grid k+q
		
		call monhkhorst_packq(q,ngrid(1),ngrid(2),ngrid(3),mshift,&
				      rlat(1,:),rlat(2,:),rlat(3,:),qpt)


		allocate(eaux(w90basis),vaux(w90basis,w90basis))

		allocate(energyq(ngkpt,nc+nv),vectorq(ngkpt,nc+nv,w90basis))

	allocate (stt(ngkpt*nc*nv,4))


	!$omp parallel do default(shared) private(i2,j,l,h,eaux,vaux)
		do i2=1,ngkpt


#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif	




			call eigsys(nthreads,dft,systype,scs+tcor,exc,nocpq(i2),&
			            ffactor,qpt(i2,1),qpt(i2,2),qpt(i2,3),w90basis,nvec,&
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
	
					if (gwin .eq. 0) then
				
					energyq(i2,j)= eaux(nocpq(i2)-nv+j)+gwcor(nocpq(i2)-nv+j)
				
					else
	
					energyq(i2,j)= eaux(nocpq(i2)-nv+j) 

					end if	

				end do

	
				do l=1,nc+nv


					do h=1,w90basis

						vectorq(i2,l,h)=vaux(nocpq(i2)-nv+l,h)

	
					end do
			

				end do



	
		end do
	!$omp end parallel do
	
	if (dft .eq. "S") then
	
		allocate(skq(ngkpt,w90basis,w90basis))
		
		 
		do i2=1,ngkpt
		
			call overlap(w90basis,nvec,rvec,ovp,qpt(i2,1),qpt(i2,2),qpt(i2,3),skq(i2,:,:))
		
		end do
		write(300,*) 'overlap matrices calculated'
		call flush(300)	
	
	
	end if	

		call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

		deallocate(eaux,vaux)



		allocate(hbse(dimbse,dimbse),W(dimbse))

		hbse=0.

	!$omp parallel do 
        !collapse(2)

		do i2=1,dimbse



			do j=i2,dimbse




hbse(i2,j)= matrizelbsekqtemp(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,q,rlat,stt(i2,:),energyq(stt(i2,4),stt(i2,3))&
          ,energy(stt(i2,4),stt(i2,2)),vectorq(stt(i2,4)&
          ,stt(i2,3),:) ,vector(stt(i2,4),stt(i2,2),:),kpt(stt(i2,4),:),stt(j,:)&
          ,energyq(stt(j,4),stt(j,3)),energy(stt(j,4),stt(j,2))&
          ,vectorq(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:),temp,dft,nvec,rvec,&
          sk(stt(i2,4),:,:),sk(stt(j,4),:,:),skq(stt(i2,4),:,:),skq(stt(j,4),:,:))

			end do



		end do

	!$omp end parallel do

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
	
	if (bsewf) then
	
	select case (bsealgo)
		
		case ("cheev")

     	
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      		IF( INFO.GT. 0 ) THEN
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
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
      		LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 1 + 5*dimbse + 2*dimbse**2, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 3 + 5*dimbse, IWORK( 1 ) )
		CALL CHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
       	        IF( INFO.GT. 0 ) THEN
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
	
	
	else

 	select case (bsealgo)
		
		case ("cheev")

     	
      		CALL CHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		CALL CHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      		IF( INFO.GT. 0 ) THEN
        	WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	STOP
      		END IF   

		case ("cheevr")
		
		CALL CHEEVR( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, ISUPPZ, WORK, LWORK,RWORK, LRWORK, IWORK, LIWORK,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 24*dimbse, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 10*dimbse, IWORK( 1 ) )	
		CALL CHEEVR( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
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
		
		CALL CHEEVX( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, WORK, LWORK,RWORK, IWORK, ifailx,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )

		CALL CHEEVX( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
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
		
		CALL CHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
      		LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 1 + 5*dimbse + 2*dimbse**2, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 3 + 5*dimbse, IWORK( 1 ) )
		CALL CHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
       	        IF( INFO.GT. 0 ) THEN
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
      	
      	end if  

		!write(*,*) W(1)


		do i2=1,dimbse

			exk(i2,i)=W(i2)

		!	write(400,*) qauxv(i,1),W(i2)
		!	call flush(400)

		end do


	if (bsewf) then
	
	      	do i2=excwf0,excwff
      	
      			call excwfi(outputfolder,ngkpt,kpt,q,nc,nv,nocpk,stt,W(i2),i2,i,hbse(:,i2))
      	
      		end do
	
	else
	
	 continue
	
	end if

		deallocate(hbse,W)
		deallocate(energyq,vectorq)
		deallocate(stt)

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
		 go to 580	
		 	
	end select
		!deallocate (IWORK)


		qpt=0.0
		nocpq=0

		write(300,*) 'progress:',i,'/',(nks/2)*nkpts
		call flush(300)
		
		if (dft .eq. "S") then
		
			deallocate(skq)
		
		end if		

	end do


	if (gwin .eq. 0) then
    	
    		deallocate(gwcor)
        end if

	do i2=1,dimbse

		do i=1,(nks/2)*nkpts

			write(400,*) real(qauxv(i,1)),exk(i2,i)
			call flush(400)

		end do

		write(400,*)

	end do

580 continue

	deallocate(energy,vector)
	deallocate(kpt)
	deallocate(qpt)
	deallocate(exk)
	deallocate(nocpq,nocpk)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ovp,ffactor)
	
		if (dft .eq. "S") then
		
			deallocate(sk)
		
		end if		

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)

	close(200)



	close(300)
	close(400)
	close(500)

	close(304)

end subroutine bsebndstemp
