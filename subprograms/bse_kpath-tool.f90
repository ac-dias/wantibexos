!gfortran -mcmodel=large bse_kpath-tool.f90 -o bse_kpath-tool.x -llapack95 -lopenblas -fopenmp 

subroutine bsebnds(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	double precision,parameter:: pi=acos(-1.)
	integer :: dimbse
	double precision,allocatable,dimension(:,:) :: exk !(dimbse,nqpts*npath)


	integer :: ncaux,nvaux

	integer,allocatable,dimension(:,:) :: stt !(ngrid*ngrid*nc*nv,4)
	double precision,allocatable,dimension(:,:) :: kpt,qpt !pontos k do grid (ngrid*ngrid,2)

	double precision,dimension(4) :: q

	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	double precision,allocatable,dimension(:,:) :: energy !variavel que guarda os autovalores para cada ponto k
	double complex,allocatable,dimension(:,:,:) :: vector !variavel que guarda os autovetores para cada ponto k

	double precision,allocatable,dimension(:,:) :: energyq !variavel que guarda os autovalores para cada ponto k
	double complex,allocatable,dimension(:,:,:) :: vectorq !variavel que guarda os autovetores para cada ponto k

	double precision,allocatable,dimension(:,:) :: qauxv !(nqpts*npath,4) 

	integer:: counter,c,v,i,j,f,l,h,erro,i2,k

	double complex:: matrizelbsekq

	double precision :: a,r0,ed

	integer :: ngkpt

	!variaveis relacionadas a marcacao do tempo

	double precision:: t0,tf
	integer,dimension(8) :: values,values2

	!variaveis kpoints

	integer :: nks
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks

	integer,allocatable,dimension(:) :: nocpk,nocpq


	!definicoes diagonalizacao 
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
	COMPLEX*16,allocatable,dimension(:,:) :: hbse

        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
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
	double precision :: ez,w1,lc
	logical :: bsewf
	integer :: excwf0,excwff	
	

	!call input_read

	! INPUT 
	OPEN(UNIT=500, FILE= kpathsbse,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse-kpath input file"
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse_kpath output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"bands_bse.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands_bse output file"


	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
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
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,*) 'ktol-coulomb:', ktol
	write(300,*)
	write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'number of kpoints in the path:','   ',(nks/2)*nkpts
	write(300,*)
	write(300,*)
	write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
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

	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)

			do j=1,nc+nv
	
				energy(i,j)= eaux(nocpk(i)-nv+j)

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


    !write(*,*) "quantum numbers for exciton basis set finished"

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


	! $omp parallel do
		do i2=1,ngkpt

            !write(*,*) i2

			call eigsys(nthreads,scs,exc,nocpq(i2),ffactor,qpt(i2,1),qpt(i2,2),qpt(i2,3),w90basis,nvec,&
			             rlat,rvec,hopmatrices,ihopmatrices,efermi,eaux,vaux)

				do j=1,nc+nv
	
					energyq(i2,j)= eaux(nocpq(i2)-nv+j) 

				end do

	
				do l=1,nc+nv


					do h=1,w90basis

						vectorq(i2,l,h)=vaux(nocpq(i2)-nv+l,h)

	
					end do
			

				end do



	
		end do
	! $omp end parallel do

		call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

		deallocate(eaux,vaux)


		allocate(hbse(dimbse,dimbse),W(dimbse))

		hbse=0.

	!$omp parallel do 
        !collapse(2)

		do i2=1,dimbse



			do j=i2,dimbse




hbse(i2,j)= matrizelbsekq(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,q,rlat,stt(i2,:),energyq(stt(i2,4),stt(i2,3))&
          ,energy(stt(i2,4),stt(i2,2)),vectorq(stt(i2,4)&
          ,stt(i2,3),:) ,vector(stt(i2,4),stt(i2,2),:),kpt(stt(i2,4),:),stt(j,:)&
          ,energyq(stt(j,4),stt(j,3)),energy(stt(j,4),stt(j,2))&
          ,vectorq(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:))


			end do



		end do

	!$omp end parallel do

    		!call LA_HEEVR( hbse, W, JOBZ='N', UPLO='U', ABSTOL=ABSTOL, INFO=INFO ) 



      	!LWORK = 2*dimbse+dimbse**2
      	!LIWORK = 3 + 5*dimbse
      	!LRWORK = 1 + 5*dimbse + 2*dimbse**2

	!allocate(WORK(LWORK))
	!allocate (IWORK(LIWORK))
	!allocate (RWORK(LRWORK))

 

      	!CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
      	!LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	!LRWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	!LIWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      	!CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,&
	!	      RWORK, LRWORK, IWORK, LIWORK,INFO )
      	!IF( INFO.GT. 0 ) THEN
        !WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        ! STOP
      	!END IF   


	allocate (RWORK(3*dimbse-2))
        LWMAX = 2*dimbse-1
	allocate(WORK(LWMAX))
	
	if (bsewf) then
	
	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF 
	
	
	else

     	LWORK = -1
      	CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF 
      	
      	end if  

		!write(*,*) W(1)


		do i2=1,dimbse

			exk(i2,i)=W(i2)

		!	write(400,*) qauxv(i,1),W(i2)
		!	call flush(400)

		end do


	if (bsewf) then
	
	      	do i2=excwf0,excwff
      	
      			call excwfi2(outputfolder,ngkpt,kpt,q,nc,nv,nocpk,stt,W(i2),i2,i,hbse(:,i2))
      	
      		end do
	
	else
	
	 continue
	
	end if

		deallocate(hbse,W)
		deallocate(energyq,vectorq)
		deallocate(stt)
		deallocate(RWORK,WORK)
		!deallocate (IWORK)


		qpt=0.0
		nocpq=0

		write(300,*) 'progress:',i,'/',(nks/2)*nkpts
		call flush(300)

	end do



	do i2=1,dimbse

		do i=1,(nks/2)*nkpts

			write(400,*) real(qauxv(i,1)),exk(i2,i)
			call flush(400)

		end do

		write(400,*)

	end do



	deallocate(energy,vector)
	deallocate(kpt)
	deallocate(qpt)
	deallocate(exk)
	deallocate(nocpq,nocpk)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)

	close(200)



	close(300)
	close(400)
	close(500)



end subroutine bsebnds
