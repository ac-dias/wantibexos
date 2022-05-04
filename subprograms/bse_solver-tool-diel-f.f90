
subroutine bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	double precision,parameter:: pi=acos(-1.)

	integer :: dimbse !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector

	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:) :: stt

	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	!double precision,allocatable,dimension(:,:) :: lcount !pertencimento a camada de um dado estado

	!double precision,allocatable,dimension(:) :: orbweight

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	double complex:: matrizelbse


	double complex,allocatable,dimension(:) :: hrx,hry,hrz,hrsp,hrsm 

	double precision,allocatable,dimension(:) :: actxx,actyy,actzz,actxy,actxz,actyz,actsp,actsm

	double complex,parameter :: imag=cmplx(0.0,1.0)

	double precision :: a,r0,ed,ec,ev,egap

	integer :: ngkpt
	real,allocatable,dimension(:,:) :: vecres
		
	!double precision,dimension(3) :: shift !shift no mhkpack

	!variaveis relacionadas a marcacao do tempo

	double precision:: t0,tf
	integer,dimension(8) :: values,values2

	integer,allocatable,dimension(:) :: nocpk

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
	logical :: bsewf
	integer :: excwf0,excwff
	double precision :: ez,w1,lc
	
	logical :: cpol,dtfull
	logical :: tmcoef	

	!fim modificacoes versao 2.1

	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada orb weight"

	!OUTPUT

	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse-diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse-diel output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bse_opt_diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel output file"
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"bse_opt_diel-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening output file"	


	OPEN(UNIT=401, FILE=trim(outputfolder)//"sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=402, FILE=trim(outputfolder)//"tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=404, FILE=trim(outputfolder)//"tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=403, FILE=trim(outputfolder)//"sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics-pol output file"		



	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)
	!call OPENBLAS_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)

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
	!write(300,*) 'constante dielétrica:',ed,'  ','r0:',r0

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

	allocate(hbse(dimbse,dimbse),W(dimbse))

	allocate (RWORK(3*dimbse-2))
        LWMAX = 2*dimbse-1
	allocate(WORK(LWMAX))

	allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	allocate(actxz(dimbse),actyz(dimbse))
	allocate(actsp(dimbse),actsm(dimbse))	

	!allocate(orbweight(w90basis))
	
	!do i=1,w90basis
	!	read(203,*) orbweight(i)
	!end do

	!allocate(lcount(ngkpt,w90basis))

	egap = 50.0

	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


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
	! $omp end parallel do

	deallocate(eaux,vaux)
	write(300,*) 'direct gap:', egap
	write(300,*) 'eigenvalues and eigenvectors calculated'
	call flush(300)


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
	write(402,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(404,*) "number of kpoints:",ngkpt
	write(404,*) "number of conduction states",nc
	write(404,*) "number of valence states",nv
	write(404,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
		      
	else
	
	 continue
	end if	

 
	!$omp parallel do &
	!$OMP DEFAULT (SHARED)	

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
 	
	call Bubblem(7,15,vecres, dimbse)
	
	do i=1,dimbse
	
		write(401,"(7F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(403,"(6F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(402,"(3F10.4,3I10.0,7F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(404,"(3F10.4,3I10.0,6F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do	

	write(300,*) 'single particle optics finished'
	call flush(300)

	deallocate(vecres)
	
	!go to 789

	!allocate(hbse(dimbse,dimbse),W(dimbse))

	hbse=0.

	!$omp parallel do 
!collapse(2)

	do i=1,dimbse



		do j=i,dimbse



  hbse(i,j)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt(i,:),eigv(stt(i,4)&
  	    ,stt(i,3)),eigv(stt(i,4),stt(i,2)),vector(stt(i,4)&
            ,stt(i,3),:) ,vector(stt(i,4),stt(i,2),:),kpt(stt(i,4),:),stt(j,:),eigv(stt(j,4),stt(j,3))&
  	    ,eigv(stt(j,4),stt(j,2)) &
            ,vector(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:))



		end do



	end do

	!$omp end parallel do

	write(300,*) 'exciton Hamiltonian matrix finished'
	call flush(300)	



      	!LWORK = 2*dimbse+dimbse**2
      	!LIWORK = 3 + 5*dimbse
      	!LRWORK = 1 + 5*dimbse + 2*dimbse**2


!	allocate(WORK(LWORK))
!	allocate (IWORK(LIWORK))
!	allocate (RWORK(LRWORK))

 !     	CALL ZHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 
!	      RWORK, LRWORK, IWORK, LIWORK,INFO )
 !     	IF( INFO.GT. 0 ) THEN
  !      WRITE(*,*)'The algorithm failed to compute eigenvalues.'
  !       STOP
  !    	END IF   

	!allocate (RWORK(3*dimbse-2))
        !LWMAX = 2*dimbse-1
	!allocate(WORK(LWMAX))

     	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF   

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
      	
      			call excwf2(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,hbse(:,i))
      	
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
	deallocate(RWORK,WORK)
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

	!$omp do ordered
	do i=1,dimbse
		!$omp ordered
		write(301,"(7F15.4)") W(i),actxx(i),actyy(i),actzz(i),actxy(i),actxz(i),actyz(i)
		call flush(301)
		
		if (cpol) then	
		write(302,"(6F15.4)") W(i),actxx(i),actyy(i),actzz(i),actsp(i),actsm(i)
		call flush(302)
		end if		


		!$omp end ordered
	end do
	!$omp end do
	
	deallocate(hrx,hry,hrz,hrsp,hrsm)
	deallocate(actxx,actxy,actxz,actyy,actyz,actzz)
	deallocate(actsp,actsm)

	deallocate(hbse,W,stt,nocpk)

	deallocate(kpt)

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


	


end subroutine bsesolver
