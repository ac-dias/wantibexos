
!gfortran tdos-tool.f90 -o tdos-tool.x -llapack95 -lopenblas -fopenmp

subroutine dostool(nthreads,outputfolder,ngrid,numdos, &
		     sme,params,orbw,exc,mag,mshift,dft,nocpf,fermishift,spintxt)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	integer :: i,j,k,erro
	real,parameter :: pi=acos(-1.)
	complex,parameter :: imag=cmplx(0.0,1.0)

	integer :: ngkpt
	!real,dimension(3) :: shift !shift no mhkpack

	complex,allocatable,dimension(:,:) :: autovetores
	real,allocatable,dimension(:) :: eigv

	real,allocatable,dimension(:,:) :: kpts

	real,allocatable,dimension(:,:) :: res
        real,allocatable,dimension(:,:,:) :: res2

	real,allocatable,dimension(:) :: orbweight

	real :: kx,ky,kz,kp


	integer  :: ndim,counter

	!variaveis para dos

	real :: tdos,wdos
	real,allocatable,dimension(:)   :: atdos,awdos

	real  :: cup,cdown
	real,allocatable,dimension(:)  :: updos,dndos,en
	
	real :: deltaen
	real :: sumdos

	real :: spinf
	
	complex :: spinz,spinx,spiny	

	!variaveis para grafico dos

	real :: at0

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
	character(len=1) :: dft
	real,dimension(3) :: ediel,mag
	integer :: nocpf
	real :: fermishift
	
	logical :: spintxt
	character(len=100) :: f1,f2,f3
	
	integer,allocatable,dimension(:) :: nocpj
	integer :: nkc,nkv,nkgap
	real :: cbm,vbm,gap
	real,allocatable,dimension(:,:) :: ebands
	complex,allocatable,dimension(:,:) :: ovptb
	
	

	!fim modificacoes versao 2.1

	!call input_read


	! INPUT : lendo os parametros do modelo de tight-binding

	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"



	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"dos.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening dos output file"
    	OPEN(UNIT=500, FILE=trim(outputfolder)//"cbm_vbm-kmesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening cbm-vbm kmesh output file"
    	OPEN(UNIT=501, FILE=trim(outputfolder)//"gap-kmesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gap kmesh output file"      	


	call OMP_SET_NUM_THREADS(nthreads)
	!call OPENBLAS_SET_NUM_THREADS(nthread)

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)

	end select

	!ediel(2) = edielh
	
	
	if ( systype .eq. "NP" ) then

	spinf = 2.0
	write(300,*) "#energy tdos sumdos awdos"

	else

	spinf = 1.0
	write(300,*) "#energy tdos sumdos updos dndos awdos"
	
	end if
	

	allocate(orbweight(w90basis))

	OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)

    	if (erro/=0) then
	
	orbweight= 1.0

	else

	do i=1,w90basis
		read(203,*) orbweight(i)
	end do

	end if




	allocate(eigv(w90basis))
	allocate(autovetores(w90basis,w90basis))

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	allocate(kpts(ngkpt,3))

	allocate(ovptb(w90basis,w90basis))
	
	!parametros do calculo

	!OPEN(UNIT=2077, FILE= trim(calcparms)//"dos_calc.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	!call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
	!	     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
	!	     exc,mshift,coultype)

	!close(2077)

	!termino parametros calculo

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpts)


	counter=0

	ndim=ngkpt*w90basis
	
	allocate(res(ndim,5))
	allocate(nocpj(ngkpt))
	allocate(ebands(w90basis,ngkpt))
	
	if (spintxt) then

	allocate(res2(ngkpt,w90basis,7))



	end if

	edos0= 0.0
	edosf= 0.0

	!$omp parallel do default(shared) private(j,i,kx,ky,kz,eigv,autovetores,ovptb,spinx,spiny,spinz,cup,cdown,tdos)
	do j=1,ngkpt

		kx= kpts(j,1)
		ky= kpts(j,2)
		kz= kpts(j,3)

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif	

#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif	

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif		

		call eigsys(nthreads,dft,systype,scs,exc,nocpj(j),ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)

#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif	

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	  	      
		       		    
		    

		if (eigv(1) .lt. edos0 ) then
		
			edos0 = eigv(1)
			
		end if
		
		
		if (eigv(w90basis) .gt. edosf ) then
		
			edosf = eigv(w90basis)
			
		end if		


		if (dft .eq. "S") then
		 call overlap(w90basis,nvec,rvec,ovp,kx,ky,kz,ovptb)
		else
		 continue
		end if
		
		

		do i=1,w90basis

			
			tdos = 1.0
			call spincont(w90basis,autovetores(i,:),cup,cdown,dft,systype,ovptb)
			call layercont(w90basis,dft,autovetores(i,:),orbweight,ovptb,wdos)
			
 
			
			!counter=counter+1

			res(i+(j-1)*w90basis,1) = eigv(i)
			res(i+(j-1)*w90basis,2) = tdos
			res(i+(j-1)*w90basis,3) = cup
			res(i+(j-1)*w90basis,4) = cdown
			res(i+(j-1)*w90basis,5) = wdos
			!res(counter,4) = dcont
			!res(counter,5) = pcont
			!res(counter,6) = 1.0
			ebands(i,j) = eigv(i)
			
			
			
			if (spintxt) then
			
			call spinvl(w90basis,autovetores(i,:),dft,systype,ovptb,spinx,spiny,spinz)
			
			res2(j,i,1) = kx
			res2(j,i,2) = ky
			res2(j,i,3) = kz					
			res2(j,i,4) = eigv(i)
			res2(j,i,5) = real(spinx)
			res2(j,i,6) = real(spiny)
			res2(j,i,7) = real(spinz)				
			
			else
			
			continue
			
			end if


		end do
		



	end do
	!$omp end parallel do
	

	call gapfinder(w90basis,ngkpt,nocpj,ebands,nkc,nkv,nkgap,gap,cbm,vbm)
	
	write(500,*) "2"
	write(500,*) "C"	
	write(500,"(3E15.8,I7,1F15.8)") kpts(nkv,2),kpts(nkv,3),kpts(nkv,4),nocpj(nkv),vbm
	write(500,"(3E15.8,I7,1F15.8)") kpts(nkc,2),kpts(nkc,3),kpts(nkc,4),nocpj(nkc)+1,cbm
	
	write(501,*) "fundamental band gap (eV):",cbm-vbm
	write(501,*) "direct band gap (eV):",gap
	write(501,*)
	write(501,*)
	write(501,*)"kpoint cbm"
	write(501,"(3E15.8)")kpts(nkc,2),kpts(nkc,3),kpts(nkc,4)
	write(501,*)	
	write(501,*)"kpoint vbm"
	write(501,"(3E15.8)")kpts(nkv,2),kpts(nkv,3),kpts(nkv,4)
	write(501,*)		 
	write(501,*)"kpoint direct band gap"
	write(501,"(3E15.8)")kpts(nkgap,2),kpts(nkgap,3),kpts(nkgap,4)
	
	deallocate(ebands)		
	
	if (spintxt) then
		
	do i=1,w90basis
	
	 write(f1,"(a18)")'spin_texture_band_'
	 write(f2,"(I0)") i
	 write(f3,"(a4)")'.dat'
	 	 	 
	 OPEN(UNIT=700+i, FILE=trim(outputfolder)//trim(f1)//trim(f2)//trim(f3),STATUS='unknown',IOSTAT=erro)
    	 if (erro/=0) stop "Error opening spin texture band bz output file"
    	 
    	 write(700+i,*) "#kx ky kz energy <sx> <sy> <sz>"
    	 
    	  do j=1,ngkpt
    	  
    	  	write(700+i,"(7F15.4)") (res2(j,i,k), k=1,7)
    	  
    	  end do
    	 
    	 close(700+i)

	end do	
	
	deallocate(res2)
	
	else
	
	continue
	
	end if	
	
	allocate(atdos(int(numdos)+1),updos(int(numdos)),dndos(int(numdos)),awdos(int(numdos)))
	allocate(en(int(numdos)))
	
	deltaen= dble(((edosf-edos0)))/(numdos-1.)
	
	sumdos = 0.0
	atdos = 0.0
	
	!$omp parallel do default(shared) private(i)
	do i=1,int(numdos)

		en(i) = edos0 + dble(((edosf-edos0)*(i-1))/(numdos-1.))


		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,2),sme,atdos(i))
		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,3),sme,updos(i))
		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,4),sme,dndos(i))
		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,5),sme,awdos(i))
		!call endos(a0,ndim,ngrid,ngrid,1,en,res(:,1),res(:,6),sme,tcont)

		!sumdos = sumdos+ deltaen*(at0)

		! $omp ordered
		!write(300,"(6F15.4)") en,atdos,sumdos,updos,-dndos,awdos
		!write(301,*) en,d2cont,p2cont,tcont
		! $omp end ordered

	end do
	!$omp end parallel do
	

	
	atdos = spinf*atdos
	awdos = spinf*awdos
	
	do i=1,int(numdos)
	
		if ( systype .eq. "NP" ) then
	
		write(300,"(4F15.4)") en(i),atdos(i),sumdos,awdos(i)
		
		else
		
		write(300,"(6F15.4)") en(i),atdos(i),sumdos,updos(i),-dndos(i),awdos(i)	

		end if
		
		if (i .eq. 1) then
		
		sumdos= sumdos 
		
		else
		
		sumdos= sumdos + (deltaen*(atdos(i-1)+atdos(i))*0.5)
		
		end if
	
	end do

	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(ovp,ovptb)

	deallocate(eigv)
	deallocate(autovetores)
	deallocate(kpts)
	deallocate(res)
	deallocate(orbweight)
	deallocate(atdos,updos,dndos,awdos)
	deallocate(nocpj)


	close(200)

	close(300)
	!close(301)
	close(500)
	close(501)





end subroutine dostool

