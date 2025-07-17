

!gfortran bands-kpath-tool.f90 -o bands_tool.x -llapack95 -lopenblas -fopenmp

subroutine bandstool(nthreads,outputfolder,params,kpaths,orbw, &
		     exc,mag,mshift,dft,nocpf,fermishift)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,k,erro,j2
	real,parameter :: pi=acos(-1.)

	!complex,allocatable,dimension(:,:) :: htb
	real,allocatable,dimension(:) :: eigv

	real,allocatable,dimension(:,:,:) :: ebands

	complex,allocatable,dimension(:,:) :: autovetores

	real,allocatable,dimension(:,:) :: kpts

	real,allocatable,dimension(:) :: orbweight

	real :: kx,ky,kz,kp

	real :: lco

	complex :: spz,spx,spy

	!variaveis kpoints

	integer :: nks
	integer :: nkpts
	real,allocatable,dimension(:,:) :: ks

	!call input_read

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
	real,dimension(3) :: ediel,mag
	character(len=1) :: dft	
	
	integer,allocatable,dimension(:) :: nocpj
	integer :: nkc,nkv,nkgap
	real :: cbm,vbm,gap
	
	integer :: nocpf
	real :: fermishift
	complex,allocatable,dimension(:,:) :: ovptb
	
	!integer :: power
	!logical :: lowdin

	!fim modificacoes versao 2.1

	! INPUT : lendo os parametros do modelo de tight-binding
	
	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file"


	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"bands.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bands_kpath_coord.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands kpath coord output file"    	
    	OPEN(UNIT=500, FILE=trim(outputfolder)//"cbm_vbm-kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening cbm-vbm kpath output file"
    	OPEN(UNIT=501, FILE=trim(outputfolder)//"gap-kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gap kpath output file"    	

	!OPEN(UNIT=500, FILE=trim(outputfolder)//"nocp.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida nocp"

	call OMP_SET_NUM_THREADS(nthreads)
	
	
	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select


	

	read(202,*) nks
	read(202,*) nkpts

	allocate(ks(nks,3))

	do i=1,nks

		read(202,*) ks(i,1),ks(i,2),ks(i,3)
	
	end do

	allocate(orbweight(w90basis))

	OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	if (erro/=0) then
	
	  orbweight= 1.0

	else

	  do i=1,w90basis
		read(203,*) orbweight(i)
	  end do

	end if



	!termino leitura parametros

	!parametros do calculo

	!OPEN(UNIT=2077, FILE= trim(calcparms)//"bands_calc.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	!call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
	!	     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
	!	     exc,mshift,coultype)

	!close(2077)

	!termino parametros calculo
 

	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	allocate(ebands(w90basis,((nks/2)*nkpts),9))
	!allocate(kpts(nkpts*(nks-1),4))
	allocate(kpts((nks/2)*nkpts,4))
	allocate(nocpj((nks/2)*nkpts))

	!definindo kpath

	call kpath(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpts)

	!termino definicao kpath

	ebands=0.0
	
	
	allocate(ovptb(w90basis,w90basis))
	
		!write(*,*) "before main loop"
	
	
	!$omp parallel do default(shared) private(j,i,kp,kx,ky,kz,eigv,autovetores,ovptb,spx,spy,spz)	
	do j=1,(nks/2)*nkpts

	        kp= kpts(j,1)
		kx= kpts(j,2)
		ky= kpts(j,3)
		kz= kpts(j,4)
		

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
		    
		!call eigsysl(nthreads,dft,systype,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		 !   ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag,lowdin,power,ovptb,smh)	    
		
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif	

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	
    
		   
		    

		!write(500,*) nocp,kx,ky,kz 
		
		
		if (dft .eq. "S") then
		 call overlap(w90basis,nvec,rvec,ovp,kx,ky,kz,ovptb)
		else
		 continue
		end if
		
		
		do i=1,w90basis

			call layercont(w90basis,dft,autovetores(i,:),orbweight,ovptb,lco)
			call spinvl(w90basis,autovetores(i,:),dft,systype,ovptb,spx,spy,spz)

			
			ebands(i,j,1) = kp
			ebands(i,j,2) = eigv(i)
			ebands(i,j,3) = real(spx)
			ebands(i,j,4) = real(spy)
			ebands(i,j,5) = real(spz)						
			ebands(i,j,6) = lco
			ebands(i,j,7) = kx
			ebands(i,j,8) = ky
			ebands(i,j,9) = kz									


		end do
		




	end do
	!$omp end parallel do
	
	
	!write(*,*) "autovetores e autovalores"
	
	write(500,*) "2"
	write(500,*) "C"
	

	
	call gapfinder(w90basis,((nks/2)*nkpts),nocpj,ebands(:,:,2),nkc,nkv,nkgap,gap,cbm,vbm)

	write(500,"(3E18.8,I7,1F18.8)") kpts(nkv,2),kpts(nkv,3),kpts(nkv,4),nocpj(nkv),vbm
	write(500,"(3E18.8,I7,1F18.8)") kpts(nkc,2),kpts(nkc,3),kpts(nkc,4),nocpj(nkc)+1,cbm
	
	
	write(501,*) "fundamental band gap (eV):",cbm-vbm
	write(501,*) "direct band gap (eV):",gap
	write(501,*)
	write(501,*)
	write(501,*)"kpoint cbm"
	write(501,"(3E18.8)")kpts(nkc,2),kpts(nkc,3),kpts(nkc,4)
	write(501,*)	
	write(501,*)"kpoint vbm"
	write(501,"(3E18.8)")kpts(nkv,2),kpts(nkv,3),kpts(nkv,4)
	write(501,*)		 
	write(501,*)"kpoint direct band gap"
	write(501,"(3E18.8)")kpts(nkgap,2),kpts(nkgap,3),kpts(nkgap,4)


	write(300,*) "#k-path energy <sx> <sy> <sz> orbital"	
	write(301,*) "#kx ky kz k-path energy "	
	
	do i=1,w90basis

			write(300,*) "#band",i
			write(301,*) "#band",i			

		
		do j=1,(nks/2)*nkpts

	               write(300,"(1E20.8,1F20.8,4F20.6)") ebands(i,j,1),ebands(i,j,2),ebands(i,j,3),&
	                            ebands(i,j,4),ebands(i,j,5),ebands(i,j,6)
	                            
	               write(301,"(4E20.8,1F20.8)") ebands(i,j,7),ebands(i,j,8),ebands(i,j,9),&
	                            ebands(i,j,1),ebands(i,j,2)	                            

		end do

			write(300,*)
			write(301,*)			

	end do



	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	
	!if (dft .eq. "S") then
	deallocate(ovp,ovptb)
	!else
	! continue
	!end if

	deallocate(eigv)
	deallocate(kpts)
	deallocate(orbweight)
	deallocate(autovetores)
	deallocate(ebands)
	deallocate(ks)
	deallocate(nocpj)


	close(200)
	close(300)
	close(301)
	close(202)
	close(203)
	close(500)
	close(501)




end subroutine bandstool
