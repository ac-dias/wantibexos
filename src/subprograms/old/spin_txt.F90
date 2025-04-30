subroutine spintexture(nthreads,outputfolder,ngrid,params,mshift,dft,nocpf,fermishift)

	use omp_lib
	!use hamiltonian_input_variables
	
	implicit none
	
	integer :: i,j,k,ngkpt,erro
	double complex :: spz,spx,spy
	
	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores
	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid	
	
	double precision,allocatable,dimension(:,:,:) :: ebnd	
	
	character(len=200) :: file1,file2,file3	
	
	
	!data from hamiltonian input variables
	
	integer :: w90basis,ntype,nocp
	double precision :: efermi,scs
	!double precision:: edielh
	integer :: nvec
	double precision,dimension(3,3) :: rlat

	integer,allocatable,dimension(:,:) :: rvec
	double precision,allocatable,dimension(:,:,:) :: hopmatrices,ihopmatrices
	integer,allocatable,dimension(:) :: ffactor

	character(len=4) :: systype

	!double precision :: lc
	
	!variables for spin polarized hamiltonian
	double precision,allocatable,dimension(:,:,:) :: hopmatricesu,ihopmatricesu
	double precision,allocatable,dimension(:,:,:) :: hopmatricesd,ihopmatricesd
	integer,allocatable,dimension(:) :: ffactoru,ffactord
	integer,allocatable,dimension(:,:) :: rvecu,rvecd
	integer :: w90basisu,w90basisd
	integer :: nvecu,nvecd		
	
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
	character(len=1) :: dft
	double precision,dimension(3) :: ediel
	integer :: nocpf
	double precision :: fermishift
	
				
	call OMP_SET_NUM_THREADS(nthreads)
	
	call hamiltonian_input_read(200,params)	
	
	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)	
	
	allocate(kpt(ngkpt,3))
	
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)
	
	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(ebnd(ngkpt,w90basis,7))
	
	do i=1,ngkpt
	
		call eigsys(nthreads,scs,exc,nocp,ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux,nocpf,fermishift)
		             
		             
		do j=1,w90basis
		
			call spinvl(w90basis,vaux(j,:),dft,systype,spx,spy,spz)
			
			ebnd(i,j,1) = kpt(i,1)
			ebnd(i,j,2) = kpt(i,2)
			ebnd(i,j,3) = kpt(i,3)						
			ebnd(i,j,4) = eaux(j)
			ebnd(i,j,5) = real(spx)
			ebnd(i,j,6) = real(spy)
			ebnd(i,j,7) = real(spz)						
		
		end do             	
	
	
	end do				

	deallocate(eaux,vaux,kpt)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	
	do i=1,w90basis
	
	 write(file1,"(a18)")'spin_texture_band_'
	 write(file2,"(I0)") i
	 write(file3,"(a4)")'.dat'
	 	 	 
	 OPEN(UNIT=700+i, FILE=trim(outputfolder)//trim(file1)//trim(file2)//trim(file3),STATUS='unknown',IOSTAT=erro)
    	 if (erro/=0) stop "Error opening spin texture band bz output file"
    	 
    	 write(700+i,*) "#kx ky kz energy <sx> <sy> <sz>"
    	 
    	  do j=1,ngkpt
    	  
    	  	write(700+i,"(7F15.4)") (ebnd(j,i,k), k=1,7)
    	  
    	  end do
    	 
    	 close(700+i)

	end do	
	
	
	deallocate(ebnd)

		
	
end subroutine spintexture	
