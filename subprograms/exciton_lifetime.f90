
subroutine lifetime(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)

	implicit none
	
	integer :: nc,nv
	double precision :: numbse
	integer,dimension(3) :: ngrid
	double precision,dimension(3,3) :: rlat
	character(len=70) :: outputfolder
	character(len=70) :: cflag
	character(len=5) :: sysdim	
	
	integer :: erro
	double precision,allocatable,dimension(:,:) :: refr
	double precision,allocatable,dimension(:,:) :: fosc
	integer :: dimbse,i,j
	
	double precision :: resx,resy,resz
	double precision :: refraux,foscaux
	double precision,dimension(4) :: lft
	
	
	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_opt_diel.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel input file"
	!OPEN(UNIT=544, FILE=trim(outputfolder)//"bse_indice_refracao.dat",STATUS='old', IOSTAT=erro)
	!if (erro/=0) stop "Erro na abertura do arquivo de entrada indice de refracao"
	
	!arquivo de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"exciton_lifetime.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening exciton_lifetime output file"
	
	
	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv	
	
	allocate(fosc(dimbse,7),refr(int(numbse),7))
	
	!lendo os inputs
	read(100,*) cflag
	do i=1,dimbse
	
	 read(100,*) fosc(i,1),fosc(i,2),fosc(i,3),fosc(i,4),fosc(i,5),fosc(i,6),fosc(i,7)

	end do
	
	!read(544,*) cflag
	
	
	!do i=1,int(numbse)
	
	 !read(544,*) refr(i,1),refr(i,2),refr(i,3),refr(i,4),refr(i,5),refr(i,6),refr(i,7)

	 
	!end do
	
	write(300,*) "#energy(eV) sum-lft(s)  x-lft(s) y-lft(s) z-lft(s)"
	
	do i=1,dimbse
	
		!call interp1D(int(numbse),7,2,refr,fosc(i,1),resx)
		!call interp1D(int(numbse),7,3,refr,fosc(i,1),resy)	
		!call interp1D(int(numbse),7,4,refr,fosc(i,1),resz)	
		

		foscaux = fosc(i,2)+fosc(i,3)+fosc(i,4)
		
		call exclft(sysdim,ngrid,rlat,foscaux,fosc(i,1),lft(1))
		call exclft(sysdim,ngrid,rlat,fosc(i,2),fosc(i,1),lft(2))
		call exclft(sysdim,ngrid,rlat,fosc(i,3),fosc(i,1),lft(3))	
		call exclft(sysdim,ngrid,rlat,fosc(i,4),fosc(i,1),lft(4))		
		!write(*,*) refraux
		
		write(300,"(1F15.4,4E15.4)") fosc(i,1),lft(1),lft(2),lft(3),lft(4)	
	
	end do
	
	
	
	
	
	
	
	
	deallocate(fosc,refr)

	close(100)
	!close(544)
	close(300)

end subroutine lifetime


subroutine lifetimepol(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)

	implicit none
	
	integer :: nc,nv
	double precision :: numbse
	integer,dimension(3) :: ngrid
	double precision,dimension(3,3) :: rlat
	character(len=70) :: outputfolder
	character(len=70) :: cflag
	character(len=5) :: sysdim	
	
	integer :: erro
	double precision,allocatable,dimension(:,:) :: refr
	double precision,allocatable,dimension(:,:) :: fosc
	integer :: dimbse,i,j
	
	double precision :: resx,resy,resz
	double precision :: refraux,foscaux
	double precision,dimension(6) :: lft
	
	
	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_opt_diel-pol.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol input file"
	!OPEN(UNIT=544, FILE=trim(outputfolder)//"bse_indice_refracao.dat",STATUS='old', IOSTAT=erro)
	!if (erro/=0) stop "Erro na abertura do arquivo de entrada indice de refracao"
	
	!arquivo de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"exciton_lifetime-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening exciton_lifetime-pol output file"
	
	
	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv	
	
	allocate(fosc(dimbse,6),refr(int(numbse),6))
	
	!lendo os inputs
	read(100,*) cflag
	do i=1,dimbse
	
	 read(100,*) fosc(i,1),fosc(i,2),fosc(i,3),fosc(i,4),fosc(i,5),fosc(i,6)

	end do
	
	!read(544,*) cflag
	
	
	!do i=1,int(numbse)
	
	 !read(544,*) refr(i,1),refr(i,2),refr(i,3),refr(i,4),refr(i,5),refr(i,6),refr(i,7)

	 
	!end do
	
	write(300,*) "#energy(eV) sum-lft(s)  x-lft(s) y-lft(s) z-lft(s) sp-lft(s) sm-lft(s)"
	
	do i=1,dimbse
	
		!call interp1D(int(numbse),7,2,refr,fosc(i,1),resx)
		!call interp1D(int(numbse),7,3,refr,fosc(i,1),resy)	
		!call interp1D(int(numbse),7,4,refr,fosc(i,1),resz)	
		

		foscaux = fosc(i,2)+fosc(i,3)+fosc(i,4)
		
		call exclft(sysdim,ngrid,rlat,foscaux,fosc(i,1),lft(1))
		call exclft(sysdim,ngrid,rlat,fosc(i,2),fosc(i,1),lft(2))
		call exclft(sysdim,ngrid,rlat,fosc(i,3),fosc(i,1),lft(3))	
		call exclft(sysdim,ngrid,rlat,fosc(i,4),fosc(i,1),lft(4))
		call exclft(sysdim,ngrid,rlat,fosc(i,5),fosc(i,1),lft(5))
		call exclft(sysdim,ngrid,rlat,fosc(i,6),fosc(i,1),lft(6))					
		!write(*,*) refraux
		
		write(300,"(1F15.4,6E15.4)") fosc(i,1),lft(1),lft(2),lft(3),lft(4),lft(5),lft(6)	
	
	end do
	
	
	
	
	
	
	
	
	deallocate(fosc,refr)

	close(100)
	!close(544)
	close(300)

end subroutine lifetimepol

