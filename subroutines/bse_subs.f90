
function matrizelbse(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype

	integer,dimension(3) :: ngrid

	integer :: w90basis
	double precision :: a,vcell1
	double precision :: ez,w

	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	double complex, dimension(w90basis) :: vbc,vbv

	double precision,dimension(3,3) :: rlat

	double precision :: tolr
	integer :: ktol
	
	double complex:: matrizelbse

	double precision,parameter:: pi=acos(-1.)

	double precision :: auxi

	double precision :: modk

	double precision,dimension(3) :: ediel

	double precision :: lc

	double complex :: vc,vv

	double precision :: vcoul1

	double precision :: vcoul,v2dk,v3diel,v2dt,v0dt,v2dt2
	double precision :: v2dohono,v2drk
	
	double precision :: r0


	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2,rlat,ngrid,tolr)

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoul1= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoul1= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	



	if (est1(1) .eq. est2(1)) then


		matrizelbse= (ec1-ev1) + vcoul1


	else

	
		call vecconjg(vbc1,w90basis,vbc)

		call vecconjg(vbv1,w90basis,vbv)

		call prodintsq(vbc,vbc2,w90basis,vc)

		call prodintsq(vbv,vbv2,w90basis,vv)
	
		matrizelbse=  vcoul1*vc*vv


	end if
		



end function matrizelbse




subroutine opticalactivity(dimse,excitonvec,hopt,activity,description)

	implicit none


	integer :: dimse
	double precision,dimension(dimse) :: activity
	double complex :: actaux
	double precision :: description(dimse)
	double complex,dimension(dimse) :: hopt
	double complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j
	

	activity=0.
	actaux=cmplx(0.,0.)

	do i=1,dimse

		do j=1,dimse

		
		actaux=actaux+(excitonvec(j,i)*hopt(j))

		end do

		activity(i)=actaux*conjg(actaux)

		if (activity(i) .gt. 0.1) then

			description(i)= 1.

		else if ( (activity(i) .lt. 0.1) .and. (activity(i) .gt. 1.0e-8)) then

			description(i)= 0.

		else

			description(i)= -1.


		end if
	
		actaux=cmplx(0.,0.)


	end do



end subroutine opticalactivity

subroutine excitonil(i,w90basis,nc,nv,ncount,stt,ndim,wf,lc,pinter,pintra)

	implicit none

	integer :: i !numero do estado excitonico

	integer :: j,k,nc,nv

	integer :: w90basis
	integer :: ncount,ndim
	integer,dimension(ncount*nc*nv,4) :: stt

	double complex,dimension(ndim) :: wf

	double precision,dimension(ncount,2*w90basis) :: lc

	double precision :: pinter,pintra

	double precision :: inter,intra

	double precision :: p1c,p1v,p2c,p2v


	pinter = 0d0
	pintra = 0d0

	do j=1,ndim

		p1c = 1d0-lc(stt(j,4),stt(j,3))
		p1v = 1d0-lc(stt(j,4),stt(j,2))

		p2c = lc(stt(j,4),stt(j,3))
		p2v = lc(stt(j,4),stt(j,2))


		intra = (p1c*p1v)+(p2c*p2v)

		inter = (p1c*p2v)+(p2c*p1v)


		pinter = pinter + Real(conjg(wf(j))*wf(j)*inter)
 
		pintra = pintra + Real(conjg(wf(j))*wf(j)*intra)

	end do


end subroutine excitonil



subroutine dielbse(dimse,excitonvec,hopt1,hopt2,activity)

	implicit none


	integer :: dimse
	double precision,dimension(dimse) :: activity, exciton
	double precision :: actaux
	double complex,dimension(dimse) :: hopt1,hopt2
	double complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k
	

	activity=0.0
	actaux=0.0

	do i=1,dimse



		do j=1,dimse
			do k=1,dimse


		actaux=actaux+(excitonvec(j,i)*hopt1(j)*conjg(excitonvec(k,i))*conjg(hopt2(k)))

		end do
			end do

		activity(i)=actaux

	
	
		actaux=0.0


	end do



end subroutine dielbse



subroutine dielbsep(nthread,dimse,excitonvec,hopt1,hopt2,activity)

	use omp_lib
	implicit none


	integer :: dimse,nthread
	double precision,dimension(dimse) :: activity, exciton
	double complex :: actaux,actaux2
	double complex,dimension(dimse) :: hopt1,hopt2
	double complex,dimension(dimse,dimse) :: excitonvec
	integer,dimension(dimse,4) :: stt

	integer :: i,j,k,no
	

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0


	! $OMP DO PRIVATE(actaux)
	!$OMP PARALLEL DO PRIVATE(actaux,actaux2)
	do i=1,dimse


		actaux=0.0
		actaux2=0.0


		do j=1,dimse

			
			actaux=actaux+(excitonvec(j,i)*hopt1(j))
			actaux2=actaux2+(excitonvec(j,i)*hopt2(j))
			

			
		end do
			

		activity(i)=actaux*conjg(actaux2)

	
	
		


	end do
	!$OMP END PARALLEL DO



end subroutine dielbsep



subroutine dielbsev(nthread,dimse,excitonvec,hopt,activity)

	implicit none


	integer :: dimse,nthread
	double precision,dimension(dimse) :: activity, exciton
	double complex :: actaux
	double complex,dimension(dimse) :: hopt
	double complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k

	call OMP_SET_NUM_THREADS(nthread)	

	activity=0.0
	actaux=0.0

	!$OMP PARALLEL DO PRIVATE(actaux)
	do i=1,dimse

		actaux=0.0

		do j=1,dimse
			


		actaux=actaux+(excitonvec(j,i)*hopt(j))
	

		end do
		

		activity(i)=actaux*conjg(actaux)

	


	end do
	!$OMP END PARALLEL DO


end subroutine dielbsev


subroutine excwf(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf)

	implicit none
	
	integer :: i,erro
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double complex,dimension(ngkpt*nc*nv) :: ewf
	
	
	WRITE (file1,'(a7,I0,a4)') ,'exc_wf_',excnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
		
 	write(700+excnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  stt(i,3),stt(i,2),real(ewf(i)),aimag(ewf(i))
 
	end do
    	
    
	close(700+excnum)

end subroutine excwf

subroutine excwfi(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf)

	implicit none
	
	integer :: i,erro
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum,qptnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double precision,dimension(4) :: qpt
	double complex,dimension(ngkpt*nc*nv) :: ewf
	
	
	WRITE (file1,'(a7,I0,a1,I0,a4)') ,'exc_wf_',excnum,"_",qptnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum*qptnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(700+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum*qptnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
		
 	write(700+excnum*qptnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  stt(i,3),stt(i,2),real(ewf(i)),aimag(ewf(i))
 
	end do
    	
    
	close(700+excnum*qptnum)

end subroutine excwfi


subroutine exclft(sysdim,ngrid,rlat,fosc,enexc,lft) 

	implicit none
	integer,dimension(3) :: ngrid
	double precision :: lft,fosc,refr,enexc
	double precision,dimension(3,3) :: rlat
	
	double precision,parameter :: cfe= 7.297E-3 !constante estrutura fina
	double precision,parameter :: lspd= 299.792E16 !vel luz Ang/s^2
	double precision,parameter :: hbar= 6.582E-16 !const planck eV.s
	
	double precision,parameter :: cic= -(0.0904756)*10**(3)
	
	double precision,parameter:: pi=acos(-1.)	
	
	double precision :: aux1,aux2,aux3,vc
	
	character(len=5) :: sysdim
	

	
	select case (sysdim)
	
	 case("3D")
	 
	 	!for 3D systems- DOI: 10.1103/PhysRevLett.95.247402
		aux1 = (4.0*cfe)*(enexc*enexc*enexc)
		aux2 = fosc/dble(ngrid(1)*ngrid(2)*ngrid(3))
		aux3 = 3.0*(lspd*lspd)*(hbar*hbar*hbar)
		lft = (aux1*aux2)/aux3
		lft = (1.0/lft)
	 
	 case("2D")
	 
	 	 !for 2D systems- DOI:10.1021/nl503799t
	 	call vcell2D(rlat,vc)
		aux1 = dble(ngrid(1)*ngrid(2)*ngrid(3))*vc*hbar
		aux2 = (8.0*pi*cfe*enexc)*(fosc)	
		lft = aux1/aux2
	 
	 case("1D")
	 
	 	!for 1D systems- DOI: 10.1103/PhysRevLett.95.247402
		aux2 = (2.0*pi*cfe*enexc*enexc)*fosc
		aux1 = dble(ngrid(1)*ngrid(2)*ngrid(3))*rlat(1,1)*hbar*hbar*lspd
		lft = aux1/aux2
	 
	 case default

		write(*,*) "Wrong System dimension"
		STOP

	end select 
	

	

	
end subroutine exclft



subroutine excwf2(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf)

	implicit none
	
	integer :: i,erro
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double complex,dimension(ngkpt*nc*nv) :: ewf
	
	
	WRITE (file1,'(a7,I0,a4)') ,'exc_wf_',excnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
		
 	write(700+excnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  nocp(stt(i,4))+stt(i,3)-nv,nocp(stt(i,4))-nv+stt(i,2),real(ewf(i)),aimag(ewf(i))
 
	end do
    	
    
	close(700+excnum)

end subroutine excwf2

subroutine excwfi2(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf)

	implicit none
	
	integer :: i,erro
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum,qptnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double precision,dimension(4) :: qpt
	double complex,dimension(ngkpt*nc*nv) :: ewf
	
	
	WRITE (file1,'(a7,I0,a1,I0,a4)') ,'exc_wf_',excnum,"_",qptnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum*qptnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(700+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum*qptnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
		
 	write(700+excnum*qptnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  nocp(stt(i,4))+stt(i,3)-nv,nocp(stt(i,4))-nv+stt(i,2),real(ewf(i)),aimag(ewf(i))
 
	end do
    	
    
	close(700+excnum*qptnum)

end subroutine excwfi2





















