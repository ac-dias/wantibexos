
function matrizelbse(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,dft,nvec,rvec,sk,skp) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	integer :: w90basis,nvec
	real :: a,vcell1
	real :: ez,w
	
	real,dimension(nvec,3) :: rvec
	!real,dimension(nvec,w90basis,w90basis) :: ovp
	
	complex,dimension(w90basis,w90basis) :: sk,skp

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	complex, dimension(w90basis) :: vbc,vbv

	real,dimension(3,3) :: rlat

	real :: tolr
	integer :: ktol
	
	complex:: matrizelbse

	real,parameter:: pi=acos(-1.)

	real :: auxi

	real :: modk

	real,dimension(3) :: ediel

	real :: lc

	complex :: vc,vv

	real :: vcoul1

	real :: vcoul,v2dk,v3diel,v2dt,v0dt,v2dt2
	real :: v2dohono,v2drk,v1dt,v2d,v2diel
	real :: v1d,v1diel,v3davg,v2dtavg
	
	real :: r0


	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2,rlat,ngrid,tolr)
		
	case("V3DA")

		vcoul1= v3davg(kpt1,kpt2,ngrid,rlat,tolr)		

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vcoul1= v2d(kpt1,kpt2,rlat,ngrid,tolr)

	case("V2DL")

		vcoul1= v2diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)		

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DTA")

		vcoul1= v2dtavg(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoul1= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoul1= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")

		vcoul1= v1d(kpt1,kpt2,ngrid,rlat,lc,tolr)

	case("V1DL")

		vcoul1= v1diel(kpt1,kpt2,ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vcoul1= v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
				

	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	



	if (est1(1) .eq. est2(1)) then


		matrizelbse= (ec1-ev1) + vcoul1


	else

			
		select case (dft)
		
		case ("S")
		

		 !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		 !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		 
		 call sandwich(w90basis,vbc1,0.5*(sk+skp),vbc2,vc)
		 call sandwich(w90basis,vbv1,0.5*(sk+skp),vbv2,vv)
		 
		 matrizelbse=  vcoul1*vc*vv		
		
		case default
	
		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)
	
	         matrizelbse=  vcoul1*vc*vv
	         
		end select
	
		


	end if
		



end function matrizelbse




subroutine opticalactivity(dimse,excitonvec,hopt,activity,description)

	implicit none


	integer :: dimse
	real,dimension(dimse) :: activity
	complex :: actaux
	real :: description(dimse)
	complex,dimension(dimse) :: hopt
	complex,dimension(dimse,dimse) :: excitonvec

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

	complex,dimension(ndim) :: wf

	real,dimension(ncount,2*w90basis) :: lc

	real :: pinter,pintra

	real :: inter,intra

	real :: p1c,p1v,p2c,p2v


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
	real,dimension(dimse) :: activity, exciton
	real :: actaux
	complex,dimension(dimse) :: hopt1,hopt2
	complex,dimension(dimse,dimse) :: excitonvec

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
	real,dimension(dimse) :: activity, exciton
	complex :: actaux,actaux2
	complex,dimension(dimse) :: hopt1,hopt2
	complex,dimension(dimse,dimse) :: excitonvec
	integer,dimension(dimse,4) :: stt

	integer :: i,j,k,no
	

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0


	! $OMP DO PRIVATE(actaux)
	! $OMP PARALLEL DO PRIVATE(actaux,actaux2)
	do i=1,dimse


		actaux=0.0
		actaux2=0.0


		do j=1,dimse

			
			actaux=actaux+(excitonvec(j,i)*hopt1(j))
			actaux2=actaux2+(excitonvec(j,i)*hopt2(j))
			

			
		end do
			

		activity(i)=actaux*conjg(actaux2)

	
	
		


	end do
	! $OMP END PARALLEL DO



end subroutine dielbsep



subroutine dielbsev(nthread,dimse,excitonvec,hopt,activity)

	implicit none


	integer :: dimse,nthread
	real,dimension(dimse) :: activity, exciton
	complex :: actaux
	complex,dimension(dimse) :: hopt
	complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k

	call OMP_SET_NUM_THREADS(nthread)	

	activity=0.0
	actaux=0.0

	! $OMP PARALLEL DO PRIVATE(actaux)
	do i=1,dimse

		actaux=0.0

		do j=1,dimse
			


		actaux=actaux+(excitonvec(j,i)*hopt(j))
	

		end do
		

		activity(i)=actaux*conjg(actaux)

	


	end do
	! $OMP END PARALLEL DO


end subroutine dielbsev



subroutine exclft(sysdim,ngrid,rlat,fosc,enexc,lft) !tempo de vida exciton

	implicit none
	integer,dimension(3) :: ngrid
	real :: lft,fosc,refr,enexc
	real,dimension(3,3) :: rlat
	
	real,parameter :: cfe= 7.297E-3 !constante estrutura fina
	real,parameter :: lspd= 299.792E16 !vel luz Ang/s^2
	real,parameter :: hbar= 6.582E-16 !const planck eV.s
	
	real,parameter :: cic= -(0.0904756)*10**(3)
	
	real,parameter:: pi=acos(-1.)	
	
	real :: aux1,aux2,aux3,vc
	
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



subroutine excwf(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf) !exciton wf direto

	implicit none
	
	integer :: i,erro,counter,j
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1,file2
	CHARACTER(LEN=30) :: Format,Format2
	integer :: ngkpt,nc,nv,excnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	real :: excenergy
	real,dimension(ngkpt,3) :: kpt
	complex,dimension(ngkpt*nc*nv) :: ewf
	real,dimension(ngkpt) :: den
	
	
	WRITE (file1,"(a7,I0,a4)") 'exc_wf_',excnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	WRITE (file2,"(a8,I0,a4)") 'exc_den_',excnum,'.dat'
	Format2 = "(3F15.4,1E15.4)"	
	
	OPEN(UNIT=700+excnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	OPEN(UNIT=800+excnum, FILE=trim(outputfolder)//trim(file2),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening den BSE output file"
    	
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
	
	write(800+excnum,*) "#"," ","excitonic state",excnum
    	write(800+excnum,*) "#"," ","exciton energy:",excenergy
    	write(800+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(800+excnum,*) "#"," ","Number of conduction states","",nc
    	write(800+excnum,*) "#"," ","Number of valence states","",nv
    	write(800+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(800+excnum,*) "#"," ","kx"," ","ky"," ","kz"," ","den"
    	
    	counter=0
    	
    	do i=1,ngkpt
    	
    		den(i) = 0.0
    	
    		do j=1,nc*nv
    		
    			counter=counter+1
    			
    			den(i) = den(i)+ewf(counter)*conjg(ewf(counter)) 
    		
    			
    		
    		end do
    		
    		 write(800+excnum,Format2) kpt(i,1),kpt(i,2),kpt(i,3),den(i)
    	
    	end do
    	
    
	close(700+excnum)
	close(800+excnum)

end subroutine excwf

subroutine excwfi(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf) !exciton wf indireto

	implicit none
	
	integer :: i,erro,counter,j
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1,file2
	CHARACTER(LEN=30) :: Format,Format2
	integer :: ngkpt,nc,nv,excnum,qptnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	real :: excenergy
	real,dimension(ngkpt,3) :: kpt
	real,dimension(4) :: qpt
	complex,dimension(ngkpt*nc*nv) :: ewf
	real,dimension(ngkpt) :: den
	
	
	WRITE (file1,"(a7,I0,a1,I0,a4)") 'exc_wf_',excnum,"_",qptnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	WRITE (file2,"(a8,I0,a1,I0,a4)") 'exc_den_',excnum,"_",qptnum,'.dat'
	Format = "(3F15.4,1E15.4)"	
	
	OPEN(UNIT=700+excnum*qptnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	OPEN(UNIT=800+excnum*qptnum, FILE=trim(outputfolder)//trim(file2),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening den BSE output file"    	
    	
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
	
	write(800+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(800+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(800+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(800+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(800+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(800+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(800+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(800+excnum*qptnum,*) "#"," ","kx"," ","ky"," ","kz"," ","den"
    	
    	counter=0
    	
    	do i=1,ngkpt
    	
    		den(i) = 0.0
    	
    		do j=1,nc*nv
    		
    			counter=counter+1
    			
    			den(i) = den(i)+ewf(counter)*conjg(ewf(counter)) 
    		
    			
    		
    		end do
    		
    		 write(800+excnum*qptnum,Format2) kpt(i,1),kpt(i,2),kpt(i,3),den(i)
    	
    	end do    	
    	
    
	close(700+excnum*qptnum)
	close(800+excnum*qptnum)

end subroutine excwfi





















