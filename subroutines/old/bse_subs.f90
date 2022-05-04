
function matrizelbse(coultype,tolr,w90basis,ediel,lc,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype

	integer,dimension(3) :: ngrid

	integer :: w90basis
	double precision :: a,vcell1

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


	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,tolr,lc)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc)

	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat)

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


subroutine optrpa2(ev,vv,ec,vc,kx,ky,kz,ffactor,sme,&
		  w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    exx,exy,exz,eyy,eyz,ezz)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	double precision :: ev,ec,erpa,sme,factor

	integer,dimension(w90basis) :: ffactor

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	double complex,dimension(w90basis,w90basis) :: hx,hy,hz

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(w90basis) :: vc,vv

	double complex :: hrx,hry,hrz

	double precision :: exx,exy,exz,eyy,eyz,ezz



	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)

	erpa = ec - ev

	
	call sandwich(w90basis,vc,hx,vv,hrx)

	call sandwich(w90basis,vc,hy,vv,hry)

	call sandwich(w90basis,vc,hz,vv,hrz)

	factor = 1.0 !/((erpa*erpa)+(sme*sme))

	exx = factor*(hrx*conjg(hrx))
	exy = factor*(hrx*conjg(hry))
	exz = factor*(hrx*conjg(hrz))

	eyy = factor*(hry*conjg(hry))
	eyz = factor*(hry*conjg(hrz))

	ezz = factor*(hrz*conjg(hrz))


end subroutine optrpa2

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

subroutine optrpabse(ev,vv,ec,vc,kx,ky,kz,ffactor,sme,&
		  w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    hrx,hry,hrz)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	double precision :: ev,ec,sme,factor

	integer,dimension(w90basis) :: ffactor

	double complex :: erpa

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	double complex,dimension(w90basis,w90basis) :: hx,hy,hz

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(w90basis) :: vc,vv

	double complex :: hrx,hry,hrz

	double precision :: exx,exy,exz,eyy,eyz,ezz



	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)

	erpa = ec - ev + cmplx(0.0,sme)

	
	call sandwich(w90basis,vc,hx,vv,hrx)

	call sandwich(w90basis,vc,hy,vv,hry)

	call sandwich(w90basis,vc,hz,vv,hrz)

	factor = 1./erpa

	hrx = factor*hrx
	hry = factor*hry
	hrz = factor*hrz



end subroutine optrpabse

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
	double precision :: actaux
	double complex,dimension(dimse) :: hopt1,hopt2
	double complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k,no
	

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0


	! $OMP DO PRIVATE(actaux)
	!$OMP PARALLEL DO PRIVATE(actaux)
	do i=1,dimse


		actaux=0.0


		do j=1,dimse
			do k=1,dimse


		actaux=actaux+(excitonvec(j,i)*hopt1(j)*conjg(excitonvec(k,i))*conjg(hopt2(k)))

		end do
			end do

		activity(i)=actaux

	
	
		


	end do
	!$OMP END PARALLEL DO



end subroutine dielbsep

subroutine dielbsev(dimse,excitonvec,hopt,activity)

	implicit none


	integer :: dimse
	double precision,dimension(dimse) :: activity, exciton
	double complex :: actaux
	double complex,dimension(dimse) :: hopt
	double complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k
	

	activity=0.0
	actaux=0.0

	!$omp do private(actaux)
	do i=1,dimse

		actaux=0.0

		do j=1,dimse
			


		actaux=actaux+(excitonvec(j,i)*hopt(j))
	

		end do
		

		activity(i)=actaux*conjg(actaux)

	


	end do
	!$omp end do


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
    	if (erro/=0) stop "Erro na abertura do arquivo de saida wf BSE"
    	
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





























