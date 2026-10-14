! Independent lattice-gauge reference for the DFT=S interband optical vertex.
program test_nonorth_optics
  use bse_q_optics, only: rmn_data, rmn_dfts_dipoles, rmn_destroy
  implicit none
  integer, parameter :: n=2, nr=3
  real :: rv(nr,3), lat(3,3), h(nr,n,n), hi(nr,n,n), s(nr,n,n)
  real :: tau(n,3), dm(nr,n,n,3), dout(nr,n,n,3), e(n), mag(3), k, gap, eta, mu
  complex :: c(n,n), sk(n,n), dh(n,n), ds(n,n), xk(n,n), dk(n,n)
  complex :: z(3), bz(3), full, mid, dcv, cond, ph, norm(n,n)
  integer :: ff(nr), i, a, b, axis, ie, nocc
  type(rmn_data) :: data
  logical :: ok
  character(len=256) :: message
  real, parameter :: tol=2.e-6, pi=acos(-1.), hbar=6.582119569e-16
  ! Sheet optics: incidence/substrate indices, thickness (m), conductivities (s^-1).
  real, parameter :: nis(4)=[1.,1.5,1.33,2.], nss(4)=[1.47,1.,1.47,3.5], thick=1.252e-9
  complex, parameter :: sigs(4)=[(0.,0.),(0.,3.e15),(0.,-3.e15),(2.e15,1.e15)]
  real :: tr, rf, tr2, rf2
  complex :: y
  integer :: j

  rv=0.; rv(2,:)=[1.,0.3,-0.2]; rv(3,:)=-rv(2,:)
  lat=0.; do i=1,3; lat(i,i)=1.; end do
  h=0.; hi=0.; s=0.; tau=0.; dm=0.; mag=0.; ff=1; k=0.73
  h(1,1,1)=-1.; h(1,2,2)=1.6; h(1,1,2)=0.35; h(1,2,1)=0.35
  h(2,1,2)=0.4; h(2,2,1)=0.1; h(3,:,:)=transpose(h(2,:,:))
  s(1,1,1)=1.; s(1,2,2)=1.; s(1,1,2)=0.08; s(1,2,1)=0.08
  s(2,1,1)=0.02; s(2,2,2)=0.015; s(2,1,2)=0.07; s(2,2,1)=0.01
  s(3,:,:)=transpose(s(2,:,:)); tau(2,:)=[0.32,-0.12,0.08]
  do axis=1,3
    dm(1,1,2,axis)=0.12*axis; dm(1,2,1,axis)=0.12*axis
    dm(1,1,1,axis)=0.03; dm(1,2,2,axis)=-0.02
  end do
  call eigsys(1,'S','NP  ',0.,0.,nocc,ff,k,0.,0.,n,nr,lat,rv,h,hi,s,0.,e,c,1,0.,mag)
  gap=e(2)-e(1)
  call overlap(n,nr,rv,s,k,0.,0.,sk)
  norm=matmul(conjg(c),matmul(sk,transpose(c)))
  do i=1,n; norm(i,i)=norm(i,i)-1.; end do
  call check(maxval(abs(norm))<tol,'generalized eigenvector normalization')

  do ie=0,1
    eta=0.08*ie
    call optsp_s(e(1),c(1,:),e(2),c(2,:),k,0.,0.,eta,n,nr,rv,h,hi,s,tau,z(1),z(2),z(3))
    call optspbz_s(e(1),c(1,:),e(2),c(2,:),k,0.,0.,n,nr,rv,h,hi,s,tau,bz(1),bz(2),bz(3))
    call check(maxval(abs(bz/cmplx(gap,eta)-z))<tol,'OPT_BZ and IPA midpoint agreement')
    do axis=1,3
      dh=0.; ds=0.; xk=0.; dk=0.
      do i=1,nr
        ph=exp(cmplx(0.,k*rv(i,1)))
        dh=dh+cmplx(0.,rv(i,axis))*ph*h(i,:,:)
        ds=ds+cmplx(0.,rv(i,axis))*ph*s(i,:,:)
        dk=dk+ph*dm(i,:,:,axis)
        do b=1,n; do a=1,n
          xk(a,b)=xk(a,b)+ph*(dm(i,a,b,axis)+ &
            0.5*(tau(a,axis)+tau(b,axis)+rv(i,axis))*s(i,a,b))
        end do; end do
      end do
      call sandwich(n,c(2,:),dk,c(1,:),dcv)
      call sandwich(n,c(2,:),dh-e(2)*ds+cmplx(0.,gap)*xk,c(1,:),full)
      mid=full-cmplx(0.,gap)*dcv
      call check(abs(z(axis)-mid/cmplx(gap,eta))<tol,'midpoint vs lattice-gauge formula')
      bz(axis)=full/cmplx(gap,eta)
    end do
    call optdip_s(gap,eta,c(1,:),c(2,:),k,0.,0.,n,nr,rv,dm,z(1),z(2),z(3))
    call check(maxval(abs(z-bz))<tol,'full position vertex vs lattice-gauge formula')

    ! A covariant change H -> H+mu*S, E -> E+mu leaves the dipoles invariant.
    mu=0.4
    call optsp_s(e(1)+mu,c(1,:),e(2)+mu,c(2,:),k,0.,0.,eta,n,nr,rv,h+mu*s,hi,s,tau, &
      z(1),z(2),z(3))
    call optdip_s(gap,eta,c(1,:),c(2,:),k,0.,0.,n,nr,rv,dm,z(1),z(2),z(3))
    call check(maxval(abs(z-bz))<tol,'covariant energy-reference invariance')
  end do

  data%num_wann=n; data%nrpts=nr
  allocate(data%rvec(3,nr),data%rmn(3,n,n,nr))
  data%rvec=int(transpose(rv)); data%rmn=0.
  ! rmn_data uses integer translations; use an identity lattice with integer R.
  rv(:,2:3)=0.
  do i=1,nr; do b=1,n; do a=1,n
    data%rmn(:,a,b,i)=cmplx(dm(i,a,b,:)+0.5*(tau(a,:)+tau(b,:)+rv(i,:))*s(i,a,b),0.)
  end do; end do; end do
  call rmn_dfts_dipoles(data,n,nr,rv,lat,s,tau,dout,ok,message)
  call check(ok,'position-matrix midpoint subtraction: '//trim(message))
  call check(maxval(abs(dout-dm))<tol,'residual position matrix recovery')
  call rmn_destroy(data)

  call optcondcalc(0.,0.,2.,0.,cond)
  call check(abs(cond)==0.,'zero Cartesian off-diagonal conductivity')
  call optcondcalc(1.,0.,2.,1.,cond)
  call check(abs(cond)==0.,'zero diagonal/circular vacuum conductivity')
  call optcondcalc(2.3,0.7,2.,0.,cond)
  call check(abs(cond/cmplx(0.,-2./(4.*pi*hbar))-cmplx(2.3,0.7))<tol, &
    'nonzero Cartesian off-diagonal conductivity')
  call optcondcalc(2.3,0.7,2.,1.,cond)
  call check(abs(cond/cmplx(0.,-2./(4.*pi*hbar))-cmplx(1.3,0.7))<tol, &
    'nonzero diagonal/circular conductivity')

  ! Normal incidence on a sheet: 1-R-T = 4 ni Re(y)/|ni+ns+y|^2, y = 4 pi d sigma/c,
  ! which vanishes for a lossless sheet, and T is the same from either side.
  do i=1,size(nis)
    do j=1,size(sigs)
      call transmittance(nis(i),nss(i),thick,sigs(j),tr)
      call reflectance(nis(i),nss(i),thick,sigs(j),rf)
      call transmittance(nss(i),nis(i),thick,sigs(j),tr2)
      call reflectance(nss(i),nis(i),thick,sigs(j),rf2)
      y=4.*pi*thick*sigs(j)/3.e8
      call check(abs(tr-4.*nis(i)*nss(i)/abs(nis(i)+nss(i)+y)**2)<tol,'sheet transmittance')
      call check(abs(rf-abs((nis(i)-nss(i)-y)/(nis(i)+nss(i)+y))**2)<tol,'sheet reflectance')
      call check(abs(1.-rf-tr-4.*nis(i)*real(y)/abs(nis(i)+nss(i)+y)**2)<tol,'sheet absorptance')
      call check(abs(tr-tr2)<tol,'sheet transmittance from either side')
      if (real(sigs(j))==0.) then
        call check(abs(1.-rf-tr)<tol,'lossless sheet energy conservation')
        call check(abs(1.-rf2-tr2)<tol,'lossless sheet energy conservation, reversed')
      end if
    end do
  end do
  print '(A)', 'nonorth_optics: PASS'
contains
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print '(2A)', 'FAIL: ',label
      error stop 1
    end if
  end subroutine check
end program test_nonorth_optics
