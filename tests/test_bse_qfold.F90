program test_bse_qfold

	! BSE_QFOLD = M: images of the direct kernel of an M x M supercell run (bse_q_images_fold, matrizelbse) on a hexagonal lattice
	! with M = 4, 3 x 3 supercell mesh (= 12 x 12 of the unit cell).
	use bse_formfactor, only: ff_data
	use bse_q_optics, only: rmn_data
	use input_variables, only: bse_qfold
	implicit none
	integer,parameter :: m=4,nmax=27*m*m
	real,parameter :: pi=acos(-1.)
	real :: lat(3,3),ulat(3,3),b1(3),b2(3),b3(3),ed(3),zero(3),k1(3),k2(3),ku(3),kpu(3),s(3),q(3)
	real :: shift(3,nmax),wimg(nmax),ushift(3,27),oshift(3,27),w1(nmax)
	integer :: nimg,unimg,onimg,i,j,i1,i2,gc(2),gp(2),g(2),ik,jk,icl,m1,m2,n1,n2,ntie,nties,npair
	real :: classw(0:m*m-1),ql,pref,expected,qd,vcoul0
	logical :: found
	integer :: checks,grid(3),st1(4),st2(4)
	real :: rvec(1,3),centers(3,2)
	complex :: c(2),v(2),sk(2,2),phase(2),got,matrizelbse
	type(ff_data) :: ff
	type(rmn_data) :: rmn

	checks=0;nties=0;npair=0
	zero=0.0;ed=(/1.0,3.5,1.0/);grid=(/3,3,1/)
	! supercell 4 x 4 of a hexagonal cell, a = 3 A, 30 A of vacuum
	lat=0.0
	lat(1,:)=(/4.0*3.0,0.0,0.0/)
	lat(2,:)=(/-0.5*4.0*3.0,0.5*sqrt(3.0)*4.0*3.0,0.0/)
	lat(3,3)=30.0
	ulat=lat;ulat(1:2,:)=lat(1:2,:)/real(m)
	call recvec(lat(1,:),lat(2,:),lat(3,:),b1,b2,b3)

	do ik=0,8
		do jk=0,8
			k1=real(mod(ik,3))/3.0*b1+real(ik/3)/3.0*b2
			k2=real(mod(jk,3))/3.0*b1+real(jk/3)/3.0*b2
			call bse_q_images(k1,k2,lat,onimg,oshift)

			! M = 1 is bse_q_images with weights 1/n
			call bse_q_images_fold(k1,k2,lat,1,nmax,nimg,shift,wimg)
			call ok(nimg .eq. onimg,'M=1 image count')
			do i=1,nimg
				call ok(all(shift(:,i) .eq. oshift(:,i)),'M=1 shifts')
				call ok(abs(wimg(i)-1.0/real(nimg)) .lt. 1.0e-7,'M=1 weights')
			end do

			call bse_q_images_fold(k1,k2,lat,m,nmax,nimg,shift,wimg)
			call ok(abs(sum(wimg(1:nimg))-real(m*m)) .lt. 1.0e-4,'weights sum to M^2')
			classw=0.0
			do i=1,nimg
				! crystal coordinates of the shift in the basis B_1, B_2, B_3: integers, none along B_3
				n1=nint(dot_product(shift(:,i),lat(1,:))/(2.0*pi))
				n2=nint(dot_product(shift(:,i),lat(2,:))/(2.0*pi))
				call ok(abs(dot_product(shift(:,i),lat(1,:))/(2.0*pi)-real(n1)) .lt. 1.0e-3,'shift in Gamma_SC (1)')
				call ok(abs(dot_product(shift(:,i),lat(2,:))/(2.0*pi)-real(n2)) .lt. 1.0e-3,'shift in Gamma_SC (2)')
				call ok(abs(dot_product(shift(:,i),lat(3,:))) .lt. 1.0e-3,'no shift along a_3')
				icl=modulo(n1,m)+m*modulo(n2,m)
				classw(icl)=classw(icl)+wimg(i)
				q=k1-k2-shift(:,i)
				ql=sqrt(dot_product(q,q))
				! shortest of its class over Gamma_UC = M Gamma_SC (brute force) ...
				do m1=-6,6
					do m2=-6,6
						s=real(m*m1)*b1+real(m*m2)*b2
						qd=sqrt(dot_product(q-s,q-s))
						call ok(ql .le. qd*(1.0+2.0e-4)+2.0e-6,'image shorter than the rest of its class')
					end do
				end do
			end do
			! ... and every class once, with weight 1
			do icl=0,m*m-1
				call ok(abs(classw(icl)-1.0) .lt. 1.0e-4,'weight 1 per class')
			end do
			if (nimg .gt. m*m) nties=nties+1
			npair=npair+1

			! the unit-cell kernel of the folded pair (K+g, K'+g'): its images are the images of the class g-g'
			do i1=0,m*m-1
				do i2=0,m*m-1
					gc=(/mod(i1,m),i1/m/);gp=(/mod(i2,m),i2/m/)
					ku=k1+real(gc(1))*b1+real(gc(2))*b2
					kpu=k2+real(gp(1))*b1+real(gp(2))*b2
					call bse_q_images(ku,kpu,ulat,unimg,ushift)
					do j=1,unimg
						! kpu+ushift(j) = k2 + shift_j with shift_j = ushift - gc + gp
						s=ushift(:,j)-(real(gc(1))*b1+real(gc(2))*b2)+(real(gp(1))*b1+real(gp(2))*b2)
						found=.false.
						do i=1,nimg
							if (sqrt(dot_product(shift(:,i)-s,shift(:,i)-s)) .lt. 1.0e-4) found=.true.
						end do
						call ok(found,'unit-cell image of a fold pair is an image of the supercell run')
					end do
					! no more of them than the class has (ties in the unit cell = ties of the class)
					g=(/modulo(gc(1)-gp(1),m),modulo(gc(2)-gp(2),m)/)
					ntie=0
					do i=1,nimg
						n1=nint(dot_product(shift(:,i),lat(1,:))/(2.0*pi))
						n2=nint(dot_product(shift(:,i),lat(2,:))/(2.0*pi))
						if (modulo(n1,m) .eq. modulo(-g(1),m) .and. modulo(n2,m) .eq. modulo(-g(2),m)) ntie=ntie+1
					end do
					call ok(ntie .eq. unimg,'ties of the class = ties of the unit-cell pair')
				end do
			end do
		end do
	end do
	write(*,'(A,I0,A,I0,A)') 'bse_q_images_fold: ',npair,' k pairs, ',nties,' of them with ties'

	! Kernel dispatch: the elements are the sums of W over the classes, weight 1 per class, the diagonal element keeps the
	! transition energy once. All vertices are 1 (one state per band, no orbital centres), W = V2DLAVG (1/|q| away from 0).
	c=(/cmplx(1.,0.),cmplx(0.,0.)/);v=(/cmplx(0.,0.),cmplx(1.,0.)/)
	sk=0.;sk(1,1)=1.;sk(2,2)=1.;phase=1.;centers=0.;rvec=0.
	vcoul0=v2dlavg_(zero,zero)
	q=0.1*b1+0.05*b2
	pref=v2dlavg_(q,zero)*sqrt(dot_product(q,q))
	do ik=0,8
		do jk=0,8
			k1=real(mod(ik,3))/3.0*b1+real(ik/3)/3.0*b2
			k2=real(mod(jk,3))/3.0*b1+real(jk/3)/3.0*b2
			st1=(/ik+1,1,2,ik+1/);st2=(/jk+1,1,2,jk+1/)
			bse_qfold=m
			got=matrizelbse('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,lat,st1,2.,-1.,c,v,k1, &
			    st2,2.,-1.,c,v,k2,'W',1,rvec,sk,sk,.false.,phase,phase,.false.,rmn,centers,ff)
			expected=0.0
			do icl=0,m*m-1
				gc=(/mod(icl,m),icl/m/)
				call class_w(k1-k2,gc,expected,pref)
			end do
			! the q = 0 image of the diagonal element (class g = 0) is ec-ev + W(0), the other classes go through the vertices
			if (ik .eq. jk) expected=expected+3.0
			call ok(abs(real(got)-expected) .le. 2.0e-5*abs(expected)+1.0e-6,'kernel element = sum over the classes')
			call ok(abs(aimag(got)) .lt. 1.0e-6,'kernel element is real')
			! the same pair again with another dielectric constant: W scales, the cached images do not leak
			ed(2)=2.0*ed(2)
			got=matrizelbse('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,lat,st1,2.,-1.,c,v,k1, &
			    st2,2.,-1.,c,v,k2,'W',1,rvec,sk,sk,.false.,phase,phase,.false.,rmn,centers,ff)
			ed(2)=0.5*ed(2)
			if (ik .ne. jk) call ok(abs(real(got)-0.5*expected) .le. 2.0e-5*abs(expected)+1.0e-6,'cache keyed on EDIEL')
			! M = 1: only the shortest image(s)
			bse_qfold=1
			got=matrizelbse('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,lat,st1,2.,-1.,c,v,k1, &
			    st2,2.,-1.,c,v,k2,'W',1,rvec,sk,sk,.false.,phase,phase,.false.,rmn,centers,ff)
			call bse_q_images(k1,k2,lat,onimg,oshift)
			expected=0.0
			do i=1,onimg
				q=k1-k2-oshift(:,i)
				if (ik .eq. jk) then
					expected=expected+3.0+vcoul0
				else
					expected=expected+pref/sqrt(dot_product(q,q))
				end if
			end do
			expected=expected/real(onimg)
			call ok(abs(real(got)-expected) .le. 2.0e-5*abs(expected)+1.0e-6,'M=1 kernel element')
		end do
	end do
	write(*,'(A,I0,A)') 'bse_qfold: PASS (',checks,' checks)'

contains

	real function v2dlavg_(qa,qb)
		real,intent(in) :: qa(3),qb(3)
		real :: v2dlavg
		v2dlavg_=v2dlavg(qa,qb,ed,grid,lat,1.e-6)
	end function v2dlavg_

	subroutine class_w(p,g0,acc,pre)
		! sum += W of the shortest image of the class (p - g0 + Gamma_UC), brute force, ties averaged (1e-4 as in the code)
		real,intent(in) :: p(3),pre
		integer,intent(in) :: g0(2)
		real,intent(inout) :: acc
		real :: ql(-4:4,-4:4),qm,dd(3),gg(3)
		integer :: a1,a2,nt
		gg=real(g0(1))*b1+real(g0(2))*b2
		do a1=-4,4
			do a2=-4,4
				dd=p-gg-real(m*a1)*b1-real(m*a2)*b2
				ql(a1,a2)=sqrt(dot_product(dd,dd))
			end do
		end do
		qm=minval(ql)
		nt=count(ql .le. qm*(1.0+1.0e-4)+1.0e-6)
		if (qm .lt. 1.0e-6) then
			acc=acc+vcoul0
			return
		end if
		do a1=-4,4
			do a2=-4,4
				if (ql(a1,a2) .le. qm*(1.0+1.0e-4)+1.0e-6) acc=acc+pre/ql(a1,a2)/real(nt)
			end do
		end do
	end subroutine class_w

	subroutine ok(cond,label)
		logical,intent(in) :: cond
		character(len=*),intent(in) :: label
		checks=checks+1
		if (.not. cond) then
			write(*,*) 'FAIL: ',label
			stop 1
		end if
	end subroutine ok

end program test_bse_qfold
