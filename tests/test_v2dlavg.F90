program test_v2dlavg

	use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
	use bse_formfactor, only: ff_data
	use bse_q_optics, only: rmn_data
	implicit none
	real :: v2dlavg,v2diel,fermidisteh
	real(kind=8) :: v2davg_cell
	complex :: matrizelbse,matrizelbsetemp,matrizelbsekq,matrizelbsekqtemp
	real :: lat(3,3),other(3,3),ed(3),bare(3),zero(3),k(3),kp(3),q(4)
	real :: values(128),v0,vq,vd,expected,temperature
	real :: centers(3,2),rvec(1,3)
	complex :: c(2),v(2),sk(2,2),phase(2),got
	type(ff_data) :: ff
	type(rmn_data) :: rmn
	integer :: grid(3),fine(3),i,st1(4),st2(4),checks
	real(kind=8) :: pi,cic,hx,hy,rectangle,reference,area

	pi=acos(-1.0d0);cic=-0.0904756d3
	zero=0.0;bare=1.0;ed=3.5;checks=0
	lat=0.0;lat(1,1)=4.0;lat(2,2)=5.0;lat(3,3)=20.0
	grid=(/12,17,1/)
	hx=pi/(4.0d0*grid(1));hy=pi/(5.0d0*grid(2))
	! Independent rectangular integral: 4 [hx asinh(hy/hx) + hy asinh(hx/hy)].
	rectangle=(hx*asinh(hy/hx)+hy*asinh(hx/hy))/(hx*hy)
	expected=real(cic*rectangle/(20.0d0*product(grid)*ed(2)))
	! First use from multiple threads also checks publication of the cache.
	!$omp parallel do num_threads(4)
	do i=1,size(values)
		values(i)=v2dlavg(zero,zero,ed,grid,lat,1.e-6)
	end do
	!$omp end parallel do
	do i=1,size(values)
		call close(real(values(i),8),real(expected,8),'threaded rectangular q=0',2.d-6)
	end do
	call close(v2davg_cell(grid,lat),rectangle,'analytic rectangular average',1.d-12)
	v0=values(1)
	k=(/.12,.09,0.0/)
	vq=v2dlavg(k,zero,ed,grid,lat,1.e-6)
	call close(real(vq,8),real(v2diel(k,zero,ed,lat,grid,1.e-6),8),'nonzero V2DL equality',1.d-7)
	call close(real(vq,8),cic/(20.0d0*product(grid)*ed(2)*.15d0),'nonzero analytic Coulomb',2.d-6)
	call close(real(v2diel(zero,zero,ed,lat,grid,1.e-6),8),0.0d0,'legacy q=0 unchanged',1.d-7)
	ed(2)=1.0
	call close(real(v2dlavg(zero,zero,ed,grid,lat,1.e-6),8),real(3.5*v0,8),'q=0 dielectric scaling',2.d-6)
	call close(real(v2dlavg(k,zero,ed,grid,lat,1.e-6),8),real(3.5*vq,8),'q!=0 dielectric scaling',2.d-6)
	ed=(/8.0,3.5,9.0/)
	call close(real(v2dlavg(zero,zero,ed,grid,lat,1.e-6),8),real(v0,8),'only central EDIEL used',1.d-7)
	other=lat;other(3,3)=80.0
	call close(real(v2dlavg(zero,zero,ed,grid,other,1.e-6),8),real(v0,8),'vacuum independence q=0',1.d-7)
	call close(real(v2dlavg(k,zero,ed,grid,other,1.e-6),8),real(vq,8),'vacuum independence q!=0',1.d-7)
	fine=grid;fine(1:2)=2*grid(1:2)
	call close(real(v2dlavg(zero,zero,ed,fine,lat,1.e-6),8),real(v0/2,8),'q=0 mesh scaling',2.d-6)
	fine=grid;fine(3)=2
	call close(real(v2dlavg(zero,zero,ed,fine,lat,1.e-6),8),real(v0/2,8),'Nk normalization',2.d-6)
	other=lat;other(1:2,:)=2*lat(1:2,:)
	call close(real(v2dlavg(zero,zero,ed,grid,other,1.e-6),8),real(v0/2,8),'lattice scale and cache',2.d-6)
	call close(real(v2dlavg(zero,zero,ed,grid,lat,1.e-6),8),real(v0,8),'cache returns to original cell',1.d-7)

	! Skew and anisotropic grids against independent numerical triangle integrals.
	do i=1,4
		other=lat;other(2,1)=real(i)-2.4
		fine=(/7*i,31-i,1/)
		reference=numerical_average(other,fine)
		call close(v2davg_cell(fine,other),reference,'skew numerical reference',2.d-8)
	end do
	! Rigid rotation puts the same rectangular sheet into the yz plane.
	other=0.0;other(1,2)=4.0;other(2,3)=5.0;other(3,1)=20.0
	call close(v2davg_cell(grid,other),rectangle,'rotation invariance',1.d-12)
	! Very unequal rectangular steps should not require more quadrature points.
	fine=(/1,10000,1/)
	hx=pi/(4.0d0*fine(1));hy=pi/(5.0d0*fine(2))
	rectangle=(hx*asinh(hy/hx)+hy*asinh(hx/hy))/(hx*hy)
	call close(v2davg_cell(fine,lat),rectangle,'anisotropic rectangle',1.d-11)

	! Linked Q=0 and finite-temperature dispatch on a normalized diagonal state.
	c=(/cmplx(1.,0.),cmplx(0.,0.)/);v=(/cmplx(0.,0.),cmplx(1.,0.)/)
	sk=0.;sk(1,1)=1.;sk(2,2)=1.;phase=1.;centers=0.;rvec=0.
	st1=(/1,1,2,1/);st2=st1
	expected=3.0+v0
	got=matrizelbse('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,lat,st1,2.,-1.,c,v,zero, &
	    st2,2.,-1.,c,v,zero,'W',1,rvec,sk,sk,.false.,phase,phase,.false.,rmn,centers,ff)
	call close(real(got,8),real(expected,8),'Q=0 kernel dispatch',2.d-6)
	temperature=3000.
	expected=3.0+v0*fermidisteh(2.,-1.,temperature)
	got=matrizelbsetemp('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,lat,st1,2.,-1.,c,v,zero, &
	    st2,2.,-1.,c,v,zero,temperature,'W',1,rvec,sk,sk,.false.,phase,phase,.false.,rmn,centers,ff)
	call close(real(got,8),real(expected,8),'finite-T Q=0 dispatch',2.d-6)

	! Finite Q: orthogonal mixed states with different orbital centres yield
	! nonzero exchange. Its coefficient must remain bare as EDIEL changes.
	c=cmplx(1./sqrt(2.),0.);v=c;v(2)=-v(2)
	centers(1,2)=.31
	k=(/.11,0.,0./);kp=0.;q=(/0.,.2,0.,0./);st2(1)=2
	do i=1,2
		ed(2)=real(i)*2.0
		vd=v2dlavg(k,kp,ed,grid,lat,1.e-6)
		vq=v2dlavg(q(2:4),zero,bare,grid,lat,1.e-6)
		expected=vd*cos(.11*.31/2)**2-vq*sin(.2*.31/2)**2
		got=matrizelbsekq('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,q,lat,st1,2.,-1.,c,v,k, &
		    st2,2.,-1.,c,v,kp,'W',1,rvec,sk,sk,sk,sk,.true.,centers)
		call close(real(got,8),real(expected,8),'finite-Q screened direct / bare exchange',2.d-6)
		expected=expected*fermidisteh(2.,-1.,temperature)
		got=matrizelbsekqtemp('V2DLAVG   ',1.e-6,2,ed,0.,1.,0.,0.,grid,q,lat,st1,2.,-1.,c,v,k, &
		    st2,2.,-1.,c,v,kp,temperature,'W',1,rvec,sk,sk,sk,sk,.true.,centers)
		call close(real(got,8),real(expected,8),'finite-T finite-Q dispatch',2.d-6)
	end do
	write(*,'(A,I0,A)') 'v2dlavg: PASS (',checks,' checks)'

contains

	subroutine close(actual,want,label,tolerance)
		real(kind=8),intent(in) :: actual,want,tolerance
		character(len=*),intent(in) :: label
		checks=checks+1
		if (.not. ieee_is_finite(actual) .or. abs(actual-want)>tolerance*max(abs(want),1.e-12_8)) then
			write(*,*) 'FAIL: ',label,actual,want
			stop 1
		end if
	end subroutine close

	function numerical_average(lattice,mesh) result(avg)
		! Numerical Duffy triangles: the radial integral of s/(s|d|)
		! is exactly 1/|d|; integrate the edge with 20000 midpoint samples.
		real,intent(in) :: lattice(3,3)
		integer,intent(in) :: mesh(3)
		real(kind=8) :: avg,b1(2),b2(2),verts(2,4),aa(2),bb(2),d(2),s,total,jac
		integer :: edge,jj,point
		area=abs(dble(lattice(1,1))*lattice(2,2)-dble(lattice(1,2))*lattice(2,1))
		b1=2*pi*(/dble(lattice(2,2)),-dble(lattice(2,1))/)/area/mesh(1)
		b2=2*pi*(/-dble(lattice(1,2)),dble(lattice(1,1))/)/area/mesh(2)
		verts(:,1)=(b1+b2)/2;verts(:,2)=(-b1+b2)/2
		verts(:,3)=-verts(:,1);verts(:,4)=-verts(:,2);total=0.
		do edge=1,4
			jj=mod(edge,4)+1;aa=verts(:,edge);bb=verts(:,jj)
			jac=abs(aa(1)*bb(2)-aa(2)*bb(1))
			do point=1,20000
				s=(point-.5d0)/20000;d=(1-s)*aa+s*bb
				total=total+jac/sqrt(sum(d*d))/20000
			end do
		end do
		avg=total/abs(b1(1)*b2(2)-b1(2)*b2(1))
	end function numerical_average

end program test_v2dlavg
