subroutine bsesolvertemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef,st,phavg,ta,temp,nocpf,fermishift,bsealgo,bsehamwrite,bsehamread,bsehamfile,dft,mag)
    implicit none
    integer :: nthreads,ngrid(3),nc,nv,excwf0,excwff,nocpf
    real :: numdos,ebse0,ebsef,numbse,sme,ktol,ediel(3),exc,mshift(3)
    real :: ez,w1,r0,lc,rk,fermishift,mag(3)
    character(len=70) :: outputfolder,calcparms,params,kpaths,kpathsbse,orbw,meshtype,bsehamfile
    character(len=12) :: bsealgo
    character(len=7) :: coultype
    character(len=1) :: dft
    logical :: bsewf,dtfull,cpol,tmcoef,bsehamwrite,bsehamread
    real :: st,phavg,temp
    character(len=2) :: ta

    ! Both temperature branches share the distributed assembly, solvers and optics.
    call bsesolver_core(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef,nocpf,fermishift,bsealgo,bsehamwrite,bsehamread,bsehamfile,dft,mag, &
        st,phavg,ta,temp)
end subroutine bsesolvertemp
