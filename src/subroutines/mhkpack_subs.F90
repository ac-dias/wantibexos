
subroutine rkmesh(rk,rlat,ngrid)

	implicit none
	
	integer,dimension(3) :: ngrid
	real :: rk
	real,dimension(3,3) :: rlat,blat
	real,dimension(3) :: vsize
	real,parameter :: pi=acos(-1.0)

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))

	call vecsize(blat(1,:),vsize(1))
	call vecsize(blat(2,:),vsize(2))
	call vecsize(blat(3,:),vsize(3))

	vsize = vsize/(2*pi)

	ngrid(1) = int(max(1.0,(rk*vsize(1))+0.5))
	ngrid(2) = int(max(1.0,(rk*vsize(2))+0.5))
	ngrid(3) = int(max(1.0,(rk*vsize(3))+0.5))



end subroutine rkmesh

subroutine rkmesh2D(rk,rlat,ngrid)

	implicit none
	
	integer,dimension(3) :: ngrid
	real :: rk
	real,dimension(3,3) :: rlat,blat
	real,dimension(3) :: vsize
	real,parameter :: pi=acos(-1.0)

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))

	call vecsize(blat(1,:),vsize(1))
	call vecsize(blat(2,:),vsize(2))
	call vecsize(blat(3,:),vsize(3))

	vsize = vsize/(2*pi)

	ngrid(1) = int(max(1.0,(rk*vsize(1))+0.5))
	ngrid(2) = int(max(1.0,(rk*vsize(2))+0.5))
	ngrid(3) = 1



end subroutine rkmesh2D


subroutine rkmesh1D(rk,rlat,ngrid)

	implicit none
	
	integer,dimension(3) :: ngrid
	real :: rk
	real,dimension(3,3) :: rlat,blat
	real,dimension(3) :: vsize
	real,parameter :: pi=acos(-1.0)

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))

	call vecsize(blat(1,:),vsize(1))
	call vecsize(blat(2,:),vsize(2))
	call vecsize(blat(3,:),vsize(3))

	vsize = vsize/(2*pi)

	ngrid(1) = int(max(1.0,(rk*vsize(1))+0.5))
	ngrid(2) = 1
	ngrid(3) = 1



end subroutine rkmesh1D

subroutine kpath(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3

	real,dimension((nks/2)*npts,4) :: kpt

	real,dimension(nks,3) :: ks
	real,dimension(nks,3) :: ksaux

	real,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	real :: kdistot !soma do caminho todo
	real :: kdisaux


	
	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS output file"


	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do



	kdistot= 0.0
	do i=1,nks/2

		call distvec(ksaux(2*i,:),ksaux(2*i-1,:),kdis(i))

		!kdis(i)=sqrt((ksaux(2*i,1)-ksaux(2*i-1,1))**2+(ksaux(2*i,2)-ksaux(2*i-1,2))**2+(ksaux(2*i,3)-ksaux(2*i-1,3))**2)
		kdistot=kdistot+kdis(i)

	end do

	kdis=kdis/kdistot


	kdisaux=0.0
	write(1733,*) real(kdisaux)

	 do j=1,nks/2
		
	  do i=1,npts

		!kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*dble((i)/(npts))
		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,2) = ksaux(2*j-1,1)+ (ksaux(2*j,1)-ksaux(2*j-1,1))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,3) = ksaux(2*j-1,2)+ (ksaux(2*j,2)-ksaux(2*j-1,2))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,4) = ksaux(2*j-1,3)+ (ksaux(2*j,3)-ksaux(2*j-1,3))*(i-1.)/(npts-1.)	
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)
		

	  end do

		kdisaux=kdisaux+kdis(j)
		write(1733,*) real(kdisaux)
	 end do


	close(1733)


end subroutine kpath

subroutine qgridbse(rlat1,rlat2,rlat3,nqgrid,shift,qpt)

	implicit none

	integer,dimension(3) :: nqgrid
	integer :: i,j,k,counter
	real,dimension(3) :: rlat1,rlat2,rlat3,shift
	real,dimension(3) :: blat1,blat2,blat3
	real,dimension(nqgrid(1)*nqgrid(2)*nqgrid(3),4) :: qpt

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	counter=1
	do i=0,nqgrid(1)-1
		do j=0,nqgrid(2)-1
			do k=0,nqgrid(3)-1
				! qpt(:,2:4) is Cartesian.  shift is in units of a grid spacing.
				qpt(counter,1)=real(counter)
				qpt(counter,2)=((real(i)+shift(1))/real(nqgrid(1)))*blat1(1) &
					+ ((real(j)+shift(2))/real(nqgrid(2)))*blat2(1) &
					+ ((real(k)+shift(3))/real(nqgrid(3)))*blat3(1)
				qpt(counter,3)=((real(i)+shift(1))/real(nqgrid(1)))*blat1(2) &
					+ ((real(j)+shift(2))/real(nqgrid(2)))*blat2(2) &
					+ ((real(k)+shift(3))/real(nqgrid(3)))*blat3(2)
				qpt(counter,4)=((real(i)+shift(1))/real(nqgrid(1)))*blat1(3) &
					+ ((real(j)+shift(2))/real(nqgrid(2)))*blat2(3) &
					+ ((real(k)+shift(3))/real(nqgrid(3)))*blat3(3)
				counter=counter+1
			end do
		end do
	end do

end subroutine qgridbse

subroutine kpathbse(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3

	real,dimension((nks/2)*npts,4) :: kpt

	real,dimension(nks,3) :: ks
	real,dimension(nks,3) :: ksaux

	real,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	real :: kdistot !soma do caminho todo
	real :: kdisaux


	
	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS-BSE.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS-BSE output file"


	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do



	kdistot= 0.0
	do i=1,nks/2

		call distvec(ksaux(2*i,:),ksaux(2*i-1,:),kdis(i))

		!kdis(i)=sqrt((ksaux(2*i,1)-ksaux(2*i-1,1))**2+(ksaux(2*i,2)-ksaux(2*i-1,2))**2+(ksaux(2*i,3)-ksaux(2*i-1,3))**2)
		kdistot=kdistot+kdis(i)

	end do

	kdis=kdis/kdistot


	kdisaux=0.0
	write(1733,*) real(kdisaux)

	 do j=1,nks/2
		
	  do i=1,npts

		!kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*dble((i)/(npts))
		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,2) = ksaux(2*j-1,1)+ (ksaux(2*j,1)-ksaux(2*j-1,1))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,3) = ksaux(2*j-1,2)+ (ksaux(2*j,2)-ksaux(2*j-1,2))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,4) = ksaux(2*j-1,3)+ (ksaux(2*j,3)-ksaux(2*j-1,3))*(i-1.)/(npts-1.)	
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)
		

	  end do

		kdisaux=kdisaux+kdis(j)
		write(1733,*) real(kdisaux)
	 end do


	close(1733)


end subroutine kpathbse



subroutine distvec(vec1,vec2,distvecout) !subroutina para calcular a distância entre dois vetores

	implicit none

	real,dimension(3) :: vec1,vec2
	real :: distvecout

	distvecout=sqrt((vec1(1)-vec2(1))**2+(vec1(2)-vec2(2))**2+(vec1(3)-vec2(3))**2)


end subroutine distvec


subroutine monhkhorst_packq(q,n1,n2,n3,shift,rlat1,rlat2,rlat3,qpt)

	implicit none
	
	integer :: n1,n2,n3
	real,dimension(3) :: shift,kshift
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(n1*n2*n3,3) :: qpt

	real,dimension(4) :: q

	integer :: i,j,k, counter
	real,dimension(3) :: blat1,blat2,blat3

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	kshift(1) = blat1(1)*shift(1)+blat2(1)*shift(2)+blat3(1)*shift(3) 
	kshift(2) = blat1(2)*shift(1)+blat2(2)*shift(2)+blat3(2)*shift(3)
	kshift(3) = blat1(3)*shift(1)+blat2(3)*shift(2)+blat3(3)*shift(3)

	counter = 1

	do i=0,n1-1
	 do j=0,n2-1
	  do k=0,n3-1

	    !qpt(counter,1) = q(2)+(blat1(1)/dble(n1))*(dble(i)+shift(1))+(blat2(1)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(1)/dble(n3))*(dble(k)+shift(3))

	    !qpt(counter,2) = q(3)+(blat1(2)/dble(n1))*(dble(i)+shift(1))+(blat2(2)/dble(n2))*(dble(j)+shift(2))&
	!		     +(blat3(2)/dble(n3))*(dble(k)+shift(3))

	    !qpt(counter,3) =q(4)+(blat1(3)/dble(n1))*(dble(i)+shift(1))+(blat2(3)/dble(n2))*(dble(j)+shift(2))&
	!		     +(blat3(3)/dble(n3))*(dble(k)+shift(3))	


	    qpt(counter,1) = q(2)+(blat1(1)/dble(n1))*(dble(i))+(blat2(1)/dble(n2))*(dble(j))&
			     +(blat3(1)/dble(n3))*(dble(k))+kshift(1)

	    qpt(counter,2) = q(3)+(blat1(2)/dble(n1))*(dble(i))+(blat2(2)/dble(n2))*(dble(j))&
			     +(blat3(2)/dble(n3))*(dble(k))+kshift(2)

	    qpt(counter,3) =q(4)+(blat1(3)/dble(n1))*(dble(i))+(blat2(3)/dble(n2))*(dble(j))&
			     +(blat3(3)/dble(n3))*(dble(k))+kshift(3)	
   
	    counter = counter+1

	  end do
	 end do
	end do

end subroutine monhkhorst_packq


subroutine monhkhorst_pack(n1,n2,n3,shift,rlat1,rlat2,rlat3,kpt)

	implicit none
	
	integer :: n1,n2,n3
	real,dimension(3) :: shift,kshift
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(n1*n2*n3,3) :: kpt

	integer :: i,j,k, counter
	real,dimension(3) :: blat1,blat2,blat3

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	counter = 1

	kshift(1) = blat1(1)*shift(1)+blat2(1)*shift(2)+blat3(1)*shift(3) 
	kshift(2) = blat1(2)*shift(1)+blat2(2)*shift(2)+blat3(2)*shift(3)
	kshift(3) = blat1(3)*shift(1)+blat2(3)*shift(2)+blat3(3)*shift(3)

	do i=0,n1-1
	 do j=0,n2-1
	  do k=0,n3-1

	    !kpt(counter,1) = (blat1(1)/dble(n1))*(dble(i)+shift(1))+(blat2(1)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(1)/dble(n3))*(dble(k)+shift(3))

	    !kpt(counter,2) = (blat1(2)/dble(n1))*(dble(i)+shift(1))+(blat2(2)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(2)/dble(n3))*(dble(k)+shift(3))

	    !kpt(counter,3) =(blat1(3)/dble(n1))*(dble(i)+shift(1))+(blat2(3)/dble(n2))*(dble(j)+shift(2))&
		!	     +(blat3(3)/dble(n3))*(dble(k)+shift(3))	

	    kpt(counter,1) = (blat1(1)/real(n1))*(real(i))+(blat2(1)/real(n2))*(real(j))&
			     +(blat3(1)/real(n3))*(real(k))+kshift(1)

	    kpt(counter,2) = (blat1(2)/real(n1))*(real(i))+(blat2(2)/real(n2))*(real(j))&
			     +(blat3(2)/real(n3))*(real(k))+kshift(2)

	    kpt(counter,3) =(blat1(3)/real(n1))*(real(i))+(blat2(3)/real(n2))*(real(j))&
			     +(blat3(3)/real(n3))*(real(k))+kshift(3)	
   
	    counter = counter+1

	  end do
	 end do
	end do

end subroutine monhkhorst_pack

subroutine recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3) !calcula os vetores da rede recíproca, a partir dos vetores da rede real

	implicit none
	
	real,parameter :: pi=acos(-1.)
	
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3

	real,dimension(3) :: v23,v31,v12
	real :: vol

	call prodvec(rlat2,rlat3,v23)
	call prodvec(rlat3,rlat1,v31)
	call prodvec(rlat1,rlat2,v12)

	vol= abs((rlat1(1)*v23(1))+(rlat1(2)*v23(2))+(rlat1(3)*v23(3)))

	blat1 = ((2.0*pi)/vol)*v23
	blat2 = ((2.0*pi)/vol)*v31
	blat3 = ((2.0*pi)/vol)*v12


end subroutine recvec


subroutine conv_kpt_d2c(rlat1,rlat2,rlat3,kpt)!convert kpoint coordinates from direct to cartesian

	implicit none
	
	real,parameter :: pi=acos(-1.)
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3
	real,dimension(3) :: kpt,kptaux
	
	
	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)
	
	kptaux(1) = kpt(1)*blat1(1)+kpt(2)*blat2(1)+kpt(3)*blat3(1)
	kptaux(2) = kpt(1)*blat1(2)+kpt(2)*blat2(2)+kpt(3)*blat3(2)
	kptaux(3) = kpt(1)*blat1(3)+kpt(2)*blat2(3)+kpt(3)*blat3(3)
	
	kpt = kptaux	 	

end subroutine conv_kpt_d2c

!subrotinas para gerar grid hexagonal

subroutine gridgen2dhexcount(ngridx,ngridy,a,counter)


	implicit none

	real,parameter :: pi=acos(-1.)

	integer :: ngridx,ngridy,i,j

	integer :: counter

	real :: a!constante da rede

	real :: a0

        real :: xmax, ymax

	real :: kx,ky


	a0 = a/sqrt(3.)

        xmax=2*pi/(3.d0*a0)
        ymax=4*pi/(3.d0*a)


	counter = 0


	do i=1,ngridx
                kx=xmax*(i-1)/(ngridx-1)
		do j=1,ngridy
                         ky=ymax*(j-1)/(ngridy-1)
                         if (ky.lt.ymax/2.d0) then
                            counter=counter+1
                            if (i.ne.1.and.j.ne.1) counter=counter+1 
                            if (i.ne.1) counter=counter+1
	                    if (j.ne.1) counter=counter+1
                          
                         else if (ky.lt.(ymax-kx*tan(pi/6.d0))) then
                            counter=counter+1
                            if (i.ne.1.and.j.ne.1) counter=counter+1 
                            if (i.ne.1) counter=counter+1
	                    if (j.ne.1) counter=counter+1
                         end if
		
		end do
	end do 


end subroutine gridgen2dhexcount


subroutine gridgen2dhex(ngridx,ngridy,a,ncounter,kg)


	implicit none

	real,parameter :: pi=acos(-1.)

	integer :: ngridx,ngridy,i,j

	integer :: counter,ncounter

	real :: a!constante da rede

	real :: a0

        real :: xmax, ymax

	real :: kx,ky

	real,dimension(ncounter,3) :: kg


	a0 = a/sqrt(3.)

        xmax=2*pi/(3.d0*a0)
        ymax=4*pi/(3.d0*a)

	counter = 0

	do i=1,ngridx
                kx=xmax*(i-1)/(ngridx-1)
		do j=1,ngridy
                         ky=ymax*(j-1)/(ngridy-1)
                         if (ky.lt.ymax/2.d0) then
                            counter=counter+1 
                            kg(counter, 2)=kx
                            kg(counter, 1)=ky
                            kg(counter, 3)=0.0
                            if (i.ne.1.and.j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=-ky
                                kg(counter, 3)=0.0
                            end if
                            if (i.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=ky
                                kg(counter, 3)=0.0
			    end if
	                    if (j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=kx
                                kg(counter, 1)=-ky
                                kg(counter, 3)=0.0
                            end if
                         else if (ky.lt.(ymax-kx*tan(pi/6.d0))) then
                            counter=counter+1 
                            kg(counter, 2)=kx
                            kg(counter, 1)=ky
                            kg(counter, 3)=0.0
                            if (i.ne.1.and.j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=-ky
                                kg(counter, 3)=0.0
                            end if
                            if (i.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=ky
                                kg(counter, 3)=0.0
			    end if
	                    if (j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=kx
                                kg(counter, 1)=-ky
                                kg(counter, 3)=0.0
                            end if
                         end if
                
		
		end do
        end do

	


end subroutine gridgen2dhex

!subrotina pra gerar um grid quadrado ou retangular
subroutine gridgen2dret(ngridx,ngridy,a,b,kpt)

	implicit none
	real,parameter :: pi=acos(-1.)
	
	integer :: i,j,ngridx,ngridy,counter
	real :: a,b
	real,dimension(ngridx*ngridy,3) :: kpt
	real :: x,y,kx,ky
	
	
	x= pi/a
	y= pi/b
	
	counter = 0
	
	do i=1,ngridx
	
	 kx = -x+ 2.d0*x*(i-1.d0)/(ngridx-1.d0)
	
		do j= 1,ngridy
		
		 	 ky = -y+ 2.d0*y*(j-1.d0)/(ngridy-1.d0)
		 	 counter=counter+1
		 	 
		 	 kpt(counter,1) = kx	
		 	 kpt(counter,2) = ky
		 	 kpt(counter,3) = 0.0	
		
		end do
	
	end do


end subroutine gridgen2dret

subroutine monhkhorst_pack_adp(n1,n2,n3,shift,rlat1,rlat2,rlat3,kpt)

	implicit none
	
	integer :: n1,n2,n3
	real,dimension(3) :: shift,kshift
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(n1*n2*n3,3) :: kpt

	integer :: i,j,k, counter
	real,dimension(3) :: blat1,blat2,blat3
	
	real :: n1aux,n2aux,n3aux

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	counter = 1

	kshift(1) = blat1(1)*shift(1)+blat2(1)*shift(2)+blat3(1)*shift(3) 
	kshift(2) = blat1(2)*shift(1)+blat2(2)*shift(2)+blat3(2)*shift(3)
	kshift(3) = blat1(3)*shift(1)+blat2(3)*shift(2)+blat3(3)*shift(3)

	do i=0,n1-1
	
		n1aux = -0.5 + (1.0)/(real(n1))*(real(i))
		
	 do j=0,n2-1
	  
	  	n2aux = -0.5 + (1.0)/(real(n2))*(real(j))
	 
	  do k=0,n3-1

	  	n3aux = -0.5 + (1.0)/(real(n3))*(real(k))

	    kpt(counter,1) = (blat1(1)/real(n1))*(real(n1aux))+(blat2(1)/real(n2))*(real(n2aux))&
			     +(blat3(1)/real(n3))*(real(n3aux))+kshift(1)

	    kpt(counter,2) = (blat1(2)/real(n1))*(real(n1aux))+(blat2(2)/real(n2))*(real(n2aux))&
			     +(blat3(2)/real(n3))*(real(n3aux))+kshift(2)

	    kpt(counter,3) =(blat1(3)/real(n1))*(real(n1aux))+(blat2(3)/real(n2))*(real(n2aux))&
			     +(blat3(3)/real(n3))*(real(n3aux))+kshift(3)	
   
	    counter = counter+1

	  end do
	 end do
	end do

end subroutine monhkhorst_pack_adp





















