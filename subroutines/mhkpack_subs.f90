
subroutine rkmesh(rk,rlat,ngrid)

	implicit none
	
	integer,dimension(3) :: ngrid
	double precision :: rk
	double precision,dimension(3,3) :: rlat,blat
	double precision,dimension(3) :: vsize
	double precision,parameter :: pi=acos(-1.0)

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
	double precision :: rk
	double precision,dimension(3,3) :: rlat,blat
	double precision,dimension(3) :: vsize
	double precision,parameter :: pi=acos(-1.0)

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat(1,:),blat(2,:),blat(3,:))

	call vecsize(blat(1,:),vsize(1))
	call vecsize(blat(2,:),vsize(2))
	call vecsize(blat(3,:),vsize(3))

	vsize = vsize/(2*pi)

	ngrid(1) = int(max(1.0,(rk*vsize(1))+0.5))
	ngrid(2) = int(max(1.0,(rk*vsize(2))+0.5))
	ngrid(3) = 1



end subroutine rkmesh2D

subroutine kpath2(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks/2)*npts,4) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux


	
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


end subroutine kpath2

subroutine kpathbse(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks/2)*npts,4) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux


	
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

subroutine kpath(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension((nks-1)*npts,4) :: kpt

	double precision,dimension(nks,3) :: ks
	double precision,dimension(nks,3) :: ksaux

	double precision,dimension(nks-1) :: kdis !distancia entre os pontos k do caminho
	double precision :: kdistot !soma do caminho todo
	double precision :: kdisaux

	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS-BSE.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS-BSE output file"

	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do


	kdistot= 0.0
	do i=1,nks-1

		call distvec(ksaux(i,:),ksaux(i+1,:),kdis(i))
		kdistot=kdistot+kdis(i)
	end do

	kdis=kdis/kdistot


	do i=1,npts

		kpt(i,1) = (kdis(1))*(i-1.)/(npts-1.)
		kpt(i,2) = ksaux(1,1)+ (ksaux(2,1)-ksaux(1,1))*(i-1.)/(npts-1.)
		kpt(i,3) = ksaux(1,2)+ (ksaux(2,2)-ksaux(1,2))*(i-1.)/(npts-1.)
		kpt(i,4) = ksaux(1,3)+ (ksaux(2,3)-ksaux(1,3))*(i-1.)/(npts-1.)

	end do

	if (nks .gt. 2) then

	kdisaux=0.0
	write(1733,*) real(kdisaux)
	 do j=2,nks-1
		kdisaux=kdisaux+kdis(j-1)
		write(1733,*) real(kdisaux)
	  do i=1,npts

		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i)/(npts)
		kpt((j-1)*npts+i,2) = ksaux(j,1)+ (ksaux(j+1,1)-ksaux(j,1))*(i)/(npts)
		kpt((j-1)*npts+i,3) = ksaux(j,2)+ (ksaux(j+1,2)-ksaux(j,2))*(i)/(npts)
		kpt((j-1)*npts+i,4) = ksaux(j,3)+ (ksaux(j+1,3)-ksaux(j,3))*(i)/(npts)
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)	

	  end do
	 end do

	else 

	write(1733,*) real(0.0000)
	write(1733,*) real(1.0000)
	 continue

	end if


end subroutine kpath

subroutine distvec(vec1,vec2,distvecout) !subroutina para calcular a distância entre dois vetores

	implicit none

	double precision,dimension(3) :: vec1,vec2
	double precision :: distvecout

	distvecout=sqrt((vec1(1)-vec2(1))**2+(vec1(2)-vec2(2))**2+(vec1(3)-vec2(3))**2)


end subroutine distvec


subroutine monhkhorst_packq(q,n1,n2,n3,shift,rlat1,rlat2,rlat3,qpt)

	implicit none
	
	integer :: n1,n2,n3
	double precision,dimension(3) :: shift,kshift
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(n1*n2*n3,3) :: qpt

	double precision,dimension(4) :: q

	integer :: i,j,k, counter
	double precision,dimension(3) :: blat1,blat2,blat3

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
	double precision,dimension(3) :: shift,kshift
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(n1*n2*n3,3) :: kpt

	integer :: i,j,k, counter
	double precision,dimension(3) :: blat1,blat2,blat3

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

	    kpt(counter,1) = (blat1(1)/dble(n1))*(dble(i))+(blat2(1)/dble(n2))*(dble(j))&
			     +(blat3(1)/dble(n3))*(dble(k))+kshift(1)

	    kpt(counter,2) = (blat1(2)/dble(n1))*(dble(i))+(blat2(2)/dble(n2))*(dble(j))&
			     +(blat3(2)/dble(n3))*(dble(k))+kshift(2)

	    kpt(counter,3) =(blat1(3)/dble(n1))*(dble(i))+(blat2(3)/dble(n2))*(dble(j))&
			     +(blat3(3)/dble(n3))*(dble(k))+kshift(3)	
   
	    counter = counter+1

	  end do
	 end do
	end do

end subroutine monhkhorst_pack

subroutine recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3) !calcula os vetores da rede recíproca, a partir dos vetores da rede real

	implicit none
	
	double precision,parameter :: pi=acos(-1.)
	
	double precision,dimension(3) :: rlat1,rlat2,rlat3
	double precision,dimension(3) :: blat1,blat2,blat3

	double precision,dimension(3) :: v23,v31,v12
	double precision :: vol

	call prodvec(rlat2,rlat3,v23)
	call prodvec(rlat3,rlat1,v31)
	call prodvec(rlat1,rlat2,v12)

	vol= abs((rlat1(1)*v23(1))+(rlat1(2)*v23(2))+(rlat1(3)*v23(3)))

	blat1 = ((2.0*pi)/vol)*v23
	blat2 = ((2.0*pi)/vol)*v31
	blat3 = ((2.0*pi)/vol)*v12


end subroutine recvec



!subrotinas para gerar grid hexagonal

subroutine gridgenhexcount(ngrid,a,counter)


	implicit none

	double precision,parameter :: pi=acos(-1.)

	integer :: ngrid,i,j

	integer :: counter

	double precision :: a!constante da rede

	double precision :: a0

        double precision :: xmax, ymax

	double precision :: kx,ky


	a0 = a/sqrt(3.)

        xmax=2*pi/(3.d0*a0)
        ymax=4*pi/(3.d0*a)


	counter = 0


	do i=1,ngrid
                kx=xmax*(i-1)/(ngrid-1)
		do j=1,ngrid
                         ky=ymax*(j-1)/(ngrid-1)
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


end subroutine gridgenhexcount


subroutine gridgenhex(ngrid,a,ncounter,kg)


	implicit none

	double precision,parameter :: pi=acos(-1.)

	integer :: ngrid,i,j

	integer :: counter,ncounter

	double precision :: a!constante da rede

	double precision :: a0

        double precision :: xmax, ymax

	double precision :: kx,ky

	double precision,dimension(ncounter,2) :: kg


	a0 = a/sqrt(3.)

        xmax=2*pi/(3.d0*a0)
        ymax=4*pi/(3.d0*a)

	counter = 0

	do i=1,ngrid
                kx=xmax*(i-1)/(ngrid-1)
		do j=1,ngrid
                         ky=ymax*(j-1)/(ngrid-1)
                         if (ky.lt.ymax/2.d0) then
                            counter=counter+1 
                            kg(counter, 2)=kx
                            kg(counter, 1)=ky
                            if (i.ne.1.and.j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=-ky
                            end if
                            if (i.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=ky
			    end if
	                    if (j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=kx
                                kg(counter, 1)=-ky
                            end if
                         else if (ky.lt.(ymax-kx*tan(pi/6.d0))) then
                            counter=counter+1 
                            kg(counter, 2)=kx
                            kg(counter, 1)=ky
                            if (i.ne.1.and.j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=-ky
                            end if
                            if (i.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=-kx
                                kg(counter, 1)=ky
			    end if
	                    if (j.ne.1) then
                                counter=counter+1 
                                kg(counter, 2)=kx
                                kg(counter, 1)=-ky
                            end if
                         end if
                
		
		end do
        end do

	


end subroutine gridgenhex
