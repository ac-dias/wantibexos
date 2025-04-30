subroutine spin_z_matrix(w90basis,dft,systype,ovptp,spmz)

	implicit none

	integer :: w90basis,w2

	complex,dimension(w90basis,w90basis) :: spmz,ovptp

	integer :: i,j
	
	character(len=1) :: dft
	character(len=4) :: systype
	
	
	spmz = 0.0
	w2=w90basis/2
	
	select case (systype)
	
	 case ("NP")
	 
	 if (dft .eq. "S") then
	 
	 	spmz= ovptp
	 
	 else
	 
	  do i=1,w90basis

		spmz(i,i) = cmplx(1.,0.)

	  end do
	 
	 end if
	 
	 case ("SP")
	 
	 if (dft .eq. "S") then
	 
	 	spmz(1:w2,1:w2) = cmplx(1.,0.)
	 	
	 	spmz(w2+1:w90basis,w2+1:w90basis) = cmplx(-1.,0.)
	 	
	 	call applyovp(w90basis,ovptp,spmz)
	 
	 else
	 
	  do i=1,w2

		spmz(i,i) = cmplx(1.,0.)
		spmz(w2+i,w2+i) = cmplx(-1.,0.)

	  end do
	 
	 end if
	 
	 case ("SOC")
	 
	  if (dft .eq. "W") then
	  
	   do i=1,w2

		spmz((2*i)-1,(2*i)-1) = cmplx(1.,0.)
		spmz(2*i,2*i) = cmplx(-1.,0.)

	   end do
	   
	   else if (dft .eq. "S") then
	   
	   
	   	spmz(1:w2,1:w2) = cmplx(1.,0.)
	 	
	 	spmz(w2+1:w90basis,w2+1:w90basis) = cmplx(-1.,0.)
	 	
	 	call applyovp(w90basis,ovptp,spmz)
	   
	  
	  else
	  	   
	   do i=1,w2

		spmz(i,i) = cmplx(1.,0.)
		spmz(w2+i,w2+i) = cmplx(-1.,0.)

	   end do
	   
	  
	  end if
	 
	end select	

end subroutine spin_z_matrix


subroutine spin_x_matrix(w90basis,dft,systype,ovptp,spmx)

	implicit none

	integer :: w90basis,w2

	complex,dimension(w90basis,w90basis) :: spmx,ovptp

	integer :: i,j
	
	character(len=1) :: dft
	character(len=4) :: systype
	
	spmx = 0.0
	w2=w90basis/2
	
	select case (systype)
	
	 case ("NP")
	 
	  spmx = 0.0
	 
	 case ("SP")
	 
	 
	 if (dft .eq. "S") then
	 
	 	spmx(1:w2,w2+1:w90basis) = cmplx(1.,0.)
	 	spmx(w2+1:w90basis,1:w2) = cmplx(1.,0.)	 
	 	
	 	call applyovp(w90basis,ovptp,spmx)	
	 
	 
	 else
	 
	  do i=1,w2
	  
	  	spmx(i,w2+i) = cmplx(1.,0.)
	  	spmx(w2+i,i) = cmplx(1.,0.) 
	  
	  end do
	 
	 end if
	 
	 case ("SOC")
	 
	  if (dft .eq. "W") then
	  
	   do i=1,w2
	   
	   	spmx((2*i)-1,2*i) = cmplx(1.,0.)
	   	spmx(2*i,(2*i)-1) = cmplx(1.,0.)	   
	   end do	  

	  else if (dft .eq. "S") then
	  
	  	spmx(1:w2,w2+1:w90basis) = cmplx(1.,0.)
	 	spmx(w2+1:w90basis,1:w2) = cmplx(1.,0.)	 
	 	
	 	call applyovp(w90basis,ovptp,spmx)
	  
	  else
	  
	  do i=1,w2
	  
	  	spmx(i,w2+i) = cmplx(1.,0.)
	  	spmx(w2+i,i) = cmplx(1.,0.) 
	  
	  end do
	  
	  
	  end if
	 
	end select	

end subroutine spin_x_matrix


subroutine spin_y_matrix(w90basis,dft,systype,ovptp,spmy)

	implicit none

	integer :: w90basis,w2

	complex,dimension(w90basis,w90basis) :: spmy,ovptp

	integer :: i,j
	
	character(len=1) :: dft
	character(len=4) :: systype
	
	spmy = 0.0
	w2=w90basis/2	
	
	select case (systype)
	
	 case ("NP")
	 
	  spmy = 0.0
	 
	 case ("SP")
	 
	 
	 if (dft .eq. "S") then
	 
	 	spmy(1:w2,w2+1:w90basis) = cmplx(0.,-1.)
	 	spmy(w2+1:w90basis,1:w2) = cmplx(0.,1.)
	 
	 	call applyovp(w90basis,ovptp,spmy)
	 
	 else
	 
	  do i=1,w2
	  
	  	spmy(i,w2+i) = cmplx(0.,-1.)
	  	spmy(w2+i,i) = cmplx(0.,1.) 
	  
	  end do
	  
	  end if
	 
	 case ("SOC")
	 
	  if (dft .eq. "W") then

	   do i=1,w2
	   
	   	spmy((2*i)-1,2*i) = cmplx(0.,-1.)
	   	spmy(2*i,(2*i)-1) = cmplx(0.,1.)	   
	   end do
	  
	  else if (dft .eq. "S") then
	  
	 	spmy(1:w2,w2+1:w90basis) = cmplx(0.,-1.)
	 	spmy(w2+1:w90basis,1:w2) = cmplx(0.,1.)
	 
	 	call applyovp(w90basis,ovptp,spmy)
	 		  
	  else
	  
	  do i=1,w2
	  
	  	spmy(i,w2+i) = cmplx(0.,-1.)
	  	spmy(w2+i,i) = cmplx(0.,1.)
	  
	  end do	  
	  
	  
	  end if
	 
	end select	

end subroutine spin_y_matrix



subroutine spinvl(w90basis,vec,dft,systype,ovptp,spx,spy,spz)


	implicit none

	integer :: w90basis,w2

	complex,dimension(w90basis) :: vec
	complex,dimension(w90basis) :: vcconj

	complex,dimension(w90basis) :: zvvaux,yvvaux,xvvaux

	complex,dimension(w90basis,w90basis) :: spmz,spmx,spmy,ovptp

	complex :: spz,spx,spy

	integer :: i,j
	
	character(len=1) :: dft
	character(len=4) :: systype	

	spmz=0.0

	w2=w90basis/2
	

	call spin_z_matrix(w90basis,dft,systype,ovptp,spmz)
	call spin_x_matrix(w90basis,dft,systype,ovptp,spmx)		
	call spin_y_matrix(w90basis,dft,systype,ovptp,spmy)		
	


	call vecconjg(vec,w90basis,vcconj)

	!calculando valor medio de sz

	call matvec(spmz,vec,w90basis,zvvaux)

	call prodintsq(vcconj,zvvaux,w90basis,spz)
	
	!calculando valor medio de sx

	call matvec(spmx,vec,w90basis,xvvaux)

	call prodintsq(vcconj,xvvaux,w90basis,spx)
	
	!calculando valor medio de sy

	call matvec(spmy,vec,w90basis,yvvaux)

	call prodintsq(vcconj,yvvaux,w90basis,spy)		


end subroutine spinvl
