subroutine second_order_mass_tensor(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
                                    ihopmatrices,hxx,hxy,hxz,hyy,hyz,hzz) !subroutine that takes second order derivative of MLWF-TB Hamiltonian for a given k-point

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	real :: kx,ky,kz

	integer,dimension(nvec) :: ffactor

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec

	complex,dimension(w90basis,w90basis) :: hxx,hxy,hxz,hyy,hyz,hzz

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices


	complex,dimension(nvec) :: kvecs

	real,dimension(nvec) :: kvecsxx,kvecsxy,kvecsxz
	real,dimension(nvec) :: kvecsyy,kvecsyz,kvecszz
	
	complex,parameter :: imag=cmplx(0.0,1.0)	


		do i=1,nvec

		kvecs(i) = cmplx(0.0,kx*rvec(i,1))+cmplx(0.0,ky*rvec(i,2))+cmplx(0.0,kz*rvec(i,3))
			   
		
		kvecsxx(i) = -rvec(i,1)*rvec(i,1)
		              		              
		kvecsxy(i) = -rvec(i,1)*rvec(i,2)		              

		kvecsxz(i) = -rvec(i,1)*rvec(i,3)		              
		              
		kvecsyy(i) = -rvec(i,2)*rvec(i,2)
		              
		kvecsyz(i) = -rvec(i,2)*rvec(i,3)		              

		kvecszz(i) = -rvec(i,3)*rvec(i,3)		              
		              
		end do


		hxx = 0.0
		hxy = 0.0
		hxz = 0.0
		hyy = 0.0
		hyz = 0.0
		hzz = 0.0


		do i=1,nvec
		

			hxx = hxx+kvecsxx(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))
			hxy = hxy+kvecsxy(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))
			hxz = hxz+kvecsxz(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))
			hyy = hyy+kvecsyy(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))
			hyz = hyz+kvecsyz(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))
			hzz = hzz+kvecszz(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))
							

		
		end do


end subroutine second_order_mass_tensor


subroutine effective_mass_calculator(dft,systype,kx,ky,kz,ffactor,w90basis,wf,nvec,rlat,rvec,hopmatrices,&
                                     ihopmatrices,ovp,exc,mag,mtensor) !calculates effective mass tensor for a given kpoint and specific electronic state

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	real :: kx,ky,kz,exc
	
	real,dimension(3) :: mag

	integer,dimension(nvec) :: ffactor

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec
	
	complex,dimension(w90basis,w90basis) :: hxx,hxy,hxz,hyy,hyz,hzz
	
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices,ovp
	
	real :: c !constant to get non-dimensional effective mass
	
	complex,dimension(w90basis) :: wf !electronic state eigenvector
	
	complex,dimension(3,3) :: mtensoraux 
	complex,dimension(3,3) :: mtensor !effective mass tensor
	
	character(len=1) :: dft
	character(len=4) :: systype 
	
	
	call second_order_mass_tensor(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
                                    ihopmatrices,hxx,hxy,hxz,hyy,hyz,hzz)
                                    
        call sandwich(w90basis,wf,hxx,wf,mtensoraux(1,1))
        call sandwich(w90basis,wf,hxy,wf,mtensoraux(1,2))
        call sandwich(w90basis,wf,hxz,wf,mtensoraux(1,3))
        
        mtensoraux(2,1) = conjg(mtensoraux(1,2))
        call sandwich(w90basis,wf,hyy,wf,mtensoraux(2,2))
        call sandwich(w90basis,wf,hyz,wf,mtensoraux(2,3))
        
        mtensoraux(3,1) = conjg(mtensoraux(1,3))
        mtensoraux(3,2) = conjg(mtensoraux(2,3))
        call sandwich(w90basis,wf,hzz,wf,mtensoraux(3,3))
        
        c= (9.109383)/((6.582119)*(6.582119)*(1.602177))       
	

	do i=1,3
	
		do j=1,3
		
			mtensor(i,j) = 1.0/(c*mtensoraux(i,j))
		
		end do
	
	end do	

end subroutine effective_mass_calculator


subroutine num_effective_mass_calculator(dft,systype,nthreads,scs,exc,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
                                     ihopmatrices,ovp,efermi,eigvout,n,h,mtensor,nocpf,fermishift,mag) !calculates effective mass tensor for a given kpoint and specific electronic
                                     
                                     
 	implicit none

	integer :: i,j,k,nthreads,nocp

	integer :: w90basis,nvec

	real :: kx,ky,kz
	
	real :: scs,exc,efermi

	integer,dimension(nvec) :: ffactor

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec
	
	real,dimension(3) :: mag
	
	complex,dimension(w90basis,w90basis) :: hxx,hxy,hxz,hyy,hyz,hzz
	
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices,ovp
	
	real :: c !constant to get non-dimensional effective mass
	

	real,dimension(3,3) :: mtensoraux 
	real,dimension(3,3) :: mtensor !effective mass tensor 
	
	real,dimension(w90basis) :: eigv,eigvout 
	
	complex,dimension(w90basis,w90basis) :: autovetores                                  
                                     
        real :: energy
        
        real :: h !kpoint step
        
        integer :: n !band number
        
        real,dimension(5) :: enxx,enyy,enzz
        real,dimension(16) :: enxy,enxz,enyz
        
        real,dimension(5,3) :: kptaa
        real,dimension(16,3) :: kptab
        
        real,dimension(5) :: factorsaa
        real,dimension(16) :: factorsab
        
	real :: caux,daux
	
	integer :: nocpf
	real :: fermishift
	
	character(len=1) :: dft
	character(len=4) :: systype
	                                      
        c= (9.109383)/((6.582119)*(6.582119)*(1.602177))    
        
        mtensoraux = 0.0  
        
        !defining factorsaa 
        
        	caux = 1.0/(12.0*h*h)
        
        	factorsaa(1) = caux*(-1.0)
        	factorsaa(2) = caux*(-1.0)
        	factorsaa(3) = caux*(16.0)
        	factorsaa(4) = caux*(16.0)        	
        	factorsaa(5) = caux*(-30.0)
        	        	        	        
        !defining factorsab 
        
        	daux = 1.0/(600.0*h*h)         
                            
        	factorsab(1) = daux*(-63.0)
        	factorsab(2) = daux*(-63.0)
        	factorsab(3) = daux*(-63.0)
        	factorsab(4) = daux*(-63.0)      	
        	factorsab(5) = daux*(63.0)
        	factorsab(6) = daux*(63.0)
        	factorsab(7) = daux*(63.0)
        	factorsab(8) = daux*(63.0)
        	factorsab(9) = daux*(44.00)        	
        	factorsab(10) = daux*(44.00)
        	factorsab(11) = daux*(-44.00)
        	factorsab(12) = daux*(-44.00)
        	factorsab(13) = daux*(74.00)
        	factorsab(14) = daux*(74.00)        	
        	factorsab(15) = daux*(-74.00)
        	factorsab(16) = daux*(-74.00)
        	        	       	                                     
        !defining kpoints for xx 
        
        kptaa(1,1) = kx-2.0*h
        kptaa(1,2) = ky
        kptaa(1,3) = kz                                   

        kptaa(2,1) = kx+2.0*h
        kptaa(2,2) = ky
        kptaa(2,3) = kz                                   
        
        kptaa(3,1) = kx-1.0*h
        kptaa(3,2) = ky
        kptaa(3,3) = kz                                                                               

        kptaa(4,1) = kx+1.0*h
        kptaa(4,2) = ky
        kptaa(4,3) = kz                                    
                                     
        kptaa(5,1) = kx
        kptaa(5,2) = ky
        kptaa(5,3) = kz
        
        do i=1,5
        	call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kptaa(i,1),kptaa(i,2),kptaa(i,3),w90basis,nvec,&
        	            rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		            
		            enxx(i) = eigv(n)
		            
		            mtensoraux(1,1) = mtensoraux(1,1)+factorsaa(i)*enxx(i)
        end do
        
        


        !defining kpoints for yy 
        
        kptaa(1,1) = kx
        kptaa(1,2) = ky-2.0*h
        kptaa(1,3) = kz                                   

        kptaa(2,1) = kx
        kptaa(2,2) = ky+2.0*h
        kptaa(2,3) = kz                                   
        
        kptaa(3,1) = kx
        kptaa(3,2) = ky-1.0*h
        kptaa(3,3) = kz                                                                               

        kptaa(4,1) = kx
        kptaa(4,2) = ky+1.0*h
        kptaa(4,3) = kz                                    
                                     
        kptaa(5,1) = kx
        kptaa(5,2) = ky
        kptaa(5,3) = kz
        
        do i=1,5
        	call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kptaa(i,1),kptaa(i,2),kptaa(i,3),w90basis,nvec,&
        	            rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		            
		            enyy(i) = eigv(n)
		            
		            mtensoraux(2,2) = mtensoraux(2,2)+factorsaa(i)*enyy(i)
        end do
        
        !defining kpoints for zz 
        
        kptaa(1,1) = kx
        kptaa(1,2) = ky
        kptaa(1,3) = kz-2.0*h                                   

        kptaa(2,1) = kx
        kptaa(2,2) = ky
        kptaa(2,3) = kz+2.0*h                                   
        
        kptaa(3,1) = kx
        kptaa(3,2) = ky
        kptaa(3,3) = kz-1.0*h                                                                               

        kptaa(4,1) = kx
        kptaa(4,2) = ky
        kptaa(4,3) = kz+1.0*h                                   
                                     
        kptaa(5,1) = kx
        kptaa(5,2) = ky
        kptaa(5,3) = kz                 

        do i=1,5
        	call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kptaa(i,1),kptaa(i,2),kptaa(i,3),w90basis,nvec,&
        	            rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		            
		            enzz(i) = eigv(n)
		            
		            mtensoraux(3,3) = mtensoraux(3,3)+factorsaa(i)*enzz(i)
		            
		 if (i .eq. 5) then
		 
		 	  eigvout = eigv
		 
		 end if
        end do
          
        !defining kpoints for xy 
        
        kptab(1,1) = kx+h  !a
        kptab(1,2) = ky-2.0*h  !a   
        kptab(1,3) = kz                                   

        kptab(2,1) = kx+2.0*h  !a
        kptab(2,2) = ky-h  !a
        kptab(2,3) = kz                                   
        
        kptab(3,1) = kx-2.0*h  !a
        kptab(3,2) = ky+h  !a
        kptab(3,3) = kz                                                                               

        kptab(4,1) = kx-h  !a
        kptab(4,2) = ky+2.0*h  !a
        kptab(4,3) = kz                                    
                                     
        kptab(5,1) = kx-h  !a
        kptab(5,2) = ky-2.0*h  !a
        kptab(5,3) = kz 
        
        kptab(6,1) = kx-2.0*h  !a
        kptab(6,2) = ky-h  !a
        kptab(6,3) = kz                                   

        kptab(7,1) = kx+h  !a
        kptab(7,2) = ky+2.0*h  !a
        kptab(7,3) = kz                                   
        
        kptab(8,1) = kx+2.0*h  !a
        kptab(8,2) = ky+h  !a
        kptab(8,3) = kz                                                                               

        kptab(9,1) = kx+2.0*h  !a
        kptab(9,2) = ky-2.0*h  !a
        kptab(9,3) = kz                                    
                                     
        kptab(10,1) = kx-2.0*h  !a
        kptab(10,2) = ky+2.0*h  !a
        kptab(10,3) = kz 
        
        kptab(11,1) = kx-2.0*h  !a
        kptab(11,2) = ky-2.0*h  !a
        kptab(11,3) = kz                                   

        kptab(12,1) = kx+2.0*h  !a
        kptab(12,2) = ky+2.0*h  !a
        kptab(12,3) = kz                                   
        
        kptab(13,1) = kx-h  !a
        kptab(13,2) = ky-h  !a
        kptab(13,3) = kz                                                                               

        kptab(14,1) = kx+h  !a
        kptab(14,2) = ky+h  !a
        kptab(14,3) = kz                                    
                                     
        kptab(15,1) = kx+h  !a
        kptab(15,2) = ky-h  !a
        kptab(15,3) = kz                 
        
        kptab(16,1) = kx-h  !a
        kptab(16,2) = ky+h  !a
        kptab(16,3) = kz          

        do i=1,16
        	call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kptab(i,1),kptab(i,2),kptab(i,3),w90basis,nvec,&
        	            rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		            
		            enxy(i) = eigv(n)
		            
		            mtensoraux(1,2) = mtensoraux(1,2)+factorsab(i)*enxy(i)
        end do
        
        !defining kpoints for xz 
        
        kptab(1,1) = kx+h  !a
        kptab(1,2) = ky
        kptab(1,3) = kz-2.0*h  !a                                  

        kptab(2,1) = kx+2.0*h  !a
        kptab(2,2) = ky
        kptab(2,3) = kz-h  !a                                   
        
        kptab(3,1) = kx-2.0*h  !a
        kptab(3,2) = ky
        kptab(3,3) = kz+h  !a                                                                               

        kptab(4,1) = kx-h  !a
        kptab(4,2) = ky
        kptab(4,3) = kz+2.0*h  !a                                    
                                     
        kptab(5,1) = kx-h  !a
        kptab(5,2) = ky
        kptab(5,3) = kz-2.0*h  !a 
        
        kptab(6,1) = kx-2.0*h  !a
        kptab(6,2) = ky
        kptab(6,3) = kz-h  !a                                   

        kptab(7,1) = kx+h  !a
        kptab(7,2) = ky
        kptab(7,3) = kz+2.0*h  !a                                  
        
        kptab(8,1) = kx+2.0*h  !a
        kptab(8,2) = ky
        kptab(8,3) = kz+h  !a                                                                               

        kptab(9,1) = kx+2.0*h  !a
        kptab(9,2) = ky
        kptab(9,3) = kz-2.0*h  !a                                    
                                     
        kptab(10,1) = kx-2.0*h  !a
        kptab(10,2) = ky
        kptab(10,3) = kz+2.0*h  !a 
        
        kptab(11,1) = kx-2.0*h  !a
        kptab(11,2) = ky
        kptab(11,3) = kz-2.0*h  !a                                   

        kptab(12,1) = kx+2.0*h  !a
        kptab(12,2) = ky
        kptab(12,3) = kz+2.0*h  !a                                  
        
        kptab(13,1) = kx-h  !a
        kptab(13,2) = ky
        kptab(13,3) = kz-h  !a                                                                               

        kptab(14,1) = kx+h  !a
        kptab(14,2) = ky
        kptab(14,3) = kz+h  !a                                    
                                     
        kptab(15,1) = kx+h  !a
        kptab(15,2) = ky
        kptab(15,3) = kz-h  !a                 
        
        kptab(16,1) = kx-h  !a
        kptab(16,2) = ky
        kptab(16,3) = kz+h  !a               

        do i=1,16
         	call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kptab(i,1),kptab(i,2),kptab(i,3),w90basis,nvec,&
        	            rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		            
		            enxz(i) = eigv(n)
		            
		            mtensoraux(1,3) = mtensoraux(1,3)+factorsab(i)*enxz(i)
        end do        
        
        !defining kpoints for yz 
        
        kptab(1,1) = kx
        kptab(1,2) = ky+h  !a
        kptab(1,3) = kz-2.0*h  !a                                  

        kptab(2,1) = kx
        kptab(2,2) = ky+2.0*h  !a
        kptab(2,3) = kz-h  !a                                  
        
        kptab(3,1) = kx
        kptab(3,2) = ky-2.0*h  !a
        kptab(3,3) = kz+h  !a                                                                               

        kptab(4,1) = kx
        kptab(4,2) = ky-h  !a
        kptab(4,3) = kz+2.0*h  !a                                    
                                     
        kptab(5,1) = kx
        kptab(5,2) = ky-h  !a
        kptab(5,3) = kz-2.0*h  !a 
        
        kptab(6,1) = kx  
        kptab(6,2) = ky-2.0*h  !a
        kptab(6,3) = kz-h  !a                                 

        kptab(7,1) = kx
        kptab(7,2) = ky+h  !a
        kptab(7,3) = kz+2.0*h  !a                                   
        
        kptab(8,1) = kx
        kptab(8,2) = ky+2.0*h  !a
        kptab(8,3) = kz+h  !a                                                                               

        kptab(9,1) = kx
        kptab(9,2) = ky+2.0*h  !a
        kptab(9,3) = kz-2.0*h  !a                                    
                                     
        kptab(10,1) = kx
        kptab(10,2) = ky-2.0*h  !a
        kptab(10,3) = kz+2.0*h  !a 
        
        kptab(11,1) = kx
        kptab(11,2) = ky-2.0*h  !a
        kptab(11,3) = kz-2.0*h  !a                                  

        kptab(12,1) = kx
        kptab(12,2) = ky+2.0*h  !a
        kptab(12,3) = kz+2.0*h  !a                                   
        
        kptab(13,1) = kx
        kptab(13,2) = ky-h  !a
        kptab(13,3) = kz-h  !a                                                                               

        kptab(14,1) = kx
        kptab(14,2) = ky+h  !a
        kptab(14,3) = kz+h  !a                                    
                                     
        kptab(15,1) = kx
        kptab(15,2) = ky+h  !a
        kptab(15,3) = kz-h  !a                 
        
        kptab(16,1) = kx
        kptab(16,2) = ky-h  !a
        kptab(16,3) = kz+h  !a               
        
        do i=1,16
        	call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kptab(i,1),kptab(i,2),kptab(i,3),w90basis,nvec,&
        	            rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		            
		            enyz(i) = eigv(n)
		            
		            mtensoraux(2,3) = mtensoraux(2,3)+factorsab(i)*enyz(i)
        end do          
        
        
        
     mtensoraux(2,1) = mtensoraux(1,2)
     mtensoraux(3,1) = mtensoraux(1,3)         
     mtensoraux(3,2) = mtensoraux(2,3)       
        
        
        
   	do i=1,3
	
		do j=1,3
		
			mtensor(i,j) = 1.0/(c*mtensoraux(i,j))
		
		end do
	
	end do	     
        
        
                                                                                                                 
end subroutine num_effective_mass_calculator                                   
