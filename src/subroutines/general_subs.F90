subroutine gapfinder(w90basis,ngkpt,nocpj,ebands,kptcbm,kptvbm,kptgap,gapdir,cbm,vbm) !subroutine to find direct and indirect bandgap

	implicit none
	integer :: i,j
	integer :: w90basis,ngkpt
	real :: ebands(w90basis,ngkpt) 
	integer :: kptcbm,kptvbm,kptgap
	real:: gapdir,cbm,vbm
	integer,dimension(ngkpt) :: nocpj
	
	
	
	!searching CBM and VBM
	cbm= 500
	vbm= -500
	gapdir= 500
	kptvbm= 0
	kptcbm= 0

	do j=1,ngkpt
	
	
		if (ebands(nocpj(j),j) .gt. vbm ) then
		
			vbm = ebands(nocpj(j),j)
			kptvbm = j
		
		else
			continue
		end if
		
		
		if (ebands(nocpj(j)+1,j) .lt. cbm ) then
		
			cbm = ebands(nocpj(j)+1,j)
			kptcbm = j
		
		else
			continue
		end if		
		
		
		if ( (ebands(nocpj(j)+1,j)-ebands(nocpj(j),j)) .lt. gapdir ) then
		
			gapdir = ebands(nocpj(j)+1,j)-ebands(nocpj(j),j)
			kptgap = j
		
		else
			continue
		end if						
	
	
	end do



end subroutine gapfinder



subroutine quantumnumbers(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

	integer :: ngkpt,nc,nv,w90basis
	integer,dimension(ngkpt) :: nocpk,nocpq
	integer,dimension(ngkpt*nc*nv,4) :: stt

	integer :: i,c,v

	integer :: counter, ncaux,nvaux


	counter=1

	do i=1,ngkpt

	ncaux=(nocpq(i)+1)+nc-1

	nvaux=nocpk(i)-nv+1


				if (nv .gt. nocpk(i)) then

					write(*,*) "Number of valence states higher then ocuppied states"
					STOP

				else if (nc .gt. (w90basis-nocpq(i))) then

					write(*,*) "Number of conduction states higher then unocuppied states"
					STOP

				else 

					continue
					
				end if


		do c=(nocpq(i)+1),ncaux


			do v=nvaux,nocpk(i)



				stt(counter,1) = counter !numero do estado

				stt(counter,2) = v   !numero da banda da valencia

				stt(counter,3) = c   !numero da banda da condução

				stt(counter,4) = i   !numero do ponto k no grid

				counter=counter+1
 			

			end do



		end do



	end do


end subroutine quantumnumbers

subroutine sandwich(w90basis,lvec,hm,rvec,res)

	implicit none

	integer :: w90basis

	complex,dimension(w90basis) :: lvec,rvec,clvec,auxvec
	complex,dimension(w90basis,w90basis) :: hm

	complex :: res

	call matvec(hm,rvec,w90basis,auxvec)
	call vecconjg(lvec,w90basis,clvec)

	call prodintsq(clvec,auxvec,w90basis,res)


end subroutine sandwich


subroutine sandwich_average(w90basis,lvec,hma,hmb,rvec,res)

	implicit none

	integer :: w90basis,i,j
	complex,dimension(w90basis) :: lvec,rvec
	complex,dimension(w90basis,w90basis) :: hma,hmb
	complex :: res,aux

	! Evaluate conjg(lvec)^T [(hma + hmb) / 2] rvec directly.  Passing
	! 0.5*(hma+hmb) to sandwich creates a full w90basis-by-w90basis
	! compiler temporary for every BSE Hamiltonian element.
	res = 0.0
	do i=1,w90basis
		aux = 0.0
		do j=1,w90basis
			aux = aux + 0.5*(hma(i,j)+hmb(i,j))*rvec(j)
		end do
		res = res + conjg(lvec(i))*aux
	end do

end subroutine sandwich_average

subroutine applyovp(w90basis,ovp,hm)

	implicit none
	
	integer :: i,j,w90basis
	complex,dimension(w90basis,w90basis) :: ovp,hm
	
	do i=1,w90basis
	
		do j=1,w90basis
		
			hm(i,j) = hm(i,j)*ovp(i,j)
		
		end do
	
	end do
	

end subroutine applyovp


subroutine prodintsq(vetora,vetorb,n,resultado)

	implicit none

	integer:: n,i
	complex :: resultado1, auxi
	complex, dimension(n) :: vetora,vetorb,vetorc
	complex :: resultado

	

	resultado1=0

	do i=1,n

		auxi=vetora(i)*vetorb(i)

		resultado1=resultado1+auxi

	end do

	resultado=resultado1

end subroutine prodintsq

subroutine vecconjg(vetor,n,vetout)

	implicit none

	integer :: i,j,n
	complex, dimension(n) :: vetor,vetout


	do i=1,n

		vetout(i)=conjg(vetor(i))

	end do



end subroutine vecconjg 

subroutine matvec(matriz,vetor,n,veout)

	implicit none

	integer :: n !ordem da matriz e tamanho do vetor
	complex, dimension(n,n) :: matriz
	complex, dimension(n):: vetor,veout
	complex :: flag
	integer :: i,j

	veout=0
	
	do i=1,n


		do j=1,n

			flag=matriz(i,j)*vetor(j)

			veout(i)=flag+veout(i)
		

		end do
	

	end do


end subroutine matvec

SUBROUTINE invers(h11inv,n1)

	integer n1
	complex h11inv(n1,n1)
        integer    lda1,lwk1,if1,ipiv1(n1)
        complex wk1(64*n1)
	lda1     =n1
	lwk1     =64*n1

        call zgetrf(n1,n1,h11inv,lda1,ipiv1,if1)
  	if (if1.eq.0) then
	call zgetri(n1,h11inv,lda1,ipiv1,wk1,lwk1,if1)
	else
	write(*,*)'error in routine invers'
	endif

   

END SUBROUTINE invers

subroutine prodvec(v1,v2,vx) !calcula o vetor oriundo do produto vetorial entre dois vetores v1 X v2

	implicit none
	
	real,dimension(3) :: v1,v2,vx

	vx(1) = (v1(2)*v2(3))-(v1(3)*v2(2))
	vx(2) = (v1(3)*v2(1))-(v1(1)*v2(3))
	vx(3) = (v1(1)*v2(2))-(v1(2)*v2(1))


end subroutine prodvec

subroutine vcell2D(rlat,vc)

	real,dimension(3,3) :: rlat
	real :: vc

	real,dimension(3) :: vx
	real :: aux

		vx(1)=rlat(1,2)*rlat(2,3)-rlat(1,3)*rlat(2,2)
		vx(2)=rlat(1,3)*rlat(2,1)-rlat(1,1)*rlat(2,3)
		vx(3)=rlat(1,1)*rlat(2,2)-rlat(1,2)*rlat(2,1)
		call vecsize(vx,vc)

end subroutine vcell2D

subroutine vcell3D(rlat,vc)

	real,dimension(3,3) :: rlat
	real :: vc

	real,dimension(3) :: vx
	real :: aux

		vx(1)=rlat(1,2)*rlat(2,3)-rlat(1,3)*rlat(2,2)
		vx(2)=rlat(1,3)*rlat(2,1)-rlat(1,1)*rlat(2,3)
		vx(3)=rlat(1,1)*rlat(2,2)-rlat(1,2)*rlat(2,1)

		vc=vx(1)*rlat(3,1)+vx(2)*rlat(3,2)+vx(3)*rlat(3,3)
		vc = abs(vc)

end subroutine vcell3D

subroutine vcell(systype,rlat,vc)

	character(len=2) :: systype
	real,dimension(3,3) :: rlat
	real :: vc

	real,dimension(3) :: vx
	real :: aux

	if (systype .eq. "2D") then

		vx(1)=rlat(1,2)*rlat(2,3)-rlat(1,3)*rlat(2,2)
		vx(2)=rlat(1,3)*rlat(2,1)-rlat(1,1)*rlat(2,3)
		vx(3)=rlat(1,1)*rlat(2,2)-rlat(1,2)*rlat(2,1)
		call vecsize(vx,vc)
		
	else
		vx(1)=rlat(1,2)*rlat(2,3)-rlat(1,3)*rlat(2,2)
		vx(2)=rlat(1,3)*rlat(2,1)-rlat(1,1)*rlat(2,3)
		vx(3)=rlat(1,1)*rlat(2,2)-rlat(1,2)*rlat(2,1)

		vc=vx(1)*rlat(3,1)+vx(2)*rlat(3,2)+vx(3)*rlat(3,3)
		vc = abs(vc)

	end if


end subroutine vcell

subroutine alat2D(rlat,a0)

	real,dimension(3,3) :: rlat
	real :: a0
	real,dimension(3) :: vsize

	vsize(1)=sqrt(rlat(1,1)**2+rlat(1,2)**2+rlat(1,3)**2)
	vsize(2)=sqrt(rlat(2,1)**2+rlat(2,2)**2+rlat(2,3)**2)
	vsize(3)=sqrt(rlat(3,1)**2+rlat(3,2)**2+rlat(3,3)**2)

		a0 = 0.5*(vsize(1)+vsize(2))

end subroutine alat2D

subroutine alat(systype,rlat,a0)

	character(len=2) :: systype
	real,dimension(3,3) :: rlat
	real :: a0

	real,dimension(3) :: vsize

	vsize(1)=sqrt(rlat(1,1)**2+rlat(1,2)**2+rlat(1,3)**2)
	vsize(2)=sqrt(rlat(2,1)**2+rlat(2,2)**2+rlat(2,3)**2)
	vsize(3)=sqrt(rlat(3,1)**2+rlat(3,2)**2+rlat(3,3)**2)


	if (systype .eq. "2D") then

		a0 = 0.5*(vsize(1)+vsize(2))
	else
		a0 = (1.0/3.0)*(vsize(1)+vsize(2)+vsize(3))

	end if
	

end subroutine alat

subroutine vecsize(vec,vsize)

	real,dimension(3) :: vec
	real :: vsize

	vsize= sqrt((vec(1)**2)+(vec(2)**2)+(vec(3)**2))


end subroutine vecsize


subroutine modvec(k,kp,modk) !subrotina que calcula o módulo da diferença de dois vetores

	implicit none

	real, dimension(3) :: k,kp
	real:: modk


	modk=sqrt((k(1)-kp(1))**2 + (k(2)-kp(2))**2 + (k(3)-kp(3))**2 )


end subroutine modvec

subroutine modvecq(q,modq) !subrotina que calcula o módulo do vetor q

	implicit none

	real, dimension(4) :: q
	real:: modq


	modq=sqrt((q(2))**2 + (q(3))**2 + (q(4))**2)


end subroutine modvecq

subroutine nelec4(nedos,efermi,en,tdos,neletrons) !simpson 3/8 - integral dos

	implicit none

	integer :: nedos,ifermi,i
	real,dimension(nedos) :: en,tdos

	real :: neletrons,hv
	real :: dx

	real :: flag1,efermi


	dx=abs(en(2)-en(1))

	neletrons=0
	
	do i=1,nedos/3
	
	flag1= (tdos(3*i-2)*hv(en(3*i-2),efermi))+3.0*(tdos(3*i-1)*hv(en(3*i-1),efermi))&
              +3.0*(tdos(3*i)*hv(en(3*i),efermi))+(tdos(3*i+1)*hv(en(3*i+1),efermi))
		
		neletrons=neletrons+flag1
	
		!write(*,*) "progresso:",i,"/",npassos

	end do
	
		neletrons=3.0*(dx/8.)*neletrons


end subroutine nelec4

subroutine gauss(nedos,en,edft,sme,gaussout) !smearing gaussiano

	implicit none
	integer :: nedos,i
	real :: en,sme
	real,dimension(nedos) :: edft,gaussout

	real,parameter :: pi=acos(-1.)
	real :: temp1,rnorm,deno

	rnorm=1.0/(sme*sqrt(2.0*pi))
	deno = 2.0*sme*sme

	do i=1,nedos
		temp1= -1.0*(edft(i)-en)*(edft(i)-en)/deno
		gaussout(i) = rnorm*exp(temp1)

	end do




end subroutine gauss


function hv(e,e0)

	implicit none

	real :: hv,e,e0

	if (e .lt. e0) then

		hv=1.0

	else

		hv=0.0
	end if


end function hv

!algoritmo buble sort modificado 
!a=coluna do array referencia para ordenação
!b=numero de colunas do array 
!n=numero de linhas do array

!return p,q in ascending order
Subroutine Orderm(a,b,p,q)

	implicit none

	integer :: a,b
	real,dimension(b) :: p,q,temp

  	if (p(a)>q(a)) then

    		temp=p

    		p=q

    		q=temp
  	end if
  	return

end subroutine Orderm

!Buuble sorting of integer array A
Subroutine Bubblem(a,b,vec, n)

	implicit none

	integer :: n,i,j,a,b
	real,dimension(n,b) :: vec


  	do i=1, n
    		do j=n, i+1, -1
      		call Orderm(a,b,vec(j-1,:), vec(j,:))
    		end do
  	end do

  	return
end subroutine Bubblem

!subrotina pra interpolação
subroutine interp1D(ndim,b,a,vec,x,res)

	implicit none

	integer :: a,b
	integer :: ndim,i,iaux
	real,dimension(ndim,b) :: vec
	real :: x,res

	do i=1,ndim

		if (x .le. vec(i,1)) then

			iaux=i

			go to 105
		else
			continue
		end if

	end do

105 continue
        
	res= vec(iaux-1,a) + (x-vec(iaux-1,1))*((vec(iaux,a)-vec(iaux-1,a))/(vec(iaux,1)-vec(iaux-1,1)))

end subroutine interp1D

function gaussian(deltaen,sme)

	implicit none
	
	real :: gaussian
	real :: deltaen,sme
	real,parameter :: pi=acos(-1.)
	
	real :: norm,deno,aux1	
	
	norm = 1.0D0/(sme*sqrt(2.0D0*pi))
	deno = 2.0D0*sme*sme
	aux1 = -1.0D0*((deltaen)*(deltaen))/deno
	
	gaussian = norm*exp(aux1)	

end function gaussian

function lorentzian(deltaen,sme)

	implicit none
	
	real :: lorentzian
	real :: deltaen,sme
	real,parameter :: pi=acos(-1.)
	
	lorentzian = ((sme)/(pi*((deltaen))**2+sme**2))

end function lorentzian

subroutine quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

	integer :: ngkpt,nc,nv,w90basis
	integer,dimension(ngkpt) :: nocpk,nocpq
	integer,dimension(ngkpt*nc*nv,4) :: stt

	integer :: i,c,v

	integer :: counter, ncaux,nvaux


	counter=1

	do i=1,ngkpt

	ncaux=(nocpq(i)+1)+nc-1

	nvaux=nocpk(i)-nv+1


				if (nv .gt. nocpk(i)) then

					write(*,*) "Number of valence states higher then ocuppied states"
					STOP

				else if (nc .gt. (w90basis-nocpq(i))) then

					write(*,*) "Number of conduction states higher then unocuppied states"
					STOP

				else 

					continue
					
				end if


		do c=nv+1,nv+nc


			do v=1,nv



				stt(counter,1) = counter !numero do estado

				stt(counter,2) = v   !numero da banda da valencia

				stt(counter,3) = c   !numero da banda da condução

				stt(counter,4) = i   !numero do ponto k no grid

				counter=counter+1
 			

			end do



		end do



	end do


end subroutine quantumnumbers2

function gapcortemp(sparam,avgphonon,temp) !doi: 10.1063/1.104723

	implicit none

	real :: gapcortemp
	real :: sparam, avgphonon, temp
	real,parameter :: kb =8.617330350E-5
	real :: aux1,aux2
	
	
	aux1 = avgphonon/(2.0*kb*temp)
	aux2 = 1/tanh(aux1)
	
	
	gapcortemp = - sparam*avgphonon*((aux2)-1d0)


end function gapcortemp

function gapcortemp2(sparam,avgphonon,temp) !(bose-einstein model)doi: https://doi.org/10.1038/s41598-020-71808-y

	implicit none

	real :: gapcortemp2
	real :: sparam, avgphonon, temp
	real,parameter :: kb =8.617330350E-5
	real :: aux1,aux2
	
	
	aux1 = 2.0*sparam
	aux2 = exp(avgphonon/temp)
	
	
	gapcortemp2 = - (aux1/(aux2-1.0))


end function gapcortemp2

function fermidisteh(ec,ev,temp) result(difference)
    implicit none
    real, intent(in) :: ec,ev,temp
    real :: difference
    real, parameter :: kb=8.617330350E-5
    ! Preserve the existing energy reference: the chemical potential is zero.
    difference = occupation(ev)-occupation(ec)
contains
    real function occupation(energy)
        real, intent(in) :: energy
        real :: decay
        if (temp <= 0.0) then
            occupation = 0.5
            if (energy < 0.0) occupation = 1.0
            if (energy > 0.0) occupation = 0.0
        else
            ! exp always has a nonpositive argument, including at very low T.
            decay = exp(-abs(energy)/(kb*temp))
            if (energy >= 0.0) then
                occupation = decay/(1.0+decay)
            else
                occupation = 1.0/(1.0+decay)
            end if
        end if
    end function occupation
end function fermidisteh


