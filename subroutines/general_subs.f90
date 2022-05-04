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

	double complex,dimension(w90basis) :: lvec,rvec,clvec,auxvec
	double complex,dimension(w90basis,w90basis) :: hm

	double complex :: res

	call matvec(hm,rvec,w90basis,auxvec)
	call vecconjg(lvec,w90basis,clvec)

	call prodintsq(clvec,auxvec,w90basis,res)


end subroutine sandwich


subroutine prodintsq(vetora,vetorb,n,resultado)

	implicit none

	integer:: n,i
	double complex :: resultado1, auxi
	double complex, dimension(n) :: vetora,vetorb,vetorc
	double complex :: resultado

	

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
	double complex, dimension(n) :: vetor,vetout


	do i=1,n

		vetout(i)=conjg(vetor(i))

	end do



end subroutine vecconjg 

subroutine matvec(matriz,vetor,n,veout)

	implicit none

	integer :: n !ordem da matriz e tamanho do vetor
	double complex, dimension(n,n) :: matriz
	double complex, dimension(n):: vetor,veout
	double complex :: flag
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
	double complex h11inv(n1,n1)
        integer    lda1,lwk1,if1,ipiv1(n1)
        double complex wk1(64*n1)
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
	
	double precision,dimension(3) :: v1,v2,vx

	vx(1) = (v1(2)*v2(3))-(v1(3)*v2(2))
	vx(2) = (v1(3)*v2(1))-(v1(1)*v2(3))
	vx(3) = (v1(1)*v2(2))-(v1(2)*v2(1))


end subroutine prodvec

subroutine vcell2D(rlat,vc)

	double precision,dimension(3,3) :: rlat
	double precision :: vc

	double precision,dimension(3) :: vx
	double precision :: aux

		call prodvec(rlat(1,:),rlat(2,:),vx)
		call vecsize(vx,vc)

end subroutine vcell2D

subroutine vcell3D(rlat,vc)

	double precision,dimension(3,3) :: rlat
	double precision :: vc

	double precision,dimension(3) :: vx
	double precision :: aux

		call prodvec(rlat(1,:),rlat(2,:),vx)

		vc=vx(1)*rlat(3,1)+vx(2)*rlat(3,2)+vx(3)*rlat(3,3)
		vc = abs(vc)

end subroutine vcell3D

subroutine vcell(systype,rlat,vc)

	character(len=2) :: systype
	double precision,dimension(3,3) :: rlat
	double precision :: vc

	double precision,dimension(3) :: vx
	double precision :: aux

	if (systype .eq. "2D") then

		call prodvec(rlat(1,:),rlat(2,:),vx)
		call vecsize(vx,vc)
		
	else
		call prodvec(rlat(1,:),rlat(2,:),vx)

		vc=vx(1)*rlat(3,1)+vx(2)*rlat(3,2)+vx(3)*rlat(3,3)
		vc = abs(vc)

	end if


end subroutine vcell

subroutine alat2D(rlat,a0)

	double precision,dimension(3,3) :: rlat
	double precision :: a0
	double precision,dimension(3) :: vsize

	call vecsize(rlat(1,:),vsize(1))
	call vecsize(rlat(2,:),vsize(2))
	call vecsize(rlat(3,:),vsize(3))

		a0 = 0.5*(vsize(1)+vsize(2))

end subroutine alat2D

subroutine alat(systype,rlat,a0)

	character(len=2) :: systype
	double precision,dimension(3,3) :: rlat
	double precision :: a0

	double precision,dimension(3) :: vsize

	call vecsize(rlat(1,:),vsize(1))
	call vecsize(rlat(2,:),vsize(2))
	call vecsize(rlat(3,:),vsize(3))


	if (systype .eq. "2D") then

		a0 = 0.5*(vsize(1)+vsize(2))
	else
		a0 = (1.0/3.0)*(vsize(1)+vsize(2)+vsize(3))

	end if
	

end subroutine alat

subroutine vecsize(vec,vsize)

	double precision,dimension(3) :: vec
	double precision :: vsize

	vsize= sqrt((vec(1)**2)+(vec(2)**2)+(vec(3)**2))


end subroutine vecsize


subroutine modvec(k,kp,modk) !subrotina que calcula o módulo da diferença de dois vetores

	implicit none

	double precision, dimension(3) :: k,kp
	double precision:: modk


	modk=sqrt((k(1)-kp(1))**2 + (k(2)-kp(2))**2 + (k(3)-kp(3))**2 )


end subroutine modvec

subroutine modvecq(q,modq) !subrotina que calcula o módulo do vetor q

	implicit none

	double precision, dimension(4) :: q
	double precision:: modq


	modq=sqrt((q(2))**2 + (q(3))**2 + (q(4))**2)


end subroutine modvecq

subroutine nelec4(nedos,efermi,en,tdos,neletrons) !simpson 3/8 - integral dos

	implicit none

	integer :: nedos,ifermi,i
	double precision,dimension(nedos) :: en,tdos

	double precision :: neletrons,hv
	double precision :: dx

	double precision :: flag1,efermi


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
	double precision :: en,sme
	double precision,dimension(nedos) :: edft,gaussout

	double precision,parameter :: pi=acos(-1.)
	double precision :: temp1,rnorm,deno

	rnorm=1.0/(sme*dsqrt(2.0*pi))
	deno = 2.0*sme*sme

	do i=1,nedos
		temp1= -1.0*(edft(i)-en)*(edft(i)-en)/deno
		gaussout(i) = rnorm*dexp(temp1)

	end do




end subroutine gauss


function hv(e,e0)

	implicit none

	double precision :: hv,e,e0

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
	double precision,dimension(ndim,b) :: vec
	double precision :: x,res

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
	
	double precision :: gaussian
	double precision :: deltaen,sme
	double precision,parameter :: pi=acos(-1.)
	
	double precision :: norm,deno,aux1	
	
	norm = 1.0D0/(sme*DSQRT(2.0D0*pi))
	deno = 2.0D0*sme*sme
	aux1 = -1.0D0*((deltaen)*(deltaen))/deno
	
	gaussian = norm*dexp(aux1)	

end function gaussian

function lorentzian(deltaen,sme)

	implicit none
	
	double precision :: lorentzian
	double precision :: deltaen,sme
	double precision,parameter :: pi=acos(-1.)
	
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

	double precision :: gapcortemp
	double precision :: sparam, avgphonon, temp
	double precision,parameter :: kb =8.617330350E-5
	double precision :: aux1,aux2
	
	
	aux1 = avgphonon/(2.0*kb*temp)
	aux2 = 1/dtanh(aux1)
	
	
	gapcortemp = - sparam*avgphonon*((aux2)-1d0)


end function gapcortemp

function gapcortemp2(sparam,avgphonon,temp) !(bose-einstein model)doi: https://doi.org/10.1038/s41598-020-71808-y

	implicit none

	double precision :: gapcortemp2
	double precision :: sparam, avgphonon, temp
	double precision,parameter :: kb =8.617330350E-5
	double precision :: aux1,aux2
	
	
	aux1 = 2.0*sparam
	aux2 = dexp(avgphonon/temp)
	
	
	gapcortemp2 = - (aux1/(aux2-1.0))


end function gapcortemp2

function fermidisteh(ec,ev,temp) !distribuição fermi-dirac -> fv-fc

	implicit none
	
	double precision :: fermidisteh
	double precision :: ec,ev,temp
	double precision,parameter :: efermi = 0d0
	double precision,parameter :: kb =8.617330350E-5
	double precision :: beta
	double precision :: fv,fc
	
	if ( temp .eq. 0d0 ) then
	
		temp = 1E-8
	
	end if	

	beta = 1.0/(kb*temp)
	
	fv = 1.0/ (dexp((ev-efermi)*beta)+1.0)
	
	fc = 1.0/ (dexp((ec-efermi)*beta)+1.0)
	
	
	fermidisteh = fv-fc



end function fermidisteh






