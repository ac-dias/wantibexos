INCLUDE "./subroutines/sistema_tb.f90"
INCLUDE "./subroutines/module_input_read.f90"
INCLUDE "./subroutines/mhkpack_subs.f90"
INCLUDE "./subroutines/general_subs.f90"
INCLUDE "./subroutines/optics.f90"
INCLUDE "./subroutines/bse_subs.f90"
INCLUDE "./subroutines/corrections_sys.f90"
INCLUDE "./subroutines/coulomb_pot.f90"


!gfortran -mcmodel=large bse_optics-tool.f90 -o bse_opt-tool-teste.x -llapack95 -lopenblas -fopenmp

program main

	use omp_lib
	use input_variables
	use hamiltonian_input_variables

	implicit none

	double precision,parameter:: pi=acos(-1.)

	integer :: dimbse !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector

	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:) :: stt

	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	!double precision,allocatable,dimension(:,:) :: lcount !pertencimento a camada de um dado estado

	!double precision,allocatable,dimension(:) :: orbweight

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	double complex:: matrizelbse

	double complex :: hxsp,hysp,hzsp

	double complex,allocatable,dimension(:) :: hoptxf,hoptyf,hoptzf,hoptspf,hoptsmf 

	double precision,allocatable,dimension(:) :: activitylinearx,activitylineary,activitylinearz,activitysigmap,activitysigmam
	double precision,allocatable :: descriptionlinearx(:),descriptionlineary(:),descriptionlinearz(:)
	double precision,allocatable :: descriptionsigmap(:),descriptionsigmam(:)

	double precision,allocatable,dimension(:) :: pinter,pintra

	double complex,parameter :: imag=cmplx(0.0,1.0)

	double precision :: a,r0,ed,egap

	integer :: ngkpt
	!double precision,dimension(3) :: shift !shift no mhkpack

	!variaveis relacionadas a marcacao do tempo

	double precision:: t0,tf
	integer,dimension(8) :: values,values2

	integer,allocatable,dimension(:) :: nocpk

	!definicoes diagonalizacao 
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
	COMPLEX*16,allocatable,dimension(:,:) :: hbse

        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        double complex,allocatable,dimension (:) :: WORK

	call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada orb weight"

	!OUTPUT
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida log BSE"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bse_opt_x.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida opt linear x"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"bse_opt_y.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida opt linear y"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"bse_opt_z.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida opt linear z"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"bse_opt_sp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida opt sigma +"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"bse_opt_sm.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida opt sigma -"
	!OPEN(UNIT=306, FILE=trim(outputfolder)//"exc_layer.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida opt sigma -"


	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)
	!call OPENBLAS_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)


	!termino leitura parametros

	!parametros do calculo


	OPEN(UNIT=2077, FILE= trim(calcparms)//"bse_opt_pol_calc.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype)

	close(2077)



	!termino parametros calculo 

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if


	!Informações para o arquivo de log do calculo

	write(300,*)
	write(300,*)
	write(300,*)
	write(300,*) 'numero de threads:', nthreads
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,*) 'ktol-coulomb:', ktol
	write(300,*)
	write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	!write(300,*) 'angulo theta:',beta,'rad'

	write(300,*)
	write(300,*) 'bandas de condução, nc:',nc,'  ','bandas de valencia, nv:',nv
	write(300,*)
	!write(300,*) 'constante dielétrica:',ed,'  ','r0:',r0

	!write(300,*)
	!write(300,*) 'intensidade termo exchange:',exc
	!write(300,*)
	!write(300,*) 'intensidade termo exchange l:',excl
	!write(300,*)
	!write(300,*) 'intensidade campo eletrico:',ez
	!write(300,*)
	!write(300,*) 'direção da magnetização do termo de exchange:',mag(1),'x',mag(2),'y',mag(3),'z'
	!write(300,*)

	!write(300,*) 'arquivo de saida 1:',logbse
	!write(300,*) 'arquivo de saida 2:',bseoptx
	!write(300,*) 'arquivo de saida 3:',bseopty
	!write(300,*) 'arquivo de saida 4:',bseoptsp
	!write(300,*) 'arquivo de saida 5:',bseoptsm

	
	write(300,*)
	write(300,*) 'inicio','  ','dia',values(3),'',values(5),'horas',values(6),'min',values(7),'seg'
	write(300,*) 

	call flush(300)

	allocate(kpt(ngkpt,3))

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	!alocacao de memoria

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,w90basis),vector(ngkpt,w90basis,w90basis))
	allocate(nocpk(ngkpt))

	allocate (stt(ngkpt*nc*nv,4))
	allocate(hoptxf(dimbse),hoptyf(dimbse),hoptzf(dimbse),hoptspf(dimbse),hoptsmf(dimbse))
	allocate(hbse(dimbse,dimbse),W(dimbse))

	allocate (RWORK(3*dimbse-2))
        LWMAX = 2*dimbse-1
	allocate(WORK(LWMAX))

	allocate(activitylinearx(dimbse),activitylineary(dimbse),activitysigmap(dimbse),activitysigmam(dimbse))

	allocate(descriptionlinearx(dimbse),descriptionlineary(dimbse),descriptionsigmap(dimbse),descriptionsigmam(dimbse))

	allocate(activitylinearz(dimbse),descriptionlinearz(dimbse))

	!allocate(orbweight(w90basis))
	
	!do i=1,w90basis
	!	read(203,*) orbweight(i)
	!end do

	!allocate(lcount(ngkpt,w90basis))
	egap = 50.0
	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


			do j=1,w90basis
	
				eigv(i,j)= eaux(j)

			end do
			
			!if (ntype .eq. 2) then

			!do kl = 1,w90basis

			!	call layercont(w90basis,vaux(kl,:),orbweight,lcount(i,kl))


			!end do

			!else 

			!	continue

			!end if

			if ((eigv(i,nocpk(i)+1)-eigv(i,nocpk(i))) .le. egap) then

				egap = eigv(i,nocpk(i)+1)-eigv(i,nocpk(i))

			else
				continue

			end if


			do l=1,w90basis


				do h=1,w90basis

					vector(i,l,h)=vaux(l,h)


				end do
			

			end do

	
	end do
	! $omp end parallel do

	deallocate(eaux,vaux)
	write(300,*) 'direct gap:', egap
	write(300,*) 'termino calculo eigsys sp'
	call flush(300)


	!definindo os numeros quanticos dos estados


	!allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)


	!allocate(hoptxf(dimbse),hoptyf(dimbse),hoptzf(dimbse),hoptspf(dimbse),hoptsmf(dimbse))

	!$omp parallel do

	do i=1,dimbse

		call optsp(vector(stt(i,4),stt(i,2),:),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
		     ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)


		hoptxf(i)= hxsp
		hoptyf(i)= hysp
		hoptzf(i)= hzsp
		hoptspf(i)= (hxsp+cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		hoptsmf(i)= (hxsp-cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))


	end do

	!$omp end parallel do

	write(300,*) 'termino optica sp'
	call flush(300)

	!go to 789

	!allocate(hbse(dimbse,dimbse),W(dimbse))

	hbse=0.

	!$omp parallel do 
!collapse(2)

	do i=1,dimbse



		do j=i,dimbse



  hbse(i,j)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ngrid,rlat,stt(i,:),eigv(stt(i,4)&
  	    ,stt(i,3)),eigv(stt(i,4),stt(i,2)),vector(stt(i,4)&
            ,stt(i,3),:) ,vector(stt(i,4),stt(i,2),:),kpt(stt(i,4),:),stt(j,:),eigv(stt(j,4),stt(j,3))&
  	    ,eigv(stt(j,4),stt(j,2)) &
            ,vector(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:))



		end do



	end do

	!$omp end parallel do

	write(300,*) 'termino hamiltoniano exciton'
	call flush(300)	



      	!LWORK = 2*dimbse+dimbse**2
      	!LIWORK = 3 + 5*dimbse
      	!LRWORK = 1 + 5*dimbse + 2*dimbse**2


!	allocate(WORK(LWORK))
!	allocate (IWORK(LIWORK))
!	allocate (RWORK(LRWORK))

 !     	CALL ZHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 
!	      RWORK, LRWORK, IWORK, LIWORK,INFO )
 !     	IF( INFO.GT. 0 ) THEN
  !      WRITE(*,*)'The algorithm failed to compute eigenvalues.'
  !       STOP
  !    	END IF   

	!allocate (RWORK(3*dimbse-2))
        !LWMAX = 2*dimbse-1
	!allocate(WORK(LWMAX))

     	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF   

	!allocate(pinter(dimbse),pintra(dimbse))

	!if (ntype .eq. 2) then

	!do i=1,dimbse

	!	call excitonil(i,w90basis,nc,nv,ngrid*ngrid,stt,dimbse,hbse(:,i),lcount,pinter(i),pintra(i))

	!	write(306,*) W(i),pinter(i),pintra(i)
	!	call flush(306)

	!end do 

	!else 

	!	continue

	!end if 

	write(300,*) 'termino diagonalizacao excitons'
	call flush(300)
	write(300,*) "exciton ground state",W(1)

	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices)
	deallocate(RWORK,WORK)
	!deallocate (IWORK)

	!allocate(activitylinearx(dimbse),activitylineary(dimbse),activitysigmap(dimbse),activitysigmam(dimbse))

	!allocate(descriptionlinearx(dimbse),descriptionlineary(dimbse),descriptionsigmap(dimbse),descriptionsigmam(dimbse))

	!allocate(activitylinearz(dimbse),descriptionlinearz(dimbse))

	call opticalactivity(dimbse,hbse,hoptxf,activitylinearx,descriptionlinearx)
	call opticalactivity(dimbse,hbse,hoptyf,activitylineary,descriptionlineary)
	call opticalactivity(dimbse,hbse,hoptzf,activitylinearz,descriptionlinearz)
	call opticalactivity(dimbse,hbse,hoptspf,activitysigmap,descriptionsigmap)
	call opticalactivity(dimbse,hbse,hoptsmf,activitysigmam,descriptionsigmam)

	!$omp do ordered
	do i=1,dimbse
		!$omp ordered
		write(301,*) W(i),activitylinearx(i)/W(i),descriptionlinearx(i)
		call flush(301)

		write(302,*) W(i),activitylineary(i)/W(i),descriptionlineary(i)
		call flush(302)

		write(303,*) W(i),activitylinearz(i)/W(i),descriptionlinearz(i)
		call flush(303)

		write(304,*) W(i),activitysigmap(i)/W(i),descriptionsigmap(i)
		call flush(304)

		write(305,*) W(i),activitysigmam(i)/W(i),descriptionsigmam(i)
		call flush(305)
		!$omp end ordered
	end do
	!$omp end do
	


	deallocate(activitylinearx,activitylineary,activitysigmap,activitysigmam)
	deallocate(descriptionlinearx,descriptionlineary,descriptionsigmap,descriptionsigmam)
	deallocate(activitylinearz,descriptionlinearz)

	deallocate(hoptxf,hoptyf,hoptzf,hoptspf,hoptsmf)
	!deallocate(lcount,pinter,pintra)

	deallocate(hbse,W,stt,nocpk)

	deallocate(kpt)

789     continue

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'fim','   ','dia',values2(3),'',values2(5),'horas',values2(6),'min',values2(7),'seg'
	write(300,*)




	close(200)
	close(203)



	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
	


end program main
