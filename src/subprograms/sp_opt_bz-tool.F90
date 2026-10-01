
!funciona apenas para semicondutores

!gfortran -mcmodel=large sp_opt_bz-tool.f90 -o sp_opt_bz-tool.x -llapack95 -lopenblas -fopenmp

subroutine spoptpolbz(nthreads,dft,outputfolder,ngrid,nc,nv, &
		     sme,params,exc,mshift,nocpf,fermishift,meshgen,mag)

	use omp_lib
	use hamiltonian_input_variables
	use input_variables, only: bsecenterfile
	use bse_q_optics, only: rmn_data, rmn_center_setup, rmn_destroy, rmn_dfts_dipoles, &
		rmn_apply_q0_optical_correction

	implicit none

	real,parameter:: pi=acos(-1.)

	integer :: dimrpa !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	real,allocatable,dimension(:,:) :: eigv
	real,allocatable,dimension(:,:) :: tau !DFT=S orbital centres
	real,allocatable,dimension(:,:,:,:) :: dmat !DFT=S orbital dipoles (BSE_CENTER_FILE)
	type(rmn_data) :: rmn
	logical :: use_dmat,use_rmn,rmn_ok
	character(len=256) :: rmn_message
	complex,allocatable,dimension(:,:,:) :: vector

	real,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:,:) :: stto

	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores


	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux,nocpf
	real :: fermishift

	complex :: hxsp,hysp,hzsp
 

	real :: auxx,auxy,auxz,auxsp,auxsm


	complex,parameter :: imag=cmplx(0.0,1.0)

	integer :: nk

	real,allocatable,dimension(:,:) :: output

	real,parameter :: norm=(1.0/sqrt(2.))

	integer :: ngkpt

	integer,allocatable,dimension(:) :: nocpk

	!modificacoes versao 2.1

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	real :: edos0,edosf,numdos
	real :: ebse0,ebsef,numbse
	real :: sme,exc,rk
	real,dimension(3) :: mshift
	real :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	!character(len=70) :: diein    !ambiente dieletrico
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=70) :: meshgen	
	character(len=5) :: coultype
	character(len=1) :: dft
	real,dimension(3) :: ediel,mag


	!fim modificacoes versao 2.1
	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"


	!OUTPUT
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bz_act_x.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_x output file"
	OPEN(UNIT=302, FILE= trim(outputfolder)//"bz_act_y.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_y output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"bz_act_z.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_z output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"bz_act_sp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_sp output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"bz_act_sm.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_sm output file"
	OPEN(UNIT=306, FILE=trim(outputfolder)//"dichroism_cp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening dichroism output file"
	OPEN(UNIT=307, FILE=trim(outputfolder)//"dichroism_xy.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening dichroism output file"    	




	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	call hamiltonian_input_read(200,params)

	allocate(ovp(nvec,w90basis,w90basis))
	ovp = 0.0

	end select
	
	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo
	dimrpa = nc*nv

	

	select case (meshgen)
	
	case ("2DHEX")
	
	 call gridgen2dhexcount(ngrid(1),ngrid(2),rlat(1,1),ngkpt)
	
	 allocate(kpt(ngkpt,3))
	
	 call gridgen2dhex(ngrid(1),ngrid(2),rlat(1,1),ngkpt,kpt)
	
	
	case ("2DRET")
	
	 ngkpt = ngrid(1)*ngrid(2)
	
	 allocate(kpt(ngkpt,3))
	
	 call gridgen2dret(ngrid(1),ngrid(2),rlat(1,1),rlat(2,2),kpt)
	
	case default

	 ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	
	 allocate(kpt(ngkpt,3))
	 
	 call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)
	
	end select



	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,w90basis),vector(ngkpt,w90basis,w90basis))
	allocate(nocpk(ngkpt))


	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux)
	do i=1,ngkpt

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,ovp,efermi,eaux,vaux,nocpf,fermishift,mag)
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif		

			do j=1,w90basis
	
				eigv(i,j)= eaux(j)

			end do
			


			do l=1,w90basis


				do h=1,w90basis

					vector(i,l,h)=vaux(l,h)


				end do
			

			end do
			
			

	
	end do
	!$omp end parallel do

	deallocate(eaux,vaux)



	!definindo os numeros quanticos dos estados




	allocate (stto(ngkpt,nc*nv,3))

do i=1,ngkpt

	counter=1
	
	ncaux=(nocpk(i)+1)+nc-1

	nvaux=nocpk(i)-nv+1

	do c=(nocpk(i)+1),ncaux


		do v=nvaux,nocpk(1)


			stto(i,counter,1) = counter !numero do estado

			stto(i,counter,2) = v   !numero da banda da valencia

			stto(i,counter,3) = c   !numero da banda da condução

			counter=counter+1
 			

		end do



	end do

end do

	counter=counter-1 !numero total de estados para equação bse

	!allocate(hoptxf(nk,dimrpa),hoptyf(nk,dimrpa),hoptspf(nk,dimrpa),hoptsmf(nk,dimrpa))

	allocate(output(ngkpt,7))
	
	output = 0.0

	! DFT=S: orbital centres for the atomic-gauge optical vertex
	allocate(tau(w90basis,3))
	tau = 0.0
	if (dft .eq. "S") call orbital_centres(w90basis,tau)
	! BSE_CENTER_FILE, the position matrix of the basis, as in bsesolver:
	! DFT=S: the dipoles between basis orbitals (siesta2wtb.py --rmatrix), for
	! optdip_s; DFT=W: the -i[A,H] term of the Wannier position matrix
	use_rmn = len_trim(bsecenterfile) > 0
	use_dmat = dft .eq. "S" .and. use_rmn
	if (use_rmn) then
		call rmn_center_setup(trim(bsecenterfile),dft,w90basis,rmn,rmn_ok,rmn_message)
		if (rmn_ok .and. use_dmat) then
			allocate(dmat(nvec,w90basis,w90basis,3))
			call rmn_dfts_dipoles(rmn,w90basis,nvec,rvec,rlat,ovp,tau,dmat,rmn_ok,rmn_message)
			call rmn_destroy(rmn)
		end if
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			if (use_dmat) stop 'Unable to build the DFT=S orbital dipoles'
			stop 'Unable to initialize the Wannier position data'
		end if
	end if

	! The optical vertex takes the DFT energies: the scissor scs, which eigsys
	! adds to the conduction bands, shifts the transition energies (BSE
	! diagonal, IPA spectra) but not the eigenvectors, so it must not change
	! r_cv. With the shifted energies <c|dH/dk|v>/(Ec-Ev) would shrink by
	! (Ec-Ev)/(Ec-Ev+scs), and for DFT=S the -(Ec+Ev)/2 dS/dk term would change.

	!$omp parallel do default(shared) private(j,hxsp,hysp,hzsp,rmn_ok,rmn_message)
	do j=1,ngkpt

	!auxx = 0.0
	!auxy = 0.0
	!auxz = 0.0
	!auxsp = 0.0
	!auxsm = 0.0

	! $omp parallel do reduction(+:auxx, auxy, auxz, auxsp, auxsm) private(i)
	do i=1,dimrpa

		if (dft .eq. "S") then
		call optspbz_s(eigv(j,stto(j,i,2)),vector(j,stto(j,i,2),:),&
		     eigv(j,stto(j,i,3))-scs,vector(j,stto(j,i,3),:),&
		     kpt(j,1),kpt(j,2),kpt(j,3),&
		     w90basis,nvec,rvec,hopmatrices,ihopmatrices,ovp,tau,&
		     hxsp,hysp,hzsp)
		else
		call optspbz(vector(j,stto(j,i,2),:),vector(j,stto(j,i,3),:),&
		     kpt(j,1),kpt(j,2),kpt(j,3),&
		     ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)
		end if
		     
		    

		hxsp= hxsp/(cmplx(eigv(j,stto(j,i,3))-scs-eigv(j,stto(j,i,2)),sme))
		hysp= hysp/(cmplx(eigv(j,stto(j,i,3))-scs-eigv(j,stto(j,i,2)),sme))
		hzsp= hzsp/(cmplx(eigv(j,stto(j,i,3))-scs-eigv(j,stto(j,i,2)),sme))
		if (use_dmat) call optdip_s(vector(j,stto(j,i,2),:),vector(j,stto(j,i,3),:),&
		     kpt(j,1),kpt(j,2),kpt(j,3),w90basis,nvec,rvec,dmat,hxsp,hysp,hzsp)
		if (use_rmn .and. .not. use_dmat) then
			call rmn_apply_q0_optical_correction(rmn,kpt(j,:),rlat,&
				eigv(j,stto(j,i,2)),vector(j,stto(j,i,2),:),&
				eigv(j,stto(j,i,3))-scs,vector(j,stto(j,i,3),:),sme,&
				hxsp,hysp,hzsp,rmn_ok,rmn_message)
			if (.not. rmn_ok) stop 'Unable to evaluate the Q=0 Wannier position matrix'
		end if
		
		output(j,1) = output(j,1)+real(hxsp*conjg(hxsp))	
		output(j,2) = output(j,2)+real(hysp*conjg(hysp))
		output(j,3) = output(j,3)+real(hzsp*conjg(hzsp))
		output(j,4) = output(j,4)&
		        +real(norm*(hxsp+cmplx(0.,1.)*hysp)*conjg(norm*(hxsp+cmplx(0.,1.)*hysp)))
	        output(j,5) = output(j,5)&
	               +real(norm*(hxsp-cmplx(0.,1.)*hysp)*conjg(norm*(hxsp-cmplx(0.,1.)*hysp)))		



	end do
	! $omp end parallel do

		!output(j,1) = auxx
		!output(j,2) = auxy
		!output(j,3) = auxz
		!output(j,4) = auxsp
		!output(j,5) = auxsm
		output(j,6) = (output(j,4)-output(j,5))/(output(j,4)+output(j,5))
		output(j,7) = (output(j,1)-output(j,2))/(output(j,1)+output(j,2))
		! $omp ordered
		
		!write(301,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxx
		!write(302,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxy
		!write(303,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxz
		!write(304,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxsp
		!write(305,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxsm
		!write(306,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),(auxsp-auxsm)/(auxsp+auxsm)
		!write(307,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),(auxx-auxy)/(auxx+auxy)
		
		!call flush(301)
		!call flush(302)
		!call flush(303)
		!call flush(304)
		!call flush(305)
		!call flush(306)
		!call flush(307)		
		! $omp end ordered
	end do
	!$omp end parallel do
	
	write(301,*) "#kx ky kz optx"
	write(302,*) "#kx ky kz opty"
	write(303,*) "#kx ky kz optz"
	write(304,*) "#kx ky kz optsp"
	write(305,*) "#kx ky kz optsm"
	write(306,*) "#kx ky kz circ dichroism-xy"
	write(307,*) "#kx ky kz linear dichroism-xy"


	do j=1,ngkpt
		write(301,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,1)
		write(302,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,2)
		write(303,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,3)
		write(304,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,4)
		write(305,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,5)
		write(306,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,6)
		write(307,"(3E20.8,1E20.6)") kpt(j,1),kpt(j,2),kpt(j,3),output(j,7)

	end do



	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(tau)
	if (use_dmat) deallocate(dmat)
	if (use_rmn .and. .not. use_dmat) call rmn_destroy(rmn)
	deallocate(ovp)
	!deallocate(hoptxf,hoptyf,hoptspf,hoptsmf)
	deallocate(stto)
	deallocate(kpt)
	deallocate(output)





	close(200)

	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
	close(306)
	close(307)

	


end subroutine spoptpolbz
