! Analytical dielectric-screening models.
module dielectric_models

	use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan,ieee_value
	implicit none
	private

	public :: trolle_pedersen_veniard_dielectric
	public :: tpv_bulk_dielectric
	public :: tpv_effective_2d_dielectric

contains

! Bulk-like non-local dielectric function from Eq. (1) of
! M. L. Trolle, T. G. Pedersen, and V. Veniard,
! Sci. Rep. 7, 39844 (2017), DOI: 10.1038/srep39844.
!
! q and q_tf are in inverse Angstrom, while hbar_omega_p is in eV.  The
! free-electron kinetic prefactor hbar**2/(2*m_e) is therefore in eV Angstrom**2.
! The optional Cappellini fitting parameter alpha defaults to the value 1.5
! used in the paper.
pure real function tpv_bulk_dielectric(q,kappa,q_tf,hbar_omega_p,alpha) result(epsilon_q)

	implicit none

	real,intent(in) :: q,kappa,q_tf,hbar_omega_p
	real,intent(in),optional :: alpha

	real(kind=8),parameter :: hbar2_over_2me=3.809982116154859d0
	real(kind=8) :: alpha_value,denominator,q_abs

	alpha_value=1.5d0
	if (present(alpha)) alpha_value=dble(alpha)

	if (kappa .lt. 1.0 .or. q_tf .le. 0.0 .or. hbar_omega_p .le. 0.0 .or. &
	    alpha_value .lt. 0.0d0) then
		epsilon_q=ieee_value(1.0,ieee_quiet_nan)
		return
	end if

	! Vacuum is an exact and useful special case; it also avoids 1/(kappa-1).
	if (abs(dble(kappa)-1.0d0) .le. epsilon(1.0d0)) then
		epsilon_q=1.0
		return
	end if

	q_abs=abs(dble(q))
	denominator=1.0d0/(dble(kappa)-1.0d0) &
	            +alpha_value*(q_abs/dble(q_tf))**2 &
	            +(hbar2_over_2me*q_abs*q_abs/dble(hbar_omega_p))**2
	epsilon_q=real(1.0d0+1.0d0/denominator)

end function tpv_bulk_dielectric


! Effective finite-thickness 2D dielectric function from Eq. (12) of the
! reference above.  epsilon_layer_q, epsilon_top_q, and epsilon_bottom_q are
! the dielectric functions of the slab and its two surrounding half-spaces at
! the same in-plane momentum q.  q is in inverse Angstrom and thickness is in
! Angstrom.
pure real function tpv_effective_2d_dielectric(q,thickness,epsilon_layer_q, &
	                                           epsilon_top_q,epsilon_bottom_q) &
	                                           result(epsilon_2d_q)

	implicit none

	real,intent(in) :: q,thickness,epsilon_layer_q
	real,intent(in) :: epsilon_top_q,epsilon_bottom_q

	real(kind=8) :: beta,delta,denominator,difference,environment_product
	real(kind=8) :: environment_sum,layer,one_minus_exp,scaled_factor

	if (thickness .lt. 0.0 .or. epsilon_layer_q .le. 0.0 .or. &
	    epsilon_top_q .le. 0.0 .or. epsilon_bottom_q .le. 0.0) then
		epsilon_2d_q=ieee_value(1.0,ieee_quiet_nan)
		return
	end if

	beta=abs(dble(q))*dble(thickness)
	layer=dble(epsilon_layer_q)
	environment_sum=dble(epsilon_top_q)+dble(epsilon_bottom_q)
	environment_product=dble(epsilon_top_q)*dble(epsilon_bottom_q)

	! The beta=0 limit of Eq. (12) is the average environmental screening.
	if (beta .le. epsilon(1.0d0)) then
		epsilon_2d_q=real(0.5d0*environment_sum)
		return
	end if

	! Stable evaluations of 1-exp(-beta) and beta-1+exp(-beta).  Their direct
	! forms lose precision in the long-wavelength limit.
	if (beta .lt. 1.0d-3) then
		one_minus_exp=beta-beta**2/2.0d0+beta**3/6.0d0-beta**4/24.0d0 &
		              +beta**5/120.0d0-beta**6/720.0d0
		delta=beta**2/2.0d0-beta**3/6.0d0+beta**4/24.0d0 &
		      -beta**5/120.0d0+beta**6/720.0d0
	else
		one_minus_exp=1.0d0-exp(-beta)
		delta=beta-one_minus_exp
	end if

	! This is Eq. (12), multiplied by 2*exp(-beta) and regrouped so that
	! neither hyperbolic overflow nor small-beta cancellation is introduced.
	difference=layer*environment_sum-environment_product-layer*layer
	scaled_factor=2.0d0*layer*environment_sum &
	              -2.0d0*difference*one_minus_exp &
	              +difference*one_minus_exp*one_minus_exp
	denominator=2.0d0*layer*environment_sum*delta &
	            +difference*beta*(-2.0d0*one_minus_exp &
	                              +one_minus_exp*one_minus_exp) &
	            +(layer*environment_sum-2.0d0*environment_product) &
	             *one_minus_exp*one_minus_exp

	if (abs(denominator) .le. tiny(denominator)) then
		epsilon_2d_q=ieee_value(1.0,ieee_quiet_nan)
	else
		epsilon_2d_q=real(layer*scaled_factor*delta/denominator)
	end if

end function tpv_effective_2d_dielectric


! Complete Trolle-Pedersen-Veniard model: Eq. (1) supplies the slab response
! used by the finite-thickness effective screening in Eq. (12).  The two
! environmental dielectric values may themselves be evaluated with
! tpv_bulk_dielectric when their q-dependence is required.
pure real function trolle_pedersen_veniard_dielectric(q,thickness,kappa,q_tf, &
	                                                  hbar_omega_p,epsilon_top_q, &
	                                                  epsilon_bottom_q,alpha) &
	                                                  result(epsilon_2d_q)

	implicit none

	real,intent(in) :: q,thickness,kappa,q_tf,hbar_omega_p
	real,intent(in) :: epsilon_top_q,epsilon_bottom_q
	real,intent(in),optional :: alpha

	real :: alpha_value,epsilon_layer_q

	alpha_value=1.5
	if (present(alpha)) alpha_value=alpha
	epsilon_layer_q=tpv_bulk_dielectric(q,kappa,q_tf,hbar_omega_p,alpha_value)
	epsilon_2d_q=tpv_effective_2d_dielectric(q,thickness,epsilon_layer_q, &
	                                        epsilon_top_q,epsilon_bottom_q)

end function trolle_pedersen_veniard_dielectric

end module dielectric_models
