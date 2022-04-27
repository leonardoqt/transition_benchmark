program test
	!
	use benchmark_system,  only: assign_model, assign_psi, rho0, evo_npi, evo_hst, evo_loc01, nT, &
	                             final_rho_hop, final_rho_hop_loc19
	use model_H,           only: sz1, sz2, H1, H2
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp)   , allocatable :: psi(:), Tv_npi(:,:,:), hop_p(:)
	complex(dp), allocatable :: Ut_npi(:,:,:), rho_f(:,:)
	integer                  :: idt
	!
	call assign_model(1.0d0,-0.5d0,2.5d0,1.0d-4,H1,sz1)
	call assign_model(1.0d0,-0.5d0,2.5d0,1.0d-4,H2,sz2)
	!write(*,'(3(ES12.4,1X))') U
	!
	allocate(psi(sz2))
	psi(1) = 1.d0
	call assign_psi(psi)
	!
	call evo_loc01(Ut_npi, Tv_npi)
	!call final_rho_hop(Ut_npi, Tv_npi, rho_f, hop_p, 1)
	call final_rho_hop_loc19(Ut_npi, Tv_npi, rho_f, hop_p, 1)
	!
	write(*,'(3(ES12.4,1X))') hop_p
	write(*,*) ''
	write(*,*) ''
	write(*,'(3(ES12.4,1X))')  dble(rho_f)
	write(*,*) ''
	write(*,'(3(ES12.4,1X))') aimag(rho_f)
end program test
