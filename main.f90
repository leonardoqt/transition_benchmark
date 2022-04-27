program test
	!
	use benchmark_system,  only: assign_model, assign_psi, rho0, evo_npi, evo_hst, evo_loc01, nT, &
	                             final_rho_hop, final_rho_hop_loc19
	use model_H,           only: sz1, sz2, H1, H2
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp)   , allocatable :: psi(:), Tv(:,:,:), hop_p_npi(:), hop_p_hst(:), hop_p_loc01(:), hop_p_loc19(:)
	complex(dp), allocatable :: Ut(:,:,:), rho_f(:,:)
	integer                  :: idt
	!
	allocate(psi(sz1))
	!
	call assign_model(1.0d0,-2.5d0,2.5d0,1.0d-4,H1,sz1)
	!call assign_model(1.0d0,-2.5d0,2.5d0,1.0d-4,H2,sz2)
	!
	psi(1) = 1.d0
	call assign_psi(psi)
	!
	call evo_npi(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f, hop_p_npi, 1)
	call evo_hst(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f, hop_p_hst, 1)
	call evo_loc01(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f, hop_p_loc01, 1)
	call final_rho_hop_loc19(Ut, Tv, rho_f, hop_p_loc19, 1)
	!
	write(*,*) 'Total hopping probability'
	write(*,'("NPI  " 2(ES12.4,1X))') hop_p_npi
	write(*,'("HST  " 2(ES12.4,1X))') hop_p_hst
	write(*,'("LOC01" 2(ES12.4,1X))') hop_p_loc01
	write(*,'("LOC19" 2(ES12.4,1X))') hop_p_loc19
	write(*,*) ''
	write(*,*) ''
	write(*,*) 'Final density matrix'
	write(*,'(2(ES12.4,1X))')  dble(rho_f)
	write(*,*) ''
	write(*,'(2(ES12.4,1X))') aimag(rho_f)
end program test
