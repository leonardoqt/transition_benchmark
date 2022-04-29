program test
	!
	use benchmark_system,  only: assign_model, assign_psi, rho0, evo_npi, evo_hst, evo_loc01, nT, &
	                             final_rho_hop, final_rho_hop_loc19, final_psi_hop_loc01
	use model_H,           only: sz1, sz2, H1, H2
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp)   , allocatable :: psi(:), Tv(:,:,:), hop_p_npi(:), hop_p_hst(:), hop_p_loc01_1(:), hop_p_loc01_2(:), hop_p_loc19(:)
	complex(dp), allocatable :: Ut(:,:,:), rho_f_npi(:,:), rho_f_hst(:,:), rho_f_loc(:,:), psi_f_loc(:,:)
	integer                  :: idt, state0 = 1
	!
	allocate(psi(sz2))
	!
	call assign_model(1.0d0,-2.5d0,2.5d0,1.0d-2,H1,sz1)
	call assign_model(1.0d0,-2.5d0,2.5d0,1.0d-3,H2,sz2)
	!
	psi = 0.d0
	psi(1) = 1.d0
	call assign_psi(psi)
	!
	call evo_npi(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_npi, hop_p_npi, state0)
	call evo_hst(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_hst, hop_p_hst, state0)
	call evo_loc01(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_loc, hop_p_loc01_1, state0)
	call final_psi_hop_loc01(Ut, Tv, psi_f_loc, hop_p_loc01_2, state0)
	call final_rho_hop_loc19(Ut, Tv, rho_f_loc, hop_p_loc19, state0)
	!
1001 format (3(ES14.6,1X))
	write(*,*) 'Total hopping probability'
	write(*,'("NPI    " 3(ES14.6,1X))') hop_p_npi
	write(*,'("HST    " 3(ES14.6,1X))') hop_p_hst
	write(*,'("LOC01-1" 3(ES14.6,1X))') hop_p_loc01_1
	write(*,'("LOC01-2" 3(ES14.6,1X))') hop_p_loc01_2
	write(*,'("LOC19  " 3(ES14.6,1X))') hop_p_loc19
	!write(*,*) ''
	!write(*,*) ''
	!write(*,*) 'Final density matrix (real)'
	!write(*,*) 'NPI'
	!write(*,1001)  dble(rho_f_npi)
	!write(*,*) 'HST'
	!write(*,1001)  dble(rho_f_hst)
	!write(*,*) 'LOC'
	!write(*,1001)  dble(rho_f_loc)
	!write(*,*) ''
	!write(*,*) 'Final density matrix (imag)'
	!write(*,*) 'NPI'
	!write(*,1001) aimag(rho_f_npi)
	!write(*,*) 'HST'
	!write(*,1001) aimag(rho_f_hst)
	!write(*,*) 'LOC'
	!write(*,1001) aimag(rho_f_loc)
	!write(*,*) ''
end program test
