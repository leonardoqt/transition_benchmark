program test
	!
	use benchmark_system,  only: assign_model, assign_psi, rho0, evo_npi, evo_hst, evo_loc01, nT, &
	                             final_rho_hop, final_rho_hop_loc19, final_psi_hop_loc01, &
	                             test_U, test_T
	use model_H,           only: sz1, sz2, sz3, szn, H1, H2, H3, Hn
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp)   , allocatable :: psi(:), Tv(:,:,:), &
	                            hop_p_npi(:), hop_p_hst(:), hop_p_loc01_1(:), hop_p_loc01_2(:), hop_p_loc19(:), &
	                            pop_p_npi(:,:), pop_p_hst(:,:), pop_p_loc01_1(:,:), pop_p_loc01_2(:,:), pop_p_loc19(:,:)
	complex(dp), allocatable :: Ut(:,:,:), rho_f_npi(:,:), rho_f_hst(:,:), rho_f_loc(:,:), psi_f_loc(:,:)
	!
	real(dp)   , allocatable :: rho_diag_npi(:), rho_diag_hst(:), rho_diag_loc(:)
	integer                  :: idt, state0 = 1, sz, istate
	character(len=100)       :: f_sz
	!
	!
	!sz = sz1
	!call assign_model(1.0d0,-2.5d0,2.5d0,1.0d-4,H1,sz)
	!sz = sz2
	!call assign_model(1.0d0,-2.5d-1,2.5d-1,1.0d-6,H2,sz)
	!sz = sz3
	!call assign_model(1.0d0,-2.5d0,2.5d0,1.0d-2,H3,sz)
	sz = szn
	call assign_model(1.0d0,-2.5d-1,2.5d-1,1.0d-5,Hn,sz)
	!
	allocate(psi(sz))
	psi = 0.d0
	psi(state0) = 1.d0
	call assign_psi(psi)
	!
	call evo_npi(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_npi, hop_p_npi, pop_p_npi, state0)
	!!call test_U()
	call evo_hst(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_hst, hop_p_hst, pop_p_hst, state0)
	call evo_loc01(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_loc, hop_p_loc01_1, pop_p_loc01_1, state0)
	call final_psi_hop_loc01(Ut, Tv, psi_f_loc, hop_p_loc01_2, pop_p_loc01_2, state0)
	call final_rho_hop_loc19(Ut, Tv, rho_f_loc, hop_p_loc19, pop_p_loc19, state0)
	!
	! use final rho of NPI as reference
	allocate(rho_diag_npi(sz))
	allocate(rho_diag_hst(sz))
	allocate(rho_diag_loc(sz))
	do istate = 1, sz
		rho_diag_npi(istate) = dble(rho_f_npi(istate,istate))
		rho_diag_hst(istate) = dble(rho_f_hst(istate,istate))
		rho_diag_loc(istate) = dble(rho_f_loc(istate,istate))
	enddo
	!
	write(f_sz,*) sz
	!write(*,*) 'Final density matrix (real)'
	!write(*,*) 'NPI'
	!write(*,'('//adjustl(f_sz)//'(ES14.6,1X))')  dble(rho_f_npi)
	!write(*,*) 'HST'
	!write(*,'('//adjustl(f_sz)//'(ES14.6,1X))')  dble(rho_f_hst)
	!write(*,*) 'LOC'
	!write(*,'('//adjustl(f_sz)//'(ES14.6,1X))')  dble(rho_f_loc)
	!write(*,*) ''
	!write(*,*) 'Final density matrix (imag)'
	!write(*,*) 'NPI'
	!write(*,'('//adjustl(f_sz)//'(ES14.6,1X))') aimag(rho_f_npi)
	!write(*,*) 'HST'
	!write(*,'('//adjustl(f_sz)//'(ES14.6,1X))') aimag(rho_f_hst)
	!write(*,*) 'LOC'
	!write(*,'('//adjustl(f_sz)//'(ES14.6,1X))') aimag(rho_f_loc)
	!write(*,*) ''
	!
	!write(*,*) 'Total hopping probability'
	!write(*,'("NPI    " '//adjustl(f_sz)//'(ES14.6,1X))') hop_p_npi
	!write(*,'("HST    " '//adjustl(f_sz)//'(ES14.6,1X))') hop_p_hst
	!write(*,'("LOC01-1" '//adjustl(f_sz)//'(ES14.6,1X))') hop_p_loc01_1
	!write(*,'("LOC01-2" '//adjustl(f_sz)//'(ES14.6,1X))') hop_p_loc01_2
	!write(*,'("LOC19  " '//adjustl(f_sz)//'(ES14.6,1X))') hop_p_loc19
	!write(*,*) ''
	!
	write(*,*) 'Final population from rho'
	write(*,'("NPI    " '//adjustl(f_sz)//'(ES14.6,1X))') rho_diag_npi
	write(*,'("HST    " '//adjustl(f_sz)//'(ES14.6,1X))') rho_diag_hst
	write(*,'("LOC    " '//adjustl(f_sz)//'(ES14.6,1X))') rho_diag_loc
	write(*,*) ''
	!
	write(*,*) 'Population from actual hopping'
	write(*,'("NPI    " '//adjustl(f_sz)//'(ES14.6,1X))') pop_p_npi
	write(*,'("HST    " '//adjustl(f_sz)//'(ES14.6,1X))') pop_p_hst
	write(*,'("LOC01-1" '//adjustl(f_sz)//'(ES14.6,1X))') pop_p_loc01_1
	write(*,'("LOC01-2" '//adjustl(f_sz)//'(ES14.6,1X))') pop_p_loc01_2
	write(*,'("LOC19  " '//adjustl(f_sz)//'(ES14.6,1X))') pop_p_loc19
end program test
