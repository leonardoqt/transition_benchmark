program test
	!
	use benchmark_system,  only: assign_model, assign_psi, rho0, Uini, nT, &
	                             evo_npi, evo_hst, evo_loc01, evo_rho_diab, &
	                             final_rho_hop, final_rho_hop_loc19, final_psi_hop_loc01, &
	                             test_U, test_T, test_E!, ZY_correct_sign_full, nstate
	use model_H,           only: sz1, sz2, sz3, szn, H1, H2, H3, Hn, H_sys_rotate, H_sys_parallel, H_sys_mix
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp)   , allocatable :: psi(:), Tv(:,:,:), &
	                            hop_p_npi(:), hop_p_hst(:), hop_p_loc01_1(:), hop_p_loc01_2(:), hop_p_loc19(:), &
	                            pop_p_npi(:,:), pop_p_hst(:,:), pop_p_loc01_1(:,:), pop_p_loc01_2(:,:), pop_p_loc19(:,:)
	complex(dp), allocatable :: Ut(:,:,:), rho_f_npi(:,:), rho_f_hst(:,:), rho_f_loc(:,:), psi_f_loc(:,:)
	!
	real(dp)   , allocatable :: rho_diag_diab(:), rho_diag_npi(:), rho_diag_hst(:), rho_diag_loc(:)
	real(dp)   , allocatable :: aux_vec(:,:)
	real(dp)                 :: dt, shift
	integer                  :: idt, imodel, sz, istate, ini_type
	character(len=100)       :: f_sz
	!
	!!!!!!!!!!!!!!!!!!!!!
	!real(dp), allocatable :: ss(:,:), uu(:,:)
	!allocate(ss(4,4))
	!allocate(uu(4,4))
	!nstate = 4
	!ss = 0.d0
	!uu = 0.d0
	!ss(1,2) = 1.d0
	!ss(2,1) = 1.d0
	!ss(3,4) = 1.d0
	!ss(4,3) = 1.d0
	!uu(1,1) = 1.d0
	!uu(2,2) = 1.d0
	!uu(3,3) = 1.d0
	!uu(4,4) = 1.d0
	!call ZY_correct_sign_full(ss,uu)
	!write(*,'(4(ES14.6,1X))') transpose(uu)
	!stop
	!!!!!!!!!!!!!!!!!!!!!
	!
	read(*,*) imodel, ini_type, sz, dt, shift
	!
	select case (imodel)
		case(1)
			call assign_model(1.0d0,-1.5d0,1.5d0,dt,H_sys_rotate,sz,shift)
		case(2)
			call assign_model(1.0d0,-1.5d0,1.5d0,dt,H_sys_parallel,sz,shift)
		case(3)
			call assign_model(1.0d0,-1.5d0,1.5d0,dt,H_sys_mix,sz,shift)
	end select
	!
	allocate(psi(sz))
	allocate(aux_vec(sz,1))
	psi = 0.d0
	aux_vec = 0.d0
	select case (ini_type)
		case(0)
			! ground diabats
			psi(1) = 1.d0
		case(1)
			! ground adiabats
			aux_vec(1,1) = 1.d0
			aux_vec = matmul(Uini, aux_vec)
			psi = aux_vec(:,1)
		case(2)
			! all adiabats
			do istate = 1, sz
				aux_vec(istate,1) = 1.d0/sqrt(sz*1.d0)
			enddo
			aux_vec = matmul(Uini, aux_vec)
			psi = aux_vec(:,1)
	end select
	call assign_psi(psi)
	call evo_rho_diab(rho_diag_diab)
	!
	call evo_npi(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_npi, hop_p_npi, pop_p_npi, 1)
	!!call test_U()
	call evo_hst(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_hst, hop_p_hst, pop_p_hst, 1)
	call evo_loc01(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_loc, hop_p_loc01_1, pop_p_loc01_1, 1)
	call final_psi_hop_loc01(Ut, Tv, psi_f_loc, hop_p_loc01_2, pop_p_loc01_2, 1)
	call final_rho_hop_loc19(Ut, Tv, rho_f_loc, hop_p_loc19, pop_p_loc19, 1)
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
	!
	write(*,*) 'Final population from rho'
	write(*,'("DIAB   " '//adjustl(f_sz)//'(ES14.6,1X))') rho_diag_diab
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
	!call test_E
end program test
