program test
	!
	use benchmark_system,  only: assign_model, assign_psi, Uini, nT, print_rho_Tvt, &
	                             evo_npi_interp, evo_npi, evo_hst, evo_loc01, evo_rho_diab, &
	                             final_rho_hop, final_rho_hop_loc19, final_psi_hop_loc01, final_psi_hop_loc01_dt, &
	                             final_rho_hop_interp, final_psi_hop_interp_dt, final_rho_hop_conditional_interp, &
	                             test_U, test_T, test_E, test_H!, ZY_correct_sign_full, nstate
	use model_H,           only: H_sys_single, H_sys_rotate, H_sys_parallel, H_sys_mix
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp)   , allocatable :: Tv(:,:,:), &
	                            pop_p_npi_dq(:,:), pop_p_npi_cdq(:,:), pop_p_npi_t_dq(:,:), pop_p_npi(:,:), pop_p_npi_t(:,:), &
	                            pop_p_hst(:,:), pop_p_loc01(:,:), pop_p_loc01_l(:,:), pop_p_loc01_t(:,:), pop_p_loc19(:,:)
	!
	complex(dp), allocatable :: psi(:), Ut(:,:,:), rho_f_npi_dq(:,:), rho_f_npi(:,:), rho_f_hst(:,:), &
	                            rho_f_loc(:,:), psi_f(:,:)
	!
	real(dp)   , allocatable :: rho_diag_diab(:), rho_diag_npi_dq(:), rho_diag_npi(:), rho_diag_hst(:), rho_diag_loc(:)
	real(dp)   , allocatable :: aux_vec(:,:)
	real(dp)                 :: dt, shift, x0, x1, threshold
	integer                  :: idt, imodel, sz, istate, ini_type, nqT, num_extra_call
	character(len=100)       :: f_sz
	!
	!
	read(*,*) imodel, ini_type, sz, dt, shift, nqT, threshold
	!
	x0 = -0.5d0
	x1 =  0.5d0
	select case (imodel)
		case(0)
			call assign_model(1.0d0,x0,x1,dt,H_sys_single,sz,shift)
		case(1)
			call assign_model(1.0d0,x0,x1,dt,H_sys_rotate,sz,shift)
		case(2)
			call assign_model(1.0d0,x0,x1,dt,H_sys_parallel,sz,shift)
		case(3)
			call assign_model(1.0d0,x0,x1,dt,H_sys_mix,sz,shift)
	end select
	!
	allocate(psi(sz))
	allocate(aux_vec(sz,1))
	psi = (0.d0, 0.d0)
	aux_vec = 0.d0
	select case (ini_type)
		case(0)
			! ground diabats
			psi(1) = (1.d0,0.d0)
		case(1)
			! ground adiabats
			aux_vec(1,1) = 1.d0
			aux_vec = matmul(Uini, aux_vec)
			psi = cmplx(aux_vec(:,1), aux_vec(:,1)*0.d0, dp)
		case(2)
			! all adiabats
			do istate = 1, sz
				aux_vec(istate,1) = 1.d0/sqrt(sz*1.d0)
			enddo
			aux_vec = matmul(Uini, aux_vec)
			psi = cmplx(aux_vec(:,1), aux_vec(:,1)*0.d0, dp)
		case(3)
			! random complex
			do istate = 1, sz
				psi(istate) = cmplx(rand(0),rand(0),dp)
			enddo
			aux_vec(1,1) = dble(dot_product(psi,psi))
			psi = psi / aux_vec(1,1)
	end select
	call assign_psi(psi)
	call print_rho_Tvt()
	stop
	!
	num_extra_call = 0
	call evo_npi_interp(Ut, Tv, nqT)
	call final_rho_hop_interp(Ut, Tv, rho_f_npi_dq, pop_p_npi_dq, nqT)
	call final_rho_hop_conditional_interp(rho_f_npi_dq, pop_p_npi_cdq, threshold, num_extra_call)
	call final_psi_hop_interp_dt(Ut, Tv, psi_f, pop_p_npi_t_dq, nqT)
	!
	!
	!-----------------------------------
	if ( nqT > 1 ) then
		select case (imodel)
			case(0)
				call assign_model(1.0d0,x0,x1,dt/nqT,H_sys_single,sz,shift)
			case(1)
				call assign_model(1.0d0,x0,x1,dt/nqT,H_sys_rotate,sz,shift)
			case(2)
				call assign_model(1.0d0,x0,x1,dt/nqT,H_sys_parallel,sz,shift)
			case(3)
				call assign_model(1.0d0,x0,x1,dt/nqT,H_sys_mix,sz,shift)
		end select
		!
		call assign_psi(psi)
	endif
	!
	call evo_rho_diab(rho_diag_diab)
	!-----------------------------------
	!
	call evo_npi(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_npi, pop_p_npi)
	call final_psi_hop_loc01_dt(Ut, Tv, psi_f, pop_p_npi_t)
	call evo_hst(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_hst, pop_p_hst)
	call evo_loc01(Ut, Tv)
	call final_rho_hop(Ut, Tv, rho_f_loc, pop_p_loc01_l)
	call final_psi_hop_loc01(Ut, Tv, psi_f, pop_p_loc01)
	call final_psi_hop_loc01_dt(Ut, Tv, psi_f, pop_p_loc01_t)
	call final_rho_hop_loc19(Ut, Tv, rho_f_loc, pop_p_loc19)
	!
	! use final rho of NPI as reference
	allocate(rho_diag_npi_dq(sz))
	allocate(rho_diag_npi(sz))
	allocate(rho_diag_hst(sz))
	allocate(rho_diag_loc(sz))
	do istate = 1, sz
		rho_diag_npi_dq(istate) = dble(rho_f_npi_dq(istate,istate))
		rho_diag_npi(istate) = dble(rho_f_npi(istate,istate))
		rho_diag_hst(istate) = dble(rho_f_hst(istate,istate))
		rho_diag_loc(istate) = dble(rho_f_loc(istate,istate))
	enddo
	!
	write(f_sz,*) sz
	!
	write(*,*) 'Final population from rho'
	write(*,'("DIAB     " '//adjustl(f_sz)//'(ES18.10,1X))') rho_diag_diab
	write(*,'("NPI-dqT  " '//adjustl(f_sz)//'(ES18.10,1X))') rho_diag_npi_dq
	write(*,'("NPI      " '//adjustl(f_sz)//'(ES18.10,1X))') rho_diag_npi
	write(*,'("HST      " '//adjustl(f_sz)//'(ES18.10,1X))') rho_diag_hst
	write(*,'("LOC      " '//adjustl(f_sz)//'(ES18.10,1X))') rho_diag_loc
	write(*,*) ''
	!
	write(*,*) 'Population from actual hopping'
	write(*,'("NPI-dqT  " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_npi_dq
	write(*,'("NPI-cdqT " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_npi_cdq
	write(*,'("NPI-t-dqT" '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_npi_t_dq
	write(*,'("NPI      " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_npi
	write(*,'("NPI-t    " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_npi_t
	write(*,'("HST      " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_hst
	write(*,'("LOC01    " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_loc01
	write(*,'("LOC01-l  " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_loc01_l
	write(*,'("LOC01-t  " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_loc01_t
	write(*,'("LOC19    " '//adjustl(f_sz)//'(ES18.10,1X))') pop_p_loc19
	!
	write(*,*) nT, nT+num_extra_call
	!call test_E
	!call test_H
end program test
