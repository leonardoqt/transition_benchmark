module benchmark_system
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	real(dp)  :: v0, dt
	integer   :: nT, nstate
	real(dp)   , allocatable :: U(:,:,:), E(:,:), H(:,:,:), Uini(:,:)
	complex(dp), allocatable :: psi0(:), rho0(:,:) ! rho is used in real calculation
	!
	contains
	!
	!
	subroutine assign_model(v_in, x_ini_in, x_fin_in, dt_in, H_diabat, nstate_in, shift_in)
		!
		implicit none
		!
		real(dp) :: v_in, x_ini_in, x_fin_in, dt_in, shift_in
		integer  :: nstate_in
		!
		! working variables
		integer  :: idt, istate, jstate, n_bad_diagonal, imax, imin
		real(dp), allocatable :: S(:,:)
		!
		interface 
			function H_diabat(x,sz,shift)
				use iso_fortran_env, only: dp=> real64
				implicit none
				real(dp) :: x, shift
				integer  :: sz
				real(dp), allocatable :: H_diabat(:,:)
			end function
		end interface
		!
		!
		v0 = v_in
		dt = dt_in
		nstate = nstate_in
		nT = ceiling( (x_fin_in - x_ini_in) / (v0*dt) )
		!
		if ( allocated(U) ) deallocate(U)
		if ( allocated(E) ) deallocate(E)
		if ( allocated(H) ) deallocate(H)
		if ( allocated(Uini) ) deallocate(Uini)
		allocate( U(nstate,nstate,nT) )
		allocate( E(nstate,nT)        )
		allocate( H(nstate,nstate,nT) )
		allocate( Uini(nstate,nstate) )
		if ( allocated(psi0) ) deallocate(psi0)
		if ( allocated(rho0) ) deallocate(rho0)
		allocate( psi0(nstate) )
		allocate( rho0(nstate,nstate) )
		!
		! diagonalize Hamiltonian to get eigenvector (U) and eigenvalue(E)
		do idt = 1, nT
			H(:,:,idt) = H_diabat(x_ini_in+v0*(idt-1)*dt,nstate,shift_in)
			call diag_real( H(:,:,idt), U(:,:,idt), E(:,idt) )
		enddo
		Uini = U(:,:,1)
		!
		! naive parallel transport of U followed by ZyZ's simplified algorithm
		allocate( S(nstate,nstate) )
		do idt = 2, nT
			! naive parallel transport
			do istate = 1, nstate
				U(:,istate,idt) = U(:,istate,idt) * sign( 1.d0, dot_product(U(:,istate,idt-1),U(:,istate,idt)) )
			enddo
			!
			! ZY's scheme
			S = matmul( transpose(U(:,:,idt-1)), U(:,:,idt) )
			call ZY_correct_sign_full(S,U(:,:,idt))
			!call correct_sign_bruteforce(S,U(:,:,idt))
			!call rand_sign(U(:,:,idt))
		enddo
		!
		deallocate(S)
		!
	end subroutine assign_model
	!
	!
	subroutine assign_psi(psi_in)
		!
		implicit none
		!
		complex(dp), dimension(:) :: psi_in
		!
		!
		psi0 = psi_in
		rho0(:,1) = psi0
		rho0 = matmul( rho0(:,1:1), conjg(transpose(rho0(:,1:1))) )
		!
	end subroutine
	!
	subroutine assign_rho(rho_in)
		!
		implicit none
		!
		complex(dp), dimension(:,:) :: rho_in
		!
		!
		rho0 = rho_in
		!
	end subroutine
	!
	!
	subroutine evo_rho_diab(rho_f)
		! it generates Ut and Tv rate, such that Ut*rho_diabat*Ut^\dagger is rho(t)_adiabats
		! and Tv*dt is the hopping rate
		!
		implicit none
		!
		real(dp)   , allocatable :: rho_f(:)
		!
		complex(dp), allocatable :: Ut(:,:), U_dt(:,:), iHTdt(:,:), rho_full(:,:)
		integer :: idt, istate
		!
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		allocate( rho_f(nstate) )
		!
		allocate(    Ut(nstate,nstate) )
		allocate(  U_dt(nstate,nstate) )
		allocate( iHTdt(nstate,nstate) )
		allocate( rho_full(nstate,nstate) )
		!
		Ut = (0.d0, 0.d0)
		do istate = 1,nstate
			Ut(istate,istate) = (1.d0, 0.d0)
		enddo
		!
		do idt = 1, nT-1
			iHTdt = cmplx( H(:,:,1)*0.d0, -dt/2*(H(:,:,idt)+H(:,:,idt+1)), dp )
			U_dt = expm(iHTdt)
			Ut = matmul(U_dt, Ut)
		enddo
		Ut = matmul( transpose( U(:,:,nT) ), Ut)
		rho_full = matmul( matmul(Ut, rho0), conjg(transpose(Ut)) )
		!
		do istate = 1,nstate
			rho_f(istate) = dble( rho_full(istate,istate) )
		enddo
		!
		deallocate(Ut,U_dt,iHTdt,rho_full)
		!
	end subroutine evo_rho_diab
	!
	!
	!
	subroutine evo_npi(Ut, Tv)
		! it generates Ut and Tv rate, such that Ut*rho_diabat*Ut^\dagger is rho(t)_adiabats
		! and Tv*dt is the hopping rate
		!
		implicit none
		!
		complex(dp), allocatable :: Ut(:,:,:)
		real(dp)   , allocatable :: Tv(:,:,:)
		!
		complex(dp), allocatable :: U_dt(:,:), iHTdt(:,:)
		integer :: idt, istate
		!
		!
		if ( allocated(Ut) ) deallocate(Ut)
		if ( allocated(Tv) ) deallocate(Tv)
		allocate( Ut(nstate,nstate,nT  ) )
		allocate( Tv(nstate,nstate,nT-1) )
		!
		allocate(  U_dt(nstate,nstate) )
		allocate( iHTdt(nstate,nstate) )
		!
		do idt = 1, nT-1
			Tv(:,:,idt) = logm( matmul(transpose(U(:,:,idt)),U(:,:,idt+1)) ) / dt
		enddo
		!
		! evolution of rho
		Ut(:,:,1) = transpose( U(:,:,1) )
		do idt = 1, nT-1
			iHTdt = (0.d0, 0.d0)
			do istate = 1, nstate
				iHTdt(istate,istate) = cmplx(0.d0, -(E(istate,idt)+E(istate,idt+1))/2*dt, dp)
			enddo
			iHTdt = iHTdt - Tv(:,:,idt)*dt
			U_dt = expm(iHTdt)
			Ut(:,:,idt+1) = matmul(U_dt, Ut(:,:,idt))
		enddo
		!
		deallocate(U_dt,iHTdt)
		!
	end subroutine evo_npi
	!
	subroutine evo_npi_interp(Ut, Tv, nqT)
		! it generates Ut and Tv rate, such that Ut*rho_diabat*Ut^\dagger is rho(t)_adiabats
		! and Tv*dt is the hopping rate
		!
		implicit none
		!
		complex(dp), allocatable :: Ut(:,:,:)
		real(dp)   , allocatable :: Tv(:,:,:)
		integer                  :: nqT ! number of quantum time steps in each ionic time step
		!
		complex(dp), allocatable :: U_dt(:,:), iHTdt(:,:)
		integer :: idt, iqt, istate
		!
		!
		if ( allocated(Ut) ) deallocate(Ut)
		if ( allocated(Tv) ) deallocate(Tv)
		allocate( Ut(nstate,nstate,(nT-1)*nqT+1 ) )
		allocate( Tv(nstate,nstate,(nT-1)*nqT   ) )
		!
		allocate(  U_dt(nstate,nstate) )
		allocate( iHTdt(nstate,nstate) )
		!
		do idt = 1, nT-1
			Tv(:,:,(idt-1)*nqT+1) = logm( matmul(transpose(U(:,:,idt)),U(:,:,idt+1)) ) / dt
			do iqt = 2, nqT
				Tv(:,:,(idt-1)*nqT+iqt) = Tv(:,:,(idt-1)*nqT+1)
			enddo
		enddo
		!
		! evolution of rho
		Ut(:,:,1) = transpose( U(:,:,1) )
		do idt = 1, nT-1
			do iqt = 1, nqT
				iHTdt = (0.d0, 0.d0)
				do istate = 1, nstate
					iHTdt(istate,istate) = cmplx(0.d0, -( E(istate,idt)+(E(istate,idt+1)-E(istate,idt))*((iqt-0.5d0)/nqT) )*(dt/nqT), dp)
				enddo
				iHTdt = iHTdt - Tv(:,:,(idt-1)*nqT+iqt)*(dt/nqT)
				U_dt = expm(iHTdt)
				Ut(:,:,(idt-1)*nqT+iqt+1) = matmul( U_dt, Ut(:,:,(idt-1)*nqT+iqt) )
			enddo
		enddo
		!
		deallocate(U_dt,iHTdt)
		!
	end subroutine evo_npi_interp
	!
	!
	subroutine evo_hst(Ut, Tv)
		! it generates Ut and Tv rate, such that Ut*rho_diabat*Ut^\dagger is rho(t)_adiabats
		! and Tv*dt is the hopping rate
		!
		implicit none
		!
		complex(dp), allocatable :: Ut(:,:,:)
		real(dp)   , allocatable :: Tv(:,:,:)
		!
		complex(dp), allocatable :: U_dt(:,:), iHTdt(:,:)
		integer :: idt, istate
		!
		!
		if ( allocated(Ut) ) deallocate(Ut)
		if ( allocated(Tv) ) deallocate(Tv)
		allocate( Ut(nstate,nstate,nT  ) )
		allocate( Tv(nstate,nstate,nT-1) )
		!
		allocate(  U_dt(nstate,nstate) )
		allocate( iHTdt(nstate,nstate) )
		!
		do idt = 1, nT-1
			Tv(:,:,idt) = ( matmul(transpose(U(:,:,idt)),U(:,:,idt+1)) - matmul(transpose(U(:,:,idt+1)),U(:,:,idt)) ) / (2*dt)
		enddo
		!
		! evolution of rho
		Ut(:,:,1) = transpose( U(:,:,1) )
		do idt = 1, nT-1
			iHTdt = (0.d0, 0.d0)
			do istate = 1, nstate
				iHTdt(istate,istate) = cmplx(0.d0, -(E(istate,idt)+E(istate,idt+1))/2*dt, dp)
			enddo
			iHTdt = iHTdt - Tv(:,:,idt)*dt
			U_dt = expm(iHTdt)
			Ut(:,:,idt+1) = matmul(U_dt, Ut(:,:,idt))
		enddo
		!
		deallocate(U_dt,iHTdt)
		!
	end subroutine evo_hst
	!
	!
	subroutine evo_loc01(Ut, Tv)
		! https://doi.org/10.1063/1.1376633
		!
		! it generates Ut and Tv rate, such that Ut*rho_diabat*Ut^\dagger is rho(t)_adiabats
		! and Tv*dt is the hopping rate
		!
		implicit none
		!
		complex(dp), allocatable :: Ut(:,:,:)
		real(dp)   , allocatable :: Tv(:,:,:)
		!
		real(dp)   , allocatable :: S(:,:),vec(:,:),ww(:,:),val(:),Hloc(:,:)
		complex(dp), allocatable :: U_dt(:,:), iHTdt(:,:)
		integer :: idt, istate
		!
		!
		if ( allocated(Ut) ) deallocate(Ut)
		if ( allocated(Tv) ) deallocate(Tv)
		allocate( Ut(nstate,nstate,nT  ) )
		allocate( Tv(nstate,nstate,nT-1) )
		!
		allocate(     S(nstate,nstate) )
		allocate(    ww(nstate,nstate) )
		allocate(   vec(nstate,nstate) )
		allocate(   val(nstate) )
		allocate(  Hloc(nstate,nstate) )
		allocate(  U_dt(nstate,nstate) )
		allocate( iHTdt(nstate,nstate) )
		!
		do idt = 1, nT-1
			S = matmul(transpose(U(:,:,idt)),U(:,:,idt+1))
			call diag_real( matmul(transpose(S),S), vec, val )
			do istate = 1, nstate
				ww(:,istate) = vec(:,istate) / sqrt(val(istate))
			enddo
			Tv(:,:,idt) = matmul( matmul(S,ww), transpose(vec) )
		enddo
		!
		! evolution of rho
		Ut(:,:,1) = transpose( U(:,:,1) )
		do idt = 1, nT-1
			do istate = 1, nstate
				Hloc(:,istate) = Tv(:,istate,idt) * E(istate,idt+1)
			enddo
			Hloc = matmul(Hloc, transpose(Tv(:,:,idt)))
			!
			iHTdt = (0.d0, 0.d0)
			do istate = 1, nstate
				iHTdt(istate,istate) = cmplx(0.d0, -E(istate,idt+1)/2*dt, dp)
			enddo
			iHTdt = iHTdt - Hloc*( (0.d0, 5.d-1) * dt )
			U_dt = expm(iHTdt)
			U_dt = matmul(transpose(Tv(:,:,idt)), U_dt)
			Ut(:,:,idt+1) = matmul(U_dt, Ut(:,:,idt))
		enddo
		!
		!
		! the following uses (19) to the first order in dt
		! TODO: the problem of (19) is how to apply it to the mixed state
		do idt = 1, nT-1
			Tv(:,:,idt) = Tv(:,:,idt) / dt
			do istate = 1, nstate
				Tv(istate,istate,idt) = 0.d0
			enddo
		enddo
		!
		deallocate(S,ww,vec,val,Hloc,U_dt,iHTdt)
		!
	end subroutine evo_loc01
	!
	!
	!
	subroutine final_rho_hop(Ut, Tv, rho_f, pop_p)
		! use given Ut, Tv and stored rho0 to calculate final rho and the 
		! also calculate final population based on hops
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		!
		complex(dp), allocatable       :: rho_f(:,:), U_dt(:,:)
		real(dp)   , allocatable       :: pop_p(:,:)
		!
		real(dp)                       :: drate
		real(dp)   , allocatable       :: T_trans(:,:), p_trans(:,:)
		integer                        :: idt, istate, jstate
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( rho_f(nstate,nstate) )
		allocate( pop_p(nstate,1) )
		!
		allocate(  U_dt(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		!
		! initial pop_p is the diagonal of rho_adiabat
		rho_f = matmul( matmul(Ut(:,:,1), rho0), transpose(conjg(Ut(:,:,1))) )
		do istate = 1, nstate
			pop_p(istate,1) = dble(rho_f(istate,istate))
		enddo
		!
		do idt = 1, nT-1
			rho_f = matmul( matmul(Ut(:,:,idt+1), rho0), transpose(conjg(Ut(:,:,idt+1))) )
			U_dt = matmul( Ut(:,:,idt+1), transpose(conjg( Ut(:,:,idt) )) )
			!
			! for calculating pop_p
			! T_trans
			T_trans = 0.d0
			do istate = 1, nstate
				do jstate = 1, nstate
					if ( jstate .ne. istate ) then
						drate = 2 * dble( Tv(istate,jstate,idt)*rho_f(jstate,istate) ) / dble( rho_f(istate,istate) + 1.d-13)
						if ( drate > 0.d0 ) T_trans(istate,jstate) = drate
					endif
				enddo
			enddo
			!
			! normalize Ti: if sum greater than 1
			do istate = 1, nstate
				if ( sum(T_trans(istate,:))*dt > 1.d0 ) T_trans(istate,:) = T_trans(istate,:) / ( sum(T_trans(istate,:))*dt )
			enddo
			!
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			!
			! naive integration for dp = p_trans*p*dt
			pop_p = pop_p + matmul(p_trans,pop_p)*dt
		enddo
		!
		deallocate(U_dt,T_trans,p_trans)
		!
	end subroutine final_rho_hop
	!
	subroutine final_rho_hop_interp(Ut, Tv, rho_f, pop_p, nqT)
		! use given Ut, Tv and stored rho0 to calculate final rho and the 
		! calculate final population based on hops
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		!
		complex(dp), allocatable       :: rho_f(:,:), U_dt(:,:)
		real(dp)   , allocatable       :: pop_p(:,:)
		integer                        :: nqT
		!
		real(dp)                       :: drate
		real(dp)   , allocatable       :: T_trans(:,:), p_trans(:,:)
		integer                        :: idt, istate, jstate
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( rho_f(nstate,nstate) )
		allocate( pop_p(nstate,1) )
		!
		allocate(  U_dt(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		!
		! initial pop_p is the diagonal of rho_adiabat
		rho_f = matmul( matmul(Ut(:,:,1), rho0), transpose(conjg(Ut(:,:,1))) )
		do istate = 1, nstate
			pop_p(istate,1) = dble(rho_f(istate,istate))
		enddo
		!
		do idt = 1, (nT-1)*nqT
			rho_f = matmul( matmul(Ut(:,:,idt+1), rho0), transpose(conjg(Ut(:,:,idt+1))) )
			U_dt = matmul( Ut(:,:,idt+1), transpose(conjg( Ut(:,:,idt) )) )
			!
			! for calculating pop_p
			! T_trans
			T_trans = 0.d0
			do istate = 1, nstate
				do jstate = 1, nstate
					if ( jstate .ne. istate ) then
						drate = 2 * dble( Tv(istate,jstate,idt)*rho_f(jstate,istate) ) / dble( rho_f(istate,istate) + 1.d-13)
						if ( drate > 0.d0 ) T_trans(istate,jstate) = drate
					endif
				enddo
			enddo
			!
			! normalize Ti: if sum greater than 1
			do istate = 1, nstate
				if ( sum(T_trans(istate,:))*dt/nqT > 1.d0 ) T_trans(istate,:)=T_trans(istate,:)/(sum(T_trans(istate,:))*dt/nqT)
			enddo
			!
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			!
			! naive integration for dp = p_trans*p*dt
			pop_p = pop_p + matmul(p_trans,pop_p)*(dt/nqT)
		enddo
		!
		deallocate(U_dt,T_trans,p_trans)
		!
	end subroutine final_rho_hop_interp
	!
	!
	subroutine final_rho_hop_conditional_interp(rho_f, pop_p, threshold, num_extra_call)
		! instead of using Ut and T, it directly use U and E to calculate Ut and T
		! therefore no evo_xxx needs to be called beforehand
		!
		implicit none
		!
		complex(dp), allocatable       :: rho_f(:,:)
		real(dp)   , allocatable       :: pop_p(:,:), Tv(:,:)
		real(dp)                       :: threshold
		integer                        :: num_extra_call
		!
		integer                        :: idt, istate
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( rho_f(nstate,nstate) )
		allocate( pop_p(nstate,1) )
		allocate( Tv(nstate,nstate) )
		!
		rho_f = matmul( matmul( transpose(U(:,:,1)),rho0 ), U(:,:,1) )
		do istate = 1, nstate
			pop_p(istate,1) = dble(rho_f(istate,istate))
		enddo
		!
		do idt = 1, nT-1
			Tv = logm( matmul(transpose(U(:,:,idt)),U(:,:,idt+1)) ) / dt
			call single_rho_hop_interp(E(:,idt),E(:,idt+1),Tv,dt,rho_f,pop_p,threshold,num_extra_call)
			!call single_rho_hop_interp_exact(E(:,idt),E(:,idt+1),Tv,dt,rho_f,pop_p,num_extra_call)
		enddo
		!
		deallocate(Tv)
		!
	end subroutine
	!
	recursive subroutine single_rho_hop_interp(E1,E2,Tv,dt_interp,rho_t,pop_t,threshold,num_extra_call)
		! calculate the evolution of rho_t and pop_t of a single time interval
		! if hop is too big, interpolate with dt_interp/10
		!
		implicit none
		!
		real(dp), dimension(:)      :: E1, E2
		real(dp), dimension(:,:)    :: Tv, pop_t
		real(dp)                    :: dt_interp,threshold
		complex(dp), dimension(:,:) :: rho_t
		integer                     :: num_extra_call
		!
		complex(dp), allocatable    :: iHTdt(:,:), U_dt(:,:), rho_new(:,:)
		real(dp)   , allocatable    :: T_trans(:,:), p_trans(:,:)
		real(dp)                    :: drate
		integer                     :: istate, jstate, t1, all_good
		!
		allocate( iHTdt  (nstate,nstate) )
		allocate( U_dt   (nstate,nstate) )
		allocate( rho_new(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		!
		iHTdt = (0.d0, 0.d0)
		do istate = 1, nstate
			iHTdt(istate,istate) = cmplx(0.d0, -(E1(istate)+E2(istate))/2*dt_interp, dp)
		enddo
		iHTdt = iHTdt - Tv * dt_interp
		U_dt = expm(iHTdt)
		rho_new = matmul( matmul(U_dt,rho_t),transpose(conjg(U_dt)) )
		!
		T_trans = 0.d0
		do istate = 1, nstate
			do jstate = 1, nstate
				if ( jstate .ne. istate ) then
					drate = 2 * dble( Tv(istate,jstate)*rho_new(jstate,istate) ) / dble( rho_t(istate,istate) + 1.d-13)
					if ( drate > 0.d0 ) T_trans(istate,jstate) = drate*dt_interp
				endif
			enddo
		enddo
		!
		! check if all good
		all_good = 1
		do istate = 1, nstate
			if (sum(T_trans(istate,:)) > threshold) then
				all_good = 0
				exit
			endif
		enddo
		!
		if ( all_good > 0 .or. dt_interp < dt / 1.d4 ) then
			! do not need interpolate
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			pop_t = pop_t + matmul(p_trans,pop_t)
			rho_t = rho_new
		else
			! need interpolate
			num_extra_call = num_extra_call+10
			do t1 = 1,10
				call single_rho_hop_interp(E1+(E2-E1)*(t1-1)/10.d0,E1+(E2-E1)*t1/10.d0,Tv,dt_interp/10,rho_t,pop_t,threshold,num_extra_call)
			enddo
		endif
		!
		deallocate(iHTdt,U_dt,rho_new,T_trans,p_trans)
		!
	end subroutine single_rho_hop_interp
	!
	!
	recursive subroutine single_rho_hop_interp_exact(E1,E2,Tv,dt_interp,rho_t,pop_t,num_extra_call)
		! calculate the evolution of rho_t and pop_t of a single time interval
		! if hop is too big, interpolate with dt_interp/10
		!
		implicit none
		!
		real(dp), dimension(:)      :: E1, E2
		real(dp), dimension(:,:)    :: Tv, pop_t
		real(dp)                    :: dt_interp
		complex(dp), dimension(:,:) :: rho_t
		integer                     :: num_extra_call
		!
		complex(dp), allocatable    :: iHTdt(:,:), U_dt(:,:), rho_new(:,:)
		real(dp)   , allocatable    :: T_trans(:,:), p_trans(:,:), WW(:,:)
		real(dp)                    :: drate
		integer                     :: istate, jstate, t1, all_good
		!
		allocate( iHTdt  (nstate,nstate) )
		allocate( U_dt   (nstate,nstate) )
		allocate( rho_new(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		allocate( WW     (nstate,nstate) )
		!
		iHTdt = (0.d0, 0.d0)
		do istate = 1, nstate
			iHTdt(istate,istate) = cmplx(0.d0, -(E1(istate)+E2(istate))/2*dt_interp, dp)
		enddo
		iHTdt = iHTdt - Tv * dt_interp
		U_dt = expm(iHTdt)
		rho_new = matmul( matmul(U_dt,rho_t),transpose(conjg(U_dt)) )
		!
		T_trans = 0.d0
		do istate = 1, nstate
			do jstate = 1, nstate
				if ( jstate .ne. istate ) then
					drate = 2 * dble( Tv(istate,jstate)*rho_new(jstate,istate) ) / dble( rho_t(istate,istate) + 1.d-13)
					if ( drate > 0.d0 ) T_trans(istate,jstate) = drate*dt_interp
				endif
			enddo
		enddo
		!
		! generate p_trans
		do istate = 1, nstate
			if ( sum(T_trans(istate,:)) > 1.d0 ) T_trans(istate,:) = T_trans(istate,:) / sum(T_trans(istate,:))
		enddo
		p_trans = transpose(T_trans)
		do istate = 1, nstate
			p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
		enddo
		!
		! find the exact WW
		call generate_exact_p_increment(rho_t,rho_new,WW,p_trans)
		! Note that W = p_trans + I
		!
		if ( minval(WW) > -1.d-7 ) then
			! do not need interpolate
			p_trans = WW
			do istate = 1,nstate
				p_trans(istate,istate) = p_trans(istate,istate) - 1.d0
			enddo
			pop_t = pop_t + matmul(p_trans,pop_t)
			rho_t = rho_new
		elseif( dt_interp < dt/1.d4 ) then
			pop_t = pop_t + matmul(p_trans,pop_t)
			rho_t = rho_new
		else
			! need interpolate
			num_extra_call = num_extra_call + 10
			do t1 = 1,10
				call single_rho_hop_interp_exact(E1+(E2-E1)*(t1-1)/10.d0,E1+(E2-E1)*t1/10.d0,Tv,dt_interp/10,rho_t,pop_t,num_extra_call)
			enddo
		endif
		!
		deallocate(iHTdt,U_dt,rho_new,T_trans,p_trans,WW)
		!
	end subroutine single_rho_hop_interp_exact
	!
	!
	subroutine final_psi_hop_interp_dt(Ut, Tv, psi_f, pop_p, nqT)
		! use given Ut, Tv and stored rho0 to calculate final rho and the 
		! calculate final population based on hops
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		!
		complex(dp), allocatable       :: psi_f(:,:), psi_l(:,:), psi00(:,:)
		real(dp)   , allocatable       :: pop_p(:,:)
		integer                        :: nqT
		!
		real(dp)                       :: drate
		real(dp)   , allocatable       :: T_trans(:,:), p_trans(:,:)
		integer                        :: idt, istate, jstate
		!
		if ( allocated(psi_f) ) deallocate(psi_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( psi_f(nstate,1) )
		allocate( pop_p(nstate,1) )
		!
		allocate( psi_l(nstate,1) )
		allocate( psi00(nstate,1) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		!
		! initial pop_p is the diagonal of rho_adiabat
		psi00(:,1) = psi0
		psi_l = matmul(Ut(:,:,1), psi00)
		do istate = 1, nstate
			pop_p(istate,1) = dble(psi_l(istate,1)*conjg(psi_l(istate,1)))
		enddo
		!
		do idt = 1, (nT-1)*nqT
			psi_f = matmul(Ut(:,:,idt+1), psi00)
			!
			! for calculating pop_p
			! T_trans
			T_trans = 0.d0
			do istate = 1, nstate
				do jstate = 1, nstate
					if ( jstate .ne. istate ) then
						drate = 2 * dble( Tv(istate,jstate,idt)*psi_f(jstate,1)*conjg(psi_l(istate,1)) ) / &
						        dble( psi_f(istate,1)*conjg(psi_f(istate,1)) + 1.d-13)
						if ( drate > 0.d0 ) T_trans(istate,jstate) = drate
					endif
				enddo
			enddo
			!
			! normalize Ti: if sum greater than 1
			do istate = 1, nstate
				if ( sum(T_trans(istate,:))*dt/nqT > 1.d0 ) T_trans(istate,:)=T_trans(istate,:)/(sum(T_trans(istate,:))*dt/nqT)
			enddo
			!
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			!
			! naive integration for dp = p_trans*p*dt
			pop_p = pop_p + matmul(p_trans,pop_p)*(dt/nqT)
			!
			psi_l = psi_f
		enddo
		!
		deallocate(T_trans,p_trans,psi_l,psi00)
		!
	end subroutine final_psi_hop_interp_dt
	!
	!
	subroutine final_rho_hop_loc19(Ut, Tv, rho_f, pop_p)
		! https://doi.org/10.1016/j.comptc.2019.02.009
		!
		! use given Ut, Tv and stored rho0 to calculate final rho and the 
		! total hopping rate from state0 to all states
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		!
		complex(dp), allocatable       :: rho_f(:,:)
		real(dp)   , allocatable       :: pop_p(:,:)
		!
		complex(dp), allocatable       :: rho_l(:,:)
		real(dp)   , allocatable       :: T_trans(:,:), p_trans(:,:)
		real(dp)                       :: drate, wk, Sk, Pdt, xkj
		integer                        :: idt, istate, jstate
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( rho_f(nstate,nstate) )
		allocate( pop_p(nstate,1) )
		!
		allocate( rho_l(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		!
		! initial pop_p is the diagonal of rho_adiabat
		rho_l = matmul( matmul(Ut(:,:,1), rho0), transpose(conjg(Ut(:,:,1))) )
		do istate = 1, nstate
			pop_p(istate,1) = dble(rho_l(istate,istate))
		enddo
		!
		!hop_p = 0.d0
		do idt = 1, nT-1
			rho_f = matmul( matmul(Ut(:,:,idt+1), rho0), transpose(conjg(Ut(:,:,idt+1))) )
			!
			! calculate T_trans and pop_p
			T_trans = 0.d0
			do istate = 1, nstate
				wk = dble( rho_l(istate,istate)-rho_f(istate,istate) ) / dble( rho_l(istate,istate) )
				if (wk > 0.d0) then
					!
					! calculate Sk
					Sk = 0.d0
					do jstate = 1, nstate
						if ( jstate .ne. istate) then
							! Tv*dt is the T in the reference
							Pdt = dble( rho_f(jstate,jstate)-rho_l(jstate,jstate) )
							if ( Pdt > 0.d0 ) then
								! dt is not necessary here, since Sk and xkj serve as a partition scheme
								Sk = Sk + abs(Tv(istate,jstate,idt))*sqrt(Pdt)
							endif
						endif
					enddo
					!
					! calculate xkj
					do jstate = 1, nstate
						if ( jstate .ne. istate) then
							Pdt = dble( rho_f(jstate,jstate)-rho_l(jstate,jstate) )
							if ( Pdt > 0.d0 ) then
								xkj = abs(Tv(istate,jstate,idt))*sqrt(Pdt) / Sk
								T_trans(istate,jstate) = xkj * wk / dt
								!if ( istate .eq. state0 ) hop_p(jstate) = hop_p(jstate) + xkj * wk
							endif
						endif
					enddo
				endif
			enddo
			!
			! normalize Ti: if sum greater than 1
			do istate = 1, nstate
				if ( sum(T_trans(istate,:))*dt > 1.d0 ) T_trans(istate,:) = T_trans(istate,:) / ( sum(T_trans(istate,:))*dt )
			enddo
			!
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			!
			! naive integration for dp = p_trans*p*dt
			pop_p = pop_p + matmul(p_trans,pop_p)*dt
			!
			rho_l = rho_f
		enddo
		!
		deallocate(rho_l,T_trans,p_trans)
		!
	end subroutine final_rho_hop_loc19
	!
	!
	subroutine final_psi_hop_loc01(Ut, Tv, psi_f, pop_p)
		! use given Ut, Tv and stored psi0 to calculate final psi and the 
		! total hopping rate from state0 to all states
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		!
		complex(dp), allocatable       :: psi_f(:,:), psi_l(:,:), U_dt(:,:)
		real(dp)   , allocatable       :: pop_p(:,:)
		!
		real(dp)   , allocatable       :: T_trans(:,:), p_trans(:,:)
		real(dp)                       :: drate
		integer                        :: idt, istate, jstate
		!
		if ( allocated(psi_f) ) deallocate(psi_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( psi_f(nstate,1) )
		allocate( pop_p(nstate,1) )
		!
		allocate( psi_l(nstate,1) )
		allocate( U_dt(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		!
		psi_l(:,1) = psi0
		psi_l = matmul(Ut(:,:,1), psi_l)
		do istate = 1, nstate
			pop_p(istate,1) = dble(psi_l(istate,1)*conjg(psi_l(istate,1)))
		enddo
		!
		!hop_p = 0.d0
		do idt = 1, nT-1
			U_dt = matmul( Ut(:,:,idt+1), transpose(conjg(Ut(:,:,idt))) )
			psi_f = matmul(U_dt, psi_l)
			!
			T_trans = 0.d0
			do istate = 1, nstate
				do jstate = 1, nstate
					if ( jstate .ne. istate ) then
						drate = -dble( U_dt(istate,jstate)*psi_l(jstate,1)*conjg(psi_f(istate,1)) ) * &
						         dble( psi_f(istate,1)*conjg(psi_f(istate,1)) - psi_l(istate,1)*conjg(psi_l(istate,1)) )
						drate = drate / ( dble( psi_f(istate,1)*conjg(psi_f(istate,1)) ) - &
						                  dble( U_dt(istate,istate)*psi_l(istate,1)*conjg(psi_f(istate,1)) )+1.d-13 )
						drate = drate/dble( psi_l(istate,1)*conjg(psi_l(istate,1)) + 1.d-13 )
						!if ( drate > 1.d0 ) drate = 1.d0
						if ( drate > 0.d0 ) then
							T_trans(istate,jstate) = drate / dt
							!if ( istate .eq. state0 ) hop_p(jstate) = hop_p(jstate) + drate
						endif
					endif
				enddo
			enddo
			!
			! normalize Ti: if sum greater than 1
			do istate = 1, nstate
				if ( sum(T_trans(istate,:))*dt > 1.d0 ) T_trans(istate,:) = T_trans(istate,:) / ( sum(T_trans(istate,:))*dt )
			enddo
			!
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			!
			! naive integration for dp = p_trans*p*dt
			pop_p = pop_p + matmul(p_trans,pop_p)*dt
			!
			psi_l = psi_f
		enddo
		!
		deallocate( psi_l,U_dt,T_trans,p_trans )
		!
	end subroutine final_psi_hop_loc01
	!
	!
	!
	subroutine final_psi_hop_loc01_dt(Ut, Tv, psi_f, pop_p)
		! use given Ut, Tv and stored psi0 to calculate final psi and the 
		! total hopping rate from state0 to all states
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		integer                        :: state0
		!
		complex(dp), allocatable       :: psi_f(:,:), psi_l(:,:), U_dt(:,:)
		real(dp)   , allocatable       :: pop_p(:,:)
		!
		real(dp)   , allocatable       :: T_trans(:,:), p_trans(:,:)
		real(dp)                       :: drate
		real(dp)   , allocatable       :: WW(:,:), p_trans_opt(:,:)
		integer                        :: idt, istate, jstate
		!
		if ( allocated(psi_f) ) deallocate(psi_f)
		if ( allocated(pop_p) ) deallocate(pop_p)
		allocate( psi_f(nstate,1) )
		allocate( pop_p(nstate,1) )
		!
		allocate( psi_l(nstate,1) )
		allocate( U_dt(nstate,nstate) )
		allocate( T_trans(nstate,nstate) )
		allocate( p_trans(nstate,nstate) )
		allocate( WW     (nstate,nstate) )
		allocate( p_trans_opt(nstate,nstate) )
		!
		psi_l(:,1) = psi0
		psi_l = matmul(Ut(:,:,1), psi_l)
		do istate = 1, nstate
			pop_p(istate,1) = dble(psi_l(istate,1)*conjg(psi_l(istate,1)))
		enddo
		!
		!hop_p = 0.d0
		do idt = 1, nT-1
			U_dt = matmul( Ut(:,:,idt+1), transpose(conjg(Ut(:,:,idt))) )
			psi_f = matmul(U_dt, psi_l)
			!
			T_trans = 0.d0
			do istate = 1, nstate
				do jstate = 1, nstate
					if ( jstate .ne. istate ) then
						drate = 2*dble( Tv(istate,jstate,idt)*psi_f(jstate,1)*conjg(psi_l(istate,1)) ) / &
						        dble(psi_f(istate,1)*conjg(psi_f(istate,1)) + 1.d-13)
						if ( drate > 0.d0 ) then
							T_trans(istate,jstate) = drate
							!if ( istate .eq. state0 ) hop_p(jstate) = hop_p(jstate) + drate
						endif
					endif
				enddo
			enddo
			!
			! normalize Ti: if sum greater than 1
			do istate = 1, nstate
				if ( sum(T_trans(istate,:))*dt > 1.d0 ) T_trans(istate,:) = T_trans(istate,:) / ( sum(T_trans(istate,:))*dt )
			enddo
			!
			p_trans = transpose(T_trans)
			do istate = 1, nstate
				p_trans(istate,istate) = p_trans(istate,istate) - sum(T_trans(istate,:))
			enddo
			!
			!========================================
			! naive integration for dp = p_trans*p*dt
			pop_p = pop_p + matmul(p_trans,pop_p)*dt
			!-------------------
			!!
			!!call generate_exact_p_increment(psi_l,psi_f,WW,p_trans*dt)
			!call generate_approx_p(psi_l,psi_f,WW,p_trans*dt)
			!! Note that W = P + I
			!!
			!! safeguard for still negative matrix elements
			!do istate = 1,nstate
			!do jstate = 1,nstate
			!	if ( WW(istate,jstate) < 0.d0) WW(istate,jstate) = 0.d0
			!enddo
			!enddo
			!!
			!p_trans_opt = WW
			!do istate = 1,nstate
			!	p_trans_opt(istate,istate) = p_trans_opt(istate,istate) - sum(p_trans_opt(:,istate))
			!	if (p_trans_opt(istate,istate) < -1.d0) p_trans_opt(:,istate) = -p_trans_opt(:,istate) / p_trans_opt(istate,istate)
			!enddo
			!!
			!pop_p = pop_p + matmul(p_trans_opt,pop_p)
			!========================================
			!
			psi_l = psi_f
		enddo
		!
		deallocate( psi_l,U_dt,T_trans,p_trans )
		deallocate( WW,p_trans_opt )
		!
	end subroutine final_psi_hop_loc01_dt
	!
	!
	!
	subroutine diag_real(H_in, vec, val)
		!
		implicit none
		!
		real(dp), dimension(:,:) :: H_in, vec
		real(dp), dimension(:)   :: val
		!
		! working variables
		real(dp), allocatable :: work(:)
		integer :: lwork, err_msg
		!
		!
		lwork = (3*nstate-1)*3
		!
		allocate( work(lwork) )
		!
		vec = H_in
		call dsyev('V','U',nstate,vec,nstate,val,work,lwork,err_msg)
		!
		deallocate(work)
		!
	end subroutine
	!
	!
	function logm(A_in)
		! A is real-orthonormal matrix
		!
		implicit none
		!
		real(dp), dimension(:,:) :: A_in
		real(dp), allocatable    :: A(:,:), logm(:,:)
		!
		! working variables
		real(dp), allocatable :: wr(:),wi(:),vl(:,:),vr(:,:),work(:),zeros(:)
		complex(dp), allocatable :: vec1(:,:), vec2(:,:), logw(:,:) ! A = vec1*w*vec2^\dagger
		integer :: lwork, err_msg, incr, istate, jstate
		integer, allocatable :: ipiv(:)
		!
		!
		lwork = (4*nstate)*3
		!
		allocate( A(nstate,nstate) )
		allocate( logm(nstate,nstate) )
		!
		allocate( wr(nstate) )
		allocate( wi(nstate) )
		allocate( vl(nstate,nstate) )
		allocate( vr(nstate,nstate) )
		allocate( work(lwork) )
		!
		A = A_in
		call dgeev('N','V',nstate,A,nstate,wr,wi,vl,nstate,vr,nstate,work,lwork,err_msg)
		!
		allocate( zeros(nstate) )
		allocate( logw(nstate,nstate) )
		allocate( vec1(nstate,nstate) )
		allocate( vec2(nstate,nstate) )
		logw = (0.d0, 0.d0)
		incr = 1
		do istate = 1,nstate
			if ( abs(wi(istate)) < 1.d-10 ) then
				vec1(:,istate) = vr(:,istate)
			else
				if (incr == 1) then
				vec1(:,istate)   = cmplx( vr(:,istate), vr(:,istate+1), dp )
				vec1(:,istate+1) = cmplx( vr(:,istate),-vr(:,istate+1), dp )
				endif
				incr = -incr
			endif
			logw(istate,istate) = log( cmplx( wr(istate), wi(istate), dp ) )
		enddo
		!
		vec2 = transpose(conjg(vec1))
		!
		logm = dble(matmul(matmul(vec1,logw),vec2))
		logm = ( logm - transpose(logm) )/2
		!
		deallocate(wr,wi,vl,vr,work,zeros,A,vec1,vec2,logw)
		!
	end function logm
	!
	!
	function expm(A_in)
		! A is anti-Hermitian
		!
		implicit none
		!
		complex(dp), dimension(:,:) :: A_in
		complex(dp), allocatable    :: expm(:,:)
		!
		! working variables
		complex(dp), allocatable    :: vec(:,:), work(:)
		real(dp)   , allocatable    :: val(:), rwork(:)
		integer :: lwork, err_msg, istate
		!
		!
		lwork = (2*nstate-1)*3
		!
		allocate( expm(nstate, nstate) )
		!
		allocate( vec(nstate, nstate) )
		allocate( val(nstate) )
		allocate( work(lwork) )
		allocate( rwork(3*nstate-2) )
		!
		vec = A_in * (0.d0, -1.d0)
		call zheev('V','U',nstate,vec,nstate,val,work,lwork,rwork,err_msg)
		!
		do istate = 1,nstate
			expm(:,istate) = vec(:,istate) * exp( cmplx(0.d0, val(istate), dp) )
		enddo
		expm = matmul( expm, transpose(conjg(vec)) )
		!
		deallocate(vec,val,work,rwork)
		!
	end function expm
	!
	!
	function sqrtm(A_in)
		! A is unitary matrix
		!
		implicit none
		!
		complex(dp), dimension(:,:) :: A_in
		complex(dp), allocatable    :: A(:,:), sqrtm(:,:)
		!
		! working variables
		complex(dp), allocatable :: ww(:),vl(:,:),vr(:,:),work(:)
		real(dp)   , allocatable :: rwork(:)
		complex(dp), allocatable :: sqrtw(:,:) ! A = vec1*w*vec2^\dagger
		integer :: lwork, err_msg, istate
		!
		!
		lwork = (2*nstate)*3
		!
		allocate( A(nstate,nstate) )
		allocate( sqrtm(nstate,nstate) )
		!
		allocate( ww(nstate) )
		allocate( vl(nstate,nstate) )
		allocate( vr(nstate,nstate) )
		allocate( work(lwork) )
		allocate( rwork(2*nstate) )
		!
		A = A_in
		call zgeev('N','V',nstate,A,nstate,ww,vl,nstate,vr,nstate,work,lwork,rwork,err_msg)
		!
		allocate( sqrtw(nstate,nstate) )
		sqrtw = ( 0.d0, 0.d0 )
		do istate = 1,nstate
			sqrtw(istate,istate) = sqrt(ww(istate))
		enddo
		!
		sqrtm = matmul( matmul(vr, sqrtw), transpose(conjg(vr)) )
		!
		deallocate(ww,vl,vr,work,rwork,A,sqrtw)
		!
	end function sqrtm
	!
	!
	function det_U(A_in)
		! A is real-orthonormal matrix
		!
		implicit none
		!
		real(dp), dimension(:,:) :: A_in
		real(dp), allocatable    :: A(:,:)
		real(dp)                 :: det_U
		!
		! working variables
		real(dp), allocatable :: wr(:),wi(:),vl(:,:),vr(:,:),work(:)
		integer :: lwork, err_msg, istate
		!
		!
		lwork = (4*nstate)*3
		!
		allocate( A(nstate,nstate) )
		!
		allocate( wr(nstate) )
		allocate( wi(nstate) )
		allocate( vl(nstate,nstate) )
		allocate( vr(nstate,nstate) )
		allocate( work(lwork) )
		!
		A = A_in
		call dgeev('N','V',nstate,A,nstate,wr,wi,vl,nstate,vr,nstate,work,lwork,err_msg)
		!
		det_U = 1.d0
		do istate = 1,nstate
			if ( abs(wi(istate)) < 1.d-10 ) det_U = det_U * wr(istate)
		enddo
		!
		deallocate(A,wr,wi,vl,vr,work)
		!
	end function det_U
	!
	!
	function TrlogU2(A_in)
		! A is real-orthonormal matrix
		!
		implicit none
		!
		real(dp), dimension(:,:) :: A_in
		real(dp), allocatable    :: A(:,:)
		real(dp)                 :: TrlogU2
		!
		! working variables
		real(dp), allocatable :: wr(:),wi(:),vl(:,:),vr(:,:),work(:)
		complex(dp)           :: wri
		integer :: lwork, err_msg, istate
		integer, allocatable :: ipiv(:)
		!
		!
		lwork = (4*nstate)*3
		!
		allocate( A(nstate,nstate) )
		!
		allocate( wr(nstate) )
		allocate( wi(nstate) )
		allocate( vl(nstate,nstate) )
		allocate( vr(nstate,nstate) )
		allocate( work(lwork) )
		!
		A = A_in
		call dgeev('N','V',nstate,A,nstate,wr,wi,vl,nstate,vr,nstate,work,lwork,err_msg)
		!
		TrlogU2 = 0.d0
		do istate = 1,nstate
			wri = cmplx( wr(istate), wi(istate), dp )
			TrlogU2 = TrlogU2 + abs(log(wri))**2
		enddo
		!
		deallocate(wr,wi,vl,vr,work,A)
		!
	end function TrlogU2
	!
	!
	recursive subroutine ZY_check_depth(max_ind, istate, jstate, numdepth, ind_viewed)
		implicit none
		!
		integer, dimension(:)     :: max_ind, ind_viewed
		integer                   :: istate, jstate, numdepth
		!
		ind_viewed(istate) = 1
		numdepth = numdepth + 1
		if ( max_ind(istate) .ne. jstate ) call ZY_check_depth(max_ind,max_ind(istate),jstate,numdepth,ind_viewed)
		!
	end subroutine
	!
	subroutine ZY_correct_sign(S, U)
		implicit none
		!
		real(dp), dimension(:,:)  :: S, U
		!
		integer , allocatable     :: max_ind(:), ind_viewed(:), max_ind1(:), max_ind2(:)
		integer                   :: istate, jstate, numdepth, t_ind
		real(dp), allocatable     :: absS(:,:), dmax(:)
		!
		allocate( max_ind   (nstate) )
		allocate( max_ind1  (nstate) )
		allocate( max_ind2  (nstate) )
		allocate( ind_viewed(nstate) )
		allocate( dmax      (nstate) )
		allocate( absS(nstate,nstate))
		!
		absS = abs(S)
		!write(*,'(8(1X,ES11.3))') absS
		do jstate = 1,nstate
			do istate = 1, nstate
				max_ind1(istate) = maxloc( absS(:,istate), 1)
				dmax(istate) = absS(max_ind1(istate),istate)
				absS(max_ind1(istate),istate) = -1.d0
				max_ind2(istate) = maxloc( absS(:,istate), 1)
				absS(max_ind1(istate),istate) = dmax(istate)
				dmax(istate) = dmax(istate) - absS(max_ind2(istate),istate)
			enddo
			t_ind = maxloc( dmax, 1)
			!write(*,'(8(1X,ES11.3))') dmax
			max_ind(t_ind) = max_ind1(t_ind)
			absS(:,t_ind) = 0.d0
			absS(max_ind(t_ind),:) = 0.d0
		enddo
		!write(*,*) max_ind
		!write(*,'(8(1X,ES11.3))') transpose(abs(S))
		!
		!absS = abs(S)
		!do istate = 1, nstate
		!	max_ind(istate) = maxloc( absS(:,istate), 1)
		!enddo
		!
		ind_viewed = 0
		do istate = 1, nstate
			if ( ind_viewed(istate) .eq. 0 ) then
				if ( max_ind(istate) .eq. istate ) then
					ind_viewed(istate) = 1
				else
					numdepth = 1
					call ZY_check_depth(max_ind,max_ind(istate),istate,numdepth,ind_viewed)
					if ( mod(numdepth,2) .eq. 0 ) then
						U(:,istate) = U(:,istate) * (-sign(1.d0, S(max_ind(istate),istate)) )
					else
						U(:,istate) = U(:,istate) *   sign(1.d0, S(max_ind(istate),istate))
					endif
				endif
			else
				U(:,istate) = U(:,istate) *   sign(1.d0, S(max_ind(istate),istate))
			endif
		enddo
		!
		deallocate(max_ind,max_ind1,max_ind2, ind_viewed,dmax,absS)
		!
	end subroutine ZY_correct_sign
	!
	!
	subroutine ZY_correct_sign_full(S, U)
		implicit none
		!
		real(dp), dimension(:,:)  :: S, U
		!
		integer                   :: istate, jstate, nochange
		real(dp)                  :: dtr
		!
		dtr = det_U(S)
		if ( dtr < 0.d0 ) then
			S(:,1) = - S(:,1)
			U(:,1) = - U(:,1)
		endif
		!
		nochange = 0
		do while ( nochange .ne. 1 )
			nochange = 1
			do istate = 1, nstate - 1
				do jstate = istate+1, nstate
					dtr = 3*( S(istate,istate)**2 + S(jstate,jstate)**2 ) + 6*S(istate,jstate)*S(jstate,istate) &
					      + 8*( S(istate,istate)+S(jstate,jstate) ) &
					      - 3*( dot_product(S(istate,:),S(:,istate))+dot_product(S(jstate,:),S(:,jstate)) )
					if ( dtr < 0.d0 ) then
						S(:,istate) = -S(:,istate)
						S(:,jstate) = -S(:,jstate)
						U(:,istate) = -U(:,istate)
						U(:,jstate) = -U(:,jstate)
						nochange = 0
					endif
				enddo
			enddo
		enddo
		!
	end subroutine ZY_correct_sign_full
	!
	!
	subroutine correct_sign_bruteforce(S, U)
		implicit none
		!
		real(dp), dimension(:,:)  :: S, U
		real(dp), allocatable     :: S_t(:,:), U_t(:,:), U_best(:,:)
		real(dp)                  :: f_t, f_best
		integer                   :: iter, jter, pow, ndigit
		!
		allocate( S_t(nstate,nstate) )
		allocate( U_t(nstate,nstate) )
		allocate( U_best(nstate,nstate) )
		!
		f_best = 1.d100
		do iter = 0, 2**nstate - 1
			S_t = S
			U_t = U
			jter = iter
			do ndigit = nstate-1,0,-1
				pow = jter / (2**ndigit)
				jter = mod(jter, 2**ndigit)
				S_t(:,ndigit+1) = ( (-1)**pow )*S_t(:,ndigit+1)
				U_t(:,ndigit+1) = ( (-1)**pow )*U_t(:,ndigit+1)
			enddo
			!
			f_t = TrlogU2(S_t)
			if ( f_t < f_best ) then
				f_best = f_t
				U_best = U_t
			endif
		enddo
		U = U_best
		!
		deallocate( S_t,U_t,U_best)
		!
	end subroutine correct_sign_bruteforce
	!
	!
	subroutine rand_sign(U)
		implicit none
		!
		real(dp), dimension(:,:)  :: U
		integer                   :: istate
		!
		do istate = 1, nstate
			if (rand(0) < 0.5d0) U(:,istate) = -U(:,istate)
		enddo
	end subroutine
	!
	!
	subroutine generate_exact_p_increment(rho_l, rho_f, WW, Q0)
		implicit none
		!
		complex(dp), dimension(:,:)  :: rho_l, rho_f
		real(dp)   , dimension(:,:)  :: Q0
		real(dp)   , allocatable     :: WW(:,:)
		!
		real(dp)   , allocatable     :: pop_f(:,:), pop_l(:,:), dpop(:,:), ppop_dl(:,:), ppop_ll(:,:)
		real(dp)   , allocatable     :: P1(:,:), P2(:,:), one(:,:), eye(:,:)
		integer                      :: istate
		!
		if (.not. allocated(WW)) allocate( WW(nstate,nstate) )
		!
		allocate( pop_f(nstate,1) )
		allocate( pop_l(nstate,1) )
		allocate( dpop (nstate,1) )
		allocate( one  (nstate,1) )
		allocate( P1     (nstate,nstate) )
		allocate( P2     (nstate,nstate) )
		allocate( eye    (nstate,nstate) )
		allocate( ppop_dl(nstate,nstate) )
		allocate( ppop_ll(nstate,nstate) )
		!
		one = 1.d0
		eye = 0.d0
		do istate = 1,nstate
			eye(istate,istate) = 1.d0
		enddo
		!
		! get the constained solution without boundaries
		do istate = 1,nstate
			pop_f(istate,1) = dble(rho_f(istate,istate))
			pop_l(istate,1) = dble(rho_l(istate,istate))
		enddo
		dpop = pop_f - pop_l
		ppop_dl = matmul(dpop ,transpose(pop_l)) / dot_product(pop_l(:,1),pop_l(:,1))
		ppop_ll = matmul(pop_l,transpose(pop_l)) / dot_product(pop_l(:,1),pop_l(:,1))
		P1 = eye-matmul(one,transpose(one))/nstate
		P2 = eye-ppop_ll
		WW = matmul(matmul(P1, Q0), P2) + ppop_dl
		!
		! applying boundaries iteratively
		! TODO: need to find a better algorithm
		do istate = 1,nstate
			WW(istate,istate) = WW(istate,istate) + 1.d0
		enddo
		!
		call remove_neg_increment(P1, P2, WW, Q0, 10)
		!
		deallocate(pop_f,pop_l,dpop,one,P1,P2,eye,ppop_dl,ppop_ll)
		!
	end subroutine generate_exact_p_increment
	!
	!
	subroutine remove_neg_increment(P1, P2, WW, TT0, niter)
		implicit none
		!
		real(dp), dimension(:,:)  :: P1, P2, WW, TT0
		integer                   :: iter, niter
		!
		real(dp), allocatable     :: TT(:,:), sigma(:,:), ss(:,:), xx(:,:)
		integer , allocatable     :: ind1(:), ind2(:), ind0(:)
		integer                   :: nnz, nnz0, t1, t2, t3
		!
		integer                   :: err_msg, lwork
		integer, allocatable      :: ipiv(:), work(:)
		!
		real(dp)                  :: err_sum
		!
		allocate( TT(nstate,nstate) )
		TT = TT0
		do t1 = 1, nstate
			TT(t1,t1) = TT(t1,t1) + 1.d0
		enddo
		!
		allocate( ind1(nstate*nstate) )
		allocate( ind2(nstate*nstate) )
		allocate( ind0(nstate*nstate) )
		!
		nnz  = 0
		nnz0 = 0
		ind1 = 0
		ind2 = 0
		ind0 = 0
		!
		! first find which one from the pair of T_ij, T_ji must be zero
		do t1 = 1,nstate
		do t2 = t1,nstate
			if ( TT(t1,t2) < 1.d-14 ) then
				nnz0 = nnz0+1
				ind1(nnz0) = t1
				ind2(nnz0) = t2
				ind0(nnz0) = t1*nstate+t2
			elseif ( t1 .ne. t2 ) then
				nnz0 = nnz0+1
				ind1(nnz0) = t2
				ind2(nnz0) = t1
				ind0(nnz0) = t2*nstate+t1
			endif
		enddo
		enddo
		!
		! start the increment loop for making all element non-zero
		nnz = nnz0
		do iter = 1, niter
			allocate(sigma(nnz,nnz))
			allocate(ss(nnz,1))
			allocate(xx(nnz,1))
			!
			do t1 = 1,nnz
				ss(t1,1) = -WW(ind1(t1),ind2(t1))
				do t2 = 1,nnz
					sigma(t1,t2) = P1(ind1(t1),ind1(t2))*P2(ind2(t2),ind2(t1))
				enddo
			enddo
			!
			lwork = 3*nnz
			allocate(ipiv(nnz))
			allocate(work(lwork))
			call dgetrf(nnz,nnz,sigma,nnz,ipiv,err_msg)
			call dgetri(nnz,sigma,nnz,ipiv,work,lwork,err_msg)
			if (err_msg .ne. 0) then
				!write(*,*) "warning! not enough dof", nnz, nnz0
				deallocate(sigma,ss,xx,ipiv,work)
				exit
			endif
			xx = matmul(sigma,ss)
			!
			do t1 = 1,nnz
				WW = WW + xx(t1,1)*matmul( P1(:,ind1(t1):ind1(t1)), P2(ind2(t1):ind2(t1),:) )
			enddo
			!
			deallocate(sigma,ss,xx,ipiv,work)
			!
			! check if need next iteration
			nnz = nnz0
			err_sum = 0.d0
			do t1 = 1,nstate
			do t2 = 1,nstate
				if (WW(t1,t2) < -1.d-13 ) then
					!
					err_sum = err_sum - WW(t1,t2)
					!
					do t3 = 1, nnz0
						if (ind0(t3) == t1*nstate+t2) exit
					enddo
					if (t3 == nnz0+1) then
						nnz = nnz + 1
						ind1(nnz) = t1
						ind2(nnz) = t2
					endif
				endif
			enddo
			enddo
			!
			if (nnz == nnz0) exit
			if (err_sum < 1.d-9) exit
		enddo
		!
		deallocate(TT,ind1,ind2,ind0)
		!
		!write(*,*) 'iter=',iter
		!do t1 = 1, nstate
		!	write(*,'(6(ES16.8,1X))')  WW(t1,:)
		!enddo
		!write(*,*) '-----'
		!do t1 = 1, nstate
		!	write(*,'(6(ES16.8,1X))')  TT(t1,:)
		!enddo
		!write(*,*) ''
	end subroutine remove_neg_increment
	!
	!
	subroutine generate_approx_p(psi_l, psi_f, WW, Q0)
		implicit none
		!
		complex(dp), dimension(:,:)  :: psi_l, psi_f
		real(dp)   , dimension(:,:)  :: Q0
		real(dp)   , allocatable     :: WW(:,:)
		!
		real(dp)   , allocatable     :: pop_f(:,:), pop_l(:,:), dpop(:,:), ppop_ll(:,:)
		real(dp)   , allocatable     :: P1(:,:), P2(:,:), Pd(:,:), one(:,:), eye(:,:), QQ(:,:)
		real(dp)                     :: alpha, beta, alpha_1, alpha_2, err_sum
		integer                      :: t1,t2, iter, niter
		!
		if (.not. allocated(WW)) allocate( WW(nstate,nstate) )
		!
		allocate( pop_f(nstate,1) )
		allocate( pop_l(nstate,1) )
		allocate( dpop (nstate,1) )
		allocate( one  (nstate,1) )
		allocate( P1     (nstate,nstate) )
		allocate( P2     (nstate,nstate) )
		allocate( Pd     (nstate,nstate) )
		allocate( eye    (nstate,nstate) )
		allocate( QQ     (nstate,nstate) )
		allocate( ppop_ll(nstate,nstate) )
		!
		one = 1.d0
		eye = 0.d0
		do t1 = 1,nstate
			eye(t1,t1) = 1.d0
		enddo
		!
		pop_f = dble(abs(psi_f)**2)
		pop_l = dble(abs(psi_l)**2)
		dpop = pop_f - pop_l
		!
		P1 = eye-matmul(one,transpose(one))/nstate
		Pd = matmul(dpop,transpose(pop_l))
		!
		alpha_1 = 1.d-7
		alpha_2 = 1.d0
		alpha = alpha_1
		niter = 20
		do iter = 1, niter
			beta = 1.d0-alpha
			ppop_ll = matmul(pop_l,transpose(pop_l))
			P2 = eye-(beta/(alpha+beta*dot_product(pop_l(:,1),pop_l(:,1))))*ppop_ll
			!
			call optimize_approx_hop(P1,P2,Pd,Q0,WW,alpha)
			!
			! check if this alpha is good
			QQ = WW + eye
			err_sum = 0.d0
			do t1 = 1,nstate
			do t2 = 1,nstate
				if (QQ(t1,t2) < 0.d0) err_sum = err_sum - QQ(t1,t2)
			enddo
			enddo
			!
			if (err_sum < 1d-14) then
				alpha_2 = alpha
				alpha = (alpha_1+alpha_2)/2.d0
			else
				alpha_1 = alpha
				alpha = (alpha_1+alpha_2)/2.d0
			endif
		enddo
		!
		WW = WW + eye
		!
		deallocate(pop_f,pop_l,dpop,one,P1,P2,Pd,eye,QQ,ppop_ll)
		!
	end subroutine generate_approx_p
	!
	!
	subroutine optimize_approx_hop(P1,P2,Pd,TT,WW,alpha)
		implicit none
		!
		real(dp), dimension(:,:)  :: P1, P2, Pd, TT
		real(dp)                  :: alpha, beta
		real(dp), allocatable     :: WW(:,:)
		!
		real(dp), allocatable     :: sigma(:,:), ss(:,:), xx(:,:)
		integer , allocatable     :: ind1(:), ind2(:)
		real(dp), allocatable     :: eye(:,:), S_mat(:,:), ss_full(:,:)
		integer                   :: nnz, t1, t2
		!
		integer                   :: err_msg, lwork
		integer, allocatable      :: ipiv(:), work(:)
		!
		beta = 1.d0 - alpha
		if (.not. allocated(WW)) allocate(WW(nstate,nstate))
		!
		allocate( ind1(nstate*nstate) )
		allocate( ind2(nstate*nstate) )
		allocate( eye    (nstate, nstate) )
		allocate( S_mat  (nstate, nstate) )
		allocate( ss_full(nstate, nstate) )
		!
		eye = 0.d0
		do t1 = 1,nstate
			eye(t1,t1) = 1.d0
		enddo
		!
		nnz = 0
		do t1 = 1,nstate
		do t2 = 1,nstate
			if (t1 == t2) then
				if (TT(t1,t2) < 1.d-14 - 1.d0) then
					nnz = nnz + 1
					ind1(nnz) = t1
					ind2(nnz) = t2
				endif
			else
				if (TT(t1,t2) < 1.d-14) then
					nnz = nnz + 1
					ind1(nnz) = t1
					ind2(nnz) = t2
				endif
			endif
		enddo
		enddo
		!
		! nnz must be > 0
		S_mat = alpha*matmul(matmul(P1,TT),P2)+beta*matmul(Pd,P2)+alpha*eye
		!
		allocate(sigma(nnz,nnz))
		allocate(ss(nnz,1))
		allocate(xx(nnz,1))
		!
		do t1 = 1,nnz
			ss(t1,1) = -S_mat(ind1(t1),ind2(t1))
			do t2 = 1,nnz
				sigma(t1,t2) = P1(ind1(t1),ind1(t2))*P2(ind2(t2),ind2(t1))
			enddo
		enddo
		!
		lwork = 3*nnz
		allocate(ipiv(nnz))
		allocate(work(lwork))
		call dgetrf(nnz,nnz,sigma,nnz,ipiv,err_msg)
		call dgetri(nnz,sigma,nnz,ipiv,work,lwork,err_msg)
		!
		xx = matmul(sigma,ss)
		!
		ss_full = 0.d0
		do t1 = 1,nnz
			ss_full(ind1(t1),ind2(t1)) = xx(t1,1)
		enddo
		WW = (matmul(matmul(P1,ss_full),P2)+S_mat)/alpha-eye
		!
		deallocate(sigma,ss,xx,ipiv,work,ind1,ind2,eye,S_mat,ss_full)
		!
	end subroutine optimize_approx_hop
	!
	!
	subroutine test_U()
		!
		implicit none
		!
		real(dp), allocatable  :: S(:,:)
		integer                :: idt
		!
		allocate(S(nstate,nstate))
		do idt = 1,nT-1
			S = matmul(transpose(U(:,:,idt)),U(:,:,idt+1))
			!write(*,'(ES16.8,3(ES16.8,1X))') idt*dt, S(1,1), S(2,2), S(3,3)
			if (S(2,2) < 0.991) then
				write(*,'(4(ES16.8,1X))') E(:,idt+0), (idt+0)*dt
				write(*,'(3(ES16.8,1X))') transpose(U(:,:,idt))
				write(*,*) ''
				write(*,'(4(ES16.8,1X))') E(:,idt+1), (idt+1)*dt
				write(*,'(3(ES16.8,1X))') transpose(U(:,:,idt+1))
				write(*,*) ''
				write(*,'(4(ES16.8,1X))') E(:,idt+2), (idt+2)*dt
				write(*,'(3(ES16.8,1X))') transpose(U(:,:,idt+2))
				write(*,*) ''
				write(*,'(4(ES16.8,1X))') E(:,idt+3), (idt+3)*dt
				write(*,'(3(ES16.8,1X))') transpose(U(:,:,idt+3))
				write(*,*) ''
				write(*,*) ''
				write(*,*) ''
			endif
		enddo
		deallocate(S)
	end subroutine
	!
	subroutine test_T(Tv)
		!
		implicit none
		!
		real(dp), dimension(:,:,:)  :: Tv
		integer                :: idt
		!
		do idt = 1,nT-1
			write(*,'(I12,3(ES16.8,1X))') idt, Tv(2,1,idt), Tv(3,1,idt), Tv(3,2,idt)
		enddo
	end subroutine
	!
	subroutine test_E()
		!
		implicit none
		!
		integer                :: idt
		character(len=100)     :: f_sz
		!
		write(f_sz,*) nstate+1
		do idt = 1,nT
			write(101,'('//adjustl(f_sz)//'(ES16.8,1X))') idt*dt, E(:,idt)
		enddo
	end subroutine
	!
	subroutine test_H()
			!
			implicit none
			!
			integer                :: idt, istate
			real(dp), allocatable  :: H_diag(:)
			character(len=100)     :: f_sz
			!
			allocate(H_diag(nstate))
			write(f_sz,*) nstate+1
			do idt = 1,nT
					do istate = 1,nstate
						H_diag(istate) = H(istate,istate,idt)
					enddo
					write(102,'('//adjustl(f_sz)//'(ES16.8,1X))') idt*dt, H_diag
			enddo
			deallocate(H_diag)
	end subroutine
end module benchmark_system
