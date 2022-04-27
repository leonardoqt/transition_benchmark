module benchmark_system
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	real(dp)  :: v0, dt
	integer   :: nT, nstate
	real(dp), allocatable :: U(:,:,:), E(:,:)
	real(dp), allocatable :: psi0(:), rho0(:,:) ! rho is used in real calculation
	!
	contains
	!
	!
	subroutine assign_model(v_in, x_ini_in, x_fin_in, dt_in, H_diabat, nstate_in)
		!
		implicit none
		!
		real(dp) :: v_in, x_ini_in, x_fin_in, dt_in
		integer  :: nstate_in
		!
		! working variables
		integer  :: idt, istate
		!
		interface 
			function H_diabat(x)
				use iso_fortran_env, only: dp=> real64
				implicit none
				real(dp) :: x
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
		allocate( U(nstate,nstate,nT) )
		allocate( E(nstate,nT)        )
		if ( allocated(psi0) ) deallocate(psi0)
		if ( allocated(rho0) ) deallocate(rho0)
		allocate( psi0(nstate) )
		allocate( rho0(nstate,nstate) )
		!
		! diagonalize Hamiltonian to get eigenvector (U) and eigenvalue(E)
		do idt = 1, nT
			call diag_real( H_diabat(x_ini_in+v0*(idt-1)*dt), U(:,:,idt), E(:,idt) )
		enddo
		!
		! naive parallel transport of U
		do idt = 2, nT
			do istate = 1, nstate
				U(:,istate,idt) = U(:,istate,idt) * sign( 1.d0, dot_product(U(:,istate,idt-1),U(:,istate,idt)) )
			enddo
		enddo
		!
	end subroutine assign_model
	!
	!
	subroutine assign_psi(psi_in)
		!
		implicit none
		!
		real(dp), dimension(:) :: psi_in
		!
		!
		psi0 = psi_in
		rho0(:,1) = psi0
		rho0 = matmul( rho0(:,1:1), transpose(rho0(:,1:1)) )
		!
	end subroutine
	!
	subroutine assign_rho(rho_in)
		!
		implicit none
		!
		real(dp), dimension(:,:) :: rho_in
		!
		!
		rho0 = rho_in
		!
	end subroutine
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
			Ut(:,:,idt+1) = matmul(U_dt, U(:,:,idt))
		enddo
		!
		deallocate(U_dt,iHTdt)
		!
	end subroutine evo_npi
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
			Ut(:,:,idt+1) = matmul(U_dt, U(:,:,idt))
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
			iHTdt = iHTdt - Tv(:,:,idt)*( (0.d0, 5.d-1) * dt )
			U_dt = expm(iHTdt)
			U_dt = matmul(transpose(Tv(:,:,idt)), U_dt)
			Ut(:,:,idt+1) = matmul(U_dt, U(:,:,idt))
		enddo
		!
		!
		! the following uses (19) to the first order in dt
		! TODO: the problem of (19) is how to apply to mixed state
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
	subroutine final_rho_hop(Ut, Tv, rho_f, hop_p, state0)
		! use given Ut, Tv and stored rho0 to calculate final rho and the 
		! total hopping rate from state0 to all states
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		integer                        :: state0
		!
		complex(dp), allocatable       :: rho_f(:,:)
		real(dp)   , allocatable       :: hop_p(:)
		!
		real(dp)                       :: drate
		integer                        :: idt, istate
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		if ( allocated(hop_p) ) deallocate(hop_p)
		allocate( rho_f(nstate,nstate) )
		allocate( hop_p(nstate) )
		!
		hop_p = 0.d0
		do idt = 1, nT-1
			rho_f = matmul( matmul(Ut(:,:,idt+1), rho0), transpose(conjg(Ut(:,:,idt+1))) )
			do istate = 1, nstate
				drate = 2 * dble( Tv(state0,istate,idt)*rho_f(istate,state0) ) / dble( rho_f(state0,state0) ) * dt
				if ( drate > 0.d0 ) hop_p(istate) = hop_p(istate) + drate
			enddo
		enddo
		!
	end subroutine final_rho_hop
	!
	!
	subroutine final_rho_hop_loc19(Ut, Tv, rho_f, hop_p, state0)
		! https://doi.org/10.1016/j.comptc.2019.02.009
		!
		! use given Ut, Tv and stored rho0 to calculate final rho and the 
		! total hopping rate from state0 to all states
		!
		implicit none
		!
		complex(dp), dimension(:,:,:)  :: Ut
		real(dp)   , dimension(:,:,:)  :: Tv
		integer                        :: state0
		!
		complex(dp), allocatable       :: rho_f(:,:)
		real(dp)   , allocatable       :: hop_p(:)
		!
		complex(dp), allocatable       :: rho_last(:,:)
		real(dp)                       :: drate, wk, Sk, dP, xkj
		integer                        :: idt, istate
		!
		if ( allocated(rho_f) ) deallocate(rho_f)
		if ( allocated(hop_p) ) deallocate(hop_p)
		allocate( rho_f(nstate,nstate) )
		allocate( hop_p(nstate) )
		!
		allocate( rho_last(nstate,nstate) )
		!
		hop_p = 0.d0
		rho_last = matmul( matmul(Ut(:,:,1), rho0), transpose(conjg(Ut(:,:,1))) )
		do idt = 1, nT-1
			rho_f = matmul( matmul(Ut(:,:,idt+1), rho0), transpose(conjg(Ut(:,:,idt+1))) )
			wk = dble( rho_last(state0,state0)-rho_f(state0,state0) ) / dble( rho_last(state0,state0) )
			if (wk > 0.d0) then
				!
				! calculate Sk
				Sk = 0.d0
				do istate = 1, nstate
					! The diagonal of Tv is set to zero, Tv*dt is the T in the reference
					dP = dble( rho_f(istate,istate)-rho_last(istate,istate) )
					if ( dP > 0.d0 ) then
						! dt is not necessary here, it is included in wk
						Sk = Sk + abs(Tv(state0,istate,idt))*sqrt(dP)
					endif
				enddo
				!
				! calculate xkj
				do istate = 1, nstate
					dP = dble( rho_f(istate,istate)-rho_last(istate,istate) )
					if ( dP > 0.d0 ) then
						xkj = abs(Tv(state0,istate,idt))*sqrt(dP) / Sk
						hop_p(istate) = hop_p(istate) + xkj * wk
					endif
				enddo
			endif
			!
			rho_last = rho_f
		enddo
		!
		deallocate(rho_last)
		!
	end subroutine final_rho_hop_loc19
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
		complex(dp), allocatable :: cwork(:), vec1(:,:), vec2(:,:), logw(:,:) ! A = vec1*w*vec2^\dagger
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
				vec1(:,istate) = cmplx( vr(:,istate), zeros, dp )
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
		vec2 = vec1
		lwork = nstate * 3
		allocate( cwork(lwork) )
		allocate( ipiv(nstate) )
		call zgetrf(nstate,nstate,vec2,nstate,ipiv,err_msg)
		call zgetri(nstate,vec2,nstate,ipiv,cwork,lwork,err_msg)
		!
		logm = dble(matmul(matmul(vec1,logw),vec2))
		logm = ( logm - transpose(logm) )/2
		!write(*,'(3(ES12.4,1X))') logm
		!
		deallocate(wr,wi,vl,vr,work,zeros,A,vec1,vec2,logw,cwork,ipiv)
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
end module benchmark_system
