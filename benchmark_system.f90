module benchmark_system
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	real(dp)  :: v0, dt
	integer   :: nT, nstate
	real(dp), allocatable :: U(:,:,:), E(:,:)
	real(dp), allocatable :: psi(:), rho(:,:) ! rho is used in real calculation
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
		integer  :: idt
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
		if ( allocated(U) ) deallocate(U)
		if ( allocated(E) ) deallocate(E)
		allocate( U(nstate,nstate,nT) )
		allocate( E(nstate,nT) )
		!
		do idt = 1,nT
			call diag_real( H_diabat(x_ini_in+v0*(idt-1)*dt), U(:,:,idt), E(:,idt) )
		enddo
		!
	end subroutine assign_model
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
	!
end module benchmark_system
