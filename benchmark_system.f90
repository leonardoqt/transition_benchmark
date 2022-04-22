module benchmark_system
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	real(dp)  :: v0, dt
	integer   :: nT, nstate
	real(dp), allocatable :: U(:,:,:), E(:,:), T_rate(:,:,:)
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
		if ( allocated(U     ) ) deallocate(U     )
		if ( allocated(E     ) ) deallocate(E     )
		if ( allocated(T_rate) ) deallocate(T_rate)
		allocate(      U(nstate,nstate,nT)   )
		allocate(      E(nstate,nT)          )
		allocate( T_rate(nstate,nstate,nT-1) )
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
	end subroutine assign_psi
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
	end subroutine assign_rho
	!
	!
	subroutine evo_npi()
		!
		implicit none
		!
		integer :: idt
		!
		do idt = 1, nT-1
			T_rate(:,:,idt) = logm(matmul(transpose(U(:,:,idt)),U(:,:,idt+1)) )
		enddo
	end subroutine evo_npi
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
		call dgeev('V','V',nstate,A,nstate,wr,wi,vl,nstate,vr,nstate,work,lwork,err_msg)
		!
		allocate( zeros(nstate) )
		allocate( logw(nstate,nstate) )
		allocate( vec1(nstate,nstate) )
		allocate( vec2(nstate,nstate) )
		logw = (0.d0, 0.d0)
		incr = 1
		do istate = 1,nstate
			if ( abs(wi(istate)) < 1e-10 ) then
				vec1(:,istate) = cmplx( vr(:,istate), zeros )
				vec2(:,istate) = cmplx( vl(:,istate), zeros )
			else
				if (incr == 1) then
				vec1(:,istate)   = cmplx( vr(:,istate), vr(:,istate+1) )
				vec1(:,istate+1) = cmplx( vr(:,istate),-vr(:,istate+1) )
				vec2(:,istate)   = cmplx( vl(:,istate),-vl(:,istate+1) ) ! conjugate is applied
				vec2(:,istate+1) = cmplx( vl(:,istate), vl(:,istate+1) ) ! conjugate is applied
				endif
				incr = -incr
			endif
			logw(istate,istate) = cmplx( wr(istate), wi(istate) )
		enddo
		write(*,'(3(ES12.4,1X))') A_in
		write(*,*) ''
		write(*,'(3(ES12.4,1X))') wi
		write(*,*) ''
		write(*,'(3(ES12.4,1X))') dble( matmul(matmul(vec1,logw),transpose(vec2)) )
		write(*,*) ''
		write(*,*) ''
		!
		deallocate(wr,wi,vl,vr,work,zeros,A,vec1,vec2,logw)
		!
	end function logm
	!
end module benchmark_system
