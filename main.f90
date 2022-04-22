program test
	!
	use benchmark_system,  only: assign_model, assign_psi, U, E, rho0
	use model_H,           only: sz1, sz2, H1, H2
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	real(dp), allocatable :: psi2(:)
	!
	call assign_model(1.0d0,-0.5d0,2.5d0,1.0d0,H1,sz1)
	call assign_model(1.0d0,-0.5d0,2.5d0,1.0d0,H2,sz2)
	write(*,'(3(ES12.4,1X))') U
	!
	allocate(psi2(sz2))
	psi2 = 1.d0
	call assign_psi(psi2)
	write(*,'(3(ES12.4,1X))') rho0
	!
end program test
