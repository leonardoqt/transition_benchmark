program test
	!
	use benchmark_system,  only: assign_model, U, E
	use model_H,           only: sz1, H1
	use iso_fortran_env,   only: dp=> real64
	!
	implicit none
	!
	!
	call assign_model(1.0d0,-0.5d0,2.5d0,1.0d0,H1,sz1)
	write(*,*) E
end program test
