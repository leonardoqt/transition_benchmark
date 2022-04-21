module model_H
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	integer  :: sz1 = 2
	contains
	!
	!
	function H1(x)
		!
		implicit none
		!
		real(dp)  :: x
		real(dp), allocatable :: H1(:,:)
		!
		allocate(H1(sz1,sz1))
		!
		H1(:,:) = 1.d0
		!
	end function H1
	!
end module model_H
