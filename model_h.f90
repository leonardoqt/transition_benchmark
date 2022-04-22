module model_H
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	integer  :: sz1 = 2, sz2 = 3
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
	function H2(x)
		!
		implicit none
		!
		real(dp)  :: x
		real(dp), allocatable :: H2(:,:)
		!
		allocate(H2(sz2,sz2))
		!
		H2(:,:) = 1.d-1 * x
		H2(1,1) = H2(1,1) + 1
		H2(2,2) = H2(2,2) + 2
		H2(3,3) = H2(3,3) + 3
		!
	end function H2

end module model_H
