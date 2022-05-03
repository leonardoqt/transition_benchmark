module model_H
	!
	use iso_fortran_env, only: dp=> real64
	!
	implicit none
	!
	integer  :: sz1 = 2, sz2 = 3, sz3 = 4, szn = 8
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
		H1(:,:) = 1.d-1
		H1(1,1) = H1(1,1) + x
		H1(2,2) = H1(2,2) - x
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
		H2(:,:) = 1.d-1
		H2(1,1) = x
		H2(2,2) = x + 1.d-5
		H2(3,3) = -x
		!
	end function H2
	!
	function H3(x)
		!
		implicit none
		!
		real(dp)  :: x
		real(dp), allocatable :: H3(:,:)
		!
		allocate(H3(sz3,sz3))
		!
		H3(:,:) = 1.d-1
		H3(1,1) = x
		H3(2,2) = x + 1.d-1
		H3(3,3) = -x
		H3(4,4) = -x + 1.d-1
		!
	end function H3
	!
	function Hn(x)
		!
		implicit none
		!
		real(dp)  :: x
		real(dp), allocatable :: Hn(:,:)
		!
		allocate(Hn(szn,szn))
		!
		Hn(:,:) = 1.d-1
		Hn(1,1) = x
		Hn(2,2) = x + 1.d-5
		Hn(3,3) = x + 2.d-5
		Hn(4,4) = x + 3.d-5
		Hn(5,5) =-x
		Hn(6,6) =-x + 1.d-5
		Hn(7,7) =-x + 2.d-5
		Hn(8,8) =-x + 3.d-5
		!
	end function Hn
	!
end module model_H
