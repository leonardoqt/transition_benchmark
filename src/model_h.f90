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
	function H1(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x,shift
		integer   :: sz
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
	function H2(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x,shift
		integer   :: sz
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
	function H3(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x,shift
		integer   :: sz
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
	function Hn(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x,shift
		integer   :: sz
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
	!
	function H_sys_single(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x, shift
		integer   :: sz
		real(dp), allocatable :: H_sys_single(:,:)
		!
		real(dp)  :: x_shift
		integer   :: t1
		!
		allocate(H_sys_single(sz,sz))
		!
		x_shift = 0.05555555555d0
		!
		H_sys_single = 0.d0
		H_sys_single(2,2) = - (x-x_shift) - shift
		do t1 = 3, sz
			H_sys_single(t1,t1) = (x-x_shift) + (t1-2)*shift
			H_sys_single(2,t1)  = 1.d-1 / sqrt(sz*1.d0)
			H_sys_single(t1,2)  = 1.d-1 / sqrt(sz*1.d0)
		enddo
		!
	end function H_sys_single
	!
	function H_sys_rotate(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x, shift
		integer   :: sz
		real(dp), allocatable :: H_sys_rotate(:,:)
		!
		real(dp)  :: ang0, dang, pi2, x_shift
		integer   :: t1
		!
		allocate(H_sys_rotate(sz,sz))
		!
		x_shift = 0.05555555555d0
		pi2 = atan(1.d0) * 2
		ang0 = pi2 / sz
		dang = ang0 * 2
		!
		H_sys_rotate(:,:) = 1.d-1 / sqrt(sz*1.d0)
		!
		do t1 = 1, sz
			H_sys_rotate(t1,t1) = tan(pi2 - ang0 - (t1-1)*dang) * (x-x_shift) + (t1-1)*shift
		enddo
		!
	end function H_sys_rotate
	!
	!
	function H_sys_parallel(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x, shift
		integer   :: sz
		real(dp), allocatable :: H_sys_parallel(:,:)
		!
		real(dp)  :: x_shift
		integer   :: t1, sz1
		!
		allocate(H_sys_parallel(sz,sz))
		!
		x_shift = 0.05555555555d0
		sz1 = sz / 2
		!
		H_sys_parallel(:,:) = 1.d-1 / sqrt(sz*1.d0)
		!
		do t1 = 1, sz1
			H_sys_parallel(t1,t1) = (x-x_shift) + (t1-1)*shift
		enddo
		do t1 = sz1+1, sz
			H_sys_parallel(t1,t1) =-(x-x_shift) + (sz-t1)*shift
		enddo
		!
	end function H_sys_parallel
	!
	!
	function H_sys_mix(x,sz,shift)
		!
		implicit none
		!
		real(dp)  :: x, shift
		integer   :: sz
		real(dp), allocatable :: H_sys_mix(:,:)
		!
		real(dp)  :: ang0, dang, pi2, x_shift
		integer   :: t1,szr,szp1
		!
		allocate(H_sys_mix(sz,sz))
		!
		szr = sz / 2
		szp1 = (sz-szr) / 2
		!
		x_shift = 0.05555555555d0
		pi2 = atan(1.d0) * 2
		ang0 = pi2 / szr
		dang = ang0 * 2
		!
		H_sys_mix(:,:) = 1.d-1 / sqrt(sz*1.d0)
		!
		do t1 = 1, szr
			H_sys_mix(t1,t1) = tan(pi2 - ang0 - (t1-1)*dang) * (x-x_shift) + (t1-1)*shift - 0.5d0
		enddo
		!
		do t1 = szr+1, szr+szp1
			H_sys_mix(t1,t1) = (x-x_shift) + (t1-szr)*shift + 0.5d0
		enddo
		do t1 = szr+szp1+1, sz
			H_sys_mix(t1,t1) =-(x-x_shift) + (sz-t1)*shift + 0.5d0
		enddo
		!
	end function H_sys_mix
	!
end module model_H
