program main

  ! No implicitly defined variables
  implicit none
  ! Declare local variables
  REAL :: alpha, dy, yn, y0, rd, bvec(0:128)
  integer :: i, n
  n = 128
  alpha = 1.05
  yn = 1
  y0 = 0
  dy = (alpha-1)*(yn-y0)/(alpha**n-1);
  bvec(0)=y0
  bvec(n)=yn

  do i=1,n-1
    bvec(i)=dy*(alpha**i-1)/(alpha-1)+y0
  end do


end program main
