PROGRAM main
IMPLICIT NONE
INTEGER :: n, i
REAL :: alpha, y0, yn, dy
!real, ALLOCATABLE :: bvec(:,:)
real, ALLOCATABLE :: bvec(:)

WRITE(*,*) 'Input (alpha, n, y0, yn)'
READ(*,*) alpha, n, y0, yn
allocate ( bvec(0:n)   )

dy = (alpha-1)*(yn-y0)/(alpha**n-1);

bvec(0)=y0
bvec(n)=yn

do i=1,n-1
  bvec(i)=dy*(alpha**i-1)/(alpha-1)+y0
end do

print*,'It works? Check your vector:', bvec


END PROGRAM main
