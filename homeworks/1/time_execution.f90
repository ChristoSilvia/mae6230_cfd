program main
  ! Get access to MPI library
  use MPI
  ! No implicitly defined variables
  implicit none
  ! Declare local variables
  integer :: iRank,nProc,comm,ierror,n,i
  real :: a,b,t1,t2
  
  ! Setup the parallel environment
  call MPI_Init(ierror)
  
  ! Rename the MPI global communicator
  comm=MPI_COMM_WORLD
  
  ! Get processor information
  call MPI_Comm_Rank(comm,iRank,ierror)
  call MPI_Comm_Size(comm,nProc,ierror) 
  
  n = 100000000 
  t1 = MPI_WTIME()
  do i=1,n/nProc,1
    call random_number(b) 
    a=exp(5.0+0.001*b)
  cycle
  enddo
  t2 = MPI_WTIME()
  
  ! Say hi
  print*,'Hello, my processor rank is',iRank,'of',nProc,'Processors taking',(t2-t1)

  ! Finalize MPI
  call MPI_FINALIZE(ierror)

end program main
