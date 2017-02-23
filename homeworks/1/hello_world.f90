program main
  ! Get access to MPI library
  use MPI
  ! No implicitly defined variables
  implicit none
  ! Declare local variables
  integer :: iRank,nProc,comm,ierror
  
  ! Setup the parallel environment
  call MPI_Init(ierror)
  
  ! Rename the MPI global communicator
  comm=MPI_COMM_WORLD
  
  ! Get processor information
  call MPI_Comm_Rank(comm,iRank,ierror)
  call MPI_Comm_Size(comm,nProc,ierror) 
  
  ! Say hi
  print*,'Hello, my processor rank is',iRank

  ! Finalize MPI
  call MPI_FINALIZE(ierror)

end program main
