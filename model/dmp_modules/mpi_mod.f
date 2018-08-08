module mpi

#ifdef MPI
  include "mpif.h"
#else
  integer :: MPI_COMM_WORLD
#endif

end module mpi
