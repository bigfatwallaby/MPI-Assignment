program test

    implicit none
    include "mpif.h"
    integer rank, size, err

    call MPI_INIT(err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
    print*, "Hello world from process ", rank, " of ", size, "!"
    call MPI_FINALIZE(err) 

    write(*,*) "successfully uploaded and written something!"


end program test