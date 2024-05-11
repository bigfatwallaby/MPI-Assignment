program MontyParallel
    implicit none
    include "mpif.h"
    
    integer :: seed_size
    integer, allocatable :: seed(:)
    real :: result
    real :: MonteCarloIntegrate
    real :: exact, error, receive, tempSend
    integer :: Nmax, N, i, Ni, counterR, counterS
    integer rank, size, err
    
    write(*,*) "hello world 1"

    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)
    
    write(*,*) "hello world 2"

    call MPI_INIT(err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)

    write(*,*) "hello world 3"

    Nmax = 100000000
    exact = 0.80656718084400884701112678335185691868951443065656
    
    open(unit = 1, file = "error vs N.dat", action = "write")
    
    N = 1
    
    counterS = 0
    counterR = 0
    tempSend = 0

    if (rank .eq. 0) then
        do while(N<Nmax)
            !send out all the tasks
            write(*,*) rank, "sending out tasks"
            do i = 1, (size - 2)
                Ni = floor(real(N)/(size-1))
                !send chunks to ppl
                write(*,*) rank, "sending ", Ni, " to ", i
                call MPI_SEND(Ni, 1, MPI_INTEGER, i, counterS, MPI_COMM_WORLD, err)


            end do
            ! send remainder to last guy
            write(*,*) rank, "sending out ", N-Ni*(size-1), " to ",  (size - 1)
            call MPI_SEND(N-Ni*(size-1), 1, MPI_INTEGER, (size - 1), counterS, MPI_COMM_WORLD, err)

            counterS = counterS + 1

            ! recieve all results
            write(*,*) rank, " waiting to receive all results"
            call MPI_REDUCE(tempSend, result, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD, err)
            write(*,*) rank, " received all results"

            result = result / real(size-1)
            error = ABS(exact - result)/ exact
            write(*,*) rank, " result wrote to file"
            write(1,*) N, error
            N = N * 2
        end do

    else
        !receive task
        write(*,*) rank, " waiting to receive task"
        call MPI_RECV(Ni, 1, MPI_INTEGER, 0, counterR, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
        counterR = counterR + 1
        write(*,*) rank, " received task"

        !perform task
        result = MonteCarloIntegrate(Ni)
        write(*,*) rank, " performed task"

        !send back
        call MPI_REDUCE(result, receive, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD, err)
        write(*,*) rank, " sent back"
    end if


    call MPI_FINALIZE(err) 
    end program MontyParallel
    
    real function rand()
        call random_number(rand)
    end function rand

    real function g(x,y,z)
        implicit none
        real,intent(in) :: x,y,z
        g = exp(-(x+y+z))
    end function g
    
    real function MonteCarloIntegrate(N)
        integer, intent(in) :: N
        real :: randX, randY, randZ, rand
        real :: bigG, g, V
    
        bigG = 0
    
        do i = 1,N
            randX = rand()*2
            randY = rand()*3
            randZ = rand()*4
    
            bigG = bigG + g(randX,randY,randZ)
        end do
        bigG = bigG/N
        V = 2 * 3 * 4
        MonteCarloIntegrate = bigG * V
    end function MonteCarloIntegrate