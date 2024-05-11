program MontyParallel
    implicit none
    include "mpif.h"
    
    integer :: seed_size
    integer, allocatable :: seed(:)
    real :: result
    real :: MonteCarloIntegrate
    real :: exact, error
    integer :: Nmax, N
    integer rank, size, err
    
    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)
    
    call MPI_INIT(err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)

    Nmax = 100000000
    exact = 0.80656718084400884701112678335185691868951443065656
    
    open(unit = 1, file = "error vs N.dat", action = "write")
    
    N = 1
    
    do while(N<Nmax)
        result = MonteCarloIntegrate(N)
        error = ABS(exact - result)/ exact
        write(1,*) N, error
        !N = ((N+1) * 1.1)
        N = N * 2
    end do

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