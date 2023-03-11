program main
    implicit none
    integer, parameter :: Ns=300, N_size=10000
    integer, dimension(:), allocatable :: Spin, index
    integer :: i, j

    open(10, file='config_temp_500.0000.dat')
    open(11, file='Spin.dat')

    allocate(Spin(Ns), index(Ns))
    do i = 1,Ns
        index(i) = i
    end do
    write(11,*) index(:)
    do i = 1,N_size
        do j = 1,Ns
            read(10,*) Spin(j)
        end do
        write(11,*) 2*Spin(:)-1
    end do

    close(10)
    close(11)
end program main