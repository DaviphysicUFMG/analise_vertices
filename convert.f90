program main
    integer :: i, j, Ns, M
    integer, dimension(:), allocatable :: S, index
    real(8), dimension(:), allocatable :: x, y, mx, my

    open(10, file="config_temp_  3.6667.dat")
    read(10,*) Ns, M

    allocate(x(Ns), y(Ns), mx(Ns), my(Ns), S(Ns), index(Ns))
    do i = 1,Ns
        read(10,*) x(i), y(i), mx(i), my(i)
        index(i) = i
    end do
    !read(10,*)
    open(11, file="S.dat")
    write(11,*) index(:)
    do j = 1,M
        do i = 1,Ns
            read(10,*) S(i)
        end do
        write(11,*) S(:)
    end do

    close(10)
    close(11)
end program main
