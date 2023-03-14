program main
	integer :: i, Ns
    integer, dimension(:), allocatable :: S, index
    real(8), dimension(:), allocatable :: x, y, mx, my
    
    open(10, file='config0.xyz')
    read(10,*) Ns
    allocate(x(Ns), y(Ns), mx(Ns), my(Ns))
    do i = 1,Ns
		read(10,*) x(i), y(i), mx(i), my(i)
	end do
	close(10)
	
	open(11, file='config0.csv')
	write(11,*) 'x ', 'y ', 'mx ', 'my '
	do i = 1,Ns
		write(11,*) x(i), y(i), mx(i), my(i)
	end do
	close(11)
end program main
	 
