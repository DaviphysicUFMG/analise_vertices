
module var_model
    integer :: Nx, Ny, Ns, N_vertex, N_ind_vert
    integer, dimension(:), allocatable :: index_vertex, index_spin
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    
    real(8) :: Lx, Ly
    real(8), dimension(:), allocatable :: x, y, mx, my
    real(8), dimension(:), allocatable :: v_x, v_y, vertex_q
    

end module var_model

subroutine calc_vertex()
    use var_model
    implicit none
    integer :: i, j, k, k_i
    real(8) :: vx0(3), vy0(3), xmin, ymin
    

    xmin = minval(x)
    ymin = minval(y)
    Lx = maxval(x) - minval(x)
    Ly = maxval(y) - minval(y) + sin(60.0d0*pi/180.0d0)
    
    vx0(1) = xmin + 0.5d0
    vx0(2) = xmin + 0.5d0
    vx0(3) = xmin + 1.5d0

    vy0(1) = ymin + 0.5d0*tan(30.0d0*pi/180.0d0)
    vy0(2) = ymin + tan(60.0d0*pi/180.0d0) - 0.5d0*tan(30.0d0*pi/180.0d0)
    vy0(3) = ymin + 0.5d0*tan(60.0d0*pi/180.0d0)

    allocate(v_x(N_vertex), v_y(N_vertex))

    k = 0
    do j = 1, Nx
        do i = 1, Ny
            do k_i = 1, 2
                k = k + 1
                v_x(k) = vx0(k_i) + 2.0d0*(i-2) + mod(j+1, 2)
                v_y(k) = vy0(k_i) + 2.0d0*sin(60.0d0*pi/180.0d0)*(j-1)
            end do
            k_i = k_i + 1
            k = k + 1
            v_x(k) = vx0(k_i) + 2.0d0*(i-2) + mod(j+1, 2)
            v_y(k) = vy0(k_i) + 2.0d0*sin(60.0d0*pi/180.0d0)*(j-1)
        end do
    end do

    return
end subroutine calc_vertex

subroutine calc_vertex_q()
    use var_model
    implicit none
    integer :: i, j, k, ni, nj
    real(8) :: dx, dy, dist, lij

    allocate(index_vertex(N_ind_vert), index_spin(N_ind_vert), vertex_q(N_ind_vert))

    k = 0
    do i = 1, N_vertex
        do ni = -1, 1
            do nj = -1, 1
                do j = 1, Ns
                    dx = v_x(i) - x(i) + real(ni*Lx, 8)
                    dy = v_y(i) - y(i) + real(nj*Ly, 8)
                    dist = sqrt(dx**2 + dy**2)
                    if ( mod(i, 3) /= 0) then
                        if (1.1d0 > dist .and. dist > 0.1d0) then
                            k = k + 1
                            dx = dx / dist
                            dy = dy / dist
                            lij = (dx*mx(j) + dy*my(j)) / dist**3
                            index_vertex(k) = i
                            index_spin(k) = j
                            vertex_q(k) = lij
                        end if
                    else
                        if (1.3d0 > dist .and. dist > 0.9d0) then
                            k = k + i 
                            dx = dx / dist
                            dy = dy / dist
                            lij = (dx*mx(j) + dy*my(j)) / dist**3
                            index_vertex(k) = i
                            index_spin(k) = j
                            vertex_q(k) = lij
                        end if
                    end if
                end do
            end do
        end do
    end do


    return
end subroutine calc_vertex_q