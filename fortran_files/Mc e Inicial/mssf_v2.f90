!gfortran -fopenmp -o MSSF mssf_v2.f90

module var_mssf
    integer :: Ns, N_sample, m
    real(8), parameter :: pi2 = 4.0d0*atan(1.0d0), eps = 1e-7
    real(8), dimension(:), allocatable :: x, y, mx, my
    real(8), dimension(:,:), allocatable :: qx, qy
    real(8), dimension(:,:,:), allocatable :: mxp, myp
    real(8), dimension(:,:,:), allocatable :: Axk, Ayk, Bxk, Byk
    real(8), dimension(:,:), allocatable :: A2, B2
    real(8), dimension(:,:), allocatable :: Intens

    integer, dimension(:), allocatable :: S

end module var_mssf

subroutine alloca
    use var_mssf
    implicit none
    
    allocate(qx(m,m), qy(m,m))
    allocate(mxp(m,m,Ns), myp(m,m,Ns))
    allocate(Axk(m,m,Ns), Ayk(m,m,Ns), Bxk(m,m,Ns), Byk(m,m,Ns))
    allocate(A2(m,m), B2(m,m))
    allocate(Intens(m,m))

end subroutine alloca

subroutine dealloca
    use var_mssf

    deallocate(x, y, mx, my, S)
    deallocate(qx, qy)
    deallocate(mxp, myp)
    deallocate(Axk, Ayk, Bxk, Byk)
    deallocate(A2, B2)
    deallocate(Intens)

end subroutine dealloca

subroutine gerar_q(m, n, qx, qy)
    use var_mssf, only : pi2
    implicit none
    integer, intent(in) :: m
    real(8), intent(in) :: n
    real(8), dimension(m,m), intent(out) :: qx, qy
    integer :: i, j
    real(8) :: qmin, qmax, dq
    
    qmin = -n*pi2
    qmax =  n*pi2
    dq = (qmax - qmin)/real(m-1, 8)

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j) SHARED(qx, qy, dq, qmin, m)
    do j = 1, m
        do i = 1, m
            qx(i, j) = qmin + (i-1)*dq
            qy(i, j) = qmin + (j-1)*dq
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine gerar_q

subroutine perp_m(n, m, mx, my, qx, qy, mxp, myp)
    use var_mssf, only : eps
    implicit none
    integer, intent(in) :: n, m
    real(8), dimension(n), intent(in) :: mx, my
    real(8), dimension(m,m), intent(in) :: qx, qy
    real(8), dimension(m,m,n), intent(out) :: mxp, myp
    integer :: i, j, k
    real(8) :: qmod, fac

    !$omp parallel do private(i,j,k,qmod,fac) shared(n,m,mx,my,qx,qy,mxp,myp)
    do k = 1, n
        do j = 1, m
            do i = 1, m
                qmod = sqrt(qx(i, j)**2 + qy(i, j)**2)
                if (qmod < eps) then
                    qmod = eps
                end if
                fac = (mx(k)*qx(i, j) + my(k)*qy(i, j))/qmod
                mxp(i, j, k) = mx(k) - fac*qx(i, j)/qmod
                myp(i, j, k) = my(k) - fac*qy(i, j)/qmod
            end do
        end do
    end do
    !$omp end parallel do

end subroutine perp_m

subroutine calc_AB(n, m, qx, qy, mxp, myp, rx, ry, Axk, Ayk, Bxk, Byk)
    implicit none
    integer, intent(in) :: n, m
    real(8), dimension(m,m), intent(in) :: qx, qy
    real(8), dimension(m,m,n), intent(in) :: mxp, myp
    real(8), dimension(n), intent(in) :: rx, ry 
    real(8), dimension(m,m,n), intent(out) :: Axk, Ayk, Bxk, Byk
    integer :: i, j, k
    real(8) :: fac

    !$omp parallel default(none) shared(n, m, qx, qy, mxp, myp, rx, ry, Axk, Ayk, Bxk, Byk) private(i,j,k,fac)
    !$omp do
    do k = 1, n
        do j = 1, m
            do i = 1, m
                fac = qx(i, j)*rx(k) + qy(i, j)*ry(k)
                Axk(i, j, k) = mxp(i, j, k)*cos(fac)
                Ayk(i, j, k) = myp(i, j, k)*cos(fac)

                Bxk(i, j, k) = mxp(i, j, k)*sin(fac)
                Byk(i, j, k) = myp(i, j, k)*sin(fac)
            end do
        end do
    end do
    !$omp end do nowait
    !$omp end parallel

end subroutine calc_AB

subroutine calc_A2B2(n, m, S, Axk, Ayk, Bxk, Byk, A2, B2)
    implicit none
    integer, intent(in) :: n, m
    integer, dimension(n), intent(in) :: S
    real(8), dimension(m,m,n), intent(in) :: Axk, Ayk, Bxk, Byk
    real(8), dimension(m,m), intent(out) :: A2, B2
    integer :: i, j, k
    real(8), dimension(m,m) :: Ax, Ay, Bx, By

    Ax = 0.0d0
    Ay = 0.0d0
    Bx = 0.0d0
    By = 0.0d0

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(n, m, S, Axk, Ayk, Bxk, Byk, Ax, Ay, Bx, By) PRIVATE(k, i, j)
    do k = 1, n
        do j = 1, m
            do i = 1, m
                Ax(i, j) = Ax(i, j) + S(k)*Axk(i, j, k)
                Ay(i, j) = Ay(i, j) + S(k)*Ayk(i, j, k)
                Bx(i, j) = Bx(i, j) + S(k)*Bxk(i, j, k)
                By(i, j) = By(i, j) + S(k)*Byk(i, j, k)
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    A2 = 0.0d0
    B2 = 0.0d0

    do j = 1, m
        do i = 1, m
            A2(i, j) = Ax(i, j)**2 + Ay(i, j)**2
            B2(i, j) = Bx(i, j)**2 + By(i, j)**2
        end do
    end do

end subroutine calc_A2B2

subroutine intensity(n, m, S, Axk, Ayk, Bxk, Byk, Intens)
    use var_mssf, only : N_sample
    implicit none
    integer, intent(in) :: n, m
    integer, dimension(n), intent(inout) :: S
    real(8), dimension(m,m), intent(in) :: Axk, Ayk, Bxk, Byk
    real(8), dimension(m,m), intent(out) :: Intens
    real(8), dimension(m,m) :: A2, B2
    integer :: i, j, k
    real(8) :: porcentagem

    A2 = 0.0d0
    B2 = 0.0d0
    Intens = 0.0d0

    do k = 1, N_sample
        call read_S(n, S)
        call calc_A2B2(n, m, S, Axk, Ayk, Bxk, Byk, A2, B2)
        do j = 1, m
            do i = 1, m
                Intens(i, j) = Intens(i, j) + (A2(i, j) + B2(i, j))
            end do
        end do

        if (mod(k, N_sample/10) == 0) then
            porcentagem = 100.0d0*real(k, 8)/real(N_sample, 8)
            print '(A, F0.2, A)', 'Concluido: ', porcentagem, '%'
        end if
    end do

    Intens = Intens / real(n*N_sample, 8)

end subroutine intensity

subroutine ler_input
    use var_mssf, only : Ns, N_sample, x, y, mx, my, S
    implicit none
    integer :: i

    open(1, file='config0.xyz')
    read(1, *) Ns
    read(1, *)
    allocate(x(Ns), y(Ns), mx(Ns), my(Ns), S(Ns))
    do i = 1, Ns
        read(1, *) x(i), y(i), mx(i), my(i)
    end do
    close(1)
    open(2, file='saida.dat')
    read(2, *) Ns, N_sample
end subroutine ler_input

subroutine read_S(n, S)
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(out) :: S
    integer :: i, j

    S = 0
    do i = 1, n
        read(2, *) j
        S(i) = 2*j - 1
    end do

end subroutine read_S

subroutine output
    use var_mssf, only : m, qx, qy, Intens
    implicit none
    integer :: i, j

    open(3, file='output.dat')
    write(3, *) 'qx   ', 'qy   ', 'Intens'
    do j = 1, m
        do i = 1, m
            write(3, *) qx(i, j), qy(i, j), Intens(i, j)
        end do
    end do
    close(3)
end subroutine output

subroutine inicializa
    use var_mssf

    m = 301

    call ler_input
    call alloca
    call gerar_q(m, 2.7d0, qx, qy)
    call perp_m(Ns, m, mx, my, qx, qy, mxp, myp)
    call calc_AB(Ns, m, qx, qy, mxp, myp, x, y, Axk, Ayk, Bxk, Byk)

end subroutine inicializa

subroutine calc_intens
    use var_mssf, only : Ns, m, S, Axk, Ayk, Bxk, Byk, Intens
    implicit none

    call inicializa
    call intensity(Ns, m, S, Axk, Ayk, Bxk, Byk, Intens)
    call output
    call dealloca

end subroutine calc_intens

program main
    use omp_lib
    implicit none

    call calc_intens

end program main
