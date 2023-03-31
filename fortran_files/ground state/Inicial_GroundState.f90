module var_inicial
    integer :: contorno
    integer :: nx,ny,Ns,N_viz     !! ny deve ser par !!
    integer :: N_k,N_T,N_svt,N_svk,N_theta
    integer, dimension(:), allocatable :: Nviz,jviz
    real(8) :: a,Lx,Ly,rc,theta,thetai,thetaf,dtheta
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    real(8), dimension(:), allocatable :: rx,ry,mx,my,Aij
    real(8), dimension(:), allocatable :: xt,yt,xk,yk
end module var_inicial

subroutine parametros_iniciais
    use var_inicial, only : contorno,nx,ny,a,thetai,thetaf,rc,N_theta,dtheta
    implicit none

    print*, "Parametros iniciais para a rede:"
    print*, "--------------------------------"
    print*, "Entre com o tamanho da rede: nx e ny"
    read(*,*) nx,ny
    print*, "Entre com o raio de corte (indicado = 6.0d0):"
    read(*,*) rc
    print*, "--------------------------------"
    print*, "A simulação é com ou sem contorno?"
    print*, "Com contorno = 1; Sem contorno = 2"
    read(*,*) contorno
    print*, "Ângulo inicial, final e número de pontos:"
    read(*,*) thetai,thetaf,N_theta
    a = 1.0d0

    dtheta = (thetaf - thetai)/real(N_theta-1,8)
  
    rc = rc*a
    return
end subroutine parametros_iniciais

subroutine kagome
    use var_inicial, only : a,pi,nx,ny,Ns,rx,ry,mx,my,Lx,Ly,theta
    implicit none
    integer :: i,j,k,is
    real(8), dimension(:), allocatable :: x0,y0,ex,ey
    character(60) :: dados

    !!Cria a célula unitária de spins na configuração Kagome!!

    Ns = 3*nx*ny
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns))
    allocate(x0(3),y0(3),ex(3),ey(3))

    x0(1) = 0.0d0
    y0(1) = 0.0d0

    x0(2) = a*cos(60.0d0*pi/180.0d0)
    y0(2) = a*sin(60.0d0*pi/180.0d0)

    x0(3) = a
    y0(3) = 0.0d0

    ex(1) = -cos((30.0d0 + theta)*pi/180.0d0)
    ey(1) = -sin((30.0d0 + theta)*pi/180.0d0)

    ex(2) = -cos((90.0d0 + theta)*pi/180.0d0) 
    ey(2) = -sin((90.0d0 + theta)*pi/180.0d0)

    ex(3) = -cos((150.0d0 + theta)*pi/180.0d0)
    ey(3) = -sin((150.0d0 + theta)*pi/180.0d0)

    is = 0
    do j = 1,ny
        do i = 1,nx
            do k = 1,3
                is = is + 1
                rx(is) = x0(k) + 2.0d0*a*(i-1) + a*mod(j+1,2)
                ry(is) = y0(k) + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
                mx(is) = ex(k)
                my(is) = ey(k)
            end do
        end do
    end do

    Lx = maxval(rx) - minval(rx)
    Ly = maxval(ry) - minval(ry) + a*sin(60.0d0*pi/180.0d0)

    rx = rx - 0.5d0*Lx
    ry = ry - 0.5d0*Ly

    write(dados,"('Conf_',I5.5,'.xyz')") int(100.0d0*theta)
    open(44,file='Conf/' // trim(dados))
    write(44,*) Ns
    write(44,*) theta
    do i = 1,Ns
        write(44,*) rx(i),ry(i),mx(i),my(i)
    end do
    close(44)

    call vertice_T
    call vertice_K

    deallocate(x0,y0,ex,ey)

    return
end subroutine kagome

subroutine quadrada
    use var_inicial, only : a,nx,ny,Ns,rx,ry,mx,my,Lx,Ly,theta,pi
    implicit none
    integer :: i,j,k,is
    real(8), dimension(:), allocatable :: x0,y0,ex,ey

    Ns = 2*nx*ny
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns))
    allocate(x0(2),y0(2),ex(2),ey(2))

    x0(1) = 0.5d0*a
    y0(1) = 0.0d0

    x0(2) = 0.0d0
    y0(2) = 0.5d0*a

    ex(1) = cos((0.0d0 + theta)*pi/180.0d0)
    ey(1) = sin((0.0d0 + theta)*pi/180.0d0)

    ex(2) = cos((90.0d0 + theta)*pi/180.0d0) 
    ey(2) = sin((90.0d0 + theta)*pi/180.0d0)

    is = 0
    do j = 1,ny
        do i = 1,nx
            do k = 1,2
                is = is + 1
                rx(is) = x0(k) + a*(i-1)
                ry(is) = y0(k) + a*(j-1)
                mx(is) = ex(k)
                my(is) = ey(k)
            end do
        end do
    end do

    Lx = maxval(rx) - minval(rx) + 0.5d0*a
    Ly = maxval(ry) - minval(ry) + 0.5d0*a

    rx = rx - 0.5d0*Lx
    ry = ry - 0.5d0*Ly

    return
end subroutine quadrada

subroutine Aij_sem_contorno_x
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,rc,theta
    implicit none
    integer :: i,j
    real(8) :: x,y,dij,D1,D2,Aijmin,Aijmax,soma,Aij
    character(60) :: A1,A2

    write(A1,"('Aij_',I5.5,'.dat')") int(100.0d0*theta)
    write(A2,"('Nviz_',I5.5,'.dat')") int(100.0d0*theta)

    open(11,file='Aij/' // trim(A1))
    open(12,file='Nviz/' // trim(A2))

    N_viz = 0
    !rc = 6.0d0*a
    Aijmin = 10000000.0d0
    Aijmax = -10000000.0d0
    do i = 1,Ns
        soma = 0.0d0
        do j = 1,Ns
            if (j.ne.i) then
                x = rx(i) - rx(j)
                y = ry(i) - ry(j)
                if ((sqrt(x*x+y*y) <= rc).and.(sqrt(x*x+y*y)>0.1d0)) then
                    N_viz = N_viz + 1
                    dij = 1.0d0/sqrt(x*x + y*y)
                    x = x*dij
                    y = y*dij
                    D1 = mx(i)*mx(j) + my(i)*my(j)
                    D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                    Aij = (D1 - D2)*dij**3
                    soma = soma + Aij
                    write(11,*) j,Aij
                end if
            end if
        end do
        if (soma > Aijmax) Aijmax = soma
        if (soma < Aijmin) Aijmin = soma
        write(12,*) N_viz
    end do

    deallocate(rx,ry,mx,my)

    close(11)
    close(12)

    return
end subroutine Aij_sem_contorno_x

subroutine Aij_com_contorno_x
    use var_inicial, only : Ns,N_viz,rx,ry,mx,my,Lx,Ly,rc,theta
    implicit none
    integer :: i,j,ni,nj,nc,itroca,jjtroca
    real(8) :: x,y,dij,D1,D2,Aijmin,Aijmax,Aij(Ns,Ns)
    real(8) :: Jnn, Jtroca, D_dip,B3
    character(60) :: A1,A2

    write(A1,"('Aij_',I5.5,'.dat')") int(100.0d0*theta)
    write(A2,"('Nviz_',I5.5,'.dat')") int(100.0d0*theta)

    open(11,file='Aij/' // trim(A1))
    open(12,file='Nviz/' // trim(A2))

    N_viz = 0
    nc = 3
    !rc = 6.0d0*a
    Aijmin = 10000000.0d0
    Aijmax = -10000000.0d0

    rc = sqrt(Lx**2 + Ly**2)
    D_dip = 1.0d0
    Jtroca = 2.0d0*(5.0d0*D_dip - 7.0d0*D_dip/4.0d0)
    Jnn = 0.5d0*Jtroca + 7.0d0*D_dip/4.0d0

    open(157, file='fator_imput.dat')

    print*, 'Com ou sem interação de troca? (0 pra sem, 1 pra com) '
    read(157,*) itroca
    if (itroca == 0) then
        Jtroca = 0.0d0
        print*, 'Qual o valor do interação dos primeiros vizinhos? '
        read(157, *) Jnn
    end if

    print*, 'D_dip: ', D_dip
    print*, 'Jtroca: ', Jtroca
    ! print*, 'Jnn: ', Jnn
    Aij = 0.0d0
    ! Troca !
    do i = 1,Ns
        do ni = -1,1
            do nj = -1,1
                do j = 1,Ns
                    x = rx(i) - rx(j) + real(ni*Lx,8)
                    y = ry(i) - ry(j) + real(nj*Ly,8)
                    dij = sqrt(x*x + y*y)
                    if (dij < 1.2d0 .and. dij > 0.8) then
                        Aij(i,j) = Aij(i,j) - Jtroca*(mx(i)*mx(j) + my(i)*my(j))
                    end if
                end do
            end do
        end do
    end do

    ! Dipolar !
    do i = 1,Ns
        do ni = -nc,nc
            do nj = -nc,nc
                do j = 1,Ns
                    x = rx(i) - rx(j) + real(ni*Lx,8)
                    y = ry(i) - ry(j) + real(nj*Ly,8)
                    dij = sqrt(x*x + y*y)
                    B3 = 0.0d0
                    if ((dij <= rc).and.(dij>0.4d0)) then
                        N_viz = N_viz + 1
                        dij = 1.0d0/dij
                        x = x*dij
                        y = y*dij
                        D1 = mx(i)*mx(j) + my(i)*my(j)
                        D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                        B3 = (D1 - D2)*dij**3
                        Aij(i,j) = Aij(i,j) + 0.5d0*D_dip*B3
                    end if
                end do
            end do
        end do
    end do

    print*, 'Calcular fator de correção de Jnn? (0 para não, 1 para sim)'
    read(157,*) jjtroca
    if (itroca == 0) then
        if (jjtroca == 1) then
            x = rx(1) - rx(2)
            y = ry(1) - ry(2)
            dij = sqrt(x**2 + y**2)
            x = x/dij
            y = y/dij
            D1 = mx(1)*mx(2) + my(1)*my(2)
            D2 = 3.0d0*(mx(1)*x + my(1)*y)*(mx(2)*x + my(2)*y)
            Aijmax = 0.5d0*(D1 - D2)/dij**3
            Aijmax = Jnn/Aijmax
            Aij = Aijmax*Aij
            print*, 'Fator de Correção de Jnn :', Aijmax
        else if (jjtroca == 0) then
            print*, 'Entre com o fator de correção de Jnn.'
            read(157,*) Aijmax
            Aij = Aijmax*Aij

            x = rx(1) - rx(2)
            y = ry(1) - ry(2)
            dij = sqrt(x**2 + y**2)
            x = x/dij
            y = y/dij
            D1 = mx(1)*mx(2) + my(1)*my(2)
            D2 = 3.0d0*(mx(1)*x + my(1)*y)*(mx(2)*x + my(2)*y)
            print*, 'Interação vizinho próximo = ', abs(Aijmax*(0.5d0*(D1 - D2)/dij**3))
        end if
    end if
    
    N_viz = 0
    do i = 1,Ns
        do j = 1,Ns
            !if (AIJ(i,j) .ne. 0.0d0) then
            N_viz = N_viz + 1
            write(11,*) j,Aij(i,j)
            !end if
        end do
        write(12,*) N_viz
    end do

    deallocate(rx, ry, mx, my)

    close(11)
    close(12)
    close(157)
    return
end subroutine Aij_com_contorno_x

subroutine vertice_T
    use var_inicial, only : Ns,nx,ny,N_T,N_svt,xT,yT,a,Lx,Ly,pi,rx,ry,mx,my,theta
    implicit none
    integer :: i,j,k,ni,nj
    real(8) :: vxu,vyu
    real(8) :: x,y,DD,Lij
    character(60) :: nome 

    N_T = nx*ny
    N_svt = 6*N_T
    allocate(xT(N_T),yT(N_T))
    
    vxu = 1.5d0*a
    vyu = 0.5d0*a*tan(60.0d0*pi/180.0d0)

    k = 0
    do j = 1,ny
        do i = 1,nx
            k = k + 1
            xT(k) = vxu + 2.0d0*a*(i-1) + a*mod(j+1,2)
            yT(k) = vyu + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
        end do
    end do

    xT = xT - 0.5d0*Lx
    yT = yT - 0.5d0*Ly

    write(nome,"('verticeT_',I5.5,'.dat')") int(100.0d0*theta)
    open(1,file='Vert/' // trim(nome))
    do i = 1,N_T
        do ni = -3,3
            do nj = -3,3
                do j = 1,Ns
                    x = xT(i) - rx(j) + real(ni*Lx,8)
                    y = yT(i) - ry(j) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<1.3d0*a .and. DD>0.9d0*a) then
                        Lij = x*mx(j) + y*my(j)
                        write(1,*) j,int(sign(1.d0,Lij))
                    end if
                end do
            end do
        end do
    end do

    do i = 1,Ns
        do ni = -3,3
            do nj = -3,3
                do j = 1,N_T
                    x = xT(j) - rx(i) + real(ni*Lx,8)
                    y = yT(j) - ry(i) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<1.3d0*a .and. DD>0.9d0*a) then
                        write(1,*) j
                    end if
                end do
            end do
        end do
    end do
    close(1)

    deallocate(xT,yT)

    return

end subroutine vertice_T

subroutine vertice_K
    use var_inicial, only : Ns,nx,ny,N_k,N_svk,xk,yk,a,Lx,Ly,pi,rx,ry,mx,my,theta
    implicit none
    integer :: i,j,k,l,ni,nj
    real(8) :: vxu(2),vyu(2)
    real(8) :: x,y,DD,Lij
    character(60) :: nome 

    N_K = 2*nx*ny
    N_svk = 3*N_k
    allocate(xk(N_k),yk(N_k))
    
    vxu(1) = 0.5d0*a
    vyu(1) = 0.5d0*a*tan(30.0d0*pi/180.0d0)

    vxu(2) = 0.5d0*a
    vyu(2) = a*tan(60.0d0*pi/180.0d0) - vyu(1)

    k = 0
    do j = 1,ny
        do i = 1,nx
            do l = 1,2
                k = k + 1
                xk(k) = vxu(l) + 2.0d0*a*(i-1) + a*mod(j+1,2)
                yk(k) = vyu(l) + 2.0d0*a*sin(60.0d0*pi/180.0d0)*(j-1)
            end do
        end do
    end do

    xk = xk - 0.5d0*Lx
    yk = yk - 0.5d0*Ly

    write(nome,"('verticeK_',I5.5,'.dat')") int(100.0d0*theta)
    open(1,file='Vert/' // trim(nome))

    do i = 1,N_K
        do ni = -3,3
            do nj = -3,3
                do j = 1,Ns
                    x = xk(i) - rx(j) + real(ni*Lx,8)
                    y = yk(i) - ry(j) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<a .and. DD>0.1d0*a) then
                        Lij = x*mx(j) + y*my(j)
                        write(1,*) j,int(sign(1.d0,Lij))
                    end if
                end do
            end do
        end do
    end do

    do i = 1,Ns
        do ni = -3,3
            do nj = -3,3
                do j = 1,N_K
                    x = xk(j) - rx(i) + real(ni*Lx,8)
                    y = yk(j) - ry(i) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<a .and. DD>0.1d0*a) then
                        write(1,*) j
                    end if
                end do
            end do
        end do
    end do
    close(1)

    deallocate(xK,yK)

    return

end subroutine vertice_K

subroutine output1
    use var_inicial
    implicit none
    integer :: i

    call system("mkdir " // trim('Aij/'))
    call system("mkdir " // trim('Nviz/'))
    call system("mkdir " // trim('Conf/'))
    call system("mkdir " // trim('Vert/'))

    if (mod(ny,2) .ne. 0) ny = ny + 1
    theta = 0.0d0
    N_svt = 6*nx*ny

    open(13,file='input.dat')
    write(13,*) N_theta,"   !N_theta"

    do i = 1,N_theta
        theta = thetai + (i-1)*dtheta
        call kagome
        if (contorno == 1) then
            call Aij_com_contorno_x
        else
            call Aij_sem_contorno_x
        end if 
        write(13,*) int(100.0d0*theta)
    end do

    
    write(13,*) 1.0d0, "! Constante Dipolar"
    write(13,*) Ns,"   !Ns"
    write(13,*) N_viz,"   !N_viz"
    write(13,*) 50000,"    !N_mc"
    write(13,*) 20,"    !N_temp"
    write(13,*) 10.0,"    !Temperatura inicial"
    write(13,*) 0.05,"    !Temperatura final"
    write(13,*) N_svt
    write(13,*) 10, "    !N_sin"
    write(13,*) 10,"    !N_kago"
    write(13,*) 10,"    !N_tri"
    write(13,*) 1,"   S ; 0-> S=1, 1->S aleatório"
    

    close(13)

    return
end subroutine output1

program inicial_v2

    call parametros_iniciais
    call output1
    
end program inicial_v2