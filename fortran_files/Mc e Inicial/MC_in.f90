module ranmod1
    integer,parameter::i4b=selected_int_kind(9)
    integer,parameter::k15=selected_int_kind(15)
    integer,parameter::dp=selected_real_kind(8)
end module ranmod1
   
module ranmod2
    use ranmod1
    real(dp)::U(97),cc,cd,cm
    integer(k15)::i97,j97
end module ranmod2
   
module ranutil
    use ranmod1
    integer(i4b)::rand_inst=0
    integer(i4b)::rand_feedback=1
    contains
    subroutine initrandom(i,i2)
        use ranmod1
        implicit none
        integer(k15),optional,intent(in)::i,i2
        integer(k15)::seed_in,kl,ij
        character(len=10)::fred
        real(dp)::klr
        
        if (present(i)) then
            seed_in = i
        else
            seed_in = -1
        end if
        if (seed_in /= -1) then
            if (present(i2)) then
                kl = i2
                if (i2 > 30081) stop 'i2 > 30081'
            else
                kl = 9373
            end if
            ij = i
        else
            call system_clock(count=ij)
            ij = mod(ij+rand_inst*100,31328)
            call date_and_time(time=fred)
            read(fred,'(e10.3)')klr
            kl = mod(int(klr*1000),30081)
        end if
        call rmarin(ij,kl)
    end subroutine initrandom

    subroutine rmarin(ij,kl)
        use ranmod1
        use ranmod2
        implicit none
        integer(k15),intent(in)::ij,kl
        integer(k15)::i,j,k,l,ii,jj,m
        real(dp)::s,t
        
        if (ij<0 .or. ij>31328 .or. kl<0 .or. kl>30081) then
            stop 'Erro ao iniciar o gerador de números aleatórios'
        end if
        i = mod(ij/177,177)+2
        j = mod(ij,177)+2
        k = mod(kl/169,178)+1
        l = mod(kl,169)
        do ii = 1,97
            s = 0.0d0
            t = 0.5d0
            do jj = 1,24
                m = mod(mod(i*j,179)*k,179)
                i = j
                j = k
                k = m
                l = mod(53*l+1,169)
                if (mod(l*m,64)>=32) then
                    s = s + t
                end if
                t = 0.5d0*t
            end do
            u(ii) = s
        end do
        cc = 362436.0d0/16777216.0d0
        cd = 7654321.0d0/16777216.0d0
        cm = 16777213.0d0/16777216.0d0
        i97 = 97
        j97 = 33
        
        return
    end subroutine rmarin

    function ranmar() result(fn_val)
        use ranmod1
        use ranmod2
        implicit none
        real(dp) :: fn_val
        real(dp) :: uni
        
        uni =u(i97) -u(j97)
        if ( uni < 0.0d0 ) uni = uni + 1.0d0
        u(i97) = uni
        i97 = i97 - 1
        if (i97 == 0) i97 = 97
        j97 = j97 - 1
        if (j97 == 0) j97 = 97
        cc = cc - cd
        if ( cc < 0.0d0 ) cc = cc + cm
        uni = uni - cc
        if ( uni < 0.0d0 ) uni = uni + 1.0d0
        fn_val = uni
        
        return
    end function ranmar

    real(dp) function rand_real(a,b)
        use ranmod1
        implicit none
        real(dp),intent(in)::a,b
        
        rand_real = (b-a)*ranmar() + a
    end function rand_real

    integer(i4b) function rand_int(i,j)
        use ranmod1
        implicit none
        integer(i4b),intent(in)::i,j
        
        rand_int = int(j*ranmar()) + i
    end function rand_int

end module ranutil
 
   !-----------------------------------------------------------------!
   !                          MARSAGLIA                              !
   ! Usage example:                                                  !
   !    use ranutil                                                  !
   !    ...                                                          !
   !    integer :: i                                                 !
   !    real(rk) :: x,g                                              !
   !    ...                                                          !                                           !
   !    call initrandom()                                            !
   !    x=ranmar()                                                   !
   !    x=rand_real(25.,59.)                                         !
   !    i=rand_int(0,100)                                            !
   !    ...                                                          !
   !                                                                 !
   !-----------------------------------------------------------------!
 

module var_inicial
    integer :: rede,contorno
    integer :: nx,ny,Ns,N_viz     !! ny deve ser par !!
    integer :: N_k,N_T,N_svt,N_svk
    real(8) :: a,Lx,Ly,rc,theta
    character(200) :: dir1,dir2
    
    real(8), dimension(:,:), allocatable :: Aij
    real(8), dimension(:), allocatable :: xt,yt,xk,yk
end module var_inicial

module var_annealing
    real(8), parameter :: pi = 4.0d0*atan(1.0d0)
    integer :: seed,N_sin,IniCon
    integer :: Ns,N_viz,N_mc,N_temp,N_skt,N_kago,N_tri,N_mssf
    integer :: ssf_acc, kag_acc, tri_acc, N_k, N_T
    integer, dimension(:), allocatable :: S,Nviz,jviz
    integer, dimension(:), allocatable :: Sk,Vk,Rk,St,Vt,Rt
    real(8), dimension(:), allocatable :: Aij,Bi
    real(8), dimension(:), allocatable :: rx,ry,mx,my
    real(8) :: Ti,Tf,dT,D
    real(8) :: E_tot,Temp,beta
    character(600) :: dir1,dir2,isimula

    integer, parameter :: iunit_conf = 10, iunit_en = 20, iunit_fin = 30
end module var_annealing

subroutine parametros_iniciais
    use var_inicial, only : contorno,nx,ny,a,theta,rc
    implicit none
    character(111) :: nada

    open(1, file="parametros_rede.dat")
    read(1, *) nada
    read(1, *) theta, nada
    read(1, *) Nx, nada
    read(1, *) Ny, nada
    read(1, *) contorno, nada

    close(1)

    a = 1.0d0
    rc = a
    return
end subroutine parametros_iniciais

subroutine unit_cel
    use var_inicial, only : nx,ny,contorno,theta
    implicit none
    ! logical :: direx
    character(20) :: tipo

    tipo = 'Kagome'
    if (mod(ny,2) .ne. 0) ny = ny + 1
    call diretorios(tipo,nx,ny,contorno,theta)
    call kagome

    if (contorno == 1) then
        call Aij_com_contorno_x
    else if (contorno == 2) then
        call Aij_sem_contorno_x
    end if

end subroutine unit_cel

subroutine diretorios(tipo,nx,ny,cont,theta)
    use var_inicial, only : dir2
    implicit none
    real(8), intent(in) :: theta
    character(8),intent(in) :: tipo
    integer,intent(in) :: nx,ny,cont
    character(10) :: BC
    logical :: direx

    if (cont == 1) then
        BC = "BC"
    else
        BC = "NBC"
    end if
    write(dir2,"(A,'_',I2.2,'x',I2.2,'_',A,'_',I4.4,'/')") trim(tipo),nx,ny,trim(BC),int(100.0d0*theta)

    inquire(file=dir2,exist=direx)
    if (direx .eqv. .false.) then
        call system("mkdir " // trim(dir2))
    end if

    return
end subroutine diretorios

subroutine kagome
    use var_inicial, only : nx,ny,Ns,Lx,Ly,theta
    use var_annealing, only : pi,rx,ry,mx,my
    implicit none
    integer :: i,j,k,is
    real(8), dimension(:), allocatable :: x0,y0,ex,ey

    !!Cria a célula unitária de spins na configuração Kagome!!
    Ns = 3*nx*ny
    allocate(rx(Ns),ry(Ns),mx(Ns),my(Ns))
    allocate(x0(3),y0(3),ex(3),ey(3))

    x0(1) = 0.0d0
    y0(1) = 0.0d0

    x0(2) = cos(60.0d0*pi/180.0d0)
    y0(2) = sin(60.0d0*pi/180.0d0)

    x0(3) = 1.0d0
    y0(3) = 0.0d0

    ex(1) = -cos((30.0d0 + theta)*pi/180.0d0)
    ey(1) = -sin((30.0d0 + theta)*pi/180.0d0)

    ex(2) = cos((90.0d0 + theta)*pi/180.0d0) 
    ey(2) = sin((90.0d0 + theta)*pi/180.0d0)

    ex(3) = -cos((150.0d0 + theta)*pi/180.0d0)
    ey(3) = -sin((150.0d0 + theta)*pi/180.0d0)

    is = 0
    do j = 1,ny
        do i = 1,nx
            do k = 1,3
                is = is + 1
                rx(is) = x0(k) + 2.0d0*(i-1) + mod(j+1,2)
                ry(is) = y0(k) + 2.0d0*sin(60.0d0*pi/180.0d0)*(j-1)
                mx(is) = ex(k)
                my(is) = ey(k)
            end do
        end do
    end do

    Lx = maxval(rx) - minval(rx)
    Ly = maxval(ry) - minval(ry) + sin(60.0d0*pi/180.0d0)

    rx = rx - 0.5d0*Lx
    ry = ry - 0.5d0*Ly

    return
end subroutine kagome

subroutine config
    use var_inicial, only : dir2,Ns,N_k,N_T,xK,yK,xT,yT
    use var_annealing, only : rx,ry,mx,my
    implicit none
    integer :: i

    open(10,file=trim(dir2) // "config0.xyz")
    write(10,*) Ns
    write(10,*) ' '
    do i = 1,Ns
        write(10,*) rx(i),ry(i),mx(i),my(i)
    end do
    close(10)

    call flush()
    
    open(10,file=trim(dir2) // "vertice_k.xyz")
    write(10,*) N_k
    write(10,*) ' '
    do i = 1,N_k
        write(10,*) xk(i),yk(i)
    end do
    close(10)

    call flush()
    
    open(10,file=trim(dir2) // "vertice_T.xyz")
    write(10,*) N_T
    write(10,*) ' '
    do i = 1,N_T
        write(10,*) xT(i),yT(i)
    end do
    close(10)

    return
end subroutine config

subroutine Aij_sem_contorno_x
    use var_inicial, only : Ns,N_viz,Lx,Ly,rc,dir2,Aij
    use var_annealing, only : rx,ry,mx,my
    implicit none
    integer :: i,j,itroca
    real(8) :: x,y,dij,D1,D2,A2
    real(8) :: Jnn, Jtroca, D_dip

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    allocate(Aij(Ns,Ns))
    Aij = 0.0d0
    N_viz = 0
    rc = sqrt(Lx**2 + Ly**2)

    itroca = 0
    if (itroca == 0) then
        Jtroca = 0.0d0
        Jnn = 5.7142857142857189d0
    else if (itroca == 1) then
        D_dip = 1.0d0
        Jtroca = 2.0d0*(5.0d0*D_dip - 7.0d0*D_dip/4.0d0)
        Jnn = 0.5d0*Jtroca + 7*D_dip/4.0d0
    end if


    ! Troca !
    do i = 1,Ns
        do j = 1,Ns
            x = rx(i) - rx(j)
            y = ry(i) - ry(j)
            dij = sqrt(x*x + y*y)
            if (dij < 1.2d0 .and. dij > 0.8) then
                Aij(i,j) = Aij(i,j) - Jtroca*(mx(i)*mx(j) + my(i)*my(j))
            end if
        end do
    end do

    ! Dipolar !
    do i = 1,Ns
        do j = 1,Ns
            x = rx(i) - rx(j)
            y = ry(i) - ry(j)
            dij = sqrt(x*x + y*y)
            A2 = 0.0d0
            if ((dij <= rc).and.(dij>0.4d0)) then
                N_viz = N_viz + 1
                dij = 1.0d0/dij
                x = x*dij
                y = y*dij
                D1 = mx(i)*mx(j) + my(i)*my(j)
                D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                A2 = (D1 - D2)*dij**3
                Aij(i,j) = Aij(i,j) + 0.5d0*D_dip*A2
            end if
        end do
    end do

    Aij = Jnn*Aij
 
    N_viz = 0
    do i = 1,Ns
        do j = 1,Ns
            N_viz = N_viz + 1
            write(11,*) j,Aij(i,j)
        end do
        write(12,*) N_viz
    end do

    close(11)
    close(12)

    return
end subroutine Aij_sem_contorno_x

subroutine Aij_com_contorno_x
    use var_inicial, only : Ns,N_viz,Lx,Ly,rc,dir2,Aij
    use var_annealing, only : rx,ry,mx,my
    implicit none
    integer :: i,j,ni,nj,nc,itroca
    real(8) :: x,y,dij,D1,D2,A2
    real(8) :: Jnn, Jtroca, D_dip

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    allocate(Aij(Ns,Ns))
    Aij = 0.0d0
    N_viz = 0
    rc = sqrt(Lx**2 + Ly**2)
    nc = 3

    itroca = 0
    if (itroca == 0) then
        Jtroca = 0.0d0
        Jnn = 5.7142857142857189d0
    else if (itroca == 1) then
        D_dip = 1.0d0
        Jtroca = 2.0d0*(5.0d0*D_dip - 7.0d0*D_dip/4.0d0)
        Jnn = 0.5d0*Jtroca + 7.0d0*D_dip/4.0d0
    end if

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
                    A2 = 0.0d0
                    if ((dij <= rc).and.(dij>0.4d0)) then
                        N_viz = N_viz + 1
                        dij = 1.0d0/dij
                        x = x*dij
                        y = y*dij
                        D1 = mx(i)*mx(j) + my(i)*my(j)
                        D2 = 3.0d0*(mx(i)*x + my(i)*y)*(mx(j)*x + my(j)*y)
                        A2 = (D1 - D2)*dij**3
                        Aij(i,j) = Aij(i,j) + 0.5d0*D_dip*A2
                    end if
                end do
            end do
        end do
    end do

    Aij = Jnn*Aij
 
    N_viz = 0
    do i = 1,Ns
        do j = 1,Ns
            N_viz = N_viz + 1
            write(11,*) j,Aij(i,j)
        end do
        write(12,*) N_viz
    end do

    close(11)
    close(12)
    return
end subroutine Aij_com_contorno_x

subroutine vertice_T
    use var_inicial, only : Ns,nx,ny,N_T,N_svt,xT,yT,Lx,Ly,dir2
    use var_annealing, only : pi,rx,ry,mx,my
    implicit none
    integer :: i,j,k,ni,nj
    real(8) :: vxu,vyu
    real(8) :: x,y,DD,Lij

    N_T = nx*ny
    N_svt = 6*N_T
    allocate(xT(N_T),yT(N_T))
    
    vxu = 1.5d0
    vyu = 0.5d0*tan(60.0d0*pi/180.0d0)

    k = 0
    do j = 1,ny
        do i = 1,nx
            k = k + 1
            xT(k) = vxu + 2.0d0*(i-1) + mod(j+1,2)
            yT(k) = vyu + 2.0d0*sin(60.0d0*pi/180.0d0)*(j-1)
        end do
    end do

    xT = xT - 0.5d0*Lx
    yT = yT - 0.5d0*Ly

    open(1,file=trim(dir2) // "vertice_T.dat")
    do i = 1,N_T
        do ni = -3,3
            do nj = -3,3
                do j = 1,Ns
                    x = xT(i) - rx(j) + real(ni*Lx,8)
                    y = yT(i) - ry(j) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<1.3d0 .and. DD>0.9d0) then
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
                    if (DD<1.3d0 .and. DD>0.9d0) then
                        write(1,*) j
                    end if
                end do
            end do
        end do
    end do
    close(1)

    return
end subroutine vertice_T

subroutine vertice_K
    use var_inicial, only : Ns,nx,ny,N_k,N_svk,xk,yk,Lx,Ly,dir2
    use var_annealing, only : pi,rx,ry,mx,my
    implicit none
    integer :: i,j,k,l,ni,nj
    real(8) :: vxu(2),vyu(2)
    real(8) :: x,y,DD,Lij

    N_K = 2*nx*ny
    N_svk = 3*N_k
    allocate(xk(N_k),yk(N_k))
    
    vxu(1) = 0.5d0
    vyu(1) = 0.5d0*tan(30.0d0*pi/180.0d0)

    vxu(2) = 0.5d0
    vyu(2) = tan(60.0d0*pi/180.0d0) - vyu(1)

    k = 0
    do j = 1,ny
        do i = 1,nx
            do l = 1,2
                k = k + 1
                xk(k) = vxu(l) + 2.0d0*(i-1) + mod(j+1,2)
                yk(k) = vyu(l) + 2.0d0*sin(60.0d0*pi/180.0d0)*(j-1)
            end do
        end do
    end do

    xk = xk - 0.5d0*Lx
    yk = yk - 0.5d0*Ly

    open(1,file=trim(dir2) // "vertice_K.dat")
    do i = 1,N_K
        do ni = -3,3
            do nj = -3,3
                do j = 1,Ns
                    x = xk(i) - rx(j) + real(ni*Lx,8)
                    y = yk(i) - ry(j) + real(nj*Ly,8)
                    DD = sqrt(x**2 + y**2)
                    if (DD<1.0d0 .and. DD>0.1d0) then
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
                    if (DD<1.0d0 .and. DD>0.1d0) then
                        write(1,*) j
                    end if
                end do
            end do
        end do
    end do
    close(1)

    return
end subroutine vertice_K

subroutine output
    use var_inicial, only: Ns,N_viz,N_svt,N_svk,N_k,N_T,dir2 
    implicit none

    call config

    open(13,file=trim(dir2) // 'input.dat')
    write(13,*) 1.0d0, "! Constante Dipolar"
    write(13,*) Ns,"   !Ns"
    write(13,*) N_viz,"   !N_viz"
    write(13,*) 10000,"    !N_mc"
    write(13,*) 30,"    !N_temp"
    write(13,*) 500.0,"    !Temperatura inicial"
    write(13,*) 0.15,"    !Temperatura final"
    write(13,*) N_svt
    write(13,*) 10, "    !N_sin"
    write(13,*) 10,"    !N_kago"
    write(13,*) 10,"    !N_tri"
    write(13,*) 500,"    !N_mssf"
    write(13,*) 1,"   S ; 1->Aleatorio, 2->Ultima config"
    close(13)
    if (N_svt .ne. N_svk) then
        print*, 'Erro em N_K e N_T',N_K,N_T
    end if
end subroutine output

subroutine inicializa_rede
    use var_inicial, only: Aij,xt,yt,xk,yk

    call parametros_iniciais
    call unit_cel
    call vertice_k
    call vertice_T
    call output  

    deallocate(Aij,xt,yt,xk,yk)

    return
end subroutine inicializa_rede

program inicial_v2
    call inicializa_rede
end program inicial_v2

