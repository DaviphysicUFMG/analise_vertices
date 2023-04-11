
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
 
module var_inicial
    integer :: rede,contorno,itroca,nc
    integer :: nx,ny,Ns,N_viz     !! ny deve ser par !!
    integer :: N_k,N_T,N_svt,N_svk
    integer :: N_mc,N_temp,N_single,N_kago,N_tri
    real(8) :: Ti,Tf
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
    character(600) :: dir1,dir3,isimula

    integer, parameter :: iunit_conf = 10, iunit_en = 20, iunit_fin = 30
end module var_annealing

subroutine parametros_iniciais
    use var_inicial, only : contorno,itroca,nc,nx,ny,a,theta,N_mc,N_temp,Ti,Tf,N_single,N_kago,N_tri
    implicit none
    character(111) :: nada

    open(1, file="parametros_rede.dat")

    read(1, *) nada
    read(1, *) theta, nada
    read(1, *) Nx, nada
    read(1, *) Ny, nada
    read(1, *) contorno, nada
    read(1, *) itroca, nada
    read(1, *) nc, nada
    read(1, *) N_mc, nada
    read(1, *) N_temp, nada
    read(1, *) Ti, nada
    read(1, *) Tf, nada
    read(1, *) N_single, nada
    read(1, *) N_kago, nada
    read(1, *) N_tri, nada

    close(1)

    a = 1.0d0
    return
end subroutine parametros_iniciais

subroutine unit_cel
    use var_inicial, only : nx,ny,contorno,theta
    implicit none
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

    return
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
    use var_inicial, only : Ns,N_viz,Lx,Ly,rc,dir2,Aij,itroca
    use var_annealing, only : rx,ry,mx,my
    implicit none
    integer :: i,j
    real(8) :: x,y,dij,D1,D2,A2
    real(8) :: Jnn, Jtroca, D_dip

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    allocate(Aij(Ns,Ns))
    Aij = 0.0d0
    N_viz = 0
    rc = sqrt(Lx**2 + Ly**2)

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
    use var_inicial, only : Ns,N_viz,Lx,Ly,rc,dir2,Aij,itroca,nc
    use var_annealing, only : rx,ry,mx,my
    implicit none
    integer :: i,j,ni,nj
    real(8) :: x,y,dij,D1,D2,A2
    real(8) :: Jnn, Jtroca, D_dip

    open(11,file=trim(dir2) // 'Aij.dat')
    open(12,file=trim(dir2) // 'Nviz.dat')

    allocate(Aij(Ns,Ns))
    Aij = 0.0d0
    N_viz = 0
    rc = sqrt(Lx**2 + Ly**2)

    if (itroca == 0) then
        D_dip = 1.0d0
        Jtroca = 0.0d0
        Jnn = 5.7142857142857189
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
    use var_inicial, only: Ns,N_viz,N_svt,N_svk,N_k,N_T,dir2,N_mc,N_temp,Ti,Tf,N_single,N_kago,N_tri
    implicit none

    call config

    open(13,file=trim(dir2) // 'input.dat')
    write(13,*) 1.0d0, "! Constante Dipolar"
    write(13,*) Ns,"   !Ns"
    write(13,*) N_viz,"   !N_viz"
    write(13,*) N_mc,"    !N_mc"
    write(13,*) N_temp,"    !N_temp"
    write(13,*) Ti,"    !Temperatura inicial"
    write(13,*) Tf,"    !Temperatura final"
    write(13,*) N_svt
    write(13,*) N_single, "    !N_sin"
    write(13,*) N_kago,"    !N_kago"
    write(13,*) N_tri,"    !N_tri"
    write(13,*) 500,"    !N_mssf"
    write(13,*) 1,"   S ; 1->Aleatorio, 2->Ultima config"
    close(13)
    if (N_svt .ne. N_svk) then
        print*, 'Erro em N_K e N_T',N_K,N_T
    end if

    return 
end subroutine output

subroutine deallocate_var
    use var_inicial
    implicit none
    rede = 0
    contorno = 0
    itroca = 0
    nc = 0
    nx = 0
    ny = 0
    Ns = 0
    N_viz = 0
    N_k = 0
    N_T = 0
    N_svt = 0
    N_svk = 0
    N_mc = 0
    N_temp = 0
    N_single = 0
    N_kago = 0
    N_tri = 0
    Ti = 0
    Tf = 0
    a = 0
    Lx = 0
    Ly = 0
    rc = 0
    theta = 0
    dir1 = ''

    deallocate(Aij,xt,yt,xk,yk)

    return
end subroutine deallocate_var

subroutine inicializa_rede
    use var_annealing, only : N_temp,Ti,Tf,dT,iunit_en

    call parametros_iniciais
    call unit_cel
    call vertice_k
    call vertice_T
    call output  
    call deallocate_var

    call ler_input(1)
    call ler_config(2)
    call ler_Aij(3)
    call ler_Nviz(4)
    call ler_Vertices(5)


    call inicia_Bi
    call diretorios_MC
    call En_save(iunit_en,0)

    dT = -N_temp/log(Tf/Ti)

    return
end subroutine inicializa_rede

!----------------------------------------------------------------------------!

subroutine ler_input(iunit)
    use ranutil, only : initrandom
    use var_inicial, only : dir2
    use var_annealing, only : Ns,N_viz,N_mc,N_temp,Ti,Tf,N_skt,N_kago,N_tri,N_mssf,N_sin,IniCon,D
    implicit none
    integer, intent(in) :: iunit

    open(unit=iunit,file=trim(dir2) // 'input.dat',status='old',action='read')
    read(iunit,*) D
    read(iunit,*) Ns
    read(iunit,*) N_viz
    read(iunit,*) N_mc
    read(iunit,*) N_temp
    read(iunit,*) Ti
    read(iunit,*) Tf
    read(iunit,*) N_skt
    read(iunit,*) N_sin
    read(iunit,*) N_kago
    read(iunit,*) N_tri
    read(iunit,*) N_mssf
    read(iunit,*) IniCon
    close(iunit)

    call initrandom()

    N_mssf = N_mc/N_mssf

    return
end subroutine ler_input
 
subroutine ler_config(iunit)
    use ranutil
    use var_inicial, only : dir2
    use var_annealing, only : Ns,S,IniCon
    implicit none
    integer, intent(in) :: iunit
    integer :: i

    allocate(S(Ns))

    if (IniCon == 0) then      !Configuração Inicial com Tudo S = 1
        S = 1
    else if (IniCon == 1) then !Configuração Inicial Aleatória
        S = 1
        do i = 1,Ns
            if (ranmar()<0.5d0) then
                S(i) = -1
            end if
        end do
    else if (IniCon == 2) then !Configuração Inicial Igual à Final Anterior
        open(unit=iunit,file=trim(dir2) // 'ConFin.dat',status='old',action='read')
        do i = 1,Ns
            read(iunit,*) S(i)
        end do
    end if

    return
end subroutine ler_config

subroutine ler_Aij(iunit)
    use var_inicial, only : dir2
    use var_annealing, only : N_viz,Aij,jviz
    implicit none
    integer, intent(in) :: iunit
    integer :: i

    allocate(Aij(N_viz),jviz(N_viz))
    open(unit=iunit,file=trim(dir2) // 'Aij.dat',status='old',action='read')
    do i = 1,N_viz
        read(iunit,*) jviz(i),Aij(i)
    end do
    close(iunit)

    return
end subroutine ler_Aij

subroutine ler_Nviz(iunit)
    use var_inicial, only : dir2
    use var_annealing, only : Ns,Nviz
    implicit none
    integer, intent(in) :: iunit
    integer :: i

    allocate(Nviz(0:Ns))
    open(unit=iunit,file=trim(dir2) // 'Nviz.dat',status='old',action='read')

    Nviz(0) = 0
    do i = 1,Ns
        read(iunit,*) Nviz(i)
    end do
    close(iunit)

    return
end subroutine ler_Nviz

subroutine ler_Vertices(iunit)
    use var_inicial, only : dir2
    use var_annealing, only : N_skt,Sk,Vk,Rk,St,Vt,Rt
    implicit none
    integer, intent(in) :: iunit
    integer :: i

    !-----------------------------------------------------------------!

    allocate(St(N_skt),Vt(N_skt),Rt(N_skt))
    open(unit=iunit,file=trim(dir2) // "vertice_T.dat",status='old',action='read')
    do i = 1,N_skt
        read(iunit,*) Vt(i),Rt(i)
    end do
    do i = 1,N_skt
        read(iunit,*) St(i)
    end do
    close(iunit)

    !-----------------------------------------------------------------!

    allocate(Sk(N_skt),Vk(N_skt),Rk(N_skt))
    open(unit=iunit,file=trim(dir2) // "vertice_K.dat",status='old',action='read')
    do i = 1,N_skt
        read(iunit,*) Vk(i),Rk(i)
    end do
    do i = 1,N_skt
        read(iunit,*) Sk(i)
    end do
    close(iunit)
    call flush()

    return
end subroutine ler_Vertices

subroutine inicia_Bi
    use var_annealing, only : Ns,Nviz,jviz,Aij,Bi,S,E_tot
    implicit none
    integer :: i,j,k

    allocate(Bi(Ns))

    E_tot = 0.0d0
    Bi(:) = 0.0d0
    do i = 1,Ns
        do k = Nviz(i-1)+1,Nviz(i)
            j = jviz(k)
            Bi(i) = Bi(i) + S(j)*Aij(k)
        end do
        E_tot = E_tot + S(i)*Bi(i)
    end do
    E_tot = 0.5d0*E_tot

    return
end subroutine inicia_Bi

subroutine diretorios_MC
    use var_inicial, only : dir2
    use var_annealing, only : dir1,dir3,isimula
    implicit none
    integer :: i
    logical :: direx

    dir1 = trim(dir2) // 'Resultados/'
    inquire(file=dir1,exist=direx)
    if (direx .eqv. .false.) then
        call system("mkdir " // trim(dir1))
    end if

    direx = .true.
    i = 0
    do while(direx .eqv. .true.)
        i = i + 1
        write(dir3,"('Resultados/Simula_',I4.4,'/')") i
        inquire(file=trim(dir2) // dir3,exist=direx)
    end do
    write(isimula,"(I4.4)") i
    call system("mkdir " // trim(dir2) // trim(dir3))

    return
end subroutine diretorios_MC

subroutine update(i,dE)
    use var_annealing, only : Nviz,jviz,Aij,Bi,S,E_tot
    implicit none
    integer, intent(in) :: i
    real(8), intent(in) :: dE
    integer :: j,k
    real(8) :: dBi

    do k = Nviz(i-1)+1,Nviz(i)
        j = jviz(k)
        dBi = -2.0d0*S(i)*Aij(k)
        Bi(j) = Bi(j) + dBi
    end do

    S(i) = -S(i)
    E_tot = E_tot + dE

    return
end subroutine update

subroutine Monte_Carlo
    use var_annealing, only : N_mc,beta,iunit_conf,iunit_en, &
                                N_sin,N_kago,N_tri,temp,ssf_acc, &
                                kag_acc,tri_acc,Ns,N_K,N_T,N_mssf
    implicit none
    integer :: imc,kmc,tmc,smc
    real(8) :: SSF, kagome_loop, triagular_loop

    beta = 1.0d0/temp
    do imc = 1,N_mc
        do smc = 1,N_sin
            call metropolis
        end do
        do kmc = 1,N_kago
            call worm_k
        end do
        do tmc = 1,N_tri
            call worm_t
        end do
    end do

    call En_save(iunit_en,1)
    call config_S(iunit_conf, 0)

    ssf_acc = 0
    kag_acc = 0
    tri_acc = 0
    N_K = 0
    N_T = 0

    do imc = 1,N_mc
        do smc = 1,N_sin
            call metropolis
        end do
        do kmc = 1,N_kago
            call worm_k
        end do
        do tmc = 1,N_tri
            call worm_t
        end do
        call En_save(iunit_en,2)
        if (mod(imc,N_mssf).eq.0) then
            call config_S(iunit_conf,1)
        end if
    end do

    SSF = real(ssf_acc, 8)/(N_mc*N_sin*Ns)
    kagome_loop = real(kag_acc, 8)/(N_mc*N_kago)
    triagular_loop = real(tri_acc, 8)/(N_mc*N_tri)

    call config_S(iunit_conf,3)
    write(113,*) temp,SSF,kagome_loop,real(N_K,8)/(N_mc*N_kago),triagular_loop,real(N_T,8)/(N_mc*N_tri)
    call flush()

    return
end subroutine Monte_Carlo

subroutine metropolis
    use var_annealing, only : Ns,Bi,S,beta,ssf_acc
    use ranutil, only : ranmar,rand_int
    implicit none
    integer :: i,im
    real(8) :: dE

    do im = 1,Ns
        i = rand_int(1,Ns)
        dE = -2.0d0*S(i)*Bi(i)
        if (ranmar() < exp(-beta*dE)) then
            call update(i,dE)
            ssf_acc = ssf_acc + 1
        end if
    end do

    return
end subroutine metropolis

subroutine config_S(iunit,flag)
    use var_inicial, only : dir2
    use var_annealing, only : Ns,N_mc,rx,ry,mx,my,S,temp,dir3
    implicit none
    integer, intent(in) :: iunit,flag
    integer :: i
    character(60) :: nome

    if (flag == 0) then
        write(nome,"('config_temp_',f8.4,'.dat')") temp
        open(unit=iunit,file=trim(dir2) // trim(dir3) // trim(nome))
        write(iunit,*) Ns, N_mc
    else if (flag == 1) then
        do i = 1,Ns
            write(iunit,*) (S(i)+1)/2
        end do
    else if (flag == 2) then
        write(iunit,*) Ns
        write(iunit,*) ' '
        do i = 1,Ns
            write(iunit,"(4(2x,f10.5))") rx(i),ry(i),S(i)*mx(i),S(i)*my(i)
        end do
    else if (flag == 3) then
        call flush()
        close(iunit)
    end if

    return
end subroutine config_S

subroutine En_save(iunit,flag)
    use var_inicial, only : dir2
    use var_annealing, only : Ns,N_mc,N_temp,temp,E_tot,dir3,S,mx,my,isimula
    implicit none
    integer, intent(in) :: iunit,flag

    if (flag == 0) then
        open(unit=iunit,file=(trim(dir2) // trim(dir3) // 'energia' // trim(isimula) // '.dat'))
        write(iunit,*) Ns,N_mc,N_temp
    else if (flag == 1) then
        write(iunit,*) temp
    else if (flag == 2) then
        write(iunit,*) E_tot,sum(S*mx),sum(S*my)
    else if (flag == 3) then
        call flush()
        close(iunit)
    end if

    return
end subroutine En_save

subroutine Annealing
    use var_inicial, only : dir2
    use var_annealing, only : N_temp,Ti,dT,temp,dir3
    implicit none
    integer :: i_temp

    open(113, file=trim(dir2) // trim(dir3) // trim('temps.dat'))

    write(113, *) 'Temp ', ' SSF ', ' Kagome ',' Kag_try ',' Triangular ',' Tri_try'
    do i_temp = 0,N_temp
        temp = Ti*exp(-i_temp/dT)
        call Monte_Carlo()
        call flush()
    end do
    close(113)

    return
end subroutine Annealing

subroutine worm_k
    use ranutil
    use var_annealing, only : Ns,N_skt,S,Sk,Vk,Rk,N_K
    implicit none
    integer :: i,j
    integer :: iworm, cont
    integer :: v_0,ivk,isk
    integer :: v_worm(N_skt), s_worm(Ns)
    integer :: v_sequ(N_skt), s_sequ(N_skt)
    integer :: pool(3)

    v_worm = 0
    s_worm = 0
    v_sequ = 0
    s_sequ = 0

    v_0 = int(N_skt*ranmar()/3.0d0)+1
    ivk = v_0
    v_worm(ivk) = 1
    v_sequ(1) = ivk
    iworm = 0

    do while (iworm <= N_skt)
        iworm = iworm + 1
        cont = 0
        pool = 0

        do i = 3*(ivk-1)+1,3*(ivk-1)+3
            j = Vk(i)
            if (Rk(i)*S(j) < 0.0d0) then
                cont = cont + 1
                pool(cont) = j
            end if
        end do

        if (cont /= 0) then
            isk = pool(int(cont*ranmar())+1)
            s_worm(isk) = 1
            s_sequ(iworm) = isk
        else
            exit
        end if

        if (Sk(2*(isk-1)+1) /= ivk) then
            ivk = Sk(2*(isk-1)+1)
        else
            ivk = Sk(2*(isk-1)+2)
        end if

        if (v_worm(ivk) == 1) then
            if (ivk == v_0) then
                N_K = N_k + 1
                call metropolis_loop(s_worm, 0)
                exit
            else
                N_K = N_k + 1
                call corta_rabo(N_skt,ivk,s_worm,s_sequ,v_sequ)
                call metropolis_loop(s_worm, 0)
                exit
            end if
        else
            v_worm(ivk) = 1
            v_sequ(iworm+1) = ivk
        end if
    end do

    return
end subroutine worm_k

subroutine worm_T
    use ranutil
    use var_annealing, only : Ns,N_skt,S,St,Vt,Rt,N_T
    implicit none
    integer :: i,j
    integer :: iworm, cont
    integer :: v_0,ivk,isk
    integer :: v_worm(N_skt), s_worm(Ns)
    integer :: v_sequ(N_skt), s_sequ(N_skt)
    integer :: pool(6)


    v_worm = 0
    s_worm = 0
    v_sequ = 0
    s_sequ = 0

    v_0 = int(N_skt*ranmar()/6.0d0)+1
    ivk = v_0
    v_worm(ivk) = 1
    v_sequ(1) = ivk
    iworm = 0

    do while (iworm <= N_skt)
        iworm = iworm + 1
        cont = 0
        pool = 0

        do i = 6*(ivk-1)+1,6*(ivk-1)+6
            j = Vt(i)
            if (Rt(i)*S(j) < 0.0d0) then
                cont = cont + 1
                pool(cont) = j
            end if
        end do

        if (cont /= 0) then
            isk = pool(int(cont*ranmar())+1)
            s_worm(isk) = 1
            s_sequ(iworm) = isk
        else
            exit
        end if

        if (St(2*(isk-1)+1) /= ivk) then
            ivk = St(2*(isk-1)+1)
        else
            ivk = St(2*(isk-1)+2)
        end if

        if (v_worm(ivk) == 1) then
            if (ivk == v_0) then
                N_T = N_T + 1
                call metropolis_loop(s_worm, 1)
                exit
            else
                N_T = N_T + 1
                call corta_rabo(N_skt,ivk,s_worm,s_sequ,v_sequ)
                call metropolis_loop(s_worm, 1)
                exit
            end if
        else
            v_worm(ivk) = 1
            v_sequ(iworm+1) = ivk
        end if
    end do

    return
end subroutine worm_T

subroutine corta_rabo(N,ivk,s_worm,s_sequ,v_sequ)
    use var_annealing, only : Ns
    implicit none
    integer, intent(in) :: N,ivk
    integer, dimension(N), intent(in) :: s_sequ,v_sequ
    integer, dimension(Ns), intent(inout) :: s_worm
    integer :: i,j

    do i = 1,N
        if (v_sequ(i) /= ivk) then
            j = s_sequ(i)
            s_worm(j) = 0
        else if (v_sequ(i) == ivk) then
            return
        end if
    end do

    return
end subroutine corta_rabo

subroutine metropolis_loop(s_worm, tipo)
    use ranutil
    use var_annealing, only : Ns,S,Bi,beta,E_tot,Aij,Nviz,jviz,kag_acc,tri_acc
    implicit none
    integer , intent(in) :: tipo
    integer, dimension(Ns), intent(in) :: s_worm
    integer, dimension(Ns) :: Sn
    integer :: i,j,k
    real(8), dimension(Ns) :: Bi_n
    real(8) :: E0,En,dE

    Sn = S
    do i = 1,Ns
        if (s_worm(i)==1) then
            Sn(i) = -Sn(i)
        end if
    end do

    E0 = E_tot
    En = 0.0d0

    Bi_n = 0.0d0
    do i = 1,Ns
        do k = Nviz(i-1)+1,Nviz(i)
            j = jviz(k)
            Bi_n(i) = Bi_n(i) + Sn(j)*Aij(k)
        end do
        En = En + Sn(i)*Bi_n(i)
    end do

    En = 0.5d0*En
    dE = En - E0

    if (ranmar() < exp(-beta*dE)) then
        S = Sn
        E_tot = En
        Bi = Bi_n
        if (tipo == 0) then
            kag_acc = kag_acc + 1
        else if (tipo == 1) then
            tri_acc = tri_acc + 1
        end if
    end if

    return
end subroutine metropolis_loop
 
program main

    call inicializa_rede
    call Annealing
    
end program main
