
   !-----------------------------------------------------------------!
   !                          MARSAGLIA                              !
   ! Usage example:                                                  !
   !    use ranutil                                                  !
   !    ...                                                          !
   !    integer :: i                                                 !
   !    real(rk) :: x,g                                              !
   !    ...                                                          !
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

subroutine bootstrap(N,B,data,bootstrap_mean,bootstrap_std)
    use ranutil
    implicit none
    integer, intent(in) :: N, B
    real(8), dimension(N), intent(in) :: data
    real(8), dimension(B), intent(out) :: bootstrap_mean, bootstrap_std
    integer :: i, j, idx
    real(8), dimension(N) :: bootstrap_sample
    real(8) :: mean_B, std_B

    do i = 1, B
        bootstrap_sample(:) = 0.0d0
        do j = 1, N
            idx = rand_int(1, N)
            bootstrap_sample(j) = data(idx)
        end do

        mean_B = sum(bootstrap_sample)/real(N, 8)
        std_B = sum(bootstrap_sample**2)/real(N, 8) - mean_B**2

        bootstrap_mean(i) = mean_B
        bootstrap_std(i) = std_B
    end do

    return
end subroutine bootstrap

subroutine jackknife(N,data,jackknife_mean,jackknife_std)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N), intent(in) :: data
    real(8), dimension(N), intent(out) :: jackknife_mean, jackknife_std
    integer :: i, j, k, JK
    real(8), dimension(N - 1) :: jackknife_sample
    real(8) :: mean, mean_JK, std_JK, erro

    JK = N - 1
    mean = sum(data) / real(N, 8)
    do i = 1, N
        jackknife_sample(:) = 0.0d0
        k = 0
        do j = 1, i-1
            k = k + 1
            jackknife_sample(k) = data(j)
        end do
        do j = i+1, N
            k = k + 1
            jackknife_sample(k) = data(j)
        end do

        mean_JK = sum(jackknife_sample)/real(JK, 8)
        std_JK = sum(jackknife_sample**2)/real(JK, 8) - mean_JK**2

        jackknife_mean(i) = mean_JK
        jackknife_std(i) = std_JK
    end do

    erro = sqrt(sum(jackknife_mean(:) - mean)**2)/real(JK - 1, 8)
    print '(a,1x,es14.7)', "Erro padrão da reamostragem (JackKnife): ", erro
    return
end subroutine jackknife

subroutine correla(N, data, corr)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N), intent(in) :: data
    real(8), dimension(0:N-1), intent(out) :: corr
    integer :: i, j
    real(8) :: mean, var, t_corr

    mean = sum(data)/real(N, 8)
    var = sum(data**2)/real(N, 8) - mean**2 
    do i = 0, N-1
        t_corr = 0.0d0
        do j = 1, N-i
            t_corr = t_corr + (data(j) - mean)*(data(j+i) - mean)
        end do
        t_corr = t_corr/real(N-i,8)
        if (t_corr < 0.0d0) then
            corr(i) = 0.0d0
        else
            corr(i) = t_corr/var
        end if
    end do

    return
end subroutine correla

subroutine measurement_bootstrap(B,Ns,temp,bootstrap_mean,bootstrap_std,E,dE,Cv,dCv,M,dM,Suc,dSuc)
    implicit none
    integer, intent(in) :: B, Ns
    real(8), intent(in) :: temp
    real(8), dimension(B), intent(in) :: bootstrap_mean, bootstrap_std
    real(8), optional, intent(out) :: E, dE, Cv, dCv, M, dM, Suc, dSuc
    real(8) :: cv1(B)
    real(8) :: suc1(B)

    ! Cálculo Calor Específico
    open(2, file='ecv.dat')
    if (present(E)) then
        E = sum(bootstrap_mean)/real(B*Ns, 8)
        dE = sum(bootstrap_std)/real(B*Ns, 8)

        cv1 = bootstrap_std/real(Ns*temp**2, 8)

        Cv = sum(cv1)/real(B, 8)
        dCv = sum(cv1**2)/real(B, 8) - (sum(cv1)/real(B, 8))**2

    ! Cálculo Susceptibilidade Magnetica
    else if (present(M)) then
        M = sum(bootstrap_mean)/real(B*Ns, 8)
        dM = sum(bootstrap_std)/real(B*Ns, 8)

        suc1 = bootstrap_std/real(Ns*temp, 8)

        Suc = sum(suc1)/real(B, 8)
        dSuc = sum(suc1**2)/real(B, 8) - (sum(suc1)/real(B, 8))**2

    end if

    return
end subroutine measurement_bootstrap

subroutine ler(N, data_E, data_M)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N), intent(out) :: data_E, data_M
    integer :: i
    real(8) :: mx, my

    do i = 1, N
        read(1, *) data_E(i), mx, my
        data_M(i) = sqrt(mx**2 + my**2)
    end do

    return
end subroutine ler

subroutine open(i)
    implicit none
    integer, intent(in) :: i
    character(60) :: nome

    write(nome,"('energia_',I6.6,'.dat')") i
    open(1, file='energias/' // trim(nome))
    
    return
end subroutine open

subroutine medias()
    use ranutil
    implicit none
    integer, parameter :: B = 10000
    integer :: Ns, N, idx_i, idx_f
    integer :: k
    real(8) :: temp
    real(8), dimension(:), allocatable :: data_E, data_M
    real(8), dimension(:), allocatable :: boot_mean, boot_std
    real(8) :: E, dE, Cv, dCv, M, dM, Suc, dSuc
    character(30) :: saida_E, saida_M

    Interface measurement_bootstrap
        subroutine measurement_bootstrap(BB,Nss,tempp,bootstrap_meann,bootstrap_stdd,EE,dEE,Cvv,dCvv,MM,dMM,Succ,dSucc)
            implicit none
            integer, intent(in) :: BB, Nss
            real(8), intent(in) :: tempp
            real(8), dimension(BB), intent(in) :: bootstrap_meann, bootstrap_stdd
            real(8), optional, intent(out) :: EE, dEE, Cvv, dCvv, MM, dMM, Succ, dSucc
        end subroutine measurement_bootstrap
    end Interface measurement_bootstrap

    call initrandom()
    allocate(boot_mean(B), boot_std(B))

    open(2, file='input_medias.dat')
    read(2, *) idx_i, idx_f
    read(2, *) saida_E
    read(2, *) saida_M
    close(2)

    open(3, file=saida_E)
    write(3, *) "Temp  ", "E  ", "dE  ", "Cv  ", "dCv  "
    open(4, file=saida_M)
    write(4, *) "Temp  ", "M  ", "dM  ", "Suc  ", "dSuc  "
    do k = idx_i, idx_f
        call open(k)
        read(1, *) Ns, N
        read(1, *) temp
        allocate(data_E(N), data_M(N))
        call ler(N, data_E, data_M)
        close(1)

        call bootstrap(N, B, data_E, boot_mean, boot_std)
        call measurement_bootstrap(B,Ns,temp,boot_mean,boot_std,E,dE,Cv,dCv)
        write(3, *) temp, E, dE, Cv, dCv

        call bootstrap(N, B, data_M, boot_mean, boot_std)
        call measurement_bootstrap(B,Ns,temp,boot_mean,boot_std,M,dM,Suc,dSuc)
        write(4, *) temp, M, dM, Suc, dSuc
        deallocate(data_E, data_M)
        call flush()
    end do
    
end subroutine medias

program main
    call medias
end program main