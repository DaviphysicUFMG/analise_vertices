
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


module var_annealing
   real(8), parameter :: pi = 4.0d0*atan(1.0d0)
   integer :: seed,N_sin,IniCon,itheta
   integer :: Ns,N_viz,N_mc,N_temp,N_skt,N_kago,N_tri,N_mssf,N_theta
   integer, dimension(:), allocatable :: S,Nviz,jviz,theta, Sgs
   integer, dimension(:), allocatable :: Sk,Vk,Rk,St,Vt,Rt
   real(8), dimension(:), allocatable :: Aij,Bi
   real(8), dimension(:), allocatable :: rx,ry,mx,my
   real(8) :: Ti,Tf,dT,D
   real(8) :: E_tot,Temp,beta, Emin
   character(600) :: dir1,dir2,isimula

   integer, parameter :: iunit_conf = 10, iunit_en = 20, iunit_fin = 30
end module var_annealing

subroutine inicial_novo
   use var_annealing
   implicit none

   call ler_input(1)
   ! dT = (log(Tf) - log(Ti))/real(N_temp-1,8)
   dT = -N_temp/log(Tf/Ti)

   allocate(S(Ns),rx(Ns),ry(Ns),mx(Ns),my(Ns), Sgs(Ns))
   allocate(Aij(N_viz),jviz(N_viz))
   allocate(Nviz(0:Ns))
   allocate(St(N_skt),Vt(N_skt),Rt(N_skt))
   allocate(Sk(N_skt),Vk(N_skt),Rk(N_skt))
   allocate(Bi(Ns))

   call diretorios
   
   open(unit=iunit_en,file=trim(dir2) // trim('Em_groundstate.dat'))
   write(iunit_en,*) "theta ","Em"

   do itheta = 1,N_theta
      call ler_config(2)
      call ler_Aij(3)
      call ler_Nviz(4)
      call ler_Vertices(5)
      ! Iniciar Campos !
      call inicia_Bi
      call MC
      print*, "Terminou theta=",theta(itheta)
   end do

   close(iunit_en)
   
   return
end subroutine inicial_novo

subroutine ler_input(iunit)
   use ranutil, only : initrandom
   use var_annealing, only : Ns,N_viz,N_mc,N_temp,N_theta,Ti,Tf,N_skt,N_kago,N_tri,N_sin,IniCon,D,theta
   implicit none
   integer, intent(in) :: iunit
   integer :: i

   open(unit=iunit,file='input.dat',status='old',action='read')
   read(iunit,*) N_theta
   allocate(theta(1:N_theta))
   do i = 1,N_theta
      read(iunit,*) theta(i)
   end do

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
   read(iunit,*) IniCon
   
   close(iunit)

   call initrandom()

   return
end subroutine ler_input

subroutine ler_config(iunit)
   use ranutil
   use var_annealing, only : Ns,S,rx,ry,mx,my,IniCon,theta,itheta
   implicit none
   integer, intent(in) :: iunit
   integer :: i
   character(60) :: dados


   write(dados,"('Conf_',I5.5,'.xyz')") int(theta(itheta))
   open(unit=iunit,file='Conf/' // trim(dados),status='old',action='read')

   read(iunit,*) i
   read(iunit,*) 
   do i = 1,Ns
      read(iunit,*) rx(i),ry(i),mx(i),my(i)
   end do
   close(iunit)

   if (IniCon == 0) then      !Configuração Inicial com Tudo S = 1
      S = 1
   else if (IniCon == 1) then !Configuração Inicial Aleatória
      S = 1
      do i = 1,Ns
         if (ranmar()<0.5d0) then
            S(i) = -1
         end if
      end do
   end if
   return
end subroutine ler_config

subroutine ler_Aij(iunit)
   use var_annealing, only : N_viz,Aij,jviz,D,theta,itheta
   implicit none
   integer, intent(in) :: iunit
   integer :: i
   character(60) :: A1

   write(A1,"('Aij_',I5.5,'.dat')") int(theta(itheta))
   open(unit=iunit,file='Aij/' // trim(A1),status='old',action='read')

   do i = 1,N_viz
      read(iunit,*) jviz(i),Aij(i)
   end do
   close(iunit)

   Aij = Aij*D

   return
end subroutine ler_Aij

subroutine ler_Nviz(iunit)
   use var_annealing, only : Ns,Nviz,theta,itheta
   implicit none
   integer, intent(in) :: iunit
   integer :: i
   character(60) :: A2

   write(A2,"('Nviz_',I5.5,'.dat')") int(theta(itheta))
   open(unit=iunit,file='Nviz/' // trim(A2),status='old',action='read')

   Nviz(0) = 0
   do i = 1,Ns
      read(iunit,*) Nviz(i)
   end do
   close(iunit)

   return
end subroutine ler_Nviz

subroutine ler_Vertices(iunit)
   use var_annealing, only : N_skt,Sk,Vk,Rk,St,Vt,Rt,theta,itheta
   implicit none
   integer, intent(in) :: iunit
   integer :: i
   character(60) :: nome

   write(nome,"('verticeT_',I5.5,'.dat')") int(theta(itheta))
   open(unit=iunit,file='Vert/' // trim(nome),status='old',action='read')

   !-----------------------------------------------------------------!

   do i = 1,N_skt
      read(iunit,*) Vt(i),Rt(i)
   end do
   do i = 1,N_skt
      read(iunit,*) St(i)
   end do
   close(iunit)

   !-----------------------------------------------------------------!
   
   write(nome,"('verticeK_',I5.5,'.dat')") int(theta(itheta))
   open(unit=iunit,file='Vert/' // trim(nome),status='old',action='read')

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

subroutine diretorios
   use var_annealing
   implicit none
   integer :: i
   logical :: direx

   dir1 = 'Resultados/'
   inquire(file=dir1,exist=direx)
   if (direx .eqv. .false.) then
      call system("mkdir " // trim(dir1))
   end if

   direx = .true.
   i = 0
   do while(direx .eqv. .true.)
      i = i + 1
      write(dir2,"('Resultados/Simula_',I4.4,'/')") i
      inquire(file=dir2,exist=direx)
   end do
   write(isimula,"(I4.4)") i
   call system("mkdir " // trim(dir2))

   return
end subroutine diretorios

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

subroutine metropolis
   use var_annealing, only : Ns,Bi,S,beta
   use ranutil, only : ranmar,rand_int
   implicit none
   integer :: i,im
   real(8) :: dE

   do i = 1,Ns
      !i = rand_int(1,Ns)
      dE = -2.0d0*S(i)*Bi(i)
      if (ranmar() .lt. exp(-beta*dE)) then
         call update(i,dE)
      end if
   end do

   return
end subroutine metropolis

subroutine termaliza
   use var_annealing, only : N_mc,N_sin,N_kago,N_tri
   implicit none
   integer :: imc,smc,kmc,tmc

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

   return
end subroutine termaliza 

subroutine MC_final
   use var_annealing!, only : N_mc,N_sin,N_kago,N_tri,E_tot,iunit_en,theta,itheta, Emin,S,Sgs,Ns,iunit_conf
   implicit none
   integer :: i,imc,smc,kmc,tmc
   integer, parameter :: N_mc2 = 1000000 
   ! real(8) :: data(N_mc2),Emed,Evar

   ! data = 0.0d0
   Emin = 0.0d0
   do imc = 1,N_mc2
      do smc = 1,N_sin
         call metropolis
      end do
      do kmc = 1,N_kago
         call worm_k
      end do
      do tmc = 1,N_tri
         call worm_t
      end do
      if (E_tot < Emin) then
         Emin = E_tot
         Sgs = S
      end if
      ! data(imc) = E_tot
   end do

   ! call medias(data,Emed,Evar)
   ! call Config_Final
   ! write(iunit_en,*) theta(itheta),Emed,Evar
   write(iunit_en, *) theta(itheta), Emin

   write(iunit_conf, *) Ns
   write(iunit_conf, *) Emin, theta(itheta)
   do i = 1,Ns
      write(iunit_conf, *) rx(i),ry(i),Sgs(i)*mx(i),Sgs(i)*my(i),atan2(Sgs(i)*my(i),Sgs(i)*mx(i)),Sgs(i)
   end do

   return
end subroutine MC_final

subroutine medias(data,ave,var)
   use var_annealing, only : N_mc
   implicit none
   real(8), dimension(N_mc), intent(in) :: data
   real(8), intent(out) :: ave,var
   real(8) :: s(N_mc)

   ave = sum(data(:))/real(N_mc,8)
   s(:) = data(:) - ave
   var = dot_product(s,s)
   var = ( var - sum(s)**2/real(N_mc,8))/real(N_mc-1,8)

   return
end subroutine medias

subroutine Config_final
   use var_annealing
   implicit none
   integer :: i

   write(iunit_conf,*) Ns
   write(iunit_conf,*) 1.0d0*theta(itheta)/100.0d0
   do i = 1,Ns
      write(iunit_conf,*) rx(i),ry(i),S(i)*mx(i),S(i)*my(i)
   end do

   return
end subroutine Config_final

subroutine worm_k
   use ranutil
   use var_annealing, only : Ns,N_skt,S,Sk,Vk,Rk
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
   !i=rand_int(0,100)
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

      if (cont .ne. 0) then
         isk = pool(int(cont*ranmar())+1)
         s_worm(isk) = 1
         s_sequ(iworm) = isk
      else
         return
      end if

      if (Sk(2*(isk-1)+1) .ne. ivk) then
         ivk = Sk(2*(isk-1)+1)
      else
         ivk = Sk(2*(isk-1)+2)
      end if

      if (v_worm(ivk) .eq. 1) then
         if (ivk .eq. v_0) then
            call metropolis_loop(s_worm)
            return
         else
            call corta_rabo(N_skt,ivk,s_worm,s_sequ,v_sequ)
            call metropolis_loop(s_worm)
            return
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
   use var_annealing, only : Ns,N_skt,S,St,Vt,Rt
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
   !i=rand_int(0,100)
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

      if (cont .ne. 0) then
         isk = pool(int(cont*ranmar())+1)
         s_worm(isk) = 1
         s_sequ(iworm) = isk
      else
         return
      end if

      if (St(2*(isk-1)+1) .ne. ivk) then
         ivk = St(2*(isk-1)+1)
      else
         ivk = St(2*(isk-1)+2)
      end if

      if (v_worm(ivk) .eq. 1) then
         if (ivk .eq. v_0) then
            call metropolis_loop(s_worm)
            return
         else
            call corta_rabo(N_skt,ivk,s_worm,s_sequ,v_sequ)
            call metropolis_loop(s_worm)
            return
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
      if (v_sequ(i) .ne. ivk) then
         j = s_sequ(i)
         s_worm(j) = 0
      else if (v_sequ(i) .eq. ivk) then
         return
      end if
   end do

   return
end subroutine corta_rabo

subroutine metropolis_loop(s_worm)
   use ranutil
   use var_annealing, only : Ns,S,Bi,beta,E_tot,Aij,Nviz,jviz
   implicit none
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

   if (ranmar() .lt. exp(-beta*dE)) then
      S = Sn
      E_tot = En
      Bi = Bi_n
   end if

   return
end subroutine metropolis_loop

subroutine MC
   use var_annealing, only : N_temp,Ti,dT,temp,beta,theta,itheta,dir2,iunit_conf
   implicit none
   integer :: i_temp
   character(60) :: nome

   do i_temp = 1,N_temp
      ! temp = exp(log(Ti) + (i_temp-1)*dT)
      temp = Ti*exp(-i_temp/dT)
      beta = 1.0d0/temp
      call termaliza
   end do

   print*, "Termalizou para theta=",theta(itheta)

   write(nome,"('ConfGS_',I5.5,'.xyz')") int(theta(itheta))
   open(unit=iunit_conf,file=trim(dir2) // trim(nome))
   call MC_final
   close(iunit_conf)

   return
end subroutine MC

program main
 
   call inicial_novo

end program main


