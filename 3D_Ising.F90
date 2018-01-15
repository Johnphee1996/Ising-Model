program monte_carlo
implicit none
integer i,j,k,i0,j0,k0,npoin,n1,n2,n3
integer i1,j1,k1,Smax
integer sam1,sam2,sam3
integer M,mediate
integer itemp,ntemp,time
integer,allocatable :: sigma(:,:,:),M_dat(:)
real*8,allocatable :: Tr(:)
real*8 p,H_in,H_se,H_to,delta_H,J_int,B
real*8 kB,judge1,judge2
real*8 T_min,T_max
real*8 M_even
real*8 t1,t2
!
npoin = 100
J_int = 1.d0
B = 0.d0
kB = 1.d0
Smax = 2*npoin**3
mediate = npoin**3
n1 = 50000
n2 = 50000
n3 = 1000
T_min = 0.1d0
T_max = 0.35d0
ntemp = 51
allocate(sigma(npoin+2,npoin+2,npoin+2))
allocate(M_dat(Smax+1))
allocate(Tr(ntemp))
!
do itemp = ntemp, 1, -1
    Tr(itemp) = (T_max - T_min)/(ntemp - 1)*(itemp - 1) + T_min
end do
!
call init_random_seed()
!=============
!Initial value
!=============
do i0 = 1, npoin+2
    do j0 = 1, npoin+2
        do k0 = 1, npoin+2
            !call random_number(judge1)
            !if (judge1 .lt. 0.5d0) then
            !    sigma(i0,j0) = -1
            !else
            !    sigma(i0,j0) = 1
            !end if
            sigma(i0,j0,k0) = 1
        end do
    end do
end do
!=============================
!The main part of this program
!=============================
do itemp = 1, ntemp
!do itemp = ntemp, 1, -1
    call cpu_time(t1)
    !initial value of the magnetic moment
    M_dat = 0
    !do i = 1, Smax+1
    !    M_dat(i) = 0
    !end do
    !n1 time loop to get the equilibrium
    do i1 = 1, n1
        call Randi(2,npoin+1,sam1)
        call Randi(2,npoin+1,sam2)
        call Randi(2,npoin+1,sam3)
        !The interaction Hamiltionian
        H_in = -J_int*sigma(sam1,sam2,sam3)*(sigma(sam1+1,sam2,sam3) + sigma(sam1-1,sam2,sam3) + &
                                           & sigma(sam1,sam2+1,sam3) + sigma(sam1,sam2-1,sam3) + &
                                           & sigma(sam1,sam2,sam3+1) + sigma(sam1,sam2,sam3-1))
        !The Hamiltonian in the magnetic field
        H_se = -B*sigma(sam1,sam2,sam3)
        !The total Hamiltonian
        H_to = H_in + H_se 
        !---------------------------------------------------------------------------------
        !If H_to is positive, the spin orientation needs to be flipped to lower the energy
        !---------------------------------------------------------------------------------
        if (H_to .gt. 0.d0) then
            sigma(sam1,sam2,sam3) = -sigma(sam1,sam2,sam3)
        else
            call random_number(judge2)
            delta_H = 2.d0*H_to !delta_H is minus
            p = exp(delta_H/Tr(itemp))
            !---------------------------------------------------------------
            !For higher energy, Boltzman probability distribution is adopted
            !---------------------------------------------------------------
            if (judge2 .lt. p) then
                sigma(sam1,sam2,sam3) = -sigma(sam1,sam2,sam3)
            end if
        end if
    end do
    !n2 times of the sampling
    do i1 = 1, n2
        !n3 steps interval of sampling
        do j1 = 1, n3
            !Periodic boundary condition
            do i = 2, npoin+1
                do j = 2, npoin+1
                    sigma(npoin+2,i,j) = sigma(2,i,j)
                    sigma(1,i,j) = sigma(npoin+1,i,j)
                    sigma(i,npoin+2,j) = sigma(i,2,j)
                    sigma(i,1,j) = sigma(i,npoin+1,j)
                    sigma(i,j,npoin+2) = sigma(i,j,2)
                    sigma(i,j,1) = sigma(i,j,npoin+1)
                end do
            end do
            !
            call Randi(2,npoin+1,sam1)
            call Randi(2,npoin+1,sam2)
            call Randi(2,npoin+1,sam3)
            !The interaction Hamiltionian
            H_in = -J_int*sigma(sam1,sam2,sam3)*(sigma(sam1+1,sam2,sam3) + sigma(sam1-1,sam2,sam3) + &
                                               & sigma(sam1,sam2+1,sam3) + sigma(sam1,sam2-1,sam3) + &
                                               & sigma(sam1,sam2,sam3+1) + sigma(sam1,sam2,sam3-1))
            !The Hamiltonian in the magnetic field
            H_se = -B*sigma(sam1,sam2,sam3)
            !The total Hamiltonian
            H_to = H_in + H_se 
            !---------------------------------------------------------------------------------
            !If H_to is positive, the spin orientation needs to be flipped to lower the energy
            !---------------------------------------------------------------------------------
            if (H_to .gt. 0.d0) then
                sigma(sam1,sam2,sam3) = -sigma(sam1,sam2,sam3)
            else
                call random_number(judge2)
                delta_H = 2.d0*H_to !delta_H is minus
                p = exp(delta_H/Tr(itemp))
                !---------------------------------------------------------------
                !For higher energy, Boltzman probability distribution is adopted
                !---------------------------------------------------------------
                if (judge2 .lt. p) then
                    sigma(sam1,sam2,sam3) = -sigma(sam1,sam2,sam3)
                end if
            end if
        end do
        !Obtain the total magnetic moment
        M = 0
        do i = 2, npoin + 1
            do j = 2, npoin + 1
                do k = 2, npoin + 1
                    M = M + sigma(i,j,k)
                end do
            end do
        end do
        !Save M data
        M_dat(M+mediate+1) = M_dat(M+mediate+1) + 1
    end do
    !Compute the even value of magnetic moment
    M_even = 0.d0
    do i = 1, Smax+1
        M_even = M_even + M_dat(i)*(i - mediate - 1)
    end do
    M_even = abs(M_even/n2)
    !output the data
    open(1,access='append',file='output_dat3D')
        !write(1,*) Tr(itemp),int(M_even)
        write(*,*) Tr(itemp),int(M_even)
    close(1)
    call cpu_time(t2)
    write(*,*) 'Time = ', t2-t1
end do
end program
!
!===============================
!Time seeds of the random number
!===============================
subroutine init_random_seed()
    implicit none
    integer :: i,n,clock
    integer,dimension(:),allocatable :: seed
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count = clock)
    seed = clock + 37*(/(i-1,i=1,n)/)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine
!
!================================
!Generate a random integer number
!================================
subroutine Randi(min,max,s)
    implicit none
    real*8 x
    integer min,max,s
    call random_number(x)
    s = floor(x*(max-min+1)) + min
end subroutine
!