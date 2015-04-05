    module datas
    real(kind=8),parameter :: PI = 4.0*atan(1.0)
    real(kind=8),parameter :: SMV = 1.0E-20
    real(kind=8),allocatable,dimension(:) :: w !solution variable
    real(kind=8),allocatable,dimension(:) :: flux !flux
    real(kind=8) :: dx !spacing in x-direction
    real(kind=8) :: dt !time step
    real(kind=8) :: cfl !cfl number
    integer :: num !number of grid points
    integer :: iter !iterations
    end module datas

    program main
    use datas
    integer :: i,nmax

    !control
    cfl = 0.9

    write(*,*) "Input max number of iteration"
    read(*,*) nmax

    !geometry
    num =400
    dx = 4.0/num

    !allocate array
    allocate(w(num))
    allocate(flux(num+1))

    !initial condition
    do i=1,num
        xpos = (-1.0+dx/2.0)+(4.0-dx)*(i-1.0)/(num-1.0+SMV)
        if (xpos >0.0 .and. xpos <1.0) then
            w(i) = sin(2*PI*xpos)
        else
            w(i) = 0.0
        end if
    end do

    iter = 1
    do while(.true.)
        call timestep()
        call calc_flux()
        call update()
        write(*,*) "iter,dt", iter, dt
        if (iter >= nmax) then
            call writeout()
            stop
        end if
        iter = iter+1
    end do
    end program main

    subroutine timestep()
    use datas
    integer :: i
    real(kind=8) :: umax

    umax = 0.0
    do i=1,num
        umax=max(umax,w(i))
    end do
    dt = cfl*dx/umax
    end subroutine timestep

    subroutine calc_flux()
    use datas
    integer :: i

    !boundary
    flux(1) = 0.5*w(i)**2
    flux(num+1) = 0.5*w(num)**2

    !inner
    do i=2,num
        if (w(i-1)>=w(i)) then !form a shock
            if (0<0.5*(w(i-1)+w(i))) then
                flux(i) = 0.5*w(i-1)**2
            else
                flux(i) = 0.5*w(i)**2
            end if
        else !form a rarefaction wave
            if (0<w(i-1)) then
                flux(i) = 0.5*w(i-1)**2
            else if (0>w(i)) then
                flux(i) = 0.5*w(i)**2
            else
                flux(i) = 0.0
            end if
        end if
    end do
    end subroutine calc_flux

    subroutine update()
    use datas
    integer :: i

    do i=1,num
        w(i) = w(i)+(flux(i)-flux(i+1))*dt/dx
    end do
    end subroutine update

    subroutine writeout()
    use datas
    integer :: i

    open(unit=10,file="out.dat")
    do i=1,num
        xpos = (-1.0+dx/2.0)+(4.0-dx)*(i-1.0)/(num-1.0+SMV)
        write(10,*) xpos,w(i)
    end do
    end subroutine writeout