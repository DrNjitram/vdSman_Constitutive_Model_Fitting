! Startup from zero stress of flow with a constant velocity gradient
! 2D flow. combined shear and planar extensional flow. 
! The model is an elastoviscoplastic model (Saramito 2009).
! The integration scheme is second-order Runge Kutta (Heun)

program startup

  use viscoelastic_models_2D_m

  implicit none

  integer, parameter :: model      = 2, &  ! 2: Oldroyd model 3: Giesekus
                        nmodes     = 1, &  ! number of modes
                        flowtype   = 0, &  ! 2D !axisymmetric flow needs to be added
                        ncompc     = 3, &  ! number of components of c
                        ncompt     = 3, &  ! number of components of tau
                        ncompg     = 4     ! number of components of gradv
                  

  real(dp) ::            modulus,       & ! modulus 
                         lambda,        & ! relaxation time
                         tau_y,         & ! yield stress
                         nexp,          & ! n power-law exponent
                         alpha,         & ! alpha parameter
                         gammadot_cr,   & !eta_eff parameter: gammadot_cr
                         eta_inf,       & ! to avoid discontinuity
                         I1c,           &

                         gamma_0,       & ! strain amplitude
                         freq             !frequency


  integer :: step, numtimesteps_percycle, elemcycle, cycles
  real(dp) :: c(1,ncompc,nmodes), rhs(1,ncompc,nmodes), tau(1,ncompt), gammadot, gammadotlist(257), timestep, sumG1, sumG2, G1, G2,&
             Lxx, Lyy, Lxy,Lyx, Ltt, Dxx, Dxy, Dyy, gamma_eff, cxx, cxy, cyy, ctt
  real(dp) :: k1(1,ncompc,nmodes), gradv(1,ncompg), L(ncompg)
  character(len=40) :: filename

  type(vemodel_t) :: vemodel

  ! namelist for input of variables; read from standard input
  namelist /comppar/ modulus, lambda, gamma_0,  freq, cycles, filename, numtimesteps_percycle, gammadotlist

  read ( unit=*, nml=comppar )


! define the model

  call create_viscoelastic_model ( model, vemodel, nmodes, flowtype)
  

! set material parameters
  vemodel%modulus = modulus
  vemodel%lambda = lambda

! initialize c
  c(1,:,1) = [ 1, 0, 1 ]

! calculate timestep
  timestep = 1 / (numtimesteps_percycle*freq)

! open output file

  open ( unit=13, file=filename, recl=300 )

! stepping using RK2


  do elemcycle = 1, cycles !for each cycle start a loop

    ! initial condition for the first step of firs cycle
    if ( elemcycle == 1) then

      step = 0
      write(13,'(1i12.4, 11e12.4)') elemcycle, step * timestep, real(0.00_dp), real(0.00_dp), real(0.00_dp) ,& 
      vonmises_2D(tau(1,:)),gammadot, real(0.00_dp), real(0.00_dp), gamma_0
      
    
    ! initial condition for the first step of remaining cycles
    else
    
      step = 0               !timestep will be one
      call stress_viscoelastic_2D ( vemodel, c, tau )
      write(13,'(1i12.4, 11e12.4)') elemcycle, step * timestep, tau, vonmises_2D(tau(1,:)),gammadot,&
      0.00_dp

    end if 

    ! initialize G1 and G2 for each cycle
    sumG1 = 0
    sumG2 = 0 

    step = 1
    do step = 1, numtimesteps_percycle
        
  !   initialize velocity gradient
      gammadot = gammadotlist(step)
      L = [ 0._dp, 0._dp, gammadot, 0._dp ]
      gradv(1,:) =  L

      print *, "      "
      print *, "entering step", step, "in cycle", elemcycle

  !  additional stuff
  !   from valocity gradient calcuate gamma_eff and eta_eff

      Lxx = gradv(1,1)
      Lxy = gradv(1,2)
      Lyx = gradv(1,3)
      Lyy = gradv(1,4)
      !if ( flowtype == 1 ) then 
        !Ltt  = gradv(1,5)  ! axisymmetric, if axisymmetric is should have 5 elements
      !end if

      Dxx = (Lxx + Lxx)/2
      Dxy = (Lxy + Lyx)/2
      Dyy = (Lyy + Lyy)/2


  !   from initial (non updated) conformation tensor calcuate G
      cxx = c(1,1,1)
      cxy = c(1,2,1)
      cyy = c(1,3,1)
      !if ( flowtype == 1 ) then 
        !ctt = c(4) ! axisymmetric, if axisymmetric it should have 4 elements at least
      !end if


  !   right-hand side 

      call rhs_viscoelastic_2D ( vemodel, gradv, c, rhs ) 
      !this calcuates rhs the change in conf tensor per unit time (delc/delt)
      !rhs is basically del(c) or \dot{c), change of c wrto time


  !   do intermediate step

      k1 = timestep * rhs !by multiplying with timestep we get the change in conformation tensor

  !   right-hand side 

      call rhs_viscoelastic_2D ( vemodel, gradv, c+k1, rhs ) !then we calcuate another step

  !   do step

      c = c + ( k1 + timestep * rhs ) / 2 !then we average the steps

  !   compute stress tensor

      call stress_viscoelastic_2D ( vemodel, c, tau )

      sumG1 = sumG1 + tau(1,2) * sin(freq*step*timestep*2*3.14)*timestep
      sumG2 = sumG2 + tau(1,2) * cos(freq*step*timestep*2*3.14)*timestep

      write(13,'(1i12.4, 11es12.4)') elemcycle, step * timestep, tau, vonmises_2D(tau(1,:)),gammadot,0.00_dp

      print *, "done with step", step

    end do

    ! finalize G1 and G2

    G1 = sumG1* (freq*2*3.14) / (gamma_0*3.14)
    G2 = sumG2* (freq*2*3.14)  / (gamma_0*3.14)

    open ( unit=10, file='strain_sweep_lin.txt', status = 'old', position = 'append' )
    write(10,'(11es12.4)') gamma_0, G1, G2
    close ( unit=10 )

  end do

  close(13)

! delete the model

  call delete ( vemodel )

end program startup
