program advectodiffusion

      ! -----------------------------------------------------------------------
      ! Computes a 2D advection-diffusion using finite volumes modeling (FVM)
      ! for a passive, atmospheric tracer. Different cases are handled (w/o
      ! vertical wind, w/o vertical turbulent diffusion, w/o chemical lifetime).
      !
      ! Command line call arguments are:
      ! - the advection scheme: [b, s] for backward or staggered
      ! - the boundaries condition: [p, n] for periodic or no-flux
      !
      ! say: "./a.out s n" would switch the staggered scheme and the no-flux
      ! boundaries conditions on
      ! -----------------------------------------------------------------------

      implicit none

      ! -----------------------------------------------------------------------

      real, parameter    :: pi = acos(-1.0)

      ! command line parameters handling
      character(len = 1) :: iname               ! holds a parameter...
      character(len = 1) :: scheme              ! char flag for the advection scheme
                                                ! - b: backward
                                                ! - s: staggered (centered)
      character(len = 1) :: bc                  ! char flag for the boundaries condition
                                                ! - n: no fluxes
                                                ! - p: periodic

      ! helpers, iterators
      integer            :: i, j, k             ! general purpose iterators

      ! physical domain and constants
      integer, parameter :: L  = 100000           ! domain length (say, 100 km)
      integer, parameter :: H  = 3000           ! domain height (say, 3 km)
                                                ! (not really a parameter, actually)
      
      ! integration domain
      integer, parameter :: Nx = 500            ! number of volume controls over L,
      integer, parameter :: Nz = 150            ! number of volume controls over Z, hence...
      real, parameter    :: dVx_step = real(L) / real(Nx)   ! the horizontal control step
      real, parameter    :: dVz_step = real(H) / real(Nz)   ! the vertical control step

      ! winds
      real, parameter     :: u_min = 5          ! horizontal wind range (min)
      real, parameter     :: u_max = 30         !                       (max)
      real, parameter     :: w_max = 5          ! vertical wind (do *not* set to 0 to disable,
                                                ! set w_enabled to .FALSE. instead)
      real, dimension(Nz) :: u                  ! horizontal wind
      real, dimension(Nx) :: w                  ! vertical wind

      ! gaussian distribution parameters
      real, parameter :: shift  = Nx / 5.0
      real, parameter :: DeltaX = Nx / 15.0
      !real, parameter :: DeltaT = 0.6 * L / u

      ! tracer
      real, dimension(Nx, Nz) :: tracer         ! the conservative, passive tracer
      real, dimension(Nx, Nz) :: tracer_old     ! prior state backup
      real, dimension(Nx, Nz) :: means          ! means over some time span
      real, dimension(Nx)     :: xmid           ! volume control midpoints on horizontal axis
      real, dimension(Nz)     :: zmid           ! volume control midpoints on vertical axis

      ! Courant–Friedrichs–Lewy conditions
      real, parameter     :: CFL_basis = 0.1    ! used to assess dt, should be, like, low
      real, dimension(Nz) :: CFL_u              ! to be computed after zmid
      real, dimension(Nx) :: CFL_w              ! to be computed after w within the time loop

      ! time stepping
      !real, parameter :: dt = CFL_u * dVx_step / abs(u)
      !real, parameter :: dt = 0.5               ! bold assumption
      real :: dt = 0                            ! will be dynamically computed later on

      ! vertical turbulence parametrization (diffusion process)
      real                :: turb_diff = 0      ! dt-step value of the turbulence
      real, parameter     :: Kz_basis  = 25     ! max Kz (m^2/s)
      integer, parameter  :: z0 = 1000           ! height scale for Kz(z, z0)
      real, dimension(Nz) :: Kz                 ! turbulent coeff.
      real, dimension(Nz) :: beta               ! time-space integrated coeff.

      ! chemistry parametrization
      real, parameter :: lifetime = 3600        ! 1/lifetime is the decrease rate, in s

      ! discrete time iterations
      integer :: it                             ! time step iterator
      integer :: nt = 2001                      ! number of steps

      ! avereged field computations
      integer :: mean_sup_limit = 2000             ! *must* be greater than 1000
                                                ! (set to 0 to *force* disabling over the flag)

      ! file handling
      character(len = 3) :: setchar             ! holds some characters
      integer            :: iunit       = 10    ! I/O basefile
      character(len = 6) :: export_file = "d_...."  ! if nt changes, edit this line

      ! figure handling
      integer :: mod_backup = 10                ! save data every mod_backup steps

      ! flags
      logical :: bc_enabled           = .TRUE. ! wether to use some boundaries conditions to begin with

      logical :: initial_gaussian     = .FALSE. ! gaussian initial tracer distribution
      logical :: initial_sinusoidal   = .FALSE. ! sinusoidal initial tracer distribution
      logical :: initial_void         = .TRUE.  ! no tracer (takes precedence over other initial states)
      
      logical :: w_enabled            = .FALSE. ! vertical wind flag
      
      logical :: force_cst_source     = .FALSE.  ! force a constant source at some location
      
      logical :: v_diff_enabled       = .TRUE. ! turbulence flag
      logical :: constant_v_diff      = .FALSE. ! whether Kz is a constant
      logical :: force_v_diff_at_ground = .TRUE. ! force a turbulent flux at some location
      
      logical :: lifetime_enabled     = .TRUE. ! lifetime flag
      
      logical :: compute_means        = .TRUE. ! whether to compute means on tracer
                                                ! (!! will be set to true if nt > mean_sup_limit;
                                                !  mean will be computed between 1000 and mean_sup_limit steps)

      ! -----------------------------------------------------------------------
      ! debug output
      print *, "L = ",      L
      print *, "H = ",      H
      print *, "shift = ",  shift
      print *, "DeltaX = ", DeltaX
      !print *, "DeltaT = ", DeltaT
      print *, "Nx = ",     Nx
      print *, "Nz = ",     Nz
      print *, "dVx_step",  dVx_step
      print *, "dVz_step",  dVz_step
      print *, "dt = ",     dt
      !print *, "alpha = ",  alpha
      !print *, "beta = ",   beta
      print *, ""
      print *, nt, "time steps"
      print *, ""
      !read(*,*)  ! uncomment to pause the program

      ! -----------------------------------------------------------------------
      ! 0. command line parameters handling

      ! default parameters
      scheme = "b"     ! backward advection scheme
      bc     = "p"     ! periodic boundaries condition

      if (iargc() .gt. 2) then
        write(*, '(/"too many arguments (2 max: [scheme (b, s) boundaries (p, n)])."/)')
        return
      else if(iargc() .ne. 0) then
        
        ! advection scheme
        call getarg(1, iname)
        scheme = iname

        if (iargc() .eq. 2) then
          ! boundaries condition
          call getarg(2, iname)
          bc = iname
        endif
      endif

      ! debug...
      print *, "* advection scheme: ",    scheme
      print *, "* boundary conditions: ", bc
      print *, ""

      ! -----------------------------------------------------------------------
      ! 1. time-stepping, initial tracer distribution, misc. initializations

      ! time-stepping
      if (w_enabled) then
        print *, "w enabled, dt"
        dt = min(CFL_basis * dVx_step / u_max, CFL_basis * dVz_step / w_max)
      else
        print *, "w disabled, dt"
        dt =  CFL_basis * dVx_step / u_max
      endif
      !real, parameter :: alpha = dt / dVx_step
      print *, dt
      !read(*,*)

      ! first, init the computation points of the (piecewise constant) 2D solution
      xmid(:) = (/ (0.5 * (2 * i + 1), i = 1, Nx) /)
      zmid(:) = (/ (0.5 * (2 * i + 1), i = 1, Nz) /)

      u       = u_min + (u_max - u_min) * zmid(:) * dVx_step / H
      CFL_u   = u(:) * dt / dVx_step
      !print *, "CFL_u: ", CFL_u
      !read(*,*)

      ! initial state
      ! TODO passer en constante les paramètres de localisation de la source
      if (initial_void) then
        tracer = 0
        if (force_cst_source) then
          ! yet force a (constant, see time loop) source at some location
          tracer(Nx/6, 3 * Nz / 4) = 1
        endif
      else
        if (initial_gaussian) then
          forall(i = 1:Nz)
            tracer(:, i) = 1 * exp(-((xmid(:) - shift) / DeltaX)**2)
          endforall
        else if (initial_sinusoidal) then
         forall(i = 1:Nz)
           tracer(:, i) = sin(4*pi*xmid(:)/L) &
                        + sin(17*pi*xmid(:)/L) &
                        + (-1 + mod(xmid(:),2.0))
          endforall
        else
          print *, "[E] unknown initial tracer distribution, aborting"
          stop
        endif
      endif

      ! space&time-integrated turbulent mixing coeff.
      if (v_diff_enabled) then
        !print *, Kz
        !read(*,*)
        if (constant_v_diff) then
          Kz(:) = Kz_basis
        else
          Kz(:) = Kz_basis * (1 - tanh(zmid(:) * dVz_step / z0))
        endif
        beta(:) = Kz(:) * dt / (dVz_step * dVz_step)
        !print *, beta
        !read(*,*)
      else
        turb_diff = 0 ! not really necessary, but...
      endif

      print *, 'exporting initial tracer distribution...'
      open(file = 'distributions.dat', unit = iunit, status = 'replace')
      do i = 1, Nx
        do j = 1, Nz
          write(iunit, fmt = "(f0.4, 1x)", advance = "no") tracer(i,j)
        enddo
        write(iunit,*)
      enddo
      print *, 'export done.'
      print *, ""
      !read(*,*)

      ! -----------------------------------------------------------------------
      ! 2. determining the advection scheme and the boundaries condition

      print *, 'advecting and diffusing tracer...'

      ! staggered scheme will output into data/s/
      ! while backward scheme will output into data/b/
      if (scheme .eq. "s") then
        print *, "  advection scheme: staggered"
      else
        print *, "  advection scheme: not staggered, falling back on backward"
        scheme = "b"
      endif

      ! boundaries condition
      if (bc .eq. "n") then
        print *, "  boundaries condition: no flux"
      else
        print *, "  boundaries condition: falling back on periodic"
        bc = "p"
      endif

      ! enable means computations?
      if (mean_sup_limit < 1001) then
        print *, "[W] mean_sup_limit must be greater than 1000. Disabling mean computation."
        compute_means = .FALSE. ! ensure that
      else
        if (nt > mean_sup_limit) then
          if (compute_means) then
          !compute_means = .TRUE.
            print *, "  means computation enabled"
          endif
        else
          compute_means = .FALSE.
          print *, "[W] not enough time steps to compute accurate mean. Mean computation disabled."
        endif
      endif
 
      ! -----------------------------------------------------------------------
      ! 3. MAIN LOOP, time-space integration using FVM
      advect: do it = 1, nt
        ! boundaries condition
        if (bc_enabled) then
          if (bc .eq. "p") then
            ! periodic boundaries condition
            print *, "periodic bc are not more available"
            ! TODO
          else if (bc .eq. "n") then
            ! no-flux boundaries condition
            ! ie. q_1,n = q_2,n and q_im,n = q_im-1,n at anytime

            forall(i = 1:Nx)
              !tracer(i, 1) = tracer(i, 2)
              !tracer(i, Nz) = tracer(i, Nz-1)
              tracer(i, 1) = 0
              tracer(i, Nz) = 0
            endforall

            forall(j = 1:Nz)
              tracer(1, j) = tracer(2, j)
              ! with a constant westward wind, let's forget about the other bc
              !tracer_old(Nx, j) = tracer_old(Nx-1, j)
            endforall
          else
            print *, "unknown boundaries conditions, assuming none"
          endif ! bc
        endif

        ! don't mess up different distribution states
        tracer_old = tracer

        if (compute_means .and. ((it > 999) .and. (it <= mean_sup_limit))) then
          ! gather data in order to compute a field mean over some time range
          !print *, "computing means for ", it, "step"
          means = means + tracer_old
        endif

        if (scheme .eq. "b") then
          ! backward advection scheme

          if (w_enabled) then
            ! compute w(x,t) and the related CFL_w(w(x,t))
            do i = 1, Nx
              w     = w_max * (0.2 * sin(30*pi*xmid(:)*dVx_step/L-(it*dt)/50.0) &
                            +  1.0 * sin(10*pi*xmid(:)*dVx_step/L-(it*dt)/66.0))
              CFL_w = w(:) * dt / dVz_step 

              !write(setchar(1:3), "(i3)") it
              !! TODO regexp qui remplace les espaces par des 0
              !open (iunit, file = "./data/wind/w/" // setchar, form = 'formatted', status = 'replace')
              !write(iunit, fmt = "(f0.4)") w
            enddo
          endif

          ! processing the inner domain
          do i = 2, Nx
            do k = 2, Nz-1

              if (v_diff_enabled .and. Kz(k) /= 0.0) then
                ! compute turbulent diffusivity if necessary
                if (force_v_diff_at_ground .and. (k == 2 .and. (i > Nx / 6 .and. i < Nx / 2))) then
                  ! set a unitary, constant turbulent flux at some part of the ground
                  !print *, it, k, i
                  !read(*,*)
                  turb_diff = 1
                else
                  ! simple parametrization of the turbulence as a diffusion process
                  turb_diff = beta(k+1) * (tracer_old(i, k+1) - tracer_old(i, k)) &
                            - beta(k)   * (tracer_old(i, k) - tracer_old(i, k-1))
                endif
              !else
                !turb_diff = 0
              endif

              tracer(i, k) = tracer_old(i, k) &
                           + CFL_u(k)*(tracer_old(i-1, k) - tracer_old(i, k))

              if (w_enabled) then
                ! the advection scheme depends on the sign of the local vertical
                ! wind ('could have done the same for the horizontal wind, actually)
                if (CFL_w(i) > 0) then ! CFL sign
                  !print *, "vertical winds are upward"
                  tracer(i, k) = tracer(i, k) &
                               + CFL_w(i)*(tracer_old(i, k-1) - tracer_old(i, k))
                else
                  !print *, "vertical winds are downward"
                  tracer(i, k) = tracer(i, k) &
                               + CFL_w(i)*(-tracer_old(i, k+1) + tracer_old(i, k))
                endif
              endif

              if (v_diff_enabled) then
                tracer(i, k) = tracer(i, k) + turb_diff
              endif
            enddo ! loop over k

            if (force_cst_source) then
              ! force a constant source at some location within the domain
              tracer(Nx/6, 3 * Nz / 4) = 1
            endif
          enddo ! loop over i

        else if (scheme .eq. "s") then
          ! TODO FIXME pas encore refait sur le modèle du schéma amont
          ! staggered (periodic boundary conditions)
          ! work for CFL_u like 0.2...

          ! TODO no flu!x bc within the staggered scheme
          !forall(j = 1:Nz)
            !! periodic boundaries conditions
            !tracer(1, j)  = tracer_old(1, j) &
                      !+ CFL_u(j)/2*(tracer_old(Nx, j) - tracer_old(2, j)) &
                      !+ beta*(tracer_old(2, j) - 2*tracer_old(1, j) + tracer_old(Nx, j))

            !! processing the inner domain
            !forall(i = 2:Nx-1)
              !tracer(i, j) = tracer_old(i, j) &
                        !+ CFL_u(j)/2*(tracer_old(i-1, j) - tracer_old(i+1, j)) &
                        !+ beta*(tracer_old(i+1, j) - 2*tracer_old(i, j) + tracer_old(i-1, j))
            !endforall

            !tracer(Nx, j)  = tracer_old(Nx, j) &
                        !+ CFL_u(j)/2*(tracer_old(Nx-1, j) - tracer_old(1, j)) &
                        !+ beta*(tracer_old(1, j) - 2*tracer_old(Nx, j) + tracer_old(Nx-1, j))
          !endforall
          print *, "the staggered scheme is not available yet, aborting"
          return
        else
          print *, "unknown advection scheme, aborting"
          return
        endif

        if (lifetime_enabled) then
          ! chemical lifetime deprecency
          tracer = tracer * exp(-dt/lifetime)
        endif

        ! export current iteration distribution state
        if (mod(it, mod_backup) == 1) then
          write(export_file(3:6),'(i4.4)') it
          open (iunit, file = "./data/" // scheme // "/" // export_file, form = 'formatted', status = 'replace')
          do i = 1, Nx
            do j = 1, Nz
              write(iunit, fmt = "(f0.4, 1x)", advance = "no") tracer(i,j)
            enddo
            write(iunit, *)
          enddo
        endif
        !read(*,*)

      enddo advect ! time loop

      if (compute_means) then
        ! hope this won't break for means got too big...
        means = means / (mean_sup_limit - 1000)

        print *, 'exporting mean tracer distribution...'
        open(file = 'mean.dat', unit = iunit, status = 'replace')
        do i = 1, Nx
          do j = 1, Nz
            write(iunit, fmt = "(f0.4, 1x)", advance = "no") means(i,j)
          enddo
          write(iunit,*)
        enddo
        print *, 'export done.'
        print *, ""
      endif

      print *, 'advecting and diffusing tracer done.'

end program

