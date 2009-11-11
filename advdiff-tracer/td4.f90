program advection

      ! computes 1D diffusion using FVM
      ! arguments:
      ! - advection scheme: [b, s] for backward or staggered
      ! - boundaries condition: [p, n] for periodic or no-flux
      ! say: ./a.out s n will compute with the staggered scheme and no-flux bc

      implicit none

      real, parameter    :: pi = acos(-1.0)

      ! command line parameters handling
      character(len = 1) :: iname, scheme, bc

      ! helpers
      integer            :: i, j, k             ! general purpose iterator

      ! physical constants
      integer, parameter :: L  = 100000        ! domain length (say, 100 km)
      integer, parameter :: H  = 3000          ! domain height (say, 3 km)
                                                ! not really a parameter, actually
      
      ! integration domain
      integer, parameter :: Nx = 100            ! number of volume controls over L
      integer, parameter :: Nz = 30             ! number of volume controls over Z
      real, parameter    :: dVx_step = real(L) / real(Nx)   ! hence, the control step
      real, parameter    :: dVz_step = real(H) / real(Nz)   ! hence, the control step

      ! winds and diffusion
      integer, parameter  :: u_min = 2          ! wind range
      integer, parameter  :: u_max = 10
      integer, parameter  :: w_max = 0          ! set to 0 to disable w wind
      real, dimension(Nz) :: u
      real, dimension(Nx) :: w

      ! tracer
      real, dimension(Nx, Nz) :: tracer         ! the conservative, passive tracer
      real, dimension(Nx, Nz) :: tracer_old     ! prior state backup
      real, dimension(Nx)     :: xmid           ! volume control midpoints
      real, dimension(Nz)     :: zmid           ! volume control midpoints

      ! gaussian distribution parameters
      !real, parameter :: shift  = Nx / 5.0
      !real, parameter :: DeltaX = Nx / 15.0
      !real, parameter :: DeltaT = 0.6 * L / u

      ! dt = CFL_u * dx / u, with 0 <  CFL_u <= 1
      ! this model features an horizontal wind only, so there is no dependence
      ! over dVz_step

      ! Courant–Friedrichs–Lewy condition
      ! it depends on z so it reproduces u(z) through alpha, though alpha is not
      ! directly used in the computations
      real, dimension(Nz) :: CFL_u                ! to be computed after zmid
      real, dimension(Nx) :: CFL_w                ! to be computed after w
                                                  ! within the time loop
      !real, parameter :: dt    = CFL_u * dVx_step / abs(u)
      real, parameter :: dt = 0.5               ! bold assumption
      !real, parameter :: alpha = dt / dVx_step

      ! turbulence
      real                :: turb_diff = 0
      real, parameter     :: Kz_basis  = 25        ! m^2/s max Kz
      real, dimension(Nz) :: Kz
      real, dimension(Nz) :: beta

      ! time iterations
      integer :: it = 1
      integer :: nt = 1001

      ! fluxes within control cells
      real, dimension(Nx, Nz) :: flux
      real, dimension(Nx, Nz) :: residual

      ! file handling
      character(len = 3) :: setchar
      integer            :: iunit       = 10        ! I/O basefile
      character(len = 6) :: export_file = "d_...."  ! if nt changes, edit this line

      ! figure handling
      integer :: mod_fig = 10     ! output each mod_fig steps

      ! debug output
      print *, "L = ",      L
      print *, "H = ",      H
      !print *, "shift = ",  shift
      !print *, "DeltaX = ", DeltaX
      !print *, "DeltaT = ", DeltaT
      print *, "Nx = ",     Nx
      print *, "Nz = ",     Nz
      print *, "dVx_step",  dVx_step
      print *, "dVz_step",  dVz_step
      !print *, "dt = ",     dt
      !print *, "alpha = ",  alpha
      !print *, "beta = ",   beta
      print *, ""
      print *, nt, "time steps"
      print *, ""
      !read(*,*)  ! uncomment to pause the program

      ! -----------------------------------------------------------------------
      ! 0. command line parameters: scheme, boundaries conditions

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

      print *, "advection scheme: ",    scheme
      print *, "boundary conditions: ", bc
      print *, ""

      ! -----------------------------------------------------------------------
      ! 1. initial tracer distribution, misceallenous initialization

      ! first, compute the mean values of the piecewise constant 2D solution
      xmid(:) = (/ (0.5 * (2 * i + 1), i = 1, Nx) /)
      zmid(:) = (/ (0.5 * (2 * i + 1), i = 1, Nz) /)

      u       = u_min + (u_max - u_min) * zmid(:) * dVx_step / H
      CFL_u   = u(:) * dt / dVx_step
      !print *, "CFL_u: ", CFL_u
      !read(*,*)

      ! initial state
      ! uniform over z
      forall(i = 1:Nz)
        !tracer(:, i) = 1 * exp(-((xmid(:) - shift) / DeltaX)**2)

        !tracer(:, i) = sin(4*pi*xmid(:)/L) + &
                       !sin(17*pi*xmid(:)/L) + &
                       !(-1 + mod(xmid(:),2.0))

        tracer(:,i) = 0
        ! yet force a (constant, see time loop) source at some location
        tracer(Nx/6, 3 * Nz / 4) = 1
        ! TODO passer en constante les paramètres de localisation de la source
      endforall

      ! space&time-integrated turbulent mixing coeff.
      !print *, Kz
      !read(*,*)
      Kz(:)   = Kz_basis * (1 - tanh(zmid(:) * dVz_step / H))
      beta(:) = Kz(:) * dt / (dVx_step * dVz_step)
      !print *, beta
      !read(*,*)

      print *, 'exporting initial tracer distribution...'
      open(file = 'distributions.dat', unit = iunit, status = 'replace')
      do i = 1, Nx
        do j = 1, Nz
          write(iunit, fmt = "(f0.4, 1x)", advance = "no") tracer(i,j)
        enddo
        write(iunit,*)
      enddo
      print *, 'export tracer done.'
      print *, ""
      !read(*,*)

      ! -----------------------------------------------------------------------
      ! 2. 1D advection, using finite volumes

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
      
      ! -----------------------------------------------------------------------
      ! 3. MAIN LOOP
      advect: do it = 1, nt
        ! compute the flux budget over the ith cell (aka. the residual),
        ! then update the passive tracer within the cell

        ! don't mess up different distribution states
        tracer_old = tracer

        if (scheme .eq. "b") then
          ! backward advection scheme

          do i = 1, Nx
            ! compute w(x,t) and CFL_w(w(x,t))
            w     = w_max * (0.2 * sin(30*pi*xmid(:)*dVx_step/L-(it*dt)/50.0) &
                           + 1.0 * sin(10*pi*xmid(:)*dVx_step/L-(it*dt)/66.0))
            CFL_w = w(:) * dt / dVz_step 

            !write(setchar(1:3), "(i3)") it
            !! TODO regexp qui remplace les espaces par des 0
            !open (iunit, file = "./data/wind/w/" // setchar, form = 'formatted', status = 'replace')
            !write(iunit, fmt = "(f0.4)") w
          enddo

          ! processing the inner domain
          do i = 2, Nx
            do k = 2, Nz-1
              ! compute turbulent diffusivity if necessary
              if (Kz(k) .ne. 0.0) then
                !print *, k
                turb_diff = beta(k+1) * (tracer_old(i, k+1) - tracer_old(i, k)) &
                          - beta(k)   * (tracer_old(i, k) - tracer_old(i, k-1))
              else
                turb_diff = 0
              endif

              if (CFL_w(i) > 0) then
                !print *, "vent positif sur la colonne"
                tracer(i, k) = tracer_old(i, k) &
                             + CFL_u(k)*(tracer_old(i-1, k) - tracer_old(i, k)) &
                             + CFL_w(i)*(tracer_old(i, k-1) - tracer_old(i, k)) &
                             + turb_diff
              else
                !print *, "vent négatif sur la colonne"
                tracer(i, k) = tracer_old(i, k) &
                             + CFL_u(k)*(tracer_old(i-1, k) - tracer_old(i, k)) &
                             + CFL_w(i)*(-tracer_old(i, k+1) + tracer_old(i, k)) &
                             + turb_diff
              endif
            enddo

            ! force a constant source at some location within the domain
            tracer(Nx/6, 3 * Nz / 4) = 1
          enddo

          ! boundaries condition
          !if (bc .eq. "p") then
            !! periodic boundaries condition
            !! TODO
          !else
            !! no flux boundaries condition
            !! ie. q_1,n = q_2,n and q_im,n = q_im-1,n

            !forall(i = 1:Nx)
              !tracer_old(i, 1) = tracer_old(i, 2)
              !tracer_old(i, Nz) = tracer_old(i, Nz-1)
            !endforall

            !forall(j = 1:Nz)
              !tracer_old(1, j) = tracer_old(2, j)
              !tracer_old(Nx, j) = tracer_old(Nx-1, j)
            !endforall
          !endif ! bc

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
          print *, "not yet staggered scheme, aborting..."
          return
        else
          print *, "unknown advection scheme, aborting..."
          return
        endif

        ! export current iteration distribution state
        if (mod(it, mod_fig) == 1) then
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

      enddo advect
      print *, 'advecting and diffusing tracer done.'

end program

