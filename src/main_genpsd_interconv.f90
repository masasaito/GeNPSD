program genpsd_interconv


  use globals
  use hparx, only  : getCmdArgs  
  use genpsd, only : GeNPSD__interconv_ArbitraryPSD, GeNPSD__interconv_GammaPSD, &
                     GeNPSD__interconv_LogNormalPSD, GeNPSD__alpha, GeNPSD__geometry, GeNPSD__sizbin
  implicit none


  integer, parameter :: NSIZMAX = 10000
  integer            :: mgeoX, mgeoY, mlogP, mlogQ, msizP, msizQ, mpsd, mnor
  integer            :: psdtyp, nsiz, sizlog
  real(R_)           :: sizmin, sizmax
  integer            :: isiz, narg, ios
  real(R_) :: alpha          
  real(R_) :: rmed, lns, mu,  gam,  xtot
  character(len=256) :: argv(19) ! argument to be red
  character(len=1536):: argmsg   ! argment message
  real(R_), allocatable :: siz(:), psdX(:), psdY(:)
  real(R_)           :: g(3), p(3)
  real(R_) :: wrk(NSIZMAX,2)
  real(R_) :: psi, rge

  ! arguments
  argmsg = 'Usage1: genpsd_interconv psdtyp mgeoX msizP mlogP mgeoY msizQ mlogQ mnor g1 g2 g3\n'&
         //'                                          nsiz sizmin sizmax sizlog mpsd p1 p2 p3 > output\n'&
         //'Usage2: genpsd_interconv psdtyp mgeoX msizP mlogP mgeoY msizQ mlogQ mnor g1 g2 g3 < input > output\n\n'&
         //'  Converting PSDs from dX/dp (dX/dlnp) to dY/dq (dY/dlnq)\n'                   &
         //'  Usage1 for psdtyp>0, Usage2 for psdtyp=-1\n'                                 &
         //'  input    : A use-specified size bin (1st column) and PSD (2nd column)\n'     &
         //'  output   : A use-specified PSD\n'                                            &
         //'  mgeoX    : Geometric variable X, 0: number N, 1: proj-area A, 2: volume V\n' &
         //'  mgeoY    : Geometric variable Y\n'                                           &
         //'  msizP    : Size descriptor p,    0: r_g, 1: r_a, 2: r_v, 3: r_e\n'           &
         //'  msizQ    : Size descriptor q,    4: D_g, 1: D_a, 6: D-v, 7: D_e\n'           &
         //'  mlogP    : 0:dX/dp, 1:dX/dlnp\n'                                             &
         //'  mlogQ    :\n'                                                                &
         //'  mnor     : PSD normalization 0:off, 1:on\n'                                  &
         //'  g1,g2,g3 : particle geometry, (g1,g2,g3) = (D_max,A,V)\n'                    &                   
         //'  nsiz     : # of size bins\n'                                                 &
         //'  sizmin   : Minimum size\n'                                                   &
         //'  sizmax   : Maximum size\n'                                                   &
         //'  sizlog   : Logarithmic interval, 0:off, 1:on\n'                              &
         //'  mpsd     : flag for PSD parameters, mpsd=0: (p1,p2,p3) = (X_tot,r_eff,v_eff)\n' &
         //'  p1,p2,p3 :                          mpsd=1: (p1,p2,p3) = (X_tot,r_med,S) for psdtyp=0\n'&
         //'           :                                  (p1,p2,p3) = (X_tot,gam,mu)  for psdtyp=1'


  ! read arguments
  narg = 11
  call getCmdArgs(narg, argv, argmsg)
  read (argv(1), *) psdtyp
  read (argv(2), *) mgeoX
  read (argv(3), *) msizP
  read (argv(4), *) mlogP
  read (argv(5), *) mgeoY
  read (argv(6), *) msizQ
  read (argv(7), *) mlogQ
  read (argv(8), *) mnor
  read (argv(9), *) g(1)
  read (argv(10),*) g(2)
  read (argv(11),*) g(3)

  ! additional arguments
  if (psdtyp == -1) then

     do isiz = 1, NSIZMAX
        read(*,*,iostat=ios) wrk(isiz,:) !siz(isiz), psdX(isiz)
        if (ios /= 0) exit
     end do
     if (ios /= 0) nsiz = isiz-1
     allocate (siz(nsiz), psdX(nsiz), psdY(nsiz) )
     siz(:)  = wrk(1:nsiz,1)
     psdX(:) = wrk(1:nsiz,2)
  else
     narg = 19
     call getCmdArgs(narg, argv, argmsg)
     read (argv(12),*) nsiz
     read (argv(13),*) sizmin
     read (argv(14),*) sizmax
     read (argv(15),*) sizlog
     read (argv(16),*) mpsd
     read (argv(17),*) p(1)
     read (argv(18),*) p(2)
     read (argv(19),*) p(3)
     allocate (siz(nsiz), psdX(nsiz), psdY(nsiz) )
     siz(:) = GeNPSD__sizbin(nsiz, sizmin, sizmax, sizlog) 
  end if


  ! particle geometry
  call GeNPSD__geometry(g(1), g(2), g(3), psi, rge)
  !psi: effective sphericty
  !rge: Vr/Ar =  r_e/r_g ratio


  if (psdtyp == -1) then ! arbitrary PSD
     psdY(:) = GeNPSD__interconv_ArbitraryPSD(mgeoX, mgeoY, msizP, msizQ, mlogP, mlogQ, mnor, &
                                              psi, rge, siz, psdX)

  else if (psdtyp == 0) then ! log-normal distribution

     ! Original PSD
     if (mpsd == 0) then
        alpha = GeNPSD__alpha(3, msizP, psi, rge)
        lns  = sqrt(log(p(3)+1.0_R_))         ! ln(S) where S is the geometric stddev  
        rmed = alpha*p(2)*exp(-2.5_R_*lns**2) ! median radius

        ! physical variable
        if      (mgeoX == 1) then ! A
           rmed = rmed * exp(2.0_R_*lns**2)
        else if (mgeoX == 2) then ! V
           rmed = rmed * exp(3.0_R_*lns**2)
        end if
     else
        rmed = p(2)
        lns  = log(p(3))
     end if
     xtot = p(1)


     ! Nonspherical PSD interconversion
     psdY(:) =  GeNPSD__interconv_LogNormalPSD(mgeoX, mgeoY, msizP, msizQ, mlogQ, mnor, &
                                               psi, rge, siz, xtot, rmed, lns)

  else if (psdtyp == 1) then ! gamma distribution

     ! Original PSD
     if (mpsd == 0) then
        alpha = GeNPSD__alpha(3, msizP, psi, rge)
        gam = alpha/(p(2)*p(3))    ! slope parameter
        mu  = 1.0_R_/p(3) - 2.0_R_ ! dispersion parameter 

        ! physical variable
        if      (mgeoX == 1) then ! A
           mu = mu + 2.0_R_
        else if (mgeoX == 2) then ! V
           mu = mu + 3.0_R_
        end if
     else
        gam = p(2)
        mu  = p(3)
     end if
     xtot = p(1)


     ! Nonspherical PSD interconversion
     psdY(:) = GeNPSD__interconv_GammaPSD(mgeoX, mgeoY, msizP, msizQ, mlogQ, mnor, &
                                          psi, rge, siz, xtot, gam, mu)

  else
     stop "psdtyp must be (-1, 0, 1)"
  end if


  ! output
  do isiz = 1, nsiz
     print*, siz(isiz),  psdY(isiz)
  end do
  deallocate ( siz, psdX, psdY )


end program genpsd_interconv

