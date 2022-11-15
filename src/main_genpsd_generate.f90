program genpsd_generate


  use globals
  use hparx, only : getCmdArgs
  use genpsd, only : GeNPSD__alpha, GeNPSD__radius2diameter, GeNPSD__geometry, &
                     GeNPSD__sizbin
  implicit none


  integer            :: psdtyp, mgeo, mlog, msiz, nsiz, sizlog, mpsd
  integer            :: isiz, narg
  character(len=256) :: argv(15) ! argument to be red
  character(len=1024):: argmsg   ! argment message
  real(R_)           :: sizmin, sizmax, p(3), g(3)
  real(R_), allocatable :: siz(:), psd(:)
  real(R_) :: psi, rge, fac, alpha         
  real(R_) :: pmed, lns, mu, gam

  ! arguments
  argmsg = ' Usage: genpsd_generate psdtyp mgeo msiz mlog nsiz sizmin sizmax sizlog mpsd p1 p2 p3 g1 g2 g3 > output\n\n'     &
         //'  output   : A use-specified PSD\n'                                            &
         //'  psdtyp   : PSD shape,            0: log-normal; 1: gamma\n'                  &
         //'  mgeo     : Geometric variable X, 0: number N, 1: proj-area A, 2: volume V\n' &
         //'  msiz     : Size descriptor p,    0: r_g, 1: r_a, 2: r_v, 3: r_e\n'           &
         //'           :                       4: D_g, 1: D_a, 6: D-v, 7: D_e\n'           &
         //'  mlog     : 0:dX/dp, 1:dX/dlnp\n'                                             &
         //'  nsiz     : # of size bins\n'                                                 &
         //'  sizmin   : Minimum size\n'                                                   &
         //'  sizmax   : Maximum size\n'                                                   &
         //'  sizlog   : Logarithmic interval, 0:off, 1:on\n'                              &
         //'  mpsd     : flag for PSD parameters, mpsd=0: (p1,p2,p3) = (X_tot,r_eff,v_eff)\n' &
         //'  p1,p2,p3 :                          mpsd=1: (p1,p2,p3) = (X_tot,r_med,S) for psdtyp=0\n'&
         //'           :                                  (p1,p2,p3) = (X_tot,gam,mu)  for psdtyp=1\n'&
         //'  g1,g2,g3 : particle geometry, (g1,g2,g3) = (D_max,A,V)'                   
  narg = 15
  call getCmdArgs(narg, argv, argmsg)
  read (argv(1), *) psdtyp
  read (argv(2), *) mgeo
  read (argv(3), *) msiz
  read (argv(4), *) mlog
  read (argv(5), *) nsiz
  read (argv(6), *) sizmin
  read (argv(7), *) sizmax
  read (argv(8), *) sizlog
  read (argv(9), *) mpsd
  read (argv(10),*) p(1)
  read (argv(11),*) p(2)
  read (argv(12),*) p(3)
  read (argv(13),*) g(1)
  read (argv(14),*) g(2)
  read (argv(15),*) g(3)
  allocate(siz(nsiz), psd(nsiz))
  

  ! particle geometry
  call GeNPSD__geometry(g(1), g(2), g(3), psi, rge)
  !psi: effective sphericty
  !rge: Vr/Ar = the r_e/r_g ratio


  ! size bins
  siz(:) = GeNPSD__sizbin(nsiz, sizmin, sizmax, sizlog)

 
  ! PSD generation
  alpha = GeNPSD__alpha(3, msiz, psi, rge)
  if (psdtyp == 0) then ! log-normal

     ! Size descriptor (p) conversion
     if (mpsd == 0) then
        lns  = sqrt(log(p(3)+1.0_R_))         ! ln(S) where S is the geometric stddev  
        pmed = alpha*p(2)*exp(-2.5_R_*lns**2) ! median radius

        ! physical variable
        if      (mgeo == 1) then ! A
           pmed = pmed * exp(2.0_R_*lns**2)
        else if (mgeo == 2) then ! V
           pmed = pmed * exp(3.0_R_*lns**2)
        end if
     else
        pmed = p(2)
        lns  = log(p(3))
     end if
    
     psd(:) = (1.0_R_/siz(:))**(1.0_R_-real(mlog)) * exp(-0.5_R_*(log(siz(:)/pmed)/lns)**2)/(sqrt(2.0_R_*PI_)*lns)
  else 

     ! Size descriptor (p) conversion
     if (mpsd == 0) then
        gam = alpha/(p(2)*p(3))    ! slope parameter
        mu  = 1.0_R_/p(3) - 2.0_R_ ! dispersion parameter 

        ! physical variable
        if      (mgeo == 1) then ! A
           mu = mu + 2.0_R_
        else if (mgeo == 2) then ! V
           mu = mu + 3.0_R_
        end if
     else
        gam = p(2)
        mu  = p(3)
     end if
     psd(:) = (siz(:)**(mu-1.0_R_+real(mlog)))*exp(-gam*siz(:))
  end if

  
  ! normalization
  if (mlog == 0) then
     fac = 0.5_R_ * sum((psd(1:nsiz-1)+psd(2:nsiz)) * abs( siz(2:nsiz)-siz(1:nsiz-1)) )
  else
     fac = 0.5_R_ * sum((psd(1:nsiz-1)+psd(2:nsiz)) * abs(log(siz(2:nsiz)/siz(1:nsiz-1))))
  end if
  psd(:) = psd(:)*p(1)/fac


  ! output
  do isiz = 1, nsiz   
     print*, siz(isiz), psd(isiz)
  end do
  deallocate(siz, psd)

end program genpsd_generate

