!+
! Generalized Nonspherical Particle Size Distribution (GeNPSD) 
! 
!   Drow a particle size distribution of nonspherical particles
!   based on analytical expressions 
!
!   Update: 
!     Masa Saito  03/17/2022 Created an initial version
!     Masa Saito  11/15/2022 Completed a published version
! 
!-
module genpsd

  use globals
  use hparx, only :  intpTab_linear_1R  
  implicit none
  private
  public :: GeNPSD__alpha
  public :: GeNPSD__interconv_ArbitraryPSD
  public :: GeNPSD__interconv_GammaPSD
  public :: GeNPSD__interconv_LogNormalPSD
  public :: GeNPSD__radius2diameter
  public :: GeNPSD__geometry
  public :: GeNPSD__sizbin


contains

  !+
  ! alpha coefficient that also includes radius-to-diameter conversion
  ! coefficients
  ! Table 1. in Saito and Yang (2022) J. Atmos. Sci.
  !- 
  function GeNPSD__alpha(msizP, msizQ, psi, rge) result(res)

    integer, intent(in)  :: msizP, msizQ
    real(R_), intent(in) :: psi
    real(R_), intent(in) :: rge
    real(R_) :: res


    ! size descriptor P (p --> e)
    if      (mod(msizP,4) == 0) then ! geometric size
       res = rge
    else if (mod(msizP,4) == 1) then ! projected area eqv size
       res = psi**(1.5_R_)
    else if (mod(msizP,4) == 2) then ! volume eqv size
       res = psi
    else if (mod(msizP,4) == 3) then ! effective size
       res = 1.0_R_
    end if

    ! size descriptor Q (e --> q)
    if      (mod(msizQ,4) == 0) then ! geometric size
       res = res/rge
    else if (mod(msizQ,4) == 1) then ! projected area eqv size
       res = res/psi**(1.5_R_)
    else if (mod(msizQ,4) == 2) then ! volume eqv size
       res = res/psi
    else if (mod(msizQ,4) == 3) then ! effective size
       res = res
    end if

    ! size descriptor radius/diameter
     res  = res * GeNPSD__radius2diameter(msizP, msizQ)   

  end function GeNPSD__alpha


  !+
  ! interconversion of arbitrary PSD
  ! Eq. (51) and Table 3 in Saito and Yang (2022) J. Atmos. Sci.
  !-
  function GeNPSD__interconv_ArbitraryPSD(mgeoX, mgeoY, msizP, msizQ, mlogP, mlogQ, mnor, psi, rge, sizP, psdX) result(psdY)
  
    integer, intent(in)  :: mgeoX, mgeoY, msizP, msizQ, mlogP, mlogQ, mnor
    real(R_), intent(in) :: psi, rge
    real(R_), intent(in) :: sizP(:)
    real(R_), intent(in) :: psdX(:)
    integer              :: nsiz
    real(R_)             :: alpha, xtot0, xtot, ytot0, ytot
    real(R_)             :: wrk1(size(sizP)), wrk2(size(sizP)), psdY(size(sizP)), sizQ(size(sizP)), dydx(size(sizP))

    ! dp/dq
    nsiz = size(sizP)
    alpha = GeNPSD__alpha(msizP, msizQ, psi, rge) !alpha_pq
    sizQ(:) = sizP(:) * alpha


    ! dy/dx|p
    dydx(:) = GeNPSD__dydx(msizP, mgeoX, mgeoY, sizP, psi, rge)


    ! PSD conversion dX/dlnp --> dX/dp
    if (mlogP == 0) then
       wrk1(:) = psdX(:)
       xtot  = 0.5_R_ * sum((wrk1(1:nsiz-1)+wrk1(2:nsiz)) * abs(    sizP(2:nsiz)-sizP(1:nsiz-1)) )
    else
       xtot0 = 0.5_R_ * sum((psdX(1:nsiz-1)+psdX(2:nsiz)) * abs(log(sizP(2:nsiz)/sizP(1:nsiz-1))))
       wrk1(:) = psdX(:) * (1.0_R_/sizP(:))**mlogP ! conversion: dX/dlnp --> dX/dp
       xtot  = 0.5_R_ * sum((wrk1(1:nsiz-1)+wrk1(2:nsiz)) * abs(    sizP(2:nsiz)-sizP(1:nsiz-1)) )
       wrk1(:) = wrk1(:)*xtot0/xtot  ! normalized to avoid numerical errors
    end if


    ! PSD conversion dX/dp --> dY/dq  [= dX/dp*(dY/dX)*(dp/dq), Eq. 41]
    wrk1(:) =  wrk1(:) * dydx(:) / alpha ! dp/dq = (alpha_pq)^1
    wrk2(:) = max(0.0_R_, intpTab_linear_1R(sizQ, wrk1, sizP, 1)) 
    !//NOTE(03/22/2022): This part may induce numerical errors when the size resolution is not sufficently fine.   


    ! PSD conversion dY/dq --> dY/dlnq
    if (mlogQ == 0) then
       psdY(:) = wrk2(:)
       ytot0 = 0.5_R_ * sum((psdY(1:nsiz-1)+psdY(2:nsiz)) * abs(    sizP(2:nsiz)-sizP(1:nsiz-1)) )
    else
       ytot0 = 0.5_R_ * sum((wrk2(1:nsiz-1)+wrk2(2:nsiz)) * abs(    sizP(2:nsiz)-sizP(1:nsiz-1)) )
       psdY(:)  = wrk2(:) * (sizP(:)*alpha)**mlogQ ! conversion: dY/dq --> dY/dlnq
       ytot  = 0.5_R_ * sum((psdY(1:nsiz-1)+psdY(2:nsiz)) * abs(log(sizP(2:nsiz)/sizP(1:nsiz-1))))
       psdY(:) = psdY(:)*ytot0/ytot ! normalized to avoid numerical errors
    end if

 
    ! normalization
    if (mnor == 1) then
       psdY(:) = psdY(:) / ytot0
    end if
 
  end function GeNPSD__interconv_ArbitraryPSD
  

  !+
  ! interconversion of Gamma PSD
  ! Eq. (42) and Table 2 in Saito and Yang (2022) J. Atmos. Sci.
  !-
  function GeNPSD__interconv_GammaPSD(mgeoX, mgeoY, msizP, msizQ, mlogQ, mnor, psi, rge, siz, &
                                      xtot, gamP, muX) result(res)
  
    integer,  intent(in)  :: mgeoX, mgeoY, msizP, msizQ, mlogQ, mnor
    real(R_), intent(in)  :: psi, rge, siz(:)
    real(R_), intent(in)  :: xtot, gamP, muX
    real(R_)              :: ytot, gamQ, muY, ytot0, wrk(size(siz)), res(size(siz))
    integer               :: nsiz, i
    real(R_)              :: alpha_pq, alpha_xy, lfac
    real(R_), parameter   :: odr(0:2)  = (/ 0.0_R_, 2.0_R_, 3.0_R_ /) 
    real(R_), parameter   :: coef(0:2) = (/ 1.0_R_, PI_, PI_*1.333333333333_R_ /)

     
    ! size descriptor conversion P --> Q
    lfac = odr(mgeoY)-odr(mgeoX)
    alpha_pq = GeNPSD__alpha(msizP, msizQ, psi, rge)
    gamQ  = gamP / alpha_pq
    nsiz = size(siz)
    

    ! physical variable conversion X --> Y
    if      (abs(lfac) == 0) then ! alpha_pp
       alpha_xy = 1.0_R_
    else if (abs(lfac) == 1) then ! alpha_pe
       alpha_xy = GeNPSD__alpha(msizP, 3, psi, rge)  
    else if (abs(lfac) == 2) then ! alpha_pa
       alpha_xy = GeNPSD__alpha(msizP, 1, psi, rge)  
    else if (abs(lfac) == 3) then ! alpha_pv
       alpha_xy = GeNPSD__alpha(msizP, 2, psi, rge)  
    end if
    muY = muX + lfac  
    ytot = 1.0_R_
    if (lfac > 0 ) then
       do i = 0, abs(int(lfac))-1
          ytot = ytot * (muX + real(i))
       end do
    else if (lfac < 0 ) then
       do i = 1, abs(int(lfac))
          ytot = ytot / (muX - real(i))
       end do
    end if
    ytot = xtot * ytot * coef(mgeoY)/coef(mgeoX) * (alpha_xy/gamP)**lfac   


    ! PSD conversion
    wrk(:) = (siz(:)**(muY-1.0_R_+real(mlogQ)))*exp(-gamQ*siz(:))
    if (mlogQ == 0) then
       ytot0 = 0.5_R_ * sum((wrk(1:nsiz-1)+wrk(2:nsiz)) * abs(    siz(2:nsiz)-siz(1:nsiz-1)) )
    else
       ytot0 = 0.5_R_ * sum((wrk(1:nsiz-1)+wrk(2:nsiz)) * abs(log(siz(2:nsiz)/siz(1:nsiz-1))))
    end if


    ! normalization
    if (mnor == 1) then 
       res(:) = wrk(:)/ytot0 
    else  
       res(:) = wrk(:)*ytot/ytot0 
    end if

  end function GeNPSD__interconv_GammaPSD


  !+
  ! interconversion of Log-Normal PSD
  ! Eq. (41) in Saito and Yang (2022) J. Atmos. Sci.
  !-
  function GeNPSD__interconv_LogNormalPSD(mgeoX, mgeoY, msizP, msizQ, mlogQ, mnor, psi, rge, &
                                          siz, xtot, pmed, lns) result(res)
                                              
    integer,  intent(in)  :: mgeoX, mgeoY, msizP, msizQ, mlogQ, mnor
    real(R_), intent(in)  :: psi, rge, siz(:)
    real(R_), intent(in)  :: xtot, pmed, lns
    real(R_)              :: ytot, qmed, ytot0, wrk(size(siz)), res(size(siz))
    integer               :: nsiz
    real(R_)              :: alpha_pq, alpha_xy, lfac
    real(R_), parameter   :: odr(0:2)  = (/ 0.0_R_, 2.0_R_, 3.0_R_ /) 
    real(R_), parameter   :: coef(0:2) = (/ 1.0_R_, PI_, PI_*1.333333333333_R_ /)
    
     
    ! size descriptor conversion P --> Q
    lfac = odr(mgeoY)-odr(mgeoX)
    alpha_pq = GeNPSD__alpha(msizP, msizQ, psi, rge)
    qmed  = pmed * alpha_pq * exp(lfac*lns**2.0_R_)
    nsiz = size(siz)


    ! physical variable conversion X --> Y
    if      (abs(lfac) == 0) then ! alpha_pp
       alpha_xy = 1.0_R_
    else if (abs(lfac) == 1) then ! alpha_pe
       alpha_xy = GeNPSD__alpha(msizP, 3, psi, rge)  
    else if (abs(lfac) == 2) then ! alpha_pa
       alpha_xy = GeNPSD__alpha(msizP, 1, psi, rge)  
    else if (abs(lfac) == 3) then ! alpha_pv
       alpha_xy = GeNPSD__alpha(msizP, 2, psi, rge)  
    end if
    ytot = xtot * (alpha_xy*pmed)**lfac * coef(mgeoY)/coef(mgeoX) * exp(0.5_R_*(lfac*lns)**2.0_R_)  


    ! PSD conversion
    wrk(:) = (1.0_R_/siz(:))**(1.0_R_-real(mlogQ)) * exp(-0.5_R_*(log(siz(:)/qmed)/lns)**2.0_R_)/(sqrt(2.0_R_*PI_)*lns)
    if (mlogQ == 0) then
       ytot0 = 0.5_R_ * sum((wrk(1:nsiz-1)+wrk(2:nsiz)) * abs(    siz(2:nsiz)-siz(1:nsiz-1)) )
    else
       ytot0 = 0.5_R_ * sum((wrk(1:nsiz-1)+wrk(2:nsiz)) * abs(log(siz(2:nsiz)/siz(1:nsiz-1))))
    end if


    ! normalization
    if (mnor == 1) then
       res(:) = wrk(:)/ytot0 
    else  
       res(:) = wrk(:)*ytot/ytot0 
    end if

  end function GeNPSD__interconv_LogNormalPSD


  !+
  ! dy/dx
  !-
  function GeNPSD__dydx(msizP, mgeoX, mgeoY, sizP, psi, rge) result(dydx)

    integer, intent(in)  :: mgeoX, mgeoY, msizP
    real(R_), intent(in) :: psi, rge
    real(R_), intent(in) :: sizP(:)
    real(R_) :: dydx(size(sizP)), alpha_xy, lfac
    real(R_), parameter   :: odr(0:2)  = (/ 0.0_R_, 2.0_R_, 3.0_R_ /) 
    real(R_), parameter   :: coef(0:2) = (/ 1.0_R_, PI_, PI_*1.333333333333_R_ /)
    
     
    ! Physical variable conversion X --> Y
    lfac = odr(mgeoY)-odr(mgeoX)
    if      (abs(lfac) == 0) then ! alpha_pp
       alpha_xy = 1.0_R_
    else if (abs(lfac) == 1) then ! alpha_pe
       alpha_xy = GeNPSD__alpha(msizP, 3, psi, rge)
    else if (abs(lfac) == 2) then ! alpha_pa
       alpha_xy = GeNPSD__alpha(msizP, 1, psi, rge)
    else if (abs(lfac) == 3) then ! alpha_pv
       alpha_xy = GeNPSD__alpha(msizP, 2, psi, rge)
    end if
    dydx(:) = (alpha_xy*sizP(:))**lfac * coef(mgeoY)/coef(mgeoX)    


  end function GeNPSD__dydx


  !+
  ! radius-to-diameter conversion, dr/dD or dD/dr
  ! Eqs. (43-45) in Saito and Yang (2022) J. Atmos. Sci.
  !- 
  function GeNPSD__radius2diameter(msizP, msizQ) result(res)

    integer, intent(in)  :: msizP, msizQ
    real(R_) :: res

    ! radius-to-diameter conversion
    ! radius   --> diameter: res = 2.0
    ! diameter --> radius:   res = 0.5
    ! otherwise:             res = 1.0
    res = (2.0_R_)**(((msizQ-mod(msizQ,4))-(msizP-mod(msizP,4)))/4)  ! diameter

  end function GeNPSD__radius2diameter


  !+
  ! sphericity and Vr/Ar ratio 
  !-
  subroutine GeNPSD__geometry(dmax, are, vol, psi, rge)
   
    real(R_), intent(in)  :: dmax, are, vol ! Maximum diameter, projected-area, volume
    real(R_), intent(out) :: psi, rge       ! Effective degree of sphericity, Vr/Ar ratio

    ! particle geometry
    psi = (PI_)**(1.0_R_/3.0_R_) * (6.0_R_*vol)**(2.0_R_/3.0_R_)/(4.0_R_*are) ! effective sphericty
    rge = 1.5_R_*vol/(are*dmax)                                               ! the Vr/Ar = r_e/r_g ratio

  end subroutine GeNPSD__geometry

  
  !+
  ! size bin generation
  !-
  function GeNPSD__sizbin(nsiz, sizmin, sizmax, sizlog) result(res)
    
    integer,  intent(in) :: nsiz, sizlog
    real(R_), intent(in) :: sizmin, sizmax
    integer  :: isiz
    real(R_) :: dsiz, res(nsiz)

    ! generate size bins
    if (sizlog == 0) then
       dsiz = (sizmax-sizmin)/real(nsiz-1)
       do isiz = 1, nsiz
          res(isiz) = sizmin + real(isiz-1)*dsiz
       end do
    else
       dsiz = log(sizmax/sizmin)/real(nsiz-1)
       do isiz = 1, nsiz
          res(isiz) = exp(log(sizmin) + real(isiz-1)*dsiz)
       end do
    end if

  end function GeNPSD__sizbin

end module genpsd

