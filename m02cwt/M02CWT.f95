!--------------------------------------------------------------!
! 'M02CWT' Module contains
! (1) 'rM02CWT_zerobc' subroutine,
!      which use the real-valued wavelet (CBSW, GDW or CWGDW) to
!      perform CWT of the signal with zero boundary condition. 
! (2) 'rM02CWT_periodbc' subroutine
!      which use the real-valued wavelet (CBSW, GDW or CWGDW) to
!      perform CWT of the signal with periodic boundary condition.
! (3) 'cM02CWT_zerobc' subroutine
!      which use the complex-valued wavelet (MW) 
!      to perform CWT of the signal with zero boundary condition.
! (4) 'cM02CWT_periodbc' subroutine
!      which use the complex-valued wavelet (MW) to perform 
!      CWT of the signal with periodic boundary condition.
!--------------------------------------------------------------!
Module M02CWT
  Use Parameters
  Use Tools
  Implicit None
  ! Output
  Real(8), Allocatable, Dimension(:) :: scales
  Real(8), Allocatable, Dimension(:,:) :: CWT
  Real(8), Allocatable, Dimension(:,:) :: realCWT
  Real(8), Allocatable, Dimension(:,:) :: imagCWT

  Contains
  !--------------------------------------------------------------!
  ! The real Continuous Wavelet Transform with 
  ! zero boundary condition
  !--------------------------------------------------------------!
  Subroutine rM02CWT_zerobc( signal, Lbox, Nsubs, wavelet_name, Nmesh )

    ! input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    ! f2py intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    ! f2py intent(in) :: Lbox
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
    ! f2py intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name
    ! f2py intent(in) :: wavelet_name

    ! Local variables
    Integer :: Nscales
    Integer :: Nd, Levl_min, Levl_max, Nlevs, l1, l2
    Integer :: b, i, j, k, l, m, startpt, endpt
    Real(8) :: cw, Dx, chi, w_0, w_min, w_max, gl, vx, N_4, N_2i, N_2i1
    Real(8) :: Nsubs_recip, scale_power, hh, xt, BScoord, BSval
    Real(8) :: ww, ss, lt, lt1, lt2 
    Real(8), Dimension(1:4) :: Cons
    Real(8), Dimension(-6:Nmesh+1) :: coeffs_zerobc
    Real(8), Allocatable, Dimension(:) :: scales_dmls
    Real(8), Allocatable, Dimension(:) :: pk
    Real(8), Allocatable, Dimension(:) :: coeffs1, coeffs4
    Character(len=80) :: CBSW_set, GDW_set, CWGDW_set
    
    ! The sampling interval of the signal
    Dx = Lbox/Real(Nmesh,8)

    ! Choose a wavelet
    CBSW_set  = "cbsw, CBSW, cubic BSW, cubic B-spline wavelet"
    GDW_set   = "gdw, GDW, Gaussian-derived wavelet, Gaussian derived wavelet"
    CWGDW_set = "cwgdw, CWGDW, cosine-weighted GDW, cosine-weighted Gaussian derived wavelet"
    If ( index(CBSW_set, wavelet_name)/= 0 ) Then 
      cw       = 0.466094761079290d0      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_cbsw
      Nd       = Nd_cbsw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk(-Nd:Nd+4) )
      pk = pcoeffs( wavelet_name, chi, Nd )
    Else if ( index(GDW_set, wavelet_name)/= 0 ) Then
      cw       = 2.0d0/Sqrt(5.0d0)      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_gdw
      Nd       = Nd_gdw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk(-Nd:Nd+4) )
      pk = pcoeffs( wavelet_name, chi, Nd )
    Else if ( index(CWGDW_set, wavelet_name)/= 0 ) Then
      cw       = 0.428218886729052d0      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_cwgdw
      Nd       = Nd_cwgdw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk(-Nd:Nd+4) )
      pk = pcoeffs( wavelet_name, chi, Nd )
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    ! Coefficients of the cubic B-spline approximation of the signal
    Call InterpCoeffs_zerobc( signal, Nmesh, coeffs_zerobc )
    Allocate( coeffs4(-6:Nmesh+1)  )
    coeffs4 = coeffs_zerobc
    Do i = 1, 4
      coeffs4 = CumSum( coeffs4 )
      Cons(i) = coeffs4(Nmesh+1)
    End Do

    ! Allocate arrays
    Allocate( CWT(0:Nscales-1,0:Nmesh-1) )

    ! Perform CWT at scales < w_0*2
    CWT(:,:) = 0.0d0
    Do b = 0, Nmesh-1
      Do j = 0, (1-Levl_min)*Nsubs-1
        ss  = scales_dmls(j)
        Do k = -Nd, Nd+4
          xt = b+(2-k)*hh/ss
          l1 = Ceiling(xt-6)
          l2 = l1+7
          vx = 0.0d0
          Do l = l1, l2
            BScoord = xt-2.0d0-l
            BSval   = BS7th( BScoord )
            If ( l<-6 ) Then 
              gl = 0.0d0
            Else if ( (l>=-6).And.(l<=Nmesh+1) ) Then
              gl = coeffs4(l)
            Else if ( l>Nmesh+1 ) Then
              lt  = Real(l-(Nmesh+1), 8)
              lt1 = lt*(lt+1)/2.0d0
              lt2 = lt*(lt+1)*(lt+2)/6.0d0
              gl  = Cons(4) + lt*Cons(3) + lt1*Cons(2) + lt2*Cons(1)
            End If
            vx = vx + gl*BSval
          End Do
          CWT(j,b) = CWT(j,b) + pk(k)*vx
        End Do
      End Do
    End Do
    Deallocate( coeffs4 )
    
    ! Allocate arrays
    N_4 = 0.25*Nmesh
    Allocate( coeffs1(Floor(-N_4)-6:Nmesh+2+Ceiling(N_4)) )
    coeffs1(Floor(-N_4)-6:-7)             = 0.0d0
    coeffs1(-6:Nmesh+1)                   = coeffs_zerobc(:)
    coeffs1(Nmesh+2:Nmesh+2+Ceiling(N_4)) = 0.0d0

    ! Perform CWT scales >= w_0*2
    Do i = 1, Levl_max
      N_2i = Real(Nmesh,8)/2**i
      N_2i1= Real(Nmesh,8)/2**(i+1)
      Do m = 0, 2**i-1
        startpt = Floor(m*N_2i-N_2i1)
        endpt   = Floor((m+1)*N_2i+N_2i1)
        Allocate( coeffs4(startpt-6:endpt+2) )
        coeffs4(:) = CumSum4( coeffs1(startpt-6:endpt+2) )
        Do b = Floor(m*N_2i), Floor((m+1)*N_2i)-1
          Do j = (i-Levl_min)*Nsubs, (i+1-Levl_min)*Nsubs-1
            ss = scales_dmls(j)
            Do k = -Nd, Nd+4
              xt = b+(2-k)*hh/ss
              l1 = Ceiling(xt-6)
              l2 = l1+7
              vx = 0.0d0
              Do l = l1, l2
                BScoord = xt-2.0d0-l
                BSval   = BS7th( BScoord )
                gl      = coeffs4(l)
                vx      = vx + gl*BSval
              End Do
              CWT(j,b) = CWT(j,b) + pk(k)*vx
            End Do
          End Do
        End Do
        Deallocate( coeffs4 )
      End Do
    End Do 

    ! Multiply scale
    Do i = 0, Nscales-1
      ww       = scales(i)
      ss       = scales_dmls(i)
      CWT(i,:) = Sqrt( ww )*(ss**3)*CWT(i,:)
    End Do
    Dx = Dx/hh**3
    CWT(:,:) = Dx*CWT(:,:)

  End Subroutine rM02CWT_zerobc

  !--------------------------------------------------------------!
  ! The real Continuous Wavelet Transform with 
  ! periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine rM02CWT_periodbc( signal, Lbox, Nsubs, wavelet_name, Nmesh )

    ! input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    ! f2py intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    ! f2py intent(in) :: Lbox
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
    ! f2py intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name
    ! f2py intent(in) :: wavelet_name

    ! Local variables
    Integer :: Nscales
    Integer :: Nd, Levl_min, Levl_max, Nlevs, l1, l2
    Integer :: b, i, j, k, l, ll, m, startpt, endpt
    Real(8) :: cw, Dx, chi, w_0, w_min, w_max, gl, vx, N_4, N_2i, N_2i1
    Real(8) :: Nsubs_recip, scale_power, hh, xt, BScoord, BSval
    Real(8) :: ww, ss
    Real(8), Dimension(0:Nmesh-1) :: coeffs_periodbc
    Real(8), Allocatable, Dimension(:) :: scales_dmls
    Real(8), Allocatable, Dimension(:) :: pk
    Real(8), Allocatable, Dimension(:) :: coeffs1, coeffs4
    Character(len=80) :: CBSW_set, GDW_set, CWGDW_set
    
    ! The sampling interval of the signal
    Dx = Lbox/Real(Nmesh,8)

    ! Choose a wavelet
    CBSW_set  = "cbsw, CBSW, cubic BSW, cubic B-spline wavelet"
    GDW_set   = "gdw, GDW, Gaussian-derived wavelet, Gaussian derived wavelet"
    CWGDW_set = "cwgdw, CWGDW, cosine-weighted GDW, cosine-weighted Gaussian derived wavelet"
    If ( index(CBSW_set, wavelet_name)/= 0 ) Then 
      cw       = 0.466094761079290d0      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_cbsw
      Nd       = Nd_cbsw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk(-Nd:Nd+4) )
      pk = pcoeffs( wavelet_name, chi, Nd )
    Else if ( index(GDW_set, wavelet_name)/= 0 ) Then
      cw       = 2.0d0/Sqrt(5.0d0)      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_gdw
      Nd       = Nd_gdw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk(-Nd:Nd+4) )
      pk = pcoeffs( wavelet_name, chi, Nd )
    Else if ( index(CWGDW_set, wavelet_name)/= 0 ) Then
      cw        = 0.428218886729052d0      ! The relationship between Fourier freq and wavelet scale
      chi       = chi_cwgdw
      Nd        = Nd_cwgdw
      hh        = chi/(Nd+2)
      w_0       = 2.0d0*chi/Lbox   
      w_min     = cw*Pi/Lbox 
      w_max     = cw*Pi*Nmesh/Lbox
      Levl_min  = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max  = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs     = Levl_max - Levl_min + 1
      Nscales   = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk(-Nd:Nd+4) )
      pk = pcoeffs( wavelet_name, chi, Nd )
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    ! Coefficients of the cubic B-spline approximation of the signal
    Call InterpCoeffs_periodbc( signal, Nmesh, coeffs_periodbc )
    Allocate( coeffs4(0:Nmesh-1)  )
    coeffs4(:) = coeffs_periodbc(:)
    Do i = 1, 4
      coeffs4 = CumSum( coeffs4 )
      coeffs4 = coeffs4 - Sum( coeffs4 )/Nmesh
    End Do

    ! Allocate arrays
    Allocate( CWT(0:Nscales-1,0:Nmesh-1) )

    ! Perform CWT at scales < w_0*2
    CWT(:,:) = 0.0d0
    Do b = 0, Nmesh-1
      Do j = 0, (1-Levl_min)*Nsubs-1
        ss  = scales_dmls(j)
        Do k = -Nd, Nd+4
          xt = b+(2-k)*hh/ss
          l1 = Ceiling(xt-6)
          l2 = l1+7
          vx = 0.0d0
          Do l = l1, l2
            BScoord = xt-2.0d0-l
            BSval   = BS7th( BScoord )
            If ( (l>=0).And.(l<=Nmesh-1) ) Then
              gl = coeffs4(l)
            Else 
              ll = Modulo(l, Nmesh)
              gl = coeffs4(ll)
            End If
            vx = vx + gl*BSval
          End Do
          CWT(j,b) = CWT(j,b) + pk(k)*vx
        End Do
      End Do
    End Do
    Deallocate( coeffs4 )
    
    ! Allocate arrays
    N_4 = 0.25*Nmesh
    Allocate( coeffs1(Floor(-N_4)-6:Nmesh+2+Ceiling(N_4)) )
    coeffs1(:-1)       = coeffs_periodbc(Nmesh+Floor(-N_4)-6:Nmesh-1)
    coeffs1(0:Nmesh-1) = coeffs_periodbc(:)
    coeffs1(Nmesh:)    = coeffs_periodbc(0:2+Ceiling(N_4))

    ! Perform CWT scales >= w_0*2
    Do i = 1, Levl_max
      N_2i = Real(Nmesh,8)/2**i
      N_2i1= Real(Nmesh,8)/2**(i+1)
      Do m = 0, 2**i-1
        startpt = Floor(m*N_2i-N_2i1)
        endpt   = Floor((m+1)*N_2i+N_2i1)
        Allocate( coeffs4(startpt-6:endpt+2) )
        coeffs4(:) = CumSum4( coeffs1(startpt-6:endpt+2) )
        Do b = Floor(m*N_2i), Floor((m+1)*N_2i)-1
          Do j = (i-Levl_min)*Nsubs, (i+1-Levl_min)*Nsubs-1
            ss = scales_dmls(j)
            Do k = -Nd, Nd+4
              xt = b+(2-k)*hh/ss
              l1 = Ceiling(xt-6)
              l2 = l1+7
              vx = 0.0d0
              Do l = l1, l2
                BScoord = xt-2.0d0-l
                BSval   = BS7th( BScoord )
                gl      = coeffs4(l)
                vx      = vx + gl*BSval
              End Do
              CWT(j,b) = CWT(j,b) + pk(k)*vx
            End Do
          End Do
        End Do
        Deallocate( coeffs4 )
      End Do
    End Do 

    ! Multiply scale
    Do i = 0, Nscales-1
      ww       = scales(i)
      ss       = scales_dmls(i)
      CWT(i,:) = Sqrt( ww )*(ss**3)*CWT(i,:)
    End Do
    Dx = Dx/hh**3
    CWT(:,:) = Dx*CWT(:,:)

  End Subroutine rM02CWT_periodbc

  !--------------------------------------------------------------!
  ! The complex Continuous Wavelet Transform with 
  ! zero boundary condition
  !--------------------------------------------------------------!
  Subroutine cM02CWT_zerobc( signal, Lbox, Nsubs, wavelet_name, Nmesh )

    ! input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    ! f2py intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    ! f2py intent(in) :: Lbox
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
    ! f2py intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name
    ! f2py intent(in) :: wavelet_name

    ! Local variables
    Integer :: Nscales
    Integer :: Nd, Levl_min, Levl_max, Nlevs, l1, l2
    Integer :: b, i, j, k, l, m, startpt, endpt
    Real(8) :: cw, Dx, chi, w_0, w_min, w_max, gl, vx, N_4, N_2i, N_2i1
    Real(8) :: Nsubs_recip, scale_power, hh, xt, BScoord, BSval
    Real(8) :: ww, ss, lt, lt1, lt2 
    Real(8), Dimension(1:4) :: Cons
    Real(8), Dimension(-6:Nmesh+1) :: coeffs_zerobc
    Real(8), Allocatable, Dimension(:) :: scales_dmls
    Real(8), Allocatable, Dimension(:) :: pk_r, pk_i
    Real(8), Allocatable, Dimension(:) :: coeffs1, coeffs4
    Character(len=80) :: MW_set
    
    ! The sampling interval of the signal
    Dx = Lbox/Real(Nmesh,8)

    ! Choose a wavelet
    MW_set  = "mw, MW, MorletW, Morlet wavelet"
    If ( index(MW_set, wavelet_name)/= 0 ) Then 
      cw       = 0.242640671273266d0      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_mw
      Nd       = Nd_mw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk_r(-Nd:Nd+4), pk_i(-Nd:Nd+4) )
      pk_r = pcoeffs( "rmw", chi, Nd )
      pk_i = pcoeffs( "imw", chi, Nd )
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    ! Coefficients of the cubic B-spline approximation of the signal
    Call InterpCoeffs_zerobc( signal, Nmesh, coeffs_zerobc )
    Allocate( coeffs4(-6:Nmesh+1)  )
    coeffs4(:) = coeffs_zerobc(:)
    Do i = 1, 4
      coeffs4 = CumSum( coeffs4 )
      Cons(i) = coeffs4(Nmesh+1)
    End Do

    ! Allocate arrays
    Allocate( realCWT(0:Nscales-1,0:Nmesh-1), imagCWT(0:Nscales-1,0:Nmesh-1) )

    ! Perform CWT at scales < w_0*2
    realCWT(:,:) = 0.0d0
    imagCWT(:,:) = 0.0d0
    Do b = 0, Nmesh-1
      Do j = 0, (1-Levl_min)*Nsubs-1
        ss       = scales_dmls(j)
        Do k = -Nd, Nd+4
          xt = b+(2-k)*hh/ss
          l1 = Ceiling(xt-6)
          l2 = l1+7
          vx = 0.0d0
          Do l = l1, l2
            BScoord = xt-2.0d0-l
            BSval   = BS7th( BScoord )
            If ( l<-6 ) Then 
              gl = 0.0d0
            Else if ( (l>=-6).And.(l<=Nmesh+1) ) Then
              gl = coeffs4(l)
            Else if ( l>Nmesh+1 ) Then
              lt  = Real(l-(Nmesh+1), 8)
              lt1 = lt*(lt+1)/2.0d0
              lt2 = lt*(lt+1)*(lt+2)/6.0d0
              gl  = Cons(4) + lt*Cons(3) + lt1*Cons(2) + lt2*Cons(1)
            End If
            vx = vx + gl*BSval
          End Do
          realCWT(j,b) = realCWT(j,b) + pk_r(k)*vx
          imagCWT(j,b) = imagCWT(j,b) + pk_i(k)*vx
        End Do
      End Do
    End Do
    Deallocate( coeffs4 )
    
    ! Allocate arrays
    N_4 = 0.25*Nmesh
    Allocate( coeffs1(Floor(-N_4)-6:Nmesh+2+Ceiling(N_4)) )
    coeffs1(Floor(-N_4)-6:-7)             = 0.0d0
    coeffs1(-6:Nmesh+1)                   = coeffs_zerobc(:)
    coeffs1(Nmesh+2:Nmesh+2+Ceiling(N_4)) = 0.0d0

    ! Perform CWT scales >= w_0*2
    Do i = 1, Levl_max
      N_2i = Real(Nmesh,8)/2**i
      N_2i1= Real(Nmesh,8)/2**(i+1)
      Do m = 0, 2**i-1
        startpt = Floor(m*N_2i-N_2i1)
        endpt   = Floor((m+1)*N_2i+N_2i1)
        Allocate( coeffs4(startpt-6:endpt+2) )
        coeffs4(:) = CumSum4( coeffs1(startpt-6:endpt+2) )
        Do b = Floor(m*N_2i), Floor((m+1)*N_2i)-1
          Do j = (i-Levl_min)*Nsubs, (i+1-Levl_min)*Nsubs-1
            ss = scales_dmls(j)
            Do k = -Nd, Nd+4
              xt = b+(2-k)*hh/ss
              l1 = Ceiling(xt-6)
              l2 = l1+7
              vx = 0.0d0
              Do l = l1, l2
                BScoord = xt-2.0d0-l
                BSval   = BS7th( BScoord )
                gl      = coeffs4(l)
                vx      = vx + gl*BSval
              End Do
              realCWT(j,b) = realCWT(j,b) + pk_r(k)*vx
              imagCWT(j,b) = imagCWT(j,b) + pk_i(k)*vx
            End Do
          End Do
        End Do
        Deallocate( coeffs4 )
      End Do
    End Do 

    ! Multiply scale
    Do i = 0, Nscales-1
      ww           = scales(i)
      ss           = scales_dmls(i)
      realCWT(i,:) = Sqrt( ww )*(ss**3)*realCWT(i,:)
      imagCWT(i,:) = Sqrt( ww )*(ss**3)*imagCWT(i,:)
    End Do
    Dx = Dx/hh**3
    realCWT(:,:) = Dx*realCWT(:,:)
    imagCWT(:,:) = Dx*imagCWT(:,:)

  End Subroutine cM02CWT_zerobc

  !--------------------------------------------------------------!
  ! The complex Continuous Wavelet Transform with 
  ! periodic boundary condition
  !--------------------------------------------------------------!
  Subroutine cM02CWT_periodbc( signal, Lbox, Nsubs, wavelet_name, Nmesh )

    ! input variables
    Integer, Intent(in) :: Nmesh
    Integer, Intent(in) :: Nsubs
    ! f2py intent(in) :: Nsubs
    Real(8), Intent(in) :: Lbox
    ! f2py intent(in) :: Lbox
    Real(8), Dimension(0:Nmesh-1), Intent(in) :: signal
    ! f2py intent(in) :: signal
    Character(len=*), Intent(in) :: wavelet_name
    ! f2py intent(in) :: wavelet_name

    ! Local variables
    Integer :: Nscales
    Integer :: Nd, Levl_min, Levl_max, Nlevs, l1, l2
    Integer :: b, i, j, k, l, ll, m, startpt, endpt
    Real(8) :: cw, Dx, chi, w_0, w_min, w_max, gl, vx, N_4, N_2i, N_2i1
    Real(8) :: Nsubs_recip, scale_power, hh, xt, BScoord, BSval
    Real(8) :: ww, ss
    Real(8), Dimension(0:Nmesh-1) :: coeffs_periodbc
    Real(8), Allocatable, Dimension(:) :: scales_dmls
    Real(8), Allocatable, Dimension(:) :: pk_r, pk_i
    Real(8), Allocatable, Dimension(:) :: coeffs1, coeffs4
    Character(len=80) :: MW_set
    
    ! The sampling interval of the signal
    Dx = Lbox/Real(Nmesh,8)

    ! Choose a wavelet
    MW_set  = "mw, MW, MorletW, Morlet wavelet"
    If ( index(MW_set, wavelet_name)/= 0 ) Then 
      cw       = 0.242640671273266d0      ! The relationship between Fourier freq and wavelet scale
      chi      = chi_mw
      Nd       = Nd_mw
      hh       = chi/(Nd+2)
      w_0      = 2.0d0*chi/Lbox   
      w_min    = cw*Pi/Lbox 
      w_max    = cw*Pi*Nmesh/Lbox
      Levl_min = Nint( Log(w_min/w_0)/Log(2.0d0) )   ! The minimum scale level, which is negative 
      Levl_max = Nint( Log(w_max/w_0)/Log(2.0d0) )-1 ! The maximum scale level, which is positive
      Nlevs    = Levl_max - Levl_min + 1
      Nscales  = Nsubs*Nlevs
      ! Set scales
      Allocate( scales(0:Nscales-1), scales_dmls(0:Nscales-1) )
      Nsubs_recip = 1.0d0/Nsubs
      Do i = 0, Nscales-1
        scale_power = Real(i + Levl_min*Nsubs, 8)*Nsubs_recip
        scales(i)   = w_0*( 2.0d0**scale_power )
      End Do
      scales_dmls(:) = scales(:)*Dx 
      ! pk coefficients
      Allocate( pk_r(-Nd:Nd+4), pk_i(-Nd:Nd+4) )
      pk_r = pcoeffs( "rmw", chi, Nd )
      pk_i = pcoeffs( "imw", chi, Nd )
    Else
      Stop "There is no wavelet function called "//wavelet_name
    End If
    
    ! Coefficients of the cubic B-spline approximation of the signal
    Call InterpCoeffs_periodbc( signal, Nmesh, coeffs_periodbc )
    Allocate( coeffs4(0:Nmesh-1)  )
    coeffs4(:) = coeffs_periodbc(:)
    Do i = 1, 4
      coeffs4 = CumSum( coeffs4 )
      coeffs4 = coeffs4 - Sum( coeffs4 )/Nmesh
    End Do

    ! Allocate arrays
    Allocate( realCWT(0:Nscales-1,0:Nmesh-1), imagCWT(0:Nscales-1,0:Nmesh-1) )

    ! Perform CWT at scales < w_0*2
    realCWT(:,:) = 0.0d0
    imagCWT(:,:) = 0.0d0
    Do b = 0, Nmesh-1
      Do j = 0, (1-Levl_min)*Nsubs-1
        ss       = scales_dmls(j)
        Do k = -Nd, Nd+4
          xt = b+(2-k)*hh/ss
          l1 = Ceiling(xt-6)
          l2 = l1+7
          vx = 0.0d0
          Do l = l1, l2
            BScoord = xt-2.0d0-l
            BSval   = BS7th( BScoord )
            If ( (l>=0).And.(l<=Nmesh-1) ) Then
              gl = coeffs4(l)
            Else 
              ll = Modulo(l, Nmesh)
              gl = coeffs4(ll)
            End If
            vx = vx + gl*BSval
          End Do
          realCWT(j,b) = realCWT(j,b) + pk_r(k)*vx
          imagCWT(j,b) = imagCWT(j,b) + pk_i(k)*vx
        End Do
      End Do
    End Do
    Deallocate( coeffs4 )
    
    ! Allocate arrays
    N_4 = 0.25*Nmesh
    Allocate( coeffs1(Floor(-N_4)-6:Nmesh+2+Ceiling(N_4)) )
    coeffs1(:-1)       = coeffs_periodbc(Nmesh+Floor(-N_4)-6:Nmesh-1)
    coeffs1(0:Nmesh-1) = coeffs_periodbc(:)
    coeffs1(Nmesh:)    = coeffs_periodbc(0:2+Ceiling(N_4))

    ! Perform CWT scales >= w_0*2
    Do i = 1, Levl_max
      N_2i = Real(Nmesh,8)/2**i
      N_2i1= Real(Nmesh,8)/2**(i+1)
      Do m = 0, 2**i-1
        startpt = Floor(m*N_2i-N_2i1)
        endpt   = Floor((m+1)*N_2i+N_2i1)
        Allocate( coeffs4(startpt-6:endpt+2) )
        coeffs4(:) = CumSum4( coeffs1(startpt-6:endpt+2) )
        Do b = Floor(m*N_2i), Floor((m+1)*N_2i)-1
          Do j = (i-Levl_min)*Nsubs, (i+1-Levl_min)*Nsubs-1
            ss = scales_dmls(j)
            Do k = -Nd, Nd+4
              xt = b+(2-k)*hh/ss
              l1 = Ceiling(xt-6)
              l2 = l1+7
              vx = 0.0d0
              Do l = l1, l2
                BScoord = xt-2.0d0-l
                BSval   = BS7th( BScoord )
                gl      = coeffs4(l)
                vx      = vx + gl*BSval
              End Do
              realCWT(j,b) = realCWT(j,b) + pk_r(k)*vx
              imagCWT(j,b) = imagCWT(j,b) + pk_i(k)*vx
            End Do
          End Do
        End Do
        Deallocate( coeffs4 )
      End Do
    End Do 

    ! Multiply scale
    Do i = 0, Nscales-1
      ww           = scales(i)
      ss           = scales_dmls(i)
      realCWT(i,:) = Sqrt( ww )*(ss**3)*realCWT(i,:)
      imagCWT(i,:) = Sqrt( ww )*(ss**3)*imagCWT(i,:)
    End Do
    Dx = Dx/hh**3
    realCWT(:,:) = Dx*realCWT(:,:)
    imagCWT(:,:) = Dx*imagCWT(:,:)

  End Subroutine cM02CWT_periodbc

End Module M02CWT