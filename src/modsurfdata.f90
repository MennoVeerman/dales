!> \file modsurfdata.f90
!! Variable definitions and auxilary routines for the surface model

!>
!! Variable definitions and auxilary routines for surface model
!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Thijs Heus, MPI-M
!!  \todo Documentation
!!  \par Revision list
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modsurfdata

! implicit none

SAVE
  integer :: isurf        = -1            !<   Flag for surface parametrization

  ! Soil properties

  ! Domain-uniform properties
  integer, parameter  :: ksoilmax = 4       !<  Number of soil layers [-]

  real              :: lambdasat          !<  heat conductivity saturated soil [W/m/K]
  real              :: Ke                 !<  Kersten number [-]

  real, allocatable :: zsoil  (:)         !<  Height of bottom soil layer from surface [m]
  real, allocatable :: zsoilc (:)        !<  Height of center soil layer from surface [m]
  real, allocatable :: dzsoil (:)         !<  Depth of soil layer [m]
  real, allocatable :: dzsoilh(:)         !<  Depth of soil layer between center of layers [m]

  ! Spatially varying properties
  real, allocatable :: lambda  (:,:,:)    !<  Heat conductivity soil layer [W/m/K]
  real, allocatable :: lambdah (:,:,:)    !<  Heat conductivity soil layer half levels [W/m/K]
  real, allocatable :: lambdas (:,:,:)    !<  Soil moisture diffusivity soil layer
  real, allocatable :: lambdash(:,:,:)    !<  Soil moisture diffusivity soil half levels
  real, allocatable :: gammas  (:,:,:)    !<  Soil moisture conductivity soil layer
  real, allocatable :: gammash (:,:,:)    !<  Soil moisture conductivity soil half levels
  real, allocatable :: Dh      (:,:,:)    !<  Heat diffusivity
  real, allocatable :: phiw    (:,:,:)    !<  Water content soil matrix [-]
  real, allocatable :: phiwm   (:,:,:)    !<  Water content soil matrix previous time step [-]
  real, allocatable :: phifrac (:,:,:)    !<  Relative water content per layer [-]
  real              :: phiwav  (ksoilmax)
  real, allocatable :: phitot  (:,:)      !<  Total soil water content [-]
  real, allocatable :: pCs     (:,:,:)    !<  Volumetric heat capacity [J/m3/K]
  real, allocatable :: rootf   (:,:,:)    !<  Root fraction per soil layer [-]
  real              :: rootfav (ksoilmax)
  real, allocatable :: tsoil   (:,:,:)    !<  Soil temperature [K]
  real, allocatable :: tsoilm  (:,:,:)    !<  Soil temperature previous time step [K]
  real              :: tsoilav (ksoilmax)
  real, allocatable :: tsoildeep (:,:)    !<  Soil temperature [K]
  real              :: tsoildeepav

  real, allocatable :: swdavn  (:,:,:)
  real, allocatable :: swuavn  (:,:,:)
  real, allocatable :: lwdavn  (:,:,:)
  real, allocatable :: lwuavn  (:,:,:)

  integer           :: nradtime  = 60

  ! Soil related properties [case specific]
  real              :: phi       = 0.472  !<  volumetric soil porosity [-]
  real              :: phifc     = 0.323  !<  volumetric moisture at field capacity [-]
  real              :: phiwp     = 0.171  !<  volumetric moisture at wilting point [-]

  ! Soil related constants [adapted from ECMWF]
  real, parameter   :: pCm       = 2.19e6 !<  Volumetric soil heat capacity [J/m3/K]
  real, parameter   :: pCw       = 4.2e6  !<  Volumetric water heat capacity [J/m3/K]

  real, parameter   :: lambdadry = 0.190  !<  Heat conductivity dry soil [W/m/K]
  real, parameter   :: lambdasm  = 3.11   !<  Heat conductivity soil matrix [W/m/K]
  real, parameter   :: lambdaw   = 0.57   !<  Heat conductivity water [W/m/K]

  real, parameter   :: bc        = 6.04     !< Clapp and Hornberger non-dimensional exponent [-]
  real, parameter   :: gammasat  = 0.57e-6  !< Hydraulic conductivity at saturation [m s-1]
  real, parameter   :: psisat    = -0.388   !< Matrix potential at saturation [m]

  ! Land surface properties

  ! Surface properties
  real, allocatable :: z0m        (:,:) !<  Roughness length for momentum [m]
  real              :: z0mav    = -1
  real, allocatable :: z0h        (:,:) !<  Roughness length for heat [m]
  real              :: z0hav    = -1
  real, allocatable :: tskin_surf      (:,:) !< Surface  skin temperature [K]
  real, allocatable :: tskinm_surf     (:,:) !< Surface skin temperature previous timestep [K]
  real, allocatable :: Wl         (:,:) !<  Liquid water reservoir [m]
  real              :: Wlav     = -1
  real, parameter   :: Wmax     = 0.0002 !<  Maximum layer of liquid water on surface [m]
  real, allocatable :: Wlm        (:,:) !<  Liquid water reservoir previous timestep [m]
  real, allocatable :: qskin_surf (:,:) !<  Surface skin specific humidity [kg/kg]
  real, allocatable :: albedo_surf(:,:) !<  Surface albedo [-]
  real              :: albedoav_surf = -1
  real, allocatable :: LAI_surf     (:,:) !<  Leaf area index understory vegetation [-]
  real              :: LAI_surfav    = -1
  real, allocatable :: cveg       (:,:) !<  Vegetation cover [-]
  real              :: cvegav   = -1
  real, allocatable :: cliq       (:,:) !<  Fraction of vegetated surface covered with liquid water [-]
  real, allocatable :: Cskin      (:,:) !<  Heat capacity skin layer [J]
  real              :: Cskinav  = -1
  real, allocatable :: lambdaskin (:,:) !<  Heat conductivity skin layer [W/m/K]
  real              :: lambdaskinav = -1
  real              :: ps       = -1    !<  Surface pressure [Pa]

  ! AGS options (require interactive landsurface: isurf=2)
  !<Namelist options
  logical           :: lrsAgs        = .false.!<  Switch to apply AGS to calculate resistances
  logical           :: lCO2Ags       = .false.!<  Switch to calculate CO2 fluxes with AGS
  integer           :: planttype     = 3      !<  Integer to switch between (C)3 and (C)4 plants
  logical           :: lrelaxgc_surf = .false.!<  Switch to delay plant response at surface. Timescale is equal to 1/kgc_surf
  real              :: kgc_surf      = 0.00113!<  Standard stomatal response rate at surface(corresponding to a time scale of 14.75 min.) [1/s]
  logical           :: gcsurf_old_set = .false.!<  Only apply relaxing function after initial gc is calculated once for surface
  logical           :: lrelaxci_surf  = .false.!<  Switch to delay internal CO2 concentration in plant leafs; Timescale equal to that for gc
  real              :: kci_surf       = 0.00113!<  Standard internal CO2 response rate (corresponding to a time scale of 14.75 min.) [1/s]
  logical           :: cisurf_old_set = .false.!<  Only apply relaxing function after initial ci is calculated once
  logical           :: tskin_set      = .false.!< Switch needed if lcanopyeb, to use surface radiative ppties in initaliation
  real, allocatable :: gc_old       (:,:)  !<  Old value for gc
  real              :: gc_inf              !<  Attractor for stomatal response rate
  real, allocatable :: ci_old       (:,:)  !<  Old value for ci
  real              :: ci_inf              !<  Attractor for ci
  real              :: wco2av     = 0.0
  real              :: Anav       = 0.0
  real              :: gcco2av    = 0.0
  real              :: alb_canav  = 0.0
  real              :: Respav     = 0.0
  !<Non namelist options
  logical           :: linags     = .false.!<  Switch to make additional initialization for AGS
  logical           :: lCHon      = .false.!<  Equal to lchem, but due to compilation has to be outside modchem.f90
  integer           :: indCO2     = -1     !<  Index of CO2 in the scalars
  integer           :: CO2loc     = -1     !<  Index of CO2 in the scalars
  real, allocatable :: CO2flux(:,:)        !<  Surface flux of CO2 as calculated by AGS

  !AGS variables
  real              :: CO2comp298 =  62.93 !<  CO2 compensation concentration (default was 68.5)
  real              :: Q10CO2     =    2.0 !<  Parameter to calculate the CO2 compensation concentration (default was 1.5)
  real              :: gm298      =   2.65 !<  Mesophyll conductance at 298 K (default was 7.0)
  real              :: Q10gm      =    2.0 !<  Parameter to calculate the mesophyll conductance
  real              :: T1gm       =  278.0 !<  Reference temperature to calculate the mesophyll conductance
  real              :: T2gm       =  311.0 !<  Reference temperature to calculate the mesophyll conductance (default was 301)
  real              :: gmin       = 2.5e-4 !<  Cuticular (minimum) conductance
  real              :: nuco2q     =    1.6 !<  Ratio molecular viscosity water to carbon dioxide
  real              :: f0         =   0.89 !<  Maximum value Cfrac
  real              :: ad         =   0.07 !<  Regression coefficient to calculate Cfrac
  real              :: Ammax298   =  0.780 !<  CO2 maximal primary productivity (default was 2.2)
  real              :: Q10am      =    2.0 !<  Parameter to calculate maximal primary productivity
  real              :: T1Am       =    281 !<  Reference temperature to calculate maximal primary productivity
  real              :: T2Am       =    311 !<  Reference temperature to calculate maximal primary productivity
  real              :: alpha0     =  0.017 !<  Initial low light conditions
  real              :: Kx         =    0.7 !<  Extinction coefficient PAR
  real              :: Cw         = 1.6e-3 !<  Constant water stress correction
  real              :: wsmax      =   0.55 !<  Upper reference value soil water
  real              :: wsmin      =  0.005 !<  Lower reference value soil water
  real              :: R10        =   0.23 !<  Respiration at 10oC (Jacobs 2007)
  real              :: Eact0      = 53.3e3 !<  Activation energy
  real              :: MW_Air     =  28.97 !<  Molecular weight of air
  real              :: MW_CO2     =     44 !<  Molecular weight of CO2
  real, allocatable :: wco2Field    (:,:)
  real, allocatable :: AnField      (:,:)
  real, allocatable :: gcco2Field   (:,:)
  real, allocatable :: rsco2Field   (:,:)
  real, allocatable :: RespField    (:,:)
  real, allocatable :: fstrField    (:,:)
  real, allocatable :: tauField     (:,:)
  real, allocatable :: ciField      (:,:)
  real, allocatable :: PARField     (:,:)
  real, allocatable :: PARdirField  (:,:)
  real, allocatable :: PARdifField  (:,:)

  !Variables for 2leaf AGS
  logical                   :: lsplitleaf =                  .false. !<  Switch to split AGS calculations over different parts of the leaf (direct & diffuse at different layers)
  logical                   :: l3leaves   =                  .true. !<  Switch to keep 3 leaf orientation angles on sunny leaves
  integer,parameter         :: nz_gauss   =                        3 !<  Amount of bins to use for Gaussian integrations
  real, dimension(nz_gauss)     ::  LAI_g      = (/0.8873,   0.5,0.1127/) !<  Ratio of integrated LAI at locations where shaded leaves are evaluated in the second Gaussian integration
  !real, dimension(nz_gauss)     ::  LAI_g      = (/0.1127,   0.5,0.8873/) !<  Ratio of integrated LAI at locations where shaded leaves are evaluated in the second Gaussian integration
  integer,parameter         :: nangle_gauss =                      3 !< Amount of bins to use for Gaussian integrations on leaf angles
  !integer,parameter         :: nband_can = 4 !<Number of spectral bands for spectral integration
  !integer,parameter         :: iband_par = 2 !<Which spectral bands belongs to the PAR spectrum?
  !integer,parameter         :: iband_uvs_s = 1 !<first spectral band of UV/VIS bands for computing albedo
  !integer,parameter         :: iband_uvs_e = 2 !<last spectral band of UV/VIS bands for computing albedo
  !integer,parameter         :: iband_nir_s = 3 !<first spectral band of NIR bands for computing albedo
  !integer,parameter         :: iband_nir_e = 4 !<last spectral band of NIR bands for computing albedo
  real, dimension(nangle_gauss) :: weight_g   = (/0.2778,0.4444,0.2778/) !<  Weights of the Gaussian bins (must add up to 1)
  real, dimension(nangle_gauss) :: angle_g    = (/0.1127,   0.5,0.8873/) !<  Sines of the leaf angles compared to the sun in the first Gaussian integration
  !real, dimension(nband_can)  :: sigma_b = (/0.157, 0.092, 0.622, 0.750/) ! UV PAR FAR-RED NIR/) !<  Scattering coefficients of leaves-controls the effective albedo
  !real, dimension(nband_can)  :: weight_b = (/0.128, 0.450, 0.055, 0.367/) !UV PAR FAR-RED NIR/) !<  Fraction of top-of-canopy irradiance per spectral band

  type :: canrad_band
    integer    :: spectral_type   ! 1 (UV), 2 (PAR), 3 (NIR)
    real       :: sigma  ! scattering coefficient
    real       :: weight ! irradiance fraction
    real       :: alb_soil ! soil albedo
  end type canrad_band
  integer :: n_bands_sw
  type(canrad_band), allocatable :: canrad_bands_sw (:)

  !real                      :: kdfbl      =                      0.8 !<  Diffuse radiation extinction coefficient for black leaves
  real, allocatable         :: gshad_old    (:,:,:)
  real, allocatable         :: gleafsun_old (:,:,:,:)
  real, allocatable         :: gsun_old     (:,:,:)

! variables for 2leaf AGS
  real,allocatable ::   PARleaf_shad   (:)    !< PAR at shaded leaves per vertical level [W/m2]
  real,allocatable ::   PARleaf_sun    (:,:)  !< PAR at sunny leaves per vertical level and leaf orientation [W/m2]
  real,allocatable ::   PARleaf_allsun (:)    !< PAR at sunny leaves per vertical level averaged over all leaf orientations [W/m2]
  real,allocatable ::   swleaf_shad   (:,:)    !< SW at shaded leaves per vertical level [W/m2]
  real,allocatable ::   swleaf_sun    (:,:,:)  !< SW at sunny leaves per vertical level and leaf orientation [W/m2]
  real,allocatable ::   swleaf_allsun (:,:)    !< SW at sunny leaves per vertical level averaged over all leaf orientations [W/m2]
  real,allocatable ::   fSL            (:)    !< Fraction of sunlit leaves per evrtical level [-]
  real,allocatable ::   albdir_lsplit   (:,:,:)  !< effective direct SW albedo at vegetation top according to canopyrad
  real,allocatable ::   albdif_lsplit   (:,:,:)  !< effective diffuse SW albedo at vegetation top according to canopyrad
  real,allocatable ::   albswd_lsplit   (:,:,:)  !< effective global SW albedo at vegetation top according to canopyrad
  real,allocatable ::   swdir_lsplit   (:,:,:,:)!< swdir profile along vegetation
  real,allocatable ::   swdif_lsplit   (:,:,:,:)!< swdir profile along vegetation
  real,allocatable ::   swu_lsplit     (:,:,:,:)!< swu profile along vegetation
  real,allocatable ::   PARdir_lsplit   (:)!< PARdir profile along vegetation
  real,allocatable ::   PARdif_lsplit   (:)!< PARdif profile along vegetation
  real,allocatable ::   PARu_lsplit     (:)!< PAR profile along vegetation

  real         ::   swdir_TOV        !< direct shortwave at top of vegetation [W/m2]
  real         ::   swdif_TOV        !< diffuse shortwave at top of vegetation [W/m2]
  integer      ::   surfrad_meth = 2  !< (method to calculate radiation and absorbed  fluxes in surface vegetation (when lcanopyeb=false)  =1  XPB2017, =2 Goudriaan,V.Laar 1996
  logical      ::   ldiscr      = .true.  !< use discrete defintion of differentiation to calculate absorbed radiation


  ! Surface energy balance
  real, allocatable :: Qnet     (:,:)   !<  Net radiation [W/m2]
  real              :: Qnetav   = -1
  real, allocatable :: LE       (:,:)   !<  Latent heat flux [W/m2]
  real, allocatable :: H        (:,:)   !<  Sensible heat flux [W/m2]
  real, allocatable :: G0       (:,:)   !<  Ground heat flux [W/m2]
  real, allocatable :: ra       (:,:)   !<  Aerodynamic resistance [s/m]
  real, allocatable :: rs       (:,:)   !<  Composite resistance [s/m]
  real, allocatable :: rsveg    (:,:)   !<  Vegetation resistance [s/m]
  real, allocatable :: rssoil   (:,:)   !<  Soil evaporation resistance [s/m]
  real              :: rsisurf2 = 0.    !<  Vegetation resistance [s/m] if isurf2 is used
  real, allocatable :: rsmin    (:,:)   !<  Minimum vegetation resistance [s/m]
  real              :: rsminav = -1
  real, allocatable :: rssoilmin(:,:)   !<  Minimum soil evaporation resistance [s/m]
  real              :: rssoilminav = -1
  real, allocatable :: tendskin_surf (:,:)   !<  Tendency of skin [W/m2]
  real, allocatable :: gD       (:,:)   !<  Response factor vegetation to vapor pressure deficit [-]
  real              :: gDav

  ! Turbulent exchange variables
  logical           :: lmostlocal  = .false.  !<  Switch to apply MOST locally to get local Obukhov length
  logical           :: lsmoothflux = .false.  !<  Create uniform sensible and latent heat flux over domain
  logical           :: lneutral    = .false.  !<  Disable stability corrections
  real, allocatable :: obl   (:,:)      !<  Obukhov length [m]
  real              :: oblav            !<  Spatially averaged obukhov length [m]
  real, allocatable :: Cm    (:,:)      !<  Drag coefficient for momentum [-]
  real, allocatable :: Cs    (:,:)      !<  Drag coefficient for scalars [-]
  real, allocatable :: ustar (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlflux (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtflux  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: svflux  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  ! Surface gradients of prognostic variables
  real, allocatable :: dudz  (:,:)      !<  U-wind gradient in surface layer [1/s]
  real, allocatable :: dvdz  (:,:)      !<  V-wind gradient in surface layer [1/s]
  real, allocatable :: dqtdz (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: dthldz(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]

  ! Surface properties in case of prescribed conditions (previous isurf 2, 3 and 4)
  real              :: thls  = -1       !<  Surface liquid water potential temperature [K]
  real              :: qts              !<  Surface specific humidity [kg/kg]
  real              :: thvs             !<  Surface virtual temperature [K]
  real, allocatable :: svs   (:)        !<  Surface scalar concentration [-]
  real              :: z0    = -1       !<  Surface roughness length [m]

  ! prescribed surface fluxes
  real              :: ustin  = -1      !<  Prescribed friction velocity [m/s]
  real              :: wtsurf = -1e20  !<  Prescribed kinematic temperature flux [K m/s]
  real              :: wqsurf = -1e20  !<  Prescribed kinematic moisture flux [kg/kg m/s]
  real              :: wsvsurf(100) = 0 !<  Prescribed surface scalar(n) flux [- m/s]

  ! Heterogeneous surfaces
  integer, parameter:: max_lands                   = 10 !<  Amount of land types that can be defined
  integer, parameter:: mpatch                      = 16 !<  Maximum amount of patches that can be defined in each direction
  logical           :: lhetero                 = .false.!<  Switch to apply heterogeneous surfaces using surface.inp.xxx
  logical           :: loldtable               = .false.!<  Switch to use surface.inp.xxx instead of updated surface.<name>.inp.xxx
  integer           :: xpatches                    =  2 !<  Amount of patches in the x-direction
  integer           :: ypatches                    =  1 !<  Amount of patches in the y-direction
  integer           :: land_use(mpatch,mpatch)     =  0 !<  Indicator for the land type
  integer           :: landtype(max_lands)         = -1 !< Type nr of the land in surface.inp.xxx
  character(len=10),dimension(max_lands) :: landname = "none" !< Name of the land type in surface.inp.xxx
  real              :: z0mav_land(max_lands)       = -1 !< Roughness length per land type for momentum [m]
  real              :: z0hav_land(max_lands)       = -1 !< Roughness length per land type for heat [m]
  real, allocatable :: z0mav_patch(:,:)                 !< Rougness length per patch
  real, allocatable :: z0hav_patch(:,:)                 !< Rougness length per patch
  real              :: thls_land(max_lands)        = -1 !< Surface liquid water potential temperature [K]
  real, allocatable :: thls_patch(:,:)                  !< Surface liquid water potential temperature [K]
  real, allocatable :: qts_patch(:,:)                   !< Surface liquid water potential temperature [K]
  real, allocatable :: thvs_patch(:,:)                  !< Surface virtual potential temperature [K]
  real              :: ps_land(max_lands)          = -1 !< Surface pressure [Pa]
  real, allocatable :: ps_patch(:,:)                    !< Surface pressure [Pa]
  real              :: ustin_land(max_lands)       = -1 !< Prescribed friction velocity [m/s]
  real, allocatable :: ustin_patch(:,:)                 !< Prescribed friction velocity [m/s]
  real              :: wt_land(max_lands)          = -1 !< Prescribed kinematic temperature flux [K m/s]
  real              :: wq_land(max_lands)          = -1 !< Prescribed kinematic moisture flux [kg/kg m/s]
  real              :: wsv_land(100,max_lands)     = -1 !< Prescribed surface scalar(n) flux [- m/s]
  real, allocatable :: wt_patch(:,:)                    !< Prescribed kinematic temperature flux [K m/s]
  real, allocatable :: wq_patch(:,:)                    !< Prescribed kinematic moisture flux [kg/kg m/s]
  real, allocatable :: wsv_patch(:,:,:)                 !< Prescribed surface scalar(n) flux [- m/s]
  real, allocatable :: rsisurf2_patch(:,:)              !< Vegetation resistance [s/m] if isurf2 is used
  real              :: rsisurf2_land(max_lands)    = 0  !< Vegetation resistance [s/m] if isurf2 is used
  real, allocatable :: albedo_patch(:,:)                !< Albedo
  real              :: albedo_land(max_lands)      = -1 !< Albedo
  real, allocatable :: tsoil_patch(:,:,:)               !< Soil temperature [K]
  real, allocatable :: tsoildeep_patch(:,:)             !< Soil temperature [K]
  real, allocatable :: phiw_patch(:,:,:)                !<
  real, allocatable :: rootf_patch(:,:,:)               !< Root fraction per soil layer [-]
  real, allocatable :: Cskin_patch(:,:)                 !< Heat capacity skin layer [J]
  real, allocatable :: lambdaskin_patch(:,:)            !< Heat conductivity skin layer [W/m/K]
  real, allocatable :: Qnet_patch(:,:)                  !< Net radiation [W/m2]
  real, allocatable :: cveg_patch(:,:)                  !< Vegetation cover [-]
  real, allocatable :: Wl_patch(:,:)                    !< Liquid water reservoir [m]
  real, allocatable :: rsmin_patch(:,:)                 !< Minimum vegetation resistance [s/m]
  real, allocatable :: LAI_patch(:,:)                   !< Leaf area index vegetation [-]
  real, allocatable :: gD_patch(:,:)                    !< Response factor vegetation to vapor pressure deficit [-]
  real              :: tsoil_land(ksoilmax,max_lands)=-1!< Soil temperature [K]
  real              :: tsoildeep_land(max_lands)   = -1 !< Soil temperature [K]
  real              :: phiw_land(ksoilmax,max_lands) =-1!<
  real              :: rootf_land(ksoilmax,max_lands)=-1!< Root fraction per soil layer [-]
  real              :: Cskin_land(max_lands)       = -1 !< Heat capacity skin layer [J]
  real              :: lambdaskin_land(max_lands)  = -1 !< Heat conductivity skin layer [W/m/K]
  real              :: Qnet_land(max_lands)        = -1 !< Net radiation [W/m2]
  real              :: cveg_land(max_lands)        = -1 !< Vegetation cover [-]
  real              :: Wl_land(max_lands)          = -1 !< Liquid water reservoir [m]
  real              :: rsmin_land(max_lands)       = -1 !< Minimum vegetation resistance [s/m]
  real              :: LAI_land(max_lands)         = -1 !< Leaf area index vegetation [-]
  real              :: gD_land(max_lands)          = -1 !< Response factor vegetation to vapor pressure deficit [-]
  real, allocatable :: oblpatch(:,:)                    !<  Obukhov length [m]


contains

function rho_c_dif(sigma)
    real :: rho_c_dif
    real, intent(in) :: sigma

    rho_c_dif = 4.636239e-03*exp(19.707778*sigma-57.752878*sigma**2+100.713885*sigma**3-88.256337*sigma**4+30.761673*sigma**5)
    return
end function

function rho_c_dir(mu, rho_c_dif)
    real :: rho_c_dir
    real :: a,b,c
    real, intent(in) :: mu, rho_c_dif
    a = 1.677900*exp(-1.804072*rho_c_dif**0.935482)+0.347942
    b = 4.206250*exp(-0.534006*rho_c_dif**0.955533)-2.443320
    c = -0.366804*exp(-0.676143*rho_c_dif**0.983803)+1.107502
    rho_c_dir = rho_c_dif * (a*exp(-b*mu**c) + 0.4)
    return
end function

! Norman 1979 radiative transfer model, following tridiagonal matrix solution of Bonan 2019 (https://doi.org/10.1017/9781107339217.015)
subroutine canopyrad_sw_N79(layers,LA,tau_dif_can, PHIdir_TOC,PHIdif_TOC,clump, & ! in
                        Hshad,Hsun, fracSL,                              & ! out needed for vegetation
                        effalb_dir,effalb_dif,effalb_phi,phidircan,phidifcan,phiucan,isoil) ! out needed for radiation if lcanopyeb
  use modraddata , only : zenith
  use modglobal  , only : xtime,rtimee,xday,xlat,xlon
  implicit none

  integer,intent(in) :: layers
  real, intent(in),dimension(layers-1) :: LA    ! local leaf area (m2 leaf/m2 ground)
  real, intent(in)   :: PHIdir_TOC                  ! Direct  radiation at vegetation top, always >0
  real, intent(in)   :: PHIdif_TOC                  ! Diffuse radition at vegetation top, always >0
  real, intent(in),dimension(layers-1)   :: tau_dif_can
  real, intent(in)   :: clump                     ! leaf clumping index

  real, intent(out),dimension(layers)                        :: fracSL ! fraction of sunlit leaves per layer
  real, intent(out),dimension(layers-1,n_bands_sw)              :: Hshad  ! Intercepted radiation by shaded leaves between the layer and the one below
  real, intent(out),dimension(layers-1,nangle_gauss,n_bands_sw) :: Hsun   ! Intercepted radiation by sunlit leaves between the layer and the one below per layer and per leaf orientation
  real, intent(out),dimension(layers,n_bands_sw) :: phidircan          ! downwards direct comp. of radiation inside canopy
  real, intent(out),dimension(layers,n_bands_sw) :: phidifcan          ! downwards diffuse comp. of radiation inside canopy
  real, intent(out),dimension(layers,n_bands_sw) :: phiucan            ! upwards (diffuse) comp. of radition inside canopy
  real, intent(out),dimension(n_bands_sw) :: isoil                                 ! absorbed radiation by soil
  real, intent(out),dimension(n_bands_sw) ::  effalb_dir                           ! effective albedo of canopy top for direct radiation
  real, intent(out),dimension(n_bands_sw) ::  effalb_dif                           ! effective albedo of canopy top for diffuse radiation
  real, intent(out),dimension(n_bands_sw) ::  effalb_phi                            ! effective albedo of canopy top for total radiation

  real, dimension(2*layers) :: m_a, m_c, m_d ! matrix coefficients (m_b==1)
  real, dimension(2*layers) :: c_prime, d_prime
  real, dimension(layers-1) :: tau_beam
  real, dimension(layers) :: T_beam
  real, dimension(layers) :: phidifcan_nd, phiucan_nd
  real                    :: cf_a, cf_b, cf_c, cf_d ! temporary coefficients
  real                    :: sigma, weight, leaf_refl_trans,fracSL_m
  real                    :: cos_sza, alb, minsinbeta = 1.e-10

  integer :: k_can, ib

  cos_sza  = max(zenith(xtime*3600 + rtimee,xday,xlat,xlon), minsinbeta)
  if (cos_sza>0.035) then ! day, same threshold as radpar

    do ib = 1, n_bands_sw
      sigma = canrad_bands_sw(ib)%sigma
      weight = canrad_bands_sw(ib)%weight


      leaf_refl_trans = sigma / 2.
      alb = canrad_bands_sw(ib)%alb_soil

      ! Direct beam transmittance
      T_beam(layers) = 1
      do k_can=layers-1, 1, -1
        tau_beam(k_can) = exp(-0.5 / cos_sza * clump * LA(k_can))
        T_beam(k_can) = T_beam(k_can+1) * tau_beam(k_can)
      end do

      if (ib == 1) fracSL(:) = clump * T_beam(:)    ! Fraction of sun-lit leaves

      m_d(1) = PHIdir_TOC * weight * alb * T_beam(1) ! bottom boundary condition
      m_c(1) = -alb
      do k_can=1, layers-1
        cf_a = (1-tau_dif_can(k_can))*leaf_refl_trans - ((tau_dif_can(k_can) + (1-tau_dif_can(k_can))*leaf_refl_trans)**2) / ((1-tau_dif_can(k_can))*leaf_refl_trans)
        cf_b = (tau_dif_can(k_can) + (1-tau_dif_can(k_can))*leaf_refl_trans) / ((1-tau_dif_can(k_can))*leaf_refl_trans)
        m_a(2*k_can) = -cf_a
        m_c(2*k_can+1) = -cf_a

        m_a(2*k_can+1) = -cf_b
        m_c(2*k_can) = -cf_b

        cf_c = PHIdir_TOC * weight * T_beam(k_can+1) * (1-tau_beam(k_can)) * leaf_refl_trans * (1-cf_b)
        cf_d = PHIdir_TOC * weight * T_beam(k_can+1)* (1-tau_beam(k_can)) * leaf_refl_trans * (1-cf_b)
        m_d(2*k_can) = cf_d
        m_d(2*k_can+1) = cf_c
      end do

      ! forward sweep
      c_prime(1) = m_c(1)
      d_prime(1) = m_d(1)
      do k_can=2,2*layers-1
          c_prime(k_can) = m_c(k_can)/(1-m_a(k_can)*c_prime(k_can-1))
          d_prime(k_can) = (m_d(k_can) - m_a(k_can)*d_prime(k_can-1)) / (1-m_a(k_can)*c_prime(k_can-1))
      end do

      ! backward sweep (fluxes)
      phidifcan(layers,ib) = PHIdif_TOC * weight
      do k_can=layers,2,-1
          phiucan(k_can,ib) = d_prime(2*k_can-1) - c_prime(2*k_can-1) * phidifcan(k_can,ib)
          phidifcan(k_can-1,ib) = d_prime(2*k_can-2) - c_prime(2*k_can-2) * phiucan(k_can,ib)
      end do
      phiucan(1,ib) = d_prime(1) - c_prime(1) * phidifcan(1,ib)

      do k_can=1,layers-1
        ! mean sunlit fraction within layer
        fracSL_m = (fracSL(k_can+1) - fracSL(k_can)) / (0.5 / cos_sza * clump * LA(k_can))

        Hshad(k_can, ib) = (phiucan(k_can, ib) + phidifcan(k_can+1,ib)) * (1-tau_dif_can(k_can)) * (1-sigma) / LA(k_can)
        Hsun  (k_can,:,ib) =  Hshad (k_can,ib) + angle_g * (PHIdir_TOC*weight * (1.0-sigma) * (1 - tau_beam(k_can)) * T_beam(k_can+1)/ (LA(k_can) * fracSL_m)) / sum(angle_g * weight_g)
      end do

      ! direct beam
      phidircan(:,ib) = PHIdir_TOC*weight * T_beam(:)

      ! soil reflecntace
      isoil(ib) = (phidircan(1,ib) + phidifcan(1,ib)) * (1-alb)

      ! second forward sweep, no without direct
      c_prime(1) = -alb
      d_prime(1) = 0
      do k_can=2,2*layers-1
          c_prime(k_can) = m_c(k_can)/(1-m_a(k_can)*c_prime(k_can-1))
          d_prime(k_can) = (0. - m_a(k_can)*d_prime(k_can-1)) / (1-m_a(k_can)*c_prime(k_can-1))
      end do

      ! the second (no direct) backward sweep
      phidifcan_nd(layers) = PHIdif_TOC * weight
      do k_can=layers,2,-1
          phiucan_nd(k_can) = d_prime(2*k_can-1) - c_prime(2*k_can-1) * phidifcan_nd(k_can)
          phidifcan_nd(k_can-1) = d_prime(2*k_can-2) - c_prime(2*k_can-2) * phiucan_nd(k_can)
      end do
      phiucan_nd(1) = d_prime(1) - c_prime(1) * phidifcan_nd(1)


      !calculate effective albedos for direct, diffuse and total SW radiation assuming PHI and sw albedos are identical
      effalb_dir(ib) = (phiucan(layers,ib) -  phiucan_nd(layers)) / (weight * PHIdir_TOC)
      effalb_dif(ib) = phiucan_nd(layers) / (weight * PHIdif_TOC)
      effalb_phi(ib)  = phiucan(layers,ib) / (weight * (PHIdir_toc + PHIdif_TOC))

    end do
  else ! small sinbeta, night
    Hshad      = 0.0
    Hsun       = 0.0
    fracSL     = 0.0
    phidircan  = 0.0
    phidifcan  = 0.0
    phiucan    = 0.0
    effalb_dir = 0.0
    effalb_dif = 0.0
    effalb_phi = 0.0
 endif ! sinbeta
end subroutine

! Canopy radiative transfer model of Goudriaan & van Laar (1994) https://doi.org/10.1007/978-94-011-0750-1_6
subroutine canopyrad_sw_GvL94(layers,LAI,LAI_can,PHIdir_TOC,PHIdif_TOC,clump, & ! in
                        Hshad,Hsun, fracSL,                                                      & ! out needed for vegetation
                        effalb_dir,effalb_dif,effalb_phi,phidircan,phidifcan,phiucan,isoil) ! out needed for radiation if lcanopyeb
  use modraddata , only : zenith
  use modglobal  , only : xtime,rtimee,xday,xlat,xlon
  implicit none

  integer,intent(in) :: layers
  real, intent(in)   :: LAI                     ! total Leaf Area Index of the whole column
  real, intent(in),dimension(layers) :: LAI_can ! array with LAI above the evaluated level. Array goes from canopy bottom to top.
  real, intent(in)   :: PHIdir_TOC                  ! Direct  radiation at vegetation top, always >0
  real, intent(in)   :: PHIdif_TOC                  ! Diffuse radition at vegetation top, always >0
  real, intent(in)   :: clump                     ! leaf clumping index

  real, intent(out),dimension(layers-1,n_bands_sw)              :: Hshad  ! Intercepted radiation by shaded leaves between the layer and the one below
  real, intent(out),dimension(layers-1,nangle_gauss,n_bands_sw) :: Hsun   ! Intercepted radiation by sunlit leaves between the layer and the one below per layer and per leaf orientation
  real, intent(out),dimension(layers)                        :: fracSL ! fraction of sunlit leaves per layer
  real, intent(out),dimension(n_bands_sw) ::  effalb_dir                           ! effective albedo of canopy top for direct radiation
  real, intent(out),dimension(n_bands_sw) ::  effalb_dif                           ! effective albedo of canopy top for diffuse radiation
  real, intent(out),dimension(n_bands_sw) ::  effalb_phi                            ! effective albedo of canopy top for total radiation
  real, intent(out),dimension(layers,n_bands_sw) :: phidircan          ! downwards direct comp. of radiation inside canopy
  real, intent(out),dimension(layers,n_bands_sw) :: phidifcan          ! downwards diffuse comp. of radiation inside canopy
  real, intent(out),dimension(layers,n_bands_sw) :: phiucan            ! upwards (diffuse) comp. of radition inside canopy
  real, intent(out),dimension(n_bands_sw) :: isoil                                 ! absorbed radiation by soil
  !real, intent(out) :: irefl                                 ! reflected radiaiton by canopy

  real :: kdrbl,kdfbl,kdf,kdr,iLAI
  real :: ref                    ! global PHI albedo of leave layer (smaller than sigma since some radiation is absorbed after the scattering by leaves as diffuse)
  real :: ref_dir                ! direct PHI albedo of leave layer
  real :: PHIdirprof             ! profile of direct PHI inside canopy
  real :: PHIdir_TOC_b           ! Direct  radiation at vegetation top, always >0, band-weighted
  real :: PHIdif_TOC_b           ! Diffuse radition at vegetation top, always >0, band-weighted

  real    :: phidirprof_b,HdfT,HdrT,dirH ! see Pedruzo-Bagazgoitia et al. (2017)
  integer :: k,ib
  real    :: sinbeta,minsinbeta = 1.e-10
  !needed for approach by goudiraan and van laar  1996
  real :: sigma,weight,rho_s,eta_dir,eta_dif,reta1,reta2,corrv1,corrv2,denom1,denom2
  real :: exp_dir_dif,denom_dir,exp_dif_dif,denom_dif
  real :: phid_10_dir,phid_20_dir,phiu_10_dir,phiu_20_dir
  real :: phid_10_dif,phid_20_dif,phiu_10_dif,phiu_20_dif
  real :: phidirup,phidirdown,phidifup,phidifdown
  real :: phidirup_b,phidirdown_b,phidifup_b,phidifdown_b
  real :: phidifnew,phiupnew
  real :: irefl_dif,irefl_dir
  real :: irefl

  Hshad       = 0.0
  Hsun        = 0.0
  sinbeta  = max(zenith(xtime*3600 + rtimee,xday,xlat,xlon), minsinbeta)
  if (sinbeta>0.035) then ! day, same threshold as radpar
    do ib = 1, n_bands_sw
      sigma = canrad_bands_sw(ib)%sigma
      weight = canrad_bands_sw(ib)%weight

      rho_s = canrad_bands_sw(ib)%alb_soil
      kdfbl    = 0.8                                     ! Diffuse radiation extinction coefficient for black leaves
      kdrbl    = 0.5 / sinbeta                           ! Direct radiation extinction coefficient for black leaves
      kdf      = clump * kdfbl * sqrt(1.0-sigma)
      kdr      = clump * kdrbl * sqrt(1.0-sigma)

      ref      = rho_c_dif(sigma) ! canopy reflectance fitted against Goudriaan (1977) iterative method
      ref_dir  = rho_c_dir(sinbeta, ref) ! canopy reflectance fitted against Goudriaan (1977) iterative method
      PHIdir_TOC_b = PHIdir_toc * weight
      PHIdif_TOC_b = PHIdif_toc * weight

      !follow goudiraan and van laar  1996

      !DIR:
      eta_dir = (ref_dir-rho_s) /(rho_s-(1./ref_dir))
      exp_dir_dif = exp(-kdf*LAI) * exp(-kdr*LAI)
      denom_dir = (1+eta_dir*exp_dir_dif)
      phid_10_dir = PHIdir_TOC_b/denom_dir
      phid_20_dir = PHIdir_TOC_b*eta_dir*exp_dir_dif / denom_dir
      phiu_10_dir = PHIdir_TOC_b*ref_dir/denom_dir
      phiu_20_dir = PHIdir_TOC_b*eta_dir*exp_dir_dif / (ref_dir * denom_dir)

      !DIF
      eta_dif = (ref-rho_s) /(rho_s-(1./ref))
      exp_dif_dif = exp(-2*kdf*LAI)
      denom_dif = (1+eta_dif*exp_dif_dif)
      phid_10_dif = PHIdif_TOC_b/denom_dif
      phid_20_dif = PHIdif_TOC_b*eta_dif*exp_dif_dif / denom_dif
      phiu_10_dif = PHIdif_TOC_b*ref/denom_dif
      phiu_20_dif = PHIdif_TOC_b*eta_dif*exp_dif_dif / (ref * denom_dif)

      !appendix 4 in Goud1996
      reta1 = (ref-rho_s)/(rho_s*ref-1)
      reta2 = (ref_dir-rho_s)/(rho_s*ref_dir-1)
      corrv1 = reta1*exp(-2*kdf*LAI)*ref
      corrv2 = reta2*exp(-kdr*LAI)*exp(-kdf*LAI)*ref_dir
      denom1 = 1+corrv1
      denom2 = 1+corrv2

      do k = 1, layers ! loop over the different LAI locations
        iLAI        = LAI_can(k)   ! Integrated LAI between here and canopy top
        if (ib==1) fracSL(k)   = clump * exp(-kdrbl * clump * iLAI)      ! Fraction of sun-lit leaves


        phidirdown    = phid_10_dir * exp(-kdr * iLAI) + phid_20_dir * exp(kdf * iLAI)    ! PHI down due to PHIdirTOC
        phidirup      = phiu_10_dir * exp(-kdr * iLAI) + phiu_20_dir * exp(kdf * iLAI) ! PHI up due to PHIdirTOC
        phidifdown    = phid_10_dif * exp(-kdf * iLAI) + phid_20_dif * exp(kdf * iLAI)    ! PHI down due to PHIdifTOC
        phidifup      = phiu_10_dif * exp(-kdf * iLAI) + phiu_20_dif * exp(kdf * iLAI)    ! PHI up due to PHIdifTOC
        phidirprof    = PHIdir_TOC_b * exp(-kdrbl * clump * iLAI)              ! Purely direct PHI (can only be downward) reaching leaves
        phidifnew     = phidifdown + phidirdown-phidirprof  ! PHI dif profile
        phiupnew      = phidirup + phidifup                        ! PHI up (dif) profile

        if (ldiscr) then
        !using the discrete definition of the differentiation
          if (k>1) then
          !shift one k level beacuse we move from half level to the full level below
          ! following goudriaan and taking discrete derivative
            HdfT =( (phidifdown - phidifdown_b) / (-(LAI_can(k)-LAI_can(k-1))) + &
                    (phidifup   - phidifup_b)   / (  LAI_can(k)-LAI_can(k-1)))
            HdrT =( (phidirdown - phidirdown_b) / (-(LAI_can(k)-LAI_can(k-1))) + &
                    (phidirup   - phidirup_b)   / (  LAI_can(k)-LAI_can(k-1)))
            dirH = (phidirprof   - phidirprof_b) / (-(LAI_can(k)-LAI_can(k-1)))

            Hshad (k-1,ib)   =  HdfT + HdrT - dirH * (1-sigma)
            Hsun  (k-1,:,ib) =  Hshad (k-1,ib) + (angle_g * (1.0-sigma) * kdrbl * PHIdir_TOC_b / sum(angle_g * weight_g))
          endif
          !save values needed for canopy level above
          phidirprof_b = phidirprof
          phidirdown_b = phidirdown
          phidifdown_b = phidifdown
          phidirup_b   = phidirup
          phidifup_b   = phidifup
        else
        ! using analytical differentiation
          HdfT = (1-ref)     * PHIdif_TOC_b * kdf *(exp(-kdf * iLAI) + exp(kdf * iLAI) * corrv1/ref)    /denom1 ! VISDF, p217 goudriaan 1996
          HdrT = (1-ref_dir) * PHIdir_TOC_b * kdr *(exp(-kdr * iLAI) + exp(kdr * iLAI) * corrv2/ref_dir)/denom2 ! VIST, p217 goudriaan 1996
          dirH   = kdrbl * PHIdirprof !eq 6.31 applied for direct
          Hshad(k,ib)  = HdfT + HdrT - dirH * (1-sigma) ! eq 6.33
          Hsun(k,:,ib) = Hshad(k,ib) + (angle_g * (1.0-sigma) * kdrbl * PHIdir_TOC_b / sum(angle_g * weight_g))
        endif
        phidircan(k,ib) = phidirprof
        phidifcan(k,ib) = phidifnew
        phiucan(k,ib)   = phiupnew
      end do

      ! PHI reflected by canopy, see appendix4 in Goud 1996
       irefl     = PHIdir_TOC_b*(ref_dir+corrv2/ref_dir)/denom2 + PHIdif_TOC_b*(ref+corrv1/ref)/denom1! Goudriaan 1996, p218
       irefl_dif = PHIdif_TOC_b*(ref+corrv1/ref)/denom1
       irefl_dir = PHIdir_TOC_b*(ref_dir+corrv2/ref_dir)/denom2
       effalb_phi(ib) = irefl/(PHIdir_TOC_b+PHIdif_TOC_b)
       effalb_dir(ib) = irefl_dir/PHIdir_TOC_b
       effalb_dif(ib) = irefl_dif/PHIdif_TOC_b

       isoil(ib) = (phidircan(1,ib) + phidifcan(1,ib)) * (1-rho_s)
  end do

  else ! small sinbeta, night
    Hshad      = 0.0
    Hsun       = 0.0
    fracSL     = 0.0
    effalb_dir = 0.0
    effalb_dif = 0.0
    effalb_phi = 0.0
    phidircan  = 0.0
    phidifcan  = 0.0
    phiucan    = 0.0
    isoil      = 0.0
  endif ! sinbeta
return
end subroutine ! canopyrad_sw_GvL94

! Norman 1979 radiative transfer model, following tridiagonal matrix solution of Bonan 2019 (https://doi.org/10.1017/9781107339217.015)
subroutine canopyrad_lw_N79(layers,LA,lwd_TOC,T_shad,T_sun, cfSL, tau_dif_can, & ! in
                        tskin_sfc, leaf_eps,                          & ! in
                        lwd_can, lwu_can, LW_in_leaf_shad, LW_in_leaf_sun, only_fluxes)    ! out needed for vegetation
  use modraddata, only : sfc_emis
  use modglobal, only : boltz
  implicit none

  integer,intent(in) :: layers
  real, intent(in),dimension(layers) :: LA    ! local leaf area (m2 leaf/m2 ground)
  real, intent(in),dimension(layers)   :: tau_dif_can
  real, intent(in),dimension(layers) :: T_shad, T_sun ! Shaded and sunlit leaf temperatures
  real, intent(in),dimension(layers) :: cfSL    ! fraction of sunlit leaves per layer
  real, intent(in)   :: lwd_TOC                  ! Downwelling longwave irradiance at vegetation top
  real, intent(in)   :: tskin_sfc, leaf_eps
  logical, intent(in)   :: only_fluxes

  real, intent(out),dimension(layers)    :: LW_in_leaf_shad, LW_in_leaf_sun  ! LW radiation going into leaves from both side [W/m2leaf]
  real, intent(out),dimension(layers+1)    :: lwd_can, lwu_can  ! LW fluxes

  real, dimension(2*(layers+1)) :: m_a, m_c, m_d ! matrix coefficients (m_b==1)
  real, dimension(2*(layers+1)) :: c_prime, d_prime
  real, dimension(layers) :: sb_layer_mean
  real                    :: cf_a, cf_b, cf_c, cf_d ! temporary coefficients
  real                    ::  leaf_refl_trans, H, lay_emis

  integer :: k_can

  leaf_refl_trans = (1-leaf_eps) / 2.

  m_d(1) = sfc_emis * boltz * tskin_sfc**4 ! bottom boundary condition
  m_c(1) = -(1-sfc_emis)
  do k_can=1, layers
    cf_a = (1-tau_dif_can(k_can))*leaf_refl_trans - ((tau_dif_can(k_can) + (1-tau_dif_can(k_can))*leaf_refl_trans)**2) / ((1-tau_dif_can(k_can))*leaf_refl_trans)
    cf_b = (tau_dif_can(k_can) + (1-tau_dif_can(k_can))*leaf_refl_trans) / ((1-tau_dif_can(k_can))*leaf_refl_trans)
    m_a(2*k_can) = -cf_a
    m_c(2*k_can+1) = -cf_a

    m_a(2*k_can+1) = -cf_b
    m_c(2*k_can) = -cf_b

    sb_layer_mean(k_can) = leaf_eps * boltz * ( (T_sun(k_can)**4) * cfSL(k_can) + (t_shad(k_can)**4) * (1.-cfSL(k_can)) )
    cf_c = sb_layer_mean(k_can) * (1-tau_dif_can(k_can)) * (1-cf_b)
    cf_d = sb_layer_mean(k_can) * (1-tau_dif_can(k_can)) * (1-cf_b)
    m_d(2*k_can) = cf_d
    m_d(2*k_can+1) = cf_c
  end do

  ! forward sweep
  c_prime(1) = m_c(1)
  d_prime(1) = m_d(1)
  do k_can=2,2*layers+1
      c_prime(k_can) = m_c(k_can)/(1-m_a(k_can)*c_prime(k_can-1))
      d_prime(k_can) = (m_d(k_can) - m_a(k_can)*d_prime(k_can-1)) / (1-m_a(k_can)*c_prime(k_can-1))
  end do

  ! backward sweep (fluxes)
  lwd_can(layers+1) = lwd_TOC
  do k_can=layers+1,2,-1
      lwu_can(k_can) = d_prime(2*k_can-1) - c_prime(2*k_can-1) * lwd_can(k_can)
      lwd_can(k_can-1) = d_prime(2*k_can-2) - c_prime(2*k_can-2) * lwu_can(k_can)
  end do
  lwu_can(1) = d_prime(1) - c_prime(1) * lwd_can(1)


  if (.not. only_fluxes) then
    ! compute LW radiation going into leaves
    do k_can = 1, layers
      ! net lw fluxes into layer
      H = (leaf_eps * (lwd_can(k_can+1) + lwu_can(k_can)) - 2 * sb_layer_mean(k_can) ) * (1-tau_dif_can(k_can))

      ! total LW emitted by leaves in layer
      lay_emis = 2 * sb_layer_mean(k_can) * LA(k_can)

      ! net lw fluxes into leaves: net flux + lay_emis
      ! Currently, no split between sunlit and shaded leaves
      LW_in_leaf_shad(k_can) = (H + lay_emis)/LA(k_can)
      LW_in_leaf_sun(k_can) = LW_in_leaf_shad(k_can)
    end do
  end if

end subroutine


subroutine f_Ags(CO2air,qtair,dens,tairk,pair,t_skin,          & ! in
                 phi_tot,Hleaf,                                & ! in
                 lrelaxgc,gc_old_set,kgc,gleaf_old,rk3coef,    & ! in
                 lrelaxci,ci_old_set,kci,ci_old,               & ! in
                 gleaf,Fleaf,ci,                               & ! out
                 fstr,Am,Rdark,alphac,co2abs,CO2comp,Ds,D0,fmin) ! additional out for 1leaf upscaling

      implicit none
      real, intent(in) ::  CO2air   ! CO2 air concentration [ppb]
      real, intent(in) ::  qtair    ! air specific humidity [kg/kg]
      real, intent(in) ::  pair     ! air pressure          [Pa]
      real, intent(in) ::  dens     ! air density           [kg/m3]
      real, intent(in) ::  tairk    ! air temperature (K)
      real, intent(in) ::  t_skin   ! surface or leaf skin temperature (K)
      real, intent(in) ::  phi_tot  ! Total soil water content [-]
      real, intent(in) ::  Hleaf    ! absorbed PAR [W/m2]
      logical, intent(in) ::  lrelaxgc ! if true, lag in plant stomatal response
      logical, intent(in) ::  gc_old_set !false at first timestep, true at any  other
      real, intent(in) ::  kgc      ! Standard stomatal response rate
      real, intent(in) ::  gleaf_old ! gleaf at previous timestep
      real, intent(in) ::  rk3coef    ! runge kutta step
      logical, intent(in) ::  lrelaxci ! if true, lag in plant internal carbon response
      logical, intent(in) ::  ci_old_set !false at first timestep, true at any  other
      real, intent(in) ::  kci           ! Standard internal co2 concentration response rate
      real, intent(in) ::  ci_old        ! ci at previous timestep


      real, intent (out) :: gleaf   ! leaf stomatal conductance for carbon(m/s)
      real, intent (out) :: Fleaf
      real, intent (out) :: ci

      real, intent (out) :: fstr
      real, intent (out) :: Am
      real, intent (out) :: Rdark
      real, intent (out) :: alphac
      real, intent (out) :: co2abs
      real, intent (out) :: CO2comp
      real, intent (out) :: Ds
      real, intent (out) :: D0
      real, intent (out) :: fmin

      real     :: CO2ags, gm, fmin0, esatsurf, cfrac !Variables for AGS
      real     :: Ammax, betaw!Variables for AGS
      real     :: Agl,gleaf_inf
      ! variables below ma not be necessaruly defined when implemented properly
      real     :: e! needed here as well
      !this line is done in do_lsm, but we put it here because it will change per vert level
      e    = qtair * pair / 0.622

      CO2ags  = CO2air/1000.0  !From ppb (usual DALES standard) to ppm

      ! Calculate surface resistances using the plant physiological (A-gs) model
      ! Calculate the CO2 compensation concentration
      CO2comp = dens * CO2comp298 * Q10CO2 ** (0.1 * ( tairk - 298.0 ) )
      !CO2comp = dens * CO2comp298 * Q10CO2 ** (0.1 * ( t_skin - 298.0 ) ) !!proposal pending

      ! Calculate the mesophyll conductance
      gm       = gm298 * Q10gm ** (0.1 * ( tairk - 298.0) ) / ( (1. + exp(0.3 * ( T1gm - tairk ))) * (1. + exp(0.3 * (tairk- T2gm))))
      gm       = gm / 1000   ! conversion from mm s-1 to m s-1

      ! calculate CO2 concentration inside the leaf (ci)
      fmin0    = gmin/nuco2q - (1.0/9.0) * gm
      fmin     = (-fmin0 + ( fmin0 ** 2.0 + 4 * gmin/nuco2q * gm ) ** (0.5)) / (2. * gm)

      esatsurf = 0.611e3 * exp(17.2694 * (t_skin - 273.16) / (t_skin - 35.86))
      Ds       = (esatsurf - e) / 1000.0 ! In kPa
      D0       = (f0 - fmin) / ad

      cfrac    = f0 * (1.0 - Ds/D0) + fmin * (Ds/D0)
      co2abs   = CO2ags * (MW_CO2/MW_Air) * dens

      ci_inf   = cfrac * (co2abs - CO2comp) + CO2comp
      if (lrelaxci .and. ci_old_set) then
        ci            = ci_old + min(kci*rk3coef, 1.0) * (ci_inf - ci_old)
      else
        ci            = ci_inf
      endif

      ! Calculate maximal gross primary production in high light conditions (Ag)
      Ammax    = Ammax298 * Q10am ** ( 0.1 * ( tairk - 298.0) ) / ( (1.0 + exp(0.3 * ( T1Am - tairk ))) * (1. + exp(0.3 * ( tairk - T2Am))) )

      ! Calculate the effect of soil moisture stress on gross assimilation rate
      betaw    = max(1.0e-3,min(1.0,(phi_tot-phiwp)/(phifc-phiwp)))

      ! Calculate stress function
      fstr     = betaw

      ! Calculate gross assimilation rate (Am)
      Am       = Ammax * (1 - exp( -(gm * (ci - CO2comp) / Ammax) ) )

      Rdark    = (1.0/9) * Am

      ! Calculate the light use efficiency
      alphac   = alpha0 * (co2abs  - CO2comp) / (co2abs + 2 * CO2comp)

      Agl    = fstr * (Am + Rdark) * (1 - exp(-alphac*Hleaf/(Am + Rdark)))
      !gleaf_inf  = gmin/nuco2q +  Agl/(co2abs-ci)
      gleaf_inf  = max(gmin/nuco2q, Agl/(co2abs-ci)) !Xabi test

      Fleaf  = -(Agl - Rdark) ! we flip sign here to be consistent with An in old DALES

      if (lrelaxgc .and. gc_old_set) then
          gleaf       = gleaf_old + min(kgc*rk3coef, 1.0) * (gleaf_inf - gleaf_old)
        else
          gleaf = gleaf_inf
      endif
   return
end subroutine !Ags

end module modsurfdata
