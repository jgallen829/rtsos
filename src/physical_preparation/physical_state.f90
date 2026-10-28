! ===============================================================
!  RTSOS – Radiative Transfer model based on Successive Orders of Scattering
!  Copyright (C) 2025  Pengwang Zhai
!
!  Licensed under the Creative Commons Attribution–NonCommercial 4.0
!  International License (CC BY-NC 4.0).
!  You may use, modify, and share this code for research and
!  educational purposes with proper attribution.
!  Commercial use requires written permission from the author.
!
!  Full license: https://creativecommons.org/licenses/by-nc/4.0/
!  Contact: Pengwang Zhai  |  [pwzhai@gmail.com]
! ===============================================================

MODULE PHYSICAL_STATE

LOGICAL :: ILS_FLAG,SNG_SCAT_EXACT,&
		   DIFFUSE_TRANSMITTANCE, ATMOS_ZERO,SUNGLINT_INPUT

!IF ATMOS_ZERO=.TRUE., NO ATMOSPHERE WILL BE PRESENT, i.e., TAUR=0, TAU_AEROSOL=0.

INTEGER :: NDET
DOUBLE PRECISION :: WV_PACE_REF
INTEGER,PARAMETER :: NUMMIEUSE=2
DOUBLE PRECISION,DIMENSION(NUMMIEUSE) :: Height_Variance_Particle,Height_Particle
DOUBLE PRECISION,DIMENSION(0:NUMMIEUSE) :: CEXT_REF,TAU_REF

INTEGER,PARAMETER :: NFINEMODE=10
INTEGER,PARAMETER ::NWV_PF_LUT=2250
DOUBLE PRECISION,DIMENSION(NWV_PF_LUT) :: WV_PF_LUT,AW_PF_LUT,BW_PF_LUT

INTEGER,PARAMETER ::NWV_AW_LEE_LUT=144
DOUBLE PRECISION,DIMENSION(NWV_AW_LEE_LUT) :: WV_AW_LEE_LUT,AW_LEE_LUT

INTEGER,PARAMETER ::NWV_AW_MCF_LUT=131
DOUBLE PRECISION,DIMENSION(NWV_AW_MCF_LUT) :: WV_AW_MCF_LUT,AW_MCF_LUT

INTEGER,PARAMETER ::NWV_AW_IOCCG_LUT=211
DOUBLE PRECISION,DIMENSION(NWV_AW_IOCCG_LUT) :: WV_AW_IOCCG_LUT,AW_IOCCG_LUT

INTEGER, PARAMETER :: N_ALT_INIT=104
DOUBLE PRECISION,DIMENSION(N_ALT_INIT) :: ALT_LYRA_INIT= &
	(/ 100.0D0, 95.00D0, 90.0D0,85.00D0,  80.0D0,75.0D0, &
        70.00D0, 65.0D0, 60.00D0,55.0D0,  50.0D0, 47.5D0,  45.0D0,42.5D0,&
        40.0D0, 37.50D0, 35.0D0, 32.50D0, 30.0D0, 27.50D0, 25.0D0,22.50D0,     &
		20.0D0, 19.75D0, 19.5D0, 19.25D0, 19.0D0, 18.75d0, 18.5d0,18.25D0,     &
		18.0D0, 17.75D0, 17.5D0, 17.25D0, 17.0D0, 16.75d0, 16.5d0,16.25D0,     &
		16.0D0, 15.75D0, 15.5D0, 15.25D0, 15.0D0, 14.75d0, 14.5d0,14.25D0,     &
		14.0D0, 13.75D0, 13.5D0, 13.25D0, 13.0D0, 12.75d0, 12.5d0,12.25D0,     &
		12.0D0, 11.75D0, 11.5D0, 11.25D0, 11.0D0, 10.75d0, 10.5d0,10.25D0,     &
		10.0D0,  9.75D0,  9.5D0,  9.25D0,  9.0D0,  8.75d0,  8.5d0, 8.25D0,     &
		 8.0D0,  7.75D0,  7.5D0,  7.25D0,  7.0D0,  6.75d0,  6.5d0, 6.25D0,     &
		 6.0D0,  5.75D0,  5.5D0,  5.25D0,  5.0D0,  4.75d0,  4.5d0, 4.25D0,     &
		 4.0D0,  3.75D0,  3.5D0,  3.25D0,  3.0D0,  2.75d0,  2.5d0, 2.25D0,     &
		 2.0D0,  1.75D0,  1.5D0,  1.25D0,  1.0D0,  0.75D0,  0.5D0, 0.25D0, 0.0D0, 0.0D0/)
DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) :: ALT_LYRA

DOUBLE PRECISION :: OZONE_COLUMN ! OZONE IN THE WHOLE COLUMN IN DOBSON UNIT
DOUBLE PRECISION :: H2O_COLUMN   ! WATER VAPOR IN THE WHOLE COLUMN IN CENTIMETERS.
DOUBLE PRECISION :: NO2_COLUMN ! NO2 IN THE WHOLE COLUMN IN DOBSON UNIT
DOUBLE PRECISION :: PRESSURE_SURFACE ! SURFACE PRESSURE IN mb.

LOGICAL :: CHLA_HOMOGENEITY
DOUBLE PRECISION :: CHLA_C0,CHLA_Slope,CHLA_Cmax, CHLA_Etamax,CHLA_DeltaEta
DOUBLE PRECISION :: PHYTOPLANKTON_INDEX_REFRACTION,PHYTOPLANKTON_SPECTRAL_SLOPE, &
          SEDIMENT_INDEX_REFRACTION,SEDIMENT_SPECTRAL_SLOPE,        &
          SEDIMENT_CONCENTRATION

DOUBLE PRECISION :: adg440_Case3,bbp660_Case3,Bp660_backscatteringfraction_Case3, &
          Sdg_Case3,Sbp_Case3,S_Bp_Case3

DOUBLE PRECISION ::WNDSPD,THETA0
!DOUBLE PRECISION,PARAMETER :: DEPOL_A=0.0284D0
DOUBLE PRECISION,PARAMETER :: DEPOL_O=0.087D0
DOUBLE PRECISION,PARAMETER :: Xi_kok=25.6d0,Alpha_kok=4.0d0, T0_kok=0.25d0,Dpol90_Kok=0.67d0

!PURE SEA WATER Depolarization ratio value from Zhang X, OE, 2009,
    ! also see Jonasz and Fournier, 2007
DOUBLE PRECISION :: nmbre_local,nmbim_local
LOGICAL :: sngmode,OCEAN_RAMAN_FLAG,OCEAN_FCHLA_FLAG,OCEAN_FCDOM_FLAG,&
           HYSPECTRAL_FLAG,OCEAN_PHMX_ONE,MONOCHROMATIC_FLAG

LOGICAL :: ABSORPTION_TAU_REDUCE_FLAG
DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE ::TAU_REDUCE_ATMOS_PRESERVE,TAU_REDUCE_OCEAN_PRESERVE
DOUBLE PRECISION :: ABSORPTION_TAU_REPLACE

!LOGICAL :: WATER_DET_FREE_LOCATION
!INTEGER,PARAMETER :: N_DET_INWATER_REPORT=21
!DOUBLE PRECISION :: DET_INWATER_REPORT_DEPTH_INCREMENT
!DOUBLE PRECISION,DIMENSION(N_DET_INWATER_REPORT):: DET_INWATER_REPORT_DEPTH ! LARGEST WATER DETECT DEPTH REPORTED

INTEGER :: OCEAN_CASE_SELECT
                              !==-203 Ross Li Surface Reflectance
                              !==-202 Snow Reflectance
                              !==-201 mRPV land Reflectance
                              !==-200 Land Lambertian Reflectance
                              !==0 no ocean water body
                              !==1 CASE 1 WATER IOPS
                              !==2 CASE 2 WATER IOPS
							  !==3 Bio-2 model in polarimeter fitting
INTEGER :: AP_SELECT       !=1 PHYTOPLANKTON ABSORPTION COEFFICIENTS OF BRICAUD
!=2 PHYTOPLANKTON ABSORPTION COEFFICIENT OF
!SF_Coff options from Ciotti, Lewis, and Cullen, LO, 47(2) 404-417, 2002
!Assessment of the relationships between dominant cell size in
!    natural phytoplankton communities and the spectral shape of
!    the absorption coefficient
INTEGER :: WV_SEG_FLAG !wv_seg_flag, 0: all; 1: seg1+3only; 2: seg2only; 3: seg1+2+3; 4: seg4only
DOUBLE PRECISION :: AirSensor_Height
CHARACTER*360 :: aux_dir='00000', land_ref_spectra_filename='00000'

INTEGER :: nwv_land_ref
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: wavelength_land_ref,fiso_input,  &
		   fvol_input,fgeo_input,Bpol_input,nmbre_input,nmbim_input

!CONTAINS
!SUBROUTINE aux_dir_readin
!IF(aux_dir=='00000')THEN
!  OPEN(UNIT=1,FILE='./auxiliary_directory',STATUS='OLD',ACTION='READ')
!  READ(1,'(A)')aux_dir
!  CLOSE(1)
!ELSE
!  RETURN
!ENDIF
!ENDSUBROUTINE

ENDMODULE PHYSICAL_STATE

