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

PROGRAM MONOCHROMATIC
USE ATMOS_CONFIGURATION_DIRECTORY
USE PHYSICAL_STATE
USE PHYSICAL_PARTICLES
USE RAMANDATA
USE FLUORESCENCEDATA
USE SAVE_OCEAN_MOD
USE IRRADSUNL
USE RTUTILITY, ONLY : PI,FACTOR,NUMMIEANGMAX,PSEUDO_SPHERICAL_SHELL,&
                      SPHERICAL_SHELL_SINGLESCATTERING_CORRECTION
USE Atmosphere_Profile
USE MONOCHROMATIC_RESULT
USE SnowBRDF
IMPLICIT NONE

! Physical records retain the public parser order. The final records give the
! explicit solar source, wavelength count, and that many wavelengths in nm.
INTEGER :: NTLYER,NTLYERA,NTLYERO,NREC,NCOLINPUT,NQUADAINPUT,NQUADOINPUT
INTEGER :: MAXMORDINPUT_SAVE,NUMMIERT,MAXLORDINPUT,MAXMORDINPUT,DELTAM_INPUT
INTEGER :: NTHETAOUT,NPHIOUT,IAEROSOL,INTTMP,nopt
INTEGER,ALLOCATABLE :: ARSL_TYPE(:),NUMMIEANGINPUT(:)
DOUBLE PRECISION :: LBDO_LD,FLAM,CHLa,pBRDFa_input,pBRDFk_input,pBRDFb_input,pBRDFe_input
DOUBLE PRECISION :: MU_IN,WAVELENGTH_MICRON,TEMPERTURE,SALINITY,RH_simu,AeroFMF
DOUBLE PRECISION :: requested_nm,esun_value,scale_value,prep_wave(3)
DOUBLE PRECISION,ALLOCATABLE :: TAUR(:),DEPOL_A(:),TAU_TG(:),TAU_TG_H2O(:),TAU_TG_CO2(:)
DOUBLE PRECISION,ALLOCATABLE :: TAU_TG_CH4(:),TAU_TG_O2(:),TAU_TG_O3(:),TAU_TG_NO2(:),PNDLY(:)
DOUBLE PRECISION,ALLOCATABLE :: TAU_ARSL_Ext_TMP(:),TAU_ARSL_Scat_TMP(:)
DOUBLE PRECISION,ALLOCATABLE :: RECDATASTREAM(:,:),PHMXDATASTREAM(:,:,:),PHIOUT(:),MUOUT(:)
DOUBLE PRECISION,ALLOCATABLE :: DIRAD(:,:),DSTOKES(:,:,:,:)
DOUBLE PRECISION,ALLOCATABLE :: gc(:,:),gh(:,:),gco(:,:),gch(:,:),go(:,:),goz(:,:),gn(:,:)
LOGICAL :: WLR_FLAG_INPUT,WCFLG_INPUT,GAS_ABS_FLAG,file_e
CHARACTER(LEN=360) :: Land_HEADER,CFILE_AP,INFILE,OUTFILE,CFILE_AEROSOLS(NUMMIEUSE)
CHARACTER(LEN=20) :: solar_policy
INTEGER :: NARGS,IARGC,nwave,iwave
DOUBLE PRECISION,ALLOCATABLE :: wavelengths(:),esun(:),irrad_all(:,:,:),stokes_all(:,:,:,:,:)
DOUBLE PRECISION,ALLOCATABLE :: taur_all(:,:),gas_all(:,:),ext_all(:,:),scat_all(:,:)
DOUBLE PRECISION :: start_time,prepared_time,finish_time
! Read inputs from the command line
NARGS = IARGC()
IF(NARGS.LT.1) THEN
    WRITE(*,*) ''
    WRITE(*,*) '***ERROR***'
    WRITE(*,*) 'TOO FEW INPUTS'
    WRITE(*,*) '***ERROR***'
    WRITE(*,*) ''
    WRITE(*,*) 'USEAGE:'
    WRITE(*,*) 'rtsos_Monochromatic input.file'
    WRITE(*,*) ''
    STOP
ENDIF

!  First argument is the input file
CALL GETARG(1,INFILE)

! Land does not serialize inactive water fields. Give all inactive state explicit values.
OCEAN_RAMAN_FLAG=.false.; OCEAN_FCHLA_FLAG=.false.; OCEAN_FCDOM_FLAG=.false.
NPQ_FLAG=.false.; OCEAN_PHMX_ONE=.true.; CHLA_HOMOGENEITY=.true.
SUNGLINT_INPUT=.true.; AP_SELECT=1; CHLa=0.0D0; WATER_DEPTH_MAX=200.0D0
LBDO_LD=0.3D0
OPEN(unit=1,file=INFILE,status='old')

READ(1,*)WNDSPD
READ(1,*)THETA0
READ(1,*)WV_PACE_REF
READ(1,*)TAU_REF(1)
READ(1,*)Height_Particle(1)
READ(1,*)Height_Variance_Particle(1)
READ(1,*)TAU_REF(2)
READ(1,*)Height_Particle(2)
READ(1,*)Height_Variance_Particle(2)
READ(1,*)IAEROSOL
READ(1,*)AeroFMF              ! user input value used only when iaerosol==-1
READ(1,*)RH_simu

READ(1,*)OCEAN_CASE_SELECT

IF(OCEAN_CASE_SELECT>=0)THEN
	READ(1,*)LBDO_LD
	READ(1,*)WATER_DEPTH_MAX
	READ(1,*)CHLa
	READ(1,*)PHYTOPLANKTON_INDEX_REFRACTION
	READ(1,*)PHYTOPLANKTON_SPECTRAL_SLOPE
	READ(1,*)SEDIMENT_INDEX_REFRACTION
	READ(1,*)SEDIMENT_SPECTRAL_SLOPE
	READ(1,*)SEDIMENT_CONCENTRATION

	READ(1,*)adg440_Case3
	READ(1,*)bbp660_Case3
	READ(1,*)Bp660_backscatteringfraction_Case3
	READ(1,*)Sdg_Case3
	READ(1,*)Sbp_Case3
	READ(1,*)S_Bp_Case3
ENDIF
READ(1,*)NCOLINPUT
READ(1,*)NQUADAINPUT
READ(1,*)NQUADOINPUT
READ(1,*)MAXMORDINPUT_SAVE
READ(1,*)NTHETAOUT
READ(1,*)NPHIOUT

IF(OCEAN_CASE_SELECT>=0)THEN
	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
		CHLA_HOMOGENEITY=.FALSE.
	ELSE
		CHLA_HOMOGENEITY=.TRUE.
	ENDIF

	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
	   OCEAN_RAMAN_FLAG=.FALSE.
	ELSE
	   OCEAN_RAMAN_FLAG=.TRUE.
	ENDIF

	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
	   OCEAN_FCHLA_FLAG=.FALSE.
	ELSE
	   OCEAN_FCHLA_FLAG=.TRUE.
	ENDIF

	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
	   OCEAN_FCDOM_FLAG=.FALSE.
	ELSE
	   OCEAN_FCDOM_FLAG=.TRUE.
	ENDIF

	READ(1,*)AP_SELECT
ENDIF

READ(1,*)INTTMP
IF(INTTMP==0)THEN
   MONOCHROMATIC_FLAG=.FALSE.
ELSE
   MONOCHROMATIC_FLAG=.TRUE.
ENDIF


READ(1,*)INTTMP
IF(INTTMP==0)THEN
   HYSPECTRAL_FLAG=.FALSE.
ELSE
   HYSPECTRAL_FLAG=.TRUE.
ENDIF

IF(OCEAN_CASE_SELECT>=0)THEN
	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
	   NPQ_FLAG=.FALSE.
	ELSE
	   NPQ_FLAG=.TRUE.
	ENDIF

	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
	   OCEAN_PHMX_ONE=.FALSE. ! if true, only use one phase matrix for ocean waters
	! However, single scattering abeldo change is still kept
	ELSE
	   OCEAN_PHMX_ONE=.TRUE.
	ENDIF
ENDIF

READ(1,*)INTTMP
IF(INTTMP==0)THEN
  ATMOS_ZERO=.FALSE.
ELSE
  ATMOS_ZERO=.TRUE. ! NO ATMOSPHERE WILL BE PRESENT
ENDIF

IF(OCEAN_CASE_SELECT>=0)THEN
	READ(1,*)INTTMP
	IF(INTTMP==0)THEN
		SUNGLINT_INPUT=.true.
	ELSE
		SUNGLINT_INPUT=.false. ! NO Sun glint WILL BE PRESENT
	ENDIF
ENDIF

READ(1,*)INTTMP
IF(INTTMP==0)THEN
   GAS_ABS_FLAG=.FALSE. ! if .false. not gas absorption is included
ELSE
   GAS_ABS_FLAG=.TRUE.  ! if .true. gas absorption is included
ENDIF

READ(1,*)OZONE_COLUMN ! OZONE IN THE WHOLE COLUMN IN DOBSON UNIT
READ(1,*)H2O_COLUMN   ! WATER VAPOR IN THE WHOLE COLUMN IN CENTIMETERS.
READ(1,*)NO2_COLUMN   ! NO2 IN THE WHOLE COLUMN IN DOBSON UNIT
READ(1,*)PRESSURE_SURFACE ! surface pressure in mb
READ(1,*)WV_SEG_FLAG
READ(1,*)AirSensor_Height
READ(1,*)INTTMP
	IF(INTTMP==1)THEN
        SPHERICAL_SHELL_SINGLESCATTERING_CORRECTION=.true.
		PSEUDO_SPHERICAL_SHELL=.TRUE.
	ELSE
        SPHERICAL_SHELL_SINGLESCATTERING_CORRECTION=.false.
		PSEUDO_SPHERICAL_SHELL=.FALSE.
	ENDIF

READ(1,'(A)')aux_dir
READ(1,'(A)')atmos_dir

READ(1,'(A)')CFILE_AP    ! Atmosphere number density profile
READ(1,'(A)')CFILE_AEROSOLS(1) ! AEROSOL PHASE MATRIX FILE
READ(1,'(A)')CFILE_AEROSOLS(2) ! AEROSOL PHASE MATRIX FILE top layer

IF(OCEAN_CASE_SELECT < 0) READ(1,'(A)')land_ref_spectra_filename  !land surface reflectance filename

READ(1,'(A)')OUTFILE
READ(1,'(A)')solar_policy
READ(1,*)nwave
IF(nwave<1) ERROR STOP 'At least one wavelength is required'
ALLOCATE(wavelengths(nwave),esun(nwave))
DO iwave=1,nwave
 READ(1,*)wavelengths(iwave)
ENDDO
requested_nm=wavelengths(1)
CLOSE(1)

IF(index(aux_dir,'#')>1) &
  aux_dir=trim(aux_dir(1:index(aux_dir,'#')-1))
IF(index(atmos_dir,'#')>1) &
  atmos_dir=trim(atmos_dir(1:index(atmos_dir,'#')-1))
IF(index(CFILE_AP,'#')>1) &
   CFILE_AP=trim(CFILE_AP(1:index(CFILE_AP,'#')-1))
DO INTTMP=1,NUMMIEUSE
   IF(index(CFILE_AEROSOLS(INTTMP),'#')>1) &
        CFILE_AEROSOLS(INTTMP)=trim(CFILE_AEROSOLS(INTTMP)(1:index(CFILE_AEROSOLS(INTTMP),'#')-1))
ENDDO
IF(index(OUTFILE,'#')>1) &
   OUTFILE=trim(OUTFILE(1:index(OUTFILE,'#')-1))

IF(OCEAN_CASE_SELECT<0) THEN
   IF(OCEAN_CASE_SELECT==-202) snow_brdf_path=trim(aux_dir)//'SnowBRDF_Data/brdfPaperMaterial/'
   IF(index(land_ref_spectra_filename,'#')>1) &
	  land_ref_spectra_filename=trim(land_ref_spectra_filename(1:index(land_ref_spectra_filename,'#')-1))
   OPEN(unit=1,file=land_ref_spectra_filename,status='old')
   READ(1,*)Land_HEADER
   READ(1,*)nwv_land_ref
   ALLOCATE(wavelength_land_ref(nwv_land_ref),fiso_input(nwv_land_ref),&
			fvol_input(nwv_land_ref),fgeo_input(nwv_land_ref),Bpol_input(nwv_land_ref),&
            nmbre_input(nwv_land_ref),nmbim_input(nwv_land_ref))
   DO INTTMP=1,nwv_land_ref
	  READ(1,*)wavelength_land_ref(INTTMP),fiso_input(INTTMP),&
			   fvol_input(INTTMP),fgeo_input(INTTMP),Bpol_input(INTTMP),&
               nmbre_input(INTTMP),nmbim_input(INTTMP)
   ENDDO
   CLOSE(1)


ENDIF


IF(Height_Particle(1)<Height_Particle(2)) &
   stop 'check, Height_Particle_toplayer<Height_Particle_lowlayer'
IF(WV_SEG_FLAG<0 .OR. WV_SEG_FLAG>4) STOP 'CHECK WV_SEG_FLAG INPUT'

IF(OCEAN_CASE_SELECT>=0)THEN
	write(*,*)'Wind speed=',WNDSPD
	write(*,*)'Chlorophyll a concentration=',CHLa
	write(*,*)'PHYTOPLANKTON_INDEX_REFRACTION=',PHYTOPLANKTON_INDEX_REFRACTION
	write(*,*)'PHYTOPLANKTON_SPECTRAL_SLOPE=',PHYTOPLANKTON_SPECTRAL_SLOPE
	WRITE(*,'(1x,"CHLA_HOMOGENEITY=",L2)')CHLA_HOMOGENEITY
	WRITE(*,'(1x,"OCEAN_RAMAN_FLAG,OCEAN_FCHLA_FLAG,OCEAN_FCDOM_FLAG=",3(L2,2x))') &
			 OCEAN_RAMAN_FLAG,OCEAN_FCHLA_FLAG,OCEAN_FCDOM_FLAG
	WRITE(*,'(1x,"HYSPECTRAL_FLAG,NPQ_FLAG,OCEAN_PHMX_ONE,GAS_ABS_FLAG=",4(L2,2x))')&
			   HYSPECTRAL_FLAG,NPQ_FLAG,OCEAN_PHMX_ONE,GAS_ABS_FLAG

	WRITE(*,'(1x,"adg440,bbp660,Bp660=",3(E12.4,2x))') &
			   adg440_Case3,bbp660_Case3,Bp660_backscatteringfraction_Case3
    WRITE(*,'(1x,"Sdg,Sbp,S_Bp=",3(E12.4,2x))') &
               Sdg_Case3,Sbp_Case3,S_Bp_Case3
ENDIF

WRITE(*,'(1x,"SPHERICAL_SHELL_SINGLESCATTERING_CORRECTION IS:",L2)') &
		   SPHERICAL_SHELL_SINGLESCATTERING_CORRECTION

IF(ATMOS_ZERO)THEN
    TAU_REF=0.0d0
    GAS_ABS_FLAG=.false.
    WRITE(*,*)'ATMOS_ZERO IS TRUE, NO GAS ABSORPTION IS INCLUDED'
ENDIF

IF(.NOT. GAS_ABS_FLAG) MONOCHROMATIC_FLAG=.TRUE.

IF(MONOCHROMATIC_FLAG)THEN
  HYSPECTRAL_FLAG=.FALSE.
  WRITE(*,*) 'HYSPECTRAL_FLAG IS RESET TO .FALSE. DUE TO MONOCHROMATIC_FLAG'
ENDIF

IF(OCEAN_CASE_SELECT==3 .AND. AP_SELECT .NE. 1)THEN
  AP_SELECT=1
  WRITE(*,*) 'AP_SELECT IS RESET TO 1, BECAUSE OCEAN_CASE_SELECT==3'
ENDIF

! check file existance

INTTMP=len_trim(OUTFILE)
IF(OUTFILE(INTTMP-2:INTTMP) .ne. '.h5') OUTFILE=TRIM(OUTFILE)//'.h5'
WRITE(*,*)'OUTPUT FILE IS:',OUTFILE

INQUIRE(FILE=OUTFILE, EXIST=file_e)
IF(file_e) then
!  write(*,*) 'warning, filename exist, a string is appended'
!  OUTFILE='new_'//TRIM(OUTFILE)
  stop 'warning, output filename exist, exiting'
ENDIF

IF(ABS(AirSensor_Height)<1.0E-6 .OR. ABS(AirSensor_Height)>=100.0D0)THEN
  WRITE(*,*)'AirSensor_Height ==0 .or. AirSensor_Height >=100 km'
  WRITE(*,*)'Default AirSensor_Height = 2 km is used'
  AirSensor_Height=2.0D0
ENDIF

NTLYERA=N_ALT_INIT
ALLOCATE(ALT_LYRA(NTLYERA))
ALT_LYRA_INIT(N_ALT_INIT)=AirSensor_Height
ALT_LYRA=ALT_LYRA_INIT
CALL qsort2a(NTLYERA,ALT_LYRA_INIT,ALT_LYRA)
CALL REMOVE_DUP(N_ALT_INIT,ALT_LYRA_INIT,NTLYERA,ALT_LYRA)
DO INTTMP=1,NTLYERA
   ALT_LYRA_INIT(INTTMP)=ALT_LYRA(NTLYERA-INTTMP+1)
ENDDO
DEALLOCATE(ALT_LYRA)

NTLYERA=NTLYERA-1
ALLOCATE(ALT_LYRA(NTLYERA+1))
ALT_LYRA(1:NTLYERA+1)=ALT_LYRA_INIT(1:NTLYERA+1)

!write(*,*)'NTLYERA,ALT_LYRA=',NTLYERA,ALT_LYRA
IF(OCEAN_CASE_SELECT<0)THEN
   NDET=3
ELSE
   NDET=4
ENDIF

ABSORPTION_TAU_REDUCE_FLAG=.true. ! IF TRUE, USE THINNER OCEAN LAYERS BELOW ABSORPTION OPTICAL DEPTH OF 10.
! THIS WOULD EXPEDITE THE SIMULATION VERY MUCH WITHOUT SACRIFY THE ACCURACY AT TOO

FLAM=1.0
TEMPERTURE=20.0D0
SALINITY=37.0D0

ILS_FLAG=.true.
SNG_SCAT_EXACT=.FALSE.
DIFFUSE_TRANSMITTANCE=.false.

MU_IN=-1.0d0*ABS(COS(THETA0/180.D0*PI))


CALL CPU_TIME(start_time)
! Scene preparation is shared across the wavelength set. Optical storage has
! only the current wavelength and physical reference slots, independent of count.
IF(OCEAN_RAMAN_FLAG .OR. OCEAN_FCHLA_FLAG .OR. OCEAN_FCDOM_FLAG .OR. NPQ_FLAG) &
 ERROR STOP 'Requested-wavelength inelastic excitation/emission preparation is not connected'
IF(IAEROSOL < -1) ERROR STOP 'Requested-wavelength external particle matrix loading is not connected'
prep_wave(1)=requested_nm
prep_wave(2)=WV_PACE_REF
nopt=2
IF(OCEAN_CASE_SELECT==2)THEN
 nopt=3
 prep_wave(3)=665.0D0
ENDIF
CALL Atmosphere_Profile_READIN(CFILE_AP)
CALL Surface_Pressure_Rescale(PRESSURE_SURFACE)
CALL GAS_COLUMN_RESCALE(OZONE_COLUMN,H2O_COLUMN)
CALL NO2_COLUMN_RESCALE(NO2_COLUMN)
ALLOCATE(gc(nwave,NPGRID),gh(nwave,NPGRID),gco(nwave,NPGRID),gch(nwave,NPGRID), &
 go(nwave,NPGRID),goz(nwave,NPGRID),gn(nwave,NPGRID))
gc=0; gh=0; gco=0; gch=0; go=0; goz=0; gn=0
IF(GAS_ABS_FLAG) CALL PREPARE_GAS_ABSORPTION(nwave,wavelengths,gc,gh,gco,gch,go,goz,gn)
CALL SUNLREADIN(aux_dir,TRIM(solar_policy))
IF(ANY(wavelengths/1000.0D0<WVSUN(1)) .OR. ANY(wavelengths/1000.0D0>WVSUN(NSUNL))) &
 ERROR STOP 'Requested wavelength outside selected solar data'
DO iwave=1,nwave
 CALL FSUNCAL(wavelengths(iwave)/1000.0D0,esun(iwave))
ENDDO
ALLOCATE(MUOUT(NTHETAOUT),PHIOUT(NPHIOUT))
CALL MUPHIOUTASS(NTHETAOUT,NPHIOUT,MUOUT,PHIOUT)
MAXLORDINPUT=MAX(NQUADAINPUT,MAXMORDINPUT_SAVE)
MAXMORDINPUT=MAXMORDINPUT_SAVE
IF(MAXLORDINPUT>3000) ERROR STOP 'Phase expansion exceeds core MAXLFULL=3000'
WLR_FLAG_INPUT=.false.; WCFLG_INPUT=.true.; DELTAM_INPUT=2
nmbre_local=1.5D0; nmbim_local=0.0D0
IF(OCEAN_CASE_SELECT<0)THEN
 NTLYERO=0
ELSEIF(OCEAN_CASE_SELECT==0)THEN
 NTLYERO=1
ELSE
 NTLYERO=NWATER_DEPTH-1
ENDIF
NTLYER=NTLYERA+NTLYERO
IF(OCEAN_CASE_SELECT<0)THEN
 NREC=NTLYER+NDET+1
 NUMMIERT=NTLYERA
ELSE
 NREC=NTLYER+NDET+2
 IF(OCEAN_PHMX_ONE)THEN
  NUMMIERT=NTLYERA+1
 ELSE
  NUMMIERT=NTLYERA+NTLYERO
 ENDIF
ENDIF
ALLOCATE(TAUR(NTLYERA),DEPOL_A(NTLYERA),TAU_TG(NTLYERA),TAU_TG_H2O(NTLYERA), &
 TAU_TG_CO2(NTLYERA),TAU_TG_CH4(NTLYERA),TAU_TG_O2(NTLYERA),TAU_TG_O3(NTLYERA),TAU_TG_NO2(NTLYERA))
ALLOCATE(TAU_ARSL_Ext_TMP(NTLYERA),TAU_ARSL_Scat_TMP(NTLYERA))
ALLOCATE(NUMMIEANGINPUT(NUMMIERT),PHMXDATASTREAM(NUMMIERT,NUMMIEANGMAX,0:6))
ALLOCATE(DIRAD(NDET,2),DSTOKES(NDET,NTHETAOUT,NPHIOUT,4))
ALLOCATE(RECDATASTREAM(NREC,7),PNDLY(NTLYERA),ARSL_TYPE(NTLYERA))
ALLOCATE(CHLA_VERTICAL_PROFILE(NWATER_DEPTH))
ALLOCATE(TAU_REDUCE_ATMOS_PRESERVE(NTLYERA+1),TAU_REDUCE_OCEAN_PRESERVE(NTLYERO+1))
CALL WATER_DEPTH_SETUP
CALL PHMXMIE_ALLO(NUMMIEUSE,OCEAN_CASE_SELECT,nopt,prep_wave(1:nopt))
CALL prepare_particles(nopt)
CALL PREPARE_PARTICLE_COEFFICIENTS(NUMMIEUSE,MAXLORDINPUT,DELTAM_INPUT,1,1)
CEXT_REF=0.0D0
CEXT_REF(1:NUMMIEUSE)=CEXT_ARSL_BAND(1:NUMMIEUSE,2)
CALL PREPARE_AEROSOL_PROFILE(IAEROSOL,NTLYERA,PNDLY,ARSL_TYPE)
IF(OCEAN_CASE_SELECT>=0)THEN
 CALL WATER_ABSORPTION_READ
 CALL WATER_ABSORPTION_LEE_READ
 CALL WATER_ABSORPTION_MCF_READ
 CALL WATER_ABSORPTION_IOCCG_READ
ENDIF
ALLOCATE(irrad_all(NDET,2,nwave),stokes_all(NDET,NTHETAOUT,NPHIOUT,4,nwave))
ALLOCATE(taur_all(NTLYERA,nwave),gas_all(NTLYERA,nwave),ext_all(NTLYERA,nwave),scat_all(NTLYERA,nwave))
CALL CPU_TIME(prepared_time)
WRITE(*,*) 'Shared scene preparation count: 1; wavelengths: ',nwave
DO iwave=1,nwave
 requested_nm=wavelengths(iwave)
 esun_value=esun(iwave)
 IF(iwave>1)THEN
  prep_wave(1)=requested_nm
  OPTICAL_WAVELENGTHS(1)=requested_nm
  ! Only the current optical slot changes. Reference cross sections persist;
  ! their phase matrices are not consumed by a requested-wavelength RT solve.
  CALL prepare_particles(1)
  CALL PREPARE_PARTICLE_COEFFICIENTS(NUMMIEUSE,MAXLORDINPUT,DELTAM_INPUT,1,1)
 ENDIF
 WAVELENGTH_MICRON=requested_nm/1000.0D0
CALL TAURCAL(WAVELENGTH_MICRON,NTLYERA,TAUR,DEPOL_A)
CALL PREPARE_GAS_LAYERS(GAS_ABS_FLAG,NTLYERA,nwave,iwave,gc,gh,gco,gch,go,goz,gn, &
 TAU_TG,TAU_TG_H2O,TAU_TG_CO2,TAU_TG_CH4,TAU_TG_O2,TAU_TG_O3,TAU_TG_NO2)
CALL PREPARE_RT(0,1,WAVELENGTH_MICRON,NREC,NTLYERA,NTLYERO,PNDLY,ARSL_TYPE, &
 TAUR,DEPOL_A,TAU_TG,TAU_ARSL_Ext_TMP,TAU_ARSL_Scat_TMP,RECDATASTREAM,NTLYER,NUMMIERT, &
 MAXLORDINPUT,NUMMIEANGINPUT,PHMXDATASTREAM,DELTAM_INPUT,LBDO_LD,FLAM,CHLa,TEMPERTURE,SALINITY, &
 pBRDFa_input,pBRDFk_input,pBRDFb_input,pBRDFe_input,1,CEXT_ARSL_BAND(:,1),CSCAT_ARSL_BAND(:,1))
CALL RTSOS(NREC,NDET,RECDATASTREAM,NCOLINPUT,NQUADAINPUT,NQUADOINPUT, &
 0.05D0,0.1D0,NUMMIERT,MAXLORDINPUT,MAXMORDINPUT,NUMMIEANGINPUT,PHMXDATASTREAM, &
 MU_IN,NTHETAOUT,NPHIOUT,PHIOUT,MUOUT,DIRAD,DSTOKES,WLR_FLAG_INPUT,WCFLG_INPUT,DELTAM_INPUT)
! RTSOS elastic output uses its unit solar beam. Convert to the selected spectral irradiance.
scale_value=esun_value*ABS(MU_IN)/DIRAD(1,1)
DIRAD=DIRAD*scale_value
DSTOKES=DSTOKES*scale_value
! Restore attenuation removed by PREPARE_RT, as in the legacy PACE wavelength loop.
IF(ABSORPTION_TAU_REDUCE_FLAG)THEN
 DIRAD(2:NDET,:)=DIRAD(2:NDET,:)*EXP(-TAU_REDUCE_ATMOS_PRESERVE(NTLYERA+1))
 DSTOKES(2:NDET,:,:,:)=DSTOKES(2:NDET,:,:,:)*EXP(-TAU_REDUCE_ATMOS_PRESERVE(NTLYERA+1))
ENDIF
 irrad_all(:,:,iwave)=DIRAD
 stokes_all(:,:,:,:,iwave)=DSTOKES
 taur_all(:,iwave)=TAUR
 gas_all(:,iwave)=TAU_TG
 ext_all(:,iwave)=TAU_ARSL_Ext_TMP
 scat_all(:,iwave)=TAU_ARSL_Scat_TMP
ENDDO
CALL CPU_TIME(finish_time)
WRITE(*,*) 'Scene preparation CPU seconds: ',prepared_time-start_time
WRITE(*,*) 'Wavelength evaluation CPU seconds: ',finish_time-prepared_time
CALL write_result(OUTFILE,nwave,wavelengths,esun,MAXMORDINPUT,MAXLORDINPUT, &
 NDET,NTHETAOUT,NPHIOUT,NTLYERA,MUOUT,PHIOUT,taur_all,gas_all,ext_all,scat_all,irrad_all,stokes_all)
CONTAINS
SUBROUTINE prepare_particles(last_slot)
 INTEGER,INTENT(IN) :: last_slot
 DOUBLE PRECISION :: rr1,vv1,rr2,vv2,nd1,nd2,mr1(nopt),mi1(nopt),mr2(nopt),mi2(nopt)
 INTEGER :: kind
 CALL PHMXMIE_INIT(1,IAEROSOL,NUMMIEUSE,AeroFMF,RH_simu,nopt,rr1,vv1,rr2,vv2, &
  WV_PACE_REF,mr1,mi1,mr2,mi2,nd1,nd2,PHYTOPLANKTON_INDEX_REFRACTION, &
  PHYTOPLANKTON_SPECTRAL_SLOPE,SEDIMENT_INDEX_REFRACTION,SEDIMENT_SPECTRAL_SLOPE, &
  prep_wave(1:nopt),1,last_slot)
 IF(OCEAN_CASE_SELECT==2)THEN
  DO kind=2,3
   CALL PHMXMIE_INIT(kind,IAEROSOL,NUMMIEUSE,AeroFMF,RH_simu,nopt,rr1,vv1,rr2,vv2, &
    WV_PACE_REF,mr1,mi1,mr2,mi2,nd1,nd2,PHYTOPLANKTON_INDEX_REFRACTION, &
    PHYTOPLANKTON_SPECTRAL_SLOPE,SEDIMENT_INDEX_REFRACTION,SEDIMENT_SPECTRAL_SLOPE, &
    prep_wave(1:nopt),1,last_slot)
  ENDDO
 ENDIF
END SUBROUTINE
END PROGRAM MONOCHROMATIC
