MODULE MONOCHROMATIC_RESULT
CONTAINS
SUBROUTINE write_result(OUTFILE,nwave,requested_nm,esun_value,MAXMORDINPUT,MAXLORDINPUT, &
 NDET,NTHETAOUT,NPHIOUT,NTLYERA,MUOUT,PHIOUT,TAUR,TAU_TG,TAU_ARSL_Ext_TMP,TAU_ARSL_Scat_TMP,DIRAD,DSTOKES, &
 observation_theta,observation_phi)
 USE RTUTILITY, ONLY : PI
 USE HDF5
 USE H5DS
 USE HDF5_UTILS, ONLY : CHECK_HDF5
 IMPLICIT NONE
 CHARACTER(LEN=*) :: OUTFILE
 INTEGER :: nwave,MAXMORDINPUT,MAXLORDINPUT,NDET,NTHETAOUT,NPHIOUT,NTLYERA
 INTEGER :: observation_theta(:),observation_phi(:),nobs
 DOUBLE PRECISION :: requested_nm(nwave),esun_value(nwave),MUOUT(NTHETAOUT),PHIOUT(NPHIOUT), &
 TAUR(NTLYERA,nwave),TAU_TG(NTLYERA,nwave),TAU_ARSL_Ext_TMP(NTLYERA,nwave),TAU_ARSL_Scat_TMP(NTLYERA,nwave), &
 DIRAD(NDET,2,nwave),DSTOKES(:,:,:,:,:)
 INTEGER(HID_T) :: file,theta_id,phi_id,stokes_id,layer_id,wave_id,observation_id
 INTEGER :: err,det
 CHARACTER(LEN=10) :: labels(4)=[CHARACTER(LEN=10)::'TOA','AirSensor','BOA','TOO']
 CALL h5open_f(err)
 CALL h5fcreate_f(TRIM(OUTFILE),H5F_ACC_TRUNC_F,file,err)
 CALL CHECK_HDF5(err,'create monochromatic result',OUTFILE)
 CALL vector('evaluation',[(DBLE(det-1),det=1,nwave)],wave_id)
 CALL spectrum('Wavelength_Evaluated',requested_nm)
 CALL spectrum('Solar_Irradiance',esun_value)
 CALL scalar('rt_fourier_order_limit',DBLE(MAXMORDINPUT))
 CALL scalar('phase_expansion_order',DBLE(MAXLORDINPUT))
 nobs=SIZE(observation_theta)
 IF(nobs>0)THEN
  CALL vector('observation',[(DBLE(det-1),det=1,nobs)],observation_id)
  CALL vector('ThetaV',ACOS(MUOUT(observation_theta))/PI*180.0D0,theta_id,.false.)
  CALL vector('PhiV',PHIOUT(observation_phi)/PI*180.0D0,phi_id,.false.)
  CALL h5dsattach_scale_f(theta_id,observation_id,1,err)
  CALL h5dsattach_scale_f(phi_id,observation_id,1,err)
 ELSE
  CALL vector('ThetaV',ACOS(MUOUT)/PI*180.0D0,theta_id)
  CALL vector('PhiV',PHIOUT/PI*180.0D0,phi_id)
 ENDIF
 CALL vector('N_Stokes_Parameter',[1.0D0,2.0D0,3.0D0,4.0D0],stokes_id)
 CALL vector('Atmosphere_Layer',[(DBLE(det),det=1,NTLYERA)],layer_id)
 CALL field('Tau_Rayleigh_Extinction',TAUR,layer_id)
 CALL field('Tau_Gas_Absorption',TAU_TG,layer_id)
 CALL field('Tau_Aerosol_Extinction',TAU_ARSL_Ext_TMP,layer_id)
 CALL field('Tau_Aerosol_Scattering',TAU_ARSL_Scat_TMP,layer_id)
 DO det=1,NDET
  CALL spectrum('Irradiance_'//TRIM(labels(det))//'_Downwelling',DIRAD(det,1,:))
  CALL spectrum('Irradiance_'//TRIM(labels(det))//'_Upwelling',DIRAD(det,2,:))
  IF(nobs>0)THEN
   CALL observed_stokes('Stokes_Vector_'//TRIM(labels(det)),DSTOKES(det,:,1,:,:))
  ELSE
   CALL stokes('Stokes_Vector_'//TRIM(labels(det)),DSTOKES(det,:,:,:,:))
  ENDIF
 ENDDO
 CALL h5dclose_f(theta_id,err); CALL h5dclose_f(phi_id,err)
 IF(nobs>0) CALL h5dclose_f(observation_id,err)
 CALL h5dclose_f(stokes_id,err); CALL h5dclose_f(layer_id,err); CALL h5dclose_f(wave_id,err)
 CALL h5fclose_f(file,err); CALL h5close_f(err)
 CONTAINS
 SUBROUTINE scalar(name,value)
  CHARACTER(LEN=*) :: name
  DOUBLE PRECISION :: value
  INTEGER(HID_T) :: space,dset
  INTEGER(HSIZE_T) :: dims(1)=[1]
  CALL h5screate_f(H5S_SCALAR_F,space,err)
  CALL h5dcreate_f(file,name,H5T_NATIVE_DOUBLE,space,dset,err)
  CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,value,dims,err)
  CALL CHECK_HDF5(err,'write scalar',name)
  CALL h5dclose_f(dset,err); CALL h5sclose_f(space,err)
 END SUBROUTINE
 SUBROUTINE spectrum(name,values)
  CHARACTER(LEN=*) :: name
  DOUBLE PRECISION :: values(:)
  INTEGER(HID_T) :: dset
  CALL vector(name,values,dset,.false.)
  CALL h5dsattach_scale_f(dset,wave_id,1,err)
  CALL h5dclose_f(dset,err)
 END SUBROUTINE
 SUBROUTINE vector(name,values,dset,scale)
  LOGICAL,OPTIONAL :: scale
  CHARACTER(LEN=*) :: name
  DOUBLE PRECISION :: values(:)
  INTEGER(HID_T) :: dset,space
  INTEGER(HSIZE_T) :: dims(1)
  dims=SIZE(values)
  CALL h5screate_simple_f(1,dims,space,err)
  CALL h5dcreate_f(file,name,H5T_NATIVE_DOUBLE,space,dset,err)
  CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,values,dims,err)
  CALL CHECK_HDF5(err,'write vector',name)
  IF(PRESENT(scale))THEN
   IF(scale) CALL h5dsset_scale_f(dset,err)
  ELSE
   CALL h5dsset_scale_f(dset,err)
  ENDIF
  CALL h5sclose_f(space,err)
 END SUBROUTINE
 SUBROUTINE field(name,values,dim)
  CHARACTER(LEN=*) :: name
  DOUBLE PRECISION :: values(:,:)
  INTEGER(HID_T) :: dim,dset,space
  INTEGER(HSIZE_T) :: dims(2)
  dims=SHAPE(values)
  CALL h5screate_simple_f(2,dims,space,err)
  CALL h5dcreate_f(file,name,H5T_NATIVE_DOUBLE,space,dset,err)
  CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,values,dims,err)
  CALL CHECK_HDF5(err,'write field',name)
  CALL h5dsattach_scale_f(dset,dim,2,err)
  CALL h5dsattach_scale_f(dset,wave_id,1,err)
  CALL h5dclose_f(dset,err); CALL h5sclose_f(space,err)
 END SUBROUTINE
 SUBROUTINE stokes(name,values)
  CHARACTER(LEN=*) :: name
  DOUBLE PRECISION :: values(:,:,:,:)
  INTEGER(HID_T) :: space,dset
  INTEGER(HSIZE_T) :: dims(4)
  dims=SHAPE(values)
  CALL h5screate_simple_f(4,dims,space,err)
  CALL h5dcreate_f(file,name,H5T_NATIVE_DOUBLE,space,dset,err)
  CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,values,dims,err)
  CALL CHECK_HDF5(err,'write Stokes',name)
  CALL h5dsattach_scale_f(dset,theta_id,4,err)
  CALL h5dsattach_scale_f(dset,phi_id,3,err)
  CALL h5dsattach_scale_f(dset,stokes_id,2,err)
  CALL h5dsattach_scale_f(dset,wave_id,1,err)
  CALL h5dclose_f(dset,err); CALL h5sclose_f(space,err)
 END SUBROUTINE
 SUBROUTINE observed_stokes(name,values)
  CHARACTER(LEN=*) :: name
  DOUBLE PRECISION :: values(:,:,:),ordered(4,nobs,nwave)
  INTEGER :: component
  INTEGER(HID_T) :: space,dset
  INTEGER(HSIZE_T) :: dims(3)
  DO component=1,4
   ordered(component,:,:)=values(:,component,:)
  ENDDO
  dims=SHAPE(ordered)
  CALL h5screate_simple_f(3,dims,space,err)
  CALL h5dcreate_f(file,name,H5T_NATIVE_DOUBLE,space,dset,err)
  CALL h5dwrite_f(dset,H5T_NATIVE_DOUBLE,ordered,dims,err)
  CALL CHECK_HDF5(err,'write observed Stokes',name)
  CALL h5dsattach_scale_f(dset,stokes_id,3,err)
  CALL h5dsattach_scale_f(dset,observation_id,2,err)
  CALL h5dsattach_scale_f(dset,wave_id,1,err)
  CALL h5dclose_f(dset,err); CALL h5sclose_f(space,err)
 END SUBROUTINE
END SUBROUTINE write_result
END MODULE MONOCHROMATIC_RESULT
