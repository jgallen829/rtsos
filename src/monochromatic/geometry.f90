MODULE MONOCHROMATIC_GEOMETRY
! The solver consumes Cartesian axes. Evaluate only unique requested axis values,
! then gather the selected pairs before storing wavelength results. This shares
! one multiple-scattering solve per wavelength; angular workspace is ntheta*nphi.
IMPLICIT NONE
CONTAINS
SUBROUTINE read_geometry(unit,ntheta,nphi,mu,phi,theta_index,phi_index)
 USE RTUTILITY, ONLY : PI,EXPLICIT_DIRECTIONS
 USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE
 INTEGER,INTENT(IN) :: unit
 INTEGER,INTENT(INOUT) :: ntheta,nphi
 DOUBLE PRECISION,ALLOCATABLE,INTENT(OUT) :: mu(:),phi(:)
 INTEGER,ALLOCATABLE,INTENT(OUT) :: theta_index(:),phi_index(:)
 CHARACTER(LEN=40) :: mode
 INTEGER :: status,nobs,i
 DOUBLE PRECISION,ALLOCATABLE :: theta_requested(:),phi_requested(:)
 READ(unit,'(A)',IOSTAT=status) mode
 EXPLICIT_DIRECTIONS=.false.
 IF(status<0)THEN
  ! Existing requested-wavelength records end after the wavelengths.
  ALLOCATE(mu(ntheta),phi(nphi),theta_index(0),phi_index(0))
  CALL MUPHIOUTASS(ntheta,nphi,mu,phi)
  RETURN
 ENDIF
 IF(status/=0 .OR. TRIM(mode)/='explicit_directions') ERROR STOP 'Unknown native geometry record'
 READ(unit,*) nobs
 IF(nobs<1) ERROR STOP 'At least one paired observation is required'
 ALLOCATE(theta_requested(nobs),phi_requested(nobs),theta_index(nobs),phi_index(nobs))
 DO i=1,nobs
  READ(unit,*) theta_requested(i),phi_requested(i)
 ENDDO
 IF(ANY(.NOT.IEEE_IS_FINITE(theta_requested)) .OR. ANY(.NOT.IEEE_IS_FINITE(phi_requested))) &
  ERROR STOP 'Viewing directions must be finite'
 IF(ANY(theta_requested<0.0D0) .OR. ANY(theta_requested>180.0D0)) &
  ERROR STOP 'Propagation zenith must be in [0,180] degrees'
 IF(ANY(theta_requested==90.0D0)) &
  ERROR STOP 'Exact horizon needs the mu=0 transport/direct-beam limit'
 theta_requested=COS(theta_requested/180.0D0*PI)
 phi_requested=MODULO(phi_requested,360.0D0)/180.0D0*PI
 CALL exact_axis(theta_requested,mu,theta_index)
 CALL exact_axis(phi_requested,phi,phi_index)
 ntheta=SIZE(mu); nphi=SIZE(phi)
 EXPLICIT_DIRECTIONS=.true.
 WRITE(*,*) 'Ordered observations: ',nobs,'; internal requested axes: ',ntheta,nphi
END SUBROUTINE

SUBROUTINE exact_axis(requested,axis,indices)
 DOUBLE PRECISION,INTENT(IN) :: requested(:)
 DOUBLE PRECISION,ALLOCATABLE,INTENT(OUT) :: axis(:)
 INTEGER,INTENT(OUT) :: indices(:)
 DOUBLE PRECISION :: sorted(SIZE(requested)),work(SIZE(requested))
 INTEGER :: i,n,found(1)
 sorted=requested; work=requested
 CALL qsort2a(SIZE(requested),work,sorted)
 n=1
 DO i=2,SIZE(sorted)
  IF(sorted(i)==sorted(n)) CYCLE
  n=n+1
  sorted(n)=sorted(i)
 ENDDO
 axis=sorted(:n)
 DO i=1,SIZE(requested)
  ! An exact member lookup, never nearest-angle selection.
  found=FINDLOC(axis,requested(i))
  indices(i)=found(1)
 ENDDO
END SUBROUTINE
END MODULE
