#ifdef DOUBLE    
#define RP 8
#else
#define RP 4
#endif

#ifdef COMPLEX
#define DATATYPE complex(RP)
#else
#define DATATYPE real(RP)
#endif


  program liVelocity
  use Mod_Velocity4GeneralAnisotropicMedium
  use Mod_Medium
  implicit none
  integer, parameter:: nc = 1
  integer, parameter:: nTheta=180*nc, nPhi = 360*nc
  real(RP), parameter:: DEG2RAD = acos(-1.0_RP)/180.0_RP 
  ! phase velocity, ray velocity, ray velocity components
  DATATYPE, allocatable:: pv(:,:,:), rv(:,:,:), vi(:,:,:,:)
  DATATYPE v(0:180,0:360)
  ! partial derivatives of group (or phase) velocity with respect to 21 elastic parameters
  DATATYPE pv_a(6,6,3), rv_a(6,6,3)
  ! elastic parameters
  DATATYPE aa(21), a(6,6)
  DATATYPE theta, phi, n(3)
  integer waveType, nPar, i, j
  real(8), allocatable:: theta_ray(:,:), phi_ray(:,:)

  allocate(pv(3,0:nTheta,0:nPhi), rv(3,0:nTheta,0:nPhi), vi(3,3,0:nTheta,0:nPhi))
  allocate(theta_ray(0:nTheta,0:nPhi), phi_ray(0:nTheta,0:nPhi))

  nPar=5
  aa(1:5)=[15.064d0, 1.639d0, 10.837d0, 3.126d0, 4.251d0] ! VTI medium

  nPar = 21
  !aa=[24.1,11.7,13.5,-0.5,0.6,-0.5,   23.5,12.6,-0.2,0.0,1.1,  26.5,0.1,-0.4,-0.2,   5.9,0.0,0.2,   5.2,-0.3,  5.5]/2.19
  !aa=[10.3,0.9,1.3,1.4,1.1,0.8,  10.6,2.1,0.2,-0.2,-0.6,  14.1,0.0,-0.5,-1.0,  5.1,0.0,0.2,  6.0,0.0,  4.9]/2.08
  !aa=[86.0,7.4,11.91,-18.04,0.0,0.0,  86.0,11.91,18.04,0.0,0.0,  105.75,0.0,0.0,0.0,  58.2,0.0,0.0,   58.2,-18.04,  39.3]/2.65 ! Li fang
  aa=[45.597,0.935,0.895,0.103,0.074,0.07,  14.442,0.887,-0.083,0.018,-0.059,  44.303,-0.049,-0.04,-0.026, 0.459,0.09,0.052,  0.374,0.101,  0.45] ! Grechka (2017)

  !nPar = 9
  !aa(1:nPar)=[8.56,4.9,5.13,11.29,5.13,12.76,2.79,2.57,2.27]
  !aa(1:nPar)=[15.1,5.56,1.64,12.1,1.04,10.84,3.13,2.5,4.25]

  waveType = qsv
  a = mediumParTrans2Aij(nPar, aa(1:nPar))

  !$OMP parallel do num_threads(5) private(i,j,theta,phi,n, pv_a, rv_a)
  do i = 0, nPhi
  phi = i*DEG2RAD/nc
    do j = 0, nTheta
       theta = j*DEG2RAD/nc
      n = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
      !call phaseAndRayVelocity_mkl(a, n, pv(:,j,i), rv(:,j,i), vi(:,:,j,i)) ! eigenvector method
      call phaseAndRayVelocity(a, n, pv(:,j,i), rv(:,j,i), vi(:,:,j,i), pv_a, rv_a)
      ! calculate ray angle
      call normal2Angle(vi(1,waveType,j,i),vi(2,waveType,j,i),vi(3,waveType,j,i),theta_ray(j,i),phi_ray(j,i))
    end do
  end do
  !$OMP end parallel do

  ! output result
  open(10,file='v.csv')
  do i = 0, nPhi
    do j = 0, nTheta
      !write(10,"(*(g0,','))") phi_ray(j,i),theta_ray(j,i), -1.0_RP/rv(waveType,j,i)/rv(waveType,j,i) *  [rv_a(1,1:,waveType),rv_a(2,2:,waveType),rv_a(3,3:,waveType),rv_a(4,4:,waveType),rv_a(5,5:,waveType),rv_a(6,6:,waveType)]
      !write(10,"(*(g0,','))") phi_ray(j,i), theta_ray(j,i), real(rv(waveType,j,i))
      !write(10,"(*(g0,','))") rv_a_TI(rv_a(:,:,waveType), rv(waveType,j,i)) ! derivations for TI media
      !write(10,"(*(g0,','))") phi_ray(j,i), theta_ray(j,i), real(rv(waveType,j,i))
      write(10,"(*(g0,','))") real(vi(:,waveType,j,i)), real(rv(waveType,j,i))
    end do
  end do
  close(10)

 
  contains
  function rv_a_TI(rv_a, rv) result(res)
  implicit none
  DATATYPE rv_a(6,6), rv, res(5)
  res(1) = rv_a(1,1) + rv_a(1,2) + rv_a(2,2)
  res(2) = rv_a(1,3) + rv_a(2,3)
  res(3) = rv_a(3,3)
  res(4) = rv_a(4,4) + rv_a(5,5)
  res(5) = rv_a(6,6) -2*rv_a(1,2)
  res = -res / (rv*rv)
  end function

  end program liVelocity

