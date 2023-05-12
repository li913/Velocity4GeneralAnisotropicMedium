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

  module Mod_Medium
  implicit none
  contains


  elemental subroutine normal2Angle(x,y,z,theta,phi)
  implicit none
  !//**********************************************************************
  !// trans vector to 2 angles
  !//**********************************************************************
  real(RP), intent(in):: x, y, z
  real(RP), intent(out):: theta
  real(RP), intent(out):: phi
  theta = atan2d(sqrt(x*x+y*y),z)
  theta = abs(theta)
  phi   = atan2d(y,x)
  if(phi<0.0_RP) phi = 360.0_RP + phi
  end subroutine

  pure function mediumParTrans2Aij(n, c, Q) result(res)
  implicit none
  !//**********************************************************************
  !// trans 1D medium parameters to secend-order matrix
  !//**********************************************************************
  integer, intent(in):: n
  DATATYPE, intent(in):: c(n)
  real(RP), optional, intent(in):: Q(n)
  complex(RP) res(6,6), a(n)
  integer i, j
  res = 0.0_RP
  a = c
  if(present(Q)) then
    do i = 1, n
      a(i)%im = -a(i)%re/Q(i)
    end do
  end if

  select case (n)
  case (2) !isotropic
    res(1,1) = a(1); res(2,2) = a(1); res(3,3) = a(1); res(4,4) = a(2); res(5,5) = a(2); res(6,6) = a(2)
    res(1,2) = a(1) - 2.0*a(2); res(1,3) = res(1,2); res(2,3) = res(1,2)
  case (5) !VTI
    res(1,1) = a(1); res(1,3) = a(2); res(3,3) = a(3); res(4,4) = a(4); res(6,6) = a(5)
    res(1,2)=res(1,1)-2*res(6,6); res(2,2)=res(1,1); res(2,3)=res(1,3); res(5,5)=res(4,4)
    !res(3,1)=res(1,3); res(3,2)=res(1,3); res(2,1)=res(1,2)
  case (9) !orthorhombic
    res(1,1) = a(1); res(1,2) = a(2); res(1,3) = a(3); res(2,2) = a(4); res(2,3) = a(5)
    res(3,3) = a(6); res(4,4) = a(7); res(5,5) = a(8); res(6,6) = a(9)
  case (13) !Monoclinic
    res(1,1:3) = a(1:3); res(1,6) = a(4); res(2,2:3) = a(5:6); res(2,6) = a(7); res(3,3) = a(8)
    res(3,6) = a(9); res(4,4:5) = a(10:11); res(5,5) = a(12); res(6,6) = a(13)
  case (21) !Triclinic
    res(1,1:6) = a(1:6); res(2,2:6) = a(7:11); res(3,3:6) = a(12:15)
    res(4,4:6) = a(16:18); res(5,5:6) = a(19:20); res(6,6) = a(21)
  end select
  do i = 1, 6
    do j = i+1, 6
      res(j,i) = res(i,j)
    end do
  end do
  end function

  end module


