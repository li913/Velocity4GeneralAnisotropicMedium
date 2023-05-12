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

  module Mod_Velocity4GeneralAnisotropicMedium
  implicit none
  enum, bind(c)
    enumerator RAY_VELOCITY, RAY_SLOWNESS, RAY_ATTENUATION, IMAGINARY_VELOCITY
    enumerator RAY_QUALITY, PHASE_QUALITY, PHASE_ATTENUATION, PHASE_VELOCITY
  end enum
  enum, bind(c)
    enumerator:: qP=1, qSV, qSH
  end enum
  enum, bind(c)
    enumerator thetaPositive, thetaNegative, phiPositive, phiNegative
  end enum
  contains
  
  pure function christoffelMatrix(a, n) result(res)
  implicit none
  !//**********************************************************************
  !// christoffel matrix
  !//**********************************************************************
  DATATYPE, intent(in):: a(6,6), n(3)
  DATATYPE nn(6), res(6)
  !integer, parameter:: IND(3,3) = reshape([1,6,5,6,2,4,5,4,3],[3,3]) !convert Aijkl into Amn
  nn = [n(1)*n(1), n(1)*n(2), n(1)*n(3), n(2)*n(2), n(2)*n(3), n(3)*n(3)]
  res(1) = a(1,1)*nn(1) +   2.0_RP*a(1,6)*nn(2) +   2.0_RP*a(1,5)*nn(3) + a(6,6)*nn(4) +   2.0_RP*a(5,6)*nn(5) + a(5,5)*nn(6)
  res(2) = a(1,6)*nn(1) + (a(1,2)+a(6,6))*nn(2) + (a(1,4)+a(5,6))*nn(3) + a(2,6)*nn(4) + (a(2,5)+a(4,6))*nn(5) + a(4,5)*nn(6)
  res(3) = a(1,5)*nn(1) + (a(1,4)+a(5,6))*nn(2) + (a(1,3)+a(5,5))*nn(3) + a(4,6)*nn(4) + (a(3,6)+a(4,5))*nn(5) + a(3,5)*nn(6)
  res(4) = a(6,6)*nn(1) +   2.0_RP*a(2,6)*nn(2) +   2.0_RP*a(4,6)*nn(3) + a(2,2)*nn(4) +   2.0_RP*a(2,4)*nn(5) + a(4,4)*nn(6)
  res(5) = a(5,6)*nn(1) + (a(2,5)+a(4,6))*nn(2) + (a(3,6)+a(4,5))*nn(3) + a(2,4)*nn(4) + (a(2,3)+a(4,4))*nn(5) + a(3,4)*nn(6)
  res(6) = a(5,5)*nn(1) +   2.0_RP*a(4,5)*nn(2) +   2.0_RP*a(3,5)*nn(3) + a(4,4)*nn(4)  +  2.0_RP*a(3,4)*nn(5) + a(3,3)*nn(6)
  end function

  pure function christoffel_n(a, n) result(res)
  implicit none
  !//**********************************************************************
  !// derivative of christoffel with respect to slowness direction n
  !//**********************************************************************
  DATATYPE, intent(in):: a(6,6), n(3)
  DATATYPE res(3,6), mat(6,3), B(6), n2(3)
  n2 = 2.0_RP * n
  mat(:,1) = [a(1,1), a(1,6), a(1,5), a(6,6), a(5,6), a(5,5)]
  mat(:,2) = [a(1,6), 0.5_RP*(a(1,2)+a(6,6)), 0.5_RP*(a(1,4)+a(5,6)), a(2,6), 0.5_RP*(a(2,5)+a(4,6)), a(4,5)]
  B(:)     = [a(1,5), 0.5_RP*(a(1,4)+a(5,6)), 0.5_RP*(a(1,3)+a(5,5)), a(4,6), 0.5_RP*(a(3,6)+a(4,5)), a(3,5)]
  mat(:,3) = B(:)
  res(1,:) = matmul(mat,n2) !n1
  mat(:,1) = mat(:,2)
  mat(:,2) = [a(6,6), a(2,6), a(4,6), a(2,2), a(2,4), a(4,4)]
  mat(:,3) = [a(5,6), 0.5_RP*(a(2,5)+a(4,6)), 0.5_RP*(a(3,6)+a(4,5)), a(2,4), 0.5_RP*(a(2,3)+a(4,4)), a(3,4)]
  res(2,:) = matmul(mat,n2) !n2
  mat(:,1) = B(:)
  mat(:,2) = mat(:,3)
  mat(:,3) = [a(5,5), a(4,5), a(3,5), a(4,4), a(3,4), a(3,3)]
  res(3,:) = matmul(mat,n2) !n3
  end function

  pure function christoffel_a(n) result(res)
  implicit none
  !//**********************************************************************
  !//derivative of christoffel with respect to parameter a
  !//**********************************************************************
  DATATYPE, intent(in):: n(3)
  DATATYPE res(6,6,6), n11, n12, n13, n22, n23, n33
  integer i, j
  n11 = n(1) * n(1); n12 = n(1) * n(2); n13 = n(1) * n(3)
  n22 = n(2) * n(2); n23 = n(2) * n(3); n33 = n(3) * n(3)
  res = 0.0_RP
  !F11
  res(1,1,1) = n11; res(1,6,1) = 2.0_RP*n12; res(1,5,1) = 2.0_RP*n13; res(6,6,1) = n22
  res(5,6,1) = 2.0_RP*n23; res(5,5,1) = n33
  !F12 
  res(1,6,2) = n11; res(1,2,2) = n12; res(6,6,2) = n12; res(1,4,2) = n13; res(5,6,2) = n13
  res(2,6,2) = n22; res(2,5,2) = n23; res(4,6,2) = n23; res(4,5,2) = n33
  !F13
  res(1,5,3) = n11; res(1,4,3) = n12; res(5,6,3) = n12; res(1,3,3) = n13; res(5,5,3) = n13
  res(4,6,3) = n22; res(3,6,3) = n23; res(4,5,3) = n23; res(3,5,3) = n33
  !F22 
  res(6,6,4) = n11; res(2,6,4) = 2.0_RP*n12; res(4,6,4) = 2.0_RP*n13; res(2,2,4) = n22
  res(2,4,4) = 2.0_RP*n23; res(4,4,4) = n33
  !F23 
  res(5,6,5) = n11; res(2,5,5) = n12; res(4,6,5) = n12; res(3,6,5) = n13; res(4,5,5) = n13
  res(2,4,5) = n22; res(2,3,5) = n23; res(4,4,5) = n23; res(3,4,5) = n33
  !F33 
  res(5,5,6) = n11; res(4,5,6) = 2.0_RP*n12; res(3,5,6) = 2.0_RP*n13; res(4,4,6) = n22
  res(3,4,6) = 2.0_RP*n23; res(3,3,6) = n33

  do i = 1, 6
    do j = i+1, 6
      res(j,i,:) = res(i,j,:)
    end do
  end do
  end function

  pure function christoffel_n_a(n) result(res)
  implicit none
  !//**********************************************************************
  !//derivative of christoffel with respect to a and n
  !//**********************************************************************
  DATATYPE, intent(in):: n(3)
  DATATYPE res(6,6,3,6), n1, n2, n3
  DATATYPE, parameter:: zero = 0.0_RP
  integer i, j
  n1 = n(1); n2 = n(2); n3 = n(3)
  res = 0.0_RP
  !F11 
  res(1,1,:,1) = [n1, zero, zero]; res(1,5,:,1) = [n3, zero, n1]; res(1,6,:,1) = [n2, n1, zero]
  res(5,5,:,1) = [zero, zero, n3]; res(5,6,:,1) = [zero, n3, n2]; res(6,6,:,1) = [zero, n2, zero]
  res(:,:,:,1) = res(:,:,:,1) * 2.0_RP
  !F12 
  res(1,2,:,2) = [n2, n1, zero]; res(1,4,:,2) = [n3, zero, n1];          res(1,6,:,2) = [n1*2.0_RP, zero, zero]
  res(2,5,:,2) = [zero, n3, n2]; res(2,6,:,2) = [zero, n2*2.0_RP, zero]; res(4,5,:,2) = [zero, zero, n3*2.0_RP]
  res(4,6,:,2) = [zero, n3, n2]; res(5,6,:,2) = [n3, zero, n1];          res(6,6,:,2) = [n2, n1, zero]
  !F13 
  res(1,3,:,3) = [n3, zero, n1];          res(1,4,:,3) = [n2, n1, zero]; res(1,5,:,3) = [n1*2.0_RP, zero, zero]
  res(3,5,:,3) = [zero, zero, n3*2.0_RP]; res(3,6,:,3) = [zero, n3, n2]; res(4,5,:,3) = [zero, n3, n2]
  res(4,6,:,3) = [zero, n2*2.0_RP, zero]; res(5,5,:,3) = [n3, zero, n1]; res(5,6,:,3) = [n2, n1, zero]
  !F22 
  res(2,2,:,4) = [zero, n2, zero]; res(2,4,:,4) = [zero, n3, n2]; res(2,6,:,4) = [n2, n1, zero]
  res(4,4,:,4) = [zero, zero, n3]; res(4,6,:,4) = [n3, zero, n1]; res(6,6,:,4) = [n1, zero, zero]
  res(:,:,:,4) = res(:,:,:,4) * 2.0_RP
  !F23 
  res(2,3,:,5) = [zero, n3, n2];          res(2,4,:,5) = [zero, n2*2.0_RP, zero]; res(2,5,:,5) = [n2, n1, zero]
  res(3,4,:,5) = [zero, zero, n3*2.0_RP]; res(3,6,:,5) = [n3, zero, n1];          res(4,4,:,5) = [zero, n3, n2]
  res(4,5,:,5) = [n3, zero, n1];          res(4,6,:,5) = [n2, n1, zero];          res(5,6,:,5) = [n1*2.0_RP, zero, zero]
  !F33 
  res(3,3,:,6) = [zero, zero, n3]; res(3,4,:,6) = [zero, n3, n2]; res(3,5,:,6) = [n3, zero, n1]
  res(4,4,:,6) = [zero, n2, zero]; res(4,5,:,6) = [n2, n1, zero]; res(5,5,:,6) = [n1, zero, zero]
  res(:,:,:,6) = res(:,:,:,6) * 2.0_RP

  do i = 1, 6
    do j = i+1, 6
      res(j,i,:,:) = res(i,j,:,:)
    end do
  end do
  end function

  pure subroutine BCD(cstl, B, C, D)
  implicit none
  !//************************************************************************
  !//calculate B, C, D
  !//************************************************************************
  DATATYPE, intent(in) :: cstl(6)
  DATATYPE, intent(out):: B, C, D
  B = -(cstl(1)+cstl(4)+cstl(6))
  C = (cstl(1)+cstl(4))*cstl(6) + cstl(1)*cstl(4) - (cstl(2)*cstl(2) + cstl(3)*cstl(3) + cstl(5)*cstl(5))
  D = cstl(1)*cstl(5)*cstl(5) + cstl(4)*cstl(3)*cstl(3) + cstl(6)*cstl(2)*cstl(2) - &
    (cstl(1)*cstl(4)*cstl(6) + 2.0_RP*cstl(2)*cstl(3)*cstl(5))
  end subroutine

  pure subroutine BCD_n(cstl, cstl_n, B_n, C_n, D_n)
  implicit none
  !//************************************************************************
  !//derivative of B, C, D with respect to n
  !//************************************************************************
  DATATYPE, intent(in) :: cstl(6), cstl_n(3,6)
  DATATYPE, intent(out):: B_n(3), C_n(3), D_n(3)
  B_n(:) = -(cstl_n(:,1) + cstl_n(:,4) + cstl_n(:,6))
  C_n(:) = (cstl(4)+cstl(6))*cstl_n(:,1) + (cstl(1)+cstl(6))*cstl_n(:,4) + (cstl(1)+cstl(4))*cstl_n(:,6) - &
    2.0_RP*(cstl(2)*cstl_n(:,2) + cstl(3)*cstl_n(:,3) + cstl(5)*cstl_n(:,5))
  D_n(:) = (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_n(:,1) + 2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_n(:,2) + &
    2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_n(:,3) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_n(:,4) + &
    2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_n(:,5) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_n(:,6)
  end subroutine

  pure subroutine BCD_a(cstl, cstl_a, B_a, C_a, D_a)
  implicit none
  !//************************************************************************
  !//derivative of B, C, D with respect to a
  !//************************************************************************
  DATATYPE, intent(in) :: cstl(6), cstl_a(6,6,6)
  DATATYPE, intent(out):: B_a(6,6), C_a(6,6), D_a(6,6)
  B_a(:,:) = -(cstl_a(:,:,1) + cstl_a(:,:,4) + cstl_a(:,:,6))
  C_a(:,:) = (cstl(4)+cstl(6))*cstl_a(:,:,1) + (cstl(1)+cstl(6))*cstl_a(:,:,4) + (cstl(1)+cstl(4))*cstl_a(:,:,6) - &
    2.0_RP*(cstl(2)*cstl_a(:,:,2) + cstl(3)*cstl_a(:,:,3) + cstl(5)*cstl_a(:,:,5))
  D_a(:,:) = (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_a(:,:,1) + 2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_a(:,:,2) + &
    2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_a(:,:,3) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_a(:,:,4) + &
    2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_a(:,:,5) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_a(:,:,6)
  end subroutine

  pure subroutine BCD_n_a(cstl, cstl_n, cstl_a, cstl_n_a, B_n_a, C_n_a, D_n_a)
  implicit none
  !//************************************************************************
  !//derivative of B, C, D with respect to n and a
  !//************************************************************************
  DATATYPE, intent(in) :: cstl(6), cstl_n(3,6), cstl_a(6,6,6), cstl_n_a(6,6,3,6)
  DATATYPE, intent(out):: B_n_a(6,6,3), C_n_a(6,6,3), D_n_a(6,6,3)
  integer i
  B_n_a(:,:,:) = -(cstl_n_a(:,:,:,1) + cstl_n_a(:,:,:,4) + cstl_n_a(:,:,:,6))
  do concurrent(i=1:3)
    C_n_a(:,:,i) = (cstl_n(i,4)+cstl_n(i,6))*cstl_a(:,:,1) + (cstl(4)+cstl(6))*cstl_n_a(:,:,i,1) + &
      (cstl_n(i,1)+cstl_n(i,6))*cstl_a(:,:,4) + (cstl(1)+cstl(6))*cstl_n_a(:,:,i,4) + &
      (cstl_n(i,1)+cstl_n(i,4))*cstl_a(:,:,6) + (cstl(1)+cstl(4))*cstl_n_a(:,:,i,6) - &
      2.0_RP*( cstl_n(i,2)*cstl_a(:,:,2) + cstl(2)*cstl_n_a(:,:,i,2) + &
      cstl_n(i,3)*cstl_a(:,:,3) + cstl(3)*cstl_n_a(:,:,i,3) + &
      cstl_n(i,5)*cstl_a(:,:,5) + cstl(5)*cstl_n_a(:,:,i,5) )
    D_n_a(:,:,i) = (2.0_RP*cstl(5)*cstl_n(i,5)-cstl(6)*cstl_n(i,4)-cstl(4)*cstl_n(i,6))*cstl_a(:,:,1) + (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_n_a(:,:,i,1) + &
      (2.0_RP*cstl(3)*cstl_n(i,3)-cstl(6)*cstl_n(i,1)-cstl(1)*cstl_n(i,6))*cstl_a(:,:,4) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_n_a(:,:,i,4) + &
      (2.0_RP*cstl(2)*cstl_n(i,2)-cstl(4)*cstl_n(i,1)-cstl(1)*cstl_n(i,4))*cstl_a(:,:,6) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_n_a(:,:,i,6) + &
      2.0_RP*(cstl(6)*cstl_n(i,2)+cstl(2)*cstl_n(i,6)-cstl(5)*cstl_n(i,3)-cstl(3)*cstl_n(i,5))*cstl_a(:,:,2) + &
      2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_n_a(:,:,i,2) + &
      2.0_RP*(cstl(4)*cstl_n(i,3)+cstl(3)*cstl_n(i,4)-cstl(5)*cstl_n(i,2)-cstl(2)*cstl_n(i,5))*cstl_a(:,:,3) + &
      2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_n_a(:,:,i,3) + &
      2.0_RP*(cstl(5)*cstl_n(i,1)+cstl(1)*cstl_n(i,5)-cstl(3)*cstl_n(i,2)-cstl(2)*cstl_n(i,3))*cstl_a(:,:,5) + &
      2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_n_a(:,:,i,5)
  end do
  end subroutine
  
  pure function unaryCubicEquation(B, C, D) result(res)
  implicit none
  !//************************************************************************
  !//solve equation x**3 + Bx**2 + Cx + D = 0
  !//************************************************************************
  complex(RP), parameter:: w = cmplx(-0.5_RP,sqrt(0.75_RP),RP), w2 = w*w
  complex(RP), intent(in):: B, C, D
  complex(RP) res(3), u, v, P, Q
  integer sig
  sig = 1
  P = B*C/6.0_RP - B**3/27.0_RP - 0.5_RP*D
  Q = (3.0_RP*C-B*B) / 9.0_RP
  v = sqrt(P*P+Q**3)
  u = P + v
  v = P - v
  if(abs(v)>abs(u)) then
    u = v
    sig = -1
  end if
  u = u**(1.0_RP/3.0_RP)
  if(abs(u)>1.0e-15_RP) then
    v = -Q/u
  else
    v = 0.0_RP
  end if
  res = -B / 3.0_RP
  res(1) = res(1) + u    + v
  res(2) = res(2) + u*w2  + v*w
  res(3) = res(3) + u*w + v*w2
  if( any(abs(res**3 + B*res*res + C*res + D)>1.0e-10) ) error stop
  end function unaryCubicEquation

  pure subroutine phaseAndRayVelocity(a, n, pv, rv, vi, pv_a, rv_a)
  implicit none
  !//**********************************************************************
  !//calculate phase and group velocity, and their derivative wrt a
  !//**********************************************************************
  DATATYPE, intent(in):: a(6,6), n(3)
  DATATYPE, intent(out):: pv(3), rv(3), vi(3,3)
  DATATYPE, intent(out), optional:: pv_a(6,6,3), rv_a(6,6,3)
  real(RP), parameter:: EPS = 1.0e-7_RP
  DATATYPE cstl(6), cstl_n(3,6), cstl_a(6,6,6), cstl_n_a(6,6,3,6), vi_a(6,6,3,3)
  DATATYPE B, C, D, B_n(3), C_n(3), D_n(3), B_a(6,6), C_a(6,6), D_a(6,6), B_n_a(6,6,3), C_n_a(6,6,3), D_n_a(6,6,3)
  DATATYPE g(3), c2(3), c4(3)
  complex(RP) B0, C0, D0
  integer nRoot
  integer i, j
  !christoffel matrix
  cstl = christoffelMatrix(a, n)
  call BCD(cstl, B, C, D)
  !wrt n
  cstl_n = christoffel_n(a, n)
  call BCD_n(cstl, cstl_n, B_n, C_n, D_n)
  !wrt a
  if(present(pv_a) .or. present(rv_a)) then
    cstl_a = christoffel_a(n)
    call BCD_a(cstl, cstl_a, B_a, C_a, D_a)
    !wrt a and n
    if(present(rv_a)) then
      cstl_n_a = christoffel_n_a(n)
      call BCD_n_a(cstl, cstl_n, cstl_a, cstl_n_a, B_n_a, C_n_a, D_n_a)
    end if
  end if

  !phase velocity
#ifdef COMPLEX
  c2(:) = unaryCubicEquation(B, C, D)
#else
  B0 = B; C0 = C; D0 = D 
  c2(:) = unaryCubicEquation(B0, C0, D0)
#endif
  !qSV,qSH
  if(abs(a(6,6)*(n(1)*n(1)+n(2)*n(2))+a(4,4)*n(3)*n(3)-c2(2))<EPS) then
    c4(1) = c2(2); c2(2) = c2(3); c2(3) = c4(1)
  end if
  c4(:) = c2(:) * c2(:) 
  pv = sqrt(c2) 
  !num of roots
  if(abs(pv(2)-pv(3))<EPS) then
    nRoot = 2
    if(abs(pv(1)-pv(3))<EPS) nRoot = 1
  else
    nRoot = 3
  end if

  !group velocity
  select case(nRoot)
  case(1)
    g(:) = 6.0_RP*pv(:) 
    do concurrent(i=1:3)
      vi(:,i) = -B_n(:) / g(i)
    end do
  case(2)
    g(1) = 2.0_RP*pv(1)*(3.0_RP*c4(1) + 2.0_RP*B*c2(1) + C)
    vi(:,1) = -(c4(1)*B_n(:) + c2(1)*C_n(:) + D_n(:)) / g(1)
    g(2) = 4.0_RP*pv(2)*(3.0_RP*c2(2) + B) 
    vi(:,2) = -(2.0_RP*c2(2)*B_n(:) + C_n(:)) / g(2)
    g(3) = g(2); vi(:,3) = vi(:,2)
  case(3)
    g(:) = 2.0_RP*pv(:)*(3.0_RP*c4(:) + 2.0_RP*B*c2(:) + C) 
    do concurrent(i=1:3)
      vi(:,i) = -(c4(i)*B_n(:) + c2(i)*C_n(:) + D_n(:)) / g(i)
    end do
  end select
  do concurrent(i=1:3)
    rv(i) = sqrt(vi(1,i)*vi(1,i)+vi(2,i)*vi(2,i)+vi(3,i)*vi(3,i))
  end do

  !phase velocity wrt a
  if((.not.present(pv_a)) .and. (.not.present(rv_a))) return
  select case(nRoot)
  case(1)
    do concurrent(i=1:3)
      pv_a(:,:,i) = -B_a(:,:) / g(i)
    end do
  case(2)
    pv_a(:,:,1) = -(c4(1)*B_a(:,:) + c2(1)*C_a(:,:) + D_a(:,:)) / g(1)
    pv_a(:,:,2) = -(2.0_RP*c2(2)*B_a(:,:) + C_a(:,:)) / g(2)
    pv_a(:,:,3) = pv_a(:,:,2)
  case(3)
    do concurrent(i=1:3)
      pv_a(:,:,i) = -(c4(i)*B_a(:,:) + c2(i)*C_a(:,:) + D_a(:,:)) / g(i)
    end do
  end select

  !group velocity wrt a
  if(.not.present(rv_a)) return
  vi_a = 0.0_RP
  select case(nRoot)
  case(1)
    do j = 1, 3 
      do i = 1, 3 
        vi_a(:,:,i,j) = -(6.0_RP*pv_a(:,:,j)*vi(i,j)+B_n_a(:,:,i)) / g(j)
      end do
    end do
  case(2)
    j = 1 !qP
    do i = 1, 3 
      vi_a(:,:,i,j) = (30.0_RP*c4(j)+12.0_RP*B*c2(j)+2.0_RP*C)*pv_a(:,:,j) + 2.0_RP*pv(j)*(2.0_RP*c2(j)*B_a(:,:)+C_a(:,:))
      vi_a(:,:,i,j) = vi_a(:,:,i,j) * vi(i,j) 
      vi_a(:,:,i,j) = vi_a(:,:,i,j) + 2.0_RP*pv(j)*pv_a(:,:,j)*(2.0_RP*c2(j)*B_n(i)+C_n(i)) + c4(j)*B_n_a(:,:,i) + &
        c2(j)*C_n_a(:,:,i) + D_n_a(:,:,i) 
    end do
    vi_a(:,:,:,j) = -vi_a(:,:,:,j) / g(j) 
    j = 2 !qS1
    do i = 1, 3 
      vi_a(:,:,i,j) = (36.0_RP*c2(j)+4.0_RP*B)*pv_a(:,:,j) + 4.0_RP*pv(j)*B_a(:,:)
      vi_a(:,:,i,j) = vi_a(:,:,i,j) * vi(i,j) 
      vi_a(:,:,i,j) = vi_a(:,:,i,j) + 4.0_RP*pv(j)*pv_a(:,:,j)*B_n(i) + 2.0_RP*c2(j)*B_n_a(:,:,i) + C_n_a(:,:,i) 
    end do
    vi_a(:,:,:,j) = -vi_a(:,:,:,j) / g(j) 
    vi_a(:,:,:,3) = vi_a(:,:,:,2) !qS2
  case(3)
    do concurrent( j = 1: 3) 
      do i = 1, 3 
        vi_a(:,:,i,j) = (30.0_RP*c4(j)+12.0_RP*B*c2(j)+2.0_RP*C)*pv_a(:,:,j) + 2.0_RP*pv(j)*(2.0_RP*c2(j)*B_a(:,:)+C_a(:,:))
        vi_a(:,:,i,j) = vi_a(:,:,i,j) * vi(i,j) 
        vi_a(:,:,i,j) = vi_a(:,:,i,j) + 2.0_RP*pv(j)*pv_a(:,:,j)*(2.0_RP*c2(j)*B_n(i)+C_n(i)) + c4(j)*B_n_a(:,:,i) + &
          c2(j)*C_n_a(:,:,i) + D_n_a(:,:,i) 
      end do
      vi_a(:,:,:,j) = -vi_a(:,:,:,j) / g(j) 
    end do
  end select
  !group velocity wrt a
  rv_a = 0.0_RP
  do concurrent( j = 1: 3) 
    do i = 1, 3 
      rv_a(:,:,j) = rv_a(:,:,j) + vi(i,j)*vi_a(:,:,i,j)
    end do
    rv_a(:,:,j) = rv_a(:,:,j) / rv(j)
  end do
  end subroutine

  pure subroutine phaseAndRayVelocity_mkl(a, n, pv, rv, vi)
  use lapack95, only: syev, geev
  implicit none
  !//**********************************************************************
  !// eigenvector method
  !//**********************************************************************
  DATATYPE, intent(in):: a(6,6), n(3)
  DATATYPE, intent(out):: pv(3), rv(3), vi(3,3)
  DATATYPE cstl(6), g(3,3), temp(3,3)
  integer i
  !christoffel
  cstl = christoffelMatrix(a, n)
  g(:,1) = cstl(1:3)
  g(2:3,2) = cstl(4:5)
  g(3,3) = cstl(6)
#ifdef COMPLEX
  g(1,2) = g(2,1); g(1:2,3) = g(3,1:2)
  call geev(g,w=pv,vr=temp)
  g = temp
  i = maxloc(abs(pv),dim=1)
  if(i/=3) then
    cstl(1) = pv(3); pv(3) = pv(i); pv(i) = cstl(1)
    cstl(1:3) = g(:,3); g(:,3) = g(:,i); g(:,i) = cstl(1:3)
  end if
  i = maxloc(abs(pv(1:2)),dim=1)
  if(i/=2) then
    cstl(1) = pv(2); pv(2) = pv(i); pv(i) = cstl(1)
    cstl(1:3) = g(:,2); g(:,2) = g(:,i); g(:,i) = cstl(1:3)
  end if
  do concurrent (i=1:3)
    g(:,i) = g(:,i) / sqrt(sum(g(:,i)*g(:,i)))
  end do
#else
  call syev(g, pv, "V", "L")
#endif
  pv = sqrt(pv)
  !group velocity
  do i=1,3
    cstl = christoffelMatrix(a, g(:,i))
    vi(1,i) = cstl(1)*n(1) + cstl(2)*n(2) + cstl(3)*n(3)
    vi(2,i) = cstl(2)*n(1) + cstl(4)*n(2) + cstl(5)*n(3)
    vi(3,i) = cstl(3)*n(1) + cstl(5)*n(2) + cstl(6)*n(3)
    vi(:,i) = vi(:,i) / pv(i)
    rv(i) = sqrt(vi(1,i)*vi(1,i)+vi(2,i)*vi(2,i)+vi(3,i)*vi(3,i))
  end do
  end subroutine



  end module

