# Velocity4GeneralAnisotropicMedium
author: LI XingWang, Chang'an University

e-mail: lixingwang@chd.edu.cn

reference: Wang D, Bai C Y, Li X W, Hu J. 2023. Analytic Formula for the Group Velocity and Its Derivatives with Respect to Elastic Moduli for a General Anisotropic Medium. Geophysics, 88(4): C111-C121. DOI: 10.1190/geo2022-0566.1

This code is used to calculate the group velocity and its derivatives with respect to 21 elastic moduli for a general anisotropic medium
1. floating point variables are double-precision, 
2. and you should compile the codes with option -DDOUBLE.
3. lapack or mkl is needed by subroutine phaseAndRayVelocity_mkl, but this sub is not necessary.
