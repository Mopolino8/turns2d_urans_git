      real sigmak1,sigmaw1,beta1,betastar,kappa,gamma1
      real sigmak2,sigmaw2,beta2,gamma2
      data sigmak1,sigmaw1,beta1,sigmak2,sigmaw2,beta2,betastar,kappa /0.85,0.5,0.075,1.0,0.856,0.0828,0.09,0.41/
!      data sigmak1,sigmaw1,beta1,sigmak2,sigmaw2,beta2,betastar,kappa /1.176,2.0,0.075,1.0,0.856,0.0828,0.09,0.41/
!      gamma1 = beta1/betastar - sigmaw1*kappa**2/sqrt(betastar)
!      gamma2 = beta2/betastar - sigmaw2*kappa**2/sqrt(betastar)
      gamma1 = 5./9
      gamma2 = 0.44
