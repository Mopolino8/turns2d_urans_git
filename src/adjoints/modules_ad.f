      module params_adjoints

        real cfx_totb,cfy_totb,cl_totb,cd_totb,cm_totb,cpower_totb

        data cfx_totb,cfy_totb,cl_totb,cd_totb/0.0,0.0,0.0,0.0/
        data cm_totb,cpower_totb/0.0,0.0/

        namelist /adinputs/ cfx_totb,cfy_totb,cl_totb,cd_totb,cm_totb,
     &                      cpower_totb

      end module params_adjoints

      module params_sensitivity
        real :: uinfb,vinfb,einfb,alfab,fsmachb,fmtipb,reyb
      end module params_sensitivity
