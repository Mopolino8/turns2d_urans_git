C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.6 (r4343) - 10 Feb 2012 10:52
C
C  Differentiation of calc_cp in reverse (adjoint) mode:
C   gradient     of useful results: fmtip fsmach cp
C   with respect to varying inputs: fmtip fsmach q
C
C
      SUBROUTINE CALC_CP_BO(jd, kd, q, qb, cp, cpb)
      USE PARAMS_GLOBAL
      USE PARAMS_SENSITIVITY
      IMPLICIT NONE
      INTEGER jd, kd
      REAL q(jd, kd, nq), cp(jd)
      REAL qb(jd, kd, nq), cpb(jd)
      INTEGER j
      INTEGER j1, j2, k
      REAL fmm, rr, rho, u, v, e, vsq, pp, pp1, cp_
      REAL fmmb, rrb, rhob, ub, vb, eb, vsqb, ppb, pp1b, cp_b
      INTEGER branch
      REAL temp0
      REAL temp0b
      REAL tempb0
      REAL tempb
      REAL temp0b0
      REAL temp
      k = kbeg
C     set normalization factor
      IF (fsmach .NE. 0) THEN
        fmm = fsmach
        CALL PUSHCONTROL2B(2)
      ELSE IF (fmtip .NE. 0) THEN
        fmm = fmtip
        CALL PUSHCONTROL2B(1)
      ELSE
        CALL PUSHCONTROL2B(0)
        fmm = 1.0
      END IF
      qb = 0.0
      fmmb = 0.0
      DO j=jtail2,jtail1,-1
        cp_b = 0.5*cpb(j)
        cpb(j) = 0.5*cpb(j)
        e = q(j, k-1, 4)*q(j, k-1, nq)
        rr = 1./q(j, k-1, 1)
        u = q(j, k-1, 2)*rr
        v = q(j, k-1, 3)*rr
        vsq = u*u + v*v
        rho = q(j, k-1, 1)*q(j, k-1, nq)
        pp = gm1*(e-rho*vsq/2.)
        pp1 = pp/pinf
        temp0 = gamma*fmm**2
        temp0b = 2.*cp_b/temp0
        pp1b = temp0b
        fmmb = fmmb - gamma*(pp1-1.)*2*fmm*temp0b/temp0
        ppb = pp1b/pinf
        temp0b0 = gm1*ppb
        eb = temp0b0
        rhob = -(vsq*temp0b0/2.)
        vsqb = -(rho*temp0b0/2.)
        ub = 2*u*vsqb
        vb = 2*v*vsqb
        qb(j, k-1, 4) = qb(j, k-1, 4) + q(j, k-1, nq)*eb
        qb(j, k-1, nq) = qb(j, k-1, nq) + q(j, k-1, 4)*eb
        qb(j, k-1, 3) = qb(j, k-1, 3) + rr*vb
        rrb = q(j, k-1, 2)*ub + q(j, k-1, 3)*vb
        qb(j, k-1, 2) = qb(j, k-1, 2) + rr*ub
        qb(j, k-1, 1) = qb(j, k-1, 1) + q(j, k-1, nq)*rhob
        qb(j, k-1, nq) = qb(j, k-1, nq) + q(j, k-1, 1)*rhob
        qb(j, k-1, 1) = qb(j, k-1, 1) - rrb/q(j, k-1, 1)**2
        cp_b = cpb(j)
        cpb(j) = 0.0
        e = q(j, k, 4)*q(j, k, nq)
        rr = 1./q(j, k, 1)
        u = q(j, k, 2)*rr
        v = q(j, k, 3)*rr
        vsq = u*u + v*v
        rho = q(j, k, 1)*q(j, k, nq)
        pp = gm1*(e-rho*vsq/2.)
        pp1 = pp/pinf
        temp = gamma*fmm**2
        tempb = 2.*cp_b/temp
        pp1b = tempb
        fmmb = fmmb - gamma*(pp1-1.)*2*fmm*tempb/temp
        ppb = pp1b/pinf
        tempb0 = gm1*ppb
        eb = tempb0
        rhob = -(vsq*tempb0/2.)
        vsqb = -(rho*tempb0/2.)
        ub = 2*u*vsqb
        vb = 2*v*vsqb
        qb(j, k, 4) = qb(j, k, 4) + q(j, k, nq)*eb
        qb(j, k, nq) = qb(j, k, nq) + q(j, k, 4)*eb
        qb(j, k, 3) = qb(j, k, 3) + rr*vb
        rrb = q(j, k, 2)*ub + q(j, k, 3)*vb
        qb(j, k, 2) = qb(j, k, 2) + rr*ub
        qb(j, k, 1) = qb(j, k, 1) + q(j, k, nq)*rhob
        qb(j, k, nq) = qb(j, k, nq) + q(j, k, 1)*rhob
        qb(j, k, 1) = qb(j, k, 1) - rrb/q(j, k, 1)**2
      ENDDO
      CALL POPCONTROL2B(branch)
      IF (branch .NE. 0) THEN
        IF (branch .EQ. 1) THEN
          fmtipb = fmtipb + fmmb
        ELSE
          fsmachb = fsmachb + fmmb
        END IF
      END IF
      END