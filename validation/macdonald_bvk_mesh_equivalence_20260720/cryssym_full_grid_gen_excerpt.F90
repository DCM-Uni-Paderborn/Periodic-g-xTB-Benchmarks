! Exact excerpt of the qualified CP2K cryssym.F full_grid_gen coordinate rule.
! Qualified CP2K source revision:
! 8520b2e592cd04d35081ab4ad46d92c606071e23
DO idim = 1, 3
   IF (gamma_mesh .AND. MOD(nk(idim), 2) == 0) THEN
      kpt_latt(idim) = REAL(2*ik(idim) - nk(idim), KIND=dp)/ &
                       (2._dp*REAL(nk(idim), KIND=dp))
   ELSE
      kpt_latt(idim) = REAL(2*ik(idim) - nk(idim) - 1, KIND=dp)/ &
                       (2._dp*REAL(nk(idim), KIND=dp))
   END IF
END DO
xkp(1:3, i) = kpt_latt(1:3)
DO i = 1, nk(1)*nk(2)*nk(3)
   xkp(1:3, i) = xkp(1:3, i) + shift(1:3)
END DO
