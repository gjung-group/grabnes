module kuboarrays

   use mio

   implicit none

   real(dp), public, pointer :: a(:)=>NULL(),b(:)=>NULL()
   complex(dp), public, pointer :: Psi(:)=>NULL(), ZUPsi(:)=>NULL(), c(:)=>NULL()
   complex(dp), public, pointer :: Psin(:)=>NULL(), Psinm1(:)=>NULL(), XpnPsi(:)=>NULL(), XpnPsim1(:)=>NULL()
   complex(dp), public, pointer :: H(:)=>NULL()
   complex(dp), public, pointer :: tempZ(:)=>NULL()
   real(dp), public, pointer :: tempD(:)=>NULL()

end module kuboarrays
