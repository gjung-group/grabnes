module constants

   use mio

   implicit none

   PRIVATE

   real(dp), parameter, public :: pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640_dp
   real(dp), parameter, public :: twopi = 2.0_dp*pi
   complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
   complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
   complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)

   real(dp), parameter, public :: hPlanck = 6.62606896d-34, hbar =  hPlanck/twopi
   real(dp), parameter, public :: qe = 1.602176565d-19, fluxq = hPlanck/qe
   real(dp), parameter, public :: uB=5.7883818066d-5

end module constants
