!> Functions to compute did beta distribution quantities 
!> required additional functions are included so no outside code is needed
!> Currently working beta distributions functions
!> betacdf => compute CDFs for beta distribution using incomplete beta function
!> betapdf => compute PDFs for beta distrubution
!> betainv => compute inverse CDFs for beta distribution


module truncnormal_dist_functions

use types_mod, only: digits12, i8, r8, PI
use utilities_mod,       only : error_handler, E_ERR, E_MSG
use 
implicit none
private


character(len=512)     :: msgstring
character(len=*), parameter :: source = 'truncnormal_dist_functions.f90'

contains








end module truncnormal_dist_functions