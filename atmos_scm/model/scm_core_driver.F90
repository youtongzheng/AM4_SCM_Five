
module scm_core_driver_mod

!-----------------------------------------------------------------------
!
!      wrapper module for initializing the b-grid dynamics core
!
!         reads namelist
!         initializes bgrid_core_mod
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------

use fms_mod  , only : error_mesg, FATAL
use    time_manager_mod, only:  time_type, operator(+), get_time
use scm_prog_var_mod, only:   prog_var_type

use  constants_mod, only: rdgas, rvgas, grav
use   scm_forc_mod, only: update_scm_forc, scm_forc_end, scm_data_read,&
                          scm_forc_init
use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index
!-----------------------------------------------------------------------

implicit none
private

public  scm_core_driver_init, &
        scm_core_driver,      &
        scm_core_driver_end

!-----------------------------------------------------------------------
character(len=128) :: version =  '$Id$'
character(len=128) :: tag =  '$Name$'
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------

subroutine scm_core_driver_init (Time_init, pdamp,elev,As, Bs &
                                 )
#include "fv_arrays.h"

TYPE(TIME_TYPE)                          :: Time_init
REAL,  INTENT (IN)                       :: pdamp
REAL,  INTENT (IN)   , DIMENSION(:)      :: As,Bs
REAL,  INTENT (INOUT), DIMENSION(:,:)    :: elev

integer :: unit, clat,clon,kdim, k, kd

!-----------------------------------------------------------------------

#include "fv_point.inc"
           call scm_forc_init(Time_init,pdamp,elev,As,Bs)

 end subroutine scm_core_driver_init

!#######################################################################

 subroutine scm_core_driver (Time, Time_step_atmos, pdamp, elev, Forc_tend, do_five) !yzheng
#include "fv_arrays.h"

!-----------------------------------------------------------------------
   type (time_type),      intent(in)    :: Time, Time_step_atmos
   real,                  intent(in)    :: pdamp
   real, dimension(:,:),  intent(in)    :: elev
   type (prog_var_type),  intent(inout) :: Forc_tend
   logical,                      intent(in)  :: do_five !yzheng

   type (time_type) :: Time_diag  
#include "fv_point.inc"

   Time_diag = Time + Time_step_atmos

   call update_scm_forc (Time, Time_diag, Time_step_atmos, pdamp, elev, do_five)

   !---assign forcing tendencies
   Forc_tend%T = t_dt
   Forc_tend%U = u_dt
   Forc_tend%V = v_dt
   Forc_tend%R = q_dt(:,:,:,1:Forc_tend%ntrace)


 end subroutine scm_core_driver

!#######################################################################

 subroutine scm_core_driver_end ()



   call scm_forc_end ()

!-----------------------------------------------------------------------

 end subroutine scm_core_driver_end

!#######################################################################

end module scm_core_driver_mod

