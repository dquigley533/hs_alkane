! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                          Q U A T E R N I O N                                !
!=============================================================================!
!                                                                             !
! $Id: quaternion.F90,v 1.1 2011/02/02 11:48:36 phseal Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! Routines to compute, manipulate and apply quaternion rotations.             !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: quaternion.F90,v $
! Revision 1.1  2011/02/02 11:48:36  phseal
! Initial revision
!
!
!=============================================================================!
module quaternion

  use constants, only : dp
  implicit none

  private                                       ! Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: quat_get_minimum_arc_q
  public :: quat_axis_angle_to_quat
  public :: quat_product
  public :: quat_inverse
  public :: quat_conjugate_q_with_v


  contains

    subroutine quat_get_minimum_arc_q(v1,v2,quat)
      !-------------------------------------------------------------------------!
      ! Returns the normalised quaternion required to rotate v1 onto v2         !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none

      real(kind=dp),dimension(3),intent(in)  :: v1,v2
      real(kind=dp),dimension(4),intent(out) :: quat

      real(kind=dp) :: m2v1,m2v2

      m2v1   = dot_product(v1,v1)
      m2v2   = dot_product(v2,v2)

      quat(1)   = sqrt(m2v2*m2v1) + dot_product(v1,v2)
      quat(2:4) = cross_product(v1,v2)
      quat      = quat/sqrt(dot_product(quat,quat))

      return

    end subroutine quat_get_minimum_arc_q

    subroutine quat_axis_angle_to_quat(axis,angle,quat)
      !-------------------------------------------------------------------------!
      ! Returns the normalised quaternion which represents rotation by angle    !
      ! about UNIT vector axis. Angle must be in radians.                       !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      
      real(kind=dp),dimension(3),intent(in)  :: axis
      real(kind=dp)             ,intent(in)  :: angle
      real(kind=dp),dimension(4),intent(out) :: quat
      real(kind=dp) :: sh,ch

#ifdef DEBUG
      if  ( sqrt(dot_product(axis,axis))> 1.0_dp ) then
         write(0,'("Warning in quat_axis_angle_to_quat : axis is not a unit vector")')
      end if
#endif 
      
      ch = cos(0.5_dp*angle)
      sh = sin(0.5_dp*angle)

      quat(1)   = ch
      quat(2:4) = axis(:)*sh

      return

    end subroutine quat_axis_angle_to_quat

    function quat_product(a,b,normalise)
      !-------------------------------------------------------------------------!
      ! Returns the quaternion product ab i.e. the quaternion representation of !
      ! the rotation b followed by a.                                           !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      
      real(kind=dp),dimension(4),intent(in) :: a,b
      logical,intent(in),optional           :: normalise

      real(kind=dp),dimension(4) :: quat_product
      real(kind=dp),dimension(4) :: tmpquat
      real(kind=dp)              :: tnorm

      tmpquat(1)   = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4) 
      tmpquat(2:4) = a(1)*b(2:4) + b(1)*a(2:4) + cross_product(a(2:4),b(2:4)) 

      if (present(normalise).and.normalise) then
         tnorm        = sqrt(dot_product(tmpquat,tmpquat))
         quat_product = tmpquat/tnorm
      else
         quat_product = tmpquat
      end if
      
    end function quat_product

    function quat_inverse(a)
      !-------------------------------------------------------------------------!
      ! Returns the inverse quaternion to a. Assumes a is normalised            !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      real(kind=dp),dimension(4) :: quat_inverse
      real(kind=dp),dimension(4),intent(in) :: a
      real(kind=dp),dimension(4) :: tmpquat
      
#ifdef DEBUG
      if  ( sqrt(dot_product(a,a))> 1.0_dp ) then
         write(0,'("Warning in quat_inverse : a is not normalised")')
      end if
#endif 

      tmpquat      = -a
      tmpquat(1)   = -tmpquat(1)
      quat_inverse = tmpquat   
      
    end function quat_inverse

    function quat_conjugate_q_with_v(a,v)
      !-------------------------------------------------------------------------!
      ! Conjugates the quaternion a with the vector v, i.e. performs the        !
      ! rotation represented by the quaternion, which is assumed to be          !
      ! appropriately normalised.                                               !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      real(kind=dp),dimension(3) :: quat_conjugate_q_with_v
      real(kind=dp),dimension(4),intent(in) :: a
      real(kind=dp),dimension(3),intent(in) :: v
      real(kind=dp),dimension(4) :: q,b
      
#ifdef DEBUG
      if  ( sqrt(dot_product(a,a))> 1.0_dp ) then
         write(0,'("Warning in quat_inverse : a is not normalised")')
      end if
#endif 

      q(1)   = 0.0_dp
      q(2:4) = v 
      
      b = quat_inverse(a)
      q = quat_product(q,b)
      q = quat_product(a,q)
      
      quat_conjugate_q_with_v = q(2:4)
                   
    end function quat_conjugate_q_with_v

    function cross_product(a,b)
      !-------------------------------------------------------------------------!
      ! Does exactly what is says on the tin.                                   !
      !-------------------------------------------------------------------------!
      ! D.Quigley January 2010                                                  !
      !-------------------------------------------------------------------------!
      implicit none
      real(kind=dp),dimension(3),intent(in) :: a,b
      real(kind=dp),dimension(3) :: cross_product

      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = -(a(1)*b(3)-a(3)*b(1))
      cross_product(3) = a(1)*b(2)-b(1)*a(2)

    end function cross_product


end module quaternion
