!------------------------------------------------------------------------------
module SF6TracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  SF6TracersMetadata_mod encapsulates the SF6 tracers metadata
!@auth NCCS ASTG
!  use sharedTracersMetadata_mod
  use TRACER_COM, only: n_SF6, n_SF6_c
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_tr_mm, set_ntm_power
  use Tracer_mod, only: Tracer
  use OldTracer_mod, only: set_has_chemistry
  implicit none

  private
  public SF6_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine SF6_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  SF6_setSpec('SF6')
    call  SF6_c_setSpec('SF6_c')

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine SF6_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_has_chemistry(n, .true.)
    end subroutine SF6_setSpec

    subroutine SF6_c_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6_c = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_has_chemistry(n, .true.)
    end subroutine SF6_c_setSpec

  end subroutine SF6_InitMetadata

end module SF6TracersMetadata_mod
