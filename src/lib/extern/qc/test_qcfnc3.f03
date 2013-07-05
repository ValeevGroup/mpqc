program test_qcfnc3
    use iso_c_binding

    implicit none

    integer, parameter :: natoms=3
    integer, parameter :: use_symmetry=1
    real*8 :: Z(natoms), xyz(3,natoms)
    integer :: nshell, nbasis
    type(c_ptr) :: cptr
    integer, pointer :: shell_to_nfunction(:), shell_to_function(:)
    integer :: ishape(1)
    integer :: s

    interface
     subroutine init_molecule (natoms, Z, xyz, use_symmetry)
       integer :: natoms, use_symmetry
       real*8 :: Z(:), xyz(:,:)
     end subroutine init_molecule
     function basis_set_nshell ()
       integer :: basis_set_nshell
     end function basis_set_nshell
     function basis_set_nbasis ()
       integer :: basis_set_nbasis
     end function basis_set_nbasis
     function basis_set_shell_to_nfunction ()
       import :: c_ptr
       type(c_ptr) :: basis_set_shell_to_nfunction
     end function basis_set_shell_to_nfunction
    end interface

    Z(1) = 8.0
    xyz(1,1) = 0.0
    xyz(2,1) = 0.0
    xyz(3,1) = 0.0
    Z(2) = 1.0
    xyz(1,2) = 0.0
    xyz(2,2) = 1.0
    xyz(3,2) = 1.0
    Z(3) = 1.0
    xyz(1,3) = 0.0
    xyz(2,3) =-1.0
    xyz(3,3) = 1.0

    call init_molecule_(natoms, Z, xyz, use_symmetry)

    call init_basis_set_('cc-pVDZ')
    nshell = basis_set_nshell_()
    nbasis = basis_set_nbasis_()
    write (*,*) 'nshell=', nshell, ' nbasis=', nbasis

    cptr = basis_set_shell_to_nfunction()
    ishape(1) = nshell
    call c_f_pointer(cptr, shell_to_nfunction, ishape)
    do s=1,nshell
      print *, shell_to_nfunction(s)
    end do

end program test_qcfnc3
