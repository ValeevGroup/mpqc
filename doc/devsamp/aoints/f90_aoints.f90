program f90_aoints
    implicit none

    ! interfaces need to be described for certain functions
    interface
     function basis_set_nshell ()
       integer :: basis_set_nshell
     end function basis_set_nshell
     function basis_set_nbasis ()
       integer :: basis_set_nbasis
     end function basis_set_nbasis
     function basis_set_max_nfunction_in_shell ()
       integer :: basis_set_max_nfunction_in_shell
     end function basis_set_max_nfunction_in_shell
    end interface

    ! not needed if allocating memory dynamically
    integer, parameter :: max_nshell = 1000
    ! not needed if allocating memory dynamically
    integer, parameter :: max_nfunction_per_shell = 30

    integer, parameter :: natoms=3
    integer, parameter :: use_symmetry=1
    real*8 :: Z(natoms), xyz(3,natoms)

    integer :: nshell, nbasis
    integer :: shell_to_nfunction(0:max_nshell-1), shell_to_function(0:max_nshell-1)
    integer :: s
    integer :: s1, s2, s3, s4
    integer :: bf1, bf2, bf3, bf4, bf12, bf1234
    integer :: bf1_offset, bf2_offset, bf3_offset, bf4_offset
    integer :: nbf1, nbf2, nbf3, nbf4

    real*8 :: oneel_ints_buffer(0:max_nfunction_per_shell**2-1)
    real*8 :: twoel_ints_buffer(0:max_nfunction_per_shell**4-1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! define molecular geometry
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! create molecule
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_molecule(%val(natoms), Z, xyz, %val(use_symmetry))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! create basis set
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_basis_set('cc-pVDZ')
    nshell = basis_set_nshell()
    nbasis = basis_set_nbasis()
    call basis_set_shell_to_nfunction_subrt(shell_to_nfunction)
    call basis_set_shell_to_function_subrt(shell_to_function)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! prepare integrals factory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_integrals()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! computing overlap integrals
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, 'overlap integrals:'
    call init_overlap_integrals()
    ! overlap integral shells will be returned in oneel_ints_buffer
    call set_overlap_integrals_buffer(oneel_ints_buffer)
    do s1=0,nshell-1
      bf1_offset = shell_to_function(s1)
      nbf1 = shell_to_nfunction(s1)

      do s2=0,nshell-1
        bf2_offset = shell_to_function(s2)
        nbf2 = shell_to_nfunction(s2)

        call compute_overlap_shell(%val(s1),%val(s2))

        bf12 = 0
        do bf1=0,nbf1-1
          do bf2=0,nbf2-1
            print *, bf1+bf1_offset, bf2+bf2_offset, oneel_ints_buffer(bf12)
            bf12 = bf12 + 1
          end do
        end do
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! computing hcore integrals
    ! (kinetic + nuclear attraction)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, 'hcore integrals:'
    call init_hcore_integrals()
    call set_hcore_integrals_buffer(oneel_ints_buffer)
    do s1=0,nshell-1
      bf1_offset = shell_to_function(s1)
      nbf1 = shell_to_nfunction(s1)

      do s2=0,nshell-1
        bf2_offset = shell_to_function(s2)
        nbf2 = shell_to_nfunction(s2)

        call compute_hcore_shell(%val(s1),%val(s2))

        bf12 = 0
        do bf1=0,nbf1-1
          do bf2=0,nbf2-1
            print *, bf1+bf1_offset, bf2+bf2_offset, oneel_ints_buffer(bf12)
            bf12 = bf12 + 1
          end do
        end do
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! computing two-electron Coulomb integrals
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, 'two-e Coulomb integrals:'
    call init_twoecoulomb_integrals()
    call set_twoecoulomb_integrals_buffer(twoel_ints_buffer)
    do s1=0,nshell-1
      bf1_offset = shell_to_function(s1)
      nbf1 = shell_to_nfunction(s1)

      do s2=0,nshell-1
        bf2_offset = shell_to_function(s2)
        nbf2 = shell_to_nfunction(s2)

        do s3=0,nshell-1
          bf3_offset = shell_to_function(s3)
          nbf3 = shell_to_nfunction(s3)

          do s4=0,nshell-1
            bf4_offset = shell_to_function(s4)
            nbf4 = shell_to_nfunction(s4)

            call compute_twoecoulomb_shell(%val(s1),%val(s2),%val(s3),%val(s4))

            bf1234 = 0
            do bf1=0,nbf1-1
              do bf2=0,nbf2-1
                do bf3=0,nbf3-1
                 do bf4=0,nbf4-1
                   print *, bf1+bf1_offset, bf2+bf2_offset, &
                            bf3+bf3_offset, bf4+bf4_offset, &
                            twoel_ints_buffer(bf1234)
                   bf1234 = bf1234 + 1
                 end do
               end do
              end do
            end do

          end do
        end do
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! clean-up
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call done_integrals()

end program f90_aoints
