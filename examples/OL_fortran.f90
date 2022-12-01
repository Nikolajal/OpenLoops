

program main
  use openloops
  implicit none
  integer :: id, error, k
  real(8) :: m2_tree, m2_loop(0:2), acc
  real(8) :: p_ex(0:3,4)
  real(8) :: mu = 100, alpha_s = 0.1, energy=1000


  call set_parameter("verbose", 1)
  call set_parameter("no_splash", 1)

  call set_parameter("yuk(3)", 0.)
  call set_parameter("yuk(4)", 1.)
  call set_parameter("yuk(5)", 0.)

  call set_parameter("order_ew", 1)

  id = register_process("-4 4 -> 21 25", 11)

  ! start
  call start()

  if (id > 0) then
    ! generate a random phase space point with rambo
    call phase_space_point(id, energy, p_ex)

    print *
    do k = 1, size(p_ex,2)
      print *, 'P[', int(k,1), '] =', p_ex(:,k)
    end do

    ! evaluate tree matrix element
    call evaluate_tree(id, p_ex, m2_tree)
    ! print tree result
    print *
    print *, "evaluate_tree"
    print *, "Tree:       ", m2_tree

    ! evaluate tree and loop matrix elements
    call evaluate_loop(id, p_ex, m2_tree, m2_loop(0:2), acc)
    ! print loop result
    print *
    print *, "evaluate_loop"
    print *, "Tree:       ", m2_tree
    print *, "Loop ep^0:  ", m2_loop(0)
    print *, "Loop ep^-1: ", m2_loop(1)
    print *, "Loop ep^-2: ", m2_loop(2)
    print *, "accuracy:   ", acc
    print *
  end if

  ! finish
  call finish()

end program main
