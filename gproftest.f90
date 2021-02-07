program main

  use iso_fortran_env

  implicit none

  integer, parameter :: n = 100000000

  real, allocatable :: matrix(:, :, :)
  real, allocatable :: vector(:, :, :)
  real, allocatable :: res_me(:, :, :)
  real, allocatable :: res_blas(:, :, :)

  integer :: seed(8)

  integer :: nthreads
  character(len=32) :: arg

  integer :: i

  call get_command_argument(1, arg)
  read(arg, *) nthreads
  print *, "using #threads: ", nthreads

  call clock(" ")  ! initialize

  allocate(matrix(4, 4, n))
  allocate(vector(4, 1, n))
  allocate(res_me(4, 1, n))
  allocate(res_blas(4, 1, n))

  call clock("memory allocation")

  seed = 0
  call random_seed(put=seed)
  call random_number(matrix)
  call random_number(vector)

  call clock("random number assignment")

  call openblas_set_num_threads(1)
  call omp_set_num_threads(nthreads)

!!$OMP PARALLEL DO
  do i = 1, n
     call multiply(matrix(:, :, n), vector(:, :, n), res_me(:, :, n))
  end do

  call clock("multiply subroutine")

!!$OMP PARALLEL DO
  do i = 1, n
     CALL SGEMV('N' ,4, 4, 1., matrix(:, :, n), 4, vector(:, :, n), 1, 0.0, res_blas(:, :, n), 1)
  end do

  call clock("blas sgemv")

  print *, ""
  print *, "difference:"
  print *, sum(res_me - res_blas)

end program main


subroutine clock(msg)
  use iso_fortran_env
  implicit none

  character(len=*) :: msg
  real, save :: tstart, tstop
  real(real64), save :: ostart, ostop
  real(real64) :: omp_get_wtime
  logical, save :: initialized = .false.

  if (initialized) then
     call cpu_time(tstop)
     ostop = omp_get_wtime()
     print *, msg
     print "(2(F12.8,1X))", tstop - tstart, ostop - ostart
  end if

  initialized = .true.

  call cpu_time(tstart)
  ostart = omp_get_wtime()

end subroutine clock


subroutine multiply(matrix, vector, new_vector)

  implicit none

  real, intent(in) :: matrix(4, 4)
  real, intent(in) :: vector(4, 1)
  real, intent(inout) :: new_vector(4, 1)
  integer :: i, j

  new_vector = 0

  do j = 1,4
     do i = 1,4
        new_vector(i, 1) = new_vector(i, 1) + matrix(i, j) * vector(j, 1)
     end do
  end do

end subroutine multiply
