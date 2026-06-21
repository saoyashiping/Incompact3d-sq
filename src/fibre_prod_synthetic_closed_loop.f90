module fibre_prod_synthetic_closed_loop
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fibre_prod_grid_adapter, only : fibre_prod_grid_type, fibre_prod_grid_init_from_coordinates, fibre_prod_grid_destroy
  use fibre_prod_ibm_force_buffer, only : fibre_prod_force_buffer_type, fibre_prod_force_buffer_allocate, &
                                          fibre_prod_force_buffer_destroy
  use fibre_prod_runtime_config, only : fibre_prod_runtime_config_type, fibre_prod_runtime_config_default
  use fibre_prod_state, only : fibre_prod_state_type, fibre_prod_state_allocate, fibre_prod_state_destroy, &
                               fibre_prod_state_attach_structure_u
  use fibre_prod_state_velocity_attachment, only : fibre_prod_state_velocity_attachment_type, &
                                                   fibre_prod_state_velocity_attachment_init, &
                                                   fibre_prod_state_velocity_attachment_set_points, &
                                                   fibre_prod_state_velocity_attachment_sample, &
                                                   fibre_prod_state_velocity_attachment_attach_to_state, &
                                                   fibre_prod_state_velocity_attachment_finalize
  use fibre_prod_hydro_input_candidate, only : fibre_prod_hydro_input_candidate_type, &
                                               fibre_prod_hydro_input_candidate_init, &
                                               fibre_prod_hydro_input_candidate_compute, &
                                               fibre_prod_hydro_input_candidate_attach_to_state, &
                                               fibre_prod_hydro_input_candidate_finalize
  use fibre_prod_structure_input_handoff, only : fibre_prod_structure_input_handoff_type, &
                                                 fibre_prod_structure_input_handoff_init, &
                                                 fibre_prod_structure_input_handoff_from_candidate, &
                                                 fibre_prod_structure_input_handoff_attach_to_state, &
                                                 fibre_prod_structure_input_handoff_finalize
  use fibre_prod_structure_dry_step, only : fibre_prod_structure_dry_step_type, &
                                            fibre_prod_structure_dry_step_init, &
                                            fibre_prod_structure_dry_step_predict, &
                                            fibre_prod_structure_dry_step_check_bounded, &
                                            fibre_prod_structure_dry_step_finalize
  use fibre_prod_structure_commit_gate, only : fibre_prod_structure_commit_gate_type, &
                                               fibre_prod_structure_commit_gate_init, &
                                               fibre_prod_structure_commit_gate_set_enabled, &
                                               fibre_prod_structure_commit_gate_evaluate, &
                                               fibre_prod_structure_commit_gate_commit_to_state, &
                                               fibre_prod_structure_commit_gate_finalize
  use fibre_prod_reaction_force_candidate, only : fibre_prod_reaction_force_candidate_type, &
                                                  fibre_prod_reaction_force_candidate_init, &
                                                  fibre_prod_reaction_force_candidate_from_structure_input, &
                                                  fibre_prod_reaction_force_candidate_finalize
  use fibre_prod_reaction_spreading_buffer, only : fibre_prod_reaction_spreading_buffer_type, &
                                                   fibre_prod_reaction_spreading_buffer_init, &
                                                   fibre_prod_reaction_spreading_buffer_apply, &
                                                   fibre_prod_reaction_spreading_buffer_finalize
  use fibre_prod_force_buffer_rhs_gate, only : fibre_prod_force_buffer_rhs_gate_type, &
                                               fibre_prod_force_buffer_rhs_gate_init, &
                                               fibre_prod_force_buffer_rhs_gate_apply, &
                                               fibre_prod_force_buffer_rhs_gate_finalize
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: fibre_prod_synthetic_closed_loop_type
  public :: fibre_prod_synthetic_closed_loop_init
  public :: fibre_prod_synthetic_closed_loop_run
  public :: fibre_prod_synthetic_closed_loop_get_signature
  public :: fibre_prod_synthetic_closed_loop_finalize
  public :: fibre_prod_synthetic_closed_loop_env_enabled
  public :: fibre_prod_synthetic_closed_loop_diagnostics_env_enabled
  public :: fibre_prod_synthetic_closed_loop_runtime_diagnostic

  type :: fibre_prod_synthetic_closed_loop_type
    logical :: initialized = .false.
    logical :: completed = .false.
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    integer :: nnode = 0
    real(dp) :: lambda_fsi = 0.0_dp
    real(dp) :: penalty_beta = 1.0_dp
    real(dp) :: beta_hydro = 1.0_dp
    real(dp) :: dt = 1.0e-4_dp
    real(dp) :: rho_eff = 1.0_dp
    real(dp) :: max_allowed_dx = 1.0e-2_dp
    real(dp) :: max_abs_sampled_u = 0.0_dp
    real(dp) :: max_abs_hydro_candidate = 0.0_dp
    real(dp) :: max_abs_structure_input = 0.0_dp
    real(dp) :: max_abs_dx_trial = 0.0_dp
    real(dp) :: max_abs_reaction_force = 0.0_dp
    real(dp) :: max_abs_force_buffer = 0.0_dp
    real(dp) :: max_abs_rhs_increment = 0.0_dp
    real(dp) :: sum_abs_rhs_increment = 0.0_dp
    real(dp) :: net_lagrangian_reaction_force(3) = 0.0_dp
    real(dp) :: net_eulerian_force_buffer(3) = 0.0_dp
    real(dp) :: signature(16) = 0.0_dp
  end type fibre_prod_synthetic_closed_loop_type

contains

  subroutine fibre_prod_synthetic_closed_loop_init(loop, nx, ny, nz, nnode, status)
    type(fibre_prod_synthetic_closed_loop_type), intent(inout) :: loop
    integer, intent(in) :: nx, ny, nz, nnode
    integer, intent(out) :: status

    status = 0
    call fibre_prod_synthetic_closed_loop_finalize(loop)
    if (nx <= 1 .or. ny <= 1 .or. nz <= 1 .or. nnode <= 0) then
      status = 1
      return
    end if
    loop%initialized = .true.
    loop%nx = nx
    loop%ny = ny
    loop%nz = nz
    loop%nnode = nnode
  end subroutine fibre_prod_synthetic_closed_loop_init

  subroutine fibre_prod_synthetic_closed_loop_run(loop, lambda_fsi, penalty_beta, beta_hydro, dt, rho_eff, status)
    type(fibre_prod_synthetic_closed_loop_type), intent(inout) :: loop
    real(dp), intent(in) :: lambda_fsi, penalty_beta, beta_hydro, dt, rho_eff
    integer, intent(out) :: status
    type(fibre_prod_grid_type) :: grid
    type(fibre_prod_state_type) :: state
    type(fibre_prod_state_velocity_attachment_type) :: attach
    type(fibre_prod_hydro_input_candidate_type) :: hydro
    type(fibre_prod_structure_input_handoff_type) :: handoff
    type(fibre_prod_structure_dry_step_type) :: dry
    type(fibre_prod_structure_commit_gate_type) :: commit_gate
    type(fibre_prod_reaction_force_candidate_type) :: reaction
    type(fibre_prod_reaction_spreading_buffer_type) :: spreading
    type(fibre_prod_force_buffer_type) :: force_buffer
    type(fibre_prod_force_buffer_rhs_gate_type) :: rhs_gate
    type(fibre_prod_runtime_config_type) :: config
    real(dp), allocatable :: xcoord(:), ycoord(:), zcoord(:)
    real(dp), allocatable :: ux(:, :, :), uy(:, :, :), uz(:, :, :)
    real(dp), allocatable :: rhs_x(:, :, :), rhs_y(:, :, :), rhs_z(:, :, :)
    real(dp), allocatable :: rhs0_x(:, :, :), rhs0_y(:, :, :), rhs0_z(:, :, :)
    real(dp), allocatable :: points(:, :), structure_u(:, :)
    real(dp), allocatable :: rhs_increment_x(:, :, :)
    real(dp) :: scale

    status = 0
    if (.not. loop%initialized) then
      status = 2
      return
    end if
    if (.not. all(ieee_is_finite([lambda_fsi, penalty_beta, beta_hydro, dt, rho_eff])) .or. &
        lambda_fsi < 0.0_dp .or. penalty_beta < 0.0_dp .or. beta_hydro < 0.0_dp .or. &
        dt <= 0.0_dp .or. rho_eff <= 0.0_dp) then
      status = 3
      return
    end if

    allocate(xcoord(loop%nx), ycoord(loop%ny), zcoord(loop%nz), &
             ux(loop%nx,loop%ny,loop%nz), uy(loop%nx,loop%ny,loop%nz), uz(loop%nx,loop%ny,loop%nz), &
             rhs_x(loop%nx,loop%ny,loop%nz), rhs_y(loop%nx,loop%ny,loop%nz), rhs_z(loop%nx,loop%ny,loop%nz), &
             rhs0_x(loop%nx,loop%ny,loop%nz), rhs0_y(loop%nx,loop%ny,loop%nz), rhs0_z(loop%nx,loop%ny,loop%nz), &
             points(loop%nnode,3), structure_u(loop%nnode,3), &
             rhs_increment_x(loop%nx,loop%ny,loop%nz), stat=status)
    if (status /= 0) then
      status = 4
      return
    end if

    call fill_unit_coordinates(xcoord)
    call fill_unit_coordinates(ycoord)
    call fill_unit_coordinates(zcoord)
    call fill_velocity_field(xcoord, ycoord, zcoord, ux, uy, uz)
    call fill_fibre_points(points)
    structure_u = 0.0_dp
    call fill_rhs_seed(rhs0_x, rhs0_y, rhs0_z)
    rhs_x = rhs0_x
    rhs_y = rhs0_y
    rhs_z = rhs0_z

    call fibre_prod_grid_init_from_coordinates(grid, xcoord, ycoord, zcoord, 1, loop%nx, 1, loop%ny, 1, loop%nz, &
                                               .false., .false., .false., status)
    if (status == 0) call fibre_prod_state_allocate(state, 1, loop%nnode, status)
    if (status == 0) then
      state%x(1,:,:) = points
      call fibre_prod_state_attach_structure_u(state, structure_u, status)
    end if
    if (status == 0) call fibre_prod_state_velocity_attachment_init(attach, loop%nnode, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_set_points(attach, points, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_sample(grid, ux, uy, uz, attach, status)
    if (status == 0) call fibre_prod_state_velocity_attachment_attach_to_state(attach, state, status)

    if (status == 0) call fibre_prod_hydro_input_candidate_init(hydro, loop%nnode, beta_hydro, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_compute(hydro, attach%sampled_u, structure_u, status)
    if (status == 0) call fibre_prod_hydro_input_candidate_attach_to_state(hydro, state, status)

    if (status == 0) call fibre_prod_structure_input_handoff_init(handoff, loop%nnode, status)
    if (status == 0) call fibre_prod_structure_input_handoff_from_candidate(handoff, hydro%candidate_force, status)
    if (status == 0) call fibre_prod_structure_input_handoff_attach_to_state(handoff, state, status)

    if (status == 0) call fibre_prod_structure_dry_step_init(dry, loop%nnode, dt, rho_eff, status)
    if (status == 0) call fibre_prod_structure_dry_step_predict(dry, points, structure_u, state%structure_input_force, status)
    if (status == 0) call fibre_prod_structure_dry_step_check_bounded(dry, loop%max_allowed_dx, status)

    if (status == 0) call fibre_prod_structure_commit_gate_init(commit_gate, loop%nnode, loop%max_allowed_dx, status)
    if (status == 0) call fibre_prod_structure_commit_gate_set_enabled(commit_gate, .true., status)
    if (status == 0) call fibre_prod_structure_commit_gate_evaluate(commit_gate, dry%x_trial, dry%u_trial, dry%dx_trial, status)
    if (status == 0) call fibre_prod_structure_commit_gate_commit_to_state(commit_gate, state, dry%x_trial, dry%u_trial, status)

    if (status == 0) call fibre_prod_reaction_force_candidate_init(reaction, loop%nnode, status)
    if (status == 0) call fibre_prod_reaction_force_candidate_from_structure_input(reaction, state%structure_input_force, status)
    if (status == 0) call fibre_prod_force_buffer_allocate(force_buffer, grid, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_init(spreading, loop%nx, loop%ny, loop%nz, loop%nnode, status)
    if (status == 0) call fibre_prod_reaction_spreading_buffer_apply(spreading, grid, points, reaction%reaction_force, &
                                                                      force_buffer, status)

    if (status == 0) then
      call fibre_prod_runtime_config_default(config)
      config%enabled = .true.
      config%lambda_fsi = lambda_fsi
      config%penalty_beta = penalty_beta
      call fibre_prod_force_buffer_rhs_gate_init(rhs_gate, loop%nx, loop%ny, loop%nz, lambda_fsi, penalty_beta, status)
    end if
    if (status == 0) call fibre_prod_force_buffer_rhs_gate_apply(rhs_gate, config, rhs_x, rhs_y, rhs_z, force_buffer, status)

    if (status == 0) then
      rhs_increment_x = rhs_x - rhs0_x
      scale = lambda_fsi * penalty_beta
      if (lambda_fsi == 0.0_dp .and. maxval(abs(rhs_increment_x)) /= 0.0_dp) status = 40
      if (status == 0 .and. maxval(abs(rhs_increment_x - scale * force_buffer%fx)) > 1.0e-10_dp) status = 42
      ! P0_12/P0_13 closure: a nonzero RHS increment is required only when the
      ! Eulerian force-buffer component being checked is nonzero.  The
      ! zero-force proxy case deliberately sets beta_hydro=0, producing a zero
      ! reaction force, zero force_buffer%fx, and therefore zero RHS increment
      ! even when lambda_fsi>0.  That is a valid no-response proof, not a fail.
      if (status == 0 .and. lambda_fsi > 0.0_dp .and. &
          maxval(abs(force_buffer%fx)) > 0.0_dp .and. &
          maxval(abs(rhs_increment_x)) <= 0.0_dp) status = 41
    end if

    if (status == 0) then
      loop%lambda_fsi = lambda_fsi
      loop%penalty_beta = penalty_beta
      loop%beta_hydro = beta_hydro
      loop%dt = dt
      loop%rho_eff = rho_eff
      loop%max_abs_sampled_u = attach%max_abs_sampled_velocity
      loop%max_abs_hydro_candidate = hydro%max_abs_candidate_force
      loop%max_abs_structure_input = handoff%max_abs_input_force
      loop%max_abs_dx_trial = dry%max_abs_dx_trial
      loop%max_abs_reaction_force = reaction%max_abs_reaction_force
      loop%max_abs_force_buffer = spreading%max_abs_force_buffer
      loop%max_abs_rhs_increment = rhs_gate%max_abs_rhs_increment
      loop%sum_abs_rhs_increment = rhs_gate%sum_abs_rhs_increment
      loop%net_lagrangian_reaction_force = reaction%net_reaction_force
      loop%net_eulerian_force_buffer = spreading%net_eulerian_force
      loop%signature = 0.0_dp
      loop%signature(1) = attach%sum_abs_sampled_velocity
      loop%signature(2) = hydro%sum_abs_candidate_force
      loop%signature(3) = handoff%sum_abs_input_force
      loop%signature(4) = dry%sum_abs_dx_trial
      loop%signature(5:7) = reaction%net_reaction_force
      loop%signature(8:10) = spreading%net_eulerian_force
      loop%signature(11) = spreading%sum_abs_force_buffer
      loop%signature(12) = rhs_gate%sum_abs_rhs_increment
      loop%signature(13) = rhs_gate%max_abs_rhs_increment
      loop%signature(14) = lambda_fsi
      loop%signature(15) = penalty_beta
      loop%signature(16) = maxval(abs(spreading%conservation_error)) + rhs_gate%measured_scale_error
      loop%completed = .true.
    end if

    call fibre_prod_force_buffer_rhs_gate_finalize(rhs_gate)
    call fibre_prod_reaction_spreading_buffer_finalize(spreading)
    call fibre_prod_force_buffer_destroy(force_buffer)
    call fibre_prod_reaction_force_candidate_finalize(reaction)
    call fibre_prod_structure_commit_gate_finalize(commit_gate)
    call fibre_prod_structure_dry_step_finalize(dry)
    call fibre_prod_structure_input_handoff_finalize(handoff)
    call fibre_prod_hydro_input_candidate_finalize(hydro)
    call fibre_prod_state_velocity_attachment_finalize(attach)
    call fibre_prod_state_destroy(state)
    call fibre_prod_grid_destroy(grid)
  end subroutine fibre_prod_synthetic_closed_loop_run

  subroutine fibre_prod_synthetic_closed_loop_get_signature(loop, signature, status)
    type(fibre_prod_synthetic_closed_loop_type), intent(in) :: loop
    real(dp), intent(out) :: signature(:)
    integer, intent(out) :: status

    status = 0
    if (.not. loop%initialized .or. .not. loop%completed) then
      status = 50
      return
    end if
    if (size(signature) < 16) then
      status = 51
      return
    end if
    signature(1:16) = loop%signature
  end subroutine fibre_prod_synthetic_closed_loop_get_signature

  subroutine fibre_prod_synthetic_closed_loop_finalize(loop)
    type(fibre_prod_synthetic_closed_loop_type), intent(inout) :: loop

    loop%initialized = .false.
    loop%completed = .false.
    loop%nx = 0
    loop%ny = 0
    loop%nz = 0
    loop%nnode = 0
    loop%lambda_fsi = 0.0_dp
    loop%penalty_beta = 1.0_dp
    loop%beta_hydro = 1.0_dp
    loop%dt = 1.0e-4_dp
    loop%rho_eff = 1.0_dp
    loop%max_allowed_dx = 1.0e-2_dp
    loop%max_abs_sampled_u = 0.0_dp
    loop%max_abs_hydro_candidate = 0.0_dp
    loop%max_abs_structure_input = 0.0_dp
    loop%max_abs_dx_trial = 0.0_dp
    loop%max_abs_reaction_force = 0.0_dp
    loop%max_abs_force_buffer = 0.0_dp
    loop%max_abs_rhs_increment = 0.0_dp
    loop%sum_abs_rhs_increment = 0.0_dp
    loop%net_lagrangian_reaction_force = 0.0_dp
    loop%net_eulerian_force_buffer = 0.0_dp
    loop%signature = 0.0_dp
  end subroutine fibre_prod_synthetic_closed_loop_finalize

  logical function fibre_prod_synthetic_closed_loop_env_enabled() result(enabled)
    enabled = read_env_logical('FIBRE_PROD_SYNTHETIC_CLOSED_LOOP_ENABLE')
  end function fibre_prod_synthetic_closed_loop_env_enabled

  logical function fibre_prod_synthetic_closed_loop_diagnostics_env_enabled() result(enabled)
    enabled = read_env_logical('FIBRE_PROD_SYNTHETIC_CLOSED_LOOP_DIAGNOSTICS')
  end function fibre_prod_synthetic_closed_loop_diagnostics_env_enabled

  subroutine fibre_prod_synthetic_closed_loop_runtime_diagnostic(status)
    integer, intent(out) :: status
    type(fibre_prod_synthetic_closed_loop_type) :: loop

    call fibre_prod_synthetic_closed_loop_init(loop, 4, 4, 4, 5, status)
    if (status == 0) call fibre_prod_synthetic_closed_loop_run(loop, 0.0_dp, 2.0_dp, 1.0_dp, 1.0e-4_dp, 1.0_dp, status)
    call fibre_prod_synthetic_closed_loop_finalize(loop)
  end subroutine fibre_prod_synthetic_closed_loop_runtime_diagnostic

  subroutine fill_unit_coordinates(coord)
    real(dp), intent(out) :: coord(:)
    integer :: i
    if (size(coord) == 1) then
      coord = 0.0_dp
    else
      do i = 1, size(coord)
        coord(i) = real(i - 1, dp) / real(size(coord) - 1, dp)
      end do
    end if
  end subroutine fill_unit_coordinates

  subroutine fill_velocity_field(xcoord, ycoord, zcoord, ux, uy, uz)
    real(dp), intent(in) :: xcoord(:), ycoord(:), zcoord(:)
    real(dp), intent(out) :: ux(:, :, :), uy(:, :, :), uz(:, :, :)
    integer :: i, j, k
    do k = 1, size(ux, 3)
      do j = 1, size(ux, 2)
        do i = 1, size(ux, 1)
          ux(i,j,k) = 0.1_dp + xcoord(i) + 0.25_dp*ycoord(j) + 0.125_dp*zcoord(k)
          uy(i,j,k) = -0.2_dp + 0.5_dp*xcoord(i) - 0.3_dp*ycoord(j) + 0.2_dp*zcoord(k)
          uz(i,j,k) = 0.3_dp - 0.125_dp*xcoord(i) + 0.4_dp*ycoord(j) + 0.6_dp*zcoord(k)
        end do
      end do
    end do
  end subroutine fill_velocity_field

  subroutine fill_fibre_points(points)
    real(dp), intent(out) :: points(:, :)
    integer :: i
    do i = 1, size(points, 1)
      points(i,1) = 0.2_dp + 0.6_dp * real(i - 1, dp) / real(max(1, size(points, 1) - 1), dp)
      points(i,2) = 0.35_dp + 0.04_dp * real(mod(i - 1, 3), dp)
      points(i,3) = 0.45_dp + 0.03_dp * real(mod(i, 2), dp)
    end do
  end subroutine fill_fibre_points

  subroutine fill_rhs_seed(rhs_x, rhs_y, rhs_z)
    real(dp), intent(out) :: rhs_x(:, :, :), rhs_y(:, :, :), rhs_z(:, :, :)
    integer :: i, j, k
    do k = 1, size(rhs_x, 3)
      do j = 1, size(rhs_x, 2)
        do i = 1, size(rhs_x, 1)
          rhs_x(i,j,k) = 1.0e-2_dp * real(i + 2*j + 3*k, dp)
          rhs_y(i,j,k) = -2.0e-2_dp * real(2*i - j + k, dp)
          rhs_z(i,j,k) = 3.0e-2_dp * real(i - 2*j + k, dp)
        end do
      end do
    end do
  end subroutine fill_rhs_seed

  logical function read_env_logical(name) result(enabled)
    character(len=*), intent(in) :: name
    character(len=64) :: raw
    integer :: length, env_status
    enabled = .false.
    call get_environment_variable(name, raw, length=length, status=env_status)
    if (env_status /= 0 .or. length <= 0) return
    select case (adjustl(raw(1:min(length, len(raw)))))
    case ('1', 'T', 't', 'TRUE', 'true', 'True', 'YES', 'yes', 'ON', 'on')
      enabled = .true.
    end select
  end function read_env_logical

end module fibre_prod_synthetic_closed_loop
