include_directories(${2decomp_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/src)

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(xcompact3d 
                    acl_controller.f90 
                    acl_elem.f90 
                    acl_farm_controller.f90 
                    acl_model.f90
                    acl_out.f90  
                    acl_source.f90
                    acl_turb.f90
                    acl_utils.f90 
                    adm.f90
                    airfoils.f90
                    Case-ABL.f90
                    Case-Cavity.f90
                    Case-Channel.f90
                    Case-Cylinder-wake.f90
                    Case-Generic.f90
                    Case-Gravity-current.f90
                    Case-Mixing-layer.f90
                    Case-Periodic-hill.f90
                    Case-Pipe-flow.f90
                    Case-Sandbox.f90
                    Case-TBL.f90
                    Case-TGV.f90
                    Case-Uniform.f90
                    Case-PTBL.f90
                    case.f90
                    derive.f90
                    dynstall.f90
                    dynstall_legacy.f90
                    filters.f90
                    forces.f90
                    genepsi3d.f90
                    ibm.f90
                    implicit.f90
                    les_models.f90
                    mhd.f90
                    module_param.f90
                    mptool.f90
                    fibre_stage9_5_projection_regression.f90
                    fibre_stage9_6_rk3_rhs_massflux_regression.f90
                    fibre_stage9_7_stats_visu_io_smoke.f90
                    fibre_stage9_8_restart_io_regression.f90
                    fibre_stage9_9_parallel_consistency.f90
                    fibre_stage10_config.f90
                    fibre_stage10_noop_hook.f90
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage12_config.f90
                    fibre_stage12_production_feedback_candidate.f90
                    fibre_stage13_config.f90
                    fibre_stage13_spreading_kernel.f90
                    fibre_stage13_volume_normalization_audit.f90
                    fibre_stage13_conservation_sign_audit.f90
                    fibre_stage13_production_force_density_candidate.f90
                    fibre_stage14_config.f90
                    fibre_stage14_production_rhs_injection.f90
                    fibre_stage15_config.f90
                    fibre_stage15_structure_state.f90
                    fibre_stage15_production_structure_hook.f90
                    fibre_prod_state.f90
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_ibm_spreading.f90
                    fibre_prod_boundary_conditions.f90
                    fibre_prod_bending_solver.f90
                    fibre_prod_tension_solver.f90
                    fibre_prod_structure_solver.f90
                    fibre_prod_fsi_config.f90
                    fibre_prod_fluid_surrogate.f90
                    fibre_prod_fsi_coupling.f90
                    fibre_prod_fsi_diagnostics.f90
                    fibre_prod_wall_geometry.f90
                    fibre_prod_wall_contact.f90
                    fibre_prod_collision_geometry.f90
                    fibre_prod_fibre_collision.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_velocity_bridge.f90
                    fibre_prod_state_velocity_attachment.f90
                    fibre_prod_hydro_input_candidate.f90
                    fibre_prod_structure_input_handoff.f90
                    navier.f90
                    parameters.f90
                    particle.f90
                    poisson.f90
                    probes.f90
                    schemes.f90
                    statistics.f90
		    stretching.f90
                    time_integrators.f90
                    tools.f90
                    transeq.f90
                    turbine.f90
                    utilities.f90
                    variables.f90
                    visu.f90
                    fibre_stage9_3_channel_init_dryrun.f90
                    fibre_stage9_4_no_fibre_dns_smoke.f90
                    xcompact3d.f90)


add_executable(fibre_prod_state_check
                    fibre_prod_config.f90
                    fibre_prod_state.f90
                    fibre_prod_diagnostics.f90
                    fibre_prod_state_check.f90)


add_executable(fibre_prod_grid_adapter_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_grid_adapter_check.f90)


add_executable(fibre_prod_ibm_interpolation_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_interpolation_check.f90)


add_executable(fibre_prod_structure_solver_check
                    fibre_prod_config.f90
                    fibre_prod_state.f90
                    fibre_prod_diagnostics.f90
                    fibre_prod_boundary_conditions.f90
                    fibre_prod_bending_solver.f90
                    fibre_prod_tension_solver.f90
                    fibre_prod_structure_solver.f90
                    fibre_prod_structure_solver_check.f90)


add_executable(fibre_prod_ibm_spreading_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_ibm_spreading.f90
                    fibre_prod_ibm_spreading_check.f90)


add_executable(fibre_prod_fsi_closed_loop_check
                    fibre_prod_config.f90
                    fibre_prod_state.f90
                    fibre_prod_diagnostics.f90
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_ibm_spreading.f90
                    fibre_prod_boundary_conditions.f90
                    fibre_prod_bending_solver.f90
                    fibre_prod_tension_solver.f90
                    fibre_prod_structure_solver.f90
                    fibre_prod_fsi_config.f90
                    fibre_prod_fluid_surrogate.f90
                    fibre_prod_fsi_coupling.f90
                    fibre_prod_fsi_diagnostics.f90
                    fibre_prod_fsi_closed_loop_check.f90)


add_executable(fibre_prod_wall_contact_check
                    fibre_prod_config.f90
                    fibre_prod_state.f90
                    fibre_prod_diagnostics.f90
                    fibre_prod_boundary_conditions.f90
                    fibre_prod_bending_solver.f90
                    fibre_prod_tension_solver.f90
                    fibre_prod_structure_solver.f90
                    fibre_prod_wall_geometry.f90
                    fibre_prod_wall_contact.f90
                    fibre_prod_wall_contact_diagnostics.f90
                    fibre_prod_wall_contact_check.f90)

add_executable(fibre_prod_fibre_collision_check
                    fibre_prod_config.f90
                    fibre_prod_state.f90
                    fibre_prod_diagnostics.f90
                    fibre_prod_boundary_conditions.f90
                    fibre_prod_bending_solver.f90
                    fibre_prod_tension_solver.f90
                    fibre_prod_structure_solver.f90
                    fibre_prod_collision_geometry.f90
                    fibre_prod_fibre_collision.f90
                    fibre_prod_collision_diagnostics.f90
                    fibre_prod_fibre_collision_check.f90)



add_executable(fibre_prod_main_hook_check
                    fibre_prod_config.f90
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_main_hook_check.f90)


add_executable(fibre_prod_runtime_bridge_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_runtime_bridge_check.f90)


add_executable(fibre_prod_velocity_bridge_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_velocity_bridge.f90
                    fibre_prod_velocity_bridge_check.f90)


add_executable(fibre_prod_state_velocity_attachment_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_velocity_bridge.f90
                    fibre_prod_state.f90
                    fibre_prod_state_velocity_attachment.f90
                    fibre_prod_state_velocity_attachment_check.f90)



add_executable(fibre_prod_hydro_input_candidate_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_velocity_bridge.f90
                    fibre_prod_state.f90
                    fibre_prod_state_velocity_attachment.f90
                    fibre_prod_hydro_input_candidate.f90
                    fibre_prod_hydro_input_candidate_check.f90)



add_executable(fibre_prod_structure_input_handoff_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_interpolation.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_runtime_bridge.f90
                    fibre_prod_velocity_bridge.f90
                    fibre_prod_state.f90
                    fibre_prod_state_velocity_attachment.f90
                    fibre_prod_hydro_input_candidate.f90
                    fibre_prod_structure_input_handoff.f90
                    fibre_prod_structure_input_handoff_check.f90)

add_executable(fibre_prod_force_buffer_rhs_path_check
                    fibre_prod_grid_adapter.f90
                    fibre_prod_ibm_delta.f90
                    fibre_prod_ibm_force_buffer.f90
                    fibre_prod_ibm_spreading.f90
                    fibre_prod_runtime_config.f90
                    fibre_prod_main_diagnostics.f90
                    fibre_prod_rhs_adapter.f90
                    fibre_prod_main_hook.f90
                    fibre_prod_force_buffer_rhs_path_check.f90)

target_link_libraries(xcompact3d PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(xcompact3d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(ADIOS2_FOUND)
  #target_link_libraries(xcompact3d ${ADIOS_LIBRARIES})
  target_link_libraries(xcompact3d PRIVATE adios2::fortran_mpi)
endif()

add_executable(fibre_stage10_config_check
                    fibre_stage10_config.f90
                    fibre_stage10_config_check.f90)

target_link_libraries(fibre_stage10_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage10_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)



add_executable(fibre_stage10_noop_hook_check
                    fibre_stage10_config.f90
                    fibre_stage10_noop_hook.f90
                    fibre_stage10_noop_hook_check.f90)

target_link_libraries(fibre_stage10_noop_hook_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage10_noop_hook_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage11_config_check
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage11_config_check.f90)

target_link_libraries(fibre_stage11_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage11_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage11_lagrangian_state_check
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage11_lagrangian_state.f90
                    fibre_stage11_lagrangian_state_check.f90)

target_link_libraries(fibre_stage11_lagrangian_state_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage11_lagrangian_state_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage11_grid_metadata_check
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage11_lagrangian_state.f90
                    fibre_stage11_grid_metadata.f90
                    fibre_stage11_grid_metadata_check.f90)

target_link_libraries(fibre_stage11_grid_metadata_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage11_grid_metadata_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage11_oneway_interpolation_check
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage11_lagrangian_state.f90
                    fibre_stage11_grid_metadata.f90
                    fibre_stage11_oneway_interpolation.f90
                    fibre_stage11_oneway_interpolation_check.f90)

target_link_libraries(fibre_stage11_oneway_interpolation_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage11_oneway_interpolation_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage11_controlled_interpolation_check
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage11_lagrangian_state.f90
                    fibre_stage11_grid_metadata.f90
                    fibre_stage11_controlled_interpolation.f90
                    fibre_stage11_controlled_interpolation_check.f90)

target_link_libraries(fibre_stage11_controlled_interpolation_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage11_controlled_interpolation_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage11_production_oneway_hook_check
                    fibre_stage11_config.f90
                    fibre_stage11_production_oneway_hook.f90
                    fibre_stage11_production_oneway_hook_check.f90)

target_link_libraries(fibre_stage11_production_oneway_hook_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage11_production_oneway_hook_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage12_config_check
                    fibre_stage12_config.f90
                    fibre_stage12_config_check.f90)

target_link_libraries(fibre_stage12_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage12_force_buffer_check
                    fibre_stage12_config.f90
                    fibre_stage12_force_buffer.f90
                    fibre_stage12_force_buffer_check.f90)

target_link_libraries(fibre_stage12_force_buffer_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_force_buffer_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage12_prescribed_velocity_check
                    fibre_stage12_config.f90
                    fibre_stage12_prescribed_velocity.f90
                    fibre_stage12_prescribed_velocity_check.f90)

target_link_libraries(fibre_stage12_prescribed_velocity_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_prescribed_velocity_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage12_feedback_formula_check
                    fibre_stage12_config.f90
                    fibre_stage12_force_buffer.f90
                    fibre_stage12_prescribed_velocity.f90
                    fibre_stage12_feedback_formula.f90
                    fibre_stage12_feedback_formula_check.f90)

target_link_libraries(fibre_stage12_feedback_formula_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_feedback_formula_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage12_sign_convention_audit_check
                    fibre_stage12_config.f90
                    fibre_stage12_sign_convention_audit.f90
                    fibre_stage12_sign_convention_audit_check.f90)

target_link_libraries(fibre_stage12_sign_convention_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_sign_convention_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage12_power_diagnostics_check
                    fibre_stage12_config.f90
                    fibre_stage12_feedback_formula.f90
                    fibre_stage12_sign_convention_audit.f90
                    fibre_stage12_power_diagnostics.f90
                    fibre_stage12_power_diagnostics_check.f90)

target_link_libraries(fibre_stage12_power_diagnostics_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_power_diagnostics_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)


add_executable(fibre_stage12_production_feedback_candidate_check
                    fibre_stage12_config.f90
                    fibre_stage12_production_feedback_candidate.f90
                    fibre_stage12_production_feedback_candidate_check.f90)

target_link_libraries(fibre_stage12_production_feedback_candidate_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage12_production_feedback_candidate_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage13_config_check
                    fibre_stage13_config.f90
                    fibre_stage13_config_check.f90)

target_link_libraries(fibre_stage13_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage13_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage13_force_density_buffer_check
                    fibre_stage13_config.f90
                    fibre_stage13_force_density_buffer.f90
                    fibre_stage13_force_density_buffer_check.f90)

target_link_libraries(fibre_stage13_force_density_buffer_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage13_force_density_buffer_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage13_spreading_kernel_check
                    fibre_stage13_config.f90
                    fibre_stage13_force_density_buffer.f90
                    fibre_stage13_spreading_kernel.f90
                    fibre_stage13_spreading_kernel_check.f90)

target_link_libraries(fibre_stage13_spreading_kernel_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage13_spreading_kernel_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage13_volume_normalization_audit_check
                    fibre_stage13_config.f90
                    fibre_stage13_force_density_buffer.f90
                    fibre_stage13_spreading_kernel.f90
                    fibre_stage13_volume_normalization_audit.f90
                    fibre_stage13_volume_normalization_audit_check.f90)

target_link_libraries(fibre_stage13_volume_normalization_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage13_volume_normalization_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage13_conservation_sign_audit_check
                    fibre_stage13_config.f90
                    fibre_stage13_force_density_buffer.f90
                    fibre_stage13_spreading_kernel.f90
                    fibre_stage13_volume_normalization_audit.f90
                    fibre_stage13_conservation_sign_audit.f90
                    fibre_stage13_conservation_sign_audit_check.f90)

target_link_libraries(fibre_stage13_conservation_sign_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage13_conservation_sign_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage13_production_force_density_candidate_check
                    fibre_stage13_config.f90
                    fibre_stage13_spreading_kernel.f90
                    fibre_stage13_volume_normalization_audit.f90
                    fibre_stage13_conservation_sign_audit.f90
                    fibre_stage13_production_force_density_candidate.f90
                    fibre_stage13_production_force_density_candidate_check.f90)

target_link_libraries(fibre_stage13_production_force_density_candidate_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage13_production_force_density_candidate_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage14_config_check
                    fibre_stage14_config.f90
                    fibre_stage14_config_check.f90)

target_link_libraries(fibre_stage14_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage14_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_config_check
                    fibre_stage15_config.f90
                    fibre_stage15_config_check.f90)

target_link_libraries(fibre_stage15_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_config_check
                    fibre_stage16_config.f90
                    fibre_stage16_config_check.f90)

target_link_libraries(fibre_stage16_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_one_fibre_case_check
                    fibre_stage16_config.f90
                    fibre_stage16_one_fibre_case.f90
                    fibre_stage16_one_fibre_case_check.f90)

target_link_libraries(fibre_stage16_one_fibre_case_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_one_fibre_case_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_force_sign_audit_check
                    fibre_stage16_force_sign_audit.f90
                    fibre_stage16_force_sign_audit_check.f90)

target_link_libraries(fibre_stage16_force_sign_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_force_sign_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_structure_force_input_check
                    fibre_stage16_structure_force_input.f90
                    fibre_stage16_structure_force_input_check.f90)

target_link_libraries(fibre_stage16_structure_force_input_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_structure_force_input_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_closed_loop_dryrun_check
                    fibre_stage16_structure_force_input.f90
                    fibre_stage16_closed_loop_dryrun.f90
                    fibre_stage16_closed_loop_dryrun_check.f90)

target_link_libraries(fibre_stage16_closed_loop_dryrun_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_closed_loop_dryrun_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_lambda0_no_contamination_check
                    fibre_stage16_structure_force_input.f90
                    fibre_stage16_lambda0_no_contamination.f90
                    fibre_stage16_lambda0_no_contamination_check.f90)

target_link_libraries(fibre_stage16_lambda0_no_contamination_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_lambda0_no_contamination_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage16_small_lambda_response_check
                    fibre_stage16_structure_force_input.f90
                    fibre_stage16_small_lambda_response.f90
                    fibre_stage16_small_lambda_response_check.f90)

target_link_libraries(fibre_stage16_small_lambda_response_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage16_small_lambda_response_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_structure_state_check
                    fibre_stage15_config.f90
                    fibre_stage15_structure_state.f90
                    fibre_stage15_structure_state_check.f90)

target_link_libraries(fibre_stage15_structure_state_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_structure_state_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_velocity_source_adapter_check
                    fibre_stage12_prescribed_velocity.f90
                    fibre_stage12_feedback_formula.f90
                    fibre_stage15_config.f90
                    fibre_stage15_structure_state.f90
                    fibre_stage15_velocity_source_adapter.f90
                    fibre_stage15_velocity_source_adapter_check.f90)

target_link_libraries(fibre_stage15_velocity_source_adapter_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_velocity_source_adapter_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_structure_advance_formula_check
                    fibre_stage15_config.f90
                    fibre_stage15_structure_advance_formula.f90
                    fibre_stage15_structure_advance_formula_check.f90)

target_link_libraries(fibre_stage15_structure_advance_formula_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_structure_advance_formula_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_production_structure_hook_check
                    fibre_stage15_config.f90
                    fibre_stage15_structure_state.f90
                    fibre_stage15_production_structure_hook.f90
                    fibre_stage15_production_structure_hook_check.f90)

target_link_libraries(fibre_stage15_production_structure_hook_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_production_structure_hook_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_controlled_structure_step_check
                    fibre_stage15_controlled_structure_step.f90
                    fibre_stage15_controlled_structure_step_check.f90)

target_link_libraries(fibre_stage15_controlled_structure_step_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_controlled_structure_step_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage15_feedback_linkage_check
                    fibre_stage12_feedback_formula.f90
                    fibre_stage15_feedback_linkage_check.f90)

target_link_libraries(fibre_stage15_feedback_linkage_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage15_feedback_linkage_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage14_rhs_accumulator_check
                    fibre_stage14_config.f90
                    fibre_stage14_rhs_accumulator.f90
                    fibre_stage14_rhs_accumulator_check.f90)

target_link_libraries(fibre_stage14_rhs_accumulator_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage14_rhs_accumulator_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage14_rhs_addition_formula_check
                    fibre_stage14_config.f90
                    fibre_stage14_rhs_accumulator.f90
                    fibre_stage14_rhs_addition_formula.f90
                    fibre_stage14_rhs_addition_formula_check.f90)

target_link_libraries(fibre_stage14_rhs_addition_formula_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage14_rhs_addition_formula_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage14_rhs_sign_scaling_audit_check
                    fibre_stage14_config.f90
                    fibre_stage14_rhs_accumulator.f90
                    fibre_stage14_rhs_addition_formula.f90
                    fibre_stage14_rhs_sign_scaling_audit.f90
                    fibre_stage14_rhs_sign_scaling_audit_check.f90)

target_link_libraries(fibre_stage14_rhs_sign_scaling_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage14_rhs_sign_scaling_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage14_rk_timing_audit_check
                    fibre_stage14_config.f90
                    fibre_stage14_rhs_accumulator.f90
                    fibre_stage14_rk_timing_audit.f90
                    fibre_stage14_rk_timing_audit_check.f90)

target_link_libraries(fibre_stage14_rk_timing_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage14_rk_timing_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_stage14_production_rhs_injection_check
                    fibre_stage14_config.f90
                    fibre_stage14_production_rhs_injection.f90
                    fibre_stage14_production_rhs_injection_check.f90)

target_link_libraries(fibre_stage14_production_rhs_injection_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage14_production_rhs_injection_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

add_executable(fibre_structure_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_diagnostics.f90
                    fibre_structure_driver.f90)

target_link_libraries(fibre_structure_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_structure_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_structure_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_freefree_boundary_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_freefree_boundary_driver.f90)

target_link_libraries(fibre_freefree_boundary_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_freefree_boundary_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_freefree_boundary_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_bending_operator_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_bending_operator_driver.f90)

target_link_libraries(fibre_bending_operator_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_bending_operator_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_bending_operator_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_energy_diagnostics_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_diagnostics.f90
                    fibre_energy_diagnostics_driver.f90)

target_link_libraries(fibre_energy_diagnostics_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_energy_diagnostics_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_energy_diagnostics_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_implicit_bending_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_implicit_bending.f90
                    fibre_implicit_bending_driver.f90)

target_link_libraries(fibre_implicit_bending_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_implicit_bending_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_implicit_bending_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_tension_solver_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_tension_solver.f90
                    fibre_tension_solver_driver.f90)

target_link_libraries(fibre_tension_solver_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_tension_solver_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_tension_solver_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_structure_advance_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_diagnostics.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_implicit_bending.f90
                    fibre_tension_solver.f90
                    fibre_structure_integrator.f90
                    fibre_structure_advance_driver.f90)

target_link_libraries(fibre_structure_advance_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_structure_advance_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_structure_advance_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_no_force_long_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_diagnostics.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_implicit_bending.f90
                    fibre_tension_solver.f90
                    fibre_structure_integrator.f90
                    fibre_no_force_long_driver.f90)

target_link_libraries(fibre_no_force_long_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_no_force_long_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_no_force_long_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_frequency_scaling_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_diagnostics.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_implicit_bending.f90
                    fibre_tension_solver.f90
                    fibre_structure_integrator.f90
                    fibre_frequency_scaling_driver.f90)

target_link_libraries(fibre_frequency_scaling_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_frequency_scaling_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_frequency_scaling_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_external_force_interface_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_external_force.f90
                    fibre_geometry.f90
                    fibre_diagnostics.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_implicit_bending.f90
                    fibre_tension_solver.f90
                    fibre_structure_integrator.f90
                    fibre_external_force_interface_driver.f90)

target_link_libraries(fibre_external_force_interface_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_external_force_interface_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_external_force_interface_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_ibm_skeleton_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_stencil.f90
                    fibre_ibm_diagnostics.f90
                    fibre_ibm_skeleton_check.f90)

target_link_libraries(fibre_ibm_skeleton_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_skeleton_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_skeleton_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_ibm_delta_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_delta_check.f90)

target_link_libraries(fibre_ibm_delta_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_delta_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_delta_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_ibm_interpolation_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_interpolation_check.f90)

target_link_libraries(fibre_ibm_interpolation_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_interpolation_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_interpolation_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_ibm_spreading_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_spreading_check.f90)

target_link_libraries(fibre_ibm_spreading_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_spreading_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_spreading_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_ibm_power_consistency_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_power_consistency_check.f90)

target_link_libraries(fibre_ibm_power_consistency_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_power_consistency_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_power_consistency_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_ibm_feedback_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_feedback_check.f90)

target_link_libraries(fibre_ibm_feedback_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_feedback_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_feedback_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)



add_executable(fibre_ibm_structure_coupling_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_structure_coupling.f90
                    fibre_ibm_structure_coupling_check.f90)

target_link_libraries(fibre_ibm_structure_coupling_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_structure_coupling_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_structure_coupling_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)



add_executable(fibre_ibm_force_buffer_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_force_buffer_check.f90)

target_link_libraries(fibre_ibm_force_buffer_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_force_buffer_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_force_buffer_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)



add_executable(fibre_ibm_boundary_safety_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_ibm_boundary_safety_check.f90)

target_link_libraries(fibre_ibm_boundary_safety_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_boundary_safety_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_boundary_safety_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)



add_executable(fibre_ibm_stage3_smoke_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_structure_coupling.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_ibm_stage3_smoke_check.f90)

target_link_libraries(fibre_ibm_stage3_smoke_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_ibm_stage3_smoke_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_ibm_stage3_smoke_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_oneway_skeleton_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_force_buffer.f90
                    fibre_stage4_config.f90
                    fibre_stage4_diagnostics.f90
                    fibre_stage4_oneway_skeleton_check.f90)

target_link_libraries(fibre_stage4_oneway_skeleton_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_oneway_skeleton_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_oneway_skeleton_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_grid_adapter_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_stage4_config.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_grid_adapter_check.f90)

target_link_libraries(fibre_stage4_grid_adapter_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_grid_adapter_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_grid_adapter_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_interpolation_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_interpolation_check.f90)

target_link_libraries(fibre_stage4_interpolation_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_interpolation_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_interpolation_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_feedback_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_feedback.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_feedback_adapter.f90
                    fibre_stage4_feedback_check.f90)

target_link_libraries(fibre_stage4_feedback_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_feedback_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_feedback_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_poiseuille_response_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_feedback_adapter.f90
                    fibre_stage4_analytic_channel.f90
                    fibre_stage4_oneway_response.f90
                    fibre_stage4_poiseuille_response_check.f90)

target_link_libraries(fibre_stage4_poiseuille_response_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_poiseuille_response_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_poiseuille_response_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_frozen_response_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_feedback_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_oneway_response.f90
                    fibre_stage4_frozen_response_check.f90)

target_link_libraries(fibre_stage4_frozen_response_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_frozen_response_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_frozen_response_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_time_order_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_feedback_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_oneway_response.f90
                    fibre_stage4_time_order.f90
                    fibre_stage4_time_order_check.f90)

target_link_libraries(fibre_stage4_time_order_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_time_order_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_time_order_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_main_contamination_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_force_buffer.f90
                    fibre_stage4_config.f90
                    fibre_stage4_fluid_rhs_safety.f90
                    fibre_stage4_main_contamination_check.f90)

target_link_libraries(fibre_stage4_main_contamination_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_main_contamination_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_main_contamination_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_boundary_safety_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage4_boundary_safety_check.f90)

target_link_libraries(fibre_stage4_boundary_safety_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_boundary_safety_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_boundary_safety_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4_smoke_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_config.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_feedback_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_oneway_response.f90
                    fibre_stage4_time_order.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage4_fluid_rhs_safety.f90
                    fibre_stage4_smoke_check.f90)

target_link_libraries(fibre_stage4_smoke_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4_smoke_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4_smoke_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage4p_benchmark_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_feedback.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_config.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_feedback_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_oneway_response.f90
                    fibre_stage4_fluid_rhs_safety.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage4p_benchmark_check.f90)

target_link_libraries(fibre_stage4p_benchmark_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage4p_benchmark_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage4p_benchmark_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_config_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage5_config.f90
                    fibre_stage5_config_check.f90)

target_link_libraries(fibre_stage5_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_config_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_rhs_audit_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_audit_check.f90)

target_link_libraries(fibre_stage5_rhs_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_rhs_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_rhs_audit_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_rhs_wrapper_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_rhs_wrapper_check.f90)

target_link_libraries(fibre_stage5_rhs_wrapper_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_rhs_wrapper_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_rhs_wrapper_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_spreading_candidate_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_spreading_candidate_check.f90)

target_link_libraries(fibre_stage5_spreading_candidate_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_spreading_candidate_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_spreading_candidate_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_rhs_injection_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_rhs_injection_check.f90)

target_link_libraries(fibre_stage5_rhs_injection_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_rhs_injection_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_rhs_injection_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_momentum_exchange_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_momentum_exchange.f90
                    fibre_stage5_momentum_exchange_check.f90)

target_link_libraries(fibre_stage5_momentum_exchange_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_momentum_exchange_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_momentum_exchange_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_coupled_order_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_momentum_exchange.f90
                    fibre_stage5_coupled_order.f90
                    fibre_stage5_coupled_order_check.f90)

target_link_libraries(fibre_stage5_coupled_order_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_coupled_order_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_coupled_order_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage5_rk_substep_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_external_force.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_rk_substep_sync.f90
                    fibre_stage5_rk_substep_check.f90)

target_link_libraries(fibre_stage5_rk_substep_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_rk_substep_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_rk_substep_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage5_closed_loop_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_momentum_exchange.f90
                    fibre_stage5_coupled_order.f90
                    fibre_stage5_rk_substep_sync.f90
                    fibre_stage5_closed_loop.f90
                    fibre_stage5_closed_loop_check.f90)

target_link_libraries(fibre_stage5_closed_loop_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_closed_loop_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_closed_loop_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage5_main_hook_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_main_hook.f90
                    fibre_stage5_main_hook_check.f90)

target_link_libraries(fibre_stage5_main_hook_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_main_hook_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_main_hook_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage5_smoke_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_diagnostics.f90
                    fibre_external_force.f90
                    fibre_tension_solver.f90
                    fibre_implicit_bending.f90
                    fibre_structure_integrator.f90
                    fibre_ibm_types.f90
                    fibre_ibm_grid.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_ibm_force_buffer.f90
                    fibre_ibm_boundary_safety.f90
                    fibre_stage4_grid_adapter.f90
                    fibre_stage4_interpolation_adapter.f90
                    fibre_stage4_frozen_channel.f90
                    fibre_stage4_boundary_policy.f90
                    fibre_stage5_config.f90
                    fibre_stage5_rhs_audit.f90
                    fibre_stage5_rhs_wrapper.f90
                    fibre_stage5_spreading_candidate.f90
                    fibre_stage5_momentum_exchange.f90
                    fibre_stage5_coupled_order.f90
                    fibre_stage5_rk_substep_sync.f90
                    fibre_stage5_closed_loop.f90
                    fibre_stage5_main_hook.f90
                    fibre_stage5_smoke.f90
                    fibre_stage5_smoke_check.f90)

target_link_libraries(fibre_stage5_smoke_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage5_smoke_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage5_smoke_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_config_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_config_check.f90)

target_link_libraries(fibre_stage6_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_config_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_rhs_audit_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_rhs_audit_check.f90)

target_link_libraries(fibre_stage6_rhs_audit_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_rhs_audit_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_rhs_audit_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_noop_hook_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_main_noop_hook.f90
                    fibre_stage6_noop_hook_check.f90)

target_link_libraries(fibre_stage6_noop_hook_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_noop_hook_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_noop_hook_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_controlled_rhs_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage6_controlled_rhs_check.f90)

target_link_libraries(fibre_stage6_controlled_rhs_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_controlled_rhs_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_controlled_rhs_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_singlephase_noop_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_main_noop_hook.f90
                    fibre_stage6_singlephase_noop.f90
                    fibre_stage6_singlephase_noop_check.f90)

target_link_libraries(fibre_stage6_singlephase_noop_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_singlephase_noop_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_singlephase_noop_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_projection_order_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage6_projection_order.f90
                    fibre_stage6_projection_order_check.f90)

target_link_libraries(fibre_stage6_projection_order_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_projection_order_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_projection_order_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_rk_rhs_sync_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage6_rk_rhs_sync.f90
                    fibre_stage6_rk_rhs_sync_check.f90)

target_link_libraries(fibre_stage6_rk_rhs_sync_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_rk_rhs_sync_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_rk_rhs_sync_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage6_layout_guard_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage6_layout_guard.f90
                    fibre_stage6_layout_guard_check.f90)

target_link_libraries(fibre_stage6_layout_guard_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_layout_guard_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_layout_guard_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_micro_smoke_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage6_projection_order.f90
                    fibre_stage6_layout_guard.f90
                    fibre_stage6_micro_smoke.f90
                    fibre_stage6_micro_smoke_check.f90)

target_link_libraries(fibre_stage6_micro_smoke_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_micro_smoke_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_micro_smoke_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_io_safety_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_io_safety.f90
                    fibre_stage6_io_safety_check.f90)

target_link_libraries(fibre_stage6_io_safety_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_io_safety_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_io_safety_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage6_total_smoke_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage6_projection_order.f90
                    fibre_stage6_layout_guard.f90
                    fibre_stage6_io_safety.f90
                    fibre_stage6_total_smoke.f90
                    fibre_stage6_total_smoke_check.f90)

target_link_libraries(fibre_stage6_total_smoke_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage6_total_smoke_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage6_total_smoke_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage7_config_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_config_check.f90)

target_link_libraries(fibre_stage7_config_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage7_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage7_config_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage7_grid_metadata_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_grid_metadata_check.f90)

target_link_libraries(fibre_stage7_grid_metadata_check PRIVATE decomp2d)

if (MPI_FOUND)
  target_link_libraries(fibre_stage7_grid_metadata_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS fibre_stage7_grid_metadata_check
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


add_executable(fibre_stage7_feedback_sign_audit
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_feedback.f90
                    fibre_stage7_feedback_sign_audit.f90)
target_link_libraries(fibre_stage7_feedback_sign_audit PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_feedback_sign_audit PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_feedback_sign_audit RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_rho_convention_audit
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage7_rho_convention_audit.f90)
target_link_libraries(fibre_stage7_rho_convention_audit PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_rho_convention_audit PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_rho_convention_audit RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_scalar_interpolation_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_scalar_interpolation_check.f90)
target_link_libraries(fibre_stage7_scalar_interpolation_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_scalar_interpolation_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_scalar_interpolation_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_scalar_interpolation_robustness_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_scalar_interpolation_robustness_check.f90)
target_link_libraries(fibre_stage7_scalar_interpolation_robustness_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_scalar_interpolation_robustness_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_scalar_interpolation_robustness_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_velocity_interpolation_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_velocity_interpolation.f90
                    fibre_stage7_velocity_interpolation_check.f90)
target_link_libraries(fibre_stage7_velocity_interpolation_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_velocity_interpolation_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_velocity_interpolation_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_force_spreading_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_force_spreading.f90
                    fibre_stage7_force_spreading_check.f90)
target_link_libraries(fibre_stage7_force_spreading_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_force_spreading_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_force_spreading_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_power_consistency_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_velocity_interpolation.f90
                    fibre_stage7_force_spreading.f90
                    fibre_stage7_power_consistency_check.f90)
target_link_libraries(fibre_stage7_power_consistency_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_power_consistency_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_power_consistency_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_boundary_safety_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_velocity_interpolation.f90
                    fibre_stage7_force_spreading.f90
                    fibre_stage7_boundary_safety.f90
                    fibre_stage7_boundary_safety_check.f90)
target_link_libraries(fibre_stage7_boundary_safety_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_boundary_safety_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_boundary_safety_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_channel_grid_adapter_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_velocity_interpolation.f90
                    fibre_stage7_force_spreading.f90
                    fibre_stage7_boundary_safety.f90
                    fibre_stage7_channel_grid_adapter.f90
                    fibre_stage7_channel_grid_adapter_check.f90)
target_link_libraries(fibre_stage7_channel_grid_adapter_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_channel_grid_adapter_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_channel_grid_adapter_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_rhs_candidate_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage6_config.f90
                    fibre_stage6_rhs_audit.f90
                    fibre_stage6_controlled_rhs_hook.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_velocity_interpolation.f90
                    fibre_stage7_force_spreading.f90
                    fibre_stage7_boundary_safety.f90
                    fibre_stage7_channel_grid_adapter.f90
                    fibre_stage7_rhs_candidate.f90
                    fibre_stage7_rhs_candidate_check.f90)
target_link_libraries(fibre_stage7_rhs_candidate_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_rhs_candidate_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_rhs_candidate_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


add_executable(fibre_stage7_total_smoke_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_total_smoke_check.f90)
target_link_libraries(fibre_stage7_total_smoke_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_total_smoke_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_total_smoke_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


add_executable(fibre_stage8_config_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage8_config.f90
                    fibre_stage8_config_check.f90)
target_link_libraries(fibre_stage8_config_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage8_config_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage8_config_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


add_executable(fibre_stage8_runtime_grid_bridge_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_stage7_config.f90
                    fibre_stage7_grid_metadata.f90
                    fibre_stage7_channel_grid_adapter.f90
                    fibre_stage7_scalar_interpolation.f90
                    fibre_stage7_velocity_interpolation.f90
                    fibre_stage7_force_spreading.f90
                    fibre_stage7_boundary_safety.f90
                    fibre_stage8_config.f90
                    fibre_stage8_runtime_grid_bridge.f90
                    fibre_stage8_runtime_grid_bridge_check.f90)
target_link_libraries(fibre_stage8_runtime_grid_bridge_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage8_runtime_grid_bridge_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage8_runtime_grid_bridge_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage7_structure_high_gamma_check
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_geometry.f90
                    fibre_boundary_freefree.f90
                    fibre_bending_operator.f90
                    fibre_implicit_bending.f90
                    fibre_tension_solver.f90
                    fibre_diagnostics.f90
                    fibre_structure_integrator.f90
                    fibre_stage7_structure_high_gamma_check.f90)
target_link_libraries(fibre_stage7_structure_high_gamma_check PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage7_structure_high_gamma_check PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage7_structure_high_gamma_check RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(fibre_stage5_closed_loop_diagnostics_audit
                    fibre_types.f90
                    fibre_parameters.f90
                    fibre_ibm_types.f90
                    fibre_ibm_delta.f90
                    fibre_ibm_interpolation.f90
                    fibre_ibm_spreading.f90
                    fibre_ibm_power_diagnostics.f90
                    fibre_stage5_closed_loop_diagnostics_audit.f90)
target_link_libraries(fibre_stage5_closed_loop_diagnostics_audit PRIVATE decomp2d)
if (MPI_FOUND)
  target_link_libraries(fibre_stage5_closed_loop_diagnostics_audit PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
install(TARGETS fibre_stage5_closed_loop_diagnostics_audit RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS xcompact3d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(fibre_stage8_lagrangian_state_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_lagrangian_state_check.f90
)

target_link_libraries(fibre_stage8_lagrangian_state_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_lagrangian_state_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_lagrangian_state_check DESTINATION bin)

add_executable(fibre_stage8_velocity_to_fibre_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_velocity_to_fibre_check.f90
)

target_link_libraries(fibre_stage8_velocity_to_fibre_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_velocity_to_fibre_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_velocity_to_fibre_check DESTINATION bin)

add_executable(fibre_stage8_feedback_candidate_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_feedback_candidate.f90
  fibre_stage8_feedback_candidate_check.f90
)

target_link_libraries(fibre_stage8_feedback_candidate_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_feedback_candidate_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_feedback_candidate_check DESTINATION bin)

add_executable(fibre_stage8_oneway_forcing_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_feedback_candidate.f90
  fibre_stage8_oneway_forcing.f90
  fibre_stage8_oneway_forcing_check.f90
)

target_link_libraries(fibre_stage8_oneway_forcing_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_oneway_forcing_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_oneway_forcing_check DESTINATION bin)

add_executable(fibre_stage8_twoway_force_density_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_feedback_candidate.f90
  fibre_stage8_oneway_forcing.f90
  fibre_stage8_twoway_force_density.f90
  fibre_stage8_twoway_force_density_check.f90
)
target_link_libraries(fibre_stage8_twoway_force_density_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_twoway_force_density_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_twoway_force_density_check DESTINATION bin)

add_executable(fibre_stage8_rk_sync_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_feedback_candidate.f90
  fibre_stage8_oneway_forcing.f90
  fibre_stage8_twoway_force_density.f90
  fibre_stage8_rk_sync.f90
  fibre_stage8_rk_sync_check.f90
)
target_link_libraries(fibre_stage8_rk_sync_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_rk_sync_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_rk_sync_check DESTINATION bin)

add_executable(fibre_stage8_boundary_safe_workflow_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_feedback_candidate.f90
  fibre_stage8_oneway_forcing.f90
  fibre_stage8_twoway_force_density.f90
  fibre_stage8_rk_sync.f90
  fibre_stage8_boundary_safe_workflow.f90
  fibre_stage8_boundary_safe_workflow_check.f90
)
target_link_libraries(fibre_stage8_boundary_safe_workflow_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_boundary_safe_workflow_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_boundary_safe_workflow_check DESTINATION bin)

add_executable(fibre_stage8_integrated_audit_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage6_config.f90
  fibre_stage6_rhs_audit.f90
  fibre_stage6_controlled_rhs_hook.f90
  fibre_stage7_config.f90
  fibre_stage7_grid_metadata.f90
  fibre_stage7_channel_grid_adapter.f90
  fibre_stage7_scalar_interpolation.f90
  fibre_stage7_velocity_interpolation.f90
  fibre_stage7_force_spreading.f90
  fibre_stage7_boundary_safety.f90
  fibre_stage8_config.f90
  fibre_stage8_runtime_grid_bridge.f90
  fibre_stage8_lagrangian_state.f90
  fibre_stage8_velocity_to_fibre.f90
  fibre_stage8_feedback_candidate.f90
  fibre_stage8_oneway_forcing.f90
  fibre_stage8_twoway_force_density.f90
  fibre_stage8_rk_sync.f90
  fibre_stage8_boundary_safe_workflow.f90
  fibre_stage8_integrated_audit.f90
  fibre_stage8_integrated_audit_check.f90
)
target_link_libraries(fibre_stage8_integrated_audit_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_integrated_audit_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_integrated_audit_check DESTINATION bin)

add_executable(fibre_stage8_total_smoke_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage8_total_smoke_check.f90
)
target_link_libraries(fibre_stage8_total_smoke_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage8_total_smoke_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage8_total_smoke_check DESTINATION bin)

add_executable(fibre_stage9_dependency_gate_check
  fibre_types.f90
  fibre_parameters.f90
  fibre_stage9_dependency_gate_check.f90
)
target_link_libraries(fibre_stage9_dependency_gate_check PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage9_dependency_gate_check PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage9_dependency_gate_check DESTINATION bin)

add_executable(fibre_stage9_2_minimal_parallel_gate
  fibre_stage9_2_minimal_parallel_gate.f90
)
target_link_libraries(fibre_stage9_2_minimal_parallel_gate PRIVATE decomp2d)
if(MPI_FOUND)
  target_link_libraries(fibre_stage9_2_minimal_parallel_gate PRIVATE MPI::MPI_Fortran)
endif()
install(TARGETS fibre_stage9_2_minimal_parallel_gate DESTINATION bin)