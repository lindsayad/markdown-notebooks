Print position in vector in LLDB:
```
p elem_num_map.__begin_[8297]
```

How to access the multi app FEProblem's time from the master app's FEProblem:
```
p _multi_apps._active_objects.__begin_[0].__begin_[0].__ptr_->_apps.__begin_[0].__ptr_->_executioner.__ptr_->_fe_problem._time
```

`b __cxa_throw` sets a breakpoint at `throw`

To apply a patch after git crashes during a rebase:

git apply patch-name

PBS run_tests command:
./run_tests -c -q -t --no-report --sep-files --queue-project neams --pbs assess
-i assessment 2>&1 >> assess_log

etags moose:

find . \( \( -iname "*build*" -o -iname "*installed*" \) -prune \) -o \( -iname "*.h" -o -iname "*.C" -o -iname "*.c" \) -print | etags -

clang++ -std=c++11 -g -O0 -Iinclude -I$PETSC_DIR/include -o export_discontinuous libmesh_export_discotinuous_error.cxx exact_solution.C -Llib -Wl,-rpath -Wl,lib -lmesh_dbg

icecream:

ssh -Y icecream.inl.gov && icemon &

Textbook that I need to order: "Nonlinear Systems and Optimization for the
Chemical Engineer"

Action order list:
  [0] = "no_action"
  [1] = "setup_oversampling"
  [2] = "deprecated_block"
  [3] = "finish_input_file_output"
  [4] = "meta_action"
  [5] = "dynamic_object_registration"
  [6] = "common_output"
  [7] = "set_global_params"
  [8] = "setup_recover_file_base"
  [9] = "check_copy_nodal_vars"
  [10] = "setup_mesh"
  [11] = "add_partitioner"
  [12] = "add_geometric_rm"
  [13] = "init_mesh"
  [14] = "prepare_mesh"
  [15] = "add_mesh_modifier"
  [16] = "execute_mesh_modifiers"
  [17] = "add_mortar_interface"
  [18] = "uniform_refine_mesh"
  [19] = "setup_mesh_complete"
  [20] = "determine_system_type"
  [21] = "create_problem"
  [22] = "setup_postprocessor_data"
  [23] = "setup_time_integrator"
  [24] = "setup_executioner"
  [25] = "check_integrity_early"
  [26] = "setup_predictor"
  [27] = "init_displaced_problem"
  [28] = "add_elemental_field_variable"
  [29] = "add_aux_variable"
  [30] = "add_external_aux_variables"
  [31] = "add_variable"
  [32] = "setup_variable_complete"
  [33] = "setup_quadrature"
  [34] = "add_function"
  [35] = "add_distribution"
  [36] = "add_sampler"
  [37] = "add_periodic_bc"
  [38] = "add_user_object"
  [39] = "setup_function_complete"
  [40] = "setup_adaptivity"
  [41] = "set_adaptivity_options"
  [42] = "add_ic"
  [43] = "add_constraint"
  [44] = "add_field_split"
  [45] = "add_preconditioning"
  [46] = "setup_time_stepper"
  [47] = "ready_to_init"
  [48] = "setup_dampers"
  [49] = "setup_residual_debug"
  [50] = "add_bounds_vectors"
  [51] = "add_multi_app"
  [52] = "add_transfer"
  [53] = "copy_nodal_vars"
  [54] = "copy_nodal_aux_vars"
  [55] = "add_material"
  [56] = "setup_material_output"
  [57] = "add_algebraic_rm"
  [58] = "init_problem"
  [59] = "setup_debug"
  [60] = "add_output"
  [61] = "add_postprocessor"
  [62] = "add_vector_postprocessor"
  [63] = "add_kernel"
  [64] = "add_nodal_kernel"
  [65] = "add_bc"
  [66] = "add_aux_kernel"
  [67] = "add_scalar_kernel"
  [68] = "add_aux_scalar_kernel"
  [69] = "add_dirac_kernel"
  [70] = "add_dg_kernel"
  [71] = "add_interface_kernel"
  [72] = "add_damper"
  [73] = "add_indicator"
  [74] = "add_marker"
  [75] = "add_control"
  [76] = "check_output"
  [77] = "check_integrity"

Problem/System/Assembly notes, e.g. how MOOSE actually puts together the
residuals/Jacobians:

- Assembly::prepare: resizes (based on dofIndices) and zeros the local residual
  and Jacobian blocks
- Assembly::prepareNeighbor: resizes (based on dofIndicesNeighbor) and zeros the local residual
  and Jacobian blocks
- Assembly::reinitNeighborAtPhysical: calls:
  - Assembly::reinitFEFaceNeighbor(neighbor_side_elem, ref_points)
  - Assembly::reinitNeighbor(neighbor_side_elem, ref_points)
  - Assembly::reinitFEFaceNeighbor(neighbor_elem, ref_points)
  - Assembly::reinitNeighbor(neighbor_elem, ref_points)
- Assembly::reinitFE(elem): does shape function initialization
- Assembly::reinitFEFaceNeighbor(elem, ref_points): does shape function initialization (for proper elem->dim)
- Assembly::reinit(elem): sets the _current_elem member, sets up quadrature, and calls reinitFE(elem)
- Assembly::reinitNeighbor(elem, ref_points): sets the _current_neighbor_elem member and sets up quadrature


- FEProblemBase::reinitNeighborPhys: calls:
  - Assembly::reinitNeighborAtPhysical
  - *System*::prepareNeighbor
  - Assembly::prepareNeighbor
  - *System*::reinitNeighborFace
- FEProblemBase::reinitNodeFace(slave_node, slave_boundary): calls:
  - Assembly::reinit(slave_node)
  - *System*::reinitNodeFace(slave_node, slave_boundary)

- *System*::prepareNeighbor: calls
  - MooseVariable::prepareNeighbor
- *System*::reinitNeighborFace: calls:
  - MooseVariable::computeNeighborValuesFace
- *System*::reinitNodeFace: calls:
  - MooseVariable::reinitNode
  - MooseVariable::compouteNodalValues

- MooseVariable::computeNodalValues: computes variable values
- MooseVariable::computeNeighborValuesFace: computes variable values
- MooseVariable::prepareNeighbor: sets neighbor dof indices based on current neighbor element
- MooseVariable::reinitNode: sets dof indices based on current node

In summary, proper computing order (generally):
- Assembly::reinit: sets up shape functions, quadrature, and sets problem _elem, _node, _neighbor data
- *System*::prepare: sets up dof indices through calls to MooseVariable
- Assembly::prepare: resizes local residuals and Jacobians
- *System*::reinit: computes variable values through calls to MooseVariable

Note that:
- FEProblemBase::prepare calls (generally):
  - Assembly::reinit
  - *System*::prepare
  - Assembly::prepare
- FEProblemBase::reinit calls:
  - *System*::reinit

FE::reinit in libMesh summary of steps:
1. this->_fe_map->template init_reference_to_physical_map
2. this->init_shape_functions:
     This calculates the reference phi and grad_phi values, e.g. things like dphidxi
3. this->_fe_map->compute_map
4. this->compute_shape_functions
     This calculates the physical phi and grad_phi values, e.g. things like
     dphidx

AFAICT in MOOSE, we use a JxW based on first order Lagrange shape functions in
all our compute object calculations even though if we have second order elements
the mapping of shpae functions will be done using second order Lagrange shape
functions.

# 9/6/17

Failed tests in MOOSE from changing wp to ds:
```
auxkernels/nodal_aux_var.init_test.......................................................... FAILED (EXODIFF)
auxkernels/error_function_aux.error_function_aux............................................ FAILED (EXODIFF)
auxkernels/grad_component.test.............................................................. FAILED (EXODIFF)
auxkernels/time_derivative.implicit_euler................................................... FAILED (EXODIFF)
auxkernels/constant_scalar_aux.test......................................................... FAILED (EXODIFF)
auxkernels/function_scalar_aux.test......................................................... FAILED (EXODIFF)
controls/time_periods/nodalkernels.test..................................................... FAILED (EXODIFF)
executioners/eigen_executioners.test_nonlinear_eigen........................................ FAILED (EXODIFF)
executioners/eigen_executioners.test_nonlinear_eigen_material............................... FAILED (EXODIFF)
executioners/eigen_executioners.test_normal_eigenkernel..................................... FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3qns................................................ FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test4tt................................................. FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3tt................................................. FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3qnns............................................... FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3ns................................................. FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3................................................... FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test4ns................................................. FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test4................................................... FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3nns................................................ FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test3q.................................................. FAILED (EXODIFF)
geomsearch/2d_moving_penetration.pl_test4nns................................................ FAILED (EXODIFF)
geomsearch/3d_moving_penetration.pl_test1q.................................................. FAILED (EXODIFF)
geomsearch/3d_moving_penetration.pl_test2q.................................................. FAILED (EXODIFF)
geomsearch/3d_moving_penetration.pl_test4q.................................................. FAILED (EXODIFF)
geomsearch/3d_moving_penetration.pl_test3q.................................................. FAILED (EXODIFF)
geomsearch/3d_moving_penetration_smoothing.pl_test4ns....................................... FAILED (EXODIFF)
indicators/analytical_indicator.test........................................................ FAILED (EXODIFF)
geomsearch/3d_moving_penetration_smoothing.pl_test3ns....................................... FAILED (EXODIFF)
geomsearch/3d_moving_penetration_smoothing.pl_test3qnstt.................................... FAILED (EXODIFF)
geomsearch/3d_moving_penetration_smoothing.pl_test3qns...................................... FAILED (EXODIFF)
geomsearch/3d_moving_penetration_smoothing.pl_test4qnstt.................................... FAILED (EXODIFF)
geomsearch/3d_moving_penetration_smoothing.pl_test4qns...................................... FAILED (EXODIFF)
markers/error_fraction_marker.test.......................................................... FAILED (EXODIFF)
markers/error_tolerance_marker.test......................................................... FAILED (EXODIFF)
materials/material_dependency.dont_reinit_mat_for_aux....................................... FAILED (EXODIFF)
mesh/adapt_weight.test...................................................................... FAILED (EXODIFF)
mesh/adapt.test............................................................................. FAILED (EXODIFF)
nodalkernels/high_order_time_integration.test............................................... FAILED (EXODIFF)
nodalkernels/constant_rate.test............................................................. FAILED (EXODIFF)
outputs/iterative.exodus_steady............................................................. FAILED (EXODIFF)
outputs/iterative.exodus_steady_sequence.................................................... FAILED (EXODIFF)
outputs/output_interface.indicators......................................................... FAILED (EXODIFF)
time_integrators/crank-nicolson.test........................................................ FAILED (EXODIFF)
time_integrators/implicit-euler.test........................................................ FAILED (EXODIFF)
variables/fe_hier.test_hier_2_1d............................................................ FAILED (EXODIFF)
variables/fe_hier.test_hier_2_3d............................................................ FAILED (EXODIFF)
variables/fe_hier.test_hier_1_2d............................................................ FAILED (EXODIFF)
time_steppers/timesequence_stepper.timesequence_failed....................................... FAILED (ERRMSG)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
postprocessors/variable_residual_norm.variable_residual_test................................ FAILED (CSVDIFF)
time_integrators/scalar.scalar.............................................................. FAILED (CSVDIFF)
vectorpostprocessors/least_squares_fit.least_squares........................................ FAILED (CSVDIFF)
```

Failed tests in modules frmo changing wp to ds:
```
combined/test:contact.pressurePenalty_test.................................................. FAILED (TIMEOUT)
combined/test:gravity.rz_quad4_sm........................................................... FAILED (TIMEOUT)
combined/test:gravity.rz_quad8_sm........................................................... FAILED (TIMEOUT)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.moving................. FAILED (TIMEOUT)
heat_conduction/test:heat_conduction/3d_quadrature_gap_heat_transfer.moving................. FAILED (TIMEOUT)
heat_conduction/test:heat_conduction/3d_quadrature_gap_heat_transfer.second................. FAILED (TIMEOUT)
combined/test:1D_axisymmetric.axisymmetric_gps_finite....................................... FAILED (EXODIFF)
chemical_reactions/test:aqueous_equilibrium.2species........................................ FAILED (EXODIFF)
chemical_reactions/test:aqueous_equilibrium.2species_without_action......................... FAILED (EXODIFF)
combined/test:combined_plasticity_temperature.temp_dep_yield................................ FAILED (EXODIFF)
combined/test:CHSplitFlux.simple_transient_diffusion_flux................................... FAILED (EXODIFF)
combined/test:elastic_patch.rz_nonlinear.................................................... FAILED (EXODIFF)
combined/test:eigenstrain.variable_finite................................................... FAILED (EXODIFF)
combined/test:elastic_thermal_patch.elastic_thermal_jacobian_rz_smp_sm...................... FAILED (EXODIFF)
combined/test:DiffuseCreep.stress_based_chem_pot............................................ FAILED (EXODIFF)
combined/test:gap_heat_transfer_htonly.RSpherical........................................... FAILED (EXODIFF)
combined/test:heat_conduction_patch.test.................................................... FAILED (EXODIFF)
combined/test:heat_conduction_patch.test_hex20.............................................. FAILED (EXODIFF)
combined/test:homogenization.heatConduction_test............................................ FAILED (EXODIFF)
combined/test:inelastic_strain/elas_plas.elastic_plastic_sm................................. FAILED (EXODIFF)
combined/test:inelastic_strain/elas_plas.elastic_plastic_cycle_sm........................... FAILED (EXODIFF)
combined/test:phase_field_fracture.crack2d.................................................. FAILED (EXODIFF)
combined/test:gap_heat_transfer_mortar.test................................................. FAILED (EXODIFF)
combined/test:simple_contact.RZ............................................................. FAILED (EXODIFF)
combined/test:simplest_contact.test......................................................... FAILED (EXODIFF)
combined/test:solid_mechanics/HHT_time_integrator.test...................................... FAILED (EXODIFF)
combined/test:solid_mechanics/pressure.RZ................................................... FAILED (EXODIFF)
combined/test:solid_mechanics/Rayleigh_damping/Newmark_time_integration.test................ FAILED (EXODIFF)
combined/test:solid_mechanics/Time_integration/Newmark_time_integration.test................ FAILED (EXODIFF)
combined/test:phase_field_fracture.void2d_nonsplit.......................................... FAILED (EXODIFF)
combined/test:solid_mechanics/beam_pbp.sm_beam_pbp.......................................... FAILED (EXODIFF)
combined/test:thermo_mech.thermomechanics_SMP............................................... FAILED (EXODIFF)
heat_conduction/test:heat_conduction/coupled_convective_flux.test........................... FAILED (EXODIFF)
combined/test:simple_contact.3D_2........................................................... FAILED (EXODIFF)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.gap_conductivity_property FAILED (EXODIFF)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.perfect................ FAILED (EXODIFF)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.perfectQ8.............. FAILED (EXODIFF)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.perfectQ9.............. FAILED (EXODIFF)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.nonmatching............ FAILED (EXODIFF)
heat_conduction/test:heat_conduction/2d_quadrature_gap_heat_transfer.second_order........... FAILED (EXODIFF)
combined/test:contact.pressureAugLag_test................................................... FAILED (EXODIFF)
phase_field/test:conserved_noise.normal..................................................... FAILED (EXODIFF)
phase_field/test:KKS_system.kks_example..................................................... FAILED (EXODIFF)
combined/test:phase_field_fracture_viscoplastic.crack2d..................................... FAILED (EXODIFF)
heat_conduction/test:heat_conduction/3d_quadrature_gap_heat_transfer.nonmatching............ FAILED (EXODIFF)
solid_mechanics/test:cracking.plane_stress.................................................. FAILED (EXODIFF)
solid_mechanics/test:cracking.cracking_function............................................. FAILED (EXODIFF)
solid_mechanics/test:combined_creep_plasticity.stress_relaxation_sm......................... FAILED (EXODIFF)
solid_mechanics/test:cracking.xyz........................................................... FAILED (EXODIFF)
solid_mechanics/test:cracking.exponential................................................... FAILED (EXODIFF)
solid_mechanics/test:interaction_integral.ii_3d_rot......................................... FAILED (EXODIFF)
solid_mechanics/test:combined_creep_plasticity.stress_prescribed_sm......................... FAILED (EXODIFF)
solid_mechanics/test:cracking.rz_exponential................................................ FAILED (EXODIFF)
tensor_mechanics/test:2D_geometries.plane_strain............................................ FAILED (EXODIFF)
tensor_mechanics/test:2D_geometries.axisym_resid............................................ FAILED (EXODIFF)
tensor_mechanics/test:auxkernels.ranktwoscalaraux........................................... FAILED (EXODIFF)
tensor_mechanics/test:combined_creep_plasticity.stress_relaxation........................... FAILED (EXODIFF)
tensor_mechanics/test:CylindricalRankTwoAux.test............................................ FAILED (EXODIFF)
tensor_mechanics/test:dynamics/prescribed_displacement.displacement_bc...................... FAILED (EXODIFF)
tensor_mechanics/test:dynamics/prescribed_displacement.displacement_bc_gravity.............. FAILED (EXODIFF)
tensor_mechanics/test:dynamics/rayleigh_damping.newmark..................................... FAILED (EXODIFF)
tensor_mechanics/test:dynamics/time_integration.hht......................................... FAILED (EXODIFF)
tensor_mechanics/test:dynamics/time_integration.newmark..................................... FAILED (EXODIFF)
tensor_mechanics/test:finite_strain_elastic.eigen_sol....................................... FAILED (EXODIFF)
tensor_mechanics/test:truss.truss_3d........................................................ FAILED (EXODIFF)
tensor_mechanics/test:combined_creep_plasticity.stress_prescribed........................... FAILED (EXODIFF)
combined/test:power_law_creep.power_law_creep_sm............................................ FAILED (EXODIFF)
tensor_mechanics/test:2D_geometries.finite_planestrain_Bbar................................. FAILED (EXODIFF)
tensor_mechanics/test:elastic_patch.elastic_patch........................................... FAILED (EXODIFF)
xfem/test:init_solution_propagation.init_solution_propagation............................... FAILED (EXODIFF)
solid_mechanics/test:t_stress.3d............................................................ FAILED (EXODIFF)
combined/test:contact_verification/patch_tests/brick_1.frictionless_aug_sm................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/brick_1.frictionless_pen_sm................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/brick_1.mu_0_2_pen_sm......................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/brick_3.frictionless_pen_sm................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/brick_3.frictionless_aug_sm................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/cyl_1.frictionless_aug_sm..................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/cyl_1.frictionless_pen_sm..................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/cyl_1.mu_0_2_pen_sm........................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/plane_1.frictionless_pen_sm................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/ring_1.frictionless_pen_sm.................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/ring_1.mu_0_2_pen_sm.......................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/ring_1.frictionless_aug_sm.................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/plane_1.mu_0_2_pen_sm......................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/single_pnt_2d.mu_0_0_pen_sm................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/ring_3.mu_0_2_pen_sm.......................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/ring_1.glued_pen_sm........................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/ring_3.glued_pen_sm........................... FAILED (ERRMSG)
combined/test:evolving_mass_density.test_shear_z_tm.......................................... FAILED (ERRMSG)
combined/test:evolving_mass_density.test_z_tm................................................ FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/plane_3.mu_0_2_pen_sm......................... FAILED (ERRMSG)
combined/test:incremental_slip.slip_tm....................................................... FAILED (ERRMSG)
combined/test:nodal_area.Hex20_3............................................................. FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/cyl_4.mu_0_2_pen_sm........................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/brick_3.mu_0_2_pen_sm......................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/brick_4.mu_0_2_pen_sm......................... FAILED (ERRMSG)
combined/test:contact_verification/overclosure_removal.test.................................. FAILED (ERRMSG)
combined/test:thermal_elastic.test........................................................... FAILED (ERRMSG)
solid_mechanics/test:rate_dep_smear_crack.iso_test........................................... FAILED (ERRMSG)
tensor_mechanics/test:smeared_cracking.single_element........................................ FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/single_pnt_2d.mu_0_2_kin...................... FAILED (ERRMSG)
combined/test:elastic_patch.rz_sm............................................................ FAILED (ERRMSG)
tensor_mechanics/test:elastic_patch.elastic_patch_quadratic.................................. FAILED (ERRMSG)
combined/test:internal_volume.hex8........................................................... FAILED (ERRMSG)
combined/test:contact_verification/patch_tests/cyl_1.glued_pen_sm........................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/cyl_2.mu_0_2_pen_sm.......................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/cyl_2.glued_pen_sm........................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/cyl_2.frictionless_aug_sm.................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/cyl_3.mu_0_2_pen_sm.......................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/cyl_4.frictionless_aug_sm.................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/cyl_4.frictionless_pen_sm.................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/plane_1.frictionless_aug_sm.................. FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/brick_4.frictionless_pen_sm.................. FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/ring_2.glued_kin_sm.......................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/single_pnt_2d.glued_pen_sm................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/single_pnt_2d.frictionless_pen_sm............ FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/single_pnt_2d.mu_0_0_kin_sm.................. FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/single_pnt_2d.glued_kin_sm................... FAILED (CSVDIFF)
combined/test:contact_verification/patch_tests/single_pnt_2d.frictionless_kin_sm............ FAILED (CSVDIFF)
phase_field/test:conserved_noise.integral_normal_masked..................................... FAILED (CSVDIFF)
combined/test:internal_volume.hex20......................................................... FAILED (CSVDIFF)
tensor_mechanics/test:orthotropic_plasticity.test........................................... FAILED (CSVDIFF)
solid_mechanics/test:j_integral_vtest.J_ellip_cm............................................ FAILED (CSVDIFF)
solid_mechanics/test:j_integral_vtest.j_ellip_cfp........................................... FAILED (CSVDIFF)
solid_mechanics/test:j_integral_vtest.j_ellip............................................... FAILED (CSVDIFF)
solid_mechanics/test:j_integral_vtest.J_ellip_cm_cfp........................................ FAILED (CSVDIFF)
-------------------------------------------------------------------------------------------------------------
Ran 1968 tests in 334.8 seconds
1842 passed, 215 skipped, 126 FAILED
```

# 9/11/17

Back-trace for `ExplicitEuler::preSolve`
```
frame #0: 0x00000001019c7ad0 libmoose-dbg.0.dylib`ExplicitEuler::preSolve(this=0x000000010e414c80) at ExplicitEuler.C:35
frame #1: 0x00000001011b1795 libmoose-dbg.0.dylib`NonlinearSystemBase::onTimestepBegin(this=0x000000010e826618) at NonlinearSystemBase.C:592
frame #2: 0x0000000100d4b8b7 libmoose-dbg.0.dylib`FEProblemBase::onTimestepBegin(this=0x000000010e824800) at FEProblemBase.C:3919
frame #3: 0x00000001013a2fd7 libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010e410430, input_dt=-1) at Transient.C:401
frame #4: 0x00000001013a2edd libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010e410430, input_dt=-1) at Transient.C:377
frame #5: 0x00000001013a2b83 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010e410430) at Transient.C:298
frame #6: 0x0000000101114a98 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e81f600) at MooseApp.C:561
frame #7: 0x00000001011171fa libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e81f600) at MooseApp.C:717
frame #8: 0x0000000100002983 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe828) at main.C:36
```

Bt for `ExplicitEuler::computeTimeDerivatives`
```
frame #0: 0x00000001019c7b70 libmoose-dbg.0.dylib`ExplicitEuler::computeTimeDerivatives(this=0x000000010e414c80) at ExplicitEuler.C:44
frame #1: 0x00000001011b4b83 libmoose-dbg.0.dylib`NonlinearSystemBase::computeTimeDerivatives(this=0x000000010e826618) at NonlinearSystemBase.C:701
frame #2: 0x0000000100d4cca2 libmoose-dbg.0.dylib`FEProblemBase::computeResidualType(this=0x000000010e824800, soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8, type=KT_ALL) at FEProblemBase.C:4035
frame #3: 0x0000000100d4bfd0 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010e824800, soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8) at FEProblemBase.C:3971
frame #4: 0x0000000100d4bf86 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010e824800, (null)=0x000000010e090210, soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8) at FEProblemBase.C:3963
frame #5: 0x000000010119971b libmoose-dbg.0.dylib`Moose::compute_residual(soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8, sys=0x000000010e090210) at NonlinearSystem.C:44
frame #6: 0x000000010473647d libmesh_dbg.0.dylib`::__libmesh_petsc_snes_residual(snes=0x0000000110803020, x=0x0000000110007220, r=0x0000000110009420, ctx=0x000000010e090e00) at petsc_nonlinear_solver.C:130
frame #7: 0x000000010bf0e09f libpetsc.3.7.dylib`SNESComputeFunction(snes=0x0000000110803020, x=0x0000000110007220, y=0x0000000110009420) at snes.c:2145 [opt]
frame #8: 0x000000010bf3f96d libpetsc.3.7.dylib`SNESSolve_KSPONLY(snes=0x0000000110803020) at ksponly.c:25 [opt]
frame #9: 0x000000010bf1326a libpetsc.3.7.dylib`SNESSolve(snes=<unavailable>, b=<unavailable>, x=0x0000000110007220) at snes.c:4005 [opt]
frame #10: 0x000000010473393c libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x000000010e090e00, jac_in=0x000000010e090d30, x_in=0x000000010e0905a0, r_in=0x000000010e090b50, (null)=0.00000001, (null)=10000) at petsc_nonlinear_solver.C:702
frame #11: 0x0000000104850eb0 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x000000010e090210) at nonlinear_implicit_system.C:181
frame #12: 0x00000001019ceffe libmoose-dbg.0.dylib`TimeIntegrator::solve(this=0x000000010e414c80) at TimeIntegrator.C:53
frame #13: 0x000000010119b2ac libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000010e826618) at NonlinearSystem.C:159
frame #14: 0x0000000100d47e02 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x000000010e824800) at FEProblemBase.C:3749
frame #15: 0x00000001019edf64 libmoose-dbg.0.dylib`TimeStepper::step(this=0x000000010e211d50) at TimeStepper.C:160
frame #16: 0x00000001013a32e2 libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010e410430, input_dt=-1) at Transient.C:444
frame #17: 0x00000001013a2edd libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010e410430, input_dt=-1) at Transient.C:377
frame #18: 0x00000001013a2b83 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010e410430) at Transient.C:298
frame #19: 0x0000000101114a98 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e81f600) at MooseApp.C:561
frame #20: 0x00000001011171fa libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e81f600) at MooseApp.C:717
frame #21: 0x0000000100002983 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe828) at main.C:36
```

Bt for `TimeDerivative::computeQpResidual`
```
frame #0: 0x00000001014e4db0 libmoose-dbg.0.dylib`TimeDerivative::computeQpResidual(this=0x0000000110024600) at TimeDerivative.C:42
frame #1: 0x00000001014e5495 libmoose-dbg.0.dylib`TimeKernel::computeResidual(this=0x0000000110024600) at TimeKernel.C:44
frame #2: 0x0000000100caa9f4 libmoose-dbg.0.dylib`ComputeResidualThread::onElement(this=0x00007fff5fbfbda0, elem=0x000000010e113fc0) at ComputeResidualThread.C:124
frame #3: 0x0000000100b2c5ad libmoose-dbg.0.dylib`ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >::operator(this=0x00007fff5fbfbda0, range=0x000000010e092790, bypass_threading=false)(libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> const&, bool) at ThreadedElementLoopBase.h:202
frame #4: 0x00000001011b940b libmoose-dbg.0.dylib`void libMesh::Threads::parallel_reduce<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*>, ComputeResidualThread>(range=0x000000010e092790, body=0x00007fff5fbfbda0) at threads_tbb.h:164
frame #5: 0x00000001011b0421 libmoose-dbg.0.dylib`NonlinearSystemBase::computeResidualInternal(this=0x000000010e826618, type=KT_ALL) at NonlinearSystemBase.C:1109
frame #6: 0x00000001011aff09 libmoose-dbg.0.dylib`NonlinearSystemBase::computeResidual(this=0x000000010e826618, residual=0x00007fff5fbfdad8, type=KT_ALL) at NonlinearSystemBase.C:551
frame #7: 0x0000000100d4ce56 libmoose-dbg.0.dylib`FEProblemBase::computeResidualType(this=0x000000010e824800, soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8, type=KT_ALL) at FEProblemBase.C:4062
frame #8: 0x0000000100d4bfd0 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010e824800, soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8) at FEProblemBase.C:3971
frame #9: 0x0000000100d4bf86 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010e824800, (null)=0x000000010e090210, soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8) at FEProblemBase.C:3963
frame #10: 0x000000010119971b libmoose-dbg.0.dylib`Moose::compute_residual(soln=0x000000010e0906c0, residual=0x00007fff5fbfdad8, sys=0x000000010e090210) at NonlinearSystem.C:44
frame #11: 0x000000010473647d libmesh_dbg.0.dylib`::__libmesh_petsc_snes_residual(snes=0x0000000110803020, x=0x0000000110007220, r=0x0000000110009420, ctx=0x000000010e090e00) at petsc_nonlinear_solver.C:130
frame #12: 0x000000010bf0e09f libpetsc.3.7.dylib`SNESComputeFunction(snes=0x0000000110803020, x=0x0000000110007220, y=0x0000000110009420) at snes.c:2145 [opt]
frame #13: 0x000000010bf3f96d libpetsc.3.7.dylib`SNESSolve_KSPONLY(snes=0x0000000110803020) at ksponly.c:25 [opt]
frame #14: 0x000000010bf1326a libpetsc.3.7.dylib`SNESSolve(snes=<unavailable>, b=<unavailable>, x=0x0000000110007220) at snes.c:4005 [opt]
frame #15: 0x000000010473393c libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x000000010e090e00, jac_in=0x000000010e090d30, x_in=0x000000010e0905a0, r_in=0x000000010e090b50, (null)=0.00000001, (null)=10000) at petsc_nonlinear_solver.C:702
frame #16: 0x0000000104850eb0 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x000000010e090210) at nonlinear_implicit_system.C:181
frame #17: 0x00000001019ceffe libmoose-dbg.0.dylib`TimeIntegrator::solve(this=0x000000010e414c80) at TimeIntegrator.C:53
frame #18: 0x000000010119b2ac libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000010e826618) at NonlinearSystem.C:159
frame #19: 0x0000000100d47e02 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x000000010e824800) at FEProblemBase.C:3749
frame #20: 0x00000001019edf64 libmoose-dbg.0.dylib`TimeStepper::step(this=0x000000010e211d50) at TimeStepper.C:160
frame #21: 0x00000001013a32e2 libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010e410430, input_dt=-1) at Transient.C:444
frame #22: 0x00000001013a2edd libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010e410430, input_dt=-1) at Transient.C:377
frame #23: 0x00000001013a2b83 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010e410430) at Transient.C:298
frame #24: 0x0000000101114a98 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e81f600) at MooseApp.C:561
frame #25: 0x00000001011171fa libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e81f600) at MooseApp.C:717
frame #26: 0x0000000100002983 moose_test-dbg`main(argc=3,
argv=0x00007fff5fbfe828) at main.C:36
```

```
  * frame #0: 0x00000001014e4e30 libmoose-dbg.0.dylib`TimeDerivative::computeQpJacobian(this=0x0000000110024600) at TimeDerivative.C:48
    frame #1: 0x00000001014e5035 libmoose-dbg.0.dylib`TimeDerivative::computeJacobian(this=0x0000000110024600) at TimeDerivative.C:61
    frame #2: 0x0000000100c8a0a3 libmoose-dbg.0.dylib`ComputeJacobianThread::computeJacobian(this=0x00007fff5fbfc000) at ComputeJacobianThread.C:100
    frame #3: 0x0000000100c8d45e libmoose-dbg.0.dylib`ComputeJacobianThread::onElement(this=0x00007fff5fbfc000, elem=0x000000010e113fc0) at ComputeJacobianThread.C:204
    frame #4: 0x0000000100b2c5ad libmoose-dbg.0.dylib`ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >::operator(this=0x00007fff5fbfc000, range=0x000000010e092790, bypass_threading=false)(libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> const&, bool) at ThreadedElementLoopBase.h:202
    frame #5: 0x00000001011ca97b libmoose-dbg.0.dylib`void libMesh::Threads::parallel_reduce<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*>, ComputeJacobianThread>(range=0x000000010e092790, body=0x00007fff5fbfc000) at threads_tbb.h:164
    frame #6: 0x00000001011c753a libmoose-dbg.0.dylib`NonlinearSystemBase::computeJacobianInternal(this=0x000000010e826618, jacobian=0x00007fff5fbfd230, kernel_type=KT_ALL) at NonlinearSystemBase.C:1877
    frame #7: 0x00000001011cc60a libmoose-dbg.0.dylib`NonlinearSystemBase::computeJacobian(this=0x000000010e826618, jacobian=0x00007fff5fbfd230, kernel_type=KT_ALL) at NonlinearSystemBase.C:2109
    frame #8: 0x0000000100d4d8c2 libmoose-dbg.0.dylib`FEProblemBase::computeJacobian(this=0x000000010e824800, soln=0x000000010e0906c0, jacobian=0x00007fff5fbfd230, kernel_type=KT_ALL) at FEProblemBase.C:4121
    frame #9: 0x0000000100d4d16e libmoose-dbg.0.dylib`FEProblemBase::computeJacobian(this=0x000000010e824800, (null)=0x000000010e090210, soln=0x000000010e0906c0, jacobian=0x00007fff5fbfd230) at FEProblemBase.C:4070
    frame #10: 0x000000010119959b libmoose-dbg.0.dylib`Moose::compute_jacobian(soln=0x000000010e0906c0, jacobian=0x00007fff5fbfd230, sys=0x000000010e090210) at NonlinearSystem.C:34
    frame #11: 0x0000000104738733 libmesh_dbg.0.dylib`::__libmesh_petsc_snes_jacobian(snes=0x0000000110803020, x=0x0000000110007220, jac=0x0000000110018620, pc=0x0000000110018620, ctx=0x000000010e090e00) at petsc_nonlinear_solver.C:221
    frame #12: 0x000000010bf0eb5e libpetsc.3.7.dylib`SNESComputeJacobian(snes=0x0000000110803020, X=0x0000000110007220, A=0x0000000110018620, B=<unavailable>) at snes.c:2312 [opt]
    frame #13: 0x000000010bf3f9f2 libpetsc.3.7.dylib`SNESSolve_KSPONLY(snes=0x0000000110803020) at ksponly.c:38 [opt]
    frame #14: 0x000000010bf1326a libpetsc.3.7.dylib`SNESSolve(snes=<unavailable>, b=<unavailable>, x=0x0000000110007220) at snes.c:4005 [opt]
    frame #15: 0x000000010473393c libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x000000010e090e00, jac_in=0x000000010e090d30, x_in=0x000000010e0905a0, r_in=0x000000010e090b50, (null)=0.00000001, (null)=10000) at petsc_nonlinear_solver.C:702
    frame #16: 0x0000000104850eb0 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x000000010e090210) at nonlinear_implicit_system.C:181
    frame #17: 0x00000001019ceffe libmoose-dbg.0.dylib`TimeIntegrator::solve(this=0x000000010e414c80) at TimeIntegrator.C:53
    frame #18: 0x000000010119b2ac libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000010e826618) at NonlinearSystem.C:159
    frame #19: 0x0000000100d47e02 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x000000010e824800) at FEProblemBase.C:3749
    frame #20: 0x00000001019edf64 libmoose-dbg.0.dylib`TimeStepper::step(this=0x000000010e211d50) at TimeStepper.C:160
    frame #21: 0x00000001013a32e2 libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010e410430, input_dt=-1) at Transient.C:444
    frame #22: 0x00000001013a2edd libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010e410430, input_dt=-1) at Transient.C:377
    frame #23: 0x00000001013a2b83 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010e410430) at Transient.C:298
    frame #24: 0x0000000101114a98 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e81f600) at MooseApp.C:561
    frame #25: 0x00000001011171fa libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e81f600) at MooseApp.C:717
    frame #26: 0x0000000100002983 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe828) at main.C:36
    frame #27: 0x00007fffc5084235 libdyld.dylib`start + 1
    frame #28: 0x00007fffc5084235 libdyld.dylib`start + 1
```

# 9/12/17

On test problem (`selective_reinit_test.i`) without any source code mods, `MooseVariable::computeElemValues()`, called
through the aux system for constant monomial takes 111 ms which is 5% of
execution time. 89 ms is spent within the actual method itself

This line:
```c++
using buildPtr = MooseObjectPtr (*)(const InputParameters & parameters);
```
means that `buildPtr` will be a pointer to "something" that takes an argument
`const InputParameters & parameters` and returns a `MooseObjectPtr`.

# 9/13/17

Using the MooseVariable computeElemValues code, computeElemValues occupies
1256ms of the total 7746ms simulation time, which is 16% of the total simulation
time.

Ok, now using the constant monomial code, computeElemValues occupies 624 ms of
the total 5766 ms simulation time, which is 11% of the total simulation time. So
essentially this method halves the computation time.

3d test:

- Old code:
    - Total time: 58.985 s
    - time in `MooseVariable::computeElemValues`: 16.406 s
    - % = 27

- First commit:
    - Total time: 38.674 s
    - time in `MooseVariableConstMonomial::computeElemValues`:  5.474 s
    - % = 14

- Second commit
    - Total time: 43.498 s
	- time in `MooseVariableConstMonomial::computeElemValues`: 8.089 s
	- % = 18.6

- Commenting out the zeroing
    - Total time: 38.107 s
	- time in `compute`: 4.483 s
	- % = 11.8

# 9/19/17

Alright, the mesh `read` method is implemented in
`libMesh::UnstructuredMesh`. So it should be equally usable by `DistributedMesh`
in theory.

# 9/26/17

Command to compile a solo file depending on MOOSE:

```
clang++ -std=c++11 -I $HOME/projects/moose/libmesh/installed/include -I $HOME/projects/moose/framework/include/base -I $HOME/projects/moose/framework/include/utils -L $HOME/projects/moose/libmesh/installed/lib -L $HOME/projects/moose/framework -L $HOME/projects/moose/framework/contrib/hit monotone_tester.cpp -lmoose-dbg -lmesh_dbg -lhit-dbg
```

# 10/16/17

The const_cast in reinitFE should be determined by the FEType family

# 10/25/17

The execution order for multiapps:

```
  // Execute Transfers _to_ MultiApps
  // Execute MultiApps
  // Execute Transfers _from_ MultiApps
```

Ok so OpenMC generates the powers and these are mapped over to Cerberus using
Functional Expansion Tallies which are expressed as a summation of products of
Zernike and Legendre polynomials. Cerberus then generates the temperature
distribution and that temperature distribution is mapped over to OpenMC also
using Zernike and Legendre polynomials I believe.

polytemps is set from either `set_fuel_temp` or `set_distrib_temp`. The latter
gets called from the OpenMC MultiApp solve. The former is declared to be an old
version of the latter in comments. The only other place where `polytemps` is
used is in `tracking.F90`.

So Matt built code into the particle transporting algorithm that allows setting
the temperature of the material based on evaluation of Zernike temperature
polynomials coming from Cerberus. I'm guessing that we probably don't want to
drill into the code of every Monte Carlo package and make changes. Review with
Mark what we want to pass between codes.

# 10/26/17

Ok, we have the following members in `AssemblyBase`:
```
  std::map<FEType, FEShapeData *> _fe_shape_data;
  /// Each dimension's actual fe objects indexed on type
  std::map<unsigned int, std::map<FEType, FEAbstract *>> _fe;
```

# 10/30/17

Assembly::copyShapes <- FEProblemBase::prepareShapes <-
ComputeJacobianThread::computeJacobian
Assembly::copyShapes <- FEProblemBase::prepareShapes <- ComputeFullJacobianThread::computeJacobian

ComputeJacobianThread -> FEProblemBase::prepareNeighborShapes -> Assembly::copyNeighborShapes(unsigned int var) -> MooseVariableField::secondPhiFaceNeighbor() ->
Assembly::feSecondPhiFaceNeighbor(FEType) -> Currently a templated function

# 10/31/17

ComputeResidualThread::onElement(const Elem *) -> FEProblemBase::reinitElem(const Elem *, THREAD_ID) -> AuxiliarySystem::reinitElem(const Elem *, THREAD_ID) -> MooseVariableField::reinitAux()

When we're computing our residual thread, the effective execution order is:

```
FEProblemBase::prepare(const Elem *)
  _assembly[tid]->reinit(elem);
    _current_elem = elem;
    reinitFE(elem);
_nl->prepare(tid);
    var->clearDofIndices()
	var->prepare()
FEProblemBase::reinitElem
  _nl->reinitElem(elem, tid);
     var->computeElemValues();
  _aux->reinitElem(elem, tid);
    For nodal_vars:
    var->computeElemValues

    For elemental_vars:
    var->reinitAux
    var->computeElemValues

FEProblemBase::reinitMaterials
Kernel::computeResidual
```

MooseVariable::_elem is a pointer reference to Assembly::_current_elem which is
set within Assembly::reinit.

Within FEProblemBase::prepare:

```
Assembly::reinit
	Assembly::reinitFE
	    FEShapeData::_phi::shallowCopy(FEBase::get_phi())
NonlinearSystemBase::prepare
AuxiliarySystem::prepare
Assembly::prepare
```

```
ComputeNodalStuffThreads:
  FEProblemBase::reinitNode
    *System*::reinitNode
	  if (var->isNodal())
      {
        var->reinitNode();
        var->computeNodalValues();
      }
```

In AuxiliarySystem::reinitElem, we have:

For nodal_vars:
var->computeElemValues

For elemental_vars:
var->reinitAux
var->computeElemValues

So reinitAux is purely for elemental variables.

Any non-zero value in C is true. Zero is false.

_phi_face_neighbor is initialized from _assembly.fePhiFaceNeighbor(_fe_type)

fePhiFaceNeighbor first calls buildFaceNeighborFE with type and then returns `_fe_shape_data_face_neighbor[type]->_phi`. buildFaceNeighborFE
builds an element of correct dim and type and stores them in _fe_face_neighbor
but when does `_fe_shape_data_face_neighbor[type]` actually get its phi??
```
  std::map<FEType, FEShapeData *> _fe_shape_data;
  std::map<FEType, FEShapeData *> _fe_shape_data_face;
  std::map<FEType, FEShapeData *> _fe_shape_data_neighbor;
  std::map<FEType, FEShapeData *> _fe_shape_data_face_neighbor;

  /// types of finite elements
  std::map<unsigned int, std::map<FEType, FEBase *>> _fe_neighbor;
  std::map<unsigned int, std::map<FEType, FEBase *>> _fe_face_neighbor;
  /// Each dimension's helper objects
  std::map<unsigned int, FEBase **> _holder_fe_neighbor_helper;
  std::map<unsigned int, FEBase **> _holder_fe_face_neighbor_helper;
```

Truly nodal methods should have no place in `MooseVariableField`. The trick is
to determine whether the method is associated with the geometric entity known as
a node (I'm determined that we are going to establish this concept), or whether
the method is associated with the nonlinear system entity known as a dof.

For instance, when var->nodalSln() is called from Coupleable::coupledValue, that
is a method call that is geometric node motivated. And as it stands now, how the
method is implemented is ideal because if someone tried to couple in a variable
to a NodalKernel (etc.) whose solution is not continuous at nodes, then they
should get an error. Ok so I like keeping the names of those methods as they
are.

So the question for tomorrow is methods like computeNodalValues. That name as it
stands I think should not be in MooseVariableField.

# 11/1/17

Ok, `computeNodalValues` is an appropriate name since its ultimately called from
nodal object compute threads. `computeNodalNeighborValues` is also appropriately
named since it is called for nodal constraints.

```
NonlinearSystemBase::computeResidualInternal
  NonlinearSystemBase::enforceNodalConstraintsResidual
    FEProblemBaes::reinitNodesNeighbor
      SystemBase::reinitNodesNeighbor
        var->reinitNodesNeighbor(nodes);
        var->computeNodalNeighborValues();
      auxiliarySystem::reinitNodesNeighbor
```

So the neighbor node is needed for constraints I guess.

Ok, `setNodalValue` is used for both nodal and elemental variables...thus not an
appropriate name!

`_has_nodal_value` is only used in `add` and `insert`. And it only refers to the
solution vector which is an array of numbers. Thus it should really read:
`_has_dof_value`

Legitimate nodal methods:
```
isNodal
node
nodalDofIndex
isNodalDefined
nodeNeighbor
nodalDofIndexNeighbor
isNodalNeighborDefined

nodalValue
nodalValueOld
nodalValueOlder
nodalValuePreviousNL
nodalValueDot
nodalValueDuDotDu
nodalValueNeighbor
nodalValueOldNeighbor
nodalValueOlderNeighbor
nodalValuePreviousNLNeighbor
nodalValueDotNeighbor
nodalValueDuDotDuNeighbor

computeNodalValues
computeNodalNeighborValues
setNodalValue
reinitNode
reinitNodeNeighbor
reinitNodes
reinitNodesNeighbor
getNodalValue
getNodalValueOld
getNodalValueOlder
```

Legitimate nodal data members:
```
_has_nodal_dof
_need_nodal_u*
```

In computeElemValues (computeValuesHelper), `_nodal_u` is only computed if the
`_need_nodal_u` flag is true (that's the only place that flag is checked before
`_nodal_u` computations are done). The only way that flag sets set to true is if
`nodalValue` gets called.

So I was thinking I might entirely drop the nodal bits from
`computeElemValues`. However, this might break nodal (e.g. Lagrange) aux
variables. Actually think for a kernel object...say a kernel depends on a nodal
aux variable. What we still want are the interpolations of the nodal aux
variable at quadrature points, which is what is done in computeElemValues using
the non-nodal logic. The nodal logic is not necessary for this case.

```
AuxiliarySystem::compute
  AuxiliarySystem::computeScalarVars
  AuxiliarySystem::computeNodalVars
    ComputeNodalAuxVarsThread::()
	  ::pre
	  for (node : nodes)
	    ::onNode(node)
		  var->prepareAux
		  _fe_problem.reinitNode
		  aux->compute() # This is an aux kernel
		  var->insert(_aux_sys.solution())
		::postNode(node)
      ::post
  AuxiliarySystem::computeElementalVars
```

How is `_dof_indices` correct before insert is called in the above code scheme?
It's done through FEProblem::reinitNode which ends up calling
MooseVariable::reinitNode, which correctly sets `_dof_indices`. Now what does
insert actually do? Insert inserts solution values for the current variable at
the current node (in this case) using the correct dof index into the global
auxiliary solution vector (which in general is a PetscVector (which is a class
defined in libmesh)).

AFAICT there really is no need to do nodal calculations in
computeElemValues. We're axing it!!

# 11/16/17

Apparently NonlinearSystem has both a `_transient_sys` member which is
initialized in the `NonlinaerSystem` constructor and a `_sys` member
which is initialized in the `NonlinearSystemBase` constructor.

So the FEProblem is created with the action `CreateProblemAction`. It's created
using the `Factory`. In the `FEProblemBase` constructor the `EquationSystems`
object is created/initialized. In the `FEProblem` constructor, the
`NonLinearSystem` object is created using `std::make_shared`. In the
`NonlinearSystem` constructor, the `TransientNonlinearImplicitSystem` object is
created, so we now have all of our important problem objects:

```
FEProblem (moose)
EquationSystems (libmesh)
NonlinearSystem (moose)
TransientNonlinearImplicitSystem (libmesh)
```

Notes:
The `TransientNonlinearImplicitSystem` object is accessed using the `system()`
method of `NonlinearSystem`. `system()` returns a reference to an object of type
`System` so note that we have implicitly downcasted
`TransientNonlinearImplicitSystem` to `System`.

MOOSE classes:
```
FEProblem

SystemBase
NonlinearSystemBase
NonlinearSystem
```

LibMesh classes:
```
System
ImplicitSystem
NonlinearImplicitSystem
typedef Transient<NonlinearImplicitSystem> TransientNonlinearImplicitSystem

EquationSystems
```

`TimeIntegrator::solve()` calls `System::solve()` which in turn calls
`NonlinearSolver::solve()` which takes residual and jacobian arguments. In our case
we actually call `PetscNonlinaerSolver::solve()` because its the derivative
class of `NonlinearSolver` that we actually use. The `PetscNonlinearSolver`
object is created in the constructor for `NonlinearImplicitSystem` using the
static `NonlinaerSolver::build` method.

The `ExplicitSystem` constructor creates the public `rhs` data member. The
`ImplicitSystem` constructor (`ImplicitSystem` is a direct descendent of
`ExplicitSystem`) creates the public `matrix` data member. The `System` constructor
creates the public `solution` data member.

A couple of notes:

- David made a couple of commits that created optimizations for explicit time
  schemes. These can be found with the hashes: `301a7d9e3' and `794c9730f0b8`.
- Importantly these commits did not create the enum `KernelType`. It just made
  use of it.

I think I'm just going to ask John and David how different `KernelType`s are
different. My question is are they treated any different for implicit schemes?
Or is the distinction only relevant for explicit schemes?

# 11/17/17

`TimeIntegrator::postStep` is called from `NonlinearSystemBase::computeResidual`
whose bt is:
```
TimeIntegrator::postStep
NonlinearSystemBase::computeResidual
FEProblemBase::computeResidualType
FEProblemBase::computeResidual
FEProblemBase::computeResidual
Moose::compute_residual
::__libmesh_petsc_snes_residual
SNESComputeFunction
SNESSolve_KSPONLY
SNESSolve
libMesh::PetscNonlinearSolver<double>::solve
libMesh::NonlinearImplicitSystem::solve
TimeIntegrator::solve
NonlinearSystem::solve
FEProblemBase::solve
TimeStepper::step
Transient::solveStep
Transient::takeStep
Transient::execute
MooseApp::executeExecutioner
MooseApp::run
```

`TimeIntegrator::postSolve` is called this way, immediately after
`TimeIntegrator::solve`:
```
TimeIntegrator::postSolve
NonlinearSystem::solve
FEProblemBase::solve
TimeStepper::step
Transient::solveStep
Transient::takeStep
Transient::execute
MooseApp::executeExecutioner
MooseApp::run
```

```
TimeIntegrator::preSolve
NonlinearSystemBase::onTimestepBegin
FEProblemBase::onTimestepBegin
Transient::solveStep
Transient::takeStep
Transient::execute
MooseApp::executeExecutioner
MooseApp::run
```

What methods do we have (for TimeIntegrator)?:
```
preSolve
preStep
solve
postStep
postSolve
```

```
Transient::takeStep
  Transient::solveStep
    FEProblemBase::onTimestepBegin
	  NonlinaerSystemBase::onTimestepBegin
	    TimeIntegrator::preSolve
    TimeStepper::step->FEProblemBase::solve->NonlinearSystem::solve
	  TimeIntegrator::solve
	    libMesh::NonlinearImplicitSystem::solve->libMesh::PetscNonlinearSolver::solve->SNESSolve->
	        SNESSolve_KSPONLY->SNESComputeFunction->::__libmesh_petsc_snes_residual->
		    Moose::compute_residual->FEProblemBase::computeResidual->FEProblemBase::computeResidual
          FEProblemBase::computeResidualType
            NonlinearSystemBase::computeTimeDerivatives
              TimeIntegrator::preStep
              TimeIntegrator::computeTimeDerivatives
            NonlinearSystemBase::computeResidual
              NonlinearSystemBase::computeResidualInternal
          	TimeIntegrator::postStep
          	NonlinearSystemBase::computeNodalBCs
      TimeIntegrator::postSolve
```

Kincaid & Cheney: "Another advantage of iterative methods is that they are
usually stable, and they will actually dampen errors (due to roundoff or minor
blunders) as the process continues."

# 11/28/17

Need to figure out why time step is lagging in multi-app.

So incrementStepOrReject is called from:

```
Transient::incrementStepOrReject
Transient::execute
```

```
Transient::incrementStepOrReject
TransientMultiApp::solveStep
```

```
Transient::incrementStepOrReject
TransientMultiApp::advanceStep
FEProblemBase::advanceMultiApps
Transient::incrementStepOrReject
```

The computeStep method of the `sub_app` time stepper is called before its
Transient executioner's `incrementStepOrReject` method is called.

The back traces:
```
* thread #1, queue = 'com.apple.main-thread', stop reason = breakpoint 4.1
  * frame #0: 0x00000001013f84b0 libmoose-dbg.0.dylib`Transient::incrementStepOrReject(this=0x000000010ee7a990) at Transient.C:352
    frame #1: 0x00000001013f8377 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010ee7a990) at Transient.C:317
    frame #2: 0x000000010115fb33 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000011080b200) at MooseApp.C:571
    frame #3: 0x000000010116229a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000011080b200) at MooseApp.C:727
    frame #4: 0x0000000100002383 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe3a8) at main.C:35
    frame #5: 0x00007fffd4cb1235 libdyld.dylib`start + 1
    frame #6: 0x00007fffd4cb1235 libdyld.dylib`start + 1

* thread #1, queue = 'com.apple.main-thread', stop reason = breakpoint 2.1
  * frame #0: 0x0000000101aab690 libmoose-dbg.0.dylib`TimeStepper::computeStep(this=0x00000001104210e0) at TimeStepper.C:75
    frame #1: 0x00000001013f8497 libmoose-dbg.0.dylib`Transient::computeDT(this=0x0000000110455770) at Transient.C:346
    frame #2: 0x0000000101745e16 libmoose-dbg.0.dylib`TransientMultiApp::computeDT(this=0x000000010ee7a640) at TransientMultiApp.C:465
    frame #3: 0x0000000100d86845 libmoose-dbg.0.dylib`FEProblemBase::computeMultiAppsDT(this=0x000000010f060600, type=EXEC_TIMESTEP_BEGIN) at FEProblemBase.C:3349
    frame #4: 0x00000001013fb719 libmoose-dbg.0.dylib`Transient::computeConstrainedDT(this=0x000000010ee7a990) at Transient.C:668
    frame #5: 0x00000001013f87de libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010ee7a990, input_dt=-1) at Transient.C:423
    frame #6: 0x00000001013f872d libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010ee7a990, input_dt=-1) at Transient.C:405
    frame #7: 0x00000001013f83d3 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010ee7a990) at Transient.C:326
    frame #8: 0x000000010115fb33 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000011080b200) at MooseApp.C:571
    frame #9: 0x000000010116229a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000011080b200) at MooseApp.C:727
    frame #10: 0x0000000100002383 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe3a8) at main.C:35
    frame #11: 0x00007fffd4cb1235 libdyld.dylib`start + 1
    frame #12: 0x00007fffd4cb1235 libdyld.dylib`start + 1

* thread #1, queue = 'com.apple.main-thread', stop reason = breakpoint 4.1
  * frame #0: 0x00000001013f84b0 libmoose-dbg.0.dylib`Transient::incrementStepOrReject(this=0x0000000110455770) at Transient.C:352
    frame #1: 0x0000000101744f69 libmoose-dbg.0.dylib`TransientMultiApp::solveStep(this=0x000000010ee7a640, dt=0.001, target_time=0.002, auto_advance=true) at TransientMultiApp.C:344
    frame #2: 0x0000000100d46f60 libmoose-dbg.0.dylib`FEProblemBase::execMultiApps(this=0x000000010f060600, type=EXEC_TIMESTEP_BEGIN, auto_advance=true) at FEProblemBase.C:3247
    frame #3: 0x00000001013f89bf libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010ee7a990, input_dt=-1) at Transient.C:446
    frame #4: 0x00000001013f872d libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010ee7a990, input_dt=-1) at Transient.C:405
    frame #5: 0x00000001013f83d3 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010ee7a990) at Transient.C:326
    frame #6: 0x000000010115fb33 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000011080b200) at MooseApp.C:571
    frame #7: 0x000000010116229a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000011080b200) at MooseApp.C:727
    frame #8: 0x0000000100002383 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe3a8) at main.C:35
    frame #9: 0x00007fffd4cb1235 libdyld.dylib`start + 1
    frame #10: 0x00007fffd4cb1235 libdyld.dylib`start + 1
```

# 11/29/17

The master app reads from `master_out_cp/0005.rd`. It looks for the parameter:
"MOOSE Problem/FEProblemBase/MOOSE Problem/time".

The `SubApp Backup` is also written to `0005.rd`

So the master restartable data holds the key:
"MOOSE Problem/FEProblemBase/MOOSE Problem/time"

as well as the key:
"MOOSE Problem/MultiApps/sub_app/backups"

which I believe itself holds its own restartable data!

See here's a backtrace:
```
  * frame #0: 0x0000000101a56725 libmoose-dbg.0.dylib`RestartableDataIO::serializeRestartableData(this=0x00007fff5fbfb220, restartable_data=size=30, stream=0x0000000112900500) at RestartableDataIO.C:66
    frame #1: 0x0000000101a5b916 libmoose-dbg.0.dylib`RestartableDataIO::createBackup(this=0x00007fff5fbfb220) at RestartableDataIO.C:311
    frame #2: 0x0000000101161ea8 libmoose-dbg.0.dylib`MooseApp::backup(this=0x0000000111801000) at MooseApp.C:661
    frame #3: 0x0000000101729d2e libmoose-dbg.0.dylib`MultiApp::backup(this=0x0000000110365740) at MultiApp.C:342
    frame #4: 0x0000000101734569 libmoose-dbg.0.dylib`void dataStore<SubAppBackups>(stream=0x00007fff5fbfb630, backups=0x00000001103681f0, context=0x0000000110365740) at MultiApp.h:395
    frame #5: 0x0000000101734535 libmoose-dbg.0.dylib`void storeHelper<SubAppBackups>(stream=0x00007fff5fbfb630, data=0x00000001103681f0, context=0x0000000110365740) at DataIO.h:530
    frame #6: 0x00000001017343c0 libmoose-dbg.0.dylib`RestartableData<SubAppBackups>::store(this=0x00000001103681c0, stream=0x00007fff5fbfb630) at RestartableData.h:152
    frame #7: 0x0000000101a57ae5 libmoose-dbg.0.dylib`RestartableDataIO::serializeRestartableData(this=0x00000001103856b0, restartable_data=size=40, stream=0x00007fff5fbfc610) at RestartableDataIO.C:101
    frame #8: 0x0000000101a56592 libmoose-dbg.0.dylib`RestartableDataIO::writeRestartableData(this=0x00000001103856b0, base_file_name="master_checkpoint_cp/0000.rd", restartable_datas=0x000000011080bb98, (null)=size=23) at RestartableDataIO.C:56
    frame #9: 0x000000010177eb5b libmoose-dbg.0.dylib`Checkpoint::output(this=0x0000000110385400, (null)=0x00007fff5fbfd3c4) at Checkpoint.C:126
    frame #10: 0x00000001017f993f libmoose-dbg.0.dylib`Output::outputStep(this=0x0000000110385400, type=0x00007fff5fbfd3c4) at Output.C:160
    frame #11: 0x00000001018044e7 libmoose-dbg.0.dylib`OutputWarehouse::outputStep(this=0x000000011080b5a8, type=EXEC_INITIAL) at OutputWarehouse.C:152
    frame #12: 0x0000000100d8fcf5 libmoose-dbg.0.dylib`FEProblemBase::outputStep(this=0x0000000111003c00, type=EXEC_INITIAL) at FEProblemBase.C:3904
    frame #13: 0x00000001013f82bf libmoose-dbg.0.dylib`Transient::init(this=0x0000000110361c30) at Transient.C:268
    frame #14: 0x000000010115fbd3 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000011080b200) at MooseApp.C:569
    frame #15: 0x000000010116236a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000011080b200) at MooseApp.C:727
    frame #16: 0x0000000100002383 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe3a8) at main.C:35
    frame #17: 0x00007fffd4cb1235 libdyld.dylib`start + 1
    frame #18: 0x00007fffd4cb1235 libdyld.dylib`start + 1
```

Alright, so sure enough the sub app time step has not been incremented by the
time the output is written to the checkpoint file.

With execute:
```
  * frame #0: 0x0000000101a5b7c5 libmoose-dbg.0.dylib`RestartableDataIO::createBackup(this=0x00007fff5fbfbe70) at RestartableDataIO.C:300
    frame #1: 0x0000000101161ea8 libmoose-dbg.0.dylib`MooseApp::backup(this=0x000000010f003200) at MooseApp.C:661
    frame #2: 0x0000000101729d2e libmoose-dbg.0.dylib`MultiApp::backup(this=0x0000000110c657c0) at MultiApp.C:342
    frame #3: 0x0000000101734569 libmoose-dbg.0.dylib`void dataStore<SubAppBackups>(stream=0x00007fff5fbfc280, backups=0x0000000110c68270, context=0x0000000110c657c0) at MultiApp.h:395
    frame #4: 0x0000000101734535 libmoose-dbg.0.dylib`void storeHelper<SubAppBackups>(stream=0x00007fff5fbfc280, data=0x0000000110c68270, context=0x0000000110c657c0) at DataIO.h:530
    frame #5: 0x00000001017343c0 libmoose-dbg.0.dylib`RestartableData<SubAppBackups>::store(this=0x0000000110c68240, stream=0x00007fff5fbfc280) at RestartableData.h:152
    frame #6: 0x0000000101a57ae5 libmoose-dbg.0.dylib`RestartableDataIO::serializeRestartableData(this=0x0000000110c85730, restartable_data=size=40, stream=0x00007fff5fbfd260) at RestartableDataIO.C:101
    frame #7: 0x0000000101a56592 libmoose-dbg.0.dylib`RestartableDataIO::writeRestartableData(this=0x0000000110c85730, base_file_name="master_checkpoint_cp/0001.rd", restartable_datas=0x000000010f867198, (null)=size=23) at RestartableDataIO.C:56
    frame #8: 0x000000010177eb5b libmoose-dbg.0.dylib`Checkpoint::output(this=0x0000000110c85480, (null)=0x00007fff5fbfe014) at Checkpoint.C:126
    frame #9: 0x00000001017f993f libmoose-dbg.0.dylib`Output::outputStep(this=0x0000000110c85480, type=0x00007fff5fbfe014) at Output.C:160
    frame #10: 0x00000001018044e7 libmoose-dbg.0.dylib`OutputWarehouse::outputStep(this=0x000000010f866ba8, type=EXEC_TIMESTEP_END) at OutputWarehouse.C:152
    frame #11: 0x0000000100d8fcf5 libmoose-dbg.0.dylib`FEProblemBase::outputStep(this=0x0000000111003c00, type=EXEC_TIMESTEP_END) at FEProblemBase.C:3904
    frame #12: 0x00000001013f9d15 libmoose-dbg.0.dylib`Transient::endStep(this=0x0000000110c61cb0, input_time=-1) at Transient.C:601
    frame #13: 0x00000001013f84bb libmoose-dbg.0.dylib`Transient::execute(this=0x0000000110c61cb0) at Transient.C:327
    frame #14: 0x000000010115fc03 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010f866800) at MooseApp.C:571
    frame #15: 0x000000010116236a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010f866800) at MooseApp.C:727
    frame #16: 0x0000000100002383 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe3a8) at main.C:35
    frame #17: 0x00007fffd4cb1235 libdyld.dylib`start + 1
    frame #18: 0x00007fffd4cb1235 libdyld.dylib`start + 1
```

Here's bt for incrementing step:
```
  * frame #0: 0x0000000101745be0 libmoose-dbg.0.dylib`TransientMultiApp::advanceStep(this=0x0000000110c657c0) at TransientMultiApp.C:423
    frame #1: 0x0000000100d86069 libmoose-dbg.0.dylib`FEProblemBase::advanceMultiApps(this=0x0000000111003c00, type=EXEC_TIMESTEP_END) at FEProblemBase.C:3289
    frame #2: 0x00000001013f8643 libmoose-dbg.0.dylib`Transient::incrementStepOrReject(this=0x0000000110c61cb0) at Transient.C:372
    frame #3: 0x00000001013f8447 libmoose-dbg.0.dylib`Transient::execute(this=0x0000000110c61cb0) at Transient.C:317
    frame #4: 0x000000010115fc03 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010f866800) at MooseApp.C:571
    frame #5: 0x000000010116236a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010f866800) at MooseApp.C:727
    frame #6: 0x0000000100002383 moose_test-dbg`main(argc=3, argv=0x00007fff5fbfe3a8) at main.C:35
    frame #7: 0x00007fffd4cb1235 libdyld.dylib`start + 1
    frame #8: 0x00007fffd4cb1235 libdyld.dylib`start + 1
```

In Transient::execute
```
  while (true)
  {
    if (_first != true)
      incrementStepOrReject();

    _first = false;

    if (!keepGoing())
      break;

    preStep();
    computeDT();
    takeStep();
    endStep();
    postStep();

    _steps_taken++;
  }
```

In Transient::incrementStepOrReject
```
#ifdef LIBMESH_ENABLE_AMR
      _problem.adaptMesh();
#endif

      _time_old = _time; // = _time_old + _dt;
      _t_step++;

      _problem.advanceState();

      _problem.advanceMultiApps(EXEC_TIMESTEP_BEGIN);
      _problem.advanceMultiApps(EXEC_TIMESTEP_END);
    }
```

In Transient::takeStep
```
solveStep
```

In Transient::solveStep
```
  _dt_old = _dt;

  if (input_dt == -1.0)
    _dt = computeConstrainedDT();
  else
    _dt = input_dt;

  Real current_dt = _dt;

  _problem.onTimestepBegin();

  // Increment time
  _time = _time_old + _dt;
```

# 11/30/17


Alright, we make just as many calls to `serializeRestartableData` in the old
restart case as in the new. The difference is that in that final time step with
the new case, the time has not yet been incremented.

# 12/5/17

```
/Users/lindad/projects/moose/test/moose_test-opt -i function_dt_master.i --recover --recoversuffix cpr
```
```
Running command: /Users/lindad/projects/moose/test/moose_test-opt -i
picard_rel_tol_master.i --half-transient Outputs/checkpoint=true
```

```
* thread #1, queue = 'com.apple.main-thread', stop reason = step over
  * frame #0: 0x0000000101176ba1 libmoose-dbg.0.dylib`MooseApp::runInputFile(this=0x0000000112000000) at MooseApp.C:557
    frame #1: 0x000000010173e069 libmoose-dbg.0.dylib`MultiApp::createApp(this=0x0000000110364030, i=0, start_time=0) at MultiApp.C:580
    frame #2: 0x000000010173ac97 libmoose-dbg.0.dylib`MultiApp::initialSetup(this=0x0000000110364030) at MultiApp.C:205
    frame #3: 0x000000010175e70f libmoose-dbg.0.dylib`TransientMultiApp::initialSetup(this=0x0000000110364030) at TransientMultiApp.C:134
    frame #4: 0x0000000100ed295a libmoose-dbg.0.dylib`MooseObjectWarehouse<MultiApp>::initialSetup(this=0x0000000111006a38, tid=0) const at MooseObjectWarehouse.h:69
    frame #5: 0x0000000100d498cb libmoose-dbg.0.dylib`FEProblemBase::initialSetup(this=0x0000000111005400) at FEProblemBase.C:666
    frame #6: 0x00000001014153e5 libmoose-dbg.0.dylib`Transient::init(this=0x0000000110365560) at Transient.C:262
    frame #7: 0x00000001011774b3 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000011080b200) at MooseApp.C:599
    frame #8: 0x0000000101179c4a libmoose-dbg.0.dylib`MooseApp::run(this=0x000000011080b200) at MooseApp.C:757
    frame #9: 0x0000000100002383 moose_test-dbg`main(argc=6, argv=0x00007fff5fbfe370) at main.C:35
    frame #10: 0x00007fffd4cb1235 libdyld.dylib`start + 1
```

Possible sources of jacobian error:

- ThermalContact
- Creep/elasticity

Radial return.

Issues with only getting one actions' built objects per syntax association,
whereas I expect 3. This comes from:
```
  syntax.registerActionSyntax("PlenumPressureUOAction", "BCs/PlenumPressure/*");
  syntax.registerActionSyntax("CavityPressureAction",   "BCs/PlenumPressure/*");
  syntax.registerActionSyntax("CavityPressurePPAction", "BCs/PlenumPressure/*");
```

This seems strange because for Burnup for example there are three actions
associated with one syntax:
```
  syntax.registerActionSyntax("BurnupAuxKernelsAction", "Burnup/*", "add_aux_kernel");
  syntax.registerActionSyntax("BurnupAuxVarsAction",    "Burnup/*", "add_aux_variable");
  syntax.registerActionSyntax("BurnupFunctionAction",   "Burnup/*", "add_function");
```

However, note that there is a difference here: the former syntax registratino
has two arguments whereas the latter has three! Could be related to what Daniel
said about appending as opposed to setting.

# 12/11/17

After linearizing with the Newton algorithm, we solve:

Ax = b

The accuracy of the Jacobian doesn't determine the number of linear iterations
required to solve this system; only application of a preconditioner will change
the effeciency of the iterative process. However, for Jacobian-Free
Newton-Krylov, the closer the preconditoning matrix is to reproducing the
Jacobian-vector products genrerated by the JFNK process, the better the
efficiency of the solve will be.

The number of non-linear iterations is determined by the accuracy of the Jacobian.

# 12/18/17

Here's the backtrace of console residual information. Looks like we just created
our own custom Petsc monitor

```
#0  Console::output (this=0xac2070, type=@0x7fffffffcfd4: EXEC_NONLINEAR) at /home/lindad/projects/moose/framework/src/outputs/Console.C:315
#1  0x00007ffff2a72cbd in Output::outputStep (this=0xac2070, type=@0x7fffffffcfd4: EXEC_NONLINEAR)
    at /home/lindad/projects/moose/framework/src/outputs/Output.C:160
#2  0x00007ffff2a3fae0 in PetscOutput::petscNonlinearOutput (its=0, norm=0.0060435771370377948, void_ptr=0xac2070)
    at /home/lindad/projects/moose/framework/src/outputs/PetscOutput.C:195
#3  0x00007fffebf164e3 in SNESMonitor (snes=0xcd8b10, iter=0, rnorm=0.0060435771370377948) at /home/lindad/petsc/src/snes/interface/snes.c:3520
#4  0x00007fffebea8eaa in SNESSolve_NEWTONLS (snes=0xcd8b10) at /home/lindad/petsc/src/snes/impls/ls/ls.c:185
#5  0x00007fffebf1b5fe in SNESSolve (snes=0xcd8b10, b=0x0, x=0x96ab00) at /home/lindad/petsc/src/snes/interface/snes.c:4122
#6  0x00007fffeeec6c92 in libMesh::PetscNonlinearSolver<double>::solve (this=0x8da920, jac_in=..., x_in=..., r_in=...)
    at ../src/solvers/petsc_nonlinear_solver.C:702
#7  0x00007fffeef471f3 in libMesh::NonlinearImplicitSystem::solve (this=0x8d9b40) at ../src/systems/nonlinear_implicit_system.C:183
#8  0x00007ffff2a234ae in TimeIntegrator::solve (this=0x904ec0) at /home/lindad/projects/moose/framework/src/timeintegrators/TimeIntegrator.C:53
#9  0x00007ffff257db3c in NonlinearSystem::solve (this=0x8d7830) at /home/lindad/projects/moose/framework/src/base/NonlinearSystem.C:161
#10 0x00007ffff21e3455 in FEProblemBase::solve (this=0x8a1db0) at /home/lindad/projects/moose/framework/src/base/FEProblemBase.C:3752
#11 0x00007ffff282a019 in TimeStepper::step (this=0xab04c0) at /home/lindad/projects/moose/framework/src/timesteppers/TimeStepper.C:160
#12 0x00007ffff2b0314c in Transient::solveStep (this=0x912e60, input_dt=-1) at /home/lindad/projects/moose/framework/src/executioners/Transient.C:495
#13 0x00007ffff2b02a54 in Transient::takeStep (this=0x912e60, input_dt=-1) at /home/lindad/projects/moose/framework/src/executioners/Transient.C:408
#14 0x00007ffff2b0266c in Transient::execute (this=0x912e60) at /home/lindad/projects/moose/framework/src/executioners/Transient.C:326
#15 0x00007ffff26481c9 in MooseApp::executeExecutioner (this=0x72b900) at /home/lindad/projects/moose/framework/src/base/MooseApp.C:601
#16 0x00007ffff2648ddd in MooseApp::run (this=0x72b900) at /home/lindad/projects/moose/framework/src/base/MooseApp.C:757
#17 0x000000000041db68 in main (argc=16, argv=0x7fffffffd8b8) at /home/lindad/projects/bison/src/main.C:49
```

# 12/21/17

Failed tests:

outputs/vtk.files_parallel............................................................ FAILED (MISSING FILES)
functions/image_function.threshold_adapt_parallel_check_files......................... FAILED (MISSING FILES)
misc/exception.parallel_error_jacobian_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
vectorpostprocessors/csv_reader.tester_fail.................................. FAILED (EXPECTED ERROR MISSING)
misc/exception.parallel_error_residual_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
restart/kernel_restartable.parallel_error2................................... FAILED (EXPECTED ERROR MISSING)
materials/material_dependency.dont_reinit_mat_for_aux....................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac...................................... FAILED (EXODIFF)
mesh/custom_partitioner.custom_linear_partitioner........................................... FAILED (EXODIFF)
multiapps/picard_failure.test............................................................... FAILED (EXODIFF)
ics/depend_on_uo.ic_depend_on_uo............................................................ FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_2D.......................................... FAILED (EXODIFF)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.......................... FAILED (EXODIFF)
mesh/subdomain_partitioner.subdomain_partitioner............................................ FAILED (EXODIFF)
misc/exception.parallel_exception_residual_transient_non_zero_rank.......................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_2D.......................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac...................................... FAILED (EXODIFF)
utils/random.test_uo_par_mesh............................................................... FAILED (EXODIFF)
mesh/centroid_partitioner.centroid_partitioner_test......................................... FAILED (EXODIFF)
utils/random.test_par_mesh.................................................................. FAILED (EXODIFF)
restart/restartable_types.first_parallel..................................................... FAILED (ERRMSG)
mesh/nemesis.nemesis_test.................................................................... FAILED (ERRMSG)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)
userobjects/setup_interface_count.GeneralUserObject......................................... FAILED
(CSVDIFF)

25 fails

outputs/vtk.files_parallel............................................................ FAILED (MISSING FILES)
misc/exception.parallel_error_jacobian_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
vectorpostprocessors/csv_reader.tester_fail.................................. FAILED (EXPECTED ERROR MISSING)
misc/exception.parallel_error_residual_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
restart/kernel_restartable.parallel_error2................................... FAILED (EXPECTED ERROR MISSING)
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac...................................... FAILED (EXODIFF)
mesh/custom_partitioner.custom_linear_partitioner........................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac...................................... FAILED (EXODIFF)
mesh/subdomain_partitioner.subdomain_partitioner............................................ FAILED (EXODIFF)
misc/exception.parallel_exception_residual_transient_non_zero_rank.......................... FAILED (EXODIFF)
ics/depend_on_uo.ic_depend_on_uo............................................................ FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_2D.......................................... FAILED (EXODIFF)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.......................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_2D.......................................... FAILED (EXODIFF)
utils/random.test_uo_par_mesh............................................................... FAILED (EXODIFF)
utils/random.test_par_mesh.................................................................. FAILED (EXODIFF)
mesh/centroid_partitioner.centroid_partitioner_test......................................... FAILED (EXODIFF)
restart/restartable_types.first_parallel..................................................... FAILED (ERRMSG)
mesh/nemesis.nemesis_test.................................................................... FAILED (ERRMSG)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)
userobjects/setup_interface_count.GeneralUserObject......................................... FAILED (CSVDIFF)
functions/image_function.threshold_adapt_parallel_check_files................................. FAILED (CRASH)

22 fails

outputs/vtk.files_parallel............................................................ FAILED (MISSING FILES)
misc/exception.parallel_error_residual_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
vectorpostprocessors/csv_reader.tester_fail.................................. FAILED (EXPECTED ERROR MISSING)
misc/exception.parallel_error_jacobian_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
restart/kernel_restartable.parallel_error2................................... FAILED (EXPECTED ERROR MISSING)
materials/material_dependency.dont_reinit_mat_for_aux....................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac...................................... FAILED (EXODIFF)
mesh/custom_partitioner.custom_linear_partitioner........................................... FAILED (EXODIFF)
multiapps/picard_failure.test............................................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac...................................... FAILED (EXODIFF)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.......................... FAILED (EXODIFF)
misc/exception.parallel_exception_residual_transient_non_zero_rank.......................... FAILED (EXODIFF)
ics/depend_on_uo.ic_depend_on_uo............................................................ FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_2D.......................................... FAILED (EXODIFF)
mesh/subdomain_partitioner.subdomain_partitioner............................................ FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_2D.......................................... FAILED (EXODIFF)
mesh/centroid_partitioner.centroid_partitioner_test......................................... FAILED (EXODIFF)
utils/random.test_uo_par_mesh............................................................... FAILED (EXODIFF)
utils/random.test_par_mesh.................................................................. FAILED (EXODIFF)
mesh/nemesis.nemesis_test.................................................................... FAILED (ERRMSG)
restart/restartable_types.first_parallel..................................................... FAILED (ERRMSG)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)
userobjects/setup_interface_count.GeneralUserObject......................................... FAILED (CSVDIFF)
functions/image_function.threshold_adapt_parallel_check_files................................. FAILED (CRASH)
-------------------------------------------------------------------------------------------------------------
Ran 1493 tests in 50.3 seconds
1468 passed, 32 skipped, 0 pending, 25 FAILED


This also failed on my projects2 directory!

misc/exception.parallel_error_jacobian_transient_non_zero_rank...................... FAILED (NO EXPECTED ERR)
vectorpostprocessors/csv_reader.tester_fail......................................... FAILED (NO EXPECTED ERR)
misc/exception.parallel_error_residual_transient_non_zero_rank...................... FAILED (NO EXPECTED ERR)
restart/kernel_restartable.parallel_error2.......................................... FAILED (NO EXPECTED ERR)
outputs/vtk.files_parallel............................................................ FAILED (MISSING FILES)
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac...................................... FAILED (EXODIFF)
mesh/custom_partitioner.custom_linear_partitioner........................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac...................................... FAILED (EXODIFF)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.......................... FAILED (EXODIFF)
misc/exception.parallel_exception_residual_transient_non_zero_rank.......................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_2D.......................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_2D.......................................... FAILED (EXODIFF)
ics/depend_on_uo.ic_depend_on_uo............................................................ FAILED (EXODIFF)
mesh/subdomain_partitioner.subdomain_partitioner............................................ FAILED (EXODIFF)
mesh/centroid_partitioner.centroid_partitioner_test......................................... FAILED (EXODIFF)
utils/random.test_par_mesh.................................................................. FAILED (EXODIFF)
utils/random.test_uo_par_mesh............................................................... FAILED (EXODIFF)
restart/restartable_types.first_parallel..................................................... FAILED (ERRMSG)
mesh/nemesis.nemesis_test.................................................................... FAILED (ERRMSG)
userobjects/setup_interface_count.GeneralUserObject......................................... FAILED (CSVDIFF)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)
functions/image_function.threshold_adapt_parallel_check_files................................. FAILED (CRASH)

Here are errors on rod which I trust quite a bit more:

dgkernels/jacobian_testing.jacobian_test............................................ FAILED (NO EXPECTED OUT)
interfacekernels/1d_interface.jacobian_test......................................... FAILED (NO EXPECTED OUT)
interfacekernels/1d_interface.single_variable_jacobian_test......................... FAILED (NO EXPECTED OUT)
interfacekernels/1d_interface.mixed_shapes_jacobian_test............................ FAILED (NO EXPECTED OUT)
misc/jacobian.offdiag............................................................... FAILED (NO EXPECTED OUT)
misc/jacobian.simple................................................................ FAILED (NO EXPECTED OUT)
misc/jacobian.med................................................................... FAILED (NO EXPECTED OUT)
interfacekernels/2d_interface.jacobian_test......................................... FAILED (NO EXPECTED OUT)
multiapps/picard_failure.test............................................................... FAILED (EXODIFF)
materials/material_dependency.dont_reinit_mat_for_aux....................................... FAILED (EXODIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)

Ok, we call `updateMesh` in both our residual and jacobian computing chains. So
that means we move our nodes around and perhaps change our penetration info.

What should be done is that for evaulating non-linear residuals we should first
make sure the mesh is updated, e.g. 1) move our nodes 2) update our geometric
and penetration info. Then we should 3) update the set of captured contact
points and then 4) finally determine our contact residuals.

Ok, we really might be able to do something with this snes update. So F, the
residual that is sent to KSPSolve is a pointer to snes->vec_func. So as long as
we have access to snes we should be able to modify F.

# 12/22/17

Here's another computeResidual call that I don't think is necessary: (this is
the first time this computeResidual call is made)

```
  * frame #0: 0x00000001033ce7dd libmoose-dbg.0.dylib`FEProblemBase::computeResidualType(this=0x0000000114005e18, soln=0x0000000112f56730, residual=0x0000000112f56bc0, type=KT_ALL) at FEProblemBase.C:4010
    frame #1: 0x00000001033ce170 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000114005e18, soln=0x0000000112f56730, residual=0x0000000112f56bc0) at FEProblemBase.C:3978
    frame #2: 0x00000001033ce126 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000114005e18, (null)=0x0000000112f56280, soln=0x0000000112f56730, residual=0x0000000112f56bc0) at FEProblemBase.C:3970
    frame #3: 0x00000001038bbc2e libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x0000000114007e18) at NonlinearSystem.C:140
    frame #4: 0x00000001033c9f32 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x0000000114005e18) at FEProblemBase.C:3756
    frame #5: 0x0000000104192914 libmoose-dbg.0.dylib`TimeStepper::step(this=0x00000001150169a8) at TimeStepper.C:160
    frame #6: 0x0000000103ace1b3 libmoose-dbg.0.dylib`Transient::solveStep(this=0x0000000112f74518, input_dt=-1) at Transient.C:495
    frame #7: 0x0000000103acd92a libmoose-dbg.0.dylib`Transient::takeStep(this=0x0000000112f74518, input_dt=-1) at Transient.C:408
    frame #8: 0x0000000103acd263 libmoose-dbg.0.dylib`Transient::execute(this=0x0000000112f74518) at Transient.C:326
    frame #9: 0x0000000103836423 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:604
    frame #10: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #11: 0x0000000100002ab5 bison-dbg`main(argc=4, argv=0x00007fff5fbfec08) at main.C:50
    frame #12: 0x00007fffa51b6235 libdyld.dylib`start + 1
    frame #13: 0x00007fffa51b6235 libdyld.dylib`start + 1
```

Second call, this calls is from the petsc solve, determining the initial
non-linear residual. It is called outside the non-linera iteration loop:
```
  * frame #0: 0x00000001033ce7dd libmoose-dbg.0.dylib`FEProblemBase::computeResidualType(this=0x0000000114005e18, soln=0x0000000112f56730, residual=0x00007fff5fbfd468, type=KT_ALL) at FEProblemBase.C:4010
    frame #1: 0x00000001033ce170 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000114005e18, soln=0x0000000112f56730, residual=0x00007fff5fbfd468) at FEProblemBase.C:3978
    frame #2: 0x00000001033ce126 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000114005e18, (null)=0x0000000112f56280, soln=0x0000000112f56730, residual=0x00007fff5fbfd468) at FEProblemBase.C:3970
    frame #3: 0x00000001038ba2eb libmoose-dbg.0.dylib`Moose::compute_residual(soln=0x0000000112f56730, residual=0x00007fff5fbfd468, sys=0x0000000112f56280) at NonlinearSystem.C:45
    frame #4: 0x0000000106fcce1d libmesh_dbg.0.dylib`::__libmesh_petsc_snes_residual(snes=0x0000000115890260, x=0x000000011403aa60, r=0x0000000115808a60, ctx=0x0000000112f56e70) at petsc_nonlinear_solver.C:130
    frame #5: 0x000000010fb4695e libpetsc.3.08.dylib`SNESComputeFunction(snes=0x0000000115890260, x=0x000000011403aa60, y=0x0000000115808a60) at snes.c:2195
    frame #6: 0x000000010fbdbdcd libpetsc.3.08.dylib`SNESSolve_NEWTONLS(snes=0x0000000115890260) at ls.c:175
    frame #7: 0x000000010fb57f88 libpetsc.3.08.dylib`SNESSolve(snes=0x0000000115890260, b=0x0000000000000000, x=0x000000011403aa60) at snes.c:4106
    frame #8: 0x0000000106fca2dc libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x0000000112f56e70, jac_in=0x0000000112f56da0, x_in=0x0000000112f56610, r_in=0x0000000112f56bc0, (null)=0.0000000001, (null)=100) at petsc_nonlinear_solver.C:702
    frame #9: 0x00000001070eac80 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x0000000112f56280) at nonlinear_implicit_system.C:181
    frame #10: 0x00000001041739ce libmoose-dbg.0.dylib`TimeIntegrator::solve(this=0x0000000112f7b0b8) at TimeIntegrator.C:53
    frame #11: 0x00000001038bbea2 libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x0000000114007e18) at NonlinearSystem.C:161
    frame #12: 0x00000001033c9f32 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x0000000114005e18) at FEProblemBase.C:3756
    frame #13: 0x0000000104192914 libmoose-dbg.0.dylib`TimeStepper::step(this=0x00000001150169a8) at TimeStepper.C:160
    frame #14: 0x0000000103ace1b3 libmoose-dbg.0.dylib`Transient::solveStep(this=0x0000000112f74518, input_dt=-1) at Transient.C:495
    frame #15: 0x0000000103acd92a libmoose-dbg.0.dylib`Transient::takeStep(this=0x0000000112f74518, input_dt=-1) at Transient.C:408
    frame #16: 0x0000000103acd263 libmoose-dbg.0.dylib`Transient::execute(this=0x0000000112f74518) at Transient.C:326
    frame #17: 0x0000000103836423 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:604
    frame #18: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #19: 0x0000000100002ab5 bison-dbg`main(argc=4, argv=0x00007fff5fbfec08) at main.C:50
    frame #20: 0x00007fffa51b6235 libdyld.dylib`start + 1
    frame #21: 0x00007fffa51b6235 libdyld.dylib`start + 1
```

- Just following residual and jacobian computing threads. The next call is to our
compute jacobian chain and it's from the petsc non-linear iteration loop
(calculating our preconditioner).

- Next call is to the compute residual thread and this is for approximating the
  jacobian action through finite differencing. This will now get called once per
  additional linear iteration.


# 1/3/17

Ok, it looks like we get the following messages for ksp->reason:

- Nonliner residual > 0: KSP_CONVERGED_RTOL
- Linear residual: KSP_CONVERGED_ITERATING
- Nolinear residual = 0: KSP_CONVERGED_ITERATING

Here's state of system at zeroth non-linear residual evaluation:

- snes->ksp->reason = KSP_CONVERGED_ITERATING
- snes->iter = 0
- snes->ksp->its = 0

And state when doing first linear residual evaluation:

- snes->ksp->reason = KSP_CONVERGED_ITERATING
- snes->iter = 0
- snes->ksp->its = 0

The same :-(

But the fix: snes->nfuncs

So to detect non-linear residual evaluations we could use the check:

if (snes->nfuncs == 0 || snes->ksp->reason == KSP_CONVERGED_RTOL)

May also need some additional comparisons for snes->ksp->reason.

KSP_CONVERGED_RTOL == 2

Final check:

if (snes->nfuncs == 0 || snes->ksp->reason > 0)

This is because enum values greater than zero indicate converged, less than zero
indicates diverged, equal to zero indicates still iterating.

I have seen cases where even though the ksp diverges (I've seen
KSP_DIVERGED_BREAKDOWN), the non-linear solve continues, e.g. another non-linear
residual is evaluated. So better would be:

if (snes->nfuncs == 0 || snes->ksp->reason != KSP_CONVERGED_ITERATING)

DisplacedProblem::updateMesh is called from FEProblemBase::computeResidualType

There are about a million calls to GeometricSearchData::update

  * frame #0: 0x0000000103b6025e libmoose-dbg.0.dylib`PenetrationLocator::detectPenetration(this=0x0000000112b1e7d0) at PenetrationLocator.C:91
    frame #1: 0x0000000103b36fbc libmoose-dbg.0.dylib`GeometricSearchData::update(this=0x00000001138106b8, type=ALL) at GeometricSearchData.C:83
    frame #2: 0x000000010334be82 libmoose-dbg.0.dylib`DisplacedProblem::updateGeomSearch(this=0x0000000113810218, type=ALL) at DisplacedProblem.C:743
    frame #3: 0x00000001033d40c6 libmoose-dbg.0.dylib`FEProblemBase::updateGeomSearch(this=0x0000000113806618, type=ALL) at FEProblemBase.C:4396
    frame #4: 0x00000001038f1901 libmoose-dbg.0.dylib`NonlinearSystemBase::augmentSparsity(this=0x000000011380c618, sparsity=0x000000011480d518, n_nz=size=20, n_oz=size=20) at NonlinearSystemBase.C:2367
    frame #5: 0x00000001039d9705 libmoose-dbg.0.dylib`extraSparsity(sparsity=0x000000011480d518, n_nz=size=20, n_oz=size=20, context=0x000000011380c618) at SystemBase.C:47
    frame #6: 0x00000001057b469c libmesh_dbg.0.dylib`libMesh::DofMap::build_sparsity(this=0x0000000112f46ff0, mesh=0x0000000112dcffb0) const at dof_map.C:113
    frame #7: 0x000000010580226e libmesh_dbg.0.dylib`libMesh::DofMap::compute_sparsity(this=0x0000000112f46ff0, mesh=0x0000000112dcffb0) at dof_map.C:1737
    frame #8: 0x00000001070cbc99 libmesh_dbg.0.dylib`libMesh::ImplicitSystem::reinit(this=0x0000000112f46a80) at implicit_system.C:176
    frame #9: 0x00000001070e9867 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::reinit(this=0x0000000112f46a80) at nonlinear_implicit_system.C:90
    frame #10: 0x0000000107206a45 libmesh_dbg.0.dylib`libMesh::TransientSystem<libMesh::NonlinearImplicitSystem>::reinit(this=0x0000000112f46a80) at transient_system.C:90
    frame #11: 0x00000001070440c4 libmesh_dbg.0.dylib`libMesh::EquationSystems::reinit(this=0x00000001138068e8) at equation_systems.C:258
    frame #12: 0x0000000103380c0f libmoose-dbg.0.dylib`FEProblemBase::reinitBecauseOfGhostingOrNewGeomObjects(this=0x0000000113806618) at FEProblemBase.C:3015
    frame #13: 0x00000001033761fd libmoose-dbg.0.dylib`FEProblemBase::initialSetup(this=0x0000000113806618) at FEProblemBase.C:626
    frame #14: 0x0000000103accf95 libmoose-dbg.0.dylib`Transient::init(this=0x0000000112f636c8) at Transient.C:262
    frame #15: 0x00000001038363f3 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:602
    frame #16: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #17: 0x0000000100002ab5 bison-dbg`main(argc=3, argv=0x00007fff5fbfec90) at main.C:50
    frame #18: 0x00007fffa51b6235 libdyld.dylib`start + 1

  * frame #0: 0x0000000103b6025e libmoose-dbg.0.dylib`PenetrationLocator::detectPenetration(this=0x0000000112b1e7d0) at PenetrationLocator.C:91
    frame #1: 0x0000000103b36fbc libmoose-dbg.0.dylib`GeometricSearchData::update(this=0x00000001138106b8, type=ALL) at GeometricSearchData.C:83
    frame #2: 0x00000001033489f8 libmoose-dbg.0.dylib`DisplacedProblem::updateMesh(this=0x0000000113810218) at DisplacedProblem.C:201
    frame #3: 0x0000000103376267 libmoose-dbg.0.dylib`FEProblemBase::initialSetup(this=0x0000000113806618) at FEProblemBase.C:630
    frame #4: 0x0000000103accf95 libmoose-dbg.0.dylib`Transient::init(this=0x0000000112f636c8) at Transient.C:262
    frame #5: 0x00000001038363f3 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:602
    frame #6: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #7: 0x0000000100002ab5 bison-dbg`main(argc=3, argv=0x00007fff5fbfec90) at main.C:50
    frame #8: 0x00007fffa51b6235 libdyld.dylib`start + 1

  * frame #0: 0x0000000103b6025e libmoose-dbg.0.dylib`PenetrationLocator::detectPenetration(this=0x0000000112b1e7d0) at PenetrationLocator.C:91
    frame #1: 0x0000000103b36fbc libmoose-dbg.0.dylib`GeometricSearchData::update(this=0x00000001138106b8, type=ALL) at GeometricSearchData.C:83
    frame #2: 0x000000010334be82 libmoose-dbg.0.dylib`DisplacedProblem::updateGeomSearch(this=0x0000000113810218, type=ALL) at DisplacedProblem.C:743
    frame #3: 0x00000001033d40c6 libmoose-dbg.0.dylib`FEProblemBase::updateGeomSearch(this=0x0000000113806618, type=ALL) at FEProblemBase.C:4396
    frame #4: 0x00000001033762a8 libmoose-dbg.0.dylib`FEProblemBase::initialSetup(this=0x0000000113806618) at FEProblemBase.C:633
    frame #5: 0x0000000103accf95 libmoose-dbg.0.dylib`Transient::init(this=0x0000000112f636c8) at Transient.C:262
    frame #6: 0x00000001038363f3 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:602
    frame #7: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #8: 0x0000000100002ab5 bison-dbg`main(argc=3, argv=0x00007fff5fbfec90) at main.C:50
    frame #9: 0x00007fffa51b6235 libdyld.dylib`start + 1

  * frame #0: 0x0000000103b6025e libmoose-dbg.0.dylib`PenetrationLocator::detectPenetration(this=0x0000000112b1e7d0) at PenetrationLocator.C:91
    frame #1: 0x0000000103b36fbc libmoose-dbg.0.dylib`GeometricSearchData::update(this=0x00000001138106b8, type=ALL) at GeometricSearchData.C:83
    frame #2: 0x00000001033489f8 libmoose-dbg.0.dylib`DisplacedProblem::updateMesh(this=0x0000000113810218) at DisplacedProblem.C:201
    frame #3: 0x00000001033cee65 libmoose-dbg.0.dylib`FEProblemBase::computeResidualType(this=0x0000000113806618, soln=0x0000000112f46f30, residual=0x0000000112f473c0, type=KT_ALL) at FEProblemBase.C:4033
    frame #4: 0x00000001033ce170 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000113806618, soln=0x0000000112f46f30, residual=0x0000000112f473c0) at FEProblemBase.C:3978
    frame #5: 0x00000001033ce126 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000113806618, (null)=0x0000000112f46a80, soln=0x0000000112f46f30, residual=0x0000000112f473c0) at FEProblemBase.C:3970
    frame #6: 0x00000001038bbc2e libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000011380c618) at NonlinearSystem.C:140
    frame #7: 0x00000001033c9f32 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x0000000113806618) at FEProblemBase.C:3756
    frame #8: 0x0000000104192914 libmoose-dbg.0.dylib`TimeStepper::step(this=0x0000000112b03028) at TimeStepper.C:160
    frame #9: 0x0000000103ace1b3 libmoose-dbg.0.dylib`Transient::solveStep(this=0x0000000112f636c8, input_dt=-1) at Transient.C:495
    frame #10: 0x0000000103acd92a libmoose-dbg.0.dylib`Transient::takeStep(this=0x0000000112f636c8, input_dt=-1) at Transient.C:408
    frame #11: 0x0000000103acd263 libmoose-dbg.0.dylib`Transient::execute(this=0x0000000112f636c8) at Transient.C:326
    frame #12: 0x0000000103836423 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:604
    frame #13: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #14: 0x0000000100002ab5 bison-dbg`main(argc=3, argv=0x00007fff5fbfec90) at main.C:50
    frame #15: 0x00007fffa51b6235 libdyld.dylib`start + 1

  * frame #0: 0x0000000103b6025e libmoose-dbg.0.dylib`PenetrationLocator::detectPenetration(this=0x0000000112b1e7d0) at PenetrationLocator.C:91
    frame #1: 0x0000000103b36fbc libmoose-dbg.0.dylib`GeometricSearchData::update(this=0x00000001138106b8, type=ALL) at GeometricSearchData.C:83
    frame #2: 0x00000001033489f8 libmoose-dbg.0.dylib`DisplacedProblem::updateMesh(this=0x0000000113810218) at DisplacedProblem.C:201
    frame #3: 0x00000001033cee65 libmoose-dbg.0.dylib`FEProblemBase::computeResidualType(this=0x0000000113806618, soln=0x0000000112f46f30, residual=0x00007fff5fbfd4f8, type=KT_ALL) at FEProblemBase.C:4033
    frame #4: 0x00000001033ce170 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000113806618, soln=0x0000000112f46f30, residual=0x00007fff5fbfd4f8) at FEProblemBase.C:3978
    frame #5: 0x00000001033ce126 libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x0000000113806618, (null)=0x0000000112f46a80, soln=0x0000000112f46f30, residual=0x00007fff5fbfd4f8) at FEProblemBase.C:3970
    frame #6: 0x00000001038ba2eb libmoose-dbg.0.dylib`Moose::compute_residual(soln=0x0000000112f46f30, residual=0x00007fff5fbfd4f8, sys=0x0000000112f46a80) at NonlinearSystem.C:45
    frame #7: 0x0000000106fcce1d libmesh_dbg.0.dylib`::__libmesh_petsc_snes_residual(snes=0x00000001150c2460, x=0x0000000115059660, r=0x0000000115001460, ctx=0x0000000112f47670) at petsc_nonlinear_solver.C:130
    frame #8: 0x000000010fb4695e libpetsc.3.08.dylib`SNESComputeFunction(snes=0x00000001150c2460, x=0x0000000115059660, y=0x0000000115001460) at snes.c:2195
    frame #9: 0x000000010fbdbdcd libpetsc.3.08.dylib`SNESSolve_NEWTONLS(snes=0x00000001150c2460) at ls.c:175
    frame #10: 0x000000010fb57f88 libpetsc.3.08.dylib`SNESSolve(snes=0x00000001150c2460, b=0x0000000000000000, x=0x0000000115059660) at snes.c:4106
    frame #11: 0x0000000106fca2dc libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x0000000112f47670, jac_in=0x0000000112f475a0, x_in=0x0000000112f46e10, r_in=0x0000000112f473c0, (null)=0.0000000001, (null)=100) at petsc_nonlinear_solver.C:702
    frame #12: 0x00000001070eac80 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x0000000112f46a80) at nonlinear_implicit_system.C:181
    frame #13: 0x00000001041739ce libmoose-dbg.0.dylib`TimeIntegrator::solve(this=0x0000000112f6a268) at TimeIntegrator.C:53
    frame #14: 0x00000001038bbea2 libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000011380c618) at NonlinearSystem.C:161
    frame #15: 0x00000001033c9f32 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x0000000113806618) at FEProblemBase.C:3756
    frame #16: 0x0000000104192914 libmoose-dbg.0.dylib`TimeStepper::step(this=0x0000000112b03028) at TimeStepper.C:160
    frame #17: 0x0000000103ace1b3 libmoose-dbg.0.dylib`Transient::solveStep(this=0x0000000112f636c8, input_dt=-1) at Transient.C:495
    frame #18: 0x0000000103acd92a libmoose-dbg.0.dylib`Transient::takeStep(this=0x0000000112f636c8, input_dt=-1) at Transient.C:408
    frame #19: 0x0000000103acd263 libmoose-dbg.0.dylib`Transient::execute(this=0x0000000112f636c8) at Transient.C:326
    frame #20: 0x0000000103836423 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x0000000113015600) at MooseApp.C:604
    frame #21: 0x0000000103838b7f libmoose-dbg.0.dylib`MooseApp::run(this=0x0000000113015600) at MooseApp.C:761
    frame #22: 0x0000000100002ab5 bison-dbg`main(argc=3, argv=0x00007fff5fbfec90) at main.C:50
    frame #23: 0x00007fffa51b6235 libdyld.dylib`start + 1

DisplacedProblem::updateMesh is called from FEProblemBase::computeResidualType

Ok, we call `updateMesh` in both our residual and jacobian computing chains. So
that means we move our nodes around and perhaps change our penetration info.

What should be done is that for evaulating non-linear residuals we should first
make sure the mesh is updated, e.g. 1) move our nodes 2) update our geometric
and penetration info. Then we should 3) update the set of captured contact
points and then 4) finally determine our contact residuals.

0 Nonlinear |R| = 1.495345e-01

distance = 3.568794e-4

Resid1: -2220.8116223172437 master
Resid2: -1347.9038896496565 master
Resid3: 3568.7155119669001 slave

matrix-free solve:

 0 Nonlinear |R| = 4.786440e-02
    0 KSP unpreconditioned resid norm 4.786440058172e-02 true resid norm 4.786440058172e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.751252206702e-02 true resid norm 1.751252206702e-02 ||r(i)||/||b|| 3.658778101090e-01
    2 KSP unpreconditioned resid norm 7.106912253736e-03 true resid norm 7.106912253468e-03 ||r(i)||/||b|| 1.484801265052e-01
    3 KSP unpreconditioned resid norm 4.838598438201e-03 true resid norm 4.838598433259e-03 ||r(i)||/||b|| 1.010897112353e-01
    4 KSP unpreconditioned resid norm 9.402733137094e-18 true resid norm 4.587004841479e-03 ||r(i)||/||b|| 9.583332885674e-02
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 9.010475e-02
    0 KSP unpreconditioned resid norm 9.010475030347e-02 true resid norm 9.010475030347e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.522682782986e-02 true resid norm 2.522682782986e-02 ||r(i)||/||b|| 2.799722294873e-01
    2 KSP unpreconditioned resid norm 9.313539181987e-03 true resid norm 9.313539173420e-03 ||r(i)||/||b|| 1.033634646570e-01
    3 KSP unpreconditioned resid norm 6.214150880452e-03 true resid norm 6.214150821783e-03 ||r(i)||/||b|| 6.896585142130e-02
    4 KSP unpreconditioned resid norm 4.339074190343e-03 true resid norm 4.339074228001e-03 ||r(i)||/||b|| 4.815588760179e-02
    5 KSP unpreconditioned resid norm 2.431284496950e-03 true resid norm 2.431284505302e-03 ||r(i)||/||b|| 2.698286713091e-02
    6 KSP unpreconditioned resid norm 8.427780129200e-05 true resid norm 8.427787032991e-05 ||r(i)||/||b|| 9.353321555863e-04
    7 KSP unpreconditioned resid norm 1.037937761218e-11 true resid norm 1.758290685527e+00 ||r(i)||/||b|| 1.951385115219e+01
  Linear solve converged due to CONVERGED_RTOL iterations 7
 2 Nonlinear |R| = 1.649952e+00
    0 KSP unpreconditioned resid norm 1.649951538460e+00 true resid norm 1.649951538460e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.291152194435e-01 true resid norm 4.291152194435e-01 ||r(i)||/||b|| 2.600774686049e-01
    2 KSP unpreconditioned resid norm 1.477253634835e-01 true resid norm 1.477253611568e-01 ||r(i)||/||b|| 8.953315155828e-02
    3 KSP unpreconditioned resid norm 8.747908695732e-02 true resid norm 8.747908778489e-02 ||r(i)||/||b|| 5.301918616745e-02
    4 KSP unpreconditioned resid norm 3.817918283169e-02 true resid norm 3.817918188206e-02 ||r(i)||/||b|| 2.313957773432e-02
    5 KSP unpreconditioned resid norm 1.713179678242e-02 true resid norm 1.713178675329e-02 ||r(i)||/||b|| 1.038320602390e-02
    6 KSP unpreconditioned resid norm 8.948665136664e-03 true resid norm 8.948662682531e-03 ||r(i)||/||b|| 5.423591223099e-03
    7 KSP unpreconditioned resid norm 4.261823750819e-03 true resid norm 4.261808756404e-03 ||r(i)||/||b|| 2.582990261872e-03
    8 KSP unpreconditioned resid norm 1.791881412849e-03 true resid norm 1.791877253530e-03 ||r(i)||/||b|| 1.086018111297e-03
    9 KSP unpreconditioned resid norm 6.062748648289e-04 true resid norm 6.062751089429e-04 ||r(i)||/||b|| 3.674502522108e-04
   10 KSP unpreconditioned resid norm 3.765308187906e-16 true resid norm 1.901382189021e+08 ||r(i)||/||b|| 1.152386688154e+08
  Linear solve converged due to CONVERGED_RTOL iterations 10


Newton solve:
0 Nonlinear |R| = 4.562724e-02
    0 KSP unpreconditioned resid norm 4.562723505280e-02 true resid norm 4.562723505280e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.291610408077e-02 true resid norm 4.291610408077e-02 ||r(i)||/||b|| 9.405808620904e-01
    2 KSP unpreconditioned resid norm 1.351032209192e-02 true resid norm 1.351032209192e-02 ||r(i)||/||b|| 2.961021432985e-01
    3 KSP unpreconditioned resid norm 9.156650185358e-03 true resid norm 9.156650185358e-03 ||r(i)||/||b|| 2.006838716999e-01
    4 KSP unpreconditioned resid norm 5.990233609685e-03 true resid norm 5.990233609685e-03 ||r(i)||/||b|| 1.312863600600e-01
    5 KSP unpreconditioned resid norm 2.202010328664e-03 true resid norm 2.202010328664e-03 ||r(i)||/||b|| 4.826087590264e-02
    6 KSP unpreconditioned resid norm 7.283034643562e-04 true resid norm 7.283034643562e-04 ||r(i)||/||b|| 1.596203371765e-02
    7 KSP unpreconditioned resid norm 1.906180772829e-04 true resid norm 1.906180772829e-04 ||r(i)||/||b|| 4.177725804825e-03
    8 KSP unpreconditioned resid norm 1.066057187355e-05 true resid norm 1.066057187355e-05 ||r(i)||/||b|| 2.336449241603e-04
    9 KSP unpreconditioned resid norm 4.508408661214e-06 true resid norm 4.508408661215e-06 ||r(i)||/||b|| 9.880959597920e-05
   10 KSP unpreconditioned resid norm 1.662560900169e-16 true resid norm 1.680738557347e-16 ||r(i)||/||b|| 3.683630084975e-15
  Linear solve converged due to CONVERGED_RTOL iterations 10
 1 Nonlinear |R| = 1.204547e-01
    0 KSP unpreconditioned resid norm 1.204546822107e-01 true resid norm 1.204546822107e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.494828504040e-02 true resid norm 1.494828504040e-02 ||r(i)||/||b|| 1.240988292531e-01
    2 KSP unpreconditioned resid norm 1.233337228461e-03 true resid norm 1.233337228461e-03 ||r(i)||/||b|| 1.023901442290e-02
    3 KSP unpreconditioned resid norm 7.152188510119e-04 true resid norm 7.152188510118e-04 ||r(i)||/||b|| 5.937659191703e-03
    4 KSP unpreconditioned resid norm 2.620317764914e-04 true resid norm 2.620317764914e-04 ||r(i)||/||b|| 2.175355674702e-03
    5 KSP unpreconditioned resid norm 1.159634940781e-04 true resid norm 1.159634940781e-04 ||r(i)||/||b|| 9.627147068912e-04
    6 KSP unpreconditioned resid norm 4.413616365984e-05 true resid norm 4.413616365983e-05 ||r(i)||/||b|| 3.664130181560e-04
    7 KSP unpreconditioned resid norm 1.648901180449e-05 true resid norm 1.648901180450e-05 ||r(i)||/||b|| 1.368897539048e-04
    8 KSP unpreconditioned resid norm 2.804134281485e-06 true resid norm 2.804134281491e-06 ||r(i)||/||b|| 2.327957892568e-05
    9 KSP unpreconditioned resid norm 1.398263856639e-06 true resid norm 1.398263856637e-06 ||r(i)||/||b|| 1.160821506458e-05
   10 KSP unpreconditioned resid norm 1.738488310090e-17 true resid norm 1.658039547668e-17 ||r(i)||/||b|| 1.376484099446e-16
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 1.199152e-04

Matrix SMP:
initial solution vector on time step 3:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0172816
-0.3
-0.0172816
-0.3
-0.00869411
-0.199999
-0.00869411
-0.199999

After the first non-linear iteration:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0172816
-0.3
-0.0172816
-0.3
-0.0171587
-0.298569
-0.0171587
-0.298569

After the second non-linear iteration:
0.00143036
-0.000304093
0.
0.
0.
0.
-0.00394194
4.44861e-18
0.
0.
0.00143036
0.000304093
-0.0172816
-0.3
-0.0172816
-0.3
-0.0122453
-0.296938
-0.0176921
-0.298973

Matrix FDP:
1 non-linear:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0172816
-0.3
-0.0172816
-0.3
-0.0172816
-0.3
-0.0172816
-0.3

2 non-linear:
0.00132888
0.000125133
0.
0.
0.
0.
-0.00408206
0.000592353
0.
0.
0.00163352
0.000754935
-0.0172816
-0.3
-0.0172816
-0.3
-0.0130558
-0.298287
-0.0167746
-0.300244


With matrix free:
0 non-linear:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0172816
-0.3
-0.0172816
-0.3
-0.00869424
-0.2
-0.00869424
-0.2

1 non-linear:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0172816
-0.3
-0.0172816
-0.3
-0.0172816
-0.3
-0.0172816
-0.3

2 non-linaer:
0.0292859
-0.00622614
0.
0.
0.
0.
-0.0807092
-1.00132e-11
0.
0.
0.0292859
0.00622614
-0.0172816
-0.3
-0.0172816
-0.3
0.0861564
-0.237035
-0.0261878
-0.278707

nodes:

1: top right left block
2: top left left block
3: middle left left block
4: middle right left block
5: bottom left left block
6: bottom right left block
7: bottom right right block
8: top right right block
9: top left right block
10: bottom left right block

matrix-free:
0 residual:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0122123
0.0199942
0.0276254
0.0312878

1 residual:
0.
0.
0.
0.
0.
0.
0.0728156
1.93855e-12
0.
0.
0.
0.
0.
0.
0.
0.
-0.0455098
-6.57476e-12
-0.0273059
-1.78586e-11

2 residual:
-0.00107939
0.00195319
0.
0.
0.
0.
-1.334
0.168409
0.
0.
-0.00107939
-0.00195319
0.
0.
0.
0.
0.73583
-0.0914974
0.598851
-0.0772759

matrix SMP:
0 residual:
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
0.
-0.0122124
0.0199944
0.0276258
0.0312882

1 residual:
0.
0.
0.
0.
0.
0.
0.0715866
-4.96732e-18
0.
0.
0.
0.
0.
0.
0.
0.
-0.0447786
0.000330886
-0.0265873
0.000403093

2 residual:
-2.65348e-06
4.80612e-06
0.
0.
0.
0.
5.95969e-05
-2.45497e-05
0.
0.
-2.65348e-06
-4.80612e-06
0.
0.
0.
0.
-3.0547e-05
1.6313e-05
-1.91304e-05
7.87762e-06

matrix FDP:
0 nonlinear:
1 nonlinear:
0.
0.
0.
0.
0.
0.
0.0728156
-5.85789e-11
0.
0.
0.
0.
0.
0.
0.
0.
-0.0455098
1.20099e-10
-0.0273059
7.26427e-11

2 nonlinaer:
-2.50763e-06
4.85422e-06
0.
0.
0.
0.
1.89873e-05
0.000322573
0.
0.
-3.03928e-06
-5.30093e-06
0.
0.
0.
0.
-4.23526e-05
-0.000199893
2.79678e-05
-0.000122915


matrix-free from solid mechanics disp_y slave element:
  _val = size=4 {
    [0] = 0.000010295307058583291
    [1] = 0.000012199548279089882
    [2] = -0.000005363167860984916
    [3] = -0.000017131687476688257
  }


matrix from solid mechanics disp_y slave element:
  _val = size=4 {
    [0] = -330.88581794684291
    [1] = -403.0925249092017
    [2] = 330.88581794684291
    [3] = 403.0925249092017
  }

The stress tensors are fundamentally different from the matrix-free and matrix
cases:


matrix:

zeroth non-linaer residual:
(const SymmTensor) $97 = (_xx = -19266.730169981209, _yy = -8257.1700728490905,
_zz = -8257.1700728490905, _xy = -64103.233949937516, _yz = 0, _zx = 0)


first non-linear residual:
(lldb) p _stress[_qp]
(const SymmTensor) $96 = (_xx = -275.75369078613789, _yy = -118.1801531940591,
_zz = -118.1801531940591, _xy = -917.47292857006209, _yz = 0, _zx = 0)


Matrix-free:

0th non-linear residual:
1 linear residual:
(const SymmTensor) $3 = (_xx = -19266.439665524453, _yy = -8257.046893022618,
_zz = -8257.0459675641214, _xy = -64102.572896076854, _yz = 0, _zx = 0)
2 linear residual:
(const SymmTensor) $3 = (_xx = -19266.439665524453, _yy = -8257.046893022618,
_zz = -8257.0459675641214, _xy = -64102.572896076854, _yz = 0, _zx = 0)
3 linear residual:
4 linear residual:
(const SymmTensor) $7 = (_xx = -19266.425004289489, _yy = -8257.0393552307814,
_zz = -8257.039307856081, _xy = -64102.572529606805, _yz = 0, _zx = 0)

1st non-linear residual:
(lldb) p _stress[_qp]
(const SymmTensor) $2 = (_xx = -0.000046917203695813368, _yy =
-0.000010535390001637372, _zz = -0.000017235778109235224, _xy =
0.000037176401113541834, _yz = 0, _zx = 0)

So the big change happens in the first non-linear residual
evaluation. Everything is equivalent heading into that, e.g. the first linear
solve appears to be the same (potentially

Matrix:
strain:
(SymmTensor) $100 = (_xx = -0.00020484559886970247, _yy = 0, _zz = 0, _xy =
-0.0011927148071410806, _yz = 0, _zx = 0)

leads to:
stress:
(SymmTensor) $101 = (_xx = -275.75369078613789, _yy = -118.1801531940591, _zz =
-118.1801531940591, _xy = -917.47292857006209, _yz = 0, _zx = 0)

Matrix-free:
strain:
(SymmTensor) $10 = (_xx = -0.000000000038585853262551595, _yy =
0.0000000000087105045398772063, _zz = 0, _xy = 0.000000000048329321447604379,
_yz = 0, _zx = 0)

leads to:
stress:
(SymmTensor) $11 = (_xx = -0.000046917203695813368, _yy =
-0.000010535390001637372, _zz = -0.000017235778109235224, _xy =
0.000037176401113541834, _yz = 0, _zx = 0)

So our problem is indeed the solution vector obtained from incrementing the
solution vector with the linear solve. So something isn't good in the solution
of Ax = b for matrix-free. b is the same, so the problem must be in the
approximation of the Jacobian. Something in the linear solve for the explicit
matrix case doesn't allow the displacements on the LHS of the right block to be
equal to the displacements on the RHS of the right block, whereas this doesn't
happen in the matrix-free case.

So if I solve with NEWTON using an FDP preconditioner, then the solve looks
almost exactly like matrix free for the first non-linear iteration, e.g. the
solution vector after the first non-linear iteration looks the exact same and
consequently the first non-linear esidual is the exact same.. So
I've been barking up the wrong tree; I've been scrutinizing the first non-linear
iteration when I should be scrutinizing whats going wrong in the second
non-linaer iteration. The residual SHOULD increase from 0th non-linear to 1st
non-linear residual evaluation because we're moving into contact. So why is the
second linear solve not giving us a good newton step??? Why is the Jacobian
action approximated with matrix-free totally wrong?

With explicit matrix:
Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 4.786440e-02
      0 Linear |R| = 4.786440e-02
      1 Linear |R| = 1.328504e-17
  Linear solve converged due to CONVERGED_RTOL iterations 1
 1 Nonlinear |R| = 9.010475e-02
      0 Linear |R| = 9.010475e-02
      1 Linear |R| = 4.665225e-19
  Linear solve converged due to CONVERGED_RTOL iterations 1
 2 Nonlinear |R| = 4.026441e-04
      0 Linear |R| = 4.026441e-04
      1 Linear |R| = 6.858344e-19
  Linear solve converged due to CONVERGED_RTOL iterations 1
 3 Nonlinear |R| = 4.471554e-05
      0 Linear |R| = 4.471554e-05
      1 Linear |R| = 1.137866e-20
  Linear solve converged due to CONVERGED_RTOL iterations 1
 4 Nonlinear |R| = 4.297539e-10
Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 4
 Solve Converged!

With matrix free:
Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 4.786440e-02
      0 Linear |R| = 4.786440e-02
      1 Linear |R| = 2.696862e-11
  Linear solve converged due to CONVERGED_RTOL iterations 1
 1 Nonlinear |R| = 9.010475e-02
      0 Linear |R| = 9.010475e-02
      1 Linear |R| = 2.257897e-03
      2 Linear |R| = 1.125954e-03
      3 Linear |R| = 1.763171e-10
  Linear solve converged due to CONVERGED_RTOL iterations 3
 2 Nonlinear |R| = 3.450321e-01
      0 Linear |R| = 3.450321e-01
      1 Linear |R| = 1.499698e-02
      2 Linear |R| = 1.316531e-02
      3 Linear |R| = 3.271140e-03
      4 Linear |R| = 1.858556e-05
      5 Linear |R| = 2.756099e-08
  Linear solve converged due to CONVERGED_RTOL iterations 5
 3 Nonlinear |R| = 4.964049e-01
      0 Linear |R| = 4.964049e-01
      1 Linear |R| = 2.858225e-01
      2 Linear |R| = 9.975366e-02
      3 Linear |R| = 2.184165e-02
      4 Linear |R| = 4.251682e-03
      5 Linear |R| = 2.748165e-03
      6 Linear |R| = 6.535339e-04
      7 Linear |R| = 4.179196e-04
      8 Linear |R| = 1.498660e-05
      9 Linear |R| = 7.088401e-08
  Linear solve converged due to CONVERGED_RTOL iterations 9
 4 Nonlinear |R| = 1.540996e+00

Here's the key again, the true residual norm doesn't match the unpreconditioned
norm. Why??? We should not be hitting any discontinuties!

Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 4.786440e-02
    0 KSP unpreconditioned resid norm 4.786440052731e-02 true resid norm 4.786440052731e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.696861560766e-11 true resid norm 4.587005503695e-03 ||r(i)||/||b|| 9.583334280094e-02
  Linear solve converged due to CONVERGED_RTOL iterations 1
 1 Nonlinear |R| = 9.010475e-02
    0 KSP unpreconditioned resid norm 9.010474988004e-02 true resid norm 9.010474988004e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.257896880153e-03 true resid norm 1.978756666570e-01 ||r(i)||/||b|| 2.196062548539e+00
    2 KSP unpreconditioned resid norm 1.125953792339e-03 true resid norm 2.813551715533e-01 ||r(i)||/||b|| 3.122534294006e+00
    3 KSP unpreconditioned resid norm 1.781365914843e-10 true resid norm 3.404908043026e-01 ||r(i)||/||b|| 3.778833022187e+00
  Linear solve converged due to CONVERGED_RTOL iterations 3
 2 Nonlinear |R| = 3.450320e-01

More obvious if we look at pc_type none

Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 4.786440e-02
    0 KSP unpreconditioned resid norm 4.786440058172e-02 true resid norm 4.786440058172e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.751252206702e-02 true resid norm 1.751252206702e-02 ||r(i)||/||b|| 3.658778101090e-01
    2 KSP unpreconditioned resid norm 7.106912253736e-03 true resid norm 7.106912253468e-03 ||r(i)||/||b|| 1.484801265052e-01
    3 KSP unpreconditioned resid norm 4.838598438201e-03 true resid norm 4.838598433259e-03 ||r(i)||/||b|| 1.010897112353e-01
    4 KSP unpreconditioned resid norm 9.402733137094e-18 true resid norm 4.587004841479e-03 ||r(i)||/||b|| 9.583332885674e-02
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 9.010475e-02
    0 KSP unpreconditioned resid norm 9.010475030347e-02 true resid norm 9.010475030347e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.203968611083e-03 true resid norm 2.055972259145e-01 ||r(i)||/||b|| 2.281757900910e+00
    2 KSP unpreconditioned resid norm 3.047039944308e-03 true resid norm 2.764593916214e-01 ||r(i)||/||b|| 3.068199963823e+00
    3 KSP unpreconditioned resid norm 1.987903717965e-03 true resid norm 2.977504181431e-01 ||r(i)||/||b|| 3.304491906811e+00
    4 KSP unpreconditioned resid norm 1.520024041598e-03 true resid norm 3.828072914150e-01 ||r(i)||/||b|| 4.248469588182e+00
    5 KSP unpreconditioned resid norm 9.584383261815e-04 true resid norm 4.217813444189e-01 ||r(i)||/||b|| 4.681011189736e+00
    6 KSP unpreconditioned resid norm 6.171790593690e-04 true resid norm 4.678348303627e-01 ||r(i)||/||b|| 5.192121711531e+00
    7 KSP unpreconditioned resid norm 9.438734016357e-05 true resid norm 4.801062384168e-01 ||r(i)||/||b|| 5.328312178879e+00
    8 KSP unpreconditioned resid norm 3.608431834620e-05 true resid norm 4.834431947757e-01 ||r(i)||/||b|| 5.365346367949e+00
    9 KSP unpreconditioned resid norm 1.250228954521e-05 true resid norm 4.805488751314e-01 ||r(i)||/||b|| 5.333224647012e+00
   10 KSP unpreconditioned resid norm 1.707134184964e-15 true resid norm 6.036211965558e-01 ||r(i)||/||b|| 6.699105147319e+00
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 6.167561e-01

The true residual track the unpreconditioned residual until we go into a contact
state. (Note that that last true resid norm is evaulated with an updated
geometry search (and contact set), so its norm is expected to be wrong).

# 1/4/17

 0 Nonlinear |R| = 5.010181e-02
    0 KSP unpreconditioned resid norm 5.010180817294e-02 true resid norm 5.010180817294e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.100952035824e-02 true resid norm 2.100952035824e-02 ||r(i)||/||b|| 4.193365693653e-01
    2 KSP unpreconditioned resid norm 7.180877101910e-03 true resid norm 7.180877079432e-03 ||r(i)||/||b|| 1.433257070213e-01
    3 KSP unpreconditioned resid norm 5.094244660769e-03 true resid norm 5.094244572447e-03 ||r(i)||/||b|| 1.016778587085e-01
    4 KSP unpreconditioned resid norm 2.134369916921e-17 true resid norm 3.682509044977e-10 ||r(i)||/||b|| 7.350052182280e-09
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 9.780921e-03
    0 KSP unpreconditioned resid norm 9.780920542649e-03 true resid norm 9.780920542649e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.707562241232e-03 true resid norm 4.707562241232e-03 ||r(i)||/||b|| 4.813005300171e-01
    2 KSP unpreconditioned resid norm 1.741931637020e-03 true resid norm 1.741931670773e-03 ||r(i)||/||b|| 1.780948595971e-01
    3 KSP unpreconditioned resid norm 8.154433806799e-04 true resid norm 8.154436929076e-04 ||r(i)||/||b|| 8.337085342345e-02
    4 KSP unpreconditioned resid norm 6.103517240207e-04 true resid norm 6.103520745013e-04 ||r(i)||/||b|| 6.240231395807e-02
    5 KSP unpreconditioned resid norm 1.544604879250e-04 true resid norm 1.544610909317e-04 ||r(i)||/||b|| 1.579208115005e-02
    6 KSP unpreconditioned resid norm 9.798299082952e-05 true resid norm 9.798299531208e-05 ||r(i)||/||b|| 1.001776825451e-02
    7 KSP unpreconditioned resid norm 2.768979118017e-05 true resid norm 2.768946052242e-05 ||r(i)||/||b|| 2.830966717466e-03
    8 KSP unpreconditioned resid norm 2.286737086037e-06 true resid norm 2.286566300184e-06 ||r(i)||/||b|| 2.337782308131e-04
    9 KSP unpreconditioned resid norm 1.268811215870e-06 true resid norm 1.268442016853e-06 ||r(i)||/||b|| 1.296853411008e-04
   10 KSP unpreconditioned resid norm 9.200118608275e-18 true resid norm 2.074396036643e-09 ||r(i)||/||b|| 2.120859716217e-07
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 9.391402e-04
    0 KSP unpreconditioned resid norm 9.391402399113e-04 true resid norm 9.391402399113e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.487863277973e-04 true resid norm 1.487863277973e-04 ||r(i)||/||b|| 1.584282319873e-01
    2 KSP unpreconditioned resid norm 9.185801815663e-05 true resid norm 9.185801130834e-05 ||r(i)||/||b|| 9.781075009310e-02
    3 KSP unpreconditioned resid norm 6.422082037682e-05 true resid norm 6.422082556362e-05 ||r(i)||/||b|| 6.838257252152e-02
    4 KSP unpreconditioned resid norm 3.136563451299e-05 true resid norm 3.136564413463e-05 ||r(i)||/||b|| 3.339825385140e-02
    5 KSP unpreconditioned resid norm 1.620454190575e-05 true resid norm 1.620455514186e-05 ||r(i)||/||b|| 1.725467023262e-02
    6 KSP unpreconditioned resid norm 7.187043841688e-06 true resid norm 7.187066273756e-06 ||r(i)||/||b|| 7.652814743019e-03
    7 KSP unpreconditioned resid norm 1.783166346704e-06 true resid norm 1.783186663035e-06 ||r(i)||/||b|| 1.898743752268e-03
    8 KSP unpreconditioned resid norm 4.500388224186e-07 true resid norm 4.500352683525e-07 ||r(i)||/||b|| 4.791992177813e-04
    9 KSP unpreconditioned resid norm 3.432940501253e-08 true resid norm 3.431756311392e-08 ||r(i)||/||b|| 3.654146809550e-05
   10 KSP unpreconditioned resid norm 3.990716204447e-18 true resid norm 3.070023559057e-10 ||r(i)||/||b|| 3.268972437330e-07
  Linear solve converged due to CONVERGED_RTOL iterations 10
 3 Nonlinear |R| = 9.920930e-06
    0 KSP unpreconditioned resid norm 9.920930228057e-06 true resid norm 9.920930228057e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.354143086409e-06 true resid norm 1.354143086409e-06 ||r(i)||/||b|| 1.364935601078e-01
    2 KSP unpreconditioned resid norm 1.058760684754e-06 true resid norm 7.308929589806e-02 ||r(i)||/||b|| 7.367181727713e+03
    3 KSP unpreconditioned resid norm 7.324971791316e-07 true resid norm 2.605845152499e-01 ||r(i)||/||b|| 2.626613727339e+04
    4 KSP unpreconditioned resid norm 4.779488414455e-07 true resid norm 4.663984139039e-01 ||r(i)||/||b|| 4.701156072895e+04
    5 KSP unpreconditioned resid norm 1.469832541207e-07 true resid norm 5.976716947478e-01 ||r(i)||/||b|| 6.024351356263e+04
    6 KSP unpreconditioned resid norm 7.727833418229e-08 true resid norm 6.144981858583e-01 ||r(i)||/||b|| 6.193957337997e+04
    7 KSP unpreconditioned resid norm 2.620220915067e-08 true resid norm 6.199283005275e-01 ||r(i)||/||b|| 6.248691264598e+04
    8 KSP unpreconditioned resid norm 1.120793084916e-08 true resid norm 6.214574004311e-01 ||r(i)||/||b|| 6.264104132832e+04
    9 KSP unpreconditioned resid norm 2.857862878784e-09 true resid norm 6.222874212729e-01 ||r(i)||/||b|| 6.272470493876e+04
   10 KSP unpreconditioned resid norm 1.209168911852e-10 true resid norm 6.222750248087e-01 ||r(i)||/||b|| 6.272345541237e+04
   11 KSP unpreconditioned resid norm 8.502459988594e-11 true resid norm 6.222761509805e-01 ||r(i)||/||b|| 6.272356892711e+04
   12 KSP unpreconditioned resid norm 6.930921403032e-11 true resid norm 6.222770511739e-01 ||r(i)||/||b|| 6.272365966390e+04
   13 KSP unpreconditioned resid norm 5.996744265582e-11 true resid norm 6.222772359739e-01 ||r(i)||/||b|| 6.272367829118e+04
   14 KSP unpreconditioned resid norm 5.361039608480e-11 true resid norm 6.222775038739e-01 ||r(i)||/||b|| 6.272370529470e+04
   15 KSP unpreconditioned resid norm 4.892110804271e-11 true resid norm 6.222775774751e-01 ||r(i)||/||b|| 6.272371271347e+04
   16 KSP unpreconditioned resid norm 4.528165853780e-11 true resid norm 6.222777046707e-01 ||r(i)||/||b|| 6.272372553442e+04
   17 KSP unpreconditioned resid norm 4.234865180221e-11 true resid norm 6.222777440228e-01 ||r(i)||/||b|| 6.272372950099e+04
   18 KSP unpreconditioned resid norm 3.992130322828e-11 true resid norm 6.222778181060e-01 ||r(i)||/||b|| 6.272373696835e+04
   19 KSP unpreconditioned resid norm 3.786796857995e-11 true resid norm 6.222778425642e-01 ||r(i)||/||b|| 6.272373943367e+04
   20 KSP unpreconditioned resid norm 3.610250327373e-11 true resid norm 6.222778910090e-01 ||r(i)||/||b|| 6.272374431676e+04
   21 KSP unpreconditioned resid norm 3.456258841617e-11 true resid norm 6.222779076730e-01 ||r(i)||/||b|| 6.272374599644e+04
   22 KSP unpreconditioned resid norm 3.320459395800e-11 true resid norm 6.222779418108e-01 ||r(i)||/||b|| 6.272374943742e+04
   23 KSP unpreconditioned resid norm 3.199478278621e-11 true resid norm 6.222779538906e-01 ||r(i)||/||b|| 6.272375065503e+04
   24 KSP unpreconditioned resid norm 3.090846786518e-11 true resid norm 6.222779792383e-01 ||r(i)||/||b|| 6.272375321000e+04
   25 KSP unpreconditioned resid norm 2.992560214105e-11 true resid norm 6.222779883954e-01 ||r(i)||/||b|| 6.272375413301e+04
   26 KSP unpreconditioned resid norm 2.903106248891e-11 true resid norm 6.222780079584e-01 ||r(i)||/||b|| 6.272375610490e+04
   27 KSP unpreconditioned resid norm 2.821209149619e-11 true resid norm 6.222780151383e-01 ||r(i)||/||b|| 6.272375682862e+04
   28 KSP unpreconditioned resid norm 2.745885423792e-11 true resid norm 6.222780306929e-01 ||r(i)||/||b|| 6.272375839648e+04
   29 KSP unpreconditioned resid norm 2.676279381804e-11 true resid norm 6.222780364733e-01 ||r(i)||/||b|| 6.272375897912e+04
   30 KSP unpreconditioned resid norm 6.222780491364e-01 true resid norm 6.222780491364e-01 ||r(i)||/||b|| 6.272376025552e+04
  Linear solve did not converge due to DIVERGED_DTOL iterations 30
 4 Nonlinear |R| = 8.581282e-04

Ok, so in evaluation of 3rd non-linear residual, slave node 3 is in contact.

Ok, so preconditioner matrix is being evaluated in a different contact state
compared to the residual. Good god, how does this happen? It happens because of
the call the SNESComputeFunction for evaluating of the preconditioner
matrix. When evaluating the preceding residual, this order can happen:

- Node is evaluated to be in contact within
  `MechanicalContactConstraint::updateContactSet`.
- Later in the residual evaluation, the quantities that contact_pressure depends
  on are evaluated to be
  such that the node should be released. Specificially, the quantity that
  contact_pressure depends on is the _contact_force which is computed in ::shouldApply

Ok, during the non-linear residual eval:
contact_pressure = 18.2
nodalArea = .6
contact_force = _coords = ([0] = 10.945734695299096, [1] = 0.012390364283081646,
[2] = -0)

Now, during the preconditioner function evaluation (FDP):
contact_pressure = -1159
contact_force =    _coords = ([0] = -695.40033521053851, [1] =
-0.010404735891320266, [2] = 0)

Sure enought the contact force changed such that now the node should be out of
contact. But the contact force is based purely on the distance! At least for the
penalty formulation. I would hope that distance should be the same between the
nonlinear residual and preconditioner. Sure enough that distance is the same:

non-linear and preconditioner: -0.000069540033528837762

Back to badness caused by the existence of the FDP preconditioner:

Time Step  2, time = 0.2
                dt = 0.1
 0 Nonlinear |R| = 5.010181e-02
    0 KSP unpreconditioned resid norm 5.010180817294e-02 true resid norm 5.010180817294e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.100952035824e-02 true resid norm 2.100952035824e-02 ||r(i)||/||b|| 4.193365693653e-01
    2 KSP unpreconditioned resid norm 7.180877101910e-03 true resid norm 7.180877079432e-03 ||r(i)||/||b|| 1.433257070213e-01
    3 KSP unpreconditioned resid norm 5.094244660769e-03 true resid norm 5.094244572447e-03 ||r(i)||/||b|| 1.016778587085e-01
    4 KSP unpreconditioned resid norm 2.134369916921e-17 true resid norm 3.682509044977e-10 ||r(i)||/||b|| 7.350052182280e-09
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 9.780921e-03
    0 KSP unpreconditioned resid norm 9.780920542649e-03 true resid norm 9.780920542649e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.069792725991e-03 true resid norm 1.989876029154e-02 ||r(i)||/||b|| 2.034446574305e+00
    2 KSP unpreconditioned resid norm 3.136387629252e-03 true resid norm 1.079215927328e-01 ||r(i)||/||b|| 1.103388911731e+01
    3 KSP unpreconditioned resid norm 1.505880335793e-03 true resid norm 1.139030311714e-01 ||r(i)||/||b|| 1.164543057831e+01
    4 KSP unpreconditioned resid norm 9.614736733078e-04 true resid norm 1.354793784446e-01 ||r(i)||/||b|| 1.385139341986e+01
    5 KSP unpreconditioned resid norm 4.852203616420e-04 true resid norm 1.390276072277e-01 ||r(i)||/||b|| 1.421416385314e+01
    6 KSP unpreconditioned resid norm 2.496552751813e-04 true resid norm 1.579813432719e-01 ||r(i)||/||b|| 1.615199127557e+01
    7 KSP unpreconditioned resid norm 1.249831980081e-04 true resid norm 1.637468437986e-01 ||r(i)||/||b|| 1.674145527352e+01
    8 KSP unpreconditioned resid norm 5.073945616711e-06 true resid norm 1.664838636836e-01 ||r(i)||/||b|| 1.702128781823e+01
    9 KSP unpreconditioned resid norm 2.804853373524e-06 true resid norm 1.663961059497e-01 ||r(i)||/||b|| 1.701231547932e+01
   10 KSP unpreconditioned resid norm 2.109066321895e-17 true resid norm 2.996012248147e-01 ||r(i)||/||b|| 3.063118890582e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 3.003215e-01

    frame #8: 0x000000010f1753fd libpetsc.3.08.dylib`MatFDColoringApply_AIJ(J=0x00000001130a0260, coloring=0x0000000113114e60, x1=0x0000000113061e60, sctx=0x00000001130fec60) at fdmpiaij.c:268
    frame #9: 0x000000010ec352fa libpetsc.3.08.dylib`MatFDColoringApply(J=0x00000001130a0260, coloring=0x0000000113114e60, x1=0x0000000113061e60, sctx=0x00000001130fec60) at fdmatrix.c:629
    frame #10: 0x000000010fb6cf4e libpetsc.3.08.dylib`SNESComputeJacobianDefaultColor(snes=0x00000001130fec60, x1=0x0000000113061e60, J=0x0000000113125c60, B=0x00000001130a0260, ctx=0x0000000113114e60) at snesj2.c:100
    frame #11: 0x000000010fb49571 libpetsc.3.08.dylib`SNESComputeJacobian(snes=0x00000001130fec60, X=0x0000000113061e60, A=0x0000000113125c60, B=0x00000001130a0260) at snes.c:2358

    frame #8: 0x000000010f1753fd libpetsc.3.08.dylib`MatFDColoringApply_AIJ(J=0x00000001130a0260, coloring=0x0000000113114e60, x1=0x0000000113061e60, sctx=0x00000001130fec60) at fdmpiaij.c:268
    frame #9: 0x000000010ec352fa libpetsc.3.08.dylib`MatFDColoringApply(J=0x00000001130a0260, coloring=0x0000000113114e60, x1=0x0000000113061e60, sctx=0x00000001130fec60) at fdmatrix.c:629
    frame #10: 0x000000010fb6cf4e libpetsc.3.08.dylib`SNESComputeJacobianDefaultColor(snes=0x00000001130fec60, x1=0x0000000113061e60, J=0x0000000113125c60, B=0x00000001130a0260, ctx=0x0000000113114e60) at snesj2.c:100
    frame #11: 0x000000010fb49571 libpetsc.3.08.dylib`SNESComputeJacobian(snes=0x00000001130fec60, X=0x0000000113061e60, A=0x0000000113125c60, B=0x00000001130a0260) at snes.c:2358

There are 14 colors in the formation of the FDP so our residual function is
going to get called 14 times.

Here's the solve with standard finite differencing:

Time Step  2, time = 0.2
                dt = 0.1
 0 Nonlinear |R| = 5.010181e-02
    0 KSP unpreconditioned resid norm 5.010180817294e-02 true resid norm 5.010180817294e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.100952035824e-02 true resid norm 2.100952035824e-02 ||r(i)||/||b|| 4.193365693653e-01
    2 KSP unpreconditioned resid norm 7.180877101910e-03 true resid norm 7.180877079432e-03 ||r(i)||/||b|| 1.433257070213e-01
    3 KSP unpreconditioned resid norm 5.094244660769e-03 true resid norm 5.094244572447e-03 ||r(i)||/||b|| 1.016778587085e-01
    4 KSP unpreconditioned resid norm 2.134369916921e-17 true resid norm 3.682509044977e-10 ||r(i)||/||b|| 7.350052182280e-09
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 9.780921e-03
    0 KSP unpreconditioned resid norm 9.780920542650e-03 true resid norm 9.780920542650e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.038107260444e-03 true resid norm 4.038107260444e-03 ||r(i)||/||b|| 4.128555428741e-01
    2 KSP unpreconditioned resid norm 1.204558650367e-03 true resid norm 1.334555352196e-03 ||r(i)||/||b|| 1.364447596089e-01
    3 KSP unpreconditioned resid norm 8.675577590973e-04 true resid norm 1.166808216219e-03 ||r(i)||/||b|| 1.192943150014e-01
    4 KSP unpreconditioned resid norm 5.875116180079e-04 true resid norm 9.726598006835e-04 ||r(i)||/||b|| 9.944460712488e-02
    5 KSP unpreconditioned resid norm 2.903884497385e-04 true resid norm 7.971099790164e-04 ||r(i)||/||b|| 8.149641698249e-02
    6 KSP unpreconditioned resid norm 2.020233495814e-04 true resid norm 8.408918041266e-04 ||r(i)||/||b|| 8.597266488975e-02
    7 KSP unpreconditioned resid norm 1.656236789021e-05 true resid norm 8.832436431052e-04 ||r(i)||/||b|| 9.030271120738e-02
    8 KSP unpreconditioned resid norm 3.092322268614e-06 true resid norm 8.860899773040e-04 ||r(i)||/||b|| 9.059372003281e-02
    9 KSP unpreconditioned resid norm 5.638608218535e-07 true resid norm 8.846764226804e-04 ||r(i)||/||b|| 9.044919839833e-02
   10 KSP unpreconditioned resid norm 3.300608212365e-18 true resid norm 2.036073739605e-01 ||r(i)||/||b|| 2.081679051299e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 1.003715e-02

Ok, it looks to me like the contact state for the bad linear solve stays
constant between the 0th linear residual and the following. So why the
difference between resid norms and true resid norms?! Do we have if statements
elsewhere??

Ok, but what does it mean that the linear solves look perfect when just doing
JFNK vs. PJFNK??? Even during bad non-linear solves!:

Time Step  2, time = 0.2
                dt = 0.1
 0 Nonlinear |R| = 5.010181e-02
    0 KSP unpreconditioned resid norm 5.010180817294e-02 true resid norm 5.010180817294e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.100952035824e-02 true resid norm 2.100952035824e-02 ||r(i)||/||b|| 4.193365693653e-01
    2 KSP unpreconditioned resid norm 7.180877101910e-03 true resid norm 7.180877079432e-03 ||r(i)||/||b|| 1.433257070213e-01
    3 KSP unpreconditioned resid norm 5.094244660769e-03 true resid norm 5.094244572447e-03 ||r(i)||/||b|| 1.016778587085e-01
    4 KSP unpreconditioned resid norm 2.134369916921e-17 true resid norm 3.682509044977e-10 ||r(i)||/||b|| 7.350052182280e-09
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 9.780921e-03
    0 KSP unpreconditioned resid norm 9.780920542650e-03 true resid norm 9.780920542650e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.967175416004e-03 true resid norm 3.967175420154e-03 ||r(i)||/||b|| 4.056034810685e-01
    2 KSP unpreconditioned resid norm 1.084957047727e-03 true resid norm 1.084957064809e-03 ||r(i)||/||b|| 1.109258642965e-01
    3 KSP unpreconditioned resid norm 7.494431598547e-04 true resid norm 7.494431783962e-04 ||r(i)||/||b|| 7.662296970191e-02
    4 KSP unpreconditioned resid norm 5.165327162543e-04 true resid norm 5.165327265916e-04 ||r(i)||/||b|| 5.281023645364e-02
    5 KSP unpreconditioned resid norm 2.786339581119e-04 true resid norm 2.786339426304e-04 ||r(i)||/||b|| 2.848749679699e-02
    6 KSP unpreconditioned resid norm 1.839544405828e-04 true resid norm 1.839544046311e-04 ||r(i)||/||b|| 1.880747357358e-02
    7 KSP unpreconditioned resid norm 1.666837168406e-05 true resid norm 1.666840863062e-05 ||r(i)||/||b|| 1.704175855221e-03
    8 KSP unpreconditioned resid norm 2.797798021039e-06 true resid norm 2.797799298996e-06 ||r(i)||/||b|| 2.860466238117e-04
    9 KSP unpreconditioned resid norm 5.736600192002e-07 true resid norm 5.736480277787e-07 ||r(i)||/||b|| 5.864969716064e-05
   10 KSP unpreconditioned resid norm 4.283855947448e-18 true resid norm 1.742252277471e-01 ||r(i)||/||b|| 1.781276383827e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 8.454889e-03
    0 KSP unpreconditioned resid norm 8.454888633012e-03 true resid norm 8.454888633012e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.214043299665e-03 true resid norm 2.214043299665e-03 ||r(i)||/||b|| 2.618654598264e-01
    2 KSP unpreconditioned resid norm 6.401108036954e-04 true resid norm 6.401107974133e-04 ||r(i)||/||b|| 7.570895669920e-02
    3 KSP unpreconditioned resid norm 3.152150315767e-04 true resid norm 3.152150093325e-04 ||r(i)||/||b|| 3.728198241450e-02
    4 KSP unpreconditioned resid norm 1.427037751408e-04 true resid norm 1.427037677954e-04 ||r(i)||/||b|| 1.687825517159e-02
    5 KSP unpreconditioned resid norm 3.588127490001e-05 true resid norm 3.588127532490e-05 ||r(i)||/||b|| 4.243849550519e-03
    6 KSP unpreconditioned resid norm 2.330956273631e-05 true resid norm 2.330956669412e-05 ||r(i)||/||b|| 2.756933616264e-03
    7 KSP unpreconditioned resid norm 4.431430052986e-06 true resid norm 4.431431351152e-06 ||r(i)||/||b|| 5.241265194020e-04
    8 KSP unpreconditioned resid norm 2.432027735857e-06 true resid norm 2.432027608964e-06 ||r(i)||/||b|| 2.876475036547e-04
    9 KSP unpreconditioned resid norm 2.868464074042e-07 true resid norm 2.868401713706e-07 ||r(i)||/||b|| 3.392595500911e-05
   10 KSP unpreconditioned resid norm 6.736546911810e-19 true resid norm 2.532819056402e-11 ||r(i)||/||b|| 2.995685888177e-09
  Linear solve converged due to CONVERGED_RTOL iterations 10
 3 Nonlinear |R| = 3.203698e-05
    0 KSP unpreconditioned resid norm 3.203697821607e-05 true resid norm 3.203697821607e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.812632942113e-05 true resid norm 1.812632942113e-05 ||r(i)||/||b|| 5.657939802835e-01
    2 KSP unpreconditioned resid norm 6.447163950844e-06 true resid norm 6.447164085252e-06 ||r(i)||/||b|| 2.012413293716e-01
    3 KSP unpreconditioned resid norm 1.617719805954e-06 true resid norm 1.617719787206e-06 ||r(i)||/||b|| 5.049539242731e-02
    4 KSP unpreconditioned resid norm 9.120607741052e-07 true resid norm 9.120607053025e-07 ||r(i)||/||b|| 2.846899914065e-02
    5 KSP unpreconditioned resid norm 5.977828948563e-07 true resid norm 5.977828937677e-07 ||r(i)||/||b|| 1.865915348620e-02
    6 KSP unpreconditioned resid norm 4.809063940468e-07 true resid norm 4.809063435320e-07 ||r(i)||/||b|| 1.501097701190e-02
    7 KSP unpreconditioned resid norm 1.730044803222e-08 true resid norm 1.730043095362e-08 ||r(i)||/||b|| 5.400144432142e-04
    8 KSP unpreconditioned resid norm 2.426843194615e-09 true resid norm 2.426837711512e-09 ||r(i)||/||b|| 7.575114279331e-05
    9 KSP unpreconditioned resid norm 1.259473025760e-09 true resid norm 1.259496591476e-09 ||r(i)||/||b|| 3.931383862052e-05
   10 KSP unpreconditioned resid norm 2.275642503188e-21 true resid norm 2.872895034521e-13 ||r(i)||/||b|| 8.967434491309e-09
  Linear solve converged due to CONVERGED_RTOL iterations 10
 4 Nonlinear |R| = 1.577954e-09
Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 4
 Solve Converged!

Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 4.786440e-02
    0 KSP unpreconditioned resid norm 4.786440044274e-02 true resid norm 4.786440044274e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.891237203690e-02 true resid norm 1.891237203690e-02 ||r(i)||/||b|| 3.951239723462e-01
    2 KSP unpreconditioned resid norm 7.073418966111e-03 true resid norm 7.073418829106e-03 ||r(i)||/||b|| 1.477803704565e-01
    3 KSP unpreconditioned resid norm 4.981544070757e-03 true resid norm 4.981544025388e-03 ||r(i)||/||b|| 1.040761814482e-01
    4 KSP unpreconditioned resid norm 6.673626083825e-10 true resid norm 7.784211968510e-10 ||r(i)||/||b|| 1.626305123747e-08
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 1.271783e-01
    0 KSP unpreconditioned resid norm 1.271783334012e-01 true resid norm 1.271783334012e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.736459063019e-02 true resid norm 3.736459063019e-02 ||r(i)||/||b|| 2.937968255357e-01
    2 KSP unpreconditioned resid norm 1.305858762840e-02 true resid norm 1.305858747798e-02 ||r(i)||/||b|| 1.026793411169e-01
    3 KSP unpreconditioned resid norm 7.836847728700e-03 true resid norm 7.836847591321e-03 ||r(i)||/||b|| 6.162093323395e-02
    4 KSP unpreconditioned resid norm 4.188613422779e-03 true resid norm 4.188613839770e-03 ||r(i)||/||b|| 3.293496405993e-02
    5 KSP unpreconditioned resid norm 1.153289470970e-03 true resid norm 1.153290305285e-03 ||r(i)||/||b|| 9.068292329692e-03
    6 KSP unpreconditioned resid norm 8.025212320793e-04 true resid norm 8.025218247299e-04 ||r(i)||/||b|| 6.310208690958e-03
    7 KSP unpreconditioned resid norm 2.973707031036e-05 true resid norm 2.973528293238e-05 ||r(i)||/||b|| 2.338077732044e-04
    8 KSP unpreconditioned resid norm 3.156779306488e-06 true resid norm 3.156709709714e-06 ||r(i)||/||b|| 2.482112813789e-05
    9 KSP unpreconditioned resid norm 1.670984904396e-06 true resid norm 1.671661179129e-06 ||r(i)||/||b|| 1.314422932290e-05
   10 KSP unpreconditioned resid norm 2.117450310570e-17 true resid norm 2.425273915917e+00 ||r(i)||/||b|| 1.906986709966e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 2 Nonlinear |R| = 1.260783e-01
    0 KSP unpreconditioned resid norm 1.260782996476e-01 true resid norm 1.260782996476e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.366955419990e-02 true resid norm 3.366955419990e-02 ||r(i)||/||b|| 2.670527306762e-01
    2 KSP unpreconditioned resid norm 1.353999079301e-02 true resid norm 1.353999078789e-02 ||r(i)||/||b|| 1.073935072549e-01
    3 KSP unpreconditioned resid norm 7.319856638937e-03 true resid norm 7.319856438378e-03 ||r(i)||/||b|| 5.805801996726e-02
    4 KSP unpreconditioned resid norm 4.958089847979e-03 true resid norm 4.958090351989e-03 ||r(i)||/||b|| 3.932548555816e-02
    5 KSP unpreconditioned resid norm 3.702115889975e-03 true resid norm 3.702116073903e-03 ||r(i)||/||b|| 2.936362628819e-02
    6 KSP unpreconditioned resid norm 3.064744255318e-03 true resid norm 3.064744167288e-03 ||r(i)||/||b|| 2.430826062736e-02
    7 KSP unpreconditioned resid norm 7.662955657718e-04 true resid norm 7.662960904803e-04 ||r(i)||/||b|| 6.077938016472e-03
    8 KSP unpreconditioned resid norm 4.222673849240e-05 true resid norm 4.222638722895e-05 ||r(i)||/||b|| 3.349219282539e-04
    9 KSP unpreconditioned resid norm 6.621945447620e-06 true resid norm 6.621740338236e-06 ||r(i)||/||b|| 5.252085693371e-05
   10 KSP unpreconditioned resid norm 8.623299053733e-18 true resid norm 1.394581427480e-09 ||r(i)||/||b|| 1.106123283212e-08
  Linear solve converged due to CONVERGED_RTOL iterations 10
 3 Nonlinear |R| = 2.357292e-01
    0 KSP unpreconditioned resid norm 2.357291789929e-01 true resid norm 2.357291789929e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.112782919703e-02 true resid norm 7.112782919703e-02 ||r(i)||/||b|| 3.017353621682e-01
    2 KSP unpreconditioned resid norm 2.822979165942e-02 true resid norm 2.822979299478e-02 ||r(i)||/||b|| 1.197551915948e-01
    3 KSP unpreconditioned resid norm 1.731713201495e-02 true resid norm 1.731713320673e-02 ||r(i)||/||b|| 7.346198413246e-02
    4 KSP unpreconditioned resid norm 7.628780368795e-03 true resid norm 7.628781836014e-03 ||r(i)||/||b|| 3.236248422282e-02
    5 KSP unpreconditioned resid norm 2.690141226625e-03 true resid norm 2.690145224721e-03 ||r(i)||/||b|| 1.141201626466e-02
    6 KSP unpreconditioned resid norm 1.353140696655e-03 true resid norm 1.353142840795e-03 ||r(i)||/||b|| 5.740243301981e-03
    7 KSP unpreconditioned resid norm 5.875240595944e-04 true resid norm 5.875234525567e-04 ||r(i)||/||b|| 2.492366261431e-03
    8 KSP unpreconditioned resid norm 9.919935394522e-05 true resid norm 9.920138608739e-05 ||r(i)||/||b|| 4.208277758028e-04
    9 KSP unpreconditioned resid norm 4.367942585256e-05 true resid norm 4.367955358103e-05 ||r(i)||/||b|| 1.852954893732e-04
   10 KSP unpreconditioned resid norm 8.401017139199e-17 true resid norm 4.459684339840e+00 ||r(i)||/||b|| 1.891867760662e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 4 Nonlinear |R| = 2.296731e-01
    0 KSP unpreconditioned resid norm 2.296731456158e-01 true resid norm 2.296731456158e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.075249684699e-02 true resid norm 6.075249684699e-02 ||r(i)||/||b|| 2.645171976206e-01
    2 KSP unpreconditioned resid norm 2.635397336729e-02 true resid norm 2.635397358140e-02 ||r(i)||/||b|| 1.147455594373e-01
    3 KSP unpreconditioned resid norm 1.504095333293e-02 true resid norm 1.504095301432e-02 ||r(i)||/||b|| 6.548851400973e-02
    4 KSP unpreconditioned resid norm 9.512856267385e-03 true resid norm 9.512855885667e-03 ||r(i)||/||b|| 4.141910391901e-02
    5 KSP unpreconditioned resid norm 6.413688421769e-03 true resid norm 6.413687891200e-03 ||r(i)||/||b|| 2.792528431656e-02
    6 KSP unpreconditioned resid norm 5.749390489947e-03 true resid norm 5.749387983987e-03 ||r(i)||/||b|| 2.503291348482e-02
    7 KSP unpreconditioned resid norm 3.457466592195e-03 true resid norm 3.457466421535e-03 ||r(i)||/||b|| 1.505385582744e-02
    8 KSP unpreconditioned resid norm 6.419644917829e-04 true resid norm 6.419647743577e-04 ||r(i)||/||b|| 2.795123359487e-03
    9 KSP unpreconditioned resid norm 2.468424133191e-04 true resid norm 2.468420839741e-04 ||r(i)||/||b|| 1.074753791142e-03
   10 KSP unpreconditioned resid norm 8.481035170210e-17 true resid norm 8.150968147192e-09 ||r(i)||/||b|| 3.548942618144e-08
  Linear solve converged due to CONVERGED_RTOL iterations 10
 5 Nonlinear |R| = 9.079311e-01
    0 KSP unpreconditioned resid norm 9.079311460436e-01 true resid norm 9.079311460436e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.395284912123e-01 true resid norm 3.395286376917e-01 ||r(i)||/||b|| 3.739585751312e-01
    2 KSP unpreconditioned resid norm 1.545336691965e-01 true resid norm 1.545336617610e-01 ||r(i)||/||b|| 1.702041640871e-01
    3 KSP unpreconditioned resid norm 1.440840119720e-01 true resid norm 1.440839007780e-01 ||r(i)||/||b|| 1.586947439857e-01
    4 KSP unpreconditioned resid norm 9.570187556482e-02 true resid norm 9.570181447828e-02 ||r(i)||/||b|| 1.054064671042e-01
    5 KSP unpreconditioned resid norm 4.951640603417e-02 true resid norm 4.951638461747e-02 ||r(i)||/||b|| 5.453759884023e-02
    6 KSP unpreconditioned resid norm 2.559995286006e-02 true resid norm 2.559994728110e-02 ||r(i)||/||b|| 2.819591264454e-02
    7 KSP unpreconditioned resid norm 1.211332766387e-02 true resid norm 1.211335855915e-02 ||r(i)||/||b|| 1.334171496587e-02
    8 KSP unpreconditioned resid norm 3.556553490575e-03 true resid norm 3.556546721049e-03 ||r(i)||/||b|| 3.917198717708e-03
    9 KSP unpreconditioned resid norm 1.322291448801e-03 true resid norm 1.322274855435e-03 ||r(i)||/||b|| 1.456360277095e-03
   10 KSP unpreconditioned resid norm 1.305827110142e-16 true resid norm 1.485246541791e+01 ||r(i)||/||b|| 1.635858124554e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 6 Nonlinear |R| = 1.548943e+00
    0 KSP unpreconditioned resid norm 1.548943149115e+00 true resid norm 1.548943149115e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.283125653736e+00 true resid norm 1.283125471681e+00 ||r(i)||/||b|| 8.283877122370e-01
    2 KSP unpreconditioned resid norm 1.262343509587e+00 true resid norm 1.262343491852e+00 ||r(i)||/||b|| 8.149708351617e-01
    3 KSP unpreconditioned resid norm 1.133410552608e+00 true resid norm 1.133410225403e+00 ||r(i)||/||b|| 7.317313266469e-01
    4 KSP unpreconditioned resid norm 1.112094550746e+00 true resid norm 1.112094593713e+00 ||r(i)||/||b|| 7.179699231364e-01
    5 KSP unpreconditioned resid norm 1.071082663950e+00 true resid norm 1.071082971713e+00 ||r(i)||/||b|| 6.914927589983e-01
    6 KSP unpreconditioned resid norm 1.044699088118e+00 true resid norm 1.044699262080e+00 ||r(i)||/||b|| 6.744593968325e-01
    7 KSP unpreconditioned resid norm 9.288207605279e-01 true resid norm 9.288207875153e-01 ||r(i)||/||b|| 5.996480813682e-01
    8 KSP unpreconditioned resid norm 3.734238266703e-01 true resid norm 3.734237453796e-01 ||r(i)||/||b|| 2.410829252145e-01
    9 KSP unpreconditioned resid norm 1.500324154127e-01 true resid norm 1.500328216301e-01 ||r(i)||/||b|| 9.686141270955e-02
   10 KSP unpreconditioned resid norm 6.849247046622e-15 true resid norm 1.330873222781e+02 ||r(i)||/||b|| 8.592137313377e+01
  Linear solve converged due to CONVERGED_RTOL iterations 10
 7 Nonlinear |R| = 1.530024e+02

So I really want to figure out wh the solve paths diverge in time step 2, second
linear solve. Where does the difference come from?

Here's the residual vector after kernels for JFNK:
  2.05164e-09
  5.28074e-10
  -4.40357e-09
  2.91356e-09
  -7.85051e-09
  7.42882e-11
  1.26407e-08
  2.27657e-10
  -4.51141e-09
  -3.07412e-09
  2.07321e-09
  -6.69461e-10
  -0.00103083
  0.00140039
  0.00118094
  -0.0015735
  -0.00127283
  -0.00309152
  0.00112272
  0.00326463

After contraint residuals?:
  2.05164e-09
  5.28074e-10
  -4.40357e-09
  2.91356e-09
  -7.85051e-09
  7.42882e-11
  0.00689478
  0.000124175
  -4.51141e-09
  -3.07412e-09
  2.07321e-09
  -6.69461e-10
  -0.00103083
  0.00140039
  0.00118094
  -0.0015735
  -0.00476468
  -0.0031544
  -0.0022802
  0.00320334

With PJFNK after kernels:
  -2.05164e-09
  -5.28074e-10
  4.40357e-09
  -2.91356e-09
  7.85051e-09
  -7.42881e-11
  -1.26407e-08
  -2.27657e-10
  4.51141e-09
  3.07412e-09
  -2.07321e-09
  6.69461e-10
  -0.00103084
  0.00140039
  0.00118094
  -0.0015735
  -0.00127282
  -0.00309151
  0.00112273
  0.00326463

After constraint residuals:
  -2.05164e-09
  -5.28074e-10
  4.40357e-09
  -2.91356e-09
  7.85051e-09
  -7.42881e-11
  0.00689467
  0.000124172
  4.51141e-09
  3.07412e-09
  -2.07321e-09
  6.69461e-10
  -0.00103084
  0.00140039
  0.00118094
  -0.0015735
  -0.00476463
  -0.0031544
  -0.00228015
  0.00320334

Ok, JFNK:
U:
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  -0.00869424
  -0.2
  -0.00869424
  -0.2
  -0.017699
  -0.21081
  -0.0034968
  -0.199386
h: 0.0000016232913134256937
a:
  0.
  0.
  0.
  0.
  0.
  0.
  0.0076723
  0.000138177
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  -0.00530199
  -0.00351012
  -0.00253733
  0.00356458

PJFNK:
U:
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  -0.00869424
  -0.2
  -0.00869424
  -0.2
  -0.017699
  -0.21081
  -0.0034968
  -0.199386
h: 0.000011487644520672412
a:
  0.
  0.
  0.
  0.
  0.
  0.
  -0.00108415
  -1.95255e-05
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.
  0.000749212
  0.000496007
  0.000358544
  -0.000503703

For GMRES JFNK:
h: 0.000000017667764620572558
a:
0.
0.
0.
0.
0.
0.
0.704921
0.0126956
0.
0.
0.
0.
0.
0.
0.
0.
-0.48714
-0.322506
-0.233127
0.327509

for GMRES PJFNK:
h: 0.000000017667764620572562
a:
0.
0.
0.
0.
0.
0.
0.704921
0.0126956
0.
0.
0.
0.
0.
0.
0.
0.
-0.48714
-0.322506
-0.233127
0.327509

So this looks the exact same. Now the difference must come from the residual evaluation
GMRES, PJFNK:
after kernels:
2.05164e-09
5.28074e-10
-4.40357e-09
2.91356e-09
-7.85051e-09
7.42882e-11
1.26407e-08
2.27657e-10
-4.51141e-09
-3.07412e-09
2.07321e-09
-6.69461e-10
-0.00103083
0.00140039
0.00118094
-0.0015735
-0.00127283
-0.00309152
0.00112272
0.00326463

after constraints:
2.05164e-09
5.28074e-10
-4.40357e-09
2.91356e-09
-7.85051e-09
7.42882e-11
0.00689469
0.000124173
-4.51141e-09
-3.07412e-09
2.07321e-09
-6.69461e-10
-0.00103083
0.00140039
0.00118094
-0.0015735
-0.00476464
-0.0031544
-0.00228015
0.00320335

after nodal bcs:
2.05164e-09
5.28074e-10
0.
0.
0.
0.
0.00689469
0.000124173
0.
0.
2.07321e-09
-6.69461e-10
0.
0.
0.
0.
-0.00476464
-0.0031544
-0.00228015
0.00320335

GMRES, JFNK:
after kernels:
2.05164e-09
5.28074e-10
-4.40357e-09
2.91356e-09
-7.85051e-09
7.42882e-11
1.26407e-08
2.27657e-10
-4.51141e-09
-3.07412e-09
2.07321e-09
-6.69461e-10
-0.00103083
0.00140039
0.00118094
-0.0015735
-0.00127283
-0.00309152
0.00112272
0.00326463

after constraints:
2.05164e-09
5.28074e-10
-4.40357e-09
2.91356e-09
-7.85051e-09
7.42882e-11
0.00689478
0.000124175
-4.51141e-09
-3.07412e-09
2.07321e-09
-6.69461e-10
-0.00103083
0.00140039
0.00118094
-0.0015735
-0.00476468
-0.0031544
-0.0022802
0.00320334

after nodal bcs:
2.05164e-09
5.28074e-10
0.
0.
0.
0.
0.00689478
0.000124175
0.
0.
2.07321e-09
-6.69461e-10
0.
0.
0.
0.
-0.00476468
-0.0031544
-0.0022802
0.00320334


JFNK GMRES:
master resid 1: -3491.8495368255708
master resid 2: -3402.9222261716877
slave resid: 6894.7717629972594

PJFNK GMRES:
master resid 1: -3491.8042367930398
master resid 2: -3402.8780798000266
slave resid: 6894.6823165930682

They are all very slightly different! Why!?

PJFNK:
contact force:
    _coords = ([0] = 6894.6823165930682, [1] = 124.17253658544793, [2] = -0)

JFNK:
contact force:
    _coords = ([0] = 6894.7717629972594, [1] = 124.1743019811705, [2] = -0)

What I thought! Now let's look at the contact force for PJFNK.

PJFNK, non-linear resid eval:
    _coords = ([0] = 6894.7717629967046, [1] = 124.1743019811605, [2] = -0)

Ok, now see that this is identical to the JFNK eval

No coloring (standard differencing):
master resid 1: -3491.8495368255708
master resid 2: -3402.9222261716877
slave resid: 6894.7717629972594

Coloring:
slave resid: 6894.6823165930682

Looking now at this `fdp_geometric_coupling.i` test in combined that's
failing. Looks PJFNK again screws with the problem state.

Purely mechanical contact, things look great:

Time Step  0, time = 0
                dt = 0

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
+----------------+------------------+----------------+----------------+


Time Step  1, time = 0.1
                dt = 0.1
 0 Nonlinear |R| = 7.810169e+04
    0 KSP unpreconditioned resid norm 7.810169371831e+04 true resid norm 7.810169371831e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.310217327351e+04 true resid norm 1.310217327351e+04 ||r(i)||/||b|| 1.677578634948e-01
    2 KSP unpreconditioned resid norm 1.331047225422e+03 true resid norm 1.331050630500e+03 ||r(i)||/||b|| 1.704253220553e-02
    3 KSP unpreconditioned resid norm 1.972062490300e+02 true resid norm 1.972078304974e+02 ||r(i)||/||b|| 2.525013493416e-03
    4 KSP unpreconditioned resid norm 3.650798933531e+01 true resid norm 3.650814110670e+01 ||r(i)||/||b|| 4.674436541463e-04
    5 KSP unpreconditioned resid norm 5.576140484585e+00 true resid norm 5.576286957170e+00 ||r(i)||/||b|| 7.139777246422e-05
    6 KSP unpreconditioned resid norm 1.250213945993e+00 true resid norm 1.249942105507e+00 ||r(i)||/||b|| 1.600403328019e-05
    7 KSP unpreconditioned resid norm 2.313040717716e-01 true resid norm 2.306906907994e-01 ||r(i)||/||b|| 2.953721997777e-06
    8 KSP unpreconditioned resid norm 3.178315530844e-02 true resid norm 3.118655227568e-02 ||r(i)||/||b|| 3.993069905522e-07
      Line search: Using full step: fnorm 7.810169371831e+04 gnorm 1.369964438153e+04
 1 Nonlinear |R| = 1.369964e+04
    0 KSP unpreconditioned resid norm 1.369964438153e+04 true resid norm 1.369964438153e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.734334955736e+02 true resid norm 4.734334955736e+02 ||r(i)||/||b|| 3.455808650128e-02
    2 KSP unpreconditioned resid norm 1.681122063352e+01 true resid norm 1.681119525283e+01 ||r(i)||/||b|| 1.227126397201e-03
    3 KSP unpreconditioned resid norm 3.597620963660e-01 true resid norm 3.597665178783e-01 ||r(i)||/||b|| 2.626101144373e-05
    4 KSP unpreconditioned resid norm 6.848687716769e-03 true resid norm 6.838312971503e-03 ||r(i)||/||b|| 4.991598891957e-07
      Line search: Using full step: fnorm 1.369964438153e+04 gnorm 3.015701557796e+02
 2 Nonlinear |R| = 3.015702e+02
    0 KSP unpreconditioned resid norm 3.015701557796e+02 true resid norm 3.015701557796e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.658342583353e-01 true resid norm 4.658342583353e-01 ||r(i)||/||b|| 1.544696149163e-03
    2 KSP unpreconditioned resid norm 3.970594044122e-04 true resid norm 3.968563599361e-04 ||r(i)||/||b|| 1.315966956048e-06
    3 KSP unpreconditioned resid norm 1.848139672787e-07 true resid norm 1.181451009799e-06 ||r(i)||/||b|| 3.917665548651e-09
      Line search: Using full step: fnorm 3.015701557796e+02 gnorm 1.689801705013e-01
 3 Nonlinear |R| = 1.689802e-01
    0 KSP unpreconditioned resid norm 1.689801705013e-01 true resid norm 1.689801705013e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.900098605253e-08 true resid norm 5.900098605321e-08 ||r(i)||/||b|| 3.491592290277e-07
      Line search: Using full step: fnorm 1.689801705013e-01 gnorm 7.142121341985e-08
 4 Nonlinear |R| = 7.142121e-08
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
+----------------+------------------+----------------+----------------+


Time Step  2, time = 0.2
                dt = 0.1
 0 Nonlinear |R| = 8.727455e+04
    0 KSP unpreconditioned resid norm 8.727454821989e+04 true resid norm 8.727454821989e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.248222334089e+04 true resid norm 1.248222334089e+04 ||r(i)||/||b|| 1.430224916139e-01
    2 KSP unpreconditioned resid norm 2.039264817611e+03 true resid norm 2.039266866368e+03 ||r(i)||/||b|| 2.336611197608e-02
    3 KSP unpreconditioned resid norm 3.346565755800e+02 true resid norm 3.346560076194e+02 ||r(i)||/||b|| 3.834520079969e-03
    4 KSP unpreconditioned resid norm 7.431809944698e+01 true resid norm 7.431741309467e+01 ||r(i)||/||b|| 8.515359243960e-04
    5 KSP unpreconditioned resid norm 1.557955688985e+01 true resid norm 1.557859851546e+01 ||r(i)||/||b|| 1.785010502284e-04
    6 KSP unpreconditioned resid norm 2.500361622146e+00 true resid norm 2.500331950144e+00 ||r(i)||/||b|| 2.864903916597e-05
    7 KSP unpreconditioned resid norm 4.558816667560e-01 true resid norm 4.536181808609e-01 ||r(i)||/||b|| 5.197599874341e-06
    8 KSP unpreconditioned resid norm 9.216191322389e-02 true resid norm 9.286827744607e-02 ||r(i)||/||b|| 1.064093476738e-06
    9 KSP unpreconditioned resid norm 1.066610374218e-02 true resid norm 1.180456652520e-02 ||r(i)||/||b|| 1.352578359439e-07
      Line search: gnorm after quadratic fit 5.420141560003e+04
      Line search: Quadratically determined step, lambda=3.8748876266726073e-01
 1 Nonlinear |R| = 5.420142e+04
    0 KSP unpreconditioned resid norm 5.420141560003e+04 true resid norm 5.420141560003e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.746566194782e+03 true resid norm 4.746566194782e+03 ||r(i)||/||b|| 8.757273481210e-02
    2 KSP unpreconditioned resid norm 3.493715508695e+02 true resid norm 3.493711233469e+02 ||r(i)||/||b|| 6.445793333609e-03
    3 KSP unpreconditioned resid norm 3.867132949240e+01 true resid norm 3.867161444075e+01 ||r(i)||/||b|| 7.134797866926e-04
    4 KSP unpreconditioned resid norm 6.329844737937e+00 true resid norm 6.329594944747e+00 ||r(i)||/||b|| 1.167791445053e-04
    5 KSP unpreconditioned resid norm 6.318734701698e-01 true resid norm 6.315739357622e-01 ||r(i)||/||b|| 1.165235130430e-05
    6 KSP unpreconditioned resid norm 6.699263719962e-02 true resid norm 6.756977708788e-02 ||r(i)||/||b|| 1.246642294115e-06
    7 KSP unpreconditioned resid norm 4.300849748738e-03 true resid norm 5.042722847290e-03 ||r(i)||/||b|| 9.303673698307e-08
      Line search: Using full step: fnorm 5.420141560003e+04 gnorm 1.756971896169e+04
 2 Nonlinear |R| = 1.756972e+04
    0 KSP unpreconditioned resid norm 1.756971896169e+04 true resid norm 1.756971896169e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.301494023041e+02 true resid norm 1.301494018305e+02 ||r(i)||/||b|| 7.407597248098e-03
    2 KSP unpreconditioned resid norm 3.087430481221e+00 true resid norm 3.087436805151e+00 ||r(i)||/||b|| 1.757248827874e-04
    3 KSP unpreconditioned resid norm 3.371776116793e-01 true resid norm 3.372089077713e-01 ||r(i)||/||b|| 1.919261819194e-05
    4 KSP unpreconditioned resid norm 1.247584341799e-02 true resid norm 1.251249000661e-02 ||r(i)||/||b|| 7.121622169310e-07
      Line search: Using full step: fnorm 1.756971896169e+04 gnorm 6.592554696635e+02
 3 Nonlinear |R| = 6.592555e+02
    0 KSP unpreconditioned resid norm 6.592554696635e+02 true resid norm 6.592554696635e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.526893735068e+00 true resid norm 1.526893735068e+00 ||r(i)||/||b|| 2.316088080161e-03
    2 KSP unpreconditioned resid norm 3.300523189947e-03 true resid norm 3.300110098435e-03 ||r(i)||/||b|| 5.005813755508e-06
    3 KSP unpreconditioned resid norm 6.696863197380e-06 true resid norm 1.371370646549e-05 ||r(i)||/||b|| 2.080180915676e-08
      Line search: Using full step: fnorm 6.592554696635e+02 gnorm 5.352465197869e-01
 4 Nonlinear |R| = 5.352465e-01
    0 KSP unpreconditioned resid norm 5.352465197869e-01 true resid norm 5.352465197869e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.740999983212e-06 true resid norm 5.740999983214e-06 ||r(i)||/||b|| 1.072589876063e-05
    2 KSP unpreconditioned resid norm 1.273218639849e-11 true resid norm 1.214703692494e-08 ||r(i)||/||b|| 2.269428473776e-08
      Line search: Using full step: fnorm 5.352465197869e-01 gnorm 1.110845676362e-06
 5 Nonlinear |R| = 1.110846e-06
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
+----------------+------------------+----------------+----------------+


Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 8.270134e+04
    0 KSP unpreconditioned resid norm 8.270134234182e+04 true resid norm 8.270134234182e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.249175980792e+04 true resid norm 1.249175980792e+04 ||r(i)||/||b|| 1.510466390773e-01
    2 KSP unpreconditioned resid norm 1.760270548763e+03 true resid norm 1.760271689411e+03 ||r(i)||/||b|| 2.128468099266e-02
    3 KSP unpreconditioned resid norm 2.646275002922e+02 true resid norm 2.646285033783e+02 ||r(i)||/||b|| 3.199809046442e-03
    4 KSP unpreconditioned resid norm 6.526242156694e+01 true resid norm 6.526006630234e+01 ||r(i)||/||b|| 7.891052848043e-04
    5 KSP unpreconditioned resid norm 1.414141196760e+01 true resid norm 1.414261211904e+01 ||r(i)||/||b|| 1.710082535369e-04
    6 KSP unpreconditioned resid norm 2.396887290035e+00 true resid norm 2.397719122870e+00 ||r(i)||/||b|| 2.899250550202e-05
    7 KSP unpreconditioned resid norm 4.930569892840e-01 true resid norm 4.925940797076e-01 ||r(i)||/||b|| 5.956300898619e-06
    8 KSP unpreconditioned resid norm 9.718220057380e-02 true resid norm 9.688329523332e-02 ||r(i)||/||b|| 1.171483950440e-06
    9 KSP unpreconditioned resid norm 9.323556486527e-03 true resid norm 9.363430919214e-03 ||r(i)||/||b|| 1.132198178902e-07
      Line search: gnorm after quadratic fit 7.443321069910e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 1 Nonlinear |R| = 7.443321e+04
    0 KSP unpreconditioned resid norm 7.443321069910e+04 true resid norm 7.443321069910e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.020684857424e+04 true resid norm 1.020684859751e+04 ||r(i)||/||b|| 1.371276141610e-01
    2 KSP unpreconditioned resid norm 1.153469966381e+03 true resid norm 1.153471975238e+03 ||r(i)||/||b|| 1.549673814154e-02
    3 KSP unpreconditioned resid norm 1.424544975707e+02 true resid norm 1.424535964887e+02 ||r(i)||/||b|| 1.913844574898e-03
    4 KSP unpreconditioned resid norm 2.987452517287e+01 true resid norm 2.987364091522e+01 ||r(i)||/||b|| 4.013482776658e-04
    5 KSP unpreconditioned resid norm 4.621415079804e+00 true resid norm 4.620967698847e+00 ||r(i)||/||b|| 6.208206868205e-05
    6 KSP unpreconditioned resid norm 5.891336009124e-01 true resid norm 5.902055614756e-01 ||r(i)||/||b|| 7.929330952302e-06
    7 KSP unpreconditioned resid norm 8.467547073249e-02 true resid norm 8.498417775228e-02 ||r(i)||/||b|| 1.141750798522e-06
    8 KSP unpreconditioned resid norm 1.137295231327e-02 true resid norm 1.346957267118e-02 ||r(i)||/||b|| 1.809618656064e-07
      Line search: gnorm after quadratic fit 6.361044952309e+04
      Line search: Quadratically determined step, lambda=1.8955718443621131e-01
 2 Nonlinear |R| = 6.361045e+04
    0 KSP unpreconditioned resid norm 6.361044952309e+04 true resid norm 6.361044952309e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.122986438101e+03 true resid norm 8.122988132306e+03 ||r(i)||/||b|| 1.276989581619e-01
    2 KSP unpreconditioned resid norm 9.098872508608e+02 true resid norm 9.098858582747e+02 ||r(i)||/||b|| 1.430403125739e-02
    3 KSP unpreconditioned resid norm 1.984543802226e+02 true resid norm 1.984527341726e+02 ||r(i)||/||b|| 3.119813421545e-03
    4 KSP unpreconditioned resid norm 2.477294275161e+01 true resid norm 2.477399743344e+01 ||r(i)||/||b|| 3.894642722882e-04
    5 KSP unpreconditioned resid norm 3.106716189530e+00 true resid norm 3.106830768973e+00 ||r(i)||/||b|| 4.884151569854e-05
    6 KSP unpreconditioned resid norm 4.129566204632e-01 true resid norm 4.164223383656e-01 ||r(i)||/||b|| 6.546445458060e-06
    7 KSP unpreconditioned resid norm 4.427983464642e-02 true resid norm 5.154792683407e-02 ||r(i)||/||b|| 8.103688500953e-07
      Line search: Using full step: fnorm 6.361044952309e+04 gnorm 1.761950911464e+04
 3 Nonlinear |R| = 1.761951e+04
    0 KSP unpreconditioned resid norm 1.761950911464e+04 true resid norm 1.761950911464e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.172292739602e+02 true resid norm 4.172292785439e+02 ||r(i)||/||b|| 2.367996042507e-02
    2 KSP unpreconditioned resid norm 8.000090014515e+01 true resid norm 8.000079752954e+01 ||r(i)||/||b|| 4.540466877313e-03
    3 KSP unpreconditioned resid norm 7.203363322524e+00 true resid norm 7.203549597096e+00 ||r(i)||/||b|| 4.088394035399e-04
    4 KSP unpreconditioned resid norm 1.895192313686e-01 true resid norm 1.894985277703e-01 ||r(i)||/||b|| 1.075504013973e-05
    5 KSP unpreconditioned resid norm 5.311541633733e-03 true resid norm 5.381344986609e-03 ||r(i)||/||b|| 3.054196885733e-07
      Line search: gnorm after quadratic fit 1.587280357347e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 4 Nonlinear |R| = 1.587280e+04
    0 KSP unpreconditioned resid norm 1.587280357347e+04 true resid norm 1.587280357347e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.722475393147e+02 true resid norm 3.722475682512e+02 ||r(i)||/||b|| 2.345191046611e-02
    2 KSP unpreconditioned resid norm 7.215227368985e+01 true resid norm 7.215228392967e+01 ||r(i)||/||b|| 4.545654685116e-03
    3 KSP unpreconditioned resid norm 4.286495074397e+00 true resid norm 4.286702737577e+00 ||r(i)||/||b|| 2.700658845638e-04
    4 KSP unpreconditioned resid norm 1.897616061218e-01 true resid norm 1.900786718139e-01 ||r(i)||/||b|| 1.197511648992e-05
    5 KSP unpreconditioned resid norm 5.959666593236e-03 true resid norm 6.055960564889e-03 ||r(i)||/||b|| 3.815306185110e-07
      Line search: gnorm after quadratic fit 1.428332331259e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 5 Nonlinear |R| = 1.428332e+04
    0 KSP unpreconditioned resid norm 1.428332331259e+04 true resid norm 1.428332331259e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.100173874067e+02 true resid norm 3.100173831546e+02 ||r(i)||/||b|| 2.170484952065e-02
    2 KSP unpreconditioned resid norm 5.410823430117e+01 true resid norm 5.410813997280e+01 ||r(i)||/||b|| 3.788203822642e-03
    3 KSP unpreconditioned resid norm 2.930102926550e+00 true resid norm 2.930106558214e+00 ||r(i)||/||b|| 2.051417932709e-04
    4 KSP unpreconditioned resid norm 1.239892113986e-01 true resid norm 1.239790677697e-01 ||r(i)||/||b|| 8.679987497057e-06
    5 KSP unpreconditioned resid norm 3.559874631023e-03 true resid norm 3.751749982187e-03 ||r(i)||/||b|| 2.626664607444e-07
      Line search: gnorm after quadratic fit 1.285371676428e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 6 Nonlinear |R| = 1.285372e+04
    0 KSP unpreconditioned resid norm 1.285371676428e+04 true resid norm 1.285371676428e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.588420631619e+02 true resid norm 2.588420122273e+02 ||r(i)||/||b|| 2.013752263055e-02
    2 KSP unpreconditioned resid norm 4.080096089572e+01 true resid norm 4.080109132256e+01 ||r(i)||/||b|| 3.174264072471e-03
    3 KSP unpreconditioned resid norm 2.016213652444e+00 true resid norm 2.016033108543e+00 ||r(i)||/||b|| 1.568443700383e-04
    4 KSP unpreconditioned resid norm 8.187186575581e-02 true resid norm 8.147396534172e-02 ||r(i)||/||b|| 6.338553029903e-06
    5 KSP unpreconditioned resid norm 2.148540305776e-03 true resid norm 2.499813699323e-03 ||r(i)||/||b|| 1.944817787078e-07
      Line search: gnorm after quadratic fit 1.156767295044e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 7 Nonlinear |R| = 1.156767e+04
    0 KSP unpreconditioned resid norm 1.156767295044e+04 true resid norm 1.156767295044e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.166534922156e+02 true resid norm 2.166535082966e+02 ||r(i)||/||b|| 1.872922144538e-02
    2 KSP unpreconditioned resid norm 3.093883746085e+01 true resid norm 3.093899103672e+01 ||r(i)||/||b|| 2.674608036488e-03
    3 KSP unpreconditioned resid norm 1.396710891305e+00 true resid norm 1.396760369849e+00 ||r(i)||/||b|| 1.207468758698e-04
    4 KSP unpreconditioned resid norm 5.463183111967e-02 true resid norm 5.462844139338e-02 ||r(i)||/||b|| 4.722509153518e-06
    5 KSP unpreconditioned resid norm 1.310818184304e-03 true resid norm 1.531115016352e-03 ||r(i)||/||b|| 1.323615409003e-07
      Line search: gnorm after quadratic fit 1.219964359908e+04
      Line search: Cubically determined step, current gnorm 1.119769081047e+04 lambda=3.3228231815844042e-02
 8 Nonlinear |R| = 1.119769e+04
    0 KSP unpreconditioned resid norm 1.119769081047e+04 true resid norm 1.119769081047e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.070371209438e+02 true resid norm 6.070367348966e+02 ||r(i)||/||b|| 5.421088554515e-02
    2 KSP unpreconditioned resid norm 1.595940617164e+01 true resid norm 1.595944994392e+01 ||r(i)||/||b|| 1.425244741442e-03
    3 KSP unpreconditioned resid norm 3.608946642401e-01 true resid norm 3.609372610789e-01 ||r(i)||/||b|| 3.223318693005e-05
    4 KSP unpreconditioned resid norm 2.555175812755e-02 true resid norm 2.582638990578e-02 ||r(i)||/||b|| 2.306403198918e-06
    5 KSP unpreconditioned resid norm 6.492168672807e-04 true resid norm 7.610264746032e-04 ||r(i)||/||b|| 6.796280478572e-08
      Line search: Using full step: fnorm 1.119769081047e+04 gnorm 3.617310770309e+02
 9 Nonlinear |R| = 3.617311e+02
    0 KSP unpreconditioned resid norm 3.617310770309e+02 true resid norm 3.617310770309e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.417447668611e+00 true resid norm 9.417432385992e+00 ||r(i)||/||b|| 2.603434701627e-02
    2 KSP unpreconditioned resid norm 4.073494631797e-02 true resid norm 4.072932462991e-02 ||r(i)||/||b|| 1.125955916318e-04
    3 KSP unpreconditioned resid norm 3.040591958843e-04 true resid norm 3.037572033925e-04 ||r(i)||/||b|| 8.397321178089e-07
      Line search: Using full step: fnorm 3.617310770309e+02 gnorm 3.350912942520e+00
10 Nonlinear |R| = 3.350913e+00
    0 KSP unpreconditioned resid norm 3.350912942520e+00 true resid norm 3.350912942520e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.439261380404e-02 true resid norm 3.439262392208e-02 ||r(i)||/||b|| 1.026365784848e-02
    2 KSP unpreconditioned resid norm 1.143545961767e-04 true resid norm 1.143019547112e-04 ||r(i)||/||b|| 3.411069063024e-05
    3 KSP unpreconditioned resid norm 1.746431377661e-06 true resid norm 1.728303449799e-06 ||r(i)||/||b|| 5.157709195809e-07
      Line search: Using full step: fnorm 3.350912942520e+00 gnorm 8.349711371726e-05
11 Nonlinear |R| = 8.349711e-05
    0 KSP unpreconditioned resid norm 8.349711371726e-05 true resid norm 8.349711371726e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.549893380050e-06 true resid norm 1.549894871763e-06 ||r(i)||/||b|| 1.856225685850e-02
    2 KSP unpreconditioned resid norm 6.503910364894e-09 true resid norm 6.502295880661e-09 ||r(i)||/||b|| 7.787449878422e-05
    3 KSP unpreconditioned resid norm 6.190126353518e-11 true resid norm 6.205964633081e-11 ||r(i)||/||b|| 7.432549889204e-07
      Line search: Using full step: fnorm 8.349711371726e-05 gnorm 8.899443319815e-10
12 Nonlinear |R| = 8.899443e-10
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
+----------------+------------------+----------------+----------------+


Time Step  4, time = 0.4
                dt = 0.1
 0 Nonlinear |R| = 7.921692e+04
    0 KSP unpreconditioned resid norm 7.921691850282e+04 true resid norm 7.921691850282e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.338784384816e+04 true resid norm 1.338784232409e+04 ||r(i)||/||b|| 1.690023113385e-01
    2 KSP unpreconditioned resid norm 1.487613343738e+03 true resid norm 1.487613967943e+03 ||r(i)||/||b|| 1.877899312494e-02
    3 KSP unpreconditioned resid norm 1.848728565726e+02 true resid norm 1.848713977371e+02 ||r(i)||/||b|| 2.333736293094e-03
    4 KSP unpreconditioned resid norm 4.400465584080e+01 true resid norm 4.400552781420e+01 ||r(i)||/||b|| 5.555066852624e-04
    5 KSP unpreconditioned resid norm 7.246166392133e+00 true resid norm 7.247899817413e+00 ||r(i)||/||b|| 9.149434179462e-05
    6 KSP unpreconditioned resid norm 1.097335800364e+00 true resid norm 1.096682593980e+00 ||r(i)||/||b|| 1.384404511949e-05
    7 KSP unpreconditioned resid norm 1.808034946183e-01 true resid norm 1.817969352404e-01 ||r(i)||/||b|| 2.294925612815e-06
    8 KSP unpreconditioned resid norm 2.631367873419e-02 true resid norm 3.380790339085e-02 ||r(i)||/||b|| 4.267763001870e-07
      Line search: Using full step: fnorm 7.921691850282e+04 gnorm 1.702544160360e+04
 1 Nonlinear |R| = 1.702544e+04
    0 KSP unpreconditioned resid norm 1.702544160360e+04 true resid norm 1.702544160360e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.936879290108e+02 true resid norm 4.936879034764e+02 ||r(i)||/||b|| 2.899706891433e-02
    2 KSP unpreconditioned resid norm 3.615326133695e+01 true resid norm 3.615328661373e+01 ||r(i)||/||b|| 2.123485983828e-03
    3 KSP unpreconditioned resid norm 1.334792009103e+00 true resid norm 1.334771957736e+00 ||r(i)||/||b|| 7.839866881655e-05
    4 KSP unpreconditioned resid norm 6.918964447437e-02 true resid norm 6.884348712891e-02 ||r(i)||/||b|| 4.043565431769e-06
    5 KSP unpreconditioned resid norm 1.618852708753e-03 true resid norm 2.068841111930e-03 ||r(i)||/||b|| 1.215146813868e-07
      Line search: gnorm after quadratic fit 1.497347543826e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 2 Nonlinear |R| = 1.497348e+04
    0 KSP unpreconditioned resid norm 1.497347543826e+04 true resid norm 1.497347543826e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.778217839284e+02 true resid norm 3.778222051681e+02 ||r(i)||/||b|| 2.523276621556e-02
    2 KSP unpreconditioned resid norm 2.294024601391e+01 true resid norm 2.294053670185e+01 ||r(i)||/||b|| 1.532078293809e-03
    3 KSP unpreconditioned resid norm 4.320502546624e-01 true resid norm 4.318240513541e-01 ||r(i)||/||b|| 2.883926668425e-05
    4 KSP unpreconditioned resid norm 2.375254205955e-02 true resid norm 2.393710946054e-02 ||r(i)||/||b|| 1.598634168750e-06
    5 KSP unpreconditioned resid norm 3.647971061914e-04 true resid norm 2.097317106711e-03 ||r(i)||/||b|| 1.400688247267e-07
      Line search: Using full step: fnorm 1.497347543826e+04 gnorm 4.484975606987e+02
 3 Nonlinear |R| = 4.484976e+02
    0 KSP unpreconditioned resid norm 4.484975606987e+02 true resid norm 4.484975606987e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.197702003556e+00 true resid norm 7.197702294893e+00 ||r(i)||/||b|| 1.604847590181e-02
    2 KSP unpreconditioned resid norm 4.439509032672e-02 true resid norm 4.439185073553e-02 ||r(i)||/||b|| 9.897902380197e-05
    3 KSP unpreconditioned resid norm 8.467085974353e-04 true resid norm 8.474710283934e-04 ||r(i)||/||b|| 1.889577787387e-06
    4 KSP unpreconditioned resid norm 6.567418701465e-06 true resid norm 1.244185148847e-05 ||r(i)||/||b|| 2.774117983849e-08
      Line search: Using full step: fnorm 4.484975606987e+02 gnorm 4.992912287557e-01
 4 Nonlinear |R| = 4.992912e-01
    0 KSP unpreconditioned resid norm 4.992912287557e-01 true resid norm 4.992912287557e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.039916560602e-02 true resid norm 1.039918033896e-02 ||r(i)||/||b|| 2.082788508998e-02
    2 KSP unpreconditioned resid norm 5.148726564771e-05 true resid norm 5.146978863360e-05 ||r(i)||/||b|| 1.030857056349e-04
    3 KSP unpreconditioned resid norm 2.809515374181e-07 true resid norm 2.776218990212e-07 ||r(i)||/||b|| 5.560319970231e-07
      Line search: Using full step: fnorm 4.992912287557e-01 gnorm 3.718894787665e-06
 5 Nonlinear |R| = 3.718895e-06
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
+----------------+------------------+----------------+----------------+


Time Step  5, time = 0.5
                dt = 0.1
 0 Nonlinear |R| = 7.837058e+04
    0 KSP unpreconditioned resid norm 7.837057971834e+04 true resid norm 7.837057971834e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.419207447533e+04 true resid norm 1.419207056344e+04 ||r(i)||/||b|| 1.810892635277e-01
    2 KSP unpreconditioned resid norm 1.309093243923e+03 true resid norm 1.309093697773e+03 ||r(i)||/||b|| 1.670389197678e-02
    3 KSP unpreconditioned resid norm 1.541929577016e+02 true resid norm 1.541919768807e+02 ||r(i)||/||b|| 1.967472710229e-03
    4 KSP unpreconditioned resid norm 3.265839993452e+01 true resid norm 3.265867801907e+01 ||r(i)||/||b|| 4.167211488858e-04
    5 KSP unpreconditioned resid norm 5.186215038185e+00 true resid norm 5.186750659629e+00 ||r(i)||/||b|| 6.618236943341e-05
    6 KSP unpreconditioned resid norm 7.964714602642e-01 true resid norm 7.977876759260e-01 ||r(i)||/||b|| 1.017968322798e-05
    7 KSP unpreconditioned resid norm 1.114365411466e-01 true resid norm 1.129494189914e-01 ||r(i)||/||b|| 1.441222195846e-06
    8 KSP unpreconditioned resid norm 1.569352686975e-02 true resid norm 1.948198623943e-02 ||r(i)||/||b|| 2.485880072529e-07
      Line search: Using full step: fnorm 7.837057971834e+04 gnorm 1.516012930827e+04
 1 Nonlinear |R| = 1.516013e+04
    0 KSP unpreconditioned resid norm 1.516012930827e+04 true resid norm 1.516012930827e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.219890784187e+02 true resid norm 5.219890923488e+02 ||r(i)||/||b|| 3.443170448844e-02
    2 KSP unpreconditioned resid norm 2.964148768127e+01 true resid norm 2.964154902785e+01 ||r(i)||/||b|| 1.955230620077e-03
    3 KSP unpreconditioned resid norm 9.150087977007e-01 true resid norm 9.148435718274e-01 ||r(i)||/||b|| 6.034536732668e-05
    4 KSP unpreconditioned resid norm 6.699842070041e-02 true resid norm 6.698891719346e-02 ||r(i)||/||b|| 4.418756320035e-06
    5 KSP unpreconditioned resid norm 1.686526837184e-03 true resid norm 1.569823843207e-03 ||r(i)||/||b|| 1.035495022032e-07
      Line search: gnorm after quadratic fit 1.520635961110e+04
      Line search: Cubically determined step, current gnorm 1.443954868546e+04 lambda=5.0000000000000003e-02
 2 Nonlinear |R| = 1.443955e+04
    0 KSP unpreconditioned resid norm 1.443954868546e+04 true resid norm 1.443954868546e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.113546389278e+02 true resid norm 4.113545155630e+02 ||r(i)||/||b|| 2.848804519613e-02
    2 KSP unpreconditioned resid norm 1.932809277486e+01 true resid norm 1.932833729396e+01 ||r(i)||/||b|| 1.338569349708e-03
    3 KSP unpreconditioned resid norm 5.228271279245e-01 true resid norm 5.227267054135e-01 ||r(i)||/||b|| 3.620104179156e-05
    4 KSP unpreconditioned resid norm 3.083797126353e-02 true resid norm 3.069392859554e-02 ||r(i)||/||b|| 2.125684760941e-06
    5 KSP unpreconditioned resid norm 7.272371338192e-04 true resid norm 1.033833152997e-03 ||r(i)||/||b|| 7.159733143449e-08
      Line search: gnorm after quadratic fit 1.296717879791e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 3 Nonlinear |R| = 1.296718e+04
    0 KSP unpreconditioned resid norm 1.296717879791e+04 true resid norm 1.296717879791e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.604015001550e+02 true resid norm 3.604014706066e+02 ||r(i)||/||b|| 2.779336016133e-02
    2 KSP unpreconditioned resid norm 1.074795690930e+01 true resid norm 1.074777244866e+01 ||r(i)||/||b|| 8.288443165753e-04
    3 KSP unpreconditioned resid norm 4.248035478774e-01 true resid norm 4.247519397758e-01 ||r(i)||/||b|| 3.275592527839e-05
    4 KSP unpreconditioned resid norm 1.204108501922e-02 true resid norm 1.226240934374e-02 ||r(i)||/||b|| 9.456497465527e-07
      Line search: Using full step: fnorm 1.296717879791e+04 gnorm 3.661496628173e+02
 4 Nonlinear |R| = 3.661497e+02
    0 KSP unpreconditioned resid norm 3.661496628173e+02 true resid norm 3.661496628173e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.756094103152e+00 true resid norm 2.756101500967e+00 ||r(i)||/||b|| 7.527253964296e-03
    2 KSP unpreconditioned resid norm 1.616858102776e-02 true resid norm 1.616897129359e-02 ||r(i)||/||b|| 4.415945973891e-05
    3 KSP unpreconditioned resid norm 1.445686883717e-04 true resid norm 1.461622802669e-04 ||r(i)||/||b|| 3.991872589537e-07
      Line search: Using full step: fnorm 3.661496628173e+02 gnorm 4.687047704526e-01
 5 Nonlinear |R| = 4.687048e-01
    0 KSP unpreconditioned resid norm 4.687047704526e-01 true resid norm 4.687047704526e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.409678264645e-03 true resid norm 2.409679386476e-03 ||r(i)||/||b|| 5.141145425400e-03
    2 KSP unpreconditioned resid norm 1.261756740388e-05 true resid norm 1.261612180677e-05 ||r(i)||/||b|| 2.691699040015e-05
    3 KSP unpreconditioned resid norm 8.598983796368e-08 true resid norm 8.707901355197e-08 ||r(i)||/||b|| 1.857864887270e-07
      Line search: Using full step: fnorm 4.687047704526e-01 gnorm 9.913219296694e-07
 6 Nonlinear |R| = 9.913219e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
+----------------+------------------+----------------+----------------+


Time Step  6, time = 0.6
                dt = 0.1
 0 Nonlinear |R| = 8.024976e+04
    0 KSP unpreconditioned resid norm 8.024976073441e+04 true resid norm 8.024976073441e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.369790790203e+04 true resid norm 1.369790551860e+04 ||r(i)||/||b|| 1.706909203622e-01
    2 KSP unpreconditioned resid norm 1.276427612946e+03 true resid norm 1.276427712415e+03 ||r(i)||/||b|| 1.590568869905e-02
    3 KSP unpreconditioned resid norm 1.523081098237e+02 true resid norm 1.523062699260e+02 ||r(i)||/||b|| 1.897903103164e-03
    4 KSP unpreconditioned resid norm 2.491934933829e+01 true resid norm 2.492145032911e+01 ||r(i)||/||b|| 3.105485935539e-04
    5 KSP unpreconditioned resid norm 2.410246196232e+00 true resid norm 2.410773665172e+00 ||r(i)||/||b|| 3.004088290245e-05
    6 KSP unpreconditioned resid norm 3.184930862734e-01 true resid norm 3.221416470892e-01 ||r(i)||/||b|| 4.014238100414e-06
    7 KSP unpreconditioned resid norm 4.227294412935e-02 true resid norm 5.066514832123e-02 ||r(i)||/||b|| 6.313432944543e-07
      Line search: Using full step: fnorm 8.024976073441e+04 gnorm 1.741941982599e+04
 1 Nonlinear |R| = 1.741942e+04
    0 KSP unpreconditioned resid norm 1.741941982599e+04 true resid norm 1.741941982599e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.799398078482e+02 true resid norm 5.799398556119e+02 ||r(i)||/||b|| 3.329271935606e-02
    2 KSP unpreconditioned resid norm 1.605572413116e+01 true resid norm 1.605577293075e+01 ||r(i)||/||b|| 9.217168591802e-04
    3 KSP unpreconditioned resid norm 1.486884259577e+00 true resid norm 1.487672234817e+00 ||r(i)||/||b|| 8.540308745518e-05
    4 KSP unpreconditioned resid norm 1.916910356106e-02 true resid norm 1.892792828166e-02 ||r(i)||/||b|| 1.086599236412e-06
    5 KSP unpreconditioned resid norm 5.702953425128e-04 true resid norm 9.526074166818e-04 ||r(i)||/||b|| 5.468651804698e-08
      Line search: Using full step: fnorm 1.741941982599e+04 gnorm 6.350512006889e+02
 2 Nonlinear |R| = 6.350512e+02
    0 KSP unpreconditioned resid norm 6.350512006889e+02 true resid norm 6.350512006889e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.593839174593e+01 true resid norm 1.593839174593e+01 ||r(i)||/||b|| 2.509780585981e-02
    2 KSP unpreconditioned resid norm 4.526852164449e-02 true resid norm 4.526026564432e-02 ||r(i)||/||b|| 7.127026229573e-05
    3 KSP unpreconditioned resid norm 1.799499414584e-04 true resid norm 1.724055447936e-04 ||r(i)||/||b|| 2.714829050108e-07
      Line search: Using full step: fnorm 6.350512006889e+02 gnorm 4.566556202128e+00
 3 Nonlinear |R| = 4.566556e+00
    0 KSP unpreconditioned resid norm 4.566556202128e+00 true resid norm 4.566556202128e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.166990938406e-04 true resid norm 2.166990938407e-04 ||r(i)||/||b|| 4.745350418325e-05
    2 KSP unpreconditioned resid norm 4.385749616971e-09 true resid norm 1.292465075562e-07 ||r(i)||/||b|| 2.830283956563e-08
      Line search: Using full step: fnorm 4.566556202128e+00 gnorm 7.451783045044e-05
 4 Nonlinear |R| = 7.451783e-05
    0 KSP unpreconditioned resid norm 7.451783045044e-05 true resid norm 7.451783045044e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.419731633876e-12 true resid norm 4.419731633424e-12 ||r(i)||/||b|| 5.931106161717e-08
      Line search: Using full step: fnorm 7.451783045044e-05 gnorm 2.330765763340e-10
 5 Nonlinear |R| = 2.330766e-10
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
+----------------+------------------+----------------+----------------+


Time Step  7, time = 0.7
                dt = 0.1
 0 Nonlinear |R| = 8.319786e+04
    0 KSP unpreconditioned resid norm 8.319785679542e+04 true resid norm 8.319785679542e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.381886543141e+04 true resid norm 1.381886543141e+04 ||r(i)||/||b|| 1.660964111779e-01
    2 KSP unpreconditioned resid norm 1.768378427967e+03 true resid norm 1.768379482862e+03 ||r(i)||/||b|| 2.125510861668e-02
    3 KSP unpreconditioned resid norm 3.559792750981e+02 true resid norm 3.559804028517e+02 ||r(i)||/||b|| 4.278720829636e-03
    4 KSP unpreconditioned resid norm 7.882892555064e+01 true resid norm 7.882873402032e+01 ||r(i)||/||b|| 9.474851523417e-04
    5 KSP unpreconditioned resid norm 1.254574232860e+01 true resid norm 1.254554874906e+01 ||r(i)||/||b|| 1.507917298868e-04
    6 KSP unpreconditioned resid norm 1.623752513235e+00 true resid norm 1.625338858802e+00 ||r(i)||/||b|| 1.953582605858e-05
    7 KSP unpreconditioned resid norm 3.101080604162e-01 true resid norm 3.085647072713e-01 ||r(i)||/||b|| 3.708805961553e-06
    8 KSP unpreconditioned resid norm 2.044473464985e-02 true resid norm 1.988040583221e-02 ||r(i)||/||b|| 2.389533408426e-07
      Line search: Using full step: fnorm 8.319785679542e+04 gnorm 1.423597289994e+04
 1 Nonlinear |R| = 1.423597e+04
    0 KSP unpreconditioned resid norm 1.423597289994e+04 true resid norm 1.423597289994e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.942911436521e+02 true resid norm 9.942911436521e+02 ||r(i)||/||b|| 6.984356816643e-02
    2 KSP unpreconditioned resid norm 3.095460003654e+01 true resid norm 3.095447567309e+01 ||r(i)||/||b|| 2.174384279225e-03
    3 KSP unpreconditioned resid norm 7.023034551838e-01 true resid norm 7.022911657293e-01 ||r(i)||/||b|| 4.933215107009e-05
    4 KSP unpreconditioned resid norm 1.738085872761e-02 true resid norm 1.750276462207e-02 ||r(i)||/||b|| 1.229474426869e-06
    5 KSP unpreconditioned resid norm 4.274075912928e-04 true resid norm 6.711934529687e-04 ||r(i)||/||b|| 4.714770516118e-08
      Line search: Using full step: fnorm 1.423597289994e+04 gnorm 2.872212947596e+02
 2 Nonlinear |R| = 2.872213e+02
    0 KSP unpreconditioned resid norm 2.872212947596e+02 true resid norm 2.872212947596e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.763175338460e-01 true resid norm 7.763175338460e-01 ||r(i)||/||b|| 2.702855073806e-03
    2 KSP unpreconditioned resid norm 6.649486266809e-04 true resid norm 6.652383856072e-04 ||r(i)||/||b|| 2.316117912372e-06
    3 KSP unpreconditioned resid norm 4.202144800944e-07 true resid norm 2.032916072589e-05 ||r(i)||/||b|| 7.077873784709e-08
      Line search: Using full step: fnorm 2.872212947596e+02 gnorm 2.648451192147e-01
 3 Nonlinear |R| = 2.648451e-01
    0 KSP unpreconditioned resid norm 2.648451192147e-01 true resid norm 2.648451192147e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.258931201014e-07 true resid norm 4.258931201082e-07 ||r(i)||/||b|| 1.608083703302e-06
    2 KSP unpreconditioned resid norm 4.465532129305e-13 true resid norm 5.847872787390e-09 ||r(i)||/||b|| 2.208034946888e-08
      Line search: Using full step: fnorm 2.648451192147e-01 gnorm 4.107131414701e-07
 4 Nonlinear |R| = 4.107131e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
+----------------+------------------+----------------+----------------+


Time Step  8, time = 0.8
                dt = 0.1
 0 Nonlinear |R| = 8.526195e+04
    0 KSP unpreconditioned resid norm 8.526194631801e+04 true resid norm 8.526194631801e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.400423270308e+04 true resid norm 1.400423270308e+04 ||r(i)||/||b|| 1.642495076390e-01
    2 KSP unpreconditioned resid norm 1.958044160186e+03 true resid norm 1.958042471184e+03 ||r(i)||/||b|| 2.296502197922e-02
    3 KSP unpreconditioned resid norm 4.025648661438e+02 true resid norm 4.025644547689e+02 ||r(i)||/||b|| 4.721502055177e-03
    4 KSP unpreconditioned resid norm 8.757201740364e+01 true resid norm 8.757210940148e+01 ||r(i)||/||b|| 1.027094890314e-03
    5 KSP unpreconditioned resid norm 1.355533304240e+01 true resid norm 1.355677588813e+01 ||r(i)||/||b|| 1.590014827666e-04
    6 KSP unpreconditioned resid norm 2.046590213300e+00 true resid norm 2.045721819646e+00 ||r(i)||/||b|| 2.399337462947e-05
    7 KSP unpreconditioned resid norm 3.558558789268e-01 true resid norm 3.544488618484e-01 ||r(i)||/||b|| 4.157175353778e-06
    8 KSP unpreconditioned resid norm 1.935259330569e-02 true resid norm 1.731626069440e-02 ||r(i)||/||b|| 2.030948323630e-07
      Line search: Using full step: fnorm 8.526194631801e+04 gnorm 1.441830092281e+04
 1 Nonlinear |R| = 1.441830e+04
    0 KSP unpreconditioned resid norm 1.441830092281e+04 true resid norm 1.441830092281e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.135401891861e+03 true resid norm 1.135401891861e+03 ||r(i)||/||b|| 7.874727389440e-02
    2 KSP unpreconditioned resid norm 3.485648122033e+01 true resid norm 3.485647176502e+01 ||r(i)||/||b|| 2.417515902298e-03
    3 KSP unpreconditioned resid norm 8.347848956003e-01 true resid norm 8.348548532849e-01 ||r(i)||/||b|| 5.790244341232e-05
    4 KSP unpreconditioned resid norm 2.088779820303e-02 true resid norm 2.091340405684e-02 ||r(i)||/||b|| 1.450476319561e-06
    5 KSP unpreconditioned resid norm 5.360067853916e-04 true resid norm 7.717407320330e-04 ||r(i)||/||b|| 5.352508150333e-08
      Line search: Using full step: fnorm 1.441830092281e+04 gnorm 2.884410160902e+02
 2 Nonlinear |R| = 2.884410e+02
    0 KSP unpreconditioned resid norm 2.884410160902e+02 true resid norm 2.884410160902e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.338397558483e-01 true resid norm 8.338397558483e-01 ||r(i)||/||b|| 2.890850154222e-03
    2 KSP unpreconditioned resid norm 7.126101678654e-04 true resid norm 7.117029948407e-04 ||r(i)||/||b|| 2.467412590926e-06
    3 KSP unpreconditioned resid norm 4.708253071569e-07 true resid norm 1.748273732198e-05 ||r(i)||/||b|| 6.061113484813e-08
      Line search: Using full step: fnorm 2.884410160902e+02 gnorm 2.853286226934e-01
 3 Nonlinear |R| = 2.853286e-01
    0 KSP unpreconditioned resid norm 2.853286226934e-01 true resid norm 2.853286226934e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.005909574462e-07 true resid norm 2.005909574406e-07 ||r(i)||/||b|| 7.030172982546e-07
      Line search: Using full step: fnorm 2.853286226934e-01 gnorm 1.698202696636e-07
 4 Nonlinear |R| = 1.698203e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
+----------------+------------------+----------------+----------------+


Time Step  9, time = 0.9
                dt = 0.1
 0 Nonlinear |R| = 8.488066e+04
    0 KSP unpreconditioned resid norm 8.488065597705e+04 true resid norm 8.488065597705e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.397131785757e+04 true resid norm 1.397131785757e+04 ||r(i)||/||b|| 1.645995509430e-01
    2 KSP unpreconditioned resid norm 1.923124553807e+03 true resid norm 1.923124855665e+03 ||r(i)||/||b|| 2.265680953485e-02
    3 KSP unpreconditioned resid norm 3.943172891187e+02 true resid norm 3.943166468551e+02 ||r(i)||/||b|| 4.645541935512e-03
    4 KSP unpreconditioned resid norm 8.610764615604e+01 true resid norm 8.610629934674e+01 ||r(i)||/||b|| 1.014439607654e-03
    5 KSP unpreconditioned resid norm 1.338571232978e+01 true resid norm 1.338803924954e+01 ||r(i)||/||b|| 1.577278014104e-04
    6 KSP unpreconditioned resid norm 1.971208516180e+00 true resid norm 1.971053306009e+00 ||r(i)||/||b|| 2.322146646159e-05
    7 KSP unpreconditioned resid norm 3.499341436901e-01 true resid norm 3.513869907515e-01 ||r(i)||/||b|| 4.139777039972e-06
    8 KSP unpreconditioned resid norm 1.954928226793e-02 true resid norm 2.181683890283e-02 ||r(i)||/||b|| 2.570295746622e-07
      Line search: Using full step: fnorm 8.488065597705e+04 gnorm 1.438512299166e+04
 1 Nonlinear |R| = 1.438512e+04
    0 KSP unpreconditioned resid norm 1.438512299166e+04 true resid norm 1.438512299166e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.110461908845e+03 true resid norm 1.110461908845e+03 ||r(i)||/||b|| 7.719516263360e-02
    2 KSP unpreconditioned resid norm 3.416332126120e+01 true resid norm 3.416328039677e+01 ||r(i)||/||b|| 2.374903601213e-03
    3 KSP unpreconditioned resid norm 8.106230833651e-01 true resid norm 8.105381085142e-01 ||r(i)||/||b|| 5.634558070753e-05
    4 KSP unpreconditioned resid norm 2.024172665048e-02 true resid norm 2.021265866273e-02 ||r(i)||/||b|| 1.405108505117e-06
    5 KSP unpreconditioned resid norm 5.148102737572e-04 true resid norm 7.574054520610e-04 ||r(i)||/||b|| 5.265199696244e-08
      Line search: Using full step: fnorm 1.438512299166e+04 gnorm 2.879370944545e+02
 2 Nonlinear |R| = 2.879371e+02
    0 KSP unpreconditioned resid norm 2.879370944545e+02 true resid norm 2.879370944545e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.244391888794e-01 true resid norm 8.244391888794e-01 ||r(i)||/||b|| 2.863261471890e-03
    2 KSP unpreconditioned resid norm 7.045923377989e-04 true resid norm 7.068015976048e-04 ||r(i)||/||b|| 2.454708376300e-06
    3 KSP unpreconditioned resid norm 4.623867745072e-07 true resid norm 2.050079839397e-05 ||r(i)||/||b|| 7.119887916077e-08
      Line search: Using full step: fnorm 2.879370944545e+02 gnorm 2.818197880216e-01
 3 Nonlinear |R| = 2.818198e-01
    0 KSP unpreconditioned resid norm 2.818197880216e-01 true resid norm 2.818197880216e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.962520839643e-07 true resid norm 1.962520839675e-07 ||r(i)||/||b|| 6.963743935271e-07
      Line search: Using full step: fnorm 2.818197880216e-01 gnorm 1.676926269811e-07
 4 Nonlinear |R| = 1.676926e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
+----------------+------------------+----------------+----------------+


Time Step 10, time = 1
                dt = 0.1
 0 Nonlinear |R| = 8.229029e+04
    0 KSP unpreconditioned resid norm 8.229028677920e+04 true resid norm 8.229028677920e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.373005385049e+04 true resid norm 1.373005385049e+04 ||r(i)||/||b|| 1.668490217725e-01
    2 KSP unpreconditioned resid norm 1.684771910299e+03 true resid norm 1.684771600110e+03 ||r(i)||/||b|| 2.047351717986e-02
    3 KSP unpreconditioned resid norm 3.336882575014e+02 true resid norm 3.336880870743e+02 ||r(i)||/||b|| 4.055011838391e-03
    4 KSP unpreconditioned resid norm 7.416116827948e+01 true resid norm 7.415980420262e+01 ||r(i)||/||b|| 9.011975423248e-04
    5 KSP unpreconditioned resid norm 1.200769712893e+01 true resid norm 1.200819775993e+01 ||r(i)||/||b|| 1.459248500635e-04
    6 KSP unpreconditioned resid norm 1.429645119798e+00 true resid norm 1.429628176593e+00 ||r(i)||/||b|| 1.737298814415e-05
    7 KSP unpreconditioned resid norm 2.742326745896e-01 true resid norm 2.743576589591e-01 ||r(i)||/||b|| 3.334022394347e-06
    8 KSP unpreconditioned resid norm 2.094736568922e-02 true resid norm 2.164134759009e-02 ||r(i)||/||b|| 2.629878742330e-07
      Line search: Using full step: fnorm 8.229028677920e+04 gnorm 1.415312760813e+04
 1 Nonlinear |R| = 1.415313e+04
    0 KSP unpreconditioned resid norm 1.415312760813e+04 true resid norm 1.415312760813e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.263779009219e+02 true resid norm 9.263779009219e+02 ||r(i)||/||b|| 6.545393545309e-02
    2 KSP unpreconditioned resid norm 2.909384677650e+01 true resid norm 2.909370470783e+01 ||r(i)||/||b|| 2.055637842982e-03
    3 KSP unpreconditioned resid norm 6.428042044885e-01 true resid norm 6.429168274963e-01 ||r(i)||/||b|| 4.542577763003e-05
    4 KSP unpreconditioned resid norm 1.582234804418e-02 true resid norm 1.584430587518e-02 ||r(i)||/||b|| 1.119491487244e-06
    5 KSP unpreconditioned resid norm 3.838640949817e-04 true resid norm 6.543254433572e-04 ||r(i)||/||b|| 4.623186206428e-08
      Line search: Using full step: fnorm 1.415312760813e+04 gnorm 2.878742924335e+02
 2 Nonlinear |R| = 2.878743e+02
    0 KSP unpreconditioned resid norm 2.878742924335e+02 true resid norm 2.878742924335e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.449059049874e-01 true resid norm 7.449059049874e-01 ||r(i)||/||b|| 2.587608287946e-03
    2 KSP unpreconditioned resid norm 6.397222539264e-04 true resid norm 6.390926120107e-04 ||r(i)||/||b|| 2.220040583021e-06
    3 KSP unpreconditioned resid norm 3.935397164482e-07 true resid norm 2.144365117343e-05 ||r(i)||/||b|| 7.448963570927e-08
      Line search: Using full step: fnorm 2.878742924335e+02 gnorm 2.544007949734e-01
 3 Nonlinear |R| = 2.544008e-01
    0 KSP unpreconditioned resid norm 2.544007949734e-01 true resid norm 2.544007949734e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.451840402603e-07 true resid norm 3.451840402456e-07 ||r(i)||/||b|| 1.356851264092e-06
    2 KSP unpreconditioned resid norm 3.373125628620e-13 true resid norm 5.602058318320e-09 ||r(i)||/||b|| 2.202060067818e-08
      Line search: Using full step: fnorm 2.544007949734e-01 gnorm 3.282040281623e-07
 4 Nonlinear |R| = 3.282040e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
+----------------+------------------+----------------+----------------+


Time Step 11, time = 1.1
                dt = 0.1
 0 Nonlinear |R| = 7.919076e+04
    0 KSP unpreconditioned resid norm 7.919076178227e+04 true resid norm 7.919076178227e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.334780499871e+04 true resid norm 1.334780499871e+04 ||r(i)||/||b|| 1.685525520692e-01
    2 KSP unpreconditioned resid norm 1.407843137136e+03 true resid norm 1.407844159614e+03 ||r(i)||/||b|| 1.777788378252e-02
    3 KSP unpreconditioned resid norm 2.428109233038e+02 true resid norm 2.428110258504e+02 ||r(i)||/||b|| 3.066153429840e-03
    4 KSP unpreconditioned resid norm 5.030284392439e+01 true resid norm 5.030166348777e+01 ||r(i)||/||b|| 6.351961056528e-04
    5 KSP unpreconditioned resid norm 8.963187432532e+00 true resid norm 8.963661968328e+00 ||r(i)||/||b|| 1.131907531458e-04
    6 KSP unpreconditioned resid norm 9.719240983516e-01 true resid norm 9.736496911421e-01 ||r(i)||/||b|| 1.229499084526e-05
    7 KSP unpreconditioned resid norm 7.777219689827e-02 true resid norm 7.774204108017e-02 ||r(i)||/||b|| 9.817059380476e-07
      Line search: Using full step: fnorm 7.919076178227e+04 gnorm 1.384347019725e+04
 1 Nonlinear |R| = 1.384347e+04
    0 KSP unpreconditioned resid norm 1.384347019725e+04 true resid norm 1.384347019725e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.392799323067e+02 true resid norm 6.392799323067e+02 ||r(i)||/||b|| 4.617916773741e-02
    2 KSP unpreconditioned resid norm 2.133729716712e+01 true resid norm 2.133726081422e+01 ||r(i)||/||b|| 1.541323129981e-03
    3 KSP unpreconditioned resid norm 4.360451956095e-01 true resid norm 4.360047677188e-01 ||r(i)||/||b|| 3.149533762174e-05
    4 KSP unpreconditioned resid norm 9.961390009428e-03 true resid norm 9.919767331832e-03 ||r(i)||/||b|| 7.165665249023e-07
      Line search: Using full step: fnorm 1.384347019725e+04 gnorm 2.959414976869e+02
 2 Nonlinear |R| = 2.959415e+02
    0 KSP unpreconditioned resid norm 2.959414976869e+02 true resid norm 2.959414976869e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.848601566778e-01 true resid norm 5.848601566778e-01 ||r(i)||/||b|| 1.976269503430e-03
    2 KSP unpreconditioned resid norm 5.077495601524e-04 true resid norm 5.125741318268e-04 ||r(i)||/||b|| 1.732011684178e-06
    3 KSP unpreconditioned resid norm 2.676286565737e-07 true resid norm 1.411357378983e-05 ||r(i)||/||b|| 4.769041820813e-08
      Line search: Using full step: fnorm 2.959414976869e+02 gnorm 2.046028642534e-01
 3 Nonlinear |R| = 2.046029e-01
    0 KSP unpreconditioned resid norm 2.046028642534e-01 true resid norm 2.046028642534e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.618685034899e-08 true resid norm 9.618685034326e-08 ||r(i)||/||b|| 4.701148769069e-07
      Line search: Using full step: fnorm 2.046028642534e-01 gnorm 1.074609500119e-07
 4 Nonlinear |R| = 1.074610e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
+----------------+------------------+----------------+----------------+


Time Step 12, time = 1.2
                dt = 0.1
 0 Nonlinear |R| = 7.791159e+04
    0 KSP unpreconditioned resid norm 7.791158715921e+04 true resid norm 7.791158715921e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.292730069248e+04 true resid norm 1.292730069248e+04 ||r(i)||/||b|| 1.659226973012e-01
    2 KSP unpreconditioned resid norm 1.339946323491e+03 true resid norm 1.339942416076e+03 ||r(i)||/||b|| 1.719824309750e-02
    3 KSP unpreconditioned resid norm 1.789938122594e+02 true resid norm 1.789932061008e+02 ||r(i)||/||b|| 2.297388779091e-03
    4 KSP unpreconditioned resid norm 3.562294613730e+01 true resid norm 3.562505582189e+01 ||r(i)||/||b|| 4.572497765844e-04
    5 KSP unpreconditioned resid norm 5.937519605700e+00 true resid norm 5.937535713607e+00 ||r(i)||/||b|| 7.620863507085e-05
    6 KSP unpreconditioned resid norm 1.101018598748e+00 true resid norm 1.103475666606e+00 ||r(i)||/||b|| 1.416317786405e-05
    7 KSP unpreconditioned resid norm 2.845603681932e-01 true resid norm 2.852335323484e-01 ||r(i)||/||b|| 3.660989882873e-06
    8 KSP unpreconditioned resid norm 5.012861420151e-02 true resid norm 5.180266649689e-02 ||r(i)||/||b|| 6.648904018735e-07
      Line search: Using full step: fnorm 7.791158715921e+04 gnorm 1.363525138428e+04
 1 Nonlinear |R| = 1.363525e+04
    0 KSP unpreconditioned resid norm 1.363525138428e+04 true resid norm 1.363525138428e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.847088104453e+02 true resid norm 3.847088104453e+02 ||r(i)||/||b|| 2.821428073477e-02
    2 KSP unpreconditioned resid norm 1.394569018400e+01 true resid norm 1.394567496088e+01 ||r(i)||/||b|| 1.022766252550e-03
    3 KSP unpreconditioned resid norm 3.224196900495e-01 true resid norm 3.224440850288e-01 ||r(i)||/||b|| 2.364782840751e-05
    4 KSP unpreconditioned resid norm 5.815196062044e-03 true resid norm 5.923866447589e-03 ||r(i)||/||b|| 4.344523089921e-07
      Line search: Using full step: fnorm 1.363525138428e+04 gnorm 3.035442044709e+02
 2 Nonlinear |R| = 3.035442e+02
    0 KSP unpreconditioned resid norm 3.035442044709e+02 true resid norm 3.035442044709e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.779148982579e-01 true resid norm 3.779148982579e-01 ||r(i)||/||b|| 1.245007786976e-03
    2 KSP unpreconditioned resid norm 3.046593781703e-04 true resid norm 2.996055181933e-04 ||r(i)||/||b|| 9.870243403775e-07
    3 KSP unpreconditioned resid norm 1.292305143693e-07 true resid norm 1.086648145085e-05 ||r(i)||/||b|| 3.579867871235e-08
      Line search: Using full step: fnorm 3.035442044709e+02 gnorm 1.429218407383e-01
 3 Nonlinear |R| = 1.429218e-01
    0 KSP unpreconditioned resid norm 1.429218407383e-01 true resid norm 1.429218407383e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.001737965862e-08 true resid norm 4.001737966689e-08 ||r(i)||/||b|| 2.799948521525e-07
      Line search: Using full step: fnorm 1.429218407383e-01 gnorm 4.740782674334e-08
 4 Nonlinear |R| = 4.740783e-08
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
+----------------+------------------+----------------+----------------+


Time Step 13, time = 1.3
                dt = 0.1
 0 Nonlinear |R| = 7.980491e+04
    0 KSP unpreconditioned resid norm 7.980490856601e+04 true resid norm 7.980490856601e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.260368915929e+04 true resid norm 1.260368915929e+04 ||r(i)||/||b|| 1.579312524224e-01
    2 KSP unpreconditioned resid norm 1.545092605684e+03 true resid norm 1.545092755067e+03 ||r(i)||/||b|| 1.936087369599e-02
    3 KSP unpreconditioned resid norm 2.092121244642e+02 true resid norm 2.092134159212e+02 ||r(i)||/||b|| 2.621560749589e-03
    4 KSP unpreconditioned resid norm 5.452936038614e+01 true resid norm 5.452674078260e+01 ||r(i)||/||b|| 6.832504636917e-04
    5 KSP unpreconditioned resid norm 1.188372131856e+01 true resid norm 1.188236104979e+01 ||r(i)||/||b|| 1.488926090300e-04
    6 KSP unpreconditioned resid norm 2.111556700798e+00 true resid norm 2.111508461840e+00 ||r(i)||/||b|| 2.645837831007e-05
    7 KSP unpreconditioned resid norm 4.805191694398e-01 true resid norm 4.807129302345e-01 ||r(i)||/||b|| 6.023601039990e-06
    8 KSP unpreconditioned resid norm 9.302281644206e-02 true resid norm 9.638352788308e-02 ||r(i)||/||b|| 1.207739343544e-06
    9 KSP unpreconditioned resid norm 8.126141323107e-03 true resid norm 9.343076879795e-03 ||r(i)||/||b|| 1.170739625880e-07
      Line search: Using full step: fnorm 7.980490856601e+04 gnorm 1.367615532933e+04
 1 Nonlinear |R| = 1.367616e+04
    0 KSP unpreconditioned resid norm 1.367615532933e+04 true resid norm 1.367615532933e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.218319963446e+02 true resid norm 4.218319963446e+02 ||r(i)||/||b|| 3.084434083897e-02
    2 KSP unpreconditioned resid norm 1.175490137507e+01 true resid norm 1.175489808763e+01 ||r(i)||/||b|| 8.595177375919e-04
    3 KSP unpreconditioned resid norm 2.052990013499e-01 true resid norm 2.051468216760e-01 ||r(i)||/||b|| 1.500032843558e-05
    4 KSP unpreconditioned resid norm 7.192705669861e-03 true resid norm 7.285501498446e-03 ||r(i)||/||b|| 5.327156150983e-07
      Line search: Using full step: fnorm 1.367615532933e+04 gnorm 2.958248476628e+02
 2 Nonlinear |R| = 2.958248e+02
    0 KSP unpreconditioned resid norm 2.958248476628e+02 true resid norm 2.958248476628e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.880807665769e-01 true resid norm 2.880807665769e-01 ||r(i)||/||b|| 9.738220736119e-04
    2 KSP unpreconditioned resid norm 2.027964098933e-04 true resid norm 2.038568247216e-04 ||r(i)||/||b|| 6.891132585115e-07
      Line search: Using full step: fnorm 2.958248476628e+02 gnorm 1.096241446680e-01
 3 Nonlinear |R| = 1.096241e-01
    0 KSP unpreconditioned resid norm 1.096241446680e-01 true resid norm 1.096241446680e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.964707654132e-08 true resid norm 2.964707653970e-08 ||r(i)||/||b|| 2.704429451148e-07
      Line search: Using full step: fnorm 1.096241446680e-01 gnorm 3.089841861309e-08
 4 Nonlinear |R| = 3.089842e-08
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
+----------------+------------------+----------------+----------------+


Time Step 14, time = 1.4
                dt = 0.1
 0 Nonlinear |R| = 8.406723e+04
    0 KSP unpreconditioned resid norm 8.406722960703e+04 true resid norm 8.406722960703e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.247443483672e+04 true resid norm 1.247443483672e+04 ||r(i)||/||b|| 1.483864151944e-01
    2 KSP unpreconditioned resid norm 1.848752741253e+03 true resid norm 1.848754728313e+03 ||r(i)||/||b|| 2.199138400248e-02
    3 KSP unpreconditioned resid norm 2.875644920770e+02 true resid norm 2.875651934921e+02 ||r(i)||/||b|| 3.420657429016e-03
    4 KSP unpreconditioned resid norm 6.849768583136e+01 true resid norm 6.849653481076e+01 ||r(i)||/||b|| 8.147828247813e-04
    5 KSP unpreconditioned resid norm 1.473485650206e+01 true resid norm 1.473502409620e+01 ||r(i)||/||b|| 1.752766704110e-04
    6 KSP unpreconditioned resid norm 2.453059687780e+00 true resid norm 2.453465495780e+00 ||r(i)||/||b|| 2.918456463058e-05
    7 KSP unpreconditioned resid norm 4.858755826419e-01 true resid norm 4.839495592073e-01 ||r(i)||/||b|| 5.756696889734e-06
    8 KSP unpreconditioned resid norm 9.650128940861e-02 true resid norm 9.606247269307e-02 ||r(i)||/||b|| 1.142686313586e-06
    9 KSP unpreconditioned resid norm 9.766417235320e-03 true resid norm 1.061441998593e-02 ||r(i)||/||b|| 1.262610893156e-07
      Line search: Using full step: fnorm 8.406722960703e+04 gnorm 1.397410995346e+04
 1 Nonlinear |R| = 1.397411e+04
    0 KSP unpreconditioned resid norm 1.397410995346e+04 true resid norm 1.397410995346e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.458818894742e+02 true resid norm 6.458818894742e+02 ||r(i)||/||b|| 4.621989462122e-02
    2 KSP unpreconditioned resid norm 1.750015489848e+01 true resid norm 1.750018374089e+01 ||r(i)||/||b|| 1.252329042721e-03
    3 KSP unpreconditioned resid norm 4.364680996955e-01 true resid norm 4.363753479090e-01 ||r(i)||/||b|| 3.122741622632e-05
    4 KSP unpreconditioned resid norm 1.161082298597e-02 true resid norm 1.150579913975e-02 ||r(i)||/||b|| 8.233654363730e-07
      Line search: Using full step: fnorm 1.397410995346e+04 gnorm 2.791831308563e+02
 2 Nonlinear |R| = 2.791831e+02
    0 KSP unpreconditioned resid norm 2.791831308563e+02 true resid norm 2.791831308563e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.058614491817e-01 true resid norm 5.058614491817e-01 ||r(i)||/||b|| 1.811934151000e-03
    2 KSP unpreconditioned resid norm 5.281297736154e-04 true resid norm 5.290469694791e-04 ||r(i)||/||b|| 1.894981863182e-06
    3 KSP unpreconditioned resid norm 1.706597078703e-07 true resid norm 8.699121438243e-06 ||r(i)||/||b|| 3.115919436666e-08
      Line search: Using full step: fnorm 2.791831308563e+02 gnorm 1.548019303692e-01
 3 Nonlinear |R| = 1.548019e-01
    0 KSP unpreconditioned resid norm 1.548019303692e-01 true resid norm 1.548019303692e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.597331911756e-08 true resid norm 4.597331911776e-08 ||r(i)||/||b|| 2.969815622332e-07
      Line search: Using full step: fnorm 1.548019303692e-01 gnorm 7.127719337559e-08
 4 Nonlinear |R| = 7.127719e-08
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
+----------------+------------------+----------------+----------------+


Time Step 15, time = 1.5
                dt = 0.1
 0 Nonlinear |R| = 8.840340e+04
    0 KSP unpreconditioned resid norm 8.840339652876e+04 true resid norm 8.840339652876e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.249698933810e+04 true resid norm 1.249698933810e+04 ||r(i)||/||b|| 1.413632261747e-01
    2 KSP unpreconditioned resid norm 2.102460257834e+03 true resid norm 2.102458128437e+03 ||r(i)||/||b|| 2.378254921181e-02
    3 KSP unpreconditioned resid norm 3.493305970481e+02 true resid norm 3.493303872244e+02 ||r(i)||/||b|| 3.951549385444e-03
    4 KSP unpreconditioned resid norm 7.602965341080e+01 true resid norm 7.603109947513e+01 ||r(i)||/||b|| 8.600472658355e-04
    5 KSP unpreconditioned resid norm 1.574530646804e+01 true resid norm 1.574592195796e+01 ||r(i)||/||b|| 1.781144455557e-04
    6 KSP unpreconditioned resid norm 2.500432928514e+00 true resid norm 2.501710106022e+00 ||r(i)||/||b|| 2.829880077298e-05
    7 KSP unpreconditioned resid norm 4.424466449341e-01 true resid norm 4.427668217785e-01 ||r(i)||/||b|| 5.008482017255e-06
    8 KSP unpreconditioned resid norm 9.000319749100e-02 true resid norm 8.991566126940e-02 ||r(i)||/||b|| 1.017106409935e-06
    9 KSP unpreconditioned resid norm 1.094153140457e-02 true resid norm 1.137706628110e-02 ||r(i)||/||b|| 1.286949000585e-07
      Line search: Using full step: fnorm 8.840339652876e+04 gnorm 1.435706610102e+04
 1 Nonlinear |R| = 1.435707e+04
    0 KSP unpreconditioned resid norm 1.435706610102e+04 true resid norm 1.435706610102e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.346735769804e+02 true resid norm 8.346735769804e+02 ||r(i)||/||b|| 5.813677885909e-02
    2 KSP unpreconditioned resid norm 2.408273258217e+01 true resid norm 2.408276156212e+01 ||r(i)||/||b|| 1.677415245752e-03
    3 KSP unpreconditioned resid norm 7.019501928333e-01 true resid norm 7.019376240455e-01 ||r(i)||/||b|| 4.889143917750e-05
    4 KSP unpreconditioned resid norm 1.794242989170e-02 true resid norm 1.799758100210e-02 ||r(i)||/||b|| 1.253569557698e-06
    5 KSP unpreconditioned resid norm 5.409265117447e-04 true resid norm 8.688216370057e-04 ||r(i)||/||b|| 6.051526341749e-08
      Line search: Using full step: fnorm 1.435706610102e+04 gnorm 2.747215431080e+02
 2 Nonlinear |R| = 2.747215e+02
    0 KSP unpreconditioned resid norm 2.747215431080e+02 true resid norm 2.747215431080e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.946123270333e-01 true resid norm 7.946123270333e-01 ||r(i)||/||b|| 2.892428158505e-03
    2 KSP unpreconditioned resid norm 9.694189464112e-04 true resid norm 9.714141255831e-04 ||r(i)||/||b|| 3.535995446856e-06
    3 KSP unpreconditioned resid norm 3.749190911127e-07 true resid norm 3.061339485378e-05 ||r(i)||/||b|| 1.114342708892e-07
      Line search: Using full step: fnorm 2.747215431080e+02 gnorm 2.321169727744e-01
 3 Nonlinear |R| = 2.321170e-01
    0 KSP unpreconditioned resid norm 2.321169727744e-01 true resid norm 2.321169727744e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.145197895123e-07 true resid norm 1.145197895239e-07 ||r(i)||/||b|| 4.933710282153e-07
      Line search: Using full step: fnorm 2.321169727744e-01 gnorm 1.474165537313e-07
 4 Nonlinear |R| = 1.474166e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
:                :                  :                :                :
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442327e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
|   1.500000e+00 |    -0.000000e+00 |   4.000000e+00 |  -4.117662e-02 |
+----------------+------------------+----------------+----------------+


Time Step 16, time = 1.6
                dt = 0.1
 0 Nonlinear |R| = 9.062128e+04
    0 KSP unpreconditioned resid norm 9.062127675539e+04 true resid norm 9.062127675539e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.254007466929e+04 true resid norm 1.254007466929e+04 ||r(i)||/||b|| 1.383789229006e-01
    2 KSP unpreconditioned resid norm 2.222862573038e+03 true resid norm 2.222860512492e+03 ||r(i)||/||b|| 2.452912375635e-02
    3 KSP unpreconditioned resid norm 3.757691674665e+02 true resid norm 3.757689298834e+02 ||r(i)||/||b|| 4.146586136694e-03
    4 KSP unpreconditioned resid norm 7.909839029692e+01 true resid norm 7.909727837602e+01 ||r(i)||/||b|| 8.728334140504e-04
    5 KSP unpreconditioned resid norm 1.591845967357e+01 true resid norm 1.591839924729e+01 ||r(i)||/||b|| 1.756585188074e-04
    6 KSP unpreconditioned resid norm 2.487295175736e+00 true resid norm 2.486908362049e+00 ||r(i)||/||b|| 2.744287490853e-05
    7 KSP unpreconditioned resid norm 4.131037600097e-01 true resid norm 4.125812946288e-01 ||r(i)||/||b|| 4.552808230042e-06
    8 KSP unpreconditioned resid norm 8.500547435371e-02 true resid norm 8.422165878181e-02 ||r(i)||/||b|| 9.293806244769e-07
      Line search: Using full step: fnorm 9.062127675539e+04 gnorm 1.457507614177e+04
 1 Nonlinear |R| = 1.457508e+04
    0 KSP unpreconditioned resid norm 1.457507614177e+04 true resid norm 1.457507614177e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.214495797611e+02 true resid norm 9.214495797611e+02 ||r(i)||/||b|| 6.322091018932e-02
    2 KSP unpreconditioned resid norm 2.757621528820e+01 true resid norm 2.757622478458e+01 ||r(i)||/||b|| 1.892012399548e-03
    3 KSP unpreconditioned resid norm 8.448218950142e-01 true resid norm 8.448192586536e-01 ||r(i)||/||b|| 5.796328269138e-05
    4 KSP unpreconditioned resid norm 2.109663319602e-02 true resid norm 2.115269925696e-02 ||r(i)||/||b|| 1.451292538798e-06
    5 KSP unpreconditioned resid norm 6.732170718851e-04 true resid norm 7.604154218082e-04 ||r(i)||/||b|| 5.217231213147e-08
      Line search: Using full step: fnorm 1.457507614177e+04 gnorm 2.796043778277e+02
 2 Nonlinear |R| = 2.796044e+02
    0 KSP unpreconditioned resid norm 2.796043778277e+02 true resid norm 2.796043778277e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.527648405815e-01 true resid norm 9.527648405815e-01 ||r(i)||/||b|| 3.407546219354e-03
    2 KSP unpreconditioned resid norm 1.232892907370e-03 true resid norm 1.238837173984e-03 ||r(i)||/||b|| 4.430678745480e-06
    3 KSP unpreconditioned resid norm 5.522871934544e-07 true resid norm 2.339868972296e-05 ||r(i)||/||b|| 8.368499057402e-08
      Line search: Using full step: fnorm 2.796043778277e+02 gnorm 2.776849383579e-01
 3 Nonlinear |R| = 2.776849e-01
    0 KSP unpreconditioned resid norm 2.776849383579e-01 true resid norm 2.776849383579e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.844500732800e-07 true resid norm 2.844500732715e-07 ||r(i)||/||b|| 1.024362628213e-06
    2 KSP unpreconditioned resid norm 1.189244195367e-13 true resid norm 7.068606712574e-09 ||r(i)||/||b|| 2.545549194844e-08
      Line search: Using full step: fnorm 2.776849383579e-01 gnorm 3.984580878345e-07
 4 Nonlinear |R| = 3.984581e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
:                :                  :                :                :
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305756e-03 |
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
|   1.500000e+00 |    -0.000000e+00 |   4.000000e+00 |  -4.117662e-02 |
|   1.600000e+00 |    -0.000000e+00 |   4.000000e+00 |  -2.533803e-02 |
+----------------+------------------+----------------+----------------+


Time Step 17, time = 1.7
                dt = 0.1
 0 Nonlinear |R| = 8.967540e+04
    0 KSP unpreconditioned resid norm 8.967539566567e+04 true resid norm 8.967539566567e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.251959579663e+04 true resid norm 1.251959579663e+04 ||r(i)||/||b|| 1.396101539748e-01
    2 KSP unpreconditioned resid norm 2.172030853241e+03 true resid norm 2.172029663036e+03 ||r(i)||/||b|| 2.422102123902e-02
    3 KSP unpreconditioned resid norm 3.648571256615e+02 true resid norm 3.648551982590e+02 ||r(i)||/||b|| 4.068621003015e-03
    4 KSP unpreconditioned resid norm 7.782957702505e+01 true resid norm 7.782767763004e+01 ||r(i)||/||b|| 8.678821771827e-04
    5 KSP unpreconditioned resid norm 1.586738065601e+01 true resid norm 1.586516557778e+01 ||r(i)||/||b|| 1.769177092558e-04
    6 KSP unpreconditioned resid norm 2.494545436994e+00 true resid norm 2.495486259417e+00 ||r(i)||/||b|| 2.782799273862e-05
    7 KSP unpreconditioned resid norm 4.260332363084e-01 true resid norm 4.255576951118e-01 ||r(i)||/||b|| 4.745534624663e-06
    8 KSP unpreconditioned resid norm 8.725450200180e-02 true resid norm 8.712659399564e-02 ||r(i)||/||b|| 9.715774694819e-07
      Line search: Using full step: fnorm 8.967539566567e+04 gnorm 1.448042983542e+04
 1 Nonlinear |R| = 1.448043e+04
    0 KSP unpreconditioned resid norm 1.448042983542e+04 true resid norm 1.448042983542e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.850526641282e+02 true resid norm 8.850526641282e+02 ||r(i)||/||b|| 6.112060720487e-02
    2 KSP unpreconditioned resid norm 2.607375549432e+01 true resid norm 2.607397080643e+01 ||r(i)||/||b|| 1.800635140170e-03
    3 KSP unpreconditioned resid norm 7.831562935550e-01 true resid norm 7.831951353443e-01 ||r(i)||/||b|| 5.408645628935e-05
    4 KSP unpreconditioned resid norm 1.974709626414e-02 true resid norm 1.979640512882e-02 ||r(i)||/||b|| 1.367114467859e-06
    5 KSP unpreconditioned resid norm 6.154601348281e-04 true resid norm 5.401270437079e-04 ||r(i)||/||b|| 3.730048416013e-08
      Line search: Using full step: fnorm 1.448042983542e+04 gnorm 2.768587267926e+02
 2 Nonlinear |R| = 2.768587e+02
    0 KSP unpreconditioned resid norm 2.768587267926e+02 true resid norm 2.768587267926e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.844830754971e-01 true resid norm 8.844830754970e-01 ||r(i)||/||b|| 3.194709033534e-03
    2 KSP unpreconditioned resid norm 1.117073985025e-03 true resid norm 1.123999200123e-03 ||r(i)||/||b|| 4.059829405213e-06
    3 KSP unpreconditioned resid norm 4.685108991060e-07 true resid norm 5.750050639264e-05 ||r(i)||/||b|| 2.076889793534e-07
      Line search: Using full step: fnorm 2.768587267926e+02 gnorm 2.577862248883e-01
 3 Nonlinear |R| = 2.577862e-01
    0 KSP unpreconditioned resid norm 2.577862248883e-01 true resid norm 2.577862248883e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.314664264830e-07 true resid norm 2.314664264757e-07 ||r(i)||/||b|| 8.979006794329e-07
      Line search: Using full step: fnorm 2.577862248883e-01 gnorm 2.432903294579e-07
 4 Nonlinear |R| = 2.432903e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
:                :                  :                :                :
|   3.000000e-01 |     8.068667e+03 |   1.200000e+01 |   1.613733e-04 |
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
|   1.500000e+00 |    -0.000000e+00 |   4.000000e+00 |  -4.117662e-02 |
|   1.600000e+00 |    -0.000000e+00 |   4.000000e+00 |  -2.533803e-02 |
|   1.700000e+00 |    -0.000000e+00 |   4.000000e+00 |  -1.023547e-02 |
+----------------+------------------+----------------+----------------+


Time Step 18, time = 1.8
                dt = 0.1
 0 Nonlinear |R| = 8.600535e+04
    0 KSP unpreconditioned resid norm 8.600535391409e+04 true resid norm 8.600535391409e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.247243680461e+04 true resid norm 1.247243680461e+04 ||r(i)||/||b|| 1.450193067872e-01
    2 KSP unpreconditioned resid norm 1.966133734836e+03 true resid norm 1.966133810758e+03 ||r(i)||/||b|| 2.286059787303e-02
    3 KSP unpreconditioned resid norm 3.170448589100e+02 true resid norm 3.170460501972e+02 ||r(i)||/||b|| 3.686352485846e-03
    4 KSP unpreconditioned resid norm 7.222196736907e+01 true resid norm 7.222068811922e+01 ||r(i)||/||b|| 8.397231664364e-04
    5 KSP unpreconditioned resid norm 1.531905943043e+01 true resid norm 1.531860151787e+01 ||r(i)||/||b|| 1.781121851225e-04
    6 KSP unpreconditioned resid norm 2.491843004124e+00 true resid norm 2.490824045365e+00 ||r(i)||/||b|| 2.896126731660e-05
    7 KSP unpreconditioned resid norm 4.694214281965e-01 true resid norm 4.699544258760e-01 ||r(i)||/||b|| 5.464246171761e-06
    8 KSP unpreconditioned resid norm 9.423534695814e-02 true resid norm 9.410833349470e-02 ||r(i)||/||b|| 1.094214827471e-06
    9 KSP unpreconditioned resid norm 1.033090186850e-02 true resid norm 1.324656506129e-02 ||r(i)||/||b|| 1.540202377926e-07
      Line search: gnorm after quadratic fit 7.139431172831e+04
      Line search: Quadratically determined step, lambda=1.7076319954242838e-01
 1 Nonlinear |R| = 7.139431e+04
    0 KSP unpreconditioned resid norm 7.139431172831e+04 true resid norm 7.139431172831e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.632306212535e+03 true resid norm 8.632306212535e+03 ||r(i)||/||b|| 1.209102798747e-01
    2 KSP unpreconditioned resid norm 9.979514258480e+02 true resid norm 9.979527164440e+02 ||r(i)||/||b|| 1.397804239981e-02
    3 KSP unpreconditioned resid norm 1.261138320708e+02 true resid norm 1.261143585445e+02 ||r(i)||/||b|| 1.766448271459e-03
    4 KSP unpreconditioned resid norm 2.445898788709e+01 true resid norm 2.445879822524e+01 ||r(i)||/||b|| 3.425874923806e-04
    5 KSP unpreconditioned resid norm 3.263023393402e+00 true resid norm 3.261932778001e+00 ||r(i)||/||b|| 4.568897295928e-05
    6 KSP unpreconditioned resid norm 3.930383022330e-01 true resid norm 3.936380636365e-01 ||r(i)||/||b|| 5.513577400038e-06
    7 KSP unpreconditioned resid norm 3.553968876840e-02 true resid norm 3.479515467203e-02 ||r(i)||/||b|| 4.873659235549e-07
      Line search: gnorm after quadratic fit 5.732002456346e+04
      Line search: Quadratically determined step, lambda=1.9923762622940991e-01
 2 Nonlinear |R| = 5.732002e+04
    0 KSP unpreconditioned resid norm 5.732002456346e+04 true resid norm 5.732002456346e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.588830101448e+03 true resid norm 5.588830101448e+03 ||r(i)||/||b|| 9.750222795632e-02
    2 KSP unpreconditioned resid norm 4.710751463640e+02 true resid norm 4.710752043601e+02 ||r(i)||/||b|| 8.218335702885e-03
    3 KSP unpreconditioned resid norm 5.047835762541e+01 true resid norm 5.047933011004e+01 ||r(i)||/||b|| 8.806578590028e-04
    4 KSP unpreconditioned resid norm 8.095692130780e+00 true resid norm 8.095546625035e+00 ||r(i)||/||b|| 1.412341792714e-04
    5 KSP unpreconditioned resid norm 8.570854322446e-01 true resid norm 8.564162517723e-01 ||r(i)||/||b|| 1.494096100437e-05
    6 KSP unpreconditioned resid norm 9.696321153735e-02 true resid norm 9.726247954690e-02 ||r(i)||/||b|| 1.696832481975e-06
    7 KSP unpreconditioned resid norm 6.906718839691e-03 true resid norm 9.779310849294e-03 ||r(i)||/||b|| 1.706089786906e-07
      Line search: gnorm after quadratic fit 4.527702247198e+04
      Line search: Quadratically determined step, lambda=2.1181548243462364e-01
 3 Nonlinear |R| = 4.527702e+04
    0 KSP unpreconditioned resid norm 4.527702247198e+04 true resid norm 4.527702247198e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.491880310341e+03 true resid norm 3.491880310341e+03 ||r(i)||/||b|| 7.712256945567e-02
    2 KSP unpreconditioned resid norm 2.190888744152e+02 true resid norm 2.190890359210e+02 ||r(i)||/||b|| 4.838856973349e-03
    3 KSP unpreconditioned resid norm 1.930579642581e+01 true resid norm 1.930475365954e+01 ||r(i)||/||b|| 4.263697700415e-04
    4 KSP unpreconditioned resid norm 2.571458545288e+00 true resid norm 2.570818974599e+00 ||r(i)||/||b|| 5.677977115633e-05
    5 KSP unpreconditioned resid norm 2.240854889628e-01 true resid norm 2.236045853024e-01 ||r(i)||/||b|| 4.938588562019e-06
    6 KSP unpreconditioned resid norm 2.129856559535e-02 true resid norm 2.121043196498e-02 ||r(i)||/||b|| 4.684590727692e-07
      Line search: gnorm after quadratic fit 3.612702807088e+04
      Line search: Quadratically determined step, lambda=2.0320633813813826e-01
 4 Nonlinear |R| = 3.612703e+04
    0 KSP unpreconditioned resid norm 3.612702807088e+04 true resid norm 3.612702807088e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.222377221516e+03 true resid norm 2.222377221516e+03 ||r(i)||/||b|| 6.151563912635e-02
    2 KSP unpreconditioned resid norm 1.084503110636e+02 true resid norm 1.084501752490e+02 ||r(i)||/||b|| 3.001912447275e-03
    3 KSP unpreconditioned resid norm 7.847445769836e+00 true resid norm 7.847353936529e+00 ||r(i)||/||b|| 2.172155960665e-04
    4 KSP unpreconditioned resid norm 8.615892995886e-01 true resid norm 8.624214277454e-01 ||r(i)||/||b|| 2.387191733716e-05
    5 KSP unpreconditioned resid norm 6.183045329631e-02 true resid norm 6.148470396071e-02 ||r(i)||/||b|| 1.701903180081e-06
    6 KSP unpreconditioned resid norm 4.528478172022e-03 true resid norm 5.630802993217e-03 ||r(i)||/||b|| 1.558612289439e-07
      Line search: gnorm after quadratic fit 2.974748555542e+04
      Line search: Quadratically determined step, lambda=1.7720336206715326e-01
 5 Nonlinear |R| = 2.974749e+04
    0 KSP unpreconditioned resid norm 2.974748555542e+04 true resid norm 2.974748555542e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.505684162371e+03 true resid norm 1.505684162371e+03 ||r(i)||/||b|| 5.061551032829e-02
    2 KSP unpreconditioned resid norm 6.005144766335e+01 true resid norm 6.005197377846e+01 ||r(i)||/||b|| 2.018724361311e-03
    3 KSP unpreconditioned resid norm 3.650625785656e+00 true resid norm 3.650760232215e+00 ||r(i)||/||b|| 1.227250022666e-04
    4 KSP unpreconditioned resid norm 3.341790444292e-01 true resid norm 3.342828651782e-01 ||r(i)||/||b|| 1.123734860062e-05
    5 KSP unpreconditioned resid norm 1.998135787689e-02 true resid norm 2.052471410325e-02 ||r(i)||/||b|| 6.899646716361e-07
      Line search: gnorm after quadratic fit 2.537529270053e+04
      Line search: Quadratically determined step, lambda=1.4730384730063392e-01
 6 Nonlinear |R| = 2.537529e+04
    0 KSP unpreconditioned resid norm 2.537529270053e+04 true resid norm 2.537529270053e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.094916610638e+03 true resid norm 1.094916610638e+03 ||r(i)||/||b|| 4.314892535666e-02
    2 KSP unpreconditioned resid norm 3.722639582102e+01 true resid norm 3.722676108503e+01 ||r(i)||/||b|| 1.467047553870e-03
    3 KSP unpreconditioned resid norm 1.957094210192e+00 true resid norm 1.957375455330e+00 ||r(i)||/||b|| 7.713705920284e-05
    4 KSP unpreconditioned resid norm 1.532412253040e-01 true resid norm 1.523581706009e-01 ||r(i)||/||b|| 6.004193622474e-06
    5 KSP unpreconditioned resid norm 7.833506393997e-03 true resid norm 7.535521299762e-03 ||r(i)||/||b|| 2.969629311746e-07
      Line search: gnorm after quadratic fit 2.228854133594e+04
      Line search: Quadratically determined step, lambda=1.2182522798241417e-01
 7 Nonlinear |R| = 2.228854e+04
    0 KSP unpreconditioned resid norm 2.228854133594e+04 true resid norm 2.228854133594e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.443362082583e+02 true resid norm 8.443362082583e+02 ||r(i)||/||b|| 3.788207561599e-02
    2 KSP unpreconditioned resid norm 2.525239775610e+01 true resid norm 2.525219955678e+01 ||r(i)||/||b|| 1.132967796150e-03
    3 KSP unpreconditioned resid norm 1.177200512837e+00 true resid norm 1.177446008425e+00 ||r(i)||/||b|| 5.282741435064e-05
    4 KSP unpreconditioned resid norm 8.089623672614e-02 true resid norm 8.090489156018e-02 ||r(i)||/||b|| 3.629887229529e-06
    5 KSP unpreconditioned resid norm 3.628131301628e-03 true resid norm 4.003666804601e-03 ||r(i)||/||b|| 1.796289287960e-07
      Line search: gnorm after quadratic fit 2.020628923179e+04
      Line search: Quadratically determined step, lambda=1.0211729412451520e-01
 8 Nonlinear |R| = 2.020629e+04
    0 KSP unpreconditioned resid norm 2.020628923179e+04 true resid norm 2.020628923179e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.015878739924e+02 true resid norm 7.015878739924e+02 ||r(i)||/||b|| 3.472126257050e-02
    2 KSP unpreconditioned resid norm 1.817339541490e+01 true resid norm 1.817319281645e+01 ||r(i)||/||b|| 8.993829895227e-04
    3 KSP unpreconditioned resid norm 8.051157939585e-01 true resid norm 8.054583572400e-01 ||r(i)||/||b|| 3.986176521579e-05
    4 KSP unpreconditioned resid norm 4.538286721187e-02 true resid norm 4.541490259973e-02 ||r(i)||/||b|| 2.247562730532e-06
    5 KSP unpreconditioned resid norm 1.915346336655e-03 true resid norm 2.601918358813e-03 ||r(i)||/||b|| 1.287677479505e-07
      Line search: gnorm after quadratic fit 1.583004100499e+04
      Line search: Quadratically determined step, lambda=2.6307058544481471e-01
 9 Nonlinear |R| = 1.583004e+04
    0 KSP unpreconditioned resid norm 1.583004100499e+04 true resid norm 1.583004100499e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.390010040894e+02 true resid norm 4.390010095403e+02 ||r(i)||/||b|| 2.773214607606e-02
    2 KSP unpreconditioned resid norm 1.170274337059e+01 true resid norm 1.170213896977e+01 ||r(i)||/||b|| 7.392361754515e-04
    3 KSP unpreconditioned resid norm 4.574584328856e-01 true resid norm 4.574172803227e-01 ||r(i)||/||b|| 2.889552087569e-05
    4 KSP unpreconditioned resid norm 1.602677846952e-02 true resid norm 1.627553462799e-02 ||r(i)||/||b|| 1.028142291158e-06
    5 KSP unpreconditioned resid norm 4.765131945974e-04 true resid norm 2.877030152893e-03 ||r(i)||/||b|| 1.817449589667e-07
      Line search: gnorm after quadratic fit 1.239862218117e+04
      Line search: Quadratically determined step, lambda=2.2056355240930337e-01
10 Nonlinear |R| = 1.239862e+04
    0 KSP unpreconditioned resid norm 1.239862218117e+04 true resid norm 1.239862218117e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.813004839912e+02 true resid norm 2.813004886947e+02 ||r(i)||/||b|| 2.268804425075e-02
    2 KSP unpreconditioned resid norm 5.223963326458e+00 true resid norm 5.224560790722e+00 ||r(i)||/||b|| 4.213823692971e-04
    3 KSP unpreconditioned resid norm 1.353335731544e-01 true resid norm 1.357924001206e-01 ||r(i)||/||b|| 1.095221695898e-05
    4 KSP unpreconditioned resid norm 3.513762338983e-03 true resid norm 3.754576794085e-03 ||r(i)||/||b|| 3.028220990383e-07
      Line search: gnorm after quadratic fit 8.855228680737e+03
      Line search: Quadratically determined step, lambda=3.4237201987922633e-01
11 Nonlinear |R| = 8.855229e+03
    0 KSP unpreconditioned resid norm 8.855228680737e+03 true resid norm 8.855228680737e+03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.652910267699e+02 true resid norm 1.652910913498e+02 ||r(i)||/||b|| 1.866593142980e-02
    2 KSP unpreconditioned resid norm 2.628943921138e+00 true resid norm 2.629454794962e+00 ||r(i)||/||b|| 2.969381017435e-04
    3 KSP unpreconditioned resid norm 1.039716264736e-01 true resid norm 1.042884768692e-01 ||r(i)||/||b|| 1.177705067020e-05
    4 KSP unpreconditioned resid norm 1.799917871652e-03 true resid norm 1.652194706899e-03 ||r(i)||/||b|| 1.865784347831e-07
      Line search: Using full step: fnorm 8.855228680737e+03 gnorm 1.859448685901e+02
12 Nonlinear |R| = 1.859449e+02
    0 KSP unpreconditioned resid norm 1.859448685901e+02 true resid norm 1.859448685901e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.903901231493e-01 true resid norm 3.903914834523e-01 ||r(i)||/||b|| 2.099501246861e-03
    2 KSP unpreconditioned resid norm 2.823533357430e-03 true resid norm 2.823755696072e-03 ||r(i)||/||b|| 1.518598344489e-05
    3 KSP unpreconditioned resid norm 5.463783613473e-06 true resid norm 8.986694924450e-06 ||r(i)||/||b|| 4.832988935155e-08
      Line search: Using full step: fnorm 1.859448685901e+02 gnorm 2.201570529929e-01
13 Nonlinear |R| = 2.201571e-01
    0 KSP unpreconditioned resid norm 2.201570529929e-01 true resid norm 2.201570529929e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.313279478239e-04 true resid norm 4.313288700641e-04 ||r(i)||/||b|| 1.959187153900e-03
    2 KSP unpreconditioned resid norm 8.911217213245e-07 true resid norm 8.906593877985e-07 ||r(i)||/||b|| 4.045563726852e-06
    3 KSP unpreconditioned resid norm 2.941567525913e-09 true resid norm 7.166598107196e-09 ||r(i)||/||b|| 3.255220766163e-08
      Line search: Using full step: fnorm 2.201570529929e-01 gnorm 3.198993224263e-07
14 Nonlinear |R| = 3.198993e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
:                :                  :                :                :
|   4.000000e-01 |     9.417028e+03 |   5.000000e+00 |   1.890876e-04 |
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
|   1.500000e+00 |    -0.000000e+00 |   4.000000e+00 |  -4.117662e-02 |
|   1.600000e+00 |    -0.000000e+00 |   4.000000e+00 |  -2.533803e-02 |
|   1.700000e+00 |    -0.000000e+00 |   4.000000e+00 |  -1.023547e-02 |
|   1.800000e+00 |     1.951875e+03 |   1.400000e+01 |   3.903751e-05 |
+----------------+------------------+----------------+----------------+


Time Step 19, time = 1.9
                dt = 0.1
 0 Nonlinear |R| = 8.147513e+04
    0 KSP unpreconditioned resid norm 8.147513143066e+04 true resid norm 8.147513143066e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.252375677446e+04 true resid norm 1.252375676290e+04 ||r(i)||/||b|| 1.537126303817e-01
    2 KSP unpreconditioned resid norm 1.582178227819e+03 true resid norm 1.582174422641e+03 ||r(i)||/||b|| 1.941910856550e-02
    3 KSP unpreconditioned resid norm 2.083397267537e+02 true resid norm 2.083399900084e+02 ||r(i)||/||b|| 2.557099158357e-03
    4 KSP unpreconditioned resid norm 4.940829848825e+01 true resid norm 4.940811272094e+01 ||r(i)||/||b|| 6.064195522408e-04
    5 KSP unpreconditioned resid norm 8.141971212385e+00 true resid norm 8.141985980654e+00 ||r(i)||/||b|| 9.993216135629e-05
    6 KSP unpreconditioned resid norm 1.179531127873e+00 true resid norm 1.179377435289e+00 ||r(i)||/||b|| 1.447530571083e-05
    7 KSP unpreconditioned resid norm 2.034415433327e-01 true resid norm 2.020550344064e-01 ||r(i)||/||b|| 2.479959600659e-06
    8 KSP unpreconditioned resid norm 2.829197325752e-02 true resid norm 3.248390967401e-02 ||r(i)||/||b|| 3.986972356301e-07
      Line search: Using full step: fnorm 8.147513143066e+04 gnorm 2.674853541110e+04
 1 Nonlinear |R| = 2.674854e+04
    0 KSP unpreconditioned resid norm 2.674853541110e+04 true resid norm 2.674853541110e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.125206393056e+02 true resid norm 5.125206886608e+02 ||r(i)||/||b|| 1.916070097985e-02
    2 KSP unpreconditioned resid norm 4.065848800285e+01 true resid norm 4.065831847718e+01 ||r(i)||/||b|| 1.520020361949e-03
    3 KSP unpreconditioned resid norm 2.081317304929e+00 true resid norm 2.081389667025e+00 ||r(i)||/||b|| 7.781321986555e-05
    4 KSP unpreconditioned resid norm 2.690592872603e-01 true resid norm 2.688429488894e-01 ||r(i)||/||b|| 1.005075398550e-05
    5 KSP unpreconditioned resid norm 1.265549253775e-02 true resid norm 1.276752357785e-02 ||r(i)||/||b|| 4.773167346034e-07
      Line search: gnorm after quadratic fit 2.384670959104e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 2 Nonlinear |R| = 2.384671e+04
    0 KSP unpreconditioned resid norm 2.384670959104e+04 true resid norm 2.384670959104e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.093694304463e+02 true resid norm 4.093695064621e+02 ||r(i)||/||b|| 1.716670825798e-02
    2 KSP unpreconditioned resid norm 1.540066538867e+01 true resid norm 1.540071573842e+01 ||r(i)||/||b|| 6.458214153038e-04
    3 KSP unpreconditioned resid norm 7.480708735492e-01 true resid norm 7.480260731359e-01 ||r(i)||/||b|| 3.136810427787e-05
    4 KSP unpreconditioned resid norm 3.237898213414e-02 true resid norm 3.224319877643e-02 ||r(i)||/||b|| 1.352102630903e-06
    5 KSP unpreconditioned resid norm 1.922552459636e-03 true resid norm 1.963393605680e-03 ||r(i)||/||b|| 8.233394205538e-08
      Line search: Using full step: fnorm 2.384670959104e+04 gnorm 2.560901093190e+02
 3 Nonlinear |R| = 2.560901e+02
    0 KSP unpreconditioned resid norm 2.560901093190e+02 true resid norm 2.560901093190e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.384681036241e+00 true resid norm 5.384686562373e+00 ||r(i)||/||b|| 2.102653076564e-02
    2 KSP unpreconditioned resid norm 2.714483311648e-02 true resid norm 2.714134526426e-02 ||r(i)||/||b|| 1.059835748301e-04
    3 KSP unpreconditioned resid norm 6.024618850248e-04 true resid norm 6.050851162736e-04 ||r(i)||/||b|| 2.362782060902e-06
    4 KSP unpreconditioned resid norm 5.213074346678e-06 true resid norm 1.635614267389e-05 ||r(i)||/||b|| 6.386870120592e-08
      Line search: Using full step: fnorm 2.560901093190e+02 gnorm 5.384885219927e-01
 4 Nonlinear |R| = 5.384885e-01
    0 KSP unpreconditioned resid norm 5.384885219927e-01 true resid norm 5.384885219927e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.004445116048e-03 true resid norm 1.004445105294e-03 ||r(i)||/||b|| 1.865304578037e-03
    2 KSP unpreconditioned resid norm 1.362056104711e-05 true resid norm 1.362053277463e-05 ||r(i)||/||b|| 2.529400761269e-05
    3 KSP unpreconditioned resid norm 1.057910733238e-07 true resid norm 1.050709353679e-07 ||r(i)||/||b|| 1.951219591072e-07
      Line search: Using full step: fnorm 5.384885219927e-01 gnorm 2.211084371459e-07
 5 Nonlinear |R| = 2.211084e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
:                :                  :                :                :
|   5.000000e-01 |     6.944159e+03 |   6.000000e+00 |   1.397839e-04 |
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
|   1.500000e+00 |    -0.000000e+00 |   4.000000e+00 |  -4.117662e-02 |
|   1.600000e+00 |    -0.000000e+00 |   4.000000e+00 |  -2.533803e-02 |
|   1.700000e+00 |    -0.000000e+00 |   4.000000e+00 |  -1.023547e-02 |
|   1.800000e+00 |     1.951875e+03 |   1.400000e+01 |   3.903751e-05 |
|   1.900000e+00 |     8.369730e+03 |   5.000000e+00 |   1.675658e-04 |
+----------------+------------------+----------------+----------------+


Time Step 20, time = 2
                dt = 0.1
 0 Nonlinear |R| = 7.865514e+04
    0 KSP unpreconditioned resid norm 7.865513915954e+04 true resid norm 7.865513915954e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.279514702153e+04 true resid norm 1.279514703151e+04 ||r(i)||/||b|| 1.626740117459e-01
    2 KSP unpreconditioned resid norm 1.432307585767e+03 true resid norm 1.432306729293e+03 ||r(i)||/||b|| 1.820995734796e-02
    3 KSP unpreconditioned resid norm 1.729723964155e+02 true resid norm 1.729733560431e+02 ||r(i)||/||b|| 2.199136100850e-03
    4 KSP unpreconditioned resid norm 4.071697589879e+01 true resid norm 4.071782558522e+01 ||r(i)||/||b|| 5.176753359069e-04
    5 KSP unpreconditioned resid norm 6.675065620597e+00 true resid norm 6.676933045346e+00 ||r(i)||/||b|| 8.488870678625e-05
    6 KSP unpreconditioned resid norm 1.022513728553e+00 true resid norm 1.021474817656e+00 ||r(i)||/||b|| 1.298675240513e-05
    7 KSP unpreconditioned resid norm 1.609975774788e-01 true resid norm 1.605269085962e-01 ||r(i)||/||b|| 2.040895360576e-06
    8 KSP unpreconditioned resid norm 2.354700614502e-02 true resid norm 2.783217200802e-02 ||r(i)||/||b|| 3.538506486088e-07
      Line search: Using full step: fnorm 7.865513915954e+04 gnorm 1.732436078119e+04
 1 Nonlinear |R| = 1.732436e+04
    0 KSP unpreconditioned resid norm 1.732436078119e+04 true resid norm 1.732436078119e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.300170760189e+02 true resid norm 5.300170732904e+02 ||r(i)||/||b|| 3.059374484199e-02
    2 KSP unpreconditioned resid norm 5.172709563772e+01 true resid norm 5.172736073262e+01 ||r(i)||/||b|| 2.985816411118e-03
    3 KSP unpreconditioned resid norm 2.475134218153e+00 true resid norm 2.475245952040e+00 ||r(i)||/||b|| 1.428766107623e-04
    4 KSP unpreconditioned resid norm 1.578148007012e-01 true resid norm 1.575353819665e-01 ||r(i)||/||b|| 9.093286843663e-06
    5 KSP unpreconditioned resid norm 3.928488332407e-03 true resid norm 4.226476805274e-03 ||r(i)||/||b|| 2.439614862941e-07
      Line search: gnorm after quadratic fit 1.960257313299e+04
      Line search: Cubically determined step, current gnorm 1.661321571460e+04 lambda=2.6520922705359200e-02
 2 Nonlinear |R| = 1.661322e+04
    0 KSP unpreconditioned resid norm 1.661321571460e+04 true resid norm 1.661321571460e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.803353999396e+02 true resid norm 4.803354972918e+02 ||r(i)||/||b|| 2.891285501515e-02
    2 KSP unpreconditioned resid norm 2.947765304170e+01 true resid norm 2.947796929048e+01 ||r(i)||/||b|| 1.774368659077e-03
    3 KSP unpreconditioned resid norm 8.901137719777e-01 true resid norm 8.900588568097e-01 ||r(i)||/||b|| 5.357535061844e-05
    4 KSP unpreconditioned resid norm 3.865717738534e-02 true resid norm 3.845736435559e-02 ||r(i)||/||b|| 2.314865768088e-06
    5 KSP unpreconditioned resid norm 1.833756485926e-03 true resid norm 2.388088896545e-03 ||r(i)||/||b|| 1.437463365052e-07
      Line search: gnorm after quadratic fit 1.944426309344e+04
      Line search: Cubically determined step, current gnorm 1.591092891617e+04 lambda=1.9398456357208514e-02
 3 Nonlinear |R| = 1.591093e+04
    0 KSP unpreconditioned resid norm 1.591092891617e+04 true resid norm 1.591092891617e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.873725469860e+02 true resid norm 3.873725636089e+02 ||r(i)||/||b|| 2.434631979376e-02
    2 KSP unpreconditioned resid norm 1.865652571698e+01 true resid norm 1.865663544878e+01 ||r(i)||/||b|| 1.172567330737e-03
    3 KSP unpreconditioned resid norm 5.245705999030e-01 true resid norm 5.245281953319e-01 ||r(i)||/||b|| 3.296653502102e-05
    4 KSP unpreconditioned resid norm 4.326827980163e-02 true resid norm 4.334447624700e-02 ||r(i)||/||b|| 2.724195204150e-06
    5 KSP unpreconditioned resid norm 1.376836543432e-03 true resid norm 1.663641698541e-03 ||r(i)||/||b|| 1.045596839321e-07
      Line search: Using full step: fnorm 1.591092891617e+04 gnorm 3.413764264440e+02
 4 Nonlinear |R| = 3.413764e+02
    0 KSP unpreconditioned resid norm 3.413764264440e+02 true resid norm 3.413764264440e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.111336444441e+00 true resid norm 4.111340040706e+00 ||r(i)||/||b|| 1.204342105145e-02
    2 KSP unpreconditioned resid norm 3.590950915354e-02 true resid norm 3.591099386441e-02 ||r(i)||/||b|| 1.051947090737e-04
    3 KSP unpreconditioned resid norm 5.853196524100e-04 true resid norm 5.887717083451e-04 ||r(i)||/||b|| 1.724699372122e-06
    4 KSP unpreconditioned resid norm 3.574166648914e-06 true resid norm 2.607699534060e-05 ||r(i)||/||b|| 7.638780337656e-08
      Line search: Using full step: fnorm 3.413764264440e+02 gnorm 1.221583642166e+00
 5 Nonlinear |R| = 1.221584e+00
    0 KSP unpreconditioned resid norm 1.221583642166e+00 true resid norm 1.221583642166e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.523947334517e-03 true resid norm 3.523946492451e-03 ||r(i)||/||b|| 2.884736149710e-03
    2 KSP unpreconditioned resid norm 8.020136997107e-05 true resid norm 8.019989296042e-05 ||r(i)||/||b|| 6.565239594909e-05
    3 KSP unpreconditioned resid norm 5.621261881987e-07 true resid norm 5.693382376485e-07 ||r(i)||/||b|| 4.660657019269e-07
      Line search: Using full step: fnorm 1.221583642166e+00 gnorm 1.365679635734e-06
 6 Nonlinear |R| = 1.365680e-06
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
:                :                  :                :                :
|   6.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -2.981473e-03 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660047e-02 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233497e-02 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770082e-02 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027210e-02 |
|   1.100000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.806408e-02 |
|   1.200000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.984658e-02 |
|   1.300000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.533819e-02 |
|   1.400000e+00 |    -0.000000e+00 |   4.000000e+00 |  -5.525067e-02 |
|   1.500000e+00 |    -0.000000e+00 |   4.000000e+00 |  -4.117662e-02 |
|   1.600000e+00 |    -0.000000e+00 |   4.000000e+00 |  -2.533803e-02 |
|   1.700000e+00 |    -0.000000e+00 |   4.000000e+00 |  -1.023547e-02 |
|   1.800000e+00 |     1.951875e+03 |   1.400000e+01 |   3.903751e-05 |
|   1.900000e+00 |     8.369730e+03 |   5.000000e+00 |   1.675658e-04 |
|   2.000000e+00 |     1.046412e+04 |   6.000000e+00 |   2.104856e-04 |
+----------------+------------------+----------------+----------------+


Look adding in temperature things, and things are starting to look tough again!
Next task!

# 1/8/18

Failing linear solve:
 0 Nonlinear |R| = 1.602900e+02
    0 KSP unpreconditioned resid norm 1.602900149382e+02 true resid norm 1.602900149382e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.594651817350e+02 true resid norm 1.594651802172e+02 ||r(i)||/||b|| 9.948541103988e-01
    2 KSP unpreconditioned resid norm 1.586635739809e+02 true resid norm 1.587997593310e+02 ||r(i)||/||b|| 9.907027545801e-01
    3 KSP unpreconditioned resid norm 1.467847192430e+02 true resid norm 5.440916350914e+02 ||r(i)||/||b|| 3.394420016126e+00
    4 KSP unpreconditioned resid norm 1.455917192916e+02 true resid norm 6.206391743621e+02 ||r(i)||/||b|| 3.871976520817e+00
    5 KSP unpreconditioned resid norm 1.223033935014e+02 true resid norm 2.419546368594e+03 ||r(i)||/||b|| 1.509480406204e+01
    6 KSP unpreconditioned resid norm 8.813516974347e+01 true resid norm 1.882281778955e+03 ||r(i)||/||b|| 1.174297587832e+01
    7 KSP unpreconditioned resid norm 8.794997713634e+01 true resid norm 1.886432145929e+03 ||r(i)||/||b|| 1.176886873868e+01
    8 KSP unpreconditioned resid norm 8.427470829218e+01 true resid norm 2.936597377116e+03 ||r(i)||/||b|| 1.832052594323e+01
    9 KSP unpreconditioned resid norm 6.459510689520e+01 true resid norm 1.133823387262e+04 ||r(i)||/||b|| 7.073574655906e+01
   10 KSP unpreconditioned resid norm 4.710048084920e+01 true resid norm 4.747571870544e+04 ||r(i)||/||b|| 2.961863764486e+02
   11 KSP unpreconditioned resid norm 4.689195102232e+01 true resid norm 4.812481827343e+04 ||r(i)||/||b|| 3.002359085934e+02
   12 KSP unpreconditioned resid norm 3.992813317864e+01 true resid norm 9.267257062334e+04 ||r(i)||/||b|| 5.781556053823e+02
   13 KSP unpreconditioned resid norm 1.990642915114e+01 true resid norm 4.527871922676e+05 ||r(i)||/||b|| 2.824799738412e+03
   14 KSP unpreconditioned resid norm 2.348494303877e+00 true resid norm 4.654803075725e+05 ||r(i)||/||b|| 2.903988172638e+03
   15 KSP unpreconditioned resid norm 5.120050709300e-01 true resid norm 4.761141534357e+05 ||r(i)||/||b|| 2.970329459507e+03
   16 KSP unpreconditioned resid norm 2.982448174591e-01 true resid norm 4.755545133475e+05 ||r(i)||/||b|| 2.966838037484e+03
   17 KSP unpreconditioned resid norm 9.272306041853e-02 true resid norm 4.758549535672e+05 ||r(i)||/||b|| 2.968712391416e+03
   18 KSP unpreconditioned resid norm 1.234915147494e-03 true resid norm 4.757988528961e+05 ||r(i)||/||b|| 2.968362396619e+03

Bad residual:
  0.
  0.
  0.
  1.46237e-08
  1.31039e-06
  1.26006e-08
  9.89361e-09
  -56.6667
  1.78473e-06
  1.09714e-08
  3.47328e-09
  -56.6667
  0.
  0.
  0.
  1.46237e-08
  0.
  0.
  0.
  1.46237e-08
  7.07007e-07
  9.38382e-09
  1.88573e-09
  -56.6667
  9.51155e-07
  1.41793e-08
  8.31493e-09
  -56.6667
  0.
  0.
  0.
  1.46237e-08
  0.
  0.
  0.
  7.96613e-09
  0.999999
  -2.75764e-08
  -1.76383e-08
  -56.6667
  0.999999
  -2.71673e-08
  -1.38036e-08
  -56.6667
  0.
  0.
  0.
  7.96613e-09
  0.
  0.
  0.
  7.96613e-09
  0.999999
  -3.02493e-08
  -1.68856e-08
  -56.6667
  0.999999
  -3.06585e-08
  -1.45562e-08
  -56.6667
  0.
  0.
  0.
  7.96613e-09

b:
  0.
  0.
  0.
  -3.33956e-15
  -3.59712e-18
  -1.19904e-18
  1.19904e-18
  -56.6667
  -3.59712e-18
  -1.19904e-18
  -1.19904e-18
  -56.6667
  0.
  0.
  0.
  1.34116e-15
  0.
  0.
  0.
  -1.46162e-15
  -3.59712e-18
  1.19904e-18
  1.19904e-18
  -56.6667
  -3.59712e-18
  1.19904e-18
  -1.19904e-18
  -56.6667
  0.
  0.
  0.
  3.46002e-15
  0.
  0.
  0.
  -3.14155e-15
  1.
  -3.9968e-20
  3.9968e-20
  -56.6667
  1.
  3.9968e-20
  3.9968e-20
  -56.6667
  0.
  0.
  0.
  1.59813e-15
  0.
  0.
  0.
  2.20813e-15
  1.
  -3.9968e-20
  -3.9968e-20
  -56.6667
  1.
  3.9968e-20
  -3.9968e-20
  -56.6667
  0.
  0.
  0.

v:
  0.
  0.
  0.
  10.3057
  923.465
  8.87995
  6.97225
  -23.6671
  1257.74
  7.73177
  2.4477
  -23.6671
  0.
  0.
  0.
  10.3056
  0.
  0.
  0.
  10.3057
  498.244
  6.61299
  1.32892
  -23.6671
  670.3
  9.99249
  5.85972
  -23.6671
  0.
  0.
  0.
  10.3056
  0.
  0.
  0.
  5.61391
  -821.952
  -19.4337
  -12.4301
  -23.6671
  -821.778
  -19.1454
  -9.72772
  -23.6671
  0.
  0.
  0.
  5.61391
  0.
  0.
  0.
  5.61391
  -822.788
  -21.3173
  -11.8996
  -23.6671
  -821.997
  -21.6057
  -10.2581
  -23.6671
  0.
  0.
  0.
  5.61391

Solve that actually makes some progress (with PJFNK, FDP, -pc_type bjacobi):
 0 Nonlinear |R| = 1.602900e+02
    0 KSP unpreconditioned resid norm 1.602900149382e+02 true resid norm 1.602900149382e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.599671499368e+02 true resid norm 1.599671499378e+02 ||r(i)||/||b|| 9.979857447728e-01
    2 KSP unpreconditioned resid norm 1.599608879658e+02 true resid norm 1.610107178943e+02 ||r(i)||/||b|| 1.004496243614e+00
    3 KSP unpreconditioned resid norm 3.181879201229e+00 true resid norm 2.229717380990e+03 ||r(i)||/||b|| 1.391051951583e+01
    4 KSP unpreconditioned resid norm 1.750714636783e-02 true resid norm 2.223510353009e+03 ||r(i)||/||b|| 1.387179578133e+01
    5 KSP unpreconditioned resid norm 4.603842515535e-04 true resid norm 2.223496256476e+03 ||r(i)||/||b|| 1.387170783740e+01
  Linear solve converged due to CONVERGED_RTOL iterations 5
      Line search: Using full step: fnorm 1.602900149382e+02 gnorm 2.848229364226e-02
 1 Nonlinear |R| = 2.848229e-02
    0 KSP unpreconditioned resid norm 2.848229364226e-02 true resid norm 2.848229364226e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.847665574424e-02 true resid norm 2.848187577217e-02 ||r(i)||/||b|| 9.999853287765e-01
    2 KSP unpreconditioned resid norm 2.680115085788e-02 true resid norm 1.872178407740e-01 ||r(i)||/||b|| 6.573130771188e+00
    3 KSP unpreconditioned resid norm 2.007162991418e-02 true resid norm 2.176327964035e+00 ||r(i)||/||b|| 7.640985629070e+01
    4 KSP unpreconditioned resid norm 1.568616426536e-04 true resid norm 4.650601808486e+00 ||r(i)||/||b|| 1.632804530034e+02
    5 KSP unpreconditioned resid norm 1.466325219419e-04 true resid norm 4.649564242550e+00 ||r(i)||/||b|| 1.632440245490e+02
    6 KSP unpreconditioned resid norm 1.195816525192e-04 true resid norm 4.649111986576e+00 ||r(i)||/||b|| 1.632281460534e+02
    7 KSP unpreconditioned resid norm 1.055346315444e-04 true resid norm 4.648784124622e+00 ||r(i)||/||b|| 1.632166349737e+02
    8 KSP unpreconditioned resid norm 9.504592910450e-05 true resid norm 4.648604238169e+00 ||r(i)||/||b|| 1.632103192445e+02
    9 KSP unpreconditioned resid norm 8.725167536409e-05 true resid norm 4.648479427200e+00 ||r(i)||/||b|| 1.632059371898e+02
   10 KSP unpreconditioned resid norm 8.108499312439e-05 true resid norm 4.648393889082e+00 ||r(i)||/||b|| 1.632029339865e+02
   11 KSP unpreconditioned resid norm 7.606446422048e-05 true resid norm 4.648333111420e+00 ||r(i)||/||b|| 1.632008001112e+02
   12 KSP unpreconditioned resid norm 7.186965005058e-05 true resid norm 4.648290322227e+00 ||r(i)||/||b|| 1.631992978027e+02
   13 KSP unpreconditioned resid norm 6.829646207508e-05 true resid norm 4.648261026596e+00 ||r(i)||/||b|| 1.631982692468e+02
   14 KSP unpreconditioned resid norm 6.520439477987e-05 true resid norm 4.648242272456e+00 ||r(i)||/||b|| 1.631976107977e+02
   15 KSP unpreconditioned resid norm 6.249397178342e-05 true resid norm 4.648231874294e+00 ||r(i)||/||b|| 1.631972457231e+02
   16 KSP unpreconditioned resid norm 6.009243004956e-05 true resid norm 4.648228062221e+00 ||r(i)||/||b|| 1.631971118830e+02
   17 KSP unpreconditioned resid norm 5.794521577011e-05 true resid norm 4.648229312725e+00 ||r(i)||/||b|| 1.631971557876e+02
   18 KSP unpreconditioned resid norm 5.601051763586e-05 true resid norm 4.648234276869e+00 ||r(i)||/||b|| 1.631973300764e+02
   19 KSP unpreconditioned resid norm 5.425566058730e-05 true resid norm 4.648241776307e+00 ||r(i)||/||b|| 1.631975933781e+02
   20 KSP unpreconditioned resid norm 5.265465719768e-05 true resid norm 4.648250843643e+00 ||r(i)||/||b|| 1.631979117281e+02
   21 KSP unpreconditioned resid norm 5.118650792365e-05 true resid norm 4.648260710488e+00 ||r(i)||/||b|| 1.631982581484e+02
   22 KSP unpreconditioned resid norm 4.983400965297e-05 true resid norm 4.648270813056e+00 ||r(i)||/||b|| 1.631986128448e+02
   23 KSP unpreconditioned resid norm 4.858290021493e-05 true resid norm 4.648280774836e+00 ||r(i)||/||b|| 1.631989625982e+02
   24 KSP unpreconditioned resid norm 4.742123501403e-05 true resid norm 4.648290356909e+00 ||r(i)||/||b|| 1.631992990203e+02
   25 KSP unpreconditioned resid norm 4.633892489790e-05 true resid norm 4.648299429171e+00 ||r(i)||/||b|| 1.631996175432e+02
   26 KSP unpreconditioned resid norm 4.532738402369e-05 true resid norm 4.648307936198e+00 ||r(i)||/||b|| 1.631999162209e+02
   27 KSP unpreconditioned resid norm 4.437925707050e-05 true resid norm 4.648315869992e+00 ||r(i)||/||b|| 1.632001947727e+02
   28 KSP unpreconditioned resid norm 4.348820600102e-05 true resid norm 4.648323248438e+00 ||r(i)||/||b|| 1.632004538265e+02
   29 KSP unpreconditioned resid norm 4.264873778831e-05 true resid norm 4.648330105297e+00 ||r(i)||/||b|| 1.632006945676e+02
   30 KSP unpreconditioned resid norm 4.648336483097e+00 true resid norm 1.814175181757e+02 ||r(i)||/||b|| 6.369484159325e+03
   31 KSP unpreconditioned resid norm 1.186350170075e+00 true resid norm 2.997115815579e+02 ||r(i)||/||b|| 1.052273336278e+04
   32 KSP unpreconditioned resid norm 1.046281089570e+00 true resid norm 1.340946202100e+03 ||r(i)||/||b|| 4.707999358977e+04
   33 KSP unpreconditioned resid norm 5.011647603945e-03 true resid norm 1.451293794028e+02 ||r(i)||/||b|| 5.095424589943e+03
   34 KSP unpreconditioned resid norm 9.048268301703e-07 true resid norm 1.567998870127e+02 ||r(i)||/||b|| 5.505170650303e+03
   35 KSP unpreconditioned resid norm 7.664751845509e-11 true resid norm 1.568027504945e+02 ||r(i)||/||b|| 5.505271185807e+03
  Linear solve converged due to CONVERGED_RTOL iterations 35
      Line search: Using full step: fnorm 2.848229364226e-02 gnorm 1.652796105268e-02
 2 Nonlinear |R| = 1.652796e-02
    0 KSP unpreconditioned resid norm 1.652796105268e-02 true resid norm 1.652796105268e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.088113791034e-03 true resid norm 3.069173917373e-02 ||r(i)||/||b|| 1.856958585267e+00
    2 KSP unpreconditioned resid norm 7.081748071372e-03 true resid norm 1.104095439862e+01 ||r(i)||/||b|| 6.680167241095e+02
    3 KSP unpreconditioned resid norm 6.999321152146e-03 true resid norm 2.865385606521e+01 ||r(i)||/||b|| 1.733659461919e+03
    4 KSP unpreconditioned resid norm 6.919709285620e-03 true resid norm 6.599985343645e+01 ||r(i)||/||b|| 3.993224162744e+03
    5 KSP unpreconditioned resid norm 6.842753614241e-03 true resid norm 1.017442513934e+02 ||r(i)||/||b|| 6.155886444139e+03
    6 KSP unpreconditioned resid norm 6.768309651503e-03 true resid norm 1.359513499281e+02 ||r(i)||/||b|| 8.225536682645e+03
    7 KSP unpreconditioned resid norm 6.696243674625e-03 true resid norm 1.687128035879e+02 ||r(i)||/||b|| 1.020772030199e+04
    8 KSP unpreconditioned resid norm 6.626431719980e-03 true resid norm 2.001167843285e+02 ||r(i)||/||b|| 1.210777201680e+04
    9 KSP unpreconditioned resid norm 6.558758686416e-03 true resid norm 2.302453358645e+02 ||r(i)||/||b|| 1.393065576151e+04
   10 KSP unpreconditioned resid norm 6.493117536794e-03 true resid norm 2.591743272638e+02 ||r(i)||/||b|| 1.568096188258e+04
   11 KSP unpreconditioned resid norm 6.429408585956e-03 true resid norm 2.869738522180e+02 ||r(i)||/||b|| 1.736293129584e+04
   12 KSP unpreconditioned resid norm 6.367538863851e-03 true resid norm 3.137087054746e+02 ||r(i)||/||b|| 1.898048431229e+04
   13 KSP unpreconditioned resid norm 6.307421545034e-03 true resid norm 3.394388514720e+02 ||r(i)||/||b|| 2.053724899218e+04
   14 KSP unpreconditioned resid norm 6.248975436114e-03 true resid norm 3.642198482306e+02 ||r(i)||/||b|| 2.203658679191e+04
   15 KSP unpreconditioned resid norm 6.192124514575e-03 true resid norm 3.881032424737e+02 ||r(i)||/||b|| 2.348161647022e+04
   16 KSP unpreconditioned resid norm 6.136797512983e-03 true resid norm 4.111369204082e+02 ||r(i)||/||b|| 2.487523531171e+04
   17 KSP unpreconditioned resid norm 6.082927543272e-03 true resid norm 4.333654269722e+02 ||r(i)||/||b|| 2.622013844242e+04
   18 KSP unpreconditioned resid norm 6.030451756714e-03 true resid norm 4.548302477683e+02 ||r(i)||/||b|| 2.751883588778e+04
   19 KSP unpreconditioned resid norm 5.979311035592e-03 true resid norm 4.755700666639e+02 ||r(i)||/||b|| 2.877366815835e+04
   20 KSP unpreconditioned resid norm 5.929449713187e-03 true resid norm 4.956210040263e+02 ||r(i)||/||b|| 2.998682066388e+04
   21 KSP unpreconditioned resid norm 5.880815319013e-03 true resid norm 5.150168179349e+02 ||r(i)||/||b|| 3.116033588737e+04
   22 KSP unpreconditioned resid norm 5.833358346667e-03 true resid norm 5.337890964632e+02 ||r(i)||/||b|| 3.229612501880e+04
   23 KSP unpreconditioned resid norm 5.787032041957e-03 true resid norm 5.519674298924e+02 ||r(i)||/||b|| 3.339597837465e+04
   24 KSP unpreconditioned resid norm 5.741792209274e-03 true resid norm 5.695795652708e+02 ||r(i)||/||b|| 3.446157474932e+04
   25 KSP unpreconditioned resid norm 5.697597034311e-03 true resid norm 5.866515504207e+02 ||r(i)||/||b|| 3.549449012802e+04
   26 KSP unpreconditioned resid norm 5.654406921633e-03 true resid norm 6.032078571232e+02 ||r(i)||/||b|| 3.649620513992e+04
   27 KSP unpreconditioned resid norm 5.612184345518e-03 true resid norm 6.192715032059e+02 ||r(i)||/||b|| 3.746811244485e+04
   28 KSP unpreconditioned resid norm 5.570893712931e-03 true resid norm 6.348641528438e+02 ||r(i)||/||b|| 3.841152280189e+04
   29 KSP unpreconditioned resid norm 5.530501237386e-03 true resid norm 6.500062239378e+02 ||r(i)||/||b|| 3.932767156614e+04
   30 KSP unpreconditioned resid norm 6.647169687370e+02 true resid norm 6.647169687370e+02 ||r(i)||/||b|| 4.021772356665e+04
  Linear solve did not converge due to DIVERGED_DTOL iterations 30
      Line search: gnorm after quadratic fit 6.899607364361e-03
      Line search: Cubic step no good, shrinking lambda, current gnorm 6.712313357669e-03 lambda=2.1133200712480757e-01
      Line search: Cubically determined step, current gnorm 6.888729499309e-03 lambda=9.2205704831343577e-02
 3 Nonlinear |R| = 6.888729e-03
    0 KSP unpreconditioned resid norm 6.888729499309e-03 true resid norm 6.888729499309e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.570336917535e-06 true resid norm 7.570336917535e-06 ||r(i)||/||b|| 1.098945301640e-03
    2 KSP unpreconditioned resid norm 7.117100872228e-10 true resid norm 9.720593837431e-10 ||r(i)||/||b|| 1.411086592732e-07
  Linear solve converged due to CONVERGED_RTOL iterations 2
      Line search: Using full step: fnorm 6.888729499309e-03 gnorm 6.683782667303e-06
 4 Nonlinear |R| = 6.683783e-06
    0 KSP unpreconditioned resid norm 6.683782667303e-06 true resid norm 6.683782667303e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.671630640268e-06 true resid norm 6.671630640276e-06 ||r(i)||/||b|| 9.981818638289e-01
    2 KSP unpreconditioned resid norm 6.504446574071e-09 true resid norm 1.689031619885e+01 ||r(i)||/||b|| 2.527059457135e+06
    3 KSP unpreconditioned resid norm 3.575903788813e-09 true resid norm 1.689075259859e+01 ||r(i)||/||b|| 2.527124749466e+06
    4 KSP unpreconditioned resid norm 3.133135406431e-09 true resid norm 1.689065477890e+01 ||r(i)||/||b|| 2.527110114087e+06
    5 KSP unpreconditioned resid norm 2.822403664791e-09 true resid norm 1.689059401062e+01 ||r(i)||/||b|| 2.527101022188e+06
    6 KSP unpreconditioned resid norm 2.588907649210e-09 true resid norm 1.689055259444e+01 ||r(i)||/||b|| 2.527094825670e+06
    7 KSP unpreconditioned resid norm 2.405176849519e-09 true resid norm 1.689052255588e+01 ||r(i)||/||b|| 2.527090331424e+06
    8 KSP unpreconditioned resid norm 2.255723178669e-09 true resid norm 1.689049977274e+01 ||r(i)||/||b|| 2.527086922704e+06
    9 KSP unpreconditioned resid norm 2.131062809744e-09 true resid norm 1.689048189961e+01 ||r(i)||/||b|| 2.527084248601e+06
   10 KSP unpreconditioned resid norm 2.025020800417e-09 true resid norm 1.689046750370e+01 ||r(i)||/||b|| 2.527082094744e+06
   11 KSP unpreconditioned resid norm 1.933380880740e-09 true resid norm 1.689045566015e+01 ||r(i)||/||b|| 2.527080322760e+06
   12 KSP unpreconditioned resid norm 1.853152831515e-09 true resid norm 1.689044574555e+01 ||r(i)||/||b|| 2.527078839379e+06
   13 KSP unpreconditioned resid norm 1.782149156690e-09 true resid norm 1.689043732406e+01 ||r(i)||/||b|| 2.527077579390e+06
   14 KSP unpreconditioned resid norm 1.718727706859e-09 true resid norm 1.689043008204e+01 ||r(i)||/||b|| 2.527076495870e+06
   15 KSP unpreconditioned resid norm 1.661628471741e-09 true resid norm 1.689042378788e+01 ||r(i)||/||b|| 2.527075554163e+06
   16 KSP unpreconditioned resid norm 1.609866357327e-09 true resid norm 1.689041826701e+01 ||r(i)||/||b|| 2.527074728153e+06
   17 KSP unpreconditioned resid norm 1.562658580057e-09 true resid norm 1.689041338511e+01 ||r(i)||/||b|| 2.527073997743e+06
   18 KSP unpreconditioned resid norm 1.519374201138e-09 true resid norm 1.689040903731e+01 ||r(i)||/||b|| 2.527073347244e+06
   19 KSP unpreconditioned resid norm 1.479498240803e-09 true resid norm 1.689040514054e+01 ||r(i)||/||b|| 2.527072764226e+06
   20 KSP unpreconditioned resid norm 1.442605641569e-09 true resid norm 1.689040162805e+01 ||r(i)||/||b|| 2.527072238700e+06
   21 KSP unpreconditioned resid norm 1.408342037918e-09 true resid norm 1.689039844569e+01 ||r(i)||/||b|| 2.527071762570e+06
   22 KSP unpreconditioned resid norm 1.376409324779e-09 true resid norm 1.689039554899e+01 ||r(i)||/||b|| 2.527071329177e+06
   23 KSP unpreconditioned resid norm 1.346554670825e-09 true resid norm 1.689039290117e+01 ||r(i)||/||b|| 2.527070933020e+06
   24 KSP unpreconditioned resid norm 1.318562045224e-09 true resid norm 1.689039047142e+01 ||r(i)||/||b|| 2.527070569491e+06
   25 KSP unpreconditioned resid norm 1.292245605243e-09 true resid norm 1.689038823392e+01 ||r(i)||/||b|| 2.527070234726e+06
   26 KSP unpreconditioned resid norm 1.267444480486e-09 true resid norm 1.689038616671e+01 ||r(i)||/||b|| 2.527069925439e+06
   27 KSP unpreconditioned resid norm 1.244018617747e-09 true resid norm 1.689038425104e+01 ||r(i)||/||b|| 2.527069638823e+06
   28 KSP unpreconditioned resid norm 1.221845440869e-09 true resid norm 1.689038247082e+01 ||r(i)||/||b|| 2.527069372475e+06
   29 KSP unpreconditioned resid norm 1.200817143075e-09 true resid norm 1.689038081224e+01 ||r(i)||/||b|| 2.527069124324e+06
   30 KSP unpreconditioned resid norm 1.689037926319e+01 true resid norm 1.689037926284e+01 ||r(i)||/||b|| 2.527068892510e+06
  Linear solve did not converge due to DIVERGED_DTOL iterations 30
      Line search: gnorm after quadratic fit 7.488119005813e-01
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488429463787e-01 lambda=1.0000000000000002e-02
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488460509579e-01 lambda=1.0000000000000002e-03
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463614158e-01 lambda=1.0000000000000003e-04
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463924616e-01 lambda=1.0000000000000004e-05
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463955662e-01 lambda=1.0000000000000004e-06
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463958767e-01 lambda=1.0000000000000005e-07
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463959077e-01 lambda=1.0000000000000005e-08
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463959108e-01 lambda=1.0000000000000005e-09
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463959111e-01 lambda=1.0000000000000006e-10
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463959112e-01 lambda=1.0000000000000006e-11
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463959112e-01 lambda=1.0000000000000006e-12
      Line search: Cubic step no good, shrinking lambda, current gnorm 7.488463959112e-01 lambda=1.0000000000000007e-13
      Line search: unable to find good step length! After 12 tries
      Line search: fnorm=6.6837826673030203e-06, gnorm=7.4884639591115720e-01, ynorm=8.4202052367694780e-06, minlambda=9.9999999999999998e-13, lambda=1.0000000000000007e-13, initial slope=-6.8044408242173388e-06
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 4

With an SMP preconditioner, we get a negative Jacobian after one non-linear iteration.

On the 4th non-linear residual:
b (residual):
  0.
  0.
  0.
  4.20442e-07
  -2.73709e-06
  -1.19884e-06
  -5.49406e-07
  -9.9476e-13
  -4.70745e-06
  -1.51275e-06
  -8.1031e-07
  -9.9476e-13
  0.
  0.
  0.
  7.13222e-07
  0.
  0.
  0.
  1.4143e-07
  -1.59327e-06
  -7.22038e-07
  -1.48509e-07
  -9.9476e-13
  -2.12521e-06
  -9.01665e-07
  -7.16287e-07
  -9.9476e-13
  0.
  0.
  0.
  3.11961e-07
  0.
  0.
  0.
  -3.76809e-07
  2.4275e-07
  1.11118e-07
  1.58137e-07
  -9.9476e-13
  5.18644e-08
  2.88308e-08
  -3.88204e-09
  -9.9476e-13
  0.
  0.
  0.
  -1.89683e-08
  0.
  0.
  0.
  -2.08802e-07
  3.32957e-07
  1.76702e-07
  8.80909e-08
  -9.9476e-13
  1.77794e-07
  1.70435e-07
  3.83617e-08
  -9.9476e-13
  0.
  0.
  0.
  3.34843e-07

1st perturbed residual:
  0.
  0.
  0.
  4.39109e-07
  -2.4425e-06
  -1.25199e-06
  -5.73765e-07
  -1.0516e-12
  -4.33451e-06
  -1.57982e-06
  -8.46236e-07
  -1.0516e-12
  0.
  0.
  0.
  7.44872e-07
  0.
  0.
  0.
  1.47684e-07
  -1.45537e-06
  -7.54051e-07
  -1.55093e-07
  -1.0516e-12
  -1.92782e-06
  -9.41643e-07
  -7.48044e-07
  -1.0516e-12
  0.
  0.
  0.
  3.25652e-07
  0.
  0.
  0.
  -3.93629e-07
  0.374506
  3.25066e-07
  6.08453e-07
  -1.0516e-12
  0.373989
  2.38843e-07
  2.14939e-07
  -1.0516e-12
  0.
  0.
  0.
  -1.97861e-08
  0.
  0.
  0.
  -2.17889e-07
  0.374907
  6.18028e-07
  5.35776e-07
  -1.0516e-12
  0.374321
  6.10809e-07
  2.59249e-07
  -1.0516e-12
  0.
  0.
  0.
  3.49551e-07

Differenced and scaled by h:
  0.
  0.
  0.
  1.00421e-14
  1.58479e-13
  -2.85934e-14
  -1.31042e-14
  -3.05796e-20
  2.0063e-13
  -3.60814e-14
  -1.93267e-14
  -3.05796e-20
  0.
  0.
  0.
  1.70264e-14
  0.
  0.
  0.
  3.36458e-15
  7.41836e-14
  -1.72219e-14
  -3.54187e-15
  -3.05796e-20
  1.06187e-13
  -2.15068e-14
  -1.70843e-14
  -3.05796e-20
  0.
  0.
  0.
  7.36505e-15
  0.
  0.
  0.
  -9.04876e-15
  2.0147e-07
  1.15096e-13
  2.42253e-13
  -3.05796e-20
  2.01192e-07
  1.12979e-13
  1.17718e-13
  -3.05796e-20
  0.
  0.
  0.
  -4.39966e-16
  0.
  0.
  0.
  -4.88858e-15
  2.01686e-07
  2.37417e-13
  2.40837e-13
  -3.05796e-20
  2.01371e-07
  2.36905e-13
  1.18829e-13
  -3.05796e-20
  0.
  0.
  0.
  7.91233e-15
 This looks fine

2nd perturbed residual:
  0.
  0.
  0.
  4.39109e-07
  -2.44251e-06
  -1.25199e-06
  -5.73765e-07
  -1.0516e-12
  -4.33453e-06
  -1.57982e-06
  -8.46236e-07
  -1.0516e-12
  0.
  0.
  0.
  7.44872e-07
  0.
  0.
  0.
  1.47684e-07
  -1.45538e-06
  -7.54051e-07
  -1.55093e-07
  -1.0516e-12
  -1.92783e-06
  -9.41643e-07
  -7.48044e-07
  -1.0516e-12
  0.
  0.
  0.
  3.25652e-07
  0.
  0.
  0.
  -3.93629e-07
  0.374506
  3.25066e-07
  6.08453e-07
  -1.0516e-12
  0.373989
  2.38843e-07
  2.14939e-07
  -1.0516e-12
  0.
  0.
  0.
  -1.97861e-08
  0.
  0.
  0.
  -2.17889e-07
  0.374907
  6.18028e-07
  5.35776e-07
  -1.0516e-12
  0.374321
  6.10809e-07
  2.59249e-07
  -1.0516e-12
  0.
  0.
  0.
  3.49551e-07

difference and scaled by h:
  0.
  0.
  0.
  4.21027e-07
  6.64402e-06
  -1.19881e-06
  -5.49406e-07
  -1.28208e-12
  8.41112e-06
  -1.51275e-06
  -8.10294e-07
  -1.28208e-12
  0.
  0.
  0.
  7.13851e-07
  0.
  0.
  0.
  1.41064e-07
  3.11004e-06
  -7.22048e-07
  -1.48497e-07
  -1.28208e-12
  4.45175e-06
  -9.01696e-07
  -7.1628e-07
  -1.28208e-12
  0.
  0.
  0.
  3.08789e-07
  0.
  0.
  0.
  -3.79379e-07
  8.44684
  4.82552e-06
  1.01567e-05
  -1.28208e-12
  8.43519
  4.73676e-06
  4.93544e-06
  -1.28208e-12
  0.
  0.
  0.
  -1.84461e-08
  0.
  0.
  0.
  -2.04959e-07
  8.4559
  9.95398e-06
  1.00974e-05
  -1.28208e-12
  8.44269
  9.93249e-06
  4.98204e-06
  -1.28208e-12
  0.
  0.
  0.
  3.31734e-07

Oh, whaddya know, a different contact state it looks like.

Solve with PJFNK, FDP, bjacobi, bt, contact and mesh updated on every ComputeFunction call:
 0 Nonlinear |R| = 1.602900e+02
    0 KSP unpreconditioned resid norm 1.602900149382e+02 true resid norm 1.602900149382e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.599671499368e+02 true resid norm 1.599671499378e+02 ||r(i)||/||b|| 9.979857447728e-01
    2 KSP unpreconditioned resid norm 1.599608879658e+02 true resid norm 1.610107178943e+02 ||r(i)||/||b|| 1.004496243614e+00
    3 KSP unpreconditioned resid norm 3.181879201229e+00 true resid norm 2.229717380990e+03 ||r(i)||/||b|| 1.391051951583e+01
    4 KSP unpreconditioned resid norm 1.750714636783e-02 true resid norm 2.223510353009e+03 ||r(i)||/||b|| 1.387179578133e+01
    5 KSP unpreconditioned resid norm 4.603842515535e-04 true resid norm 2.223496256476e+03 ||r(i)||/||b|| 1.387170783740e+01
  Linear solve converged due to CONVERGED_RTOL iterations 5
      Line search: Using full step: fnorm 1.602900149382e+02 gnorm 2.848229364226e-02
 1 Nonlinear |R| = 2.848229e-02
    0 KSP unpreconditioned resid norm 2.848229364226e-02 true resid norm 2.848229364226e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.747173824518e-02 true resid norm 2.747173824123e-02 ||r(i)||/||b|| 9.645198728120e-01
    2 KSP unpreconditioned resid norm 2.388908220279e-05 true resid norm 3.371861780260e-01 ||r(i)||/||b|| 1.183844890657e+01
    3 KSP unpreconditioned resid norm 2.109140044462e-10 true resid norm 3.372577125111e-01 ||r(i)||/||b|| 1.184096044887e+01
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 2.848229364226e-02 gnorm 6.717873036118e-03
 2 Nonlinear |R| = 6.717873e-03
    0 KSP unpreconditioned resid norm 6.717873036118e-03 true resid norm 6.717873036118e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.133871915899e-03 true resid norm 2.133879423990e-03 ||r(i)||/||b|| 3.176421186463e-01
    2 KSP unpreconditioned resid norm 1.406269506134e-04 true resid norm 1.569573768382e-04 ||r(i)||/||b|| 2.336414755003e-02
    3 KSP unpreconditioned resid norm 3.330168642524e-09 true resid norm 5.392536312154e-05 ||r(i)||/||b|| 8.027148300007e-03
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 6.717873036118e-03 gnorm 2.945500798521e-05
 3 Nonlinear |R| = 2.945501e-05
    0 KSP unpreconditioned resid norm 2.945500798521e-05 true resid norm 2.945500798521e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.330423370058e-06 true resid norm 5.330437707762e-06 ||r(i)||/||b|| 1.809688087825e-01
    2 KSP unpreconditioned resid norm 1.373665285484e-06 true resid norm 1.369016702223e-06 ||r(i)||/||b|| 4.647823225546e-02
    3 KSP unpreconditioned resid norm 4.536520157338e-10 true resid norm 1.701702143766e-07 ||r(i)||/||b|| 5.777293099428e-03
    4 KSP unpreconditioned resid norm 6.231283464025e-14 true resid norm 1.701657910816e-07 ||r(i)||/||b|| 5.777142928190e-03
  Linear solve converged due to CONVERGED_RTOL iterations 4
      Line search: Using full step: fnorm 2.945500798521e-05 gnorm 4.401605929508e-07
 4 Nonlinear |R| = 4.401606e-07
    0 KSP unpreconditioned resid norm 4.401605929508e-07 true resid norm 4.401605929508e-07 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.936054969289e-08 true resid norm 2.936055024530e-08 ||r(i)||/||b|| 6.670417732870e-02
    2 KSP unpreconditioned resid norm 2.720237618859e-11 true resid norm 2.805685436947e-11 ||r(i)||/||b|| 6.374231318932e-05
    3 KSP unpreconditioned resid norm 3.472597388030e-16 true resid norm 6.793963493759e-12 ||r(i)||/||b|| 1.543519252420e-05
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 4.401605929508e-07 gnorm 7.008175587175e-11
 5 Nonlinear |R| = 7.008176e-11

Solve with PJFNK, FDP, bjacobi, bt, mesh updated on every compute function call and contact updated on:
if (snes->nfuncs == 0 || snes->ksp->reason != KSP_CONVERGED_ITERATING)

 0 Nonlinear |R| = 1.602900e+02
    0 KSP unpreconditioned resid norm 1.602900149382e+02 true resid norm 1.602900149382e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.460704081733e-02 true resid norm 2.460704072104e-02 ||r(i)||/||b|| 1.535157428897e-04
    2 KSP unpreconditioned resid norm 4.573493857903e-04 true resid norm 2.223477574342e+03 ||r(i)||/||b|| 1.387159128533e+01
  Linear solve converged due to CONVERGED_RTOL iterations 2
      Line search: Using full step: fnorm 1.602900149382e+02 gnorm 4.749061917961e-04
 1 Nonlinear |R| = 4.749062e-04
    0 KSP unpreconditioned resid norm 4.749061917961e-04 true resid norm 4.749061917961e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.617651332863e-04 true resid norm 4.617651332860e-04 ||r(i)||/||b|| 9.723291489202e-01
    2 KSP unpreconditioned resid norm 2.274102981444e-07 true resid norm 6.269124971280e-03 ||r(i)||/||b|| 1.320076486594e+01
    3 KSP unpreconditioned resid norm 1.660709058726e-11 true resid norm 6.268298672453e-03 ||r(i)||/||b|| 1.319902494585e+01
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: gnorm after quadratic fit 3.169237301824e-04
      Line search: Quadratically determined step, lambda=4.8874447856545611e-01
 2 Nonlinear |R| = 3.169237e-04
    0 KSP unpreconditioned resid norm 3.169237301824e-04 true resid norm 3.169237301824e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.140232101163e-04 true resid norm 3.139914185823e-04 ||r(i)||/||b|| 9.907475795568e-01
    2 KSP unpreconditioned resid norm 3.627751562608e-07 true resid norm 4.044078427143e-07 ||r(i)||/||b|| 1.276041533657e-03
    3 KSP unpreconditioned resid norm 1.885770844844e-11 true resid norm 1.790817502879e-07 ||r(i)||/||b|| 5.650626104421e-04
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 3.169237301824e-04 gnorm 5.467802405183e-07
 3 Nonlinear |R| = 5.467802e-07
    0 KSP unpreconditioned resid norm 5.467802405183e-07 true resid norm 5.467802405183e-07 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.816898341196e-08 true resid norm 2.816898341098e-08 ||r(i)||/||b|| 5.151792497892e-02
    2 KSP unpreconditioned resid norm 5.271431571858e-10 true resid norm 5.272333439791e-10 ||r(i)||/||b|| 9.642509090660e-04
    3 KSP unpreconditioned resid norm 1.301644928614e-13 true resid norm 5.014438606918e-12 ||r(i)||/||b|| 9.170848240904e-06
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 5.467802405183e-07 gnorm 3.042719445816e-10
 4 Nonlinear |R| = 3.042719e-10
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 4

Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 1.602775e+02
    0 KSP unpreconditioned resid norm 1.602775370798e+02 true resid norm 1.602775370798e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.580389240284e-05 true resid norm 9.580389240284e-05 ||r(i)||/||b|| 5.977374880371e-07
  Linear solve converged due to CONVERGED_RTOL iterations 1
      Line search: Using full step: fnorm 1.602775370798e+02 gnorm 1.324754493598e-04
 1 Nonlinear |R| = 1.324754e-04
    0 KSP unpreconditioned resid norm 1.324754493598e-04 true resid norm 1.324754493598e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.324081791880e-04 true resid norm 1.326770793026e-04 ||r(i)||/||b|| 1.001522017429e+00
    2 KSP unpreconditioned resid norm 1.297541569330e-04 true resid norm 1.831064988986e+01 ||r(i)||/||b|| 1.382191944119e+05
    3 KSP unpreconditioned resid norm 1.272597144252e-04 true resid norm 3.520594167363e+01 ||r(i)||/||b|| 2.657544612512e+05
    4 KSP unpreconditioned resid norm 1.249045333166e-04 true resid norm 5.085831127783e+01 ||r(i)||/||b|| 3.839074449161e+05
    5 KSP unpreconditioned resid norm 1.226754528939e-04 true resid norm 6.540330546020e+01 ||r(i)||/||b|| 4.937013293880e+05
    6 KSP unpreconditioned resid norm 1.205616076718e-04 true resid norm 7.895442616988e+01 ||r(i)||/||b|| 5.959928919015e+05
    7 KSP unpreconditioned resid norm 1.185533992703e-04 true resid norm 9.161017640998e+01 ||r(i)||/||b|| 6.915256891198e+05
    8 KSP unpreconditioned resid norm 1.166423125469e-04 true resid norm 1.034564587350e+02 ||r(i)||/||b|| 7.809481623574e+05
    9 KSP unpreconditioned resid norm 1.148207633217e-04 true resid norm 1.145685274071e+02 ||r(i)||/||b|| 8.648283735646e+05
   10 KSP unpreconditioned resid norm 1.130819714735e-04 true resid norm 1.250125888367e+02 ||r(i)||/||b|| 9.436660863641e+05
   11 KSP unpreconditioned resid norm 1.114198545580e-04 true resid norm 1.348471216171e+02 ||r(i)||/||b|| 1.017902730421e+06
   12 KSP unpreconditioned resid norm 1.098289381379e-04 true resid norm 1.441239721065e+02 ||r(i)||/||b|| 1.087929671521e+06
   13 KSP unpreconditioned resid norm 1.083042798084e-04 true resid norm 1.528892689132e+02 ||r(i)||/||b|| 1.154095114620e+06
   14 KSP unpreconditioned resid norm 1.068414045132e-04 true resid norm 1.611841893416e+02 ||r(i)||/||b|| 1.216709889421e+06
   15 KSP unpreconditioned resid norm 1.054362492186e-04 true resid norm 1.690456074068e+02 ||r(i)||/||b|| 1.276052341953e+06
   16 KSP unpreconditioned resid norm 1.040851153891e-04 true resid norm 1.765066405153e+02 ||r(i)||/||b|| 1.332372461224e+06
   17 KSP unpreconditioned resid norm 1.027846279945e-04 true resid norm 1.835971167480e+02 ||r(i)||/||b|| 1.385895406547e+06
   18 KSP unpreconditioned resid norm 1.015317000147e-04 true resid norm 1.903439732405e+02 ||r(i)||/||b|| 1.436824514734e+06
   19 KSP unpreconditioned resid norm 1.003235015919e-04 true resid norm 1.967715976818e+02 ||r(i)||/||b|| 1.485343877924e+06
   20 KSP unpreconditioned resid norm 9.915743312632e-05 true resid norm 2.029021235844e+02 ||r(i)||/||b|| 1.531620572453e+06
   21 KSP unpreconditioned resid norm 9.803110173477e-05 true resid norm 2.087556844365e+02 ||r(i)||/||b|| 1.575806577335e+06
   22 KSP unpreconditioned resid norm 9.694230058530e-05 true resid norm 2.143506346257e+02 ||r(i)||/||b|| 1.618040441920e+06
   23 KSP unpreconditioned resid norm 9.588899070209e-05 true resid norm 2.197037419773e+02 ||r(i)||/||b|| 1.658448739287e+06
   24 KSP unpreconditioned resid norm 9.486928489902e-05 true resid norm 2.248303553470e+02 ||r(i)||/||b|| 1.697147331324e+06
   25 KSP unpreconditioned resid norm 9.388143355353e-05 true resid norm 2.297445515297e+02 ||r(i)||/||b|| 1.734242477681e+06
   26 KSP unpreconditioned resid norm 9.292381197705e-05 true resid norm 2.344592647515e+02 ||r(i)||/||b|| 1.769831813250e+06
   27 KSP unpreconditioned resid norm 9.199490917412e-05 true resid norm 2.389863996808e+02 ||r(i)||/||b|| 1.804005201234e+06
   28 KSP unpreconditioned resid norm 9.109331781356e-05 true resid norm 2.433369325496e+02 ||r(i)||/||b|| 1.836845496471e+06
   29 KSP unpreconditioned resid norm 9.021772526004e-05 true resid norm 2.475209993138e+02 ||r(i)||/||b|| 1.868429210922e+06
   30 KSP unpreconditioned resid norm 2.515479748848e+02 true resid norm 2.515479748848e+02 ||r(i)||/||b|| 1.898827111744e+06
  Linear solve did not converge due to DIVERGED_DTOL iterations 30
      Line search: gnorm after quadratic fit 5.069995073375e-01
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069749809461e-01 lambda=1.0000000000000002e-02
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069725283216e-01 lambda=1.0000000000000002e-03
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722830593e-01 lambda=1.0000000000000003e-04
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722585331e-01 lambda=1.0000000000000004e-05
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722560804e-01 lambda=1.0000000000000004e-06
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558352e-01 lambda=1.0000000000000005e-07
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558106e-01 lambda=1.0000000000000005e-08
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558082e-01 lambda=1.0000000000000005e-09
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558079e-01 lambda=1.0000000000000006e-10
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558079e-01 lambda=1.0000000000000006e-11
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558079e-01 lambda=1.0000000000000006e-12
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.069722558079e-01 lambda=1.0000000000000007e-13
      Line search: unable to find good step length! After 12 tries
      Line search: fnorm=1.3247544935978778e-04, gnorm=5.0697225580792205e-01, ynorm=2.1610622831049657e-04, minlambda=9.9999999999999998e-13, lambda=1.0000000000000007e-13, initial slope=-1.0618809050736214e-03
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1

Now without the P in PJFNK, (e.g. removing the calls to updateContactSet and updateMesh from evaluation of the preconditioner):

 0 Nonlinear |R| = 1.602900e+02
    0 KSP unpreconditioned resid norm 1.602900149382e+02 true resid norm 1.602900149382e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.432659069023e+02 true resid norm 1.432659069023e+02 ||r(i)||/||b|| 8.937918369869e-01
    2 KSP unpreconditioned resid norm 1.350717791678e+02 true resid norm 1.350717816667e+02 ||r(i)||/||b|| 8.426712151645e-01
    3 KSP unpreconditioned resid norm 2.481742604667e+01 true resid norm 2.481742530837e+01 ||r(i)||/||b|| 1.548282674871e-01
    4 KSP unpreconditioned resid norm 1.010570069660e-01 true resid norm 1.010561048509e-01 ||r(i)||/||b|| 6.304578915280e-04
    5 KSP unpreconditioned resid norm 1.168146195546e-02 true resid norm 1.168120646278e-02 ||r(i)||/||b|| 7.287544684103e-05
    6 KSP unpreconditioned resid norm 1.150471172040e-02 true resid norm 1.150439382258e-02 ||r(i)||/||b|| 7.177236727451e-05
    7 KSP unpreconditioned resid norm 7.032277403283e-05 true resid norm 2.223475802277e+03 ||r(i)||/||b|| 1.387158022996e+01
  Linear solve converged due to CONVERGED_RTOL iterations 7
      Line search: Using full step: fnorm 1.602900149382e+02 gnorm 1.451539676472e-04
 1 Nonlinear |R| = 1.451540e-04
    0 KSP unpreconditioned resid norm 1.451539676472e-04 true resid norm 1.451539676472e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.436997115652e-04 true resid norm 1.436997115652e-04 ||r(i)||/||b|| 9.899812860400e-01
    2 KSP unpreconditioned resid norm 6.270235252351e-05 true resid norm 6.270235346347e-05 ||r(i)||/||b|| 4.319713369177e-01
    3 KSP unpreconditioned resid norm 6.061528069963e-05 true resid norm 6.061528263372e-05 ||r(i)||/||b|| 4.175930125523e-01
    4 KSP unpreconditioned resid norm 2.275531832661e-05 true resid norm 2.275532093491e-05 ||r(i)||/||b|| 1.567667856673e-01
    5 KSP unpreconditioned resid norm 1.294015142481e-05 true resid norm 1.294014345528e-05 ||r(i)||/||b|| 8.914770753444e-02
    6 KSP unpreconditioned resid norm 8.356953526925e-06 true resid norm 8.356941407308e-06 ||r(i)||/||b|| 5.757294507870e-02
    7 KSP unpreconditioned resid norm 5.534397472664e-07 true resid norm 5.534364807846e-07 ||r(i)||/||b|| 3.812754757967e-03
    8 KSP unpreconditioned resid norm 1.519307654842e-07 true resid norm 1.519339209827e-07 ||r(i)||/||b|| 1.046708701425e-03
    9 KSP unpreconditioned resid norm 9.089211828186e-08 true resid norm 9.089710130385e-08 ||r(i)||/||b|| 6.262116205102e-04
   10 KSP unpreconditioned resid norm 8.241968203378e-08 true resid norm 8.242204566905e-08 ||r(i)||/||b|| 5.678249585941e-04
   11 KSP unpreconditioned resid norm 4.800756870898e-08 true resid norm 4.799079869050e-08 ||r(i)||/||b|| 3.306199580238e-04
   12 KSP unpreconditioned resid norm 3.120668708371e-08 true resid norm 3.121308010534e-08 ||r(i)||/||b|| 2.150342881512e-04
   13 KSP unpreconditioned resid norm 1.098696203325e-08 true resid norm 1.098307244510e-08 ||r(i)||/||b|| 7.566498266027e-05
   14 KSP unpreconditioned resid norm 3.595916012503e-09 true resid norm 3.598525839150e-09 ||r(i)||/||b|| 2.479109525891e-05
   15 KSP unpreconditioned resid norm 9.416274057666e-10 true resid norm 1.601892025902e-03 ||r(i)||/||b|| 1.103581288109e+01
  Linear solve converged due to CONVERGED_RTOL iterations 15
      Line search: Using full step: fnorm 1.451539676472e-04 gnorm 2.927073831720e-06
 2 Nonlinear |R| = 2.927074e-06
    0 KSP unpreconditioned resid norm 2.927073831720e-06 true resid norm 2.927073831720e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.118848425608e-06 true resid norm 2.118848425608e-06 ||r(i)||/||b|| 7.238793919873e-01
    2 KSP unpreconditioned resid norm 5.693785885366e-07 true resid norm 5.693785771458e-07 ||r(i)||/||b|| 1.945214264757e-01
    3 KSP unpreconditioned resid norm 3.770601627407e-07 true resid norm 3.770601688226e-07 ||r(i)||/||b|| 1.288181270785e-01
    4 KSP unpreconditioned resid norm 2.456806842833e-07 true resid norm 2.456807409286e-07 ||r(i)||/||b|| 8.393390636963e-02
    5 KSP unpreconditioned resid norm 1.051857137614e-07 true resid norm 1.051858022500e-07 ||r(i)||/||b|| 3.593547969650e-02
    6 KSP unpreconditioned resid norm 1.041036484092e-07 true resid norm 1.041037186232e-07 ||r(i)||/||b|| 3.556579868094e-02
    7 KSP unpreconditioned resid norm 3.304649232858e-08 true resid norm 3.304672434933e-08 ||r(i)||/||b|| 1.129002076791e-02
    8 KSP unpreconditioned resid norm 2.266309888039e-08 true resid norm 2.266335204796e-08 ||r(i)||/||b|| 7.742664979052e-03
    9 KSP unpreconditioned resid norm 4.468913468206e-09 true resid norm 4.469055663584e-09 ||r(i)||/||b|| 1.526799773601e-03
   10 KSP unpreconditioned resid norm 1.667607266347e-09 true resid norm 1.667783783290e-09 ||r(i)||/||b|| 5.697785157368e-04
   11 KSP unpreconditioned resid norm 9.448521149244e-10 true resid norm 9.448837183297e-10 ||r(i)||/||b|| 3.228082968356e-04
   12 KSP unpreconditioned resid norm 5.678332222152e-10 true resid norm 5.682936554131e-10 ||r(i)||/||b|| 1.941507758549e-04
   13 KSP unpreconditioned resid norm 3.087508763292e-10 true resid norm 3.089761356878e-10 ||r(i)||/||b|| 1.055580260189e-04
   14 KSP unpreconditioned resid norm 1.602551439193e-10 true resid norm 1.604274410806e-10 ||r(i)||/||b|| 5.480812931403e-05
   15 KSP unpreconditioned resid norm 8.016393595624e-11 true resid norm 8.044215268461e-11 ||r(i)||/||b|| 2.748210578526e-05
   16 KSP unpreconditioned resid norm 4.520483409059e-11 true resid norm 4.560297199625e-11 ||r(i)||/||b|| 1.557971360410e-05
   17 KSP unpreconditioned resid norm 4.032688992341e-11 true resid norm 4.070345861596e-11 ||r(i)||/||b|| 1.390585306556e-05
   18 KSP unpreconditioned resid norm 1.649767825400e-11 true resid norm 3.079588765112e-10 ||r(i)||/||b|| 1.052104915065e-04
  Linear solve converged due to CONVERGED_RTOL iterations 18
      Line search: Using full step: fnorm 2.927073831720e-06 gnorm 7.545248011953e-11
 3 Nonlinear |R| = 7.545248e-11
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 3
 Solve Converged!

Outlier Variable Residual Norms:
  disp_x: 7.366649e-11

Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 1.602775e+02
    0 KSP unpreconditioned resid norm 1.602775370771e+02 true resid norm 1.602775370771e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.358824389986e+02 true resid norm 1.358824389986e+02 ||r(i)||/||b|| 8.477946534280e-01
    2 KSP unpreconditioned resid norm 2.480351611143e+01 true resid norm 2.480351605196e+01 ||r(i)||/||b|| 1.547535387946e-01
    3 KSP unpreconditioned resid norm 2.417605029291e+01 true resid norm 2.417605111222e+01 ||r(i)||/||b|| 1.508386736726e-01
    4 KSP unpreconditioned resid norm 1.015426056541e-01 true resid norm 1.015424120891e-01 ||r(i)||/||b|| 6.335411308463e-04
    5 KSP unpreconditioned resid norm 1.189753655170e-02 true resid norm 1.189761842347e-02 ||r(i)||/||b|| 7.423135294217e-05
    6 KSP unpreconditioned resid norm 1.358934838564e-04 true resid norm 2.231244551229e-03 ||r(i)||/||b|| 1.392113075806e-05
  Linear solve converged due to CONVERGED_RTOL iterations 6
      Line search: Using full step: fnorm 1.602775370771e+02 gnorm 1.916425318899e-04
 1 Nonlinear |R| = 1.916425e-04
    0 KSP unpreconditioned resid norm 1.916425318899e-04 true resid norm 1.916425318899e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.914394735040e-04 true resid norm 1.914394735040e-04 ||r(i)||/||b|| 9.989404315216e-01
    2 KSP unpreconditioned resid norm 1.314667684415e-04 true resid norm 1.314667667304e-04 ||r(i)||/||b|| 6.859999470576e-01
    3 KSP unpreconditioned resid norm 7.600486096655e-05 true resid norm 7.600487679062e-05 ||r(i)||/||b|| 3.965971229930e-01
    4 KSP unpreconditioned resid norm 5.427240213506e-05 true resid norm 5.427240212538e-05 ||r(i)||/||b|| 2.831960191205e-01
    5 KSP unpreconditioned resid norm 4.458204765489e-05 true resid norm 4.458205515014e-05 ||r(i)||/||b|| 2.326313199397e-01
    6 KSP unpreconditioned resid norm 1.641319431135e-05 true resid norm 1.641326556128e-05 ||r(i)||/||b|| 8.564521350985e-02
    7 KSP unpreconditioned resid norm 9.106064828496e-06 true resid norm 9.106088209795e-06 ||r(i)||/||b|| 4.751600868552e-02
    8 KSP unpreconditioned resid norm 4.122412916563e-06 true resid norm 4.122404555388e-06 ||r(i)||/||b|| 2.151090634597e-02
    9 KSP unpreconditioned resid norm 1.996874667020e-06 true resid norm 1.996865809031e-06 ||r(i)||/||b|| 1.041974236793e-02
   10 KSP unpreconditioned resid norm 8.544809322928e-07 true resid norm 8.544957019823e-07 ||r(i)||/||b|| 4.458799899770e-03
   11 KSP unpreconditioned resid norm 5.279300791751e-07 true resid norm 5.279367383073e-07 ||r(i)||/||b|| 2.754799433618e-03
   12 KSP unpreconditioned resid norm 1.190232449944e-07 true resid norm 1.189726958033e-07 ||r(i)||/||b|| 6.208052806966e-04
   13 KSP unpreconditioned resid norm 1.954185811965e-08 true resid norm 1.955040222974e-08 ||r(i)||/||b|| 1.020149443704e-04
   14 KSP unpreconditioned resid norm 1.305865018090e-08 true resid norm 1.306007174271e-08 ||r(i)||/||b|| 6.814808599074e-05
   15 KSP unpreconditioned resid norm 3.063329560261e-09 true resid norm 3.061003785955e-09 ||r(i)||/||b|| 1.597246579749e-05
   16 KSP unpreconditioned resid norm 3.757687762177e-10 true resid norm 1.405491286854e+03 ||r(i)||/||b|| 7.333921509977e+06
  Linear solve converged due to CONVERGED_RTOL iterations 16
      Line search: gnorm after quadratic fit 4.750906988607e-01
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.707949438106e-01 lambda=1.0000000000000002e-02
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703724421707e-01 lambda=1.0000000000000002e-03
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703302634560e-01 lambda=1.0000000000000003e-04
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703260462997e-01 lambda=1.0000000000000004e-05
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703256245912e-01 lambda=1.0000000000000004e-06
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255824205e-01 lambda=1.0000000000000005e-07
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255782034e-01 lambda=1.0000000000000005e-08
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255777817e-01 lambda=1.0000000000000005e-09
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255777395e-01 lambda=1.0000000000000006e-10
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255777353e-01 lambda=1.0000000000000006e-11
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255777349e-01 lambda=1.0000000000000006e-12
      Line search: Cubic step no good, shrinking lambda, current gnorm 4.703255777348e-01 lambda=1.0000000000000007e-13
      Line search: unable to find good step length! After 12 tries
      Line search: fnorm=1.9164253188987080e-04, gnorm=4.7032557773482497e-01, ynorm=1.3013523173137343e-03, minlambda=9.9999999999999998e-13, lambda=1.0000000000000007e-13, initial slope=-3.6726457059243798e-08
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1

This definitely looks to be the best yet.

Now with no line search:
 0 Nonlinear |R| = 1.602900e+02
    0 KSP unpreconditioned resid norm 1.602900149382e+02 true resid norm 1.602900149382e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.432659069023e+02 true resid norm 1.432659069023e+02 ||r(i)||/||b|| 8.937918369869e-01
    2 KSP unpreconditioned resid norm 1.350717791678e+02 true resid norm 1.350717816667e+02 ||r(i)||/||b|| 8.426712151645e-01
    3 KSP unpreconditioned resid norm 2.481742604667e+01 true resid norm 2.481742530837e+01 ||r(i)||/||b|| 1.548282674871e-01
    4 KSP unpreconditioned resid norm 1.010570069660e-01 true resid norm 1.010561048509e-01 ||r(i)||/||b|| 6.304578915280e-04
    5 KSP unpreconditioned resid norm 1.168146195546e-02 true resid norm 1.168120646278e-02 ||r(i)||/||b|| 7.287544684103e-05
    6 KSP unpreconditioned resid norm 1.150471172040e-02 true resid norm 1.150439382258e-02 ||r(i)||/||b|| 7.177236727451e-05
    7 KSP unpreconditioned resid norm 7.032277403283e-05 true resid norm 2.223475802277e+03 ||r(i)||/||b|| 1.387158022996e+01
  Linear solve converged due to CONVERGED_RTOL iterations 7
 1 Nonlinear |R| = 1.451540e-04
    0 KSP unpreconditioned resid norm 1.451539676472e-04 true resid norm 1.451539676472e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.436997115652e-04 true resid norm 1.436997115652e-04 ||r(i)||/||b|| 9.899812860401e-01
    2 KSP unpreconditioned resid norm 6.270235252351e-05 true resid norm 6.270235346241e-05 ||r(i)||/||b|| 4.319713369104e-01
    3 KSP unpreconditioned resid norm 6.061528069918e-05 true resid norm 6.061528264470e-05 ||r(i)||/||b|| 4.175930126280e-01
    4 KSP unpreconditioned resid norm 2.275531832558e-05 true resid norm 2.275532094223e-05 ||r(i)||/||b|| 1.567667857178e-01
    5 KSP unpreconditioned resid norm 1.294015145083e-05 true resid norm 1.294014653208e-05 ||r(i)||/||b|| 8.914772873125e-02
    6 KSP unpreconditioned resid norm 8.356953630210e-06 true resid norm 8.356948186098e-06 ||r(i)||/||b|| 5.757299177939e-02
    7 KSP unpreconditioned resid norm 5.534384347226e-07 true resid norm 5.534538961578e-07 ||r(i)||/||b|| 3.812874736590e-03
    8 KSP unpreconditioned resid norm 1.519310137152e-07 true resid norm 1.519299596478e-07 ||r(i)||/||b|| 1.046681410853e-03
    9 KSP unpreconditioned resid norm 9.089732006464e-08 true resid norm 9.089290654531e-08 ||r(i)||/||b|| 6.261827218270e-04
   10 KSP unpreconditioned resid norm 8.366728648700e-08 true resid norm 8.367901341938e-08 ||r(i)||/||b|| 5.764845065948e-04
   11 KSP unpreconditioned resid norm 4.800852290198e-08 true resid norm 4.799638783283e-08 ||r(i)||/||b|| 3.306584629466e-04
   12 KSP unpreconditioned resid norm 3.120680812598e-08 true resid norm 3.121193234027e-08 ||r(i)||/||b|| 2.150263809263e-04
   13 KSP unpreconditioned resid norm 1.098696135347e-08 true resid norm 1.098350865942e-08 ||r(i)||/||b|| 7.566798784387e-05
   14 KSP unpreconditioned resid norm 3.595917199957e-09 true resid norm 3.593440478279e-09 ||r(i)||/||b|| 2.475606100560e-05
   15 KSP unpreconditioned resid norm 9.416260286095e-10 true resid norm 1.601892024197e-03 ||r(i)||/||b|| 1.103581286934e+01
  Linear solve converged due to CONVERGED_RTOL iterations 15
 2 Nonlinear |R| = 2.927074e-06
    0 KSP unpreconditioned resid norm 2.927073815329e-06 true resid norm 2.927073815329e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.118848439229e-06 true resid norm 2.118848439229e-06 ||r(i)||/||b|| 7.238794006944e-01
    2 KSP unpreconditioned resid norm 5.693785743520e-07 true resid norm 5.693785642801e-07 ||r(i)||/||b|| 1.945214231695e-01
    3 KSP unpreconditioned resid norm 3.770601778621e-07 true resid norm 3.770601832994e-07 ||r(i)||/||b|| 1.288181327457e-01
    4 KSP unpreconditioned resid norm 2.456806629715e-07 true resid norm 2.456807188981e-07 ||r(i)||/||b|| 8.393389931318e-02
    5 KSP unpreconditioned resid norm 1.051857099077e-07 true resid norm 1.051857980821e-07 ||r(i)||/||b|| 3.593547847384e-02
    6 KSP unpreconditioned resid norm 1.041036432009e-07 true resid norm 1.041037121660e-07 ||r(i)||/||b|| 3.556579667407e-02
    7 KSP unpreconditioned resid norm 3.304649190305e-08 true resid norm 3.304672156718e-08 ||r(i)||/||b|| 1.129001988065e-02
    8 KSP unpreconditioned resid norm 2.266309599872e-08 true resid norm 2.266335056541e-08 ||r(i)||/||b|| 7.742664515917e-03
    9 KSP unpreconditioned resid norm 4.480163787678e-09 true resid norm 4.480319017262e-09 ||r(i)||/||b|| 1.530647773144e-03
   10 KSP unpreconditioned resid norm 1.669011427119e-09 true resid norm 1.669198753969e-09 ||r(i)||/||b|| 5.702619268524e-04
   11 KSP unpreconditioned resid norm 1.009174409193e-09 true resid norm 1.009246866312e-09 ||r(i)||/||b|| 3.447972036191e-04
   12 KSP unpreconditioned resid norm 6.807042441187e-10 true resid norm 6.809063040681e-10 ||r(i)||/||b|| 2.326235506950e-04
   13 KSP unpreconditioned resid norm 3.117391627444e-10 true resid norm 3.118024608062e-10 ||r(i)||/||b|| 1.065236070144e-04
   14 KSP unpreconditioned resid norm 1.602670897926e-10 true resid norm 1.604313093561e-10 ||r(i)||/||b|| 5.480945117130e-05
   15 KSP unpreconditioned resid norm 8.016556482573e-11 true resid norm 8.045376559510e-11 ||r(i)||/||b|| 2.748607335209e-05
   16 KSP unpreconditioned resid norm 4.519620882604e-11 true resid norm 4.556534750297e-11 ||r(i)||/||b|| 1.556685973013e-05
   17 KSP unpreconditioned resid norm 4.032472560574e-11 true resid norm 4.070931248260e-11 ||r(i)||/||b|| 1.390785304744e-05
   18 KSP unpreconditioned resid norm 1.649102536810e-11 true resid norm 3.079635551631e-10 ||r(i)||/||b|| 1.052120905016e-04
  Linear solve converged due to CONVERGED_RTOL iterations 18
 3 Nonlinear |R| = 7.545244e-11
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 3
 Solve Converged!

Outlier Variable Residual Norms:
  disp_x: 7.366619e-11

Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 1.602775e+02
    0 KSP unpreconditioned resid norm 1.602775370771e+02 true resid norm 1.602775370771e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.358824369198e+02 true resid norm 1.358824369198e+02 ||r(i)||/||b|| 8.477946404582e-01
    2 KSP unpreconditioned resid norm 2.480351366302e+01 true resid norm 2.480351547695e+01 ||r(i)||/||b|| 1.547535352070e-01
    3 KSP unpreconditioned resid norm 2.417605432595e+01 true resid norm 2.417605453121e+01 ||r(i)||/||b|| 1.508386950043e-01
    4 KSP unpreconditioned resid norm 1.015425861477e-01 true resid norm 1.015423814253e-01 ||r(i)||/||b|| 6.335409395294e-04
    5 KSP unpreconditioned resid norm 1.189753752713e-02 true resid norm 1.189746827201e-02 ||r(i)||/||b|| 7.423041612059e-05
    6 KSP unpreconditioned resid norm 1.359628794207e-04 true resid norm 2.231495014745e-03 ||r(i)||/||b|| 1.392269344438e-05
  Linear solve converged due to CONVERGED_RTOL iterations 6
 1 Nonlinear |R| = 1.900552e-04
    0 KSP unpreconditioned resid norm 1.900552381709e-04 true resid norm 1.900552381709e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.898702884725e-04 true resid norm 1.898702884725e-04 ||r(i)||/||b|| 9.990268634523e-01
    2 KSP unpreconditioned resid norm 1.314227422923e-04 true resid norm 1.314227395516e-04 ||r(i)||/||b|| 6.914975920497e-01
    3 KSP unpreconditioned resid norm 7.704532466503e-05 true resid norm 7.704531236459e-05 ||r(i)||/||b|| 4.053837879244e-01
    4 KSP unpreconditioned resid norm 5.509323253559e-05 true resid norm 5.509326398325e-05 ||r(i)||/||b|| 2.898802711963e-01
    5 KSP unpreconditioned resid norm 4.501216544552e-05 true resid norm 4.501216182516e-05 ||r(i)||/||b|| 2.368372598323e-01
    6 KSP unpreconditioned resid norm 1.941475764040e-05 true resid norm 1.941479742364e-05 ||r(i)||/||b|| 1.021534455482e-01
    7 KSP unpreconditioned resid norm 9.419548371442e-06 true resid norm 9.419507975058e-06 ||r(i)||/||b|| 4.956194875612e-02
    8 KSP unpreconditioned resid norm 3.438814420353e-06 true resid norm 3.438896000206e-06 ||r(i)||/||b|| 1.809419215856e-02
    9 KSP unpreconditioned resid norm 2.032790145627e-06 true resid norm 2.032794823724e-06 ||r(i)||/||b|| 1.069581056164e-02
   10 KSP unpreconditioned resid norm 8.916617451537e-07 true resid norm 8.916400002516e-07 ||r(i)||/||b|| 4.691478166205e-03
   11 KSP unpreconditioned resid norm 5.482537600575e-07 true resid norm 5.482281981120e-07 ||r(i)||/||b|| 2.884572945151e-03
   12 KSP unpreconditioned resid norm 1.255813891636e-07 true resid norm 1.256315190284e-07 ||r(i)||/||b|| 6.610263428541e-04
   13 KSP unpreconditioned resid norm 1.903434285886e-08 true resid norm 1.905638194985e-08 ||r(i)||/||b|| 1.002675965853e-04
   14 KSP unpreconditioned resid norm 9.830412365161e-09 true resid norm 9.823309558905e-09 ||r(i)||/||b|| 5.168660255537e-05
   15 KSP unpreconditioned resid norm 4.191324767477e-09 true resid norm 4.212480463306e-09 ||r(i)||/||b|| 2.216450598177e-05
   16 KSP unpreconditioned resid norm 6.928401977611e-10 true resid norm 1.405659696443e+03 ||r(i)||/||b|| 7.396058693101e+06
  Linear solve converged due to CONVERGED_RTOL iterations 16
 2 Nonlinear |R| = 4.047539e-08
    0 KSP unpreconditioned resid norm 4.047539349968e-08 true resid norm 4.047539349968e-08 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.027389740640e-08 true resid norm 4.027389740635e-08 ||r(i)||/||b|| 9.950217632021e-01
    2 KSP unpreconditioned resid norm 2.330944962215e-08 true resid norm 2.330944990101e-08 ||r(i)||/||b|| 5.758918662816e-01
    3 KSP unpreconditioned resid norm 1.508484090262e-08 true resid norm 1.508483781262e-08 ||r(i)||/||b|| 3.726915666118e-01
    4 KSP unpreconditioned resid norm 1.204155654751e-08 true resid norm 1.204155452527e-08 ||r(i)||/||b|| 2.975030872859e-01
    5 KSP unpreconditioned resid norm 8.106561437173e-09 true resid norm 8.106588851792e-09 ||r(i)||/||b|| 2.002843740569e-01
    6 KSP unpreconditioned resid norm 4.133308962119e-09 true resid norm 4.133292826637e-09 ||r(i)||/||b|| 1.021186570223e-01
    7 KSP unpreconditioned resid norm 2.119336320758e-09 true resid norm 2.119329973652e-09 ||r(i)||/||b|| 5.236094798358e-02
    8 KSP unpreconditioned resid norm 4.964851436983e-10 true resid norm 4.965117184690e-10 ||r(i)||/||b|| 1.226700164071e-02
    9 KSP unpreconditioned resid norm 2.490792441481e-10 true resid norm 2.490765817510e-10 ||r(i)||/||b|| 6.153777893547e-03
   10 KSP unpreconditioned resid norm 1.157373799621e-10 true resid norm 1.157143358966e-10 ||r(i)||/||b|| 2.858881060601e-03
   11 KSP unpreconditioned resid norm 4.316193659546e-11 true resid norm 4.320288705135e-11 ||r(i)||/||b|| 1.067386461646e-03
   12 KSP unpreconditioned resid norm 3.771815501907e-11 true resid norm 3.772774667611e-11 ||r(i)||/||b|| 9.321156241854e-04
   13 KSP unpreconditioned resid norm 2.628320854599e-11 true resid norm 2.628007916162e-11 ||r(i)||/||b|| 6.492853284261e-04
   14 KSP unpreconditioned resid norm 1.834519283302e-11 true resid norm 1.834387519976e-11 ||r(i)||/||b|| 4.532105463015e-04
   15 KSP unpreconditioned resid norm 1.408027599157e-11 true resid norm 1.408495690258e-11 ||r(i)||/||b|| 3.479881400706e-04
   16 KSP unpreconditioned resid norm 6.953505432275e-12 true resid norm 6.952460932156e-12 ||r(i)||/||b|| 1.717700640072e-04
   17 KSP unpreconditioned resid norm 1.786998910535e-12 true resid norm 1.782495219068e-12 ||r(i)||/||b|| 4.403898430493e-05
   18 KSP unpreconditioned resid norm 3.275393770191e-13 true resid norm 2.812006995075e-12 ||r(i)||/||b|| 6.947448194907e-05
  Linear solve converged due to CONVERGED_RTOL iterations 18
 3 Nonlinear |R| = 2.637458e-12
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 3
 Solve Converged!

This definitely looks great. I'm curious though why the basic line search succeeded when the bt line search did not.

For the Kyle Gamble problem, evaluating on every residual works optimally. For the thermal expansion contact problem, evaluating on non-linear residuals works optimally! Wtf!!! Why do the linear solves look perfect for Kyle Gamble when evaulating on every residual while the true residual norm looks totally wrong in the thermal problem???

I need to figure out how the mesh distances are actually evaulated.

It looks to me like the displaced nodes are moved every time UpdateDisplacedMeshThread::onNode is called. `_mesh` and `_ref_mesh` are two different objects. Good. The nodes are updated in this way:

      displaced_node(direction) =
          reference_node(direction) +
          (*_nl_ghosted_soln)(reference_node.dof_number(_nonlinear_system_number, _var_nums[i], 0));

So it depends ONLY on the solution being passed in.

It's clear looking at the 1D problem that the solves diverge on the first linear iteration of the second non-linear iteration. Consequently, since the b vectors are the same in that problem, it has to be the action of the Jacobians that are fairly dramatically different. That's the next bit to zone in on in the debugger.

# 1/9/18

Ok, the residual difference appears to come only on the node 2 displacement. It may also be present on the node 3 displacement (rhs block) but the residual is fairly large.

form 1:

1.59774e-05
-3.54754e-08
8.40226e-05
3.54754e-08
37.6692
3.56562e-08
-4.39588e-07
-3.56562e-08

form 2:

1.59774e-05
-3.54754e-08
2.16891e-05
3.54754e-08
37.6692
3.56562e-08
-4.39588e-07
-3.56562e-08

So the change is that the kinematic formulation uses the residual vector in computing its contact force. Consequently, if the contact force is only evaluated on non-linear residuals vs. every residual, then the contact force is going to be different from case to case with my code changes (where the contact force is evaluated when updating the contact set). There is a question about **what** residual is being used, e.g. is it an old residual? Is it the current residual?

As I suspected, in my current code, the residual being used in `_residual_copy` is an old residual. So if we're going to use this residual formulation, it's clear that the contact force must be evaluated at the end of the residual computing chain.

Ok, it looks like having the contact set update on every residual evaluation really is a bad idea. We're finite differencing over a non-smooth point. That's a fundamental problem. The contact state really probably should only be updated on non-linear residuals. Because look what happens if we go back to the bad set-up of having contact updates on Jacobian setup instead of residual:

 0 Nonlinear |R| = 2.043146e+02
b    0 KSP unpreconditioned resid norm 2.043145722763e+02 true resid norm 2.043145722763e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.931136770970e+02 true resid norm 1.931123251790e+02 ||r(i)||/||b|| 9.451715706202e-01
    2 KSP unpreconditioned resid norm 1.040113083427e+01 true resid norm 1.055288817760e+01 ||r(i)||/||b|| 5.165019831933e-02
    3 KSP unpreconditioned resid norm 2.358473647754e-01 true resid norm 8.010636481672e-01 ||r(i)||/||b|| 3.920736730829e-03
    4 KSP unpreconditioned resid norm 3.021798412022e-03 true resid norm 6.058872614985e-01 ||r(i)||/||b|| 2.965462789796e-03
    5 KSP unpreconditioned resid norm 4.153875121281e-05 true resid norm 4.881126090654e-01 ||r(i)||/||b|| 2.389024941429e-03
 1 Nonlinear |R| = 1.239123e+02
^R
    0 KSP unpreconditioned resid norm 1.239123135569e+02 true resid norm 1.239123135569e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.999177731873e-03 true resid norm 4.999177731875e-03 ||r(i)||/||b|| 4.034447899788e-05
    2 KSP unpreconditioned resid norm 7.990795445406e-06 true resid norm 7.984625352344e-06 ||r(i)||/||b|| 6.443770698121e-08
 2 Nonlinear |R| = 9.670852e-04
    0 KSP unpreconditioned resid norm 9.670852423165e-04 true resid norm 9.670852423165e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.326058041359e-07 true resid norm 1.326058041359e-07 ||r(i)||/||b|| 1.371190442512e-04
    2 KSP unpreconditioned resid norm 4.957370676372e-10 true resid norm 4.960800126583e-10 ||r(i)||/||b|| 5.129641017684e-07
 3 Nonlinear |R| = 6.909912e-10

Boom beautiful. Need to fix.

# 1/10/18

So the reason PJFNK with -pc_type none and JFNK ran with different solve histories is because of these lines in the constraint:

```
  // Release
  if (_model != CM_GLUED && -contact_pressure >= _tension_release && pinfo->_locked_this_step < 2)
  {
    pinfo->release();
    pinfo->_contact_force.zero();
  }
  // Capture
  else
  {
    pinfo->capture();
    if (_formulation == CF_KINEMATIC || _formulation == CF_TANGENTIAL_PENALTY)
      ++pinfo->_locked_this_step;


```
The key is the locking which prevents things from pinballing around. With PJFNK there are many more calls to updateContactSet. If I delete the locking lines, the JFNK solve (and the PJFNK solve with -pc_type none) no longer converge.

This is just more motivation to take the contact update out of every linear residual...kind of like they had in the beginning. Remember that there's usually a reason for the way things are!!!

MatMFFDSetFunction expects a residual computing function with three arguments: (void*,Vec,Vec)

However, SNESSetFunction expects a function with four arguments!: (SNES,Vec,Vec,void*)

What about MatFDColoringSetFunction? It only expects one argument! (void)

SNESComputeFunction takes three arguments so it works to pass it to MatMFFDSetFunction!

frame #1: 0x00000001038be412 libmoose-dbg.0.dylib`NonlinearSystem::setupColoringFiniteDifferencedPreconditioner(this=0x0000000113808c18) at NonlinearSystem.C:357
   354
   355 	  MatFDColoringCreate(petsc_mat->mat(), iscoloring, &_fdcoloring);
   356 	  MatFDColoringSetFromOptions(_fdcoloring);
-> 357 	  MatFDColoringSetFunction(_fdcoloring,
   358 	                           (PetscErrorCode(*)(void)) & libMesh::__libmesh_petsc_snes_fd_residual,
   359 	                           &petsc_nonlinear_solver);

frame #1: 0x0000000106fcf3ee libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x0000000112ca2b90, jac_in=0x0000000112ca2ac0, x_in=0x0000000112ca2330, r_in=0x0000000112ca28e0, (null)=0.0000000001, (null)=200) at petsc_nonlinear_solver.C:727
   724 	  // Should actually be a PetscReal, but I don't know which version of PETSc first introduced PetscReal
   725 	  Real final_residual_norm=0.;
   726
-> 727 	  ierr = SNESSetFunction (_snes, r->vec(), __libmesh_petsc_snes_residual, this);

(lldb) up
frame #1: 0x000000010fae9239 libpetsc.3.08.dylib`MatCreateSNESMF(snes=0x0000000115135c60, J=0x00007fff5fbfd5e0) at snesmfj.c:160
   157 	  if (snes->npc && snes->npcside== PC_LEFT) {
   158 	    ierr = MatMFFDSetFunction(*J,(PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunctionDefaultNPC,snes);CHKERRQ(ierr);
   159 	  } else {
-> 160 	    ierr = MatMFFDSetFunction(*J,(PetscErrorCode (*)(void*,Vec,Vec))SNESComputeFunction,snes);CHKERRQ(ierr);


I believe there has to be an intermediate step that transitions between the interface and the implementation. SNESComputeFunction is an example of that intermediate layer. It is the interface. It takes three arguments and then calls the implementation with four arguments.

# 1/12/18

Changing the locking doesn't affect the outcome of the sliding block test failure with line-search = none.

# 1/15/18

Here are the failing tests on rod (9 failures):

preconditioners/pbp.check_petsc_options_test................................ FAILED (EXPECTED OUTPUT MISSING)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
postprocessors/scalar_coupled_postprocessor.test.............................................. FAILED (CRASH)
parser/parser_dependency.testpbp.............................................................. FAILED (CRASH)
preconditioners/pbp.test...................................................................... FAILED (CRASH)
preconditioners/pbp.pbp_adapt_test............................................................ FAILED (CRASH)
preconditioners/pbp.lots_of_variables......................................................... FAILED (CRASH)

So looks to pretty much just be PBP.

Here are the tests that fail on my Mac (29 failures):

preconditioners/pbp.check_petsc_options_test................................ FAILED (EXPECTED OUTPUT MISSING)
misc/exception.parallel_error_jacobian_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
vectorpostprocessors/csv_reader.tester_fail.................................. FAILED (EXPECTED ERROR MISSING)
misc/exception.parallel_error_residual_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
restart/kernel_restartable.parallel_error2................................... FAILED (EXPECTED ERROR MISSING)
mesh/ghost_functors.geometric_edge_neighbor_two_2D.......................................... FAILED (EXODIFF)
mesh/custom_partitioner.custom_linear_partitioner........................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_2D.......................................... FAILED (EXODIFF)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.......................... FAILED (EXODIFF)
ics/depend_on_uo.ic_depend_on_uo............................................................ FAILED (EXODIFF)
misc/exception.parallel_exception_residual_transient_non_zero_rank.......................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac...................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac...................................... FAILED (EXODIFF)
mesh/subdomain_partitioner.subdomain_partitioner............................................ FAILED (EXODIFF)
utils/random.test_uo_par_mesh............................................................... FAILED (EXODIFF)
mesh/centroid_partitioner.centroid_partitioner_test......................................... FAILED (EXODIFF)
utils/random.test_par_mesh.................................................................. FAILED (EXODIFF)
mesh/nemesis.nemesis_test.................................................................... FAILED (ERRMSG)
restart/restartable_types.first_parallel..................................................... FAILED (ERRMSG)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
userobjects/setup_interface_count.GeneralUserObject......................................... FAILED (CSVDIFF)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)
parser/parser_dependency.testpbp.............................................................. FAILED (CRASH)
postprocessors/scalar_coupled_postprocessor.test.............................................. FAILED (CRASH)
preconditioners/pbp.test...................................................................... FAILED (CRASH)
preconditioners/pbp.pbp_adapt_test............................................................ FAILED (CRASH)
restart/restart.test_xda_restart_part_2....................................................... FAILED (CRASH)

This looks similar to the issue I had before. A lot of mesh stuff.

Trying a git clean -fxd in moose. Still 29 failures.

Now git clean -fxd in libmesh and git clean -fxd in MOOSE. Still 29 fails

Unload ccache, and cleaning. Still 29 fails

After commenting out these lines from my .bash_profile:

# export PATH=/Applications/Emacs.app/Contents/MacOS/bin:$PATH
# export PATH=/Applications/Emacs.app/Contents/MacOS:$PATH
# export PATH=/Applications/ParaView-5.4.1.app/Contents/MacOS:$PATH
# export PATH=/Applications/Cubit-15.2/Cubit.app/Contents/MacOS:$PATH
# export PATH=$HOME/gmsh/opencc_build:$PATH
# export PATH=$HOME/bash:$PATH

Then I'm down to the correct number of failures which is 9.

It's definitely a run-time issue (which is good) because when I run with the same compiled executable with those lines uncommented, then I'm back to 29 failures.

After commenting out my bash directory from PATH, down to 22 failures:

misc/exception.parallel_exception_residual_transient_non_zero_rank.............. skipped (skipped dependency)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.............. skipped (skipped dependency)
restart/kernel_restartable.parallel_error2...................................... skipped (skipped dependency)
restart/restart.test_xda_restart_part_2......................................... skipped (skipped dependency)
restart/restartable_types.first_parallel........................................ skipped (skipped dependency)
utils/random.test_par_mesh...................................................... skipped (skipped dependency)
utils/random.test_uo_par_mesh................................................... skipped (skipped dependency)
utils/random.threads_verification_uo............................................ skipped (skipped dependency)
preconditioners/pbp.check_petsc_options_test................................ FAILED (EXPECTED OUTPUT MISSING)
misc/exception.parallel_error_residual_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
misc/exception.parallel_error_jacobian_transient_non_zero_rank............... FAILED (EXPECTED ERROR MISSING)
vectorpostprocessors/csv_reader.tester_fail.................................. FAILED (EXPECTED ERROR MISSING)
ics/depend_on_uo.ic_depend_on_uo............................................................ FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_2D.......................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_2D.......................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac...................................... FAILED (EXODIFF)
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac...................................... FAILED (EXODIFF)
mesh/custom_partitioner.custom_linear_partitioner........................................... FAILED (EXODIFF)
mesh/centroid_partitioner.centroid_partitioner_test......................................... FAILED (EXODIFF)
mesh/subdomain_partitioner.subdomain_partitioner............................................ FAILED (EXODIFF)
mesh/nemesis.nemesis_test.................................................................... FAILED (ERRMSG)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
postprocessors/all_print_pps.test........................................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
userobjects/setup_interface_count.GeneralUserObject......................................... FAILED (CSVDIFF)
parser/parser_dependency.testpbp.............................................................. FAILED (CRASH)
postprocessors/scalar_coupled_postprocessor.test.............................................. FAILED (CRASH)
preconditioners/pbp.test...................................................................... FAILED (CRASH)
preconditioners/pbp.pbp_adapt_test............................................................ FAILED (CRASH)

After commenting out the cubit directory, down to 8 failures:

misc/exception.parallel_exception_residual_transient_non_zero_rank.............. skipped (skipped dependency)
misc/exception.parallel_exception_jacobian_transient_non_zero_rank.............. skipped (skipped dependency)
restart/kernel_restartable.parallel_error2...................................... skipped (skipped dependency)
restart/restart.test_xda_restart_part_2......................................... skipped (skipped dependency)
restart/restartable_types.first_parallel........................................ skipped (skipped dependency)
utils/random.threads_verification_uo............................................ skipped (skipped dependency)
utils/random.test_par_mesh...................................................... skipped (skipped dependency)
utils/random.test_uo_par_mesh................................................... skipped (skipped dependency)
mesh/centroid_partitioner.centroid_partitioner_test.......................................... [MIN_CPUS=4] OK
mesh/subdomain_partitioner.subdomain_partitioner............................................. [MIN_CPUS=4] OK
mesh/nemesis.nemesis_test.................................................................... [MIN_CPUS=4] OK
mesh/ghost_functors.geometric_edge_neighbor_one_2D........................................... [MIN_CPUS=3] OK
mesh/ghost_functors.geometric_edge_neighbor_two_2D........................................... [MIN_CPUS=3] OK
mesh/ghost_functors.geometric_edge_neighbor_one_3D_Mac....................................... [MIN_CPUS=3] OK
mesh/ghost_functors.geometric_edge_neighbor_two_3D_Mac....................................... [MIN_CPUS=3] OK
ics/depend_on_uo.ic_depend_on_uo............................................................. [MIN_CPUS=2] OK
mesh/custom_partitioner.custom_linear_partitioner............................................ [MIN_CPUS=2] OK
misc/exception.parallel_error_jacobian_transient_non_zero_rank............................... [MIN_CPUS=2] OK
misc/exception.parallel_error_residual_transient_non_zero_rank............................... [MIN_CPUS=2] OK
postprocessors/all_print_pps.test............................................................ [MIN_CPUS=2] OK
userobjects/setup_interface_count.GeneralUserObject.......................................... [MIN_CPUS=2] OK
vectorpostprocessors/csv_reader.tester_fail.................................................. [MIN_CPUS=2] OK
preconditioners/pbp.check_petsc_options_test................................ FAILED (EXPECTED OUTPUT MISSING)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
parser/parser_dependency.testpbp.............................................................. FAILED (CRASH)
postprocessors/scalar_coupled_postprocessor.test.............................................. FAILED (CRASH)
preconditioners/pbp.test...................................................................... FAILED (CRASH)
preconditioners/pbp.pbp_adapt_test............................................................ FAILED (CRASH)

The difference from linux is that on linux, pbp.lots_of_variables fails. This is because lots_of_variables only runs when gcc is the compiler. Ok, glad I figured out what the failure was!!

This is the additional test fail that comes from having my bash directory:
restart/restart.test_xda_restart_part_2....................................................... FAILED (CRASH)

Wait there's actually more than 8 failures on Mac, there's 9 (this is after appending to my path instead of prepending):

preconditioners/pbp.check_petsc_options_test................................ FAILED (EXPECTED OUTPUT MISSING)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
parser/parser_dependency.testpbp.............................................................. FAILED (CRASH)
postprocessors/scalar_coupled_postprocessor.test.............................................. FAILED (CRASH)
preconditioners/pbp.pbp_adapt_test............................................................ FAILED (CRASH)
preconditioners/pbp.test...................................................................... FAILED (CRASH)
restart/restart.test_xda_restart_part_2....................................................... FAILED (CRASH)

So the additional one is the xda_restart_test. Nope, wrong again. This is the full litany of failing tests on Linux:

preconditioners/pbp.check_petsc_options_test................................ FAILED (EXPECTED OUTPUT MISSING)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
postprocessors/scalar_coupled_postprocessor.test.............................................. FAILED (CRASH)
parser/parser_dependency.testpbp.............................................................. FAILED (CRASH)
restart/restart.test_xda_restart_part_2....................................................... FAILED (CRASH)
preconditioners/pbp.pbp_adapt_test............................................................ FAILED (CRASH)
preconditioners/pbp.test...................................................................... FAILED (CRASH)

That's the exact list of failed tests from mac which is good. Now that lots_of_variables test isn't running on rod because the Petsc version it claims is not a release. I don't know how in the hell it ran before and not now, but not that important.

Here are failing tests on MOOSE HEAD with PETSC HEAD:

outputs/format.json........................................................................... FAILED (CRASH)
-------------------------------------------------------------------------------------------------------------
Ran 1481 tests in 262.6 seconds
1480 passed, 63 skipped, 0 pending, 1 FAILED

And then my feature with PETSC HEAD:

misc/exception.parallel_exception_residual_transient......................................... FAILED (ERRMSG)
executioners/nullspace.test_singular........................................................ FAILED (CSVDIFF)
executioners/nullspace.test_singular_contaminated........................................... FAILED (CSVDIFF)
outputs/iterative.csv....................................................................... FAILED (CSVDIFF)
outputs/format.json........................................................................... FAILED (CRASH)
-------------------------------------------------------------------------------------------------------------
Ran 1478 tests in 321.2 seconds
1473 passed, 64 skipped, 0 pending, 5 FAILED

# 1/22/18

Want to edit these commits:

792a31c2dd99ab -> Undesirably restores to templated getVariable for a coupled
varialbe in GapValueAux

faf098727e3b -> restores AuxKernel.C getVar method

c7e8e80f35be -> Adds stupid methods in FEProblemBase

# 1/23/18

Warehouses:
0x000000010d87ebc0 (feproblem warehouse?)
0x000000010d887180 (displaced problem warehouse)

All of a sudden with scalar_strain_zz we see a different address!! This is
becaue there are aux and nonlinear systems for both feproblem and displaced
problem -> 4 systems, 4 warehouses, 4 sets of variables

0x000000010d87d3a0
0x000000010d847440

# 1/24/18

dont_reinit_mat_for_aux: additional call for four qp problem increases average
material property value by four
picard_failures: counts residual evaluations
all_print_pps: counts residual evaluations
iterative: counts non-linear iterations

AFAICT the computing_initial_residual infrastructure is not used for anything.

I'm interested in erroring solves and preferably as small as possible:

generalized_plane_strain_tm_contact: 55 DOFs

It looks like all the single_pnt tests have 22 dofs:
combined/test:contact_verification/patch_tests/single_pnt_2d.frictionless_kin_sm: 22 dofs

Well it looks in fact like that call to the residual at the end of the Jacobian evaluation really does impact the b vector in the linear solve! I want to look at this in the debugger and see the 1-to-1 correlation between norms.

The result I reported to the Bison team may have been a red herring. That bad behavior in the true residual norm may have been because I had deleted that residual call at the end of the jacobian evaluation. I need to check that tomorrow as well. In Gamble's contact examples and in the single_point glued contact test, I don't see any bad linear solves.

What's interesting is that I don't see any difference between the non-linaer residual and the initial residual in Gamble's sliding block examples even when we see a jump in the value of the non-linear residual. The only place I see it in MOOSE head is the single_point glued contact test. Aww wait on a second:

Here's the single point test:
 0 Nonlinear |R| = 5.471838e+03




*** Warning ***
Warning in PenetrationLocator. Penetration is not detected for one or more slave nodes. This could be because those slave nodes simply do not project to faces on the master surface. However, this could also be because contact should be enforced on those nodes, but the faces that they project to are outside the contact patch, which will give an erroneous result. Use appropriate options for 'patch_size' and 'patch_update_strategy' in the Mesh block to avoid this issue. Setting 'patch_update_strategy=iteration' is recommended because it completely avoids this potential issue. Also note that this warning is printed only once, so a similar situation could occur multiple times during the simulation but this warning is printed only at the first occurrence.

    0 KSP unpreconditioned resid norm 5.471838102659e+03 true resid norm 5.471838102659e+03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.678501661859e+02 true resid norm 4.678501661859e+02 ||r(i)||/||b|| 8.550146356826e-02
    2 KSP unpreconditioned resid norm 1.002293429030e+02 true resid norm 1.002293653303e+02 ||r(i)||/||b|| 1.831731192514e-02
    3 KSP unpreconditioned resid norm 1.177693095297e-01 true resid norm 1.177744290695e-01 ||r(i)||/||b|| 2.152374153984e-05
 1 Nonlinear |R| = 1.001988e+02
    0 KSP unpreconditioned resid norm 1.001988279477e+02 true resid norm 1.001988279477e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.758050384494e+01 true resid norm 1.758050384494e+01 ||r(i)||/||b|| 1.754561825226e-01
    2 KSP unpreconditioned resid norm 3.936647766797e+00 true resid norm 3.936652416703e+00 ||r(i)||/||b|| 3.928840783207e-02
    3 KSP unpreconditioned resid norm 5.889641465133e-06 true resid norm 1.071253408567e-05 ||r(i)||/||b|| 1.069127683935e-07
 2 Nonlinear |R| = 5.952291e-02
    0 KSP unpreconditioned resid norm 4.853631832011e+05 true resid norm 4.853631832011e+05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.128301534340e-01 true resid norm 2.128301534416e-01 ||r(i)||/||b|| 4.384966985710e-07
 3 Nonlinear |R| = 2.128295e-01
    0 KSP unpreconditioned resid norm 2.128294741483e-01 true resid norm 2.128294741483e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.795356458619e-02 true resid norm 8.795356458619e-02 ||r(i)||/||b|| 4.132583841508e-01
    2 KSP unpreconditioned resid norm 1.467021623138e-02 true resid norm 1.469159462127e-02 ||r(i)||/||b|| 6.902988733148e-02
    3 KSP unpreconditioned resid norm 2.151601016215e-03 true resid norm 2.146514383639e-03 ||r(i)||/||b|| 1.008560676208e-02
    4 KSP unpreconditioned resid norm 1.887306724317e-05 true resid norm 7.797521920927e-05 ||r(i)||/||b|| 3.663741571571e-04
 4 Nonlinear |R| = 4.791060e-04
    0 KSP unpreconditioned resid norm 4.791059721812e-04 true resid norm 4.791059721812e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.189691087299e-06 true resid norm 6.189691087297e-06 ||r(i)||/||b|| 1.291925262196e-02
    2 KSP unpreconditioned resid norm 5.542834546634e-07 true resid norm 5.542457829175e-07 ||r(i)||/||b|| 1.156833383634e-03
    3 KSP unpreconditioned resid norm 1.069465737333e-07 true resid norm 1.069637665094e-07 ||r(i)||/||b|| 2.232570093468e-04
    4 KSP unpreconditioned resid norm 3.167781033623e-08 true resid norm 3.169346953559e-08 ||r(i)||/||b|| 6.615127211064e-05
 5 Nonlinear |R| = 3.920824e-05

The difference in the nonlinear and linear appears BEFORE the non-linear residual jump. Yea ok, I see it in Gamble's frictionless example now:

 0 Nonlinear |R| = 1.156441e+04
    0 KSP unpreconditioned resid norm 1.156441422999e+04 true resid norm 1.156441422999e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.896997695292e+03 true resid norm 1.896997695292e+03 ||r(i)||/||b|| 1.640375083049e-01
    2 KSP unpreconditioned resid norm 4.635223532666e+02 true resid norm 4.635223584851e+02 ||r(i)||/||b|| 4.008178445244e-02
    3 KSP unpreconditioned resid norm 1.344740928203e+02 true resid norm 1.344740498281e+02 ||r(i)||/||b|| 1.162826297587e-02
    4 KSP unpreconditioned resid norm 3.914128450523e+01 true resid norm 3.914129070054e+01 ||r(i)||/||b|| 3.384632366337e-03
    5 KSP unpreconditioned resid norm 2.093572162613e+01 true resid norm 2.093571798586e+01 ||r(i)||/||b|| 1.810356976973e-03
    6 KSP unpreconditioned resid norm 7.948710377665e+00 true resid norm 7.948713576559e+00 ||r(i)||/||b|| 6.873425163161e-04
    7 KSP unpreconditioned resid norm 2.443833325685e+00 true resid norm 2.443856329107e+00 ||r(i)||/||b|| 2.113255613735e-04
    8 KSP unpreconditioned resid norm 7.000606319562e-01 true resid norm 7.000716210966e-01 ||r(i)||/||b|| 6.053671264052e-05
    9 KSP unpreconditioned resid norm 1.284004334384e-01 true resid norm 1.283851336321e-01 ||r(i)||/||b|| 1.110174117589e-05
   10 KSP unpreconditioned resid norm 3.216353937173e-02 true resid norm 3.214272227083e-02 ||r(i)||/||b|| 2.779450963239e-06
   11 KSP unpreconditioned resid norm 1.193724034104e-02 true resid norm 1.197067846507e-02 ||r(i)||/||b|| 1.035130550238e-06
   12 KSP unpreconditioned resid norm 5.540335325403e-03 true resid norm 5.527515018594e-03 ||r(i)||/||b|| 4.779762215938e-07
 1 Nonlinear |R| = 4.159738e+01
    0 KSP unpreconditioned resid norm 1.626348731637e+04 true resid norm 1.626348731637e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.272855560343e+03 true resid norm 2.272855560345e+03 ||r(i)||/||b|| 1.397520418673e-01
    2 KSP unpreconditioned resid norm 6.001377005647e+02 true resid norm 6.001377019633e+02 ||r(i)||/||b|| 3.690092354050e-02
    3 KSP unpreconditioned resid norm 1.608064348230e+02 true resid norm 1.608064438551e+02 ||r(i)||/||b|| 9.887574585141e-03
    4 KSP unpreconditioned resid norm 7.621451694899e+01 true resid norm 7.621452253879e+01 ||r(i)||/||b|| 4.686234941880e-03
    5 KSP unpreconditioned resid norm 2.279501219400e+01 true resid norm 2.279503347433e+01 ||r(i)||/||b|| 1.401607971950e-03
    6 KSP unpreconditioned resid norm 7.514032870307e+00 true resid norm 7.514029137087e+00 ||r(i)||/||b|| 4.620183230643e-04
    7 KSP unpreconditioned resid norm 3.547806137979e+00 true resid norm 3.547829780153e+00 ||r(i)||/||b|| 2.181469269867e-04
    8 KSP unpreconditioned resid norm 1.318728592824e+00 true resid norm 1.318735784576e+00 ||r(i)||/||b|| 8.108567116778e-05
    9 KSP unpreconditioned resid norm 4.605053015476e-01 true resid norm 4.604968123355e-01 ||r(i)||/||b|| 2.831476443997e-05
   10 KSP unpreconditioned resid norm 1.577946682316e-01 true resid norm 1.578018661445e-01 ||r(i)||/||b|| 9.702830830484e-06
   11 KSP unpreconditioned resid norm 4.567006979969e-02 true resid norm 4.567037691790e-02 ||r(i)||/||b|| 2.808153997324e-06
   12 KSP unpreconditioned resid norm 1.392599441462e-02 true resid norm 1.393336072663e-02 ||r(i)||/||b|| 8.567265098556e-07
 2 Nonlinear |R| = 7.051010e+01
    0 KSP unpreconditioned resid norm 7.051009946857e+01 true resid norm 7.051009946857e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.762990670694e+01 true resid norm 4.762990670197e+01 ||r(i)||/||b|| 6.755047441565e-01
    2 KSP unpreconditioned resid norm 1.868190829187e+01 true resid norm 1.868190841005e+01 ||r(i)||/||b|| 2.649536527512e-01
    3 KSP unpreconditioned resid norm 5.552020468490e+00 true resid norm 5.552020550776e+00 ||r(i)||/||b|| 7.874078454890e-02
    4 KSP unpreconditioned resid norm 1.450845044848e+00 true resid norm 1.450845093875e+00 ||r(i)||/||b|| 2.057641536191e-02
    5 KSP unpreconditioned resid norm 3.782391419499e-01 true resid norm 3.782402065688e-01 ||r(i)||/||b|| 5.364340845065e-03
    6 KSP unpreconditioned resid norm 1.353022251310e-01 true resid norm 1.353025559080e-01 ||r(i)||/||b|| 1.918910296933e-03
    7 KSP unpreconditioned resid norm 6.591479020965e-02 true resid norm 6.591418659061e-02 ||r(i)||/||b|| 9.348190838958e-04
    8 KSP unpreconditioned resid norm 2.135320277714e-02 true resid norm 2.135282139184e-02 ||r(i)||/||b|| 3.028335167979e-04
    9 KSP unpreconditioned resid norm 8.603092060783e-03 true resid norm 8.603683929768e-03 ||r(i)||/||b|| 1.220205898816e-04
   10 KSP unpreconditioned resid norm 2.414808727770e-03 true resid norm 2.414551185262e-03 ||r(i)||/||b|| 3.424404735578e-05
   11 KSP unpreconditioned resid norm 5.878644910138e-04 true resid norm 5.872373510015e-04 ||r(i)||/||b|| 8.328414729628e-06
   12 KSP unpreconditioned resid norm 1.787192897719e-04 true resid norm 1.785460763661e-04 ||r(i)||/||b|| 2.532205708286e-06
   13 KSP unpreconditioned resid norm 8.064898654405e-05 true resid norm 8.036627284175e-05 ||r(i)||/||b|| 1.139783852916e-06
   14 KSP unpreconditioned resid norm 2.916973178606e-05 true resid norm 2.949701270160e-05 ||r(i)||/||b|| 4.183374144118e-07
 3 Nonlinear |R| = 4.138528e-02

And in the Gamble penalty example:

 0 Nonlinear |R| = 1.156441e+04
    0 KSP unpreconditioned resid norm 1.156441422995e+04 true resid norm 1.156441422995e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.127321108483e+01 true resid norm 4.127321108483e+01 ||r(i)||/||b|| 3.568984149491e-03
    2 KSP unpreconditioned resid norm 5.942924286114e-01 true resid norm 5.942678624979e-01 ||r(i)||/||b|| 5.138763197871e-05
    3 KSP unpreconditioned resid norm 8.196865808799e-03 true resid norm 8.143681645983e-03 ||r(i)||/||b|| 7.042018284760e-07
 1 Nonlinear |R| = 4.159728e+01
    0 KSP unpreconditioned resid norm 2.489921379810e+05 true resid norm 2.489921379810e+05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.682255939819e+03 true resid norm 1.682255440180e+03 ||r(i)||/||b|| 6.756259269150e-03
    2 KSP unpreconditioned resid norm 3.475328940525e+02 true resid norm 3.475329851678e+02 ||r(i)||/||b|| 1.395758870083e-03
    3 KSP unpreconditioned resid norm 2.802467958181e+02 true resid norm 2.802470187309e+02 ||r(i)||/||b|| 1.125525572829e-03
    4 KSP unpreconditioned resid norm 2.618527049893e+02 true resid norm 2.618536041161e+02 ||r(i)||/||b|| 1.051654105384e-03
    5 KSP unpreconditioned resid norm 2.261539117820e+02 true resid norm 2.261552807474e+02 ||r(i)||/||b|| 9.082828180088e-04
    6 KSP unpreconditioned resid norm 1.861661249814e+02 true resid norm 1.861667621989e+02 ||r(i)||/||b|| 7.476812870818e-04
    7 KSP unpreconditioned resid norm 3.471621561051e+01 true resid norm 3.471608650988e+01 ||r(i)||/||b|| 1.394264364786e-04
    8 KSP unpreconditioned resid norm 2.787056786458e-02 true resid norm 2.843882954681e-02 ||r(i)||/||b|| 1.142157731461e-07
 2 Nonlinear |R| = 2.573203e+03
    0 KSP unpreconditioned resid norm 2.573203431198e+03 true resid norm 2.573203431198e+03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.696535044472e+02 true resid norm 1.696534382603e+02 ||r(i)||/||b|| 6.593083010981e-02
    2 KSP unpreconditioned resid norm 2.175582749333e+00 true resid norm 2.175437199198e+00 ||r(i)||/||b|| 8.454198268283e-04
    3 KSP unpreconditioned resid norm 5.426933937643e-03 true resid norm 5.442604806649e-03 ||r(i)||/||b|| 2.115108638773e-06
    4 KSP unpreconditioned resid norm 3.203630197162e-05 true resid norm 2.932633631562e-04 ||r(i)||/||b|| 1.139681999490e-07
 3 Nonlinear |R| = 5.835372e+02

Yea, I think I may have raised a red herring! Gotta do those two checks tomorrow!:

1) See that initial residual norm matches the residual norm that appears from that extraneous FEProblemBase::computeResidual call
2) Run the Gamble example with the new feature and see whether the bad linear solve appears. We already see the bad linear solve appears in the single_point test.

Both the hypotheses were confirmed!

# 1/25/18

New moose formulation:

Time Step  1, time = 0.001
                dt = 0.001




*** Warning ***
Unchanged contact state. 0 nodes in contact.


 0 Nonlinear |R| = 5.471838e+03




*** Warning ***
Warning in PenetrationLocator. Penetration is not detected for one or more slave nodes. This could be because those slave nodes simply do not project to faces on the master surface. However, this could also be because contact should be enforced on those nodes, but the faces that they project to are outside the contact patch, which will give an erroneous result. Use appropriate options for 'patch_size' and 'patch_update_strategy' in the Mesh block to avoid this issue. Setting 'patch_update_strategy=iteration' is recommended because it completely avoids this potential issue. Also note that this warning is printed only once, so a similar situation could occur multiple times during the simulation but this warning is printed only at the first occurrence.

    0 KSP unpreconditioned resid norm 5.471838102659e+03 true resid norm 5.471838102659e+03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.678501661859e+02 true resid norm 4.678501661859e+02 ||r(i)||/||b|| 8.550146356826e-02
    2 KSP unpreconditioned resid norm 1.002293429030e+02 true resid norm 1.002293653303e+02 ||r(i)||/||b|| 1.831731192514e-02
    3 KSP unpreconditioned resid norm 1.177693095297e-01 true resid norm 1.177744290694e-01 ||r(i)||/||b|| 2.152374153983e-05




*** Warning ***
Unchanged contact state. 0 nodes in contact.


 1 Nonlinear |R| = 1.001988e+02
    0 KSP unpreconditioned resid norm 1.001988279477e+02 true resid norm 1.001988279477e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.758050384721e+01 true resid norm 1.758050384721e+01 ||r(i)||/||b|| 1.754561825453e-01
    2 KSP unpreconditioned resid norm 3.936647788769e+00 true resid norm 3.936652454164e+00 ||r(i)||/||b|| 3.928840820593e-02
    3 KSP unpreconditioned resid norm 5.890404501950e-06 true resid norm 1.075864981744e-05 ||r(i)||/||b|| 1.073730106210e-07




*** Warning ***
Changed contact state!!! 1 nodes in contact.


 2 Nonlinear |R| = 6.864072e+04
 Solve Did NOT Converge!

Old moose formulation:

 0 Nonlinear |R| = 5.471838e+03




*** Warning ***
Warning in PenetrationLocator. Penetration is not detected for one or more slave nodes. This could be because those slave nodes simply do not project to faces on the master surface. However, this could also be because contact should be enforced on those nodes, but the faces that they project to are outside the contact patch, which will give an erroneous result. Use appropriate options for 'patch_size' and 'patch_update_strategy' in the Mesh block to avoid this issue. Setting 'patch_update_strategy=iteration' is recommended because it completely avoids this potential issue. Also note that this warning is printed only once, so a similar situation could occur multiple times during the simulation but this warning is printed only at the first occurrence.

    0 KSP unpreconditioned resid norm 5.471838102659e+03 true resid norm 5.471838102659e+03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.678501661859e+02 true resid norm 4.678501661859e+02 ||r(i)||/||b|| 8.550146356826e-02
    2 KSP unpreconditioned resid norm 1.002293429030e+02 true resid norm 1.002293653303e+02 ||r(i)||/||b|| 1.831731192514e-02
    3 KSP unpreconditioned resid norm 1.177693095297e-01 true resid norm 1.177744290695e-01 ||r(i)||/||b|| 2.152374153984e-05
 1 Nonlinear |R| = 1.001988e+02
    0 KSP unpreconditioned resid norm 1.001988279477e+02 true resid norm 1.001988279477e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.758050384494e+01 true resid norm 1.758050384494e+01 ||r(i)||/||b|| 1.754561825226e-01
    2 KSP unpreconditioned resid norm 3.936647766797e+00 true resid norm 3.936652416703e+00 ||r(i)||/||b|| 3.928840783207e-02
    3 KSP unpreconditioned resid norm 5.889641465133e-06 true resid norm 1.071253408567e-05 ||r(i)||/||b|| 1.069127683935e-07
 2 Nonlinear |R| = 5.952291e-02
    0 KSP unpreconditioned resid norm 6.864072015663e+04 true resid norm 6.864072015663e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.128325520873e-01 true resid norm 2.128325520873e-01 ||r(i)||/||b|| 3.100674812293e-06
 3 Nonlinear |R| = 2.128319e-01
    0 KSP unpreconditioned resid norm 2.128318588363e-01 true resid norm 2.128318588363e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.200018646607e-02 true resid norm 8.200018646607e-02 ||r(i)||/||b|| 3.852815406228e-01
    2 KSP unpreconditioned resid norm 1.454185608691e-02 true resid norm 1.454184957463e-02 ||r(i)||/||b|| 6.832553008810e-02
    3 KSP unpreconditioned resid norm 1.937920074349e-03 true resid norm 1.938542579923e-03 ||r(i)||/||b|| 9.108328943433e-03
    4 KSP unpreconditioned resid norm 1.592210243984e-05 true resid norm 1.778242896603e-05 ||r(i)||/||b|| 8.355153717711e-05
 4 Nonlinear |R| = 6.479389e-05
 Solve Converged!

Ok, here's the new solve after removing the divergence threshold criteria:

 0 Nonlinear |R| = 5.471838e+03

*** Warning ***
Warning in PenetrationLocator. Penetration is not detected for one or more slave nodes. This could be because those slave nodes simply do not project to faces on the master surface. However, this could also be because contact should be enforced on those nodes, but the faces that they project to are outside the contact patch, which will give an erroneous result. Use appropriate options for 'patch_size' and 'patch_update_strategy' in the Mesh block to avoid this issue. Setting 'patch_update_strategy=iteration' is recommended because it completely avoids this potential issue. Also note that this warning is printed only once, so a similar situation could occur multiple times during the simulation but this warning is printed only at the first occurrence.

    0 KSP unpreconditioned resid norm 5.471838102659e+03 true resid norm 5.471838102659e+03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.678501661859e+02 true resid norm 4.678501661859e+02 ||r(i)||/||b|| 8.550146356826e-02
    2 KSP unpreconditioned resid norm 1.002293429030e+02 true resid norm 1.002293653303e+02 ||r(i)||/||b|| 1.831731192514e-02
    3 KSP unpreconditioned resid norm 1.177693095297e-01 true resid norm 1.177744290694e-01 ||r(i)||/||b|| 2.152374153983e-05

*** Warning ***
Unchanged contact state. 0 nodes in contact.

1 Nonlinear |R| = 1.001988e+02
    0 KSP unpreconditioned resid norm 1.001988279477e+02 true resid norm 1.001988279477e+02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.758050384721e+01 true resid norm 1.758050384721e+01 ||r(i)||/||b|| 1.754561825453e-01
    2 KSP unpreconditioned resid norm 3.936647788769e+00 true resid norm 3.936652454164e+00 ||r(i)||/||b|| 3.928840820593e-02
    3 KSP unpreconditioned resid norm 5.890404501950e-06 true resid norm 1.075864981744e-05 ||r(i)||/||b|| 1.073730106210e-07

*** Warning ***
Changed contact state!!! 1 nodes in contact.
 2 Nonlinear |R| = 6.864072e+04
    0 KSP unpreconditioned resid norm 6.864071920856e+04 true resid norm 6.864071920856e+04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.128325469796e-01 true resid norm 2.128325469796e-01 ||r(i)||/||b|| 3.100674780707e-06




*** Warning ***
Unchanged contact state. 1 nodes in contact.


 3 Nonlinear |R| = 2.128319e-01
    0 KSP unpreconditioned resid norm 2.128318538346e-01 true resid norm 2.128318538346e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.200018517062e-02 true resid norm 8.200018517062e-02 ||r(i)||/||b|| 3.852815435905e-01
    2 KSP unpreconditioned resid norm 1.454185596452e-02 true resid norm 1.454184940907e-02 ||r(i)||/||b|| 6.832553091594e-02
    3 KSP unpreconditioned resid norm 1.937920511222e-03 true resid norm 1.938141877382e-03 ||r(i)||/||b|| 9.106446438644e-03
    4 KSP unpreconditioned resid norm 1.592182732237e-05 true resid norm 5.384435098530e-05 ||r(i)||/||b|| 2.529900953038e-04




*** Warning ***
Unchanged contact state. 1 nodes in contact.


 4 Nonlinear |R| = 6.479389e-05
 Solve Converged!

And then after removing the residual call, the same result!

# 2/8/17

What does our residual object expect?

It expects the solution, residual vector to fill, and the
`NonlinearImplicitSystem & sys`

# 2/9/18

good commit: 550de36fb
timing: 49.643s

bad commit: 2ac908c3f
timing: 1m59.923s

good commit: 70e48bcf8
timing = 41.7s

--with-thread-model=none

libmesh update commit cdeee5648 --with-thread-model=none
timing = 1m27.849s

2ac908c3f --with-thread-model=none
timing = 1m18.523s
Pyc tests:
test:pyc_irr_strain.test_dense ................................................................... OK [10.15s]
test:pyc_irr_strain.test_buffer .................................................................. OK [14.84s]

550de36fb --with-thread-model=none
timing = 1m27.232s
Pyc tests:
test:pyc_irr_strain.test_dense ................................................................... OK [11.39s]
test:pyc_irr_strain.test_buffer .................................................................. OK [15.22s]

835d924d4 --with-thread-model=none
timing = 1m28.172s
Pyc tests:
test:pyc_irr_strain.test_dense ................................................................... OK [11.36s]
test:pyc_irr_strain.test_buffer .................................................................. OK [14.36s]

835d924d4 --with-thread-model=
Pyc tests:
test:pyc_irr_strain.test_dense ................................................................... OK [1.296s]
test:pyc_irr_strain.test_buffer .................................................................. OK [1.749s]


2ac908c3f --with-thread-model=
Pyc tests:
test:pyc_irr_strain.test_dense ................................................................... OK [13.87s]
test:pyc_irr_strain.test_buffer .................................................................. OK [17.77s]

# 2/15/18

Commit message for app fixes:

Updates to MooseVariableInterface and variable getter APIs to incorporate vector finite elements

Needs to wait for https://github.com/idaholab/moose/pull/10238

issue message:

Update MooseVariableInterface and other APIs for incorporation of vector finite elements

issue content:

The [MOOSE vector PR](https://github.com/idaholab/moose/pull/10238) is about ready to go and it has some API changes that will require some changes in apps.

# 2/22/18

Ok, currently I am able to run all tests flawlessly on cone both with gcc and
clang compilations. Currently building a gcc version on Mac and going to test
that. I'm beginning to hypothesize that it could be an environment issue...well
no it can't be that because I just ran with gcc build on my Mac, and everything
passed! Wtf!!!

# 3/5/18

Total non-linear iterations with `num_cuts = 2` and `contact_ltol = ltol`: 149
Active time = 5.78 seconds

Total non-linear iterations with `num_cuts = 0` and `contact_ltol = .5`: 132
Active time = 4.47 seconds

sliding block problem:

2 lambda-cuts, 10 steps, .5 contact tolerance, asm
NEWTON, FD: 10 steps, 59 non-linear its
PJFNK, FD: 10 steps, 60 non-linear its
PJFNK, SMP: 10 steps, 65 non-linear its
NEWTON, SMP: 10 steps, 60 non-linear its

2 lambda-cuts, 10 steps, 1e-6 contact tolerance, asm
NEWTON, FD: 10 steps, 48 non-linear its
PJFNK, FD: 10 steps, 48 non-linear its
PJFNK, SMP: 10 steps, 51 non-linear its
NEWTON, SMP: 10 steps, 53 non-linear its

2 lambda-cuts, 30 steps, 1e-6 contact tolerance, asm
PJFNK, FD: 30 steps, 134 non-linear its
PJFNK, SMP: 30 steps, 149 non-linear its, active 9.18
NEWTON, FD: 30 steps, 134 non-linear its

2 lambda-cuts, 30 steps, .5 contact tolerance, asm
PJFNK, FD: 30 steps, 156 non-linear its
PJFNK, SMP: 30 steps, 166 non-linear its, active 8.66

0 lambda-cuts, 30 steps, 1e-6 contact tolerance, asm
PJFNK, SMP: 35 steps, 151 non-linear its, active 23.28
PJFNK, FD: witness jamming! 35 steps, 147 non-linear its

0 lambda-cuts, 30 steps, .5 contact tolerance, asm
PJFNK, SMP: 33 steps, 166 non-linear its, active 9.64

2 lambda-cuts, 30 steps, 1e-6 contact tolerance, lu
NEWTON, FD: 30 steps, 134 non-linear its

Verified through `-snes_view` and through `-mat_mffd_type` that NEWTON and PJFNK
with `-snes_fd` are really truly different

I think the clearest test of whether a linear system is well scaled is to
use PJFNK with a finite differenced preconditioner (most accurately with full
differencing, less accurately with coloring) (and perhaps to use a line
search because then the solve is likely to demonstrate bad behavior earlier) and
then observe whether the true residual matches the unpreconditioned residual. So
Zapdos I believe is poorly scaled and eventually I want to develop an
algorithm/diagnostic that reveals and suggests how to scale one's variables to
create the best conditioned system.

Another clear way to tell whether a system is well scaled is to run with
`-pc_type svd -pc_svd_monitor`!! I tested that when I run with an `asm`
preconditioner, then the Zapdos problem converges. However, as soon as I turn to
the `svd` "preconditioner", then the problem no longer converges and I see
absolutely horrendous condition numbers!:

Time Step  1, time = 1e-11
                dt = 1e-11
    |residual|_2 of individual variables:
                  potential: 88500
                  em:        0.676043
                  emliq:     860.265
                  Arp:       0.427591
                  mean_en:   40.3082
                  OHm:       4274.12

 0 Nonlinear |R| = 8.860733e+04
      SVD: condition number 2.861482635029e+18, 0 of 1039 singular values are (nearly) zero
      SVD: smallest singular values: 6.605894949746e-08 3.360421490264e-07 9.160910535765e-07 1.566311893390e-06 2.211462057380e-06
      SVD: largest singular values : 5.152086167118e+10 6.803285617521e+10 9.166842564772e+10 1.275351988992e+11 1.890265368752e+11
    0 KSP unpreconditioned resid norm 8.860733466940e+04 true resid norm 8.860733466940e+04 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 8.860733466940e+04 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.555721572811e-04 true resid norm 8.776899730270e-04 ||r(i)||/||b|| 9.905387362138e-09
    1 KSP Residual norm 9.555721572811e-04 % max 1.000000007606e+00 min 1.000000007606e+00 max/min 1.000000000000e+00
      Line search: Using full step: fnorm 8.860733466940e+04 gnorm 3.006870853434e+04
    |residual|_2 of individual variables:
                  potential: 1.8996
                  em:        247.828
                  emliq:     11356.5
                  Arp:       247.942
                  mean_en:   27839.4
                  OHm:       0.000686369

 1 Nonlinear |R| = 3.006871e+04
      SVD: condition number 5.058219757856e+19, 0 of 1039 singular values are (nearly) zero
      SVD: smallest singular values: 3.737016564828e-09 3.269299464171e-08 2.863796025981e-07 5.908120410986e-07 1.047076750051e-06
      SVD: largest singular values : 5.152084352599e+10 6.803283729982e+10 9.166840542000e+10 1.275351763691e+11 1.890265102365e+11
    0 KSP unpreconditioned resid norm 3.006870853434e+04 true resid norm 3.006870853434e+04 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 3.006870853434e+04 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.742836466488e-02 true resid norm 1.739977094151e-02 ||r(i)||/||b|| 5.786670525487e-07
    1 KSP Residual norm 1.742836466488e-02 % max 9.999999983198e-01 min 9.999999983198e-01 max/min 1.000000000000e+00
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1

Yes, definitely need a good diagostic/algorithm for creating a good matrix!!!

# 3/6/18

So even for an ideal preconditioner and approximation of the Jacobian, which
means that the linear solve is optimal, a line-search can be absolutely
necessary in contact to avoid jamming.

# 3/7/18

So even with diagonal scaling kind of like how I think it should be, the initial condition number
is still terrible:

 0 Nonlinear |R| = 8.624066e+02
The minimum diagonal entry for variable 0 is 0.374576
The maximum diagonal entry for variable 0 is 16159.2
The minimum diagonal occurs at node id 102 which has coordinates (0.5,0,0,)
The maximum diagonal occurs at node id 210 which has coordinates (0.999999,0,0,)
The minimum diagonal entry for variable 1 is 0.1416
The maximum diagonal entry for variable 1 is 11736.7
The maximum diagonal occurs at node id 212 which has coordinates (1,0,0,)
The minimum diagonal occurs at node id 270 which has coordinates (1.6,0,0,)
The minimum diagonal entry for variable 2 is 0.376852
The maximum diagonal entry for variable 2 is 616.21
The minimum diagonal occurs at node id 76 which has coordinates (0.0378194,0,0,)
The maximum diagonal occurs at node id 210 which has coordinates (0.999999,0,0,)
The minimum diagonal entry for variable 3 is 0.380591
The maximum diagonal entry for variable 3 is 565.643
The maximum diagonal occurs at node id 212 which has coordinates (1,0,0,)
The minimum diagonal occurs at node id 252 which has coordinates (1.01711,0,0,)
The minimum diagonal entry for variable 4 is 0.129553
The maximum diagonal entry for variable 4 is 247.431
The maximum diagonal occurs at node id 102 which has coordinates (0.5,0,0,)
The minimum diagonal occurs at node id 186 which has coordinates (0.99989,0,0,)
The minimum diagonal entry for variable 5 is 0.210796
The maximum diagonal entry for variable 5 is 279.172
The minimum diagonal occurs at node id 131 which has coordinates (0.971904,0,0,)
The maximum diagonal occurs at node id 210 which has coordinates (0.999999,0,0,)
The minimum diagonal entry for variable 6 is 0.895257
The maximum diagonal entry for variable 6 is 1466.66
The maximum diagonal occurs at node id 212 which has coordinates (1,0,0,)
The minimum diagonal occurs at node id 253 which has coordinates (1.02087,0,0,)
      SVD: condition number 1.799242736290e+13, 0 of 1040 singular values are (nearly) zero
      SVD: smallest singular values: 5.283618013495e-05 2.022141576985e-04 2.097912129974e-04 3.519080816676e-04 5.001484338878e-04
      SVD: largest singular values : 3.499345949675e+07 4.710631597814e+07 6.545521519716e+07 9.683473798936e+07 9.506511332114e+08
    0 KSP unpreconditioned resid norm 8.624065720871e+02 true resid norm 8.624065720871e+02 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 8.624065720871e+02 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.741274620680e-06 true resid norm 1.736671811022e-06 ||r(i)||/||b|| 2.013750668457e-09
    1 KSP Residual norm 1.741274620680e-06 % max 1.000000000973e+00 min 1.000000000973e+00 max/min 1.000000000000e+00
      Line search: Using full step: fnorm 8.624065720871e+02 gnorm 7.307267181252e+02
    |residual|_2 of individual variables:
               potential:    16.5556
               potentialliq: 2.43247e-09
               em:           0.665784
               emliq:        727.695
               Arp:          63.7654
               mean_en:      8.97804
               OHm:          0.0102765

This is undoubtedly why IMO, `-snes_fd` with NEWTON fails. `-snes_fd` with PJFNK
succeeds for the above conditioning (includes split potential variables and
manual scaling to get all the minimum diagonal entries between .1 and
1). `-snes_fd` with PJFNK for non-split potential variables and no manual
scaling fails. Its condition number is 2.861482635029e+18

# 3/12/18

Ok, it's not really that there is a large non-uniformity in the mesh...it's that
the mesh is just very fine. See here for a uniformly fine mesh:

Time Step  2, time = 2.2
                dt = 1.2
    |residual|_2 of individual variables:
                  potential: 3.94887

 0 Nonlinear |R| = 3.948868e+00
    0 KSP unpreconditioned resid norm 3.948868221119e+00 true resid norm 3.948868221119e+00 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 3.948868221119e+00 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.990356128632e-01 true resid norm 1.990356128632e-01 ||r(i)||/||b|| 5.040320459384e-02
    1 KSP Residual norm 1.990356128632e-01 % max 1.001232121095e+00 min 1.001232121095e+00 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 7.494869387932e-06 true resid norm 1.940997537159e-01 ||r(i)||/||b|| 4.915326185813e-02
    2 KSP Residual norm 7.494869387932e-06 % max 1.025538261374e+00 min 9.750729425046e-01 max/min 1.051755429435e+00
  Linear solve converged due to CONVERGED_RTOL iterations 2
      Line search: Using full step: fnorm 3.948868221119e+00 gnorm 1.992859727569e-01
    |residual|_2 of individual variables:
                  potential: 0.199286

 1 Nonlinear |R| = 1.992860e-01
    0 KSP unpreconditioned resid norm 1.992859727569e-01 true resid norm 1.992859727569e-01 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 1.992859727569e-01 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.988633256152e-07 true resid norm 4.988633256139e-07 ||r(i)||/||b|| 2.503253584348e-06
    1 KSP Residual norm 4.988633256152e-07 % max 9.999999877377e-01 min 9.999999877377e-01 max/min 1.000000000000e+00
  Linear solve converged due to CONVERGED_RTOL iterations 1
      Line search: Using full step: fnorm 1.992859727569e-01 gnorm 2.439376445365e-09
    |residual|_2 of individual variables:
                  potential: 2.43938e-09

 2 Nonlinear |R| = 2.439376e-09
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 2
 Solve Converged!

Time Step  3, time = 3.64
                dt = 1.44
    |residual|_2 of individual variables:
                  potential: 1.30316

 0 Nonlinear |R| = 1.303162e+00
    0 KSP unpreconditioned resid norm 1.303161773291e+00 true resid norm 1.303161773291e+00 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 1.303161773291e+00 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.693458093212e-01 true resid norm 1.693458093212e-01 ||r(i)||/||b|| 1.299499515655e-01
    1 KSP Residual norm 1.693458093212e-01 % max 1.008552587316e+00 min 1.008552587316e+00 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 4.582570360308e-07 true resid norm 1.866187828918e-01 ||r(i)||/||b|| 1.432046171985e-01
    2 KSP Residual norm 4.582570360308e-07 % max 1.067676229802e+00 min 9.366148698752e-01 max/min 1.139930898112e+00
  Linear solve converged due to CONVERGED_RTOL iterations 2
      Line search: Using full step: fnorm 1.303161773291e+00 gnorm 1.707939094624e-01
    |residual|_2 of individual variables:
                  potential: 0.170794

This small grid phenomenon appears regardless of whether the SMP or FDP
preconditioner is being used. If time is brought closer to the scale of
position, then the solve phenonemon improves.

h = .05, dt = 100
 0 Nonlinear |R| = 6.233634e-01
    0 KSP unpreconditioned resid norm 6.233634244705e-01 true resid norm 6.233634244705e-01 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 6.233634244705e-01 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.220398834842e-01 true resid norm 6.220398834842e-01 ||r(i)||/||b|| 9.978767747122e-01
    1 KSP Residual norm 6.220398834842e-01 % max 1.282619683911e+01 min 1.282619683911e+01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 4.462370692421e-01 true resid norm 4.464352564690e-01 ||r(i)||/||b|| 7.161717209318e-01
    2 KSP Residual norm 4.462370692421e-01 % max 1.285706033192e+01 min 8.185998316272e-02 max/min 1.570616048914e+02
    3 KSP unpreconditioned resid norm 6.852284705053e-02 true resid norm 6.891985496232e-02 ||r(i)||/||b|| 1.105612749430e-01
    3 KSP Residual norm 6.852284705053e-02 % max 1.285710477534e+01 min 5.701834422349e-02 max/min 2.254906723518e+02
    4 KSP unpreconditioned resid norm 8.408012471952e-03 true resid norm 1.640573974848e-02 ||r(i)||/||b|| 2.631809808607e-02
    4 KSP Residual norm 8.408012471952e-03 % max 1.285711729096e+01 min 5.666838090950e-02 max/min 2.268834416055e+02
    5 KSP unpreconditioned resid norm 7.383323102654e-04 true resid norm 1.200371623552e-02 ||r(i)||/||b|| 1.925636918097e-02
    5 KSP Residual norm 7.383323102654e-04 % max 1.285721649731e+01 min 5.666604583012e-02 max/min 2.268945416776e+02
    6 KSP unpreconditioned resid norm 7.227659694176e-05 true resid norm 1.293304604821e-02 ||r(i)||/||b|| 2.074720065457e-02
    6 KSP Residual norm 7.227659694176e-05 % max 1.285789863589e+01 min 5.666258319832e-02 max/min 2.269204457355e+02
    7 KSP unpreconditioned resid norm 1.001713886164e-05 true resid norm 1.450804227625e-02 ||r(i)||/||b|| 2.327381060024e-02
    7 KSP Residual norm 1.001713886164e-05 % max 1.285791493188e+01 min 5.666254244926e-02 max/min 2.269208965234e+02
    8 KSP unpreconditioned resid norm 1.045057641350e-06 true resid norm 1.383017353522e-02 ||r(i)||/||b|| 2.218637313693e-02
    8 KSP Residual norm 1.045057641350e-06 % max 1.285829527674e+01 min 5.666064236423e-02 max/min 2.269352188789e+02
  Linear solve converged due to CONVERGED_RTOL iterations 8
      Line search: Using full step: fnorm 6.233634244705e-01 gnorm 1.330874176494e-02
    |residual|_2 of individual variables:
                  potential: 0.0133087

 1 Nonlinear |R| = 1.330874e-02

h = .05, dt = 10
 0 Nonlinear |R| = 3.485273e+00
    0 KSP unpreconditioned resid norm 3.485273240218e+00 true resid norm 3.485273240218e+00 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 3.485273240218e+00 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.411066372603e+00 true resid norm 3.411066372603e+00 ||r(i)||/||b|| 9.787084505286e-01
    1 KSP Residual norm 3.411066372603e+00 % max 4.302042325679e+00 min 4.302042325679e+00 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 1.019511060834e+00 true resid norm 1.019211844756e+00 ||r(i)||/||b|| 2.924338422007e-01
    2 KSP Residual norm 1.019511060834e+00 % max 4.387670156309e+00 min 1.834817359803e-01 max/min 2.391338916033e+01
    3 KSP unpreconditioned resid norm 9.135445787215e-02 true resid norm 9.261009656829e-02 ||r(i)||/||b|| 2.657183244620e-02
    3 KSP Residual norm 9.135445787215e-02 % max 4.388377605805e+00 min 1.736889283692e-01 max/min 2.526573021671e+01
    4 KSP unpreconditioned resid norm 9.410242848801e-03 true resid norm 1.192180389879e-02 ||r(i)||/||b|| 3.420622452559e-03
    4 KSP Residual norm 9.410242848801e-03 % max 4.388487439697e+00 min 1.735380239321e-01 max/min 2.528833358972e+01
    5 KSP unpreconditioned resid norm 6.158872503091e-04 true resid norm 5.113914772719e-03 ||r(i)||/||b|| 1.467292352780e-03
    5 KSP Residual norm 6.158872503091e-04 % max 4.388833902130e+00 min 1.735323419247e-01 max/min 2.529115814061e+01
    6 KSP unpreconditioned resid norm 4.125683048570e-05 true resid norm 5.103559871032e-03 ||r(i)||/||b|| 1.464321308338e-03
    6 KSP Residual norm 4.125683048570e-05 % max 4.388900275474e+00 min 1.735294838547e-01 max/min 2.529195718204e+01
    7 KSP unpreconditioned resid norm 5.291593437160e-06 true resid norm 5.034838734685e-03 ||r(i)||/||b|| 1.444603733385e-03
    7 KSP Residual norm 5.291593437160e-06 % max 4.390209653942e+00 min 1.734776325484e-01 max/min 2.530706460222e+01
  Linear solve converged due to CONVERGED_RTOL iterations 7
      Line search: Using full step: fnorm 3.485273240218e+00 gnorm 5.286777749286e-03
    |residual|_2 of individual variables:
                  potential: 0.00528678

 1 Nonlinear |R| = 5.286778e-03

h = .05, dt = 1
 0 Nonlinear |R| = 1.928002e+01
    0 KSP unpreconditioned resid norm 1.928002299202e+01 true resid norm 1.928002299202e+01 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 1.928002299202e+01 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.501195750502e+01 true resid norm 1.501195750502e+01 ||r(i)||/||b|| 7.786275727592e-01
    1 KSP Residual norm 1.501195750502e+01 % max 1.490951317992e+00 min 1.490951317992e+00 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 8.214128356460e-01 true resid norm 8.216029587267e-01 ||r(i)||/||b|| 4.261421052593e-02
    2 KSP Residual norm 8.214128356460e-01 % max 1.669252069200e+00 min 5.049013893626e-01 max/min 3.306095218530e+00
    3 KSP unpreconditioned resid norm 4.946714533380e-02 true resid norm 4.865541106670e-02 ||r(i)||/||b|| 2.523617896454e-03
    3 KSP Residual norm 4.946714533380e-02 % max 1.669786057563e+00 min 5.032973535462e-01 max/min 3.317692902213e+00
    4 KSP unpreconditioned resid norm 2.919334775995e-03 true resid norm 3.619462313189e-03 ||r(i)||/||b|| 1.877312239040e-04
    4 KSP Residual norm 2.919334775995e-03 % max 1.670892671226e+00 min 5.026078528146e-01 max/min 3.324446010679e+00
    5 KSP unpreconditioned resid norm 1.171107057589e-04 true resid norm 2.417266600423e-03 ||r(i)||/||b|| 1.253767488464e-04
    5 KSP Residual norm 1.171107057589e-04 % max 1.671131162382e+00 min 5.025577227676e-01 max/min 3.325252178355e+00
  Linear solve converged due to CONVERGED_RTOL iterations 5
      Line search: Using full step: fnorm 1.928002299202e+01 gnorm 2.232151562069e-03
    |residual|_2 of individual variables:
                  potential: 0.00223215

 1 Nonlinear |R| = 2.232152e-03

h = .05, dt = .1
 0 Nonlinear |R| = 1.028467e+02
    0 KSP unpreconditioned resid norm 1.028466871259e+02 true resid norm 1.028466871259e+02 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 1.028466871259e+02 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.827037910846e+01 true resid norm 1.827037910846e+01 ||r(i)||/||b|| 1.776467440910e-01
    1 KSP Residual norm 1.827037910846e+01 % max 9.944447062505e-01 min 9.944447062505e-01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 3.974700420275e-01 true resid norm 3.970332708070e-01 ||r(i)||/||b|| 3.860438113296e-03
    2 KSP Residual norm 3.974700420275e-01 % max 1.054373521860e+00 min 8.667049344528e-01 max/min 1.216531116817e+00
    3 KSP unpreconditioned resid norm 8.424162168393e-03 true resid norm 8.735224475124e-03 ||r(i)||/||b|| 8.493442734264e-05
    3 KSP Residual norm 8.424162168393e-03 % max 1.056500933233e+00 min 8.628067988078e-01 max/min 1.224493055332e+00
    4 KSP unpreconditioned resid norm 1.494087515474e-04 true resid norm 5.023211037557e-04 ||r(i)||/||b|| 4.884173888274e-06
    4 KSP Residual norm 1.494087515474e-04 % max 1.056513117891e+00 min 8.626251719351e-01 max/min 1.224764999056e+00
  Linear solve converged due to CONVERGED_RTOL iterations 4
      Line search: Using full step: fnorm 1.028466871259e+02 gnorm 4.011681697041e-04
    |residual|_2 of individual variables:
                  potential: 0.000401168

 1 Nonlinear |R| = 4.011682e-04

h = .05, dt = .01
 0 Nonlinear |R| = 4.933394e+02
    0 KSP unpreconditioned resid norm 4.933393830700e+02 true resid norm 4.933393830700e+02 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 4.933393830700e+02 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.642657860606e+01 true resid norm 1.642657860606e+01 ||r(i)||/||b|| 3.329671047918e-02
    1 KSP Residual norm 1.642657860606e+01 % max 9.876947204243e-01 min 9.876947204243e-01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 4.485368903029e-01 true resid norm 4.486181492960e-01 ||r(i)||/||b|| 9.093499620976e-04
    2 KSP Residual norm 4.485368903029e-01 % max 9.940462814408e-01 min 9.201436575439e-01 max/min 1.080316397653e+00
    3 KSP unpreconditioned resid norm 9.104949798762e-03 true resid norm 8.967250763816e-03 ||r(i)||/||b|| 1.817663675666e-05
    3 KSP Residual norm 9.104949798762e-03 % max 9.987911020783e-01 min 9.133743829176e-01 max/min 1.093517752149e+00
    4 KSP unpreconditioned resid norm 1.260918451614e-04 true resid norm 2.558054085775e-04 ||r(i)||/||b|| 5.185181182691e-07
    4 KSP Residual norm 1.260918451614e-04 % max 1.000406867225e+00 min 9.110163641554e-01 max/min 1.098121731493e+00
  Linear solve converged due to CONVERGED_RTOL iterations 4
      Line search: Using full step: fnorm 4.933393830700e+02 gnorm 2.083801676191e-04
    |residual|_2 of individual variables:
                  potential: 0.00020838

 1 Nonlinear |R| = 2.083802e-04

h = .05, dt = .001
 0 Nonlinear |R| = 1.773842e+03
    0 KSP unpreconditioned resid norm 1.773841528311e+03 true resid norm 1.773841528311e+03 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 1.773841528311e+03 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.362685165812e+01 true resid norm 1.362685165812e+01 ||r(i)||/||b|| 7.682113334608e-03
    1 KSP Residual norm 1.362685165812e+01 % max 9.968904682400e-01 min 9.968904682400e-01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 1.363411691300e-01 true resid norm 1.362772801451e-01 ||r(i)||/||b|| 7.682607378959e-05
    2 KSP Residual norm 1.363411691300e-01 % max 9.994158763718e-01 min 9.710264833812e-01 max/min 1.029236476529e+00
    3 KSP unpreconditioned resid norm 6.796909201131e-04 true resid norm 8.925158560244e-04 ||r(i)||/||b|| 5.031542230686e-07
    3 KSP Residual norm 6.796909201131e-04 % max 9.998650545169e-01 min 9.605771610787e-01 max/min 1.040900299351e+00
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 1.773841528311e+03 gnorm 8.872656183974e-04
    |residual|_2 of individual variables:
                  potential: 0.000887266

 1 Nonlinear |R| = 8.872656e-04

h = .05, dt = 1e-4
 0 Nonlinear |R| = 3.504529e+03
    0 KSP unpreconditioned resid norm 3.504528942674e+03 true resid norm 3.504528942674e+03 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 3.504528942674e+03 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.815073345721e+01 true resid norm 3.815073345721e+01 ||r(i)||/||b|| 1.088612309422e-02
    1 KSP Residual norm 3.815073345721e+01 % max 9.898843692129e-01 min 9.898843692129e-01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 3.415607223931e-01 true resid norm 3.415685009317e-01 ||r(i)||/||b|| 9.746488230488e-05
    2 KSP Residual norm 3.415607223931e-01 % max 9.986383243515e-01 min 9.787991006767e-01 max/min 1.020268943505e+00
    3 KSP unpreconditioned resid norm 1.003061571183e-03 true resid norm 1.082585661617e-03 ||r(i)||/||b|| 3.089104639526e-07
    3 KSP Residual norm 1.003061571183e-03 % max 9.999214765593e-01 min 9.687189057860e-01 max/min 1.032210139171e+00
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 3.504528942674e+03 gnorm 1.042671052763e-03
    |residual|_2 of individual variables:
                  potential: 0.00104267

 1 Nonlinear |R| = 1.042671e-03

h = .05, dt = 1e-5
 0 Nonlinear |R| = 4.087688e+03
    0 KSP unpreconditioned resid norm 4.087687928606e+03 true resid norm 4.087687928606e+03 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 4.087687928606e+03 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.002964948858e+02 true resid norm 1.002964948858e+02 ||r(i)||/||b|| 2.453624069096e-02
    1 KSP Residual norm 1.002964948858e+02 % max 9.833720656738e-01 min 9.833720656738e-01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 2.558637581460e+00 true resid norm 2.558671805504e+00 ||r(i)||/||b|| 6.259459748867e-04
    2 KSP Residual norm 2.558637581460e+00 % max 9.957345640558e-01 min 9.500398298468e-01 max/min 1.048097703668e+00
    3 KSP unpreconditioned resid norm 2.932580896138e-02 true resid norm 2.934446149590e-02 ||r(i)||/||b|| 7.178743095954e-06
    3 KSP Residual norm 2.932580896138e-02 % max 9.997654351205e-01 min 9.349428929475e-01 max/min 1.069333156776e+00
  Linear solve converged due to CONVERGED_RTOL iterations 3
      Line search: Using full step: fnorm 4.087687928606e+03 gnorm 2.933907024986e-02
    |residual|_2 of individual variables:
                  potential: 0.0293391

 1 Nonlinear |R| = 2.933907e-02

To summarize:
h = .05, `-pc_type amg`, SMP preconditioner, converged norms:

dt = 1e2:  8 KSP unpreconditioned resid norm 1.045057641350e-06 true resid norm 1.383017353522e-02 ||r(i)||/||b|| 2.218637313693e-02
dt = 1e1:  7 KSP unpreconditioned resid norm 5.291593437160e-06 true resid norm 5.034838734685e-03 ||r(i)||/||b|| 1.444603733385e-03
dt = 1e0:  5 KSP unpreconditioned resid norm 1.171107057589e-04 true resid norm 2.417266600423e-03 ||r(i)||/||b|| 1.253767488464e-04
dt = 1e-1: 4 KSP unpreconditioned resid norm 1.494087515474e-04 true resid norm 5.023211037557e-04 ||r(i)||/||b|| 4.884173888274e-06
dt = 1e-2: 4 KSP unpreconditioned resid norm 1.260918451614e-04 true resid norm 2.558054085775e-04 ||r(i)||/||b|| 5.185181182691e-07
dt = 1e-3: 3 KSP unpreconditioned resid norm 6.796909201131e-04 true resid norm 8.925158560244e-04 ||r(i)||/||b|| 5.031542230686e-07
dt = 1e-4: 3 KSP unpreconditioned resid norm 1.003061571183e-03 true resid norm 1.082585661617e-03 ||r(i)||/||b|| 3.089104639526e-07
dt = 1e-5: 3 KSP unpreconditioned resid norm 2.932580896138e-02 true resid norm 2.934446149590e-02 ||r(i)||/||b|| 7.178743095954e-06

Analagous qualitative trends with `-pc_type lu` although one has to go to larger time steps in order to be able to see it.

# 3/13/18

h = .01, D = 1

dt = 1e-9, lin = 17 (11, 6)
dt = 1e-8, lin = 16 (10, 6)
dt = 1e-7, lin = 15 (9, 6)
dt = 1e-6, lin = 15 (9, 6)
dt = 1e-5, lin = 10 (6, 4)
dt = 1e-4, lin = 13 (7, 6)
dt = 1e-3, lin = 18 (9, 9)
dt = 1e-2, lin = 20 (10, 10)
dt = 1e-1, lin = 21 (11, 10)
dt = 1e0, lin = 23 (12, 11)
dt = 1e1, lin = 24 (12, 12)
dt = 1e2, lin = 23 (12, 11)
dt = 1e3, lin = 24 (12, 12)
dt = 1e4, lin = 24 (12, 12)
dt = 1e7, lin = 23 (12, 11)

mesh size, finish time
50, 419973

# 3/14/18

Need to consider the DirichletBCs when performing `automatic_scaling`.

# 3/20/18

I'm able to get Zapdos `mean_en.i` to run successfully with FDP and PJFNK by specifying `-mat_mffd_err 1e-6`. Boom! Nice!

Testing on the `pjfnk-linear-fail` branch:

mffd-error total-nonlinear total-linear
default    7               51
1e-9       8               54
1e-8       7               51
1e-7       6               46
1e-6       6               47
1e-5       6               47
1e-4       6               47
1e-3       6               47
1e-2       8               63
1e-1       9               71


# 3/23/18

hdf5 build process:

```
wget https://www.hdfgroup.org/package/source-gzip/?wpdmdl=4301&refresh=5ab56771a3dad1521837937
CC=mpicc CXX=mpicxx FC=mpif90 ./configure --enable-fortran --enable-parallel --prefix=$HOME/hdf5/installed
make -j
make install
```

Necessary materials:
U235, U238, O16, O17, H1, H2

# 3/26/18

Trying to understand the transfer process in Okapi. These objects require the parameter `zernike_function`:
```
functions/ZernikeLegendreReconstruction
userobject/ZernikeLegendreDeconstruction
userobject/ZLDeconstruction
```
These objects require the parameter `legendre_function`:
```
functions/FourierLegendreReconstruction
functions/ZernikeLegendreReconstruction
userobject/FLDeconstruction
userobject/FourierLegendreReconstruction
userobject/ZernikeLegendreDeconstruction
userobject/ZLDeconstruction
```
These objects require the parameter `fourier_function`:
```
functions/FourierLegendreReconstruction
userobject/FLDeconstruction
userobject/FourierLegendreDeconstruction
```

# 3/29/18

Using the functional expansion module, the first set of evaluated basis functions is:
  [0] = 1
  [1] = -0.071650234569886687
  [2] = 0.5123229830564977
  [3] = -0.073416123823084295
  [4] = -0.46478280983634346
  [5] = 0.25734108285398871
  [6] = -0.056051316512444702
  [7] = 0.08577781330818092
  [8] = -0.61333986494135351
  [9] = 0.126581468737609
  [10] = -0.037785969607150048
  [11] = 0.14166122845865634
  [12] = -0.17596540952094997
  [13] = -0.49655650600439694
  [14] = 0.060834505687257026
  [15] = -0.023717427269252349
  [16] = 0.14920619572786775
  [17] = -0.036171937424850226
  [18] = 0.25864137076001681
  [19] = -0.33695442988910645
  [20] = 0.028459541840660088

and the coefficients are:
  [0] = 51690270.747073814
  [1] = -998847.10401362366
  [2] = -486590.34052204806
  [3] = 672861.49891506473
  [4] = 1959873.3578166838
  [5] = -92574.385490837026
  [6] = -756018.90552009712
  [7] = -1178478.169628079
  [8] = -121615.76601095256
  [9] = -257326.41343374335
  [10] = -904769.95741901267
  [11] = 50737.805621445528
  [12] = -800034.26501532458
  [13] = -930427.19604526995
  [14] = -230141.46631410075
  [15] = 1220594.4844509223
  [16] = 551218.29566548322
  [17] = 306692.33506197657
  [18] = -854245.07362143102
  [19] = -575518.87827935524
  [20] = -282761.30873469013

I need to generate the corresponding basis and coefficients for April's work and compare

By comparison, the numbers I get with Aprils reconstruction are:
    [0] = 0.3989422804014327
    [1] = -0.057168615941216326
    [2] = 0.40877459832524754
    [3] = -0.071742605028269568
    [4] = -0.32115948314660991
    [5] = 0.25147499899655262
    [6] = -0.063247137840919315
    [7] = 0.096789897535959568
    [8] = -0.6920799259489987
    [9] = 0.14283189226386994
    [10] = -0.04766950439593283
    [11] = 0.17871502631668182
    [12] = -0.15697206536751265
    [13] = -0.6264389346601964
    [14] = 0.076746759880269327
    [15] = -0.032776929448154406
    [16] = 0.20619947075541317
    [17] = -0.049988771021315881
    [18] = 0.35743632163532868
    [19] = -0.46566313666054238
    [20] = 0.039330420810330981

These are most certainly not the same

With Bryson's orthonormal series:
  [0] = 0.15915494309189535
  [1] = -0.045613956021965073
  [2] = 0.32615494085210783
  [3] = -0.070107234054543138
  [4] = -0.22191744494878338
  [5] = 0.24574263238099978
  [6] = -0.071366752718111603
  [7] = 0.10921570396488606
  [8] = -0.78092857040585517
  [9] = 0.16116853162738151
  [10] = -0.060138238424980534
  [11] = 0.22546084753665432
  [12] = -0.14002882369224426
  [13] = -0.79029422455040166
  [14] = 0.096821122906789763
  [15] = -0.045296949447888285
  [16] = 0.28496284308032388
  [17] = -0.069083311708507636
  [18] = 0.49396863173424338
  [19] = -0.64353555736275336
  [20] = 0.054353721144860041

With sqrt_mu and FX module:
  [0] = 0.39894228040143265
  [1] = -0.057168615941216312
  [2] = 0.40877459832524748
  [3] = -0.071742605028269568
  [4] = -0.32115948314660986
  [5] = 0.25147499899655262
  [6] = -0.063247137840919301
  [7] = 0.096789897535959554
  [8] = -0.69207992594899859
  [9] = 0.14283189226386991
  [10] = -0.047669504395932809
  [11] = 0.17871502631668179
  [12] = -0.15697206536751265
  [13] = -0.62643893466019629
  [14] = 0.0767467598802693
  [15] = -0.032776929448154392
  [16] = 0.20619947075541309
  [17] = -0.035347398972354109
  [18] = 0.35743632163532935
  [19] = -0.46566313666054227
  [20] = 0.039330420810330967

# 4/6/18

Methods needed on MooseVariableFE:

isNodal
number
nodalDofIndex
reinitNode
computeElemValues
computeElemValuesFace
computeNeighborValuesFace
computeNeighborValues
reinitNodes
computeNodalValues
reinitNodesNeighbor
computeNodalNeighborValues
prepare
prepareNeighbor
clearDofIndices
feType
sys
kind
allDofIndices
scalingFactor
dofIndices
dofIndicesNeighbor
name

# 4/9/18

Certainly for the Gill example 2 function with an initial guess of -750,000, the
Jorge method computes a larger differencing parameter `h` than the
Walker-Pernice routine. So Jorge's method definitely seems like it's doing
things. However, it performs horribly on my simple transient diffusion test
problem. Boo.

# 4/20/18

The first set of coefficients coming from patched-openmc into Okapi is:

$2 = std::__debug::vector of length 21, capacity 21 = {51,137,266.097774886, -717987.76675515505, 186120.28773150223, 196983.29861418696,
  2035817.5976047709, 608066.28785666183, -1295173.49210542, -530779.43502793019, 181706.76884885886, -926015.03226180433, -640897.5619337773,
  686638.73304699478, -340371.58757007762, -1152052.1166714183, -330127.3761021496, 1368573.0340459668, -244161.41407632269, -4945.0136016870729,
  -921542.76396300574, -785260.50748972595, 385366.29423748387}

Now I need to check and see what master openmc would produce

Ok, but here is what it is on my Mac for patched-openc:
  [0] = 51690270.747073814
  [1] = -998847.10401362366
  [2] = -486590.34052204806
  [3] = 672861.49891506473
  [4] = 1959873.3578166838
  [5] = -92574.385490837026
  [6] = -756018.90552009712
  [7] = -1178478.169628079
  [8] = -121615.76601095256
  [9] = -257326.41343374335
  [10] = -904769.95741901267
  [11] = 50737.805621445528
  [12] = -800034.26501532458
  [13] = -930427.19604526995
  [14] = -230141.46631410075
  [15] = 1220594.4844509223
  [16] = 551218.29566548322
  [17] = 306692.33506197657
  [18] = -854245.07362143102
  [19] = -575518.87827935524
  [20] = -282761.30873469013

On Mac, with master openmc with Zernike radius 1:

0     9.079856e+07
1    -4.888190e+05
2     2.938444e+04
3     2.254495e+04
4    -1.171366e+08
5    -6.171393e+04
6    -1.040070e+05
7     8.818931e+05
8    -1.840301e+04
9    -4.345738e+04
10   -5.012342e+04
11   -5.916656e+04
12    7.373350e+07
13    1.511400e+05
14   -5.045058e+04
15   -1.228643e+03
16    3.934398e+05
17   -5.371033e+05
18   -1.185331e+05
19    1.588188e+05
20    8.675914e+03

With Zernike radius 0.5:

0     9.079856e+07
1    -9.776380e+05
2     5.876888e+04
3     9.017979e+04
4     3.256926e+06
5    -2.468557e+05
6    -8.320559e+05
7    -1.240388e+06
8     3.514464e+05
9    -3.476590e+05
10   -8.019747e+05
11    1.011295e+05
12   -6.125908e+05
13   -4.499649e+05
14   -8.072093e+05
15   -3.931658e+04
16    3.613994e+05
17    1.424739e+05
18   -7.654785e+05
19   -2.732200e+04
20    2.776293e+05

# 5/24/18

Trying to figure out what the hell is going on in material properties.

We have MaterialPropertyStorage and MaterialData. There's a `nullptr` in
`material_data` corresponding to the stateful property.

There are two MaterialPropertyStorage objects owned by
`ComputeMaterialObjectThread`, `_material_props` and `_bnd_material_props`. The
`initStatefulProps` method is called on `_bnd_material_props` for both face and
neighbor materials in `ComputeMaterialsObjectThread::onInternalSide`. The
difference is we're passing in a `_bnd_material_data` in the former case, and a
`_neighbor_material_data` in the latter case.

First call from `ComputeMaterialsObjectThread::onInternalSide` through to
`initiStatefulProps` occurs on line 176. This resizes the props to `nullptr`
first and then sets it equal to an init'd value coming from `material_data`.

First call
this=0x000000010e942eb0, material_data=0x000000010e988f78,
elem=0x000000010e90c420, side=1, n_qpoints=1

Second call
this=0x000000010e942eb0, material_data=0x000000010e988f78,
elem=0x000000010e90c420, side=1, n_qpoints=1

This is fine.

Third call (second meanginful call):
this=0x000000010e942eb0, material_data=0x000000010e989018,
elem=0x000000010e9011f0, side=0, n_qpoints=1


Second meaningful call comes for the neighbor evaluation. We see we're using the
neighboring element instead of the master element, and the left side as opposed
to the right side.

So the difference is that we're calling a different material_data object (the
neighbor one). That's
when we get `EXEC_BAD_ACCESS`.

MaterialPropertyStorage::initProps(this=0x000000010f14ae90,
material_data=0x000000010f190f58, elem=0x000000010f1007d0, side=1, n_qpoints=1)

Here is the result from `material_data.props()` for object 0x000000010f190f58
(MaterialProperties) $13 = {
  std::__1::vector<PropertyValue *, std::__1::allocator<PropertyValue *> > = size=2 {
    [0] = 0x000000010f1b53d0
    [1] = 0x000000010f1bb300
  }
}

MaterialPropertyStorage::initProps(this=0x000000010f14ae90,
material_data=0x000000010f190ff8, elem=0x000000010f100860, side=0, n_qpoints=1)

Here is the result from `material_data.props()` for object 0x000000010f190ff8
(MaterialProperties) $14 = {
  std::__1::vector<PropertyValue *, std::__1::allocator<PropertyValue *> > = size=2 {
    [0] = 0x0000000000000000
    [1] = 0x000000010f1b58d0
  }
}

Neighbor materials are really NeighborFace materials.

# 6/4/18

Solve with PJFNK with wrong sign:
 0 Nonlinear |R| = 2.121320e+00
      0 Linear |R| = 2.121320e+00
      1 Linear |R| = 3.330669e-16
    |residual|_2 of individual variables:
                     disp_x: 1.81165e-16
                     disp_y: 5.2659e-16
                     lm:     1.41421

 1 Nonlinear |R| = 1.414214e+00
      0 Linear |R| = 1.414214e+00
      1 Linear |R| = 8.164966e-01
      2 Linear |R| = 1.058549e-01
      3 Linear |R| = 2.393219e-04
    |residual|_2 of individual variables:
                     disp_x: 0.00476403
                     disp_y: 0.0189226
                     lm:     0.205581

 2 Nonlinear |R| = 2.065047e-01
      0 Linear |R| = 2.065047e-01
      1 Linear |R| = 2.012658e-01
      2 Linear |R| = 2.953744e-02
      3 Linear |R| = 1.421092e-03
      4 Linear |R| = 4.484617e-04
      5 Linear |R| = 5.106881e-10
    |residual|_2 of individual variables:
                     disp_x: 0.00201239
                     disp_y: 0.00861964
                     lm:     0.033176

 3 Nonlinear |R| = 3.433648e-02
      0 Linear |R| = 3.433648e-02
      1 Linear |R| = 2.462811e-02
      2 Linear |R| = 5.488473e-03
      3 Linear |R| = 1.774935e-03
      4 Linear |R| = 1.374019e-04
      5 Linear |R| = 2.911598e-10
    |residual|_2 of individual variables:
                     disp_x: 0.000436574
                     disp_y: 0.00143109
                     lm:     0.00108736

 4 Nonlinear |R| = 1.849584e-03
      0 Linear |R| = 1.849584e-03
      1 Linear |R| = 1.336393e-03
      2 Linear |R| = 1.203622e-03
      3 Linear |R| = 4.773999e-04
      4 Linear |R| = 4.699264e-05
      5 Linear |R| = 1.557491e-10
    |residual|_2 of individual variables:
                     disp_x: 6.27536e-06
                     disp_y: 1.93043e-05
                     lm:     1.33325e-06

 5 Nonlinear |R| = 2.034246e-05
      0 Linear |R| = 2.034246e-05
      1 Linear |R| = 1.895166e-05
      2 Linear |R| = 1.663447e-05
      3 Linear |R| = 4.713097e-06
      4 Linear |R| = 3.139378e-07
      5 Linear |R| = 3.001299e-10
    |residual|_2 of individual variables:
                     disp_x: 5.14075e-10
                     disp_y: 1.54392e-09
                     lm:     2.81494e-10

 6 Nonlinear |R| = 1.651426e-09
 Solve Converged!

With -snes_fd:
 0 Nonlinear |R| = 2.121320e+00
      0 Linear |R| = 2.121320e+00
      1 Linear |R| = 2.613001e-08
    |residual|_2 of individual variables:
                     disp_x: 7.78044e-10
                     disp_y: 2.58583e-08
                     lm:     1.41421

 1 Nonlinear |R| = 1.414214e+00
      0 Linear |R| = 1.414214e+00
      1 Linear |R| = 1.020505e-01
      2 Linear |R| = 1.458798e-02
      3 Linear |R| = 9.973429e-04
    |residual|_2 of individual variables:
                     disp_x: 0.00464774
                     disp_y: 0.0188165
                     lm:     0.205624

 2 Nonlinear |R| = 2.065358e-01
      0 Linear |R| = 2.065358e-01
      1 Linear |R| = 1.879913e-01
      2 Linear |R| = 1.702897e-01
      3 Linear |R| = 6.156904e-02
      4 Linear |R| = 1.541015e-03
      5 Linear |R| = 2.374264e-09
    |residual|_2 of individual variables:
                     disp_x: 0.00201201
                     disp_y: 0.00863122
                     lm:     0.0331786

 3 Nonlinear |R| = 3.434190e-02
      0 Linear |R| = 3.434190e-02
      1 Linear |R| = 1.559931e-02
      2 Linear |R| = 2.896319e-03
      3 Linear |R| = 2.799470e-03
      4 Linear |R| = 7.727156e-04
      5 Linear |R| = 2.568599e-10
    |residual|_2 of individual variables:
                     disp_x: 0.000436757
                     disp_y: 0.00143275
                     lm:     0.00108749

 4 Nonlinear |R| = 1.850994e-03
      0 Linear |R| = 1.850994e-03
      1 Linear |R| = 9.000313e-04
      2 Linear |R| = 7.190927e-04
      3 Linear |R| = 6.105587e-04
      4 Linear |R| = 1.609641e-04
      5 Linear |R| = 1.862600e-11
    |residual|_2 of individual variables:
                     disp_x: 6.28907e-06
                     disp_y: 1.93488e-05
                     lm:     1.33384e-06

 5 Nonlinear |R| = 2.038890e-05
      0 Linear |R| = 2.038890e-05
      1 Linear |R| = 1.331898e-05
      2 Linear |R| = 1.325214e-05
      3 Linear |R| = 2.585380e-06
      4 Linear |R| = 1.684043e-13
    |residual|_2 of individual variables:
                     disp_x: 5.05942e-10
                     disp_y: 1.54793e-09
                     lm:     9.80355e-13

 6 Nonlinear |R| = 1.628518e-09
 Solve Converged!

With correct sign and perhaps correct off-diag jacobian:
 0 Nonlinear |R| = 2.121320e+00
      0 Linear |R| = 2.121320e+00
      1 Linear |R| = 3.342200e-16
    |residual|_2 of individual variables:
                     disp_x: 5.89558e-17
                     disp_y: 9.54538e-16
                     lm:     1.41421

 1 Nonlinear |R| = 1.414214e+00
      0 Linear |R| = 1.414214e+00
      1 Linear |R| = 1.330578e-01
      2 Linear |R| = 6.544887e-02
      3 Linear |R| = 8.415802e-03
      4 Linear |R| = 5.543569e-10
    |residual|_2 of individual variables:
                     disp_x: 0.00476117
                     disp_y: 0.0189006
                     lm:     0.205573

 2 Nonlinear |R| = 2.064947e-01
      0 Linear |R| = 2.064947e-01
      1 Linear |R| = 3.276384e-02
      2 Linear |R| = 2.590318e-02
      3 Linear |R| = 5.089517e-03
      4 Linear |R| = 4.455669e-03
      5 Linear |R| = 2.010387e-03
      6 Linear |R| = 2.190189e-04
      7 Linear |R| = 7.043690e-06
    |residual|_2 of individual variables:
                     disp_x: 0.00200682
                     disp_y: 0.00860163
                     lm:     0.0331745

 3 Nonlinear |R| = 3.433019e-02
      0 Linear |R| = 3.433019e-02
      1 Linear |R| = 8.430001e-03
      2 Linear |R| = 6.183527e-03
      3 Linear |R| = 1.918139e-03
      4 Linear |R| = 1.028477e-03
      5 Linear |R| = 7.967518e-04
      6 Linear |R| = 1.091485e-04
      7 Linear |R| = 1.416016e-05
    |residual|_2 of individual variables:
                     disp_x: 0.000438092
                     disp_y: 0.00142606
                     lm:     0.0010876

 4 Nonlinear |R| = 1.846193e-03
      0 Linear |R| = 1.846193e-03
      1 Linear |R| = 5.952425e-04
      2 Linear |R| = 5.796741e-04
      3 Linear |R| = 2.742078e-04
      4 Linear |R| = 1.343655e-04
      5 Linear |R| = 2.750242e-05
      6 Linear |R| = 5.150347e-06
      7 Linear |R| = 2.802799e-06
      8 Linear |R| = 3.644887e-07
    |residual|_2 of individual variables:
                     disp_x: 6.10769e-06
                     disp_y: 1.93172e-05
                     lm:     1.34181e-06

 5 Nonlinear |R| = 2.030420e-05
      0 Linear |R| = 2.030420e-05
      1 Linear |R| = 1.840838e-05
      2 Linear |R| = 4.641810e-06
      3 Linear |R| = 3.711030e-06
      4 Linear |R| = 1.727123e-06
      5 Linear |R| = 8.006059e-07
      6 Linear |R| = 2.412153e-07
      7 Linear |R| = 3.725646e-08
      8 Linear |R| = 3.455903e-09
    |residual|_2 of individual variables:
                     disp_x: 2.95903e-09
                     disp_y: 2.57325e-09
                     lm:     4.68566e-10

 6 Nonlinear |R| = 3.949308e-09
 Solve Converged!

The solve doesn't coverge at all unless full = true when using the incorrect sign in LM::computeQpJacobian

This is the jacobian difference:

row 0: (1, 0.142941)  (3, -0.142941)
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2: (1, 0.142941)  (3, -0.142941)
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (1, -0.142941)  (3, 0.142941)
row 13:
row 14: (1, -0.142941)  (3, 0.142941)
row 15:
row 16:
row 17:

And here is the difference including the LM constraint:
row 0: (1, 0.142941)  (3, -0.142941)
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2: (1, 0.142941)  (3, -0.142941)
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (1, -0.142941)  (3, 0.142941)
row 13:
row 14: (1, -0.142941)  (3, 0.142941)
row 15:
row 16:
row 17:

Clearly the difference is the same (so the LM constraint appears to be perfect)

Un-perturbed:
_lm = 0.30970422574365242
nodalArea = 0.57860651781035366
normal(_component) = 0
_test_master[_i][_qp] = 0.038459288269898577

Other _test_master[1] = 0.96154071173010114

Ok now we've perturbed things:

_test_master[_i][_qp] = 0.038459287378721223


slave disp_x as it depends on perturbations to master disp_y:

_lm = 0.30970422379987939
nodalArea = 0.57860651471059887
normal = -0.0000000080678085828944393
test_slave = 1
From node 1: -0.0000000014457261387872412
From node 2: -0.0000000014457261693219925

The un-perturbed state:
_lm = 0.30970422379987939
nodalArea = 0.57860651471059887
normal = -0.000000031239774451023539
From node 1: -0.0000000055980701611365319
From node 2: -0.0000000055980702793717053

The difference is in the normal vectors!!!!!!!!!

Swapped vars direction:
row 0: (1, 0.142941)  (3, -0.142941)
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2: (1, 0.142941)  (3, -0.142941)
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (1, -0.142941)  (3, 0.142941)
row 13: (0, 0.00689179)  (2, 0.172305)
row 14: (1, -0.142941)  (3, 0.142941)
row 15: (0, 0.172305)  (2, 0.00689178)
row 16:
row 17:


In coupled_var direction:
row 0: (1, 0.142941)  (3, -0.142941)
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2: (1, 0.142941)  (3, -0.142941)
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (1, -0.136049)  (3, 0.315246)
row 13:
row 14: (1, 0.029364)  (3, 0.149833)
row 15:
row 16:
row 17:

Switching normal direction:
row 0: (1, 0.142941)  (3, -0.142941)
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2: (1, 0.142941)  (3, -0.142941)
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (1, -0.149833)  (3, -0.029364)
row 13:
row 14: (1, -0.315246)  (3, 0.136049)
row 15:
row 16:
row 17:

-0.0068917891353286942
-0.1723050923886707

Implementing proper geometric ideas:
row 0: (1, 0.428823)  (3, -0.428823)
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2: (1, 0.428823)  (3, -0.428823)
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4: (1, 0.285882)  (3, -0.285882)
row 5:
row 6: (1, 0.285882)  (3, -0.285882)
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

First returned value:
0.14294110877336755

After correct normals implementation:
row 0:
row 1: (0, 0.142941)  (2, 0.142941)  (12, -0.142941)  (14, -0.142941)
row 2:
row 3: (0, -0.142941)  (2, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

After test perturbations:
row 0:
row 1: (0, 0.285882)  (12, -0.142941)  (14, -0.142941)
row 2:
row 3: (2, -0.285882)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:


Test perturbations #2:
row 0:
row 1: (2, 0.285882)  (4, -0.142941)  (6, -0.142941)  (12, -0.142941)  (14, -0.142941)
row 2:
row 3: (0, -0.285882)  (4, -0.142941)  (6, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

FD:
row 0: (0, 0)  (1, -0.142941)  (3, 0.142941)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (8, 0)  (9, 0)  (10, 0)  (11, 0)  (12, 0)  (13, 0)  (14, 0)  (15, 0)  (16, 0)  (17, 0)
row 1: (0, -0.142941)  (1, 0)  (2, -0.142941)  (3, 0)  (7, 0)  (12, 0.142941)  (13, 0)  (14, 0.142941)  (15, 0)  (16, 0.0222528)  (17, 0.556354)
row 2: (1, -0.142941)  (2, 0)  (3, 0.142941)  (5, 0)  (6, 0)  (7, 0)  (12, 0)  (13, 0)  (14, 0)  (15, 0)  (16, 0)  (17, 0)
row 3: (0, 0.142941)  (1, 0)  (2, 0.142941)  (3, 0)  (5, 0)  (6, 0)  (7, 0)  (8, 0)  (9, 0)  (10, 0)  (11, 0)  (12, -0.142941)  (13, 0)  (14, -0.142941)  (15, 0)  (16, 0.556354)  (17, 0.0222528)
row 4: (4, 0)  (5, 0)  (6, 0)  (7, 0)  (8, 0)  (9, 0)  (10, 0)  (11, 0)  (13, 0)  (15, 0)
row 5: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (8, 0)  (9, 0)  (10, 0)  (11, 0)  (12, 0)  (13, 0)  (14, 0)  (15, 0)  (16, 0)  (17, 0)
row 6: (4, 0)  (5, 0)  (6, 0)  (7, 0)  (8, 0)  (9, 0)  (10, 0)  (11, 0)  (13, 0)  (15, 0)
row 7: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (8, 0)  (9, 0)  (10, 0)  (11, 0)  (12, 0)  (13, 0)  (14, 0)  (15, 0)  (16, 0)  (17, 0)
row 8: (8, 0.)
row 9: (9, 0.)
row 10: (10, 0.)
row 11: (11, 0.)
row 12: (0, 0)  (1, 0.142941)  (2, 0)  (3, -0.142941)  (12, 0.)  (16, 0)
row 13: (1, 0)  (3, 0)  (13, 0.)  (16, -0.578607)
row 14: (0, 0)  (1, 0.142941)  (2, 0)  (3, -0.142941)  (14, 0.)  (17, 0)
row 15: (1, 0)  (3, 0)  (15, 0.)  (17, -0.578607)
row 16: (16, 0.)
row 17: (17, 0.)

Solve with normals part:
 0 Nonlinear |R| = 2.121320e+00
      0 Linear |R| = 2.121320e+00
      1 Linear |R| = 3.342200e-16
    |residual|_2 of individual variables:
                     disp_x: 5.89558e-17
                     disp_y: 9.54538e-16
                     lm:     1.41421

 1 Nonlinear |R| = 1.414214e+00
      0 Linear |R| = 1.414214e+00
      1 Linear |R| = 1.062445e-01
      2 Linear |R| = 1.223260e-09
    |residual|_2 of individual variables:
                     disp_x: 0.00476117
                     disp_y: 0.0189006
                     lm:     0.205573

 2 Nonlinear |R| = 2.064947e-01
      0 Linear |R| = 2.064947e-01
      1 Linear |R| = 7.324059e-02
      2 Linear |R| = 1.705770e-02
      3 Linear |R| = 1.010074e-03
      4 Linear |R| = 8.385945e-05
    |residual|_2 of individual variables:
                     disp_x: 0.00199506
                     disp_y: 0.0086059
                     lm:     0.0331749

 3 Nonlinear |R| = 3.433094e-02
      0 Linear |R| = 3.433094e-02
      1 Linear |R| = 9.273598e-03
      2 Linear |R| = 6.048436e-03
      3 Linear |R| = 1.244520e-03
      4 Linear |R| = 2.125472e-04
      5 Linear |R| = 9.697987e-05
      6 Linear |R| = 1.465704e-05
    |residual|_2 of individual variables:
                     disp_x: 0.000436559
                     disp_y: 0.00142764
                     lm:     0.00108736

 4 Nonlinear |R| = 1.846911e-03
      0 Linear |R| = 1.846911e-03
      1 Linear |R| = 1.133112e-03
      2 Linear |R| = 5.517015e-04
      3 Linear |R| = 1.958470e-04
      4 Linear |R| = 3.341504e-05
      5 Linear |R| = 9.795720e-06
      6 Linear |R| = 1.568954e-06
    |residual|_2 of individual variables:
                     disp_x: 6.13252e-06
                     disp_y: 1.926e-05
                     lm:     1.37969e-06

 5 Nonlinear |R| = 2.025979e-05
      0 Linear |R| = 2.025979e-05
      1 Linear |R| = 1.511227e-05
      2 Linear |R| = 7.079576e-06
      3 Linear |R| = 1.854992e-06
      4 Linear |R| = 5.078927e-07
      5 Linear |R| = 9.815806e-08
      6 Linear |R| = 2.262294e-08
      7 Linear |R| = 4.314503e-09
    |residual|_2 of individual variables:
                     disp_x: 2.09696e-09
                     disp_y: 4.11604e-09
                     lm:     8.13139e-10

 6 Nonlinear |R| = 4.690443e-09
 Solve Converged!

For rotated problem, here's the finite difference jacobian:

row 0: (0, -0.142941)  (1, 0)  (2, 0)  (3, 0.142941)  (5, 0)  (7, 0)  (12, 0.0714706)  (13, -0.0714706)  (14, 0.0714706)  (15, -0.0714706)  (16, 0.0157351)  (17, 0.393401)
row 1: (0, 0)  (1, 0.142941)  (2, -0.142941)  (3, 0)  (4, 0)  (7, 0)  (12, 0.0714706)  (13, -0.0714706)  (14, 0.0714706)  (15, -0.0714706)  (16, 0.0157351)  (17, 0.393401)
row 2: (0, 0)  (1, -0.142941)  (2, 0.142941)  (3, 0)  (5, 0)  (7, 0)  (12, -0.0714706)  (13, 0.0714706)  (14, -0.0714706)  (15, 0.0714706)  (16, 0.393401)  (17, 0.0157351)
row 3: (0, 0.142941)  (1, 0)  (2, 0)  (3, -0.142941)  (5, 0)  (7, 0)  (12, -0.0714706)  (13, 0.0714706)  (14, -0.0714706)  (15, 0.0714706)  (16, 0.393401)  (17, 0.0157351)
row 4: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 5: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 6: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 7: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 8: (8, 0.)
row 9: (9, 0.)
row 10: (10, 0.)
row 11: (11, 0.)
row 12: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (12, 0.)  (16, -0.409137)
row 13: (0, -0.0714706)  (1, -0.0714706)  (2, 0.0714706)  (3, 0.0714706)  (13, 0.)  (16, -0.409137)
row 14: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (14, 0.)  (17, -0.409137)
row 15: (0, -0.0714706)  (1, -0.0714706)  (2, 0.0714706)  (3, 0.0714706)  (15, 0.)  (17, -0.409137)
row 16: (0, -0.0271948)  (1, -0.0271948)  (2, -0.679912)  (3, -0.679912)  (12, 0.707107)  (13, 0.707107)  (16, 0.)
row 17: (0, -0.679912)  (1, -0.679912)  (2, -0.0271948)  (3, -0.0271948)  (14, 0.707107)  (15, 0.707107)  (17, 0.)

Here's the current hand coded:
row 0: (0, -0.172545)  (1, 0.0296041)  (2, -0.0296041)  (3, 0.172545)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.101075)  (13, -0.101075)  (14, 0.101075)  (15, -0.101075)  (16, 0.0157351)  (17, 0.393401)
row 1: (0, -0.172545)  (1, 0.0296041)  (2, -0.0296041)  (3, 0.172545)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.101075)  (13, -0.101075)  (14, 0.101075)  (15, -0.101075)  (16, 0.0157351)  (17, 0.393401)
row 2: (0, 0.0296041)  (1, -0.172545)  (2, 0.172545)  (3, -0.0296041)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -0.101075)  (13, 0.101075)  (14, -0.101075)  (15, 0.101075)  (16, 0.393401)  (17, 0.0157351)
row 3: (0, 0.0296041)  (1, -0.172545)  (2, 0.172545)  (3, -0.0296041)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -0.101075)  (13, 0.101075)  (14, -0.101075)  (15, 0.101075)  (16, 0.393401)  (17, 0.0157351)
row 4: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 5: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 6: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 7: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 8: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 9: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 10: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 11: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 12: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, -0.409137)  (17, 0.)
row 13: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, -0.409137)  (17, 0.)
row 14: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, -0.409137)
row 15: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, -0.409137)
row 16: (0, -0.0543896)  (1, -0.0271948)  (2, -1.35982)  (3, -0.679912)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (12, 1.41421)  (13, 0.707107)  (16, 0.)
row 17: (0, -1.35982)  (1, -0.679912)  (2, -0.0543897)  (3, -0.0271948)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (14, 1.41421)  (15, 0.707106)  (17, 0.)

And difference:
row 0: (0, -0.0296041)  (1, 0.0296041)  (2, -0.0296041)  (3, 0.0296041)  (12, 0.0296041)  (13, -0.0296041)  (14, 0.0296041)  (15, -0.0296041)
row 1: (0, -0.172545)  (1, -0.113337)  (2, 0.113337)  (3, 0.172545)  (12, 0.0296041)  (13, -0.0296041)  (14, 0.0296041)  (15, -0.0296041)
row 2: (0, 0.0296041)  (1, -0.0296041)  (2, 0.0296041)  (3, -0.0296041)  (12, -0.0296041)  (13, 0.0296041)  (14, -0.0296041)  (15, 0.0296041)
row 3: (0, -0.113337)  (1, -0.172545)  (2, 0.172545)  (3, 0.113337)  (12, -0.0296041)  (13, 0.0296041)  (14, -0.0296041)  (15, 0.0296041)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13: (0, 0.142941)  (1, 0.142941)  (2, -0.142941)  (3, -0.142941)
row 14:
row 15: (0, 0.142941)  (1, 0.142941)  (2, -0.142941)  (3, -0.142941)
row 16: (0, -0.0271948)  (2, -0.679912)  (12, 0.707107)
row 17: (0, -0.679913)  (2, -0.0271949)  (14, 0.707108)


For (0, 0) Jacobian entry:
Node 1
normals: -0.0027487059762367604
perturbations: -0.0038872572574903683

Node 2
normals: -0.068721848891130505
perturbations: -0.0971873704061803

For (0, 1) Jacobian entry:
Node 1
normals: -0.0027487059577393615

after intro of _normal(comp) in testPerturbations:
Finite-difference:
row 0: (0, -0.142941)  (1, 0)  (2, 0)  (3, 0.142941)  (5, 0)  (7, 0)  (12, 0.0714706)  (13, -0.0714706)  (14, 0.0714706)  (15, -0.0714706)  (16, 0.0157351)  (17, 0.393401)
row 1: (0, 0)  (1, 0.142941)  (2, -0.142941)  (3, 0)  (4, 0)  (7, 0)  (12, 0.0714706)  (13, -0.0714706)  (14, 0.0714706)  (15, -0.0714706)  (16, 0.0157351)  (17, 0.393401)
row 2: (0, 0)  (1, -0.142941)  (2, 0.142941)  (3, 0)  (5, 0)  (7, 0)  (12, -0.0714706)  (13, 0.0714706)  (14, -0.0714706)  (15, 0.0714706)  (16, 0.393401)  (17, 0.0157351)
row 3: (0, 0.142941)  (1, 0)  (2, 0)  (3, -0.142941)  (5, 0)  (7, 0)  (12, -0.0714706)  (13, 0.0714706)  (14, -0.0714706)  (15, 0.0714706)  (16, 0.393401)  (17, 0.0157351)
row 4: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 5: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 6: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 7: (0, 0)  (1, 0)  (2, 0)  (3, 0)  (4, 0)  (5, 0)  (6, 0)  (7, 0)  (13, 0)  (14, 0)
row 8: (8, 0.)
row 9: (9, 0.)
row 10: (10, 0.)
row 11: (11, 0.)
row 12: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (12, 0.)  (16, -0.409137)
row 13: (0, -0.0714706)  (1, -0.0714706)  (2, 0.0714706)  (3, 0.0714706)  (13, 0.)  (16, -0.409137)
row 14: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (14, 0.)  (17, -0.409137)
row 15: (0, -0.0714706)  (1, -0.0714706)  (2, 0.0714706)  (3, 0.0714706)  (15, 0.)  (17, -0.409137)
row 16: (0, -0.0271948)  (1, -0.0271948)  (2, -0.679912)  (3, -0.679912)  (12, 0.707107)  (13, 0.707107)  (16, 0.)
row 17: (0, -0.679912)  (1, -0.679912)  (2, -0.0271948)  (3, -0.0271948)  (14, 0.707107)  (15, 0.707107)  (17, 0.)
Hand-coded:
row 0: (0, -0.142941)  (1, 0)  (2, 0)  (3, 0.142941)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.0714706)  (13, -0.0714706)  (14, 0.0714706)  (15, -0.0714706)  (16, 0.0157351)  (17, 0.393401)
row 1: (0, -0.142941)  (1, -0)  (2, 0)  (3, 0.142941)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.0714706)  (13, -0.0714706)  (14, 0.0714706)  (15, -0.0714706)  (16, 0.0157351)  (17, 0.393401)
row 2: (0, -0)  (1, -0.142941)  (2, 0.142941)  (3, -0)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -0.0714706)  (13, 0.0714706)  (14, -0.0714706)  (15, 0.0714706)  (16, 0.393401)  (17, 0.0157351)
row 3: (0, -0)  (1, -0.142941)  (2, 0.142941)  (3, 0)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -0.0714706)  (13, 0.0714706)  (14, -0.0714706)  (15, 0.0714706)  (16, 0.393401)  (17, 0.0157351)
row 4: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 5: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 6: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 7: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 8: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 9: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 10: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 11: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)
row 12: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, -0.409137)  (17, 0.)
row 13: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, -0.409137)  (17, 0.)
row 14: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, -0.409137)
row 15: (0, 0.0714706)  (1, 0.0714706)  (2, -0.0714706)  (3, -0.0714706)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, -0.409137)
row 16: (0, -0.0543896)  (1, -0.0271948)  (2, -1.35982)  (3, -0.679912)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (12, 1.41421)  (13, 0.707107)  (16, 0.)
row 17: (0, -1.35982)  (1, -0.679912)  (2, -0.0543897)  (3, -0.0271948)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (14, 1.41421)  (15, 0.707106)  (17, 0.)
Difference:
row 0:
row 1: (0, -0.142941)  (1, -0.142941)  (2, 0.142941)  (3, 0.142941)
row 2:
row 3: (0, -0.142941)  (1, -0.142941)  (2, 0.142941)  (3, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13: (0, 0.142941)  (1, 0.142941)  (2, -0.142941)  (3, -0.142941)
row 14:
row 15: (0, 0.142941)  (1, 0.142941)  (2, -0.142941)  (3, -0.142941)
row 16: (0, -0.0271948)  (2, -0.679912)  (12, 0.707107)
row 17: (0, -0.679913)  (2, -0.0271949)  (14, 0.707108)

Rotation into second quadrant:
Difference:
row 0: (0, -0.142941)  (1, -0.142941)  (2, -0.142941)  (3, -0.142941)  (12, 0.142941)  (13, 0.142941)  (14, 0.142941)  (15, 0.142941)
row 1:
row 2: (0, 0.142941)  (1, 0.142941)  (2, 0.142941)  (3, 0.142941)  (12, -0.142941)  (13, -0.142941)  (14, -0.142941)  (15, -0.142941)
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

After geometric thinking:
row 0: (12, 0.142941)  (14, 0.142941)
row 1: (0, 0.142941)  (1, 0.142941)  (2, 0.142941)  (3, 0.142941)  (13, -0.142941)  (15, -0.142941)
row 2: (13, -0.142941)  (15, -0.142941)
row 3: (0, -0.142941)  (1, -0.142941)  (2, -0.142941)  (3, -0.142941)  (12, 0.142941)  (14, 0.142941)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

Un-perturbed shape function: 0.03845932241566441
0.038459322415664354
Perturbed shape function:
0.038459321524486216

Perturbed system:
first resid:
-0.0068917908025728199
second resid:
-0.17230509652570158
total resid:
-0.1791968873282744

Base system:
first resid:
-0.0068917909622691802
second resid:
-0.17230510051834935
total resid:
-0.17919689148061854

finite difference jacobian:
row 0: (0, 3.82185e-09)  (1, 0.142941)  (2, -4.77732e-09)  (3, 0.142941)  (5, 2.86639e-09)  (6, 1.91093e-09)  (7, -9.55463e-10)  (8, 1.91093e-09)  (9, 1.91093e-09)  (10, 1.91093e-09)  (11, 1.91093e-09)  (12, 9.55463e-10)  (13, -0.142941)  (14, 9.55463e-10)  (15, -0.142941)  (16, -0.0222528)  (17, -0.556354)
row 1: (0, 0.142941)  (1, 2.97842e-10)  (2, -0.142941)  (5, 1.77969e-18)  (6, 1.77969e-18)  (7, -8.89845e-19)  (8, 1.77969e-18)  (9, 1.77969e-18)  (10, 1.77969e-18)  (11, 1.77969e-18)  (13, -1.48921e-10)  (15, -1.48921e-10)  (16, -2.31838e-11)  (17, -5.79629e-10)
row 2: (0, -2.86639e-09)  (1, -0.142941)  (2, 3.82185e-09)  (3, -0.142941)  (4, -9.55463e-10)  (5, -1.91093e-09)  (6, -2.86639e-09)  (7, -9.55463e-10)  (8, -1.91093e-09)  (9, -1.91093e-09)  (10, -1.91093e-09)  (11, -1.91093e-09)  (13, 0.142941)  (15, 0.142941)  (16, -0.556354)  (17, -0.0222528)
row 3: (0, 0.142941)  (2, -0.142941)  (3, -2.97842e-10)  (4, -8.89845e-19)  (5, -1.77969e-18)  (6, -2.66953e-18)  (7, -8.89845e-19)  (8, -1.77969e-18)  (9, -1.77969e-18)  (10, -1.77969e-18)  (11, -1.77969e-18)  (13, 1.48921e-10)  (15, 1.48921e-10)  (16, -5.79629e-10)  (17, -2.31837e-11)
row 4: (0, 6.58484e-11)  (1, -1.6068e-09)  (2, 1.18527e-10)  (3, -2.26533e-09)  (4, -1.2512e-09)  (5, -9.87787e-10)  (6, 1.31697e-11)  (7, 3.9509e-11)  (12, -1.64631e-09)  (13, -1.31705e-09)  (14, -1.64631e-09)  (15, -1.58046e-09)  (16, -1.64631e-09)  (17, -1.64631e-09)
row 5: (4, -1.30355e-18)  (5, -1.02911e-18)  (6, 1.37206e-20)  (7, 4.11619e-20)
row 6: (0, 1.64631e-09)  (1, 9.21939e-10)  (2, 2.96336e-09)  (3, 8.95599e-10)  (4, 1.59363e-09)  (5, -3.95091e-11)  (6, 3.29262e-10)  (7, 9.87787e-10)  (12, -6.58485e-11)  (13, -5.26788e-11)  (14, -6.58485e-11)  (15, 1.58046e-09)  (16, -6.58485e-11)  (17, -6.58486e-11)
row 7: (4, 1.6603e-18)  (5, -4.1162e-20)  (6, 3.43037e-19)  (7, 1.02911e-18)
row 8: (8, 0.)
row 9: (9, 0.)
row 10: (10, 0.)
row 11: (11, 0.)
row 12: (0, -9.55463e-10)  (2, -1.91093e-09)  (12, 0.)  (16, 0.578607)
row 13: (0, -0.142941)  (1, -1.48921e-10)  (2, 0.142941)  (3, 1.48921e-10)  (13, 0.)  (16, 6.02812e-10)
row 14: (0, -9.55463e-10)  (2, -1.91093e-09)  (14, 0.)  (17, 0.578607)
row 15: (0, -0.142941)  (1, -1.48921e-10)  (2, 0.142941)  (3, 1.48921e-10)  (15, 0.)  (17, 6.02812e-10)
row 16: (0, 0.0384593)  (1, 3.82185e-09)  (2, 0.961541)  (3, 3.82185e-09)  (12, -1.)  (14, 3.82185e-09)  (15, 3.82185e-09)  (16, 3.82185e-09)  (17, 3.82185e-09)
row 17: (0, 0.961541)  (2, 0.0384592)  (14, -1.)  (15, -3.82185e-09)  (17, 0.)

Here is the error at 315 degree rotation:
row 0:
row 1:
row 2:
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16: (2, -1.06608e-05)  (3, 1.07029e-05)  (12, 1.11425e-05)  (13, -1.10976e-05)
row 17:

And the error at 45 degree rotation:
row 0:
row 1:
row 2:
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16: (2, -1.12856e-05)  (3, -1.12644e-05)  (12, 1.17731e-05)  (13, 1.17702e-05)
row 17:

The error is not observed at any other rotation angle. I'm inclined to write this off as actually an error in the finite difference jacobian.

# 6/7/18

Important notes: The current lagrange multiplier formulation (e.g. residual statements) seems to really challenge a direct solver. Running with the `frictionless_lagrange` sliding block problem with `-snes_type fd` and NEWTON, `-pc_type lu` fails to solve. However, `-pc_type ilu` works very well. Additionally, line searches also do not work with high youngs modulus values. Ok this is weird...the problem converges marvelously with a Youngs modulus of 1e8. It also converges fairly well with Youngs modulus of 1e12. With `-snes_type fd -snes_mf_operator` the solve fails. So it looks like PJFNK with the lagrange multiplier formulation *may* struggle.

# 6/13/18

two-body-gravity.i:

  Testing hand-coded Jacobian, if (for double precision runs) ||J - Jfd||_F/||J||_F is
    O(1.e-8), the hand-coded Jacobian is probably correct.
  ||J - Jfd||_F/||J||_F = 0.725346, ||J - Jfd||_F = 2.70292
  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 0.669697)  (1, 0.)  (2, -0.119697)  (3, 0.)  (4, -0.334848)  (5, 0.)  (6, -0.215152)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 1: (0, 0.)  (1, 0.669697)  (2, 0.)  (3, -0.119697)  (4, 0.)  (5, -0.334848)  (6, 0.)  (7, -0.215152)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.0227273)  (17, 0.477273)
row 2: (0, -0.119697)  (1, 0.)  (2, 0.669697)  (3, 0.)  (4, -0.215152)  (5, 0.)  (6, -0.334848)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 3: (0, 0.)  (1, -0.119697)  (2, 0.)  (3, 0.669697)  (4, 0.)  (5, -0.215152)  (6, 0.)  (7, -0.334848)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.477273)  (17, 0.0227273)
row 4: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 1.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 5: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 1.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 6: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 1.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 7: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 1.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 8: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 1.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 9: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 1.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 10: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 1.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 11: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 1.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 12: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.333333)  (9, 0.)  (10, -0.166667)  (11, 0.)  (12, 0.666667)  (13, 0.)  (14, -0.166667)  (15, 0.)  (16, 0.)  (17, 0.)
row 13: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.333333)  (10, 0.)  (11, -0.166667)  (12, 0.)  (13, 0.666667)  (14, 0.)  (15, -0.166667)  (16, -0.5)  (17, 0.)
row 14: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.166667)  (9, 0.)  (10, -0.333333)  (11, 0.)  (12, -0.166667)  (13, 0.)  (14, 0.666667)  (15, 0.)  (16, 0.)  (17, 0.)
row 15: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.166667)  (10, 0.)  (11, -0.333333)  (12, 0.)  (13, -0.166667)  (14, 0.)  (15, 0.666667)  (16, 0.)  (17, -0.5)
row 16: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (12, 0.)  (13, 0.)  (16, 0.)
row 17: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (14, 0.)  (15, 0.)  (17, 0.)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 0.669697)  (2, -0.119697)  (4, -0.334848)  (6, -0.215152)
row 1: (1, 0.669697)  (3, -0.119697)  (5, -0.334848)  (7, -0.215152)  (16, 0.0227273)  (17, 0.477273)
row 2: (0, -0.119697)  (2, 0.669697)  (4, -0.215152)  (6, -0.334848)
row 3: (1, -0.119697)  (3, 0.669697)  (5, -0.215152)  (7, -0.334848)  (16, 0.477273)  (17, 0.0227273)
row 4: (4, 1.)
row 5: (5, 1.)
row 6: (6, 1.)
row 7: (7, 1.)
row 8: (8, 1.)
row 9: (9, 1.)
row 10: (10, 1.)
row 11: (11, 1.)
row 12: (8, -0.333333)  (10, -0.166667)  (12, 0.666667)  (14, -0.166667)
row 13: (9, -0.333333)  (11, -0.166667)  (13, 0.666667)  (15, -0.166667)  (16, -0.5)
row 14: (8, -0.166667)  (10, -0.333333)  (12, -0.166667)  (14, 0.666667)
row 15: (9, -0.166667)  (11, -0.333333)  (13, -0.166667)  (15, 0.666667)  (17, -0.5)
row 16: (1, -0.0909091)  (3, -1.90909)  (16, 0.)
row 17: (0, -8.60319e-09)  (1, -1.90909)  (3, -0.0909091)  (17, 0.)
  Hand-coded minus finite-difference Jacobian with tolerance 1e-05 ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0:
row 1:
row 2:
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16: (1, 0.0909091)  (3, 1.90909)
row 17: (1, 1.90909)  (3, 0.0909091)

  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 0.669697)  (2, -0.119697)  (4, -0.334848)  (6, -0.215152)  (17, 0.)
row 1: (1, 0.669697)  (3, -0.119697)  (5, -0.334848)  (7, -0.215152)  (16, 0.0227273)  (17, 0.477273)
row 2: (0, -0.119697)  (2, 0.669697)  (4, -0.215152)  (6, -0.334848)  (17, 0.)
row 3: (1, -0.119697)  (3, 0.669697)  (5, -0.215152)  (7, -0.334848)  (16, 0.477273)  (17, 0.0227273)
row 4: (4, 1.)  (17, 0.)
row 5: (5, 1.)  (17, 0.)
row 6: (6, 1.)  (17, 0.)
row 7: (7, 1.)  (17, 0.)
row 8: (8, 1.)  (15, 0.)
row 9: (9, 1.)  (15, 0.)
row 10: (10, 1.)  (15, 0.)
row 11: (11, 1.)  (15, 0.)
row 12: (8, -0.333333)  (10, -0.166667)  (12, 0.666667)  (14, -0.166667)  (17, 0.)
row 13: (9, -0.333333)  (11, -0.166667)  (13, 0.666667)  (15, -0.166667)  (16, -0.5)  (17, 0.)
row 14: (8, -0.166667)  (10, -0.333333)  (12, -0.166667)  (14, 0.666667)  (17, 0.)
row 15: (9, -0.166667)  (11, -0.333333)  (13, -0.166667)  (15, 0.666667)  (17, -0.5)
row 16: (16, 0.)
row 17: (17, 0.)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 0.669697)  (2, -0.119697)  (4, -0.334848)  (6, -0.215152)
row 1: (1, 0.669697)  (3, -0.119697)  (5, -0.334848)  (7, -0.215152)  (16, 0.0227273)  (17, 0.477273)
row 2: (0, -0.119697)  (2, 0.669697)  (4, -0.215152)  (6, -0.334848)
row 3: (1, -0.119697)  (3, 0.669697)  (5, -0.215152)  (7, -0.334848)  (16, 0.477273)  (17, 0.0227273)
row 4: (4, 1.)
row 5: (5, 1.)
row 6: (6, 1.)
row 7: (7, 1.)
row 8: (8, 1.)
row 9: (9, 1.)
row 10: (10, 1.)
row 11: (11, 1.)
row 12: (8, -0.333333)  (10, -0.166667)  (12, 0.666667)  (14, -0.166667)
row 13: (9, -0.333333)  (11, -0.166667)  (13, 0.666667)  (15, -0.166667)  (16, -0.5)
row 14: (8, -0.166667)  (10, -0.333333)  (12, -0.166667)  (14, 0.666667)
row 15: (9, -0.166667)  (11, -0.333333)  (13, -0.166667)  (15, 0.666667)  (17, -0.5)
row 16: (1, -0.0909091)  (3, -1.90909)  (16, 0.)
row 17: (0, -8.60319e-09)  (1, -1.90909)  (3, -0.0909091)  (17, 0.)
  Hand-coded minus finite-difference Jacobian with tolerance 1e-05 ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0:
row 1:
row 2:
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16: (1, 0.0909091)  (3, 1.90909)
row 17: (1, 1.90909)  (3, 0.0909091)

# 6/14/18

 2 Nonlinear |R| = 1.814029e-01
  ---------- Testing Jacobian -------------
  ||J - Jfd||_F/||J||_F = 0.00757875, ||J - Jfd||_F = 0.0355999
  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 0.669697)  (1, -0.0943396)  (2, -0.119697)  (3, 0.0943396)  (4, -0.334848)  (5, 0.)  (6, -0.215152)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -9.52164e-18)  (13, -9.61013e-34)  (14, -9.52164e-18)  (15, -9.61013e-34)  (16, -2.29385e-18)  (17, -4.81708e-17)
row 1: (0, -0.0943396)  (1, 0.669697)  (2, -0.0943396)  (3, -0.119697)  (4, 0.)  (5, -0.334848)  (6, 0.)  (7, -0.215152)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.0943396)  (13, 9.52164e-18)  (14, 0.0943396)  (15, 9.52164e-18)  (16, 0.0227273)  (17, 0.477273)
row 2: (0, -0.119697)  (1, -0.0943396)  (2, 0.669697)  (3, 0.0943396)  (4, -0.215152)  (5, 0.)  (6, -0.334848)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 9.52164e-18)  (13, 9.61013e-34)  (14, 9.52164e-18)  (15, 9.61013e-34)  (16, -4.81708e-17)  (17, -2.29385e-18)
row 3: (0, 0.0943396)  (1, -0.119697)  (2, 0.0943396)  (3, 0.669697)  (4, 0.)  (5, -0.215152)  (6, 0.)  (7, -0.334848)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -0.0943396)  (13, -9.52164e-18)  (14, -0.0943396)  (15, -9.52164e-18)  (16, 0.477273)  (17, 0.0227273)
row 4: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 1.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 5: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 1.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 6: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 1.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 7: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 1.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)
row 8: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 1.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 9: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 1.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 10: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 1.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 11: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 1.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)
row 12: (0, -9.52164e-18)  (1, 0.0943396)  (2, 9.52164e-18)  (3, -0.0943396)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.333333)  (9, 0.)  (10, -0.166667)  (11, 0.)  (12, 0.666667)  (13, 0.)  (14, -0.166667)  (15, 0.)  (16, 5.04647e-17)  (17, 0.)
row 13: (0, -1.15212e-17)  (1, 9.52164e-18)  (2, 1.15212e-17)  (3, -9.52164e-18)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.333333)  (10, 0.)  (11, -0.166667)  (12, 0.)  (13, 0.666667)  (14, 0.)  (15, -0.166667)  (16, -0.5)  (17, 0.)
row 14: (0, -9.52164e-18)  (1, 0.0943396)  (2, 9.52164e-18)  (3, -0.0943396)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.166667)  (9, 0.)  (10, -0.333333)  (11, 0.)  (12, -0.166667)  (13, 0.)  (14, 0.666667)  (15, 0.)  (16, 0.)  (17, 5.04647e-17)
row 15: (0, -1.15212e-17)  (1, 9.52164e-18)  (2, 1.15212e-17)  (3, -9.52164e-18)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.166667)  (10, 0.)  (11, -0.333333)  (12, 0.)  (13, -0.166667)  (14, 0.)  (15, 0.666667)  (16, 0.)  (17, -0.5)
row 16: (0, 6.63938e-18)  (1, -0.0657824)  (2, 1.39427e-16)  (3, -1.38143)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (12, -1.46066e-16)  (13, 1.44721)  (16, 0.105573)
row 17: (0, 1.39427e-16)  (1, -1.38143)  (2, 6.63938e-18)  (3, -0.0657824)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (14, -1.46066e-16)  (15, 1.44721)  (17, 0.105573)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 0.669697)  (1, -0.0943396)  (2, -0.119697)  (3, 0.0943396)  (4, -0.334848)  (6, -0.215152)
row 1: (0, -0.0943396)  (1, 0.651897)  (2, -0.0943396)  (3, -0.101897)  (4, -9.52002e-10)  (5, -0.334848)  (6, -9.52002e-10)  (7, -0.215152)  (12, 0.0943396)  (14, 0.0943396)  (16, 0.0227273)  (17, 0.477273)
row 2: (0, -0.119697)  (1, -0.0943396)  (2, 0.669697)  (3, 0.0943396)  (4, -0.215152)  (6, -0.334848)
row 3: (0, 0.0943396)  (1, -0.101897)  (2, 0.0943396)  (3, 0.651897)  (4, 3.86751e-10)  (5, -0.215152)  (6, 8.62752e-10)  (7, -0.334848)  (12, -0.0943396)  (14, -0.0943396)  (16, 0.477273)  (17, 0.0227273)
row 4: (4, 1.)
row 5: (5, 1.)
row 6: (6, 1.)
row 7: (7, 1.)
row 8: (8, 1.)
row 9: (9, 1.)
row 10: (10, 1.)
row 11: (11, 1.)
row 12: (1, 0.0943396)  (3, -0.0943396)  (8, -0.333333)  (10, -0.166667)  (12, 0.666667)  (14, -0.166667)
row 13: (1, 1.428e-09)  (3, 1.428e-09)  (9, -0.333333)  (11, -0.166667)  (13, 0.666667)  (15, -0.166667)  (16, -0.5)
row 14: (1, 0.0943396)  (3, -0.0943396)  (8, -0.166667)  (10, -0.333333)  (12, -0.166667)  (14, 0.666667)
row 15: (1, 1.428e-09)  (3, 1.428e-09)  (9, -0.166667)  (11, -0.333333)  (13, -0.166667)  (15, 0.666667)  (17, -0.5)
row 16: (0, -4.76001e-09)  (1, -0.0657824)  (3, -1.38143)  (13, 1.44721)  (16, 0.105573)
row 17: (1, -1.38143)  (3, -0.0657824)  (15, 1.44721)  (17, 0.105573)
  Hand-coded minus finite-difference Jacobian with tolerance 1e-05 ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0:
row 1: (1, 0.0177999)  (3, -0.0177999)
row 2:
row 3: (1, -0.0177999)  (3, 0.0177999)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

At second non-linear iteration:

Based:
lm = 0.20754716981132079
nodal_area = 0.5
normal(component) = -1
test_master = 0.045454545454545581
pinfo->_distance = 0.10377358490566035

Perturbed:
lm = 0.20754716981132079
nodal_area = 0.5
normal(component) = 0.99999999999999955
test_master = 0.045454542954120096
pinfo->_distance = 0.10377358623088576

You can see that the test function is indeed perturbed

At 5th non-linear iteration:

Base:
test_master = 0.04545454545454547
pinfo->_distance = 0.00000079725455170986947

Perturbed:
test_master = 0.045454545454526318
pinfo->_distance = 0.00000079857416701578664

lm = 0.26190434429523474
nodal_area = 0.5
normal = 1
wscale = 34445299.412807018

So the difference in the Jacobian I believe is correlated with the distance

Problems start with 2nd non-linear iteration:
  Hand-coded minus finite-difference Jacobian with tolerance 1e-05 ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0:
row 1: (12, -0.119047)  (14, -0.119047)
row 2:
row 3: (12, 0.119047)  (14, 0.119047)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17:

Base:
lm = 0.20754716981132079
nodal_area = 0.5
normal(component) = -1
test_master = 0.045454545454545581

Perturbed:
lm = 0.20754716981132079
nodal_area = 0.5
normal(component) = 0.99999999999999955
test_master = 0.045454571959054939

Important notes: The current lagrange multiplier formulation (e.g. residual statements) seems to really challenge a direct solver. Running with the `frictionless_lagrange` sliding block problem with `-snes_type fd` and NEWTON, `-pc_type lu` fails to solve. However, `-pc_type ilu` works very well. Additionally, line searches also do not work with high youngs modulus values. Ok this is weird...the problem converges marvelously with a Youngs modulus of 1e8. It also converges fairly well with Youngs modulus of 1e12. With `-snes_type fd -snes_mf_operator` the solve fails. So it looks like PJFNK with the lagrange multiplier formulation *may* struggle. `-pc_type lu -pc_factor_mat_solver_package mumps` actually works very well for a youngs modulus of 1e5 but fails for 1e6.

(MooseVariableFE<double>) $15 = {
  MooseVariableFEBase = {
    MooseVariableBase = {
      _var_num = 0
      _fe_type = {
        order = (_order = 2)
        family = LAGRANGE
      }
      _index = 0
      _var_kind = VAR_NONLINEAR
      _subproblem = 0x0000000113848218
      _sys = 0x0000000113876c18
      _variable = 0x00000001141ef920
      _dof_map = 0x0000000114142460
      _dof_indices = size=0 {}
      _mesh = 0x0000000113800018
      _scaling_factor = 1
    }
  }

# 6/15/18

Errors:
nx = 20
+----------------+----------------+----------------+----------------+----------------+
| time           | my_u1          | my_u2          | yaqi_u1        | yaqi_u2        |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   4.992761e-02 |   5.548607e-02 |   2.144726e-03 |   8.211431e-03 |
+----------------+----------------+----------------+----------------+----------------+

nx = 40
+----------------+----------------+----------------+----------------+----------------+
| time           | my_u1          | my_u2          | yaqi_u1        | yaqi_u2        |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   2.450870e-02 |   2.588145e-02 |   5.830046e-04 |   1.962396e-03 |
+----------------+----------------+----------------+----------------+----------------+

nx = 80
+----------------+----------------+----------------+----------------+----------------+
| time           | my_u1          | my_u2          | yaqi_u1        | yaqi_u2        |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   1.214107e-02 |   1.248190e-02 |   1.945367e-04 |   4.048233e-04 |
+----------------+----------------+----------------+----------------+----------------+

nx = 160
+----------------+----------------+----------------+----------------+----------------+
| time           | my_u1          | my_u2          | yaqi_u1        | yaqi_u2        |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   6.042278e-03 |   6.127175e-03 |   1.000789e-04 |   7.100799e-05 |
+----------------+----------------+----------------+----------------+----------------+

# 5/18/18

Currently the tangential multiplier is coming out to be .3 * the normal multiplier regardless
of the applied force on the LHS (I've applied both .4 and .2 of the normal force). The block
does get displaced twice as much when the .4*normal_force is applied.

Vec Object: 1 MPI processes
  type: mpi
Process [0]
-2.43198e-06
-4.05329e-06
-2.43174e-06
-4.05291e-06
2.43188e-06
4.05313e-06
2.43184e-06
4.05307e-06
-2.11758e-22
-8.47033e-22
-2.11758e-22
0.
2.43197e-06
4.05328e-06
2.43175e-06
4.05292e-06
0.
0.
0.
0.
0.
2.11758e-22
-8.47033e-22
0.
0.
0.
0.
-2.11758e-22

k, MatDiffusion for disp_x:
 -0.0000019774636337158934,
 -0.00000097829716956289661,
 0.0000019796167382512884,
 0.00000097614406502750118,
,
CoupledTimeDerivative for disp_x:,
 0.0000019774636337158794,
 0.0000019782971695628961,
 0.0000019807612006925181,
 0.0000019799258705917215,
,
NeumannBC for disp_x:,
 0,
 -0.00000099999999999999995,
 -0.00000099999999999999995,
 0,
,
MatDiffusion for disp_y:,
 -0.0000049610516070032426,
 -0.000004961053100969721,
 0.0000049646401144912286,
 0.0000049574645934817349,
    ,
CoupledTimeDerivative for disp_y:,
 -0.000000038948392996757606,
 -0.000000038946899030276354,
 -0.000000030676983228845547,
 -0.000000030681467618233468,
,
NeumannBC for disp_y:,
 0.0000049999999999999996,
 0.0000049999999999999996,
 0,
 0,

 -0.0000019774636337158934,
 -0.00000097829716956289661,
 0.0000019796167382512884,
 0.00000097614406502750118


 0.0000019774636337158794,
 0.0000019782971695628961,
 0.0000019807612006925181,
 0.0000019799258705917215


 0,
 -0.00000099999999999999995,
 -0.00000099999999999999995,
 0


 -0.0000049610516070032426,
 -0.000004961053100969721,
 0.0000049646401144912286,
 0.0000049574645934817349


 -0.000000038948392996757606,
 -0.000000038946899030276354,
 -0.000000030676983228845547,
 -0.000000030681467618233468


 0.0000049999999999999996,
 0.0000049999999999999996,
 0,
 0

Node coordinates:
0: _coords = ([0] = 0.5, [1] = 1, [2] = 0)
1: _coords = ([0] = -0.5, [1] = 1, [2] = 0)
2: _coords = ([0] = -0.5, [1] = 0, [2] = 0)
3: _coords = ([0] = 0.5, [1] = 0, [2] = 0)

diffusion
    [0] = 0.0000009783099597082335
    [1] = 0.0000019775371028971805
    [2] = -0.00000097723695429118657
    [3] = -0.000001978610108314227

time
    [0] = -0.00000097834501954952698
    [1] = -0.00000097751187447154328
    [2] = -0.00000097997458478785166
    [3] = -0.00000098080861825284306

bc
    [0] = 0
    [1] = -0.00000099999008053285336
    [2] = -0.00000099999008053285336
    [3] = 0

diffusion
    [0] = 0.0000009783099428078857
    [1] = 0.0000019775371033832525
    [2] = -0.00000097723695345141354
    [3] = -0.0000019786100927397242

time
    [0] = -0.00000097834501930010012
    [1] = -0.00000097751187451817753
    [2] = -0.00000097997458471741304
    [3] = -0.00000098080861806376074

bc
    [0] = 0
    [1] = -0.00000099999008053285251
    [2] = -0.00000099999008053285251
    [3] = 0

diffusion
    [0] = -0.0000019775016451033381
    [1] = -0.00000097834538095566629
    [2] = 0.0000019796523756325295
    [3] = 0.0000009761946504264753

time
    [0] = 0.0000019775018079805203
    [1] = 0.0000019783353079571464
    [2] = 0.0000019807993614216109
    [3] = 0.0000019799640636032788

bc
    [0] = 0
    [1] = -0.00000099999007247054262
    [2] = -0.00000099999007247054262
    [3] = 0

contact
0.0000029604617152777595
0.0000029561584508784956

Here's the MatDiffusion disp_x residuals in 4 consecutive iterations:
    [0] = 0.0000009783099561873058
    [1] = 0.000001977537095757813
    [2] = -0.00000097723694759681266
    [3] = -0.0000019786101043483057

    [0] = -0.0000019775016492473375
    [1] = -0.00000097834537469308097
    [2] = 0.0000019796523768443749
    [3] = 0.00000097619464709604382

    [0] = 0.00000097830995726694162
    [1] = 0.000001977537090583882
    [2] = -0.0000009772369458009243
    [3] = -0.0000019786101020498987

    [0] = -0.0000019775016445223984
    [1] = -0.00000097834538006335753
    [2] = 0.0000019796523794872235
    [3] = 0.00000097619464509853261

Perfect ping-pong!

The threshold for poor convergence behavior corresponding to stick is about a
LHS force of 6.1-6.2. The problem converges beautifully with slip. Only a
backtracking line search appears to work for this dumb example which is a bit
troubling. Basic ping pongs back and forth between states as shown above.

row 0: (17, -0.3)  (19, -0.1)
row 1:
row 2: (17, -0.7)  (19, -0.9)
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (19, 1.)
row 13:
row 14: (17, 1.)
row 15:
row 16: (16, 0.00745019)
row 17: (16, 0.00148978)  (17, 0.145842)
row 18: (18, 0.00745019)
row 19: (18, -0.0985102)
row 20:
row 21:
row 22:
row 23:
row 24:
row 25:
row 26:
row 27:
row 28:
row 29:
row 30:
row 31:
row 32:
row 33:
row 34:
row 35:

row 0: (0, -3.66402e-05)  (12, 4.00002e-05)
row 1: (1, 3.66402e-05)  (3, -3.66402e-05)  (33, 0.000200001)
row 2: (0, 3.66402e-05)  (12, -4.00002e-05)
row 3: (33, -0.000200001)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13: (1, -4.00002e-05)  (3, 4.00002e-05)
row 14:
row 15:
row 16: (13, 1.12751e-05)
row 17:
row 18:
row 19:
row 20:
row 21:
row 22:
row 23:
row 24:
row 25:
row 26:
row 27:
row 28:
row 29:
row 30:
row 31:
row 32:
row 33:
row 34:
row 35:

Currently, the normal plus tangential contact problem fails with anything but LU
and SVD preconditioners.

1st nonlinear:

Hand-coded:
row 16: (1, -0.2)   (3, -1.8)  (13, 2.)  (16, 1.)
row 17: (16, 0.1)  (17, -1.)
row 18: (1, -0.6) (3, -1.4)  (15, 2.)  (18, 1.)
row 19:          (19, 1.)


Finite-difference:
row 16: (1, -0.2)  (3, -1.8)  (13, 2.)  (16, 0.99255)
row 17: (16, 0.0985102)  (17, -1.14584)
row 18: (1, -0.6)  (3, -1.4)  (15, 2.)  (18, 0.99255)
row 19: (18, 0.0985102)  (19, 1.)

row 16: (16, 0.00745019)
row 17: (16, 0.00148978)  (17, 0.145842)
row 18: (18, 0.00745019)
row 19: (18, -0.0985102)

Extra terms are appearing in the stick case that don't appear in the FDP:
row 0:
row 1: (17, -1.16376e-05)  (19, -3.56294e-05)
row 2:
row 3: (17, -0.000108325)  (19, -8.43329e-05)
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13: (17, 0.000119962)
row 14:
row 15: (19, 0.000119962)
row 16:
row 17:
row 18:
row 19:
row 20:
row 21:
row 22:
row 23:
row 24:
row 25:
row 26:
row 27:
row 28:
row 29:
row 30:
row 31:
row 32:
row 33:
row 34:
row 35:

row 0:
row 1:
row 2:
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12:
row 13:
row 14:
row 15:
row 16:
row 17: (16, -0.0999998)
row 18:
row 19: (18, -0.0999998)
row 20:
row 21:
row 22:
row 23:
row 24:
row 25:
row 26:
row 27:
row 28:
row 29:
row 30:
row 31:
row 32:
row 33:
row 34:
row 35:

The alternative FB formulation doesn't appear to work any better for stick.

Regularization in the MechanicalContactConstraint between velocity and force
shows some real promise. Need to investigate what happens when the velocity is
small because the block is starting to slide to the left...

Problems show up on the 4th nonlinear iteration

 4 Nonlinear |R| = 3.620750e-07
  ---------- Testing Jacobian -------------
  ||J - Jfd||_F/||J||_F = 3.44142, ||J - Jfd||_F = 24.0096
  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 1.73333)  (1, -4.1544e-06)  (2, 0.766667)  (3, 4.1544e-06)  (4, -0.866667)  (5, 0.)  (6, -1.63333)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -2.66296e-08)  (13, -3.16334e-13)  (14, -7.09351e-09)  (15, -8.42642e-14)  (16, -5.93952e-07)  (17, 0.1)  (18, -1.78186e-06)  (19, -0.3)  (20, 0.0555556)  (21, 0.)  (22, 0.0277778)  (23, 0.)  (24, 0.0138889)  (25, 0.)  (26, 0.0277778)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 1: (0, -4.1544e-06)  (1, 1.73333)  (2, -1.56705e-05)  (3, 0.766667)  (4, 0.)  (5, -0.866667)  (6, 0.)  (7, -1.63333)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 8.96532e-06)  (13, 1.06499e-10)  (14, 1.08596e-05)  (15, 1.29001e-10)  (16, 0.05)  (17, 0.)  (18, 0.15)  (19, 0.)  (20, 0.)  (21, 0.0555556)  (22, 0.)  (23, 0.0277778)  (24, 0.)  (25, 0.0138889)  (26, 0.)  (27, 0.0277778)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 2: (0, 0.766667)  (1, -1.56705e-05)  (2, 1.73333)  (3, 1.56705e-05)  (4, -1.63333)  (5, 0.)  (6, -0.866667)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 2.66296e-08)  (13, 3.16334e-13)  (14, 7.09351e-09)  (15, 8.42642e-14)  (16, -5.34557e-06)  (17, 0.9)  (18, -4.15767e-06)  (19, -0.7)  (20, 0.0277778)  (21, 0.)  (22, 0.0555556)  (23, 0.)  (24, 0.0277778)  (25, 0.)  (26, 0.0138889)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 3: (0, 4.1544e-06)  (1, 0.766667)  (2, 1.56705e-05)  (3, 1.73333)  (4, 0.)  (5, -1.63333)  (6, 0.)  (7, -0.866667)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, -8.96532e-06)  (13, -1.06499e-10)  (14, -1.08596e-05)  (15, -1.29001e-10)  (16, 0.45)  (17, 0.)  (18, 0.35)  (19, 0.)  (20, 0.)  (21, 0.0277778)  (22, 0.)  (23, 0.0555556)  (24, 0.)  (25, 0.0277778)  (26, 0.)  (27, 0.0138889)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 4: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 1.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 5: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 1.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 6: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 1.)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 7: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 1.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 8: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.666667)  (9, 0.)  (10, -0.166667)  (11, 0.)  (12, -0.333333)  (13, 0.)  (14, -0.166667)  (15, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.0111111)  (29, 0.)  (30, 0.00555556)  (31, 0.)  (32, 0.00277778)  (33, 0.)  (34, 0.00555556)  (35, 0.)
row 9: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.666667)  (10, 0.)  (11, -0.166667)  (12, 0.)  (13, -0.333333)  (14, 0.)  (15, -0.166667)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.0111111)  (30, 0.)  (31, 0.00555556)  (32, 0.)  (33, 0.00277778)  (34, 0.)  (35, 0.00555556)
row 10: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.166667)  (9, 0.)  (10, 0.666667)  (11, 0.)  (12, -0.166667)  (13, 0.)  (14, -0.333333)  (15, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.00555556)  (29, 0.)  (30, 0.0111111)  (31, 0.)  (32, 0.00555556)  (33, 0.)  (34, 0.00277778)  (35, 0.)
row 11: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.166667)  (10, 0.)  (11, 0.666667)  (12, 0.)  (13, -0.166667)  (14, 0.)  (15, -0.333333)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.00555556)  (30, 0.)  (31, 0.0111111)  (32, 0.)  (33, 0.00555556)  (34, 0.)  (35, 0.00277778)
row 12: (0, -1.06499e-10)  (1, 8.96532e-06)  (2, 1.06499e-10)  (3, -8.96532e-06)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.333333)  (9, 0.)  (10, -0.166667)  (11, 0.)  (12, 0.666667)  (13, 0.)  (14, -0.166667)  (15, 0.)  (16, 5.93952e-06)  (17, -1.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.00277778)  (29, 0.)  (30, 0.00555556)  (31, 0.)  (32, 0.0111111)  (33, 0.)  (34, 0.00555556)  (35, 0.)
row 13: (0, -1.26511e-15)  (1, 1.06499e-10)  (2, 1.26511e-15)  (3, -1.06499e-10)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.333333)  (10, 0.)  (11, -0.166667)  (12, 0.)  (13, 0.666667)  (14, 0.)  (15, -0.166667)  (16, -0.5)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.00277778)  (30, 0.)  (31, 0.00555556)  (32, 0.)  (33, 0.0111111)  (34, 0.)  (35, 0.00555556)
row 14: (0, -1.29001e-10)  (1, 1.08596e-05)  (2, 1.29001e-10)  (3, -1.08596e-05)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, -0.166667)  (9, 0.)  (10, -0.333333)  (11, 0.)  (12, -0.166667)  (13, 0.)  (14, 0.666667)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 5.93952e-06)  (19, 1.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.00555556)  (29, 0.)  (30, 0.00277778)  (31, 0.)  (32, 0.00555556)  (33, 0.)  (34, 0.0111111)  (35, 0.)
row 15: (0, -1.53241e-15)  (1, 1.29001e-10)  (2, 1.53241e-15)  (3, -1.29001e-10)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, -0.166667)  (10, 0.)  (11, -0.333333)  (12, 0.)  (13, -0.166667)  (14, 0.)  (15, 0.666667)  (16, 0.)  (17, 0.)  (18, -0.5)  (19, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, 0.00555556)  (30, 0.)  (31, 0.00277778)  (32, 0.)  (33, 0.00555556)  (34, 0.)  (35, 0.0111111)
row 16: (0, 1.19087e-06)  (1, -0.100249)  (2, 1.07178e-05)  (3, -0.902244)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (12, -1.19087e-05)  (13, 1.00249)  (16, 3.10791e-06)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)
row 17: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (16, 1.2805e-11)  (17, -1.2805e-10)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (32, -9.09793e-06)  (33, 0.)
row 18: (0, 3.57259e-06)  (1, -0.300747)  (2, 8.33604e-06)  (3, -0.701744)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (14, -1.19086e-05)  (15, 1.00249)  (18, 3.10211e-06)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)
row 19: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (18, 1.57796e-09)  (19, -1.57796e-08)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (34, 1.08247e-05)  (35, 0.)
row 20: (0, 0.0555556)  (1, 0.)  (2, 0.0277778)  (3, 0.)  (4, 0.0138889)  (5, 0.)  (6, 0.0277778)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, -0.555556)  (21, 0.)  (22, -0.277778)  (23, 0.)  (24, -0.138889)  (25, 0.)  (26, -0.277778)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 21: (0, 0.)  (1, 0.0555556)  (2, 0.)  (3, 0.0277778)  (4, 0.)  (5, 0.0138889)  (6, 0.)  (7, 0.0277778)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, -0.555556)  (22, 0.)  (23, -0.277778)  (24, 0.)  (25, -0.138889)  (26, 0.)  (27, -0.277778)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 22: (0, 0.0277778)  (1, 0.)  (2, 0.0555556)  (3, 0.)  (4, 0.0277778)  (5, 0.)  (6, 0.0138889)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, -0.277778)  (21, 0.)  (22, -0.555556)  (23, 0.)  (24, -0.277778)  (25, 0.)  (26, -0.138889)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 23: (0, 0.)  (1, 0.0277778)  (2, 0.)  (3, 0.0555556)  (4, 0.)  (5, 0.0277778)  (6, 0.)  (7, 0.0138889)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, -0.277778)  (22, 0.)  (23, -0.555556)  (24, 0.)  (25, -0.277778)  (26, 0.)  (27, -0.138889)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 24: (0, 0.0138889)  (1, 0.)  (2, 0.0277778)  (3, 0.)  (4, 0.0555556)  (5, 0.)  (6, 0.0277778)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, -0.138889)  (21, 0.)  (22, -0.277778)  (23, 0.)  (24, -0.555556)  (25, 0.)  (26, -0.277778)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 25: (0, 0.)  (1, 0.0138889)  (2, 0.)  (3, 0.0277778)  (4, 0.)  (5, 0.0555556)  (6, 0.)  (7, 0.0277778)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, -0.138889)  (22, 0.)  (23, -0.277778)  (24, 0.)  (25, -0.555556)  (26, 0.)  (27, -0.277778)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 26: (0, 0.0277778)  (1, 0.)  (2, 0.0138889)  (3, 0.)  (4, 0.0277778)  (5, 0.)  (6, 0.0555556)  (7, 0.)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, -0.277778)  (21, 0.)  (22, -0.138889)  (23, 0.)  (24, -0.277778)  (25, 0.)  (26, -0.555556)  (27, 0.)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 27: (0, 0.)  (1, 0.0277778)  (2, 0.)  (3, 0.0138889)  (4, 0.)  (5, 0.0277778)  (6, 0.)  (7, 0.0555556)  (8, 0.)  (9, 0.)  (10, 0.)  (11, 0.)  (12, 0.)  (13, 0.)  (14, 0.)  (15, 0.)  (16, 0.)  (17, 0.)  (18, 0.)  (19, 0.)  (20, 0.)  (21, -0.277778)  (22, 0.)  (23, -0.138889)  (24, 0.)  (25, -0.277778)  (26, 0.)  (27, -0.555556)  (28, 0.)  (29, 0.)  (30, 0.)  (31, 0.)  (32, 0.)  (33, 0.)  (34, 0.)  (35, 0.)
row 28: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.0111111)  (9, 0.)  (10, 0.00555556)  (11, 0.)  (12, 0.00277778)  (13, 0.)  (14, 0.00555556)  (15, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, -0.111111)  (29, 0.)  (30, -0.0555556)  (31, 0.)  (32, -0.0277778)  (33, 0.)  (34, -0.0555556)  (35, 0.)
row 29: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.0111111)  (10, 0.)  (11, 0.00555556)  (12, 0.)  (13, 0.00277778)  (14, 0.)  (15, 0.00555556)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, -0.111111)  (30, 0.)  (31, -0.0555556)  (32, 0.)  (33, -0.0277778)  (34, 0.)  (35, -0.0555556)
row 30: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.00555556)  (9, 0.)  (10, 0.0111111)  (11, 0.)  (12, 0.00555556)  (13, 0.)  (14, 0.00277778)  (15, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, -0.0555556)  (29, 0.)  (30, -0.111111)  (31, 0.)  (32, -0.0555556)  (33, 0.)  (34, -0.0277778)  (35, 0.)
row 31: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.00555556)  (10, 0.)  (11, 0.0111111)  (12, 0.)  (13, 0.00555556)  (14, 0.)  (15, 0.00277778)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, -0.0555556)  (30, 0.)  (31, -0.111111)  (32, 0.)  (33, -0.0555556)  (34, 0.)  (35, -0.0277778)
row 32: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.00277778)  (9, 0.)  (10, 0.00555556)  (11, 0.)  (12, 0.0111111)  (13, 0.)  (14, 0.00555556)  (15, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, -0.0277778)  (29, 0.)  (30, -0.0555556)  (31, 0.)  (32, -0.111111)  (33, 0.)  (34, -0.0555556)  (35, 0.)
row 33: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.00277778)  (10, 0.)  (11, 0.00555556)  (12, 0.)  (13, 0.0111111)  (14, 0.)  (15, 0.00555556)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, -0.0277778)  (30, 0.)  (31, -0.0555556)  (32, 0.)  (33, -0.111111)  (34, 0.)  (35, -0.0555556)
row 34: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.00555556)  (9, 0.)  (10, 0.00277778)  (11, 0.)  (12, 0.00555556)  (13, 0.)  (14, 0.0111111)  (15, 0.)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, -0.0555556)  (29, 0.)  (30, -0.0277778)  (31, 0.)  (32, -0.0555556)  (33, 0.)  (34, -0.111111)  (35, 0.)
row 35: (0, 0.)  (1, 0.)  (2, 0.)  (3, 0.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 0.)  (8, 0.)  (9, 0.00555556)  (10, 0.)  (11, 0.00277778)  (12, 0.)  (13, 0.00555556)  (14, 0.)  (15, 0.0111111)  (20, 0.)  (21, 0.)  (22, 0.)  (23, 0.)  (24, 0.)  (25, 0.)  (26, 0.)  (27, 0.)  (28, 0.)  (29, -0.0555556)  (30, 0.)  (31, -0.0277778)  (32, 0.)  (33, -0.0555556)  (34, 0.)  (35, -0.111111)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 1.73333)  (1, -4.1544e-06)  (2, 0.766667)  (3, 4.1544e-06)  (4, -0.866667)  (6, -1.63333)  (12, -2.66296e-08)  (13, -3.1593e-13)  (14, -7.09351e-09)  (15, -8.43663e-14)  (16, -5.93952e-07)  (17, 0.1)  (18, -1.78186e-06)  (19, -0.3)  (20, 0.0555556)  (22, 0.0277778)  (24, 0.0138889)  (26, 0.0277778)  (32, 1.77971)
row 1: (0, -4.1544e-06)  (1, 1.73333)  (2, -1.56705e-05)  (3, 0.766667)  (4, 4.5469e-13)  (5, -0.866667)  (7, -1.63333)  (10, 2.27345e-13)  (12, 8.96532e-06)  (13, 1.06625e-10)  (14, 1.08596e-05)  (15, 1.29587e-10)  (16, 0.05)  (18, 0.15)  (21, 0.0555556)  (23, 0.0277778)  (24, 2.27345e-13)  (25, 0.0138889)  (27, 0.0277778)  (30, 2.27345e-13)  (33, 2.27345e-13)
row 2: (0, 0.766667)  (1, -1.56705e-05)  (2, 1.73333)  (3, 1.56705e-05)  (4, -1.63333)  (6, -0.866667)  (12, 2.66296e-08)  (13, 3.16152e-13)  (14, 7.09351e-09)  (15, 8.43663e-14)  (16, -5.34557e-06)  (17, 0.9)  (18, -4.15767e-06)  (19, -0.7)  (20, 0.0277778)  (22, 0.0555556)  (24, 0.0277778)  (26, 0.0138889)  (32, 16.0174)
row 3: (0, 4.1544e-06)  (1, 0.766667)  (2, 1.56705e-05)  (3, 1.73333)  (4, -4.5469e-13)  (5, -1.63333)  (7, -0.866667)  (10, -4.5469e-13)  (12, -8.96532e-06)  (13, -1.06397e-10)  (14, -1.08596e-05)  (15, -1.29587e-10)  (16, 0.45)  (18, 0.35)  (21, 0.0277778)  (23, 0.0555556)  (24, -4.5469e-13)  (25, 0.0277778)  (27, 0.0138889)  (30, -4.5469e-13)  (33, -4.5469e-13)
row 4: (4, 1.)
row 5: (5, 1.)
row 6: (6, 1.)
row 7: (7, 1.)
row 8: (8, 0.666667)  (10, -0.166667)  (12, -0.333333)  (14, -0.166667)  (28, 0.0111111)  (30, 0.00555556)  (32, 0.00277778)  (34, 0.00555556)
row 9: (9, 0.666667)  (11, -0.166667)  (13, -0.333333)  (15, -0.166667)  (29, 0.0111111)  (31, 0.00555556)  (33, 0.00277778)  (35, 0.00555556)
row 10: (8, -0.166667)  (10, 0.666667)  (12, -0.166667)  (14, -0.333333)  (28, 0.00555556)  (30, 0.0111111)  (32, 0.00555556)  (34, 0.00277778)
row 11: (9, -0.166667)  (11, 0.666667)  (13, -0.166667)  (15, -0.333333)  (29, 0.00555556)  (31, 0.0111111)  (33, 0.00555556)  (35, 0.00277778)
row 12: (0, -1.06499e-10)  (1, 8.96532e-06)  (2, 1.065e-10)  (3, -8.96532e-06)  (8, -0.333333)  (10, -0.166667)  (12, 0.666667)  (14, -0.166667)  (16, 5.93952e-06)  (17, -1.)  (28, 0.00277778)  (30, 0.00555556)  (32, -17.786)  (34, 0.00555556)
row 13: (0, -9.0938e-13)  (1, 1.05943e-10)  (2, -9.0938e-13)  (3, -1.06852e-10)  (9, -0.333333)  (11, -0.166667)  (13, 0.666667)  (15, -0.166667)  (16, -0.5)  (29, 0.00277778)  (31, 0.00555556)  (33, 0.0111111)  (35, 0.00555556)
row 14: (0, -1.29001e-10)  (1, 1.08596e-05)  (2, 1.29001e-10)  (3, -1.08596e-05)  (8, -0.166667)  (10, -0.333333)  (12, -0.166667)  (14, 0.666667)  (18, 5.93952e-06)  (19, 1.)  (28, 0.00555556)  (30, 0.00277778)  (32, 0.00555556)  (34, 0.0111111)
row 15: (0, -9.0938e-13)  (1, 1.28677e-10)  (2, -9.0938e-13)  (3, -1.29587e-10)  (9, -0.166667)  (11, -0.333333)  (13, -0.166667)  (15, 0.666667)  (18, -0.5)  (29, 0.00555556)  (31, 0.00277778)  (33, 0.00555556)  (35, 0.0111111)
row 16: (0, 1.19087e-06)  (1, -0.10025)  (2, 1.07178e-05)  (3, -0.902311)  (12, -1.19087e-05)  (13, 1.00241)  (16, 3.1074e-06)
row 17: (16, 1.2805e-11)  (17, -1.2805e-10)  (32, 8.94159e-06)
row 18: (0, 3.57259e-06)  (1, -0.300753)  (2, 8.33604e-06)  (3, -0.701777)  (14, -1.19086e-05)  (15, 1.00242)  (18, 3.10168e-06)
row 19: (18, 1.57796e-09)  (19, -1.57796e-08)  (34, 1.08247e-05)
row 20: (0, 0.0555556)  (2, 0.0277778)  (4, 0.0138889)  (6, 0.0277778)  (20, -0.555556)  (22, -0.277778)  (24, -0.138889)  (26, -0.277778)
row 21: (1, 0.0555556)  (3, 0.0277778)  (5, 0.0138889)  (7, 0.0277778)  (21, -0.555556)  (23, -0.277778)  (25, -0.138889)  (27, -0.277778)
row 22: (0, 0.0277778)  (2, 0.0555556)  (4, 0.0277778)  (6, 0.0138889)  (20, -0.277778)  (22, -0.555556)  (24, -0.277778)  (26, -0.138889)
row 23: (1, 0.0277778)  (3, 0.0555556)  (5, 0.0277778)  (7, 0.0138889)  (21, -0.277778)  (23, -0.555556)  (25, -0.277778)  (27, -0.138889)
row 24: (0, 0.0138889)  (2, 0.0277778)  (4, 0.0555556)  (6, 0.0277778)  (20, -0.138889)  (22, -0.277778)  (24, -0.555556)  (26, -0.277778)
row 25: (1, 0.0138889)  (3, 0.0277778)  (5, 0.0555556)  (7, 0.0277778)  (21, -0.138889)  (23, -0.277778)  (25, -0.555556)  (27, -0.277778)
row 26: (0, 0.0277778)  (2, 0.0138889)  (4, 0.0277778)  (6, 0.0555556)  (20, -0.277778)  (22, -0.138889)  (24, -0.277778)  (26, -0.555556)
row 27: (1, 0.0277778)  (3, 0.0138889)  (5, 0.0277778)  (7, 0.0555556)  (21, -0.277778)  (23, -0.138889)  (25, -0.277778)  (27, -0.555556)
row 28: (8, 0.0111111)  (10, 0.00555556)  (12, 0.00277778)  (14, 0.00555556)  (28, -0.111111)  (30, -0.0555556)  (32, -0.0277778)  (34, -0.0555556)
row 29: (9, 0.0111111)  (11, 0.00555556)  (13, 0.00277778)  (15, 0.00555556)  (29, -0.111111)  (31, -0.0555556)  (33, -0.0277778)  (35, -0.0555556)
row 30: (8, 0.00555556)  (10, 0.0111111)  (12, 0.00555556)  (14, 0.00277778)  (28, -0.0555556)  (30, -0.111111)  (32, -0.0555556)  (34, -0.0277778)
row 31: (9, 0.00555556)  (11, 0.0111111)  (13, 0.00555556)  (15, 0.00277778)  (29, -0.0555556)  (31, -0.111111)  (33, -0.0555556)  (35, -0.0277778)
row 32: (8, 0.00277778)  (10, 0.00555556)  (12, 0.0111111)  (14, 0.00555556)  (28, -0.0277778)  (30, -0.0555556)  (32, -0.111111)  (34, -0.0555556)
row 33: (9, 0.00277778)  (11, 0.00555556)  (13, 0.0111111)  (15, 0.00555556)  (29, -0.0277778)  (31, -0.0555556)  (33, -0.111111)  (35, -0.0555556)
row 34: (8, 0.00555556)  (10, 0.00277778)  (12, 0.00555556)  (14, 0.0111111)  (28, -0.0555556)  (30, -0.0277778)  (32, -0.0555556)  (34, -0.111111)
row 35: (9, 0.00555556)  (11, 0.00277778)  (13, 0.00555556)  (15, 0.0111111)  (29, -0.0555556)  (31, -0.0277778)  (33, -0.0555556)  (35, -0.111111)
  Hand-coded minus finite-difference Jacobian with tolerance 1e-05 ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (32, -1.77971)
row 1:
row 2: (32, -16.0174)
row 3:
row 4:
row 5:
row 6:
row 7:
row 8:
row 9:
row 10:
row 11:
row 12: (32, 17.7971)
row 13:
row 14:
row 15:
row 16: (3, 6.73225e-05)  (13, 8.31142e-05)
row 17: (32, -1.80395e-05)
row 18: (3, 3.36221e-05)  (15, 6.86166e-05)
row 19:
row 20:
row 21:
row 22:
row 23:
row 24:
row 25:
row 26:
row 27:
row 28:
row 29:
row 30:
row 31:
row 32:
row 33:
row 34:
row 35:
      SVD: condition number 5.682428267895e+06, 0 of 36 singular values are (nearly) zero
      SVD: smallest singular values: 6.608677756041e-07 8.065014380065e-07 2.778209012446e-02 2.784264392489e-02 8.336675255921e-02
      SVD: largest singular values : 1.394268532912e+00 1.615051296920e+00 1.656164775515e+00 3.716926338142e+00 3.755333729433e+00
      0 Linear |R| = 3.620750e-07
      1 Linear |R| = 1.081871e-21
  Linear solve converged due to CONVERGED_RTOL iterations 1
    |residual|_2 of individual variables:
                 disp_x:     7.36262e-08
                 disp_y:     9.00024e-16
                 lm:         4.36115e-10
                 tangent_lm: 1.49361e-16
                 vel_x:      2.022e-22
                 vel_y:      1.80306e-21

What's interesting is that with a different line search I can change between sticking and slipping

Overall bt seems more robust, but there are times when bt doesn't work and basic does. Ok for a
regularization of 1e-3 and 256 elem blocks, the threshold for slipping is a force of .6 * \mu\lambda_N
(going in increments of 0.1).

The same threshold is observed with 1024 elem blocks. The threshold is .5 * \mu\lambda_N with 64 elem
blocks. If I increase D (pseudo youngs modulus) to 10, then the threshold changes to .3 * \mu\lambda_N.

So for regularization 1e-3, applied force of 2e-5, D=1e-1, lambda=0.1 (for FB function), 64 element
blocks, I get sticking with bt and slipping with basic. Same with F=1e-5. Only at F=0.9e-5 does basic
line search start sticking. And actually changing to lambda=1, basic slips again at 0.9e-5. It starts
sticking at 0.8e-5 = 8e-6. Slipping I believe should only start happening at 1e-4.

We absolutely have to have the absolute value sign on the tangential LM variable in
MechanicalContactConstraint when the tangential velocity is non-zero. E.g. we have to strongly enforce
the fact that the force is in the opposite direction of the slip velocity, otherwise slipping will occur
no matter what.

Correct Jacobian:

row 17: (16, 0.0915653)  (17, -0.91566)  (32, -1.99644)

Hand-coded:

row 17: (16, 0.0915653)  (17, 0.)  (32, -1.99644)

# 6/26/18

LM: num_steps = 30, total nonlinear iterations = 59, total simulation time = 25.3 seconds
Penalty: num_steps = 31, total nonlinear iterations = 89, total simulation time = 50.5 seconds

# 6/28/18

address of material data that (first) StatefulMaterial uses = 0x000000010e28f398

We have:

```
_material_data = 0x000000010e200018
_bnd_material_data = 0x000000010e300018
_neighbor_material_data = 0x000000010e400018
```

- StatefulMaterial adds diffusivity to _bnd_material_data
- The "general" GenericConstantMaterial adds its single property to _material_data (for volume), then
to _bnd_material_data (for face...hmm), then to _neighbor_material_data (for neighbor)
- The "boundary" GenericConstantMaterial add its two properties to
_bnd_material_data

There are two MaterialPropertyStorage objects, `_material_props` and
`_bnd_material_props`. `_material_data` is initialized with `_material_props`
while the other two material data objects are intialized with
`_bnd_material_props`.

(MaterialPropertyStorage *) $3 = 0x000000010ce45890 (this is the address for the
boundary)
"diffusivity", 0

0x000000010ce456d0 this is address for non-boundary MaterialPropertyStorage

So `_prop_ids` in `MaterialPropertyStorage` is static so that explains the
nullptr in the volume material

The problem is with `_stateful_prop_id_to_prop_id` in `MaterialPropertyStorage`

`_prop_ids` is static, but `_stateful_prop_id_to_prop_id` is not.

Ok lets now add a volume stateful material.

_props_ids:
"diffusivity", 0
"diffusivity2", 1
"dummy", 2
"D", 3
"D_neighbor", 4

_bnd_material_props (MaterialPropertyStorage): _stateful_prop_id_to_prop_id
[0] : 0
[1] : 1

_material_props (MaterialPropertyStorage): _stateful_prop_id_to_prop_id
[0] : 1

_material_data (uses _material_props)
[0] : "diffusivity" : null
[1] : "diffusivity2" : non-null
[2] : "dummy" : non-null

_bnd_material_data (uses _bnd_material_props)
[0] : "diffusivity" : non-null
[1] : "diffusivity2" : non-null
[2] : "dummy" : non-null
[3] : "D" : non-null
[4] : "D_neighbor" : non-null

_neighbor_material_data (uses _bnd_material_props)
[0] : "diffusivity" : null
[1] : "diffusivity2" : non-null
[2] : "dummy" : non-null


Without the volume stateful material:

static MaterialPropertyStorage member _props_ids:
"diffusivity", 0
"dummy", 1
"D", 2
"D_neighbor", 3

_bnd_material_props (MaterialPropertyStorage): _stateful_prop_id_to_prop_id
[0] : 0

_material_props (MaterialPropertyStorage): _stateful_prop_id_to_prop_id
Empty

_material_data (uses _material_props)
[0] : "diffusivity" : null
[1] : "dummy" : non-null

_bnd_material_data (uses _bnd_material_props)
[0] : "diffusivity" : non-null
[1] : "dummy" : non-null
[2] : "D" : non-null
[3] : "D_neighbor" : non-null

_neighbor_material_data (uses _bnd_material_props)
[0] : "diffusivity" : null
[1] : "dummy" : non-null

 0 Nonlinear |R| = 2.185303e+00
  ---------- Testing Jacobian -------------
  Testing hand-coded Jacobian, if (for double precision runs) ||J - Jfd||_F/||J||_F is
    O(1.e-8), the hand-coded Jacobian is probably correct.
  ||J - Jfd||_F/||J||_F = 9.71572e-09, ||J - Jfd||_F = 1.55451e-07
  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 0.)  (1, 0.)  (2, 0.)  (3, -4.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, 4.)
row 1: (0, 0.)  (1, 0.)  (2, -4.)  (3, -8.)  (4, 0.)  (5, 0.)  (6, 4.)  (7, 0.)
row 2:
row 3:
row 4: (0, 0.)  (1, 0.)  (2, 0.)  (3, 4.)  (4, 0.)  (5, 0.)  (6, 0.)  (7, -4.)
row 5: (0, 0.)  (1, 0.)  (2, 4.)  (3, 8.)  (4, 0.)  (5, 0.)  (6, -4.)  (7, 0.)
row 6:
row 7:
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes

New Petsc:
  * frame #0: 0x0000000101d20e70 libmoose-dbg.0.dylib`NonlinearSystemBase::computeResidualTags(this=0x000000010f005618, tags=size=2) at NonlinearSystemBase.C:542
    frame #1: 0x00000001019fad87 libmoose-dbg.0.dylib`FEProblemBase::computeResidualTags(this=0x000000010f002018, tags=size=2) at FEProblemBase.C:4227
    frame #2: 0x00000001019f9cb1 libmoose-dbg.0.dylib`FEProblemBase::computeResidualInternal(this=0x000000010f002018, soln=0x000000010e56f2e0, residual=0x00007fff5fbfec00, tags=size=2) at FEProblemBase.C:4109
    frame #3: 0x00000001019f96ff libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010f002018, soln=0x000000010e56f2e0, residual=0x00007fff5fbfec00) at FEProblemBase.C:4066
    frame #4: 0x00000001019f90a6 libmoose-dbg.0.dylib`FEProblemBase::computeResidualSys(this=0x000000010f002018, (null)=0x000000010e56ecf0, soln=0x000000010e56f2e0, residual=0x00007fff5fbfec00) at FEProblemBase.C:4043
    frame #5: 0x0000000101cfde59 libmoose-dbg.0.dylib`ComputeResidualFunctor::residual(this=0x000000010f006688, soln=0x000000010e56f2e0, residual=0x00007fff5fbfec00, sys=0x000000010e56ecf0) at ComputeResidualFunctor.C:23
    frame #6: 0x00000001058356f2 libmesh_dbg.0.dylib`::libmesh_petsc_snes_residual(snes=0x000000010f022c20, x=0x000000010f00be20, r=0x000000010f00e620, ctx=0x000000010e56fa20) at petsc_nonlinear_solver.C:173
    frame #7: 0x000000010b47189f libpetsc.3.09.dylib`SNESComputeFunction(snes=0x000000010f022c20, x=0x000000010f00be20, y=0x000000010f00e620) at snes.c:2209 [opt]
    frame #8: 0x000000010b4a606d libpetsc.3.09.dylib`SNESSolve_KSPONLY(snes=0x000000010f022c20) at ksponly.c:23 [opt]
    frame #9: 0x000000010b478650 libpetsc.3.09.dylib`SNESSolve(snes=0x000000010f022c20, b=0x0000000000000000, x=0x000000010f00be20) at snes.c:4348 [opt]
    frame #10: 0x0000000105833acd libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x000000010e56fa20, pre_in=0x000000010e56f950, x_in=0x000000010e56f1c0, r_in=0x000000010e56f770, (null)=0.00000001, (null)=10000) at petsc_nonlinear_solver.C:934
    frame #11: 0x000000010595e809 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x000000010e56ecf0) at nonlinear_implicit_system.C:181
    frame #12: 0x0000000101d07ab8 libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000010f005618) at NonlinearSystem.C:191
    frame #13: 0x00000001019f4ee2 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x000000010f002018) at FEProblemBase.C:3827
    frame #14: 0x0000000100f7ea0e libmoose-dbg.0.dylib`Steady::execute(this=0x000000010e580938) at Steady.C:84
    frame #15: 0x000000010219b823 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e80a018) at MooseApp.C:801
    frame #16: 0x000000010219c01b libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e80a018) at MooseApp.C:904
    frame #17: 0x0000000100002373 moose_test-dbg`main(argc=3, argv=0x00007fff5fbff850) at main.C:31
    frame #18: 0x00007fffb5456235 libdyld.dylib`start + 1

Old Petsc:
  * frame #0: 0x0000000101d25600 libmoose-dbg.0.dylib`NonlinearSystemBase::computeResidualTags(this=0x000000010f005618, tags=size=2) at NonlinearSystemBase.C:542
    frame #1: 0x00000001019f9ba7 libmoose-dbg.0.dylib`FEProblemBase::computeResidualTags(this=0x000000010f002018, tags=size=2) at FEProblemBase.C:4227
    frame #2: 0x00000001019f8ad1 libmoose-dbg.0.dylib`FEProblemBase::computeResidualInternal(this=0x000000010f002018, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20, tags=size=2) at FEProblemBase.C:4109
    frame #3: 0x00000001019f851f libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010f002018, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20) at FEProblemBase.C:4066
    frame #4: 0x00000001019f7ec6 libmoose-dbg.0.dylib`FEProblemBase::computeResidualSys(this=0x000000010f002018, (null)=0x000000010e26ece0, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20) at FEProblemBase.C:4043
    frame #5: 0x0000000101cfdca9 libmoose-dbg.0.dylib`ComputeResidualFunctor::residual(this=0x000000010f006688, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20, sys=0x000000010e26ece0) at ComputeResidualFunctor.C:23
    frame #6: 0x00000001056e4fe2 libmesh_dbg.0.dylib`::libmesh_petsc_snes_residual(snes=0x000000010f800e20, x=0x000000010f00b820, r=0x000000010f00d620, ctx=0x000000010e26fa10) at petsc_nonlinear_solver.C:173
    frame #7: 0x000000010b72302f libpetsc.3.7.dylib`SNESComputeFunction(snes=0x000000010f800e20, x=0x000000010f00b820, y=0x000000010f00d620) at snes.c:2145 [opt]
    frame #8: 0x000000010b754f3d libpetsc.3.7.dylib`SNESSolve_KSPONLY(snes=0x000000010f800e20) at ksponly.c:25 [opt]
    frame #9: 0x000000010b72828a libpetsc.3.7.dylib`SNESSolve(snes=0x000000010f800e20, b=<unavailable>, x=0x000000010f00b820) at snes.c:4005 [opt]
    frame #10: 0x00000001056e33c7 libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x000000010e26fa10, pre_in=0x000000010e26f940, x_in=0x000000010e26f1b0, r_in=0x000000010e26f760, (null)=0.00000001, (null)=10000) at petsc_nonlinear_solver.C:930
    frame #11: 0x000000010581ef59 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x000000010e26ece0) at nonlinear_implicit_system.C:181
    frame #12: 0x0000000101d0d0b8 libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000010f005618) at NonlinearSystem.C:191
    frame #13: 0x00000001019f3d02 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x000000010f002018) at FEProblemBase.C:3827
    frame #14: 0x0000000100f7d5ce libmoose-dbg.0.dylib`Steady::execute(this=0x000000010e280928) at Steady.C:84
    frame #15: 0x00000001021a6253 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e809a18) at MooseApp.C:801
    frame #16: 0x00000001021a6a4b libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e809a18) at MooseApp.C:904
    frame #17: 0x0000000100001b63 moose_test-dbg`main(argc=3, argv=0x00007fff5fbff850) at main.C:31
    frame #18: 0x00007fffb5456235 libdyld.dylib`start + 1

  * frame #0: 0x0000000101d25600 libmoose-dbg.0.dylib`NonlinearSystemBase::computeResidualTags(this=0x000000010f005618, tags=size=2) at NonlinearSystemBase.C:542
    frame #1: 0x00000001019f9ba7 libmoose-dbg.0.dylib`FEProblemBase::computeResidualTags(this=0x000000010f002018, tags=size=2) at FEProblemBase.C:4227
    frame #2: 0x00000001019f8ad1 libmoose-dbg.0.dylib`FEProblemBase::computeResidualInternal(this=0x000000010f002018, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20, tags=size=2) at FEProblemBase.C:4109
    frame #3: 0x00000001019f851f libmoose-dbg.0.dylib`FEProblemBase::computeResidual(this=0x000000010f002018, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20) at FEProblemBase.C:4066
    frame #4: 0x00000001019f7ec6 libmoose-dbg.0.dylib`FEProblemBase::computeResidualSys(this=0x000000010f002018, (null)=0x000000010e26ece0, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20) at FEProblemBase.C:4043
    frame #5: 0x0000000101cfdca9 libmoose-dbg.0.dylib`ComputeResidualFunctor::residual(this=0x000000010f006688, soln=0x000000010e26f2d0, residual=0x00007fff5fbfec20, sys=0x000000010e26ece0) at ComputeResidualFunctor.C:23
    frame #6: 0x00000001056e4fe2 libmesh_dbg.0.dylib`::libmesh_petsc_snes_residual(snes=0x000000010f800e20, x=0x000000010f00b820, r=0x000000010f00d620, ctx=0x000000010e26fa10) at petsc_nonlinear_solver.C:173
    frame #7: 0x000000010b72302f libpetsc.3.7.dylib`SNESComputeFunction(snes=0x000000010f800e20, x=0x000000010f00b820, y=0x000000010f00d620) at snes.c:2145 [opt]
    frame #8: 0x000000010b7550c7 libpetsc.3.7.dylib`SNESSolve_KSPONLY(snes=0x000000010f800e20) at ksponly.c:51 [opt]
    frame #9: 0x000000010b72828a libpetsc.3.7.dylib`SNESSolve(snes=0x000000010f800e20, b=<unavailable>, x=0x000000010f00b820) at snes.c:4005 [opt]
    frame #10: 0x00000001056e33c7 libmesh_dbg.0.dylib`libMesh::PetscNonlinearSolver<double>::solve(this=0x000000010e26fa10, pre_in=0x000000010e26f940, x_in=0x000000010e26f1b0, r_in=0x000000010e26f760, (null)=0.00000001, (null)=10000) at petsc_nonlinear_solver.C:930
    frame #11: 0x000000010581ef59 libmesh_dbg.0.dylib`libMesh::NonlinearImplicitSystem::solve(this=0x000000010e26ece0) at nonlinear_implicit_system.C:181
    frame #12: 0x0000000101d0d0b8 libmoose-dbg.0.dylib`NonlinearSystem::solve(this=0x000000010f005618) at NonlinearSystem.C:191
    frame #13: 0x00000001019f3d02 libmoose-dbg.0.dylib`FEProblemBase::solve(this=0x000000010f002018) at FEProblemBase.C:3827
    frame #14: 0x0000000100f7d5ce libmoose-dbg.0.dylib`Steady::execute(this=0x000000010e280928) at Steady.C:84
    frame #15: 0x00000001021a6253 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010e809a18) at MooseApp.C:801
    frame #16: 0x00000001021a6a4b libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010e809a18) at MooseApp.C:904
    frame #17: 0x0000000100001b63 moose_test-dbg`main(argc=3, argv=0x00007fff5fbff850) at main.C:31
    frame #18: 0x00007fffb5456235 libdyld.dylib`start + 1

# 7/2/18

With FD:
+----------------+----------------+----------------+
| time           | cumulative     | num_nl         |
+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   6.000000e+00 |   6.000000e+00 |
|   1.500000e+01 |   1.200000e+01 |   6.000000e+00 |
|   2.500000e+01 |   1.900000e+01 |   7.000000e+00 |
|   3.500000e+01 |   2.100000e+01 |   2.000000e+00 |
|   4.500000e+01 |   2.300000e+01 |   2.000000e+00 |
|   5.500000e+01 |   2.500000e+01 |   2.000000e+00 |
|   6.500000e+01 |   2.700000e+01 |   2.000000e+00 |
|   7.500000e+01 |   2.900000e+01 |   2.000000e+00 |
|   8.500000e+01 |   3.100000e+01 |   2.000000e+00 |
|   9.500000e+01 |   3.300000e+01 |   2.000000e+00 |
|   1.000000e+02 |   3.500000e+01 |   2.000000e+00 |
+----------------+----------------+----------------+

With NEWTON:
+----------------+----------------+----------------+
| time           | cumulative     | num_nl         |
+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   6.000000e+00 |   6.000000e+00 |
|   2.000000e+01 |   1.400000e+01 |   8.000000e+00 |
|   3.000000e+01 |   1.600000e+01 |   2.000000e+00 |
|   4.000000e+01 |   1.800000e+01 |   2.000000e+00 |
|   5.000000e+01 |   2.000000e+01 |   2.000000e+00 |
|   6.000000e+01 |   2.200000e+01 |   2.000000e+00 |
|   7.000000e+01 |   2.400000e+01 |   2.000000e+00 |
|   8.000000e+01 |   2.600000e+01 |   2.000000e+00 |
|   9.000000e+01 |   2.800000e+01 |   2.000000e+00 |
|   1.000000e+02 |   3.000000e+01 |   2.000000e+00 |
+----------------+----------------+----------------+

With PJFNK:
+----------------+----------------+----------------+
| time           | cumulative     | num_nl         |
+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   6.000000e+00 |   6.000000e+00 |
|   2.000000e+01 |   1.400000e+01 |   8.000000e+00 |
|   3.000000e+01 |   1.600000e+01 |   2.000000e+00 |
|   4.000000e+01 |   1.800000e+01 |   2.000000e+00 |
|   5.000000e+01 |   2.000000e+01 |   2.000000e+00 |
|   6.000000e+01 |   2.200000e+01 |   2.000000e+00 |
|   7.000000e+01 |   2.400000e+01 |   2.000000e+00 |
|   8.000000e+01 |   2.600000e+01 |   2.000000e+00 |
|   9.000000e+01 |   2.800000e+01 |   2.000000e+00 |
|   1.000000e+02 |   3.000000e+01 |   2.000000e+00 |
+----------------+----------------+----------------+

# 7/3/18

With small two-block problem:
Omitting normals and perturbations from pmat:
PJFNK:
+----------------+----------------+----------------+----------------+----------------+
| time           | linear_its     | nl_its         | total_linear   | total_nl       |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   9.000000e+00 |   6.000000e+00 |   9.000000e+00 |   6.000000e+00 |
+----------------+----------------+----------------+----------------+----------------+
NEWTON:
+----------------+----------------+----------------+----------------+----------------+
| time           | linear_its     | nl_its         | total_linear   | total_nl       |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   8.000000e+00 |   8.000000e+00 |   8.000000e+00 |   8.000000e+00 |
+----------------+----------------+----------------+----------------+----------------+

With normals and perturbations in pmat:
PJFNK
+----------------+----------------+----------------+----------------+----------------+
| time           | linear_its     | nl_its         | total_linear   | total_nl       |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   1.100000e+01 |   6.000000e+00 |   1.100000e+01 |   6.000000e+00 |
+----------------+----------------+----------------+----------------+----------------+
NEWTON
+----------------+----------------+----------------+----------------+----------------+
| time           | linear_its     | nl_its         | total_linear   | total_nl       |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+00 |   8.000000e+00 |   8.000000e+00 |   8.000000e+00 |   8.000000e+00 |
+----------------+----------------+----------------+----------------+----------------+

With sliding block problem:
Omitting normals and perturbations from pmat:
PJFNK:
timesteps = 33
+----------------+----------------+----------------+----------------+-------------------+
| time           | linear_its     | nonlinear_its  | tot_linear_its | tot_nonlinear_its |
+----------------+----------------+----------------+----------------+-------------------+
:                :                :                :                :                   :
|   1.775000e+00 |   1.000000e+00 |   1.000000e+00 |   8.500000e+01 |      3.600000e+01 |
|   1.800000e+00 |   3.000000e+01 |   7.000000e+00 |   1.150000e+02 |      4.300000e+01 |
|   1.850000e+00 |   9.000000e+00 |   3.000000e+00 |   1.240000e+02 |      4.600000e+01 |
|   1.950000e+00 |   1.100000e+01 |   3.000000e+00 |   1.350000e+02 |      4.900000e+01 |
|   2.050000e+00 |   1.100000e+01 |   3.000000e+00 |   1.460000e+02 |      5.200000e+01 |
|   2.150000e+00 |   9.000000e+00 |   4.000000e+00 |   1.550000e+02 |      5.600000e+01 |
|   2.250000e+00 |   1.000000e+00 |   1.000000e+00 |   1.560000e+02 |      5.700000e+01 |
|   2.350000e+00 |   1.000000e+00 |   1.000000e+00 |   1.570000e+02 |      5.800000e+01 |
|   2.450000e+00 |   1.000000e+00 |   1.000000e+00 |   1.580000e+02 |      5.900000e+01 |
|   2.550000e+00 |   1.000000e+00 |   1.000000e+00 |   1.590000e+02 |      6.000000e+01 |
|   2.650000e+00 |   1.000000e+00 |   1.000000e+00 |   1.600000e+02 |      6.100000e+01 |
|   2.750000e+00 |   1.000000e+00 |   1.000000e+00 |   1.610000e+02 |      6.200000e+01 |
|   2.850000e+00 |   1.000000e+00 |   1.000000e+00 |   1.620000e+02 |      6.300000e+01 |
|   2.950000e+00 |   1.000000e+00 |   1.000000e+00 |   1.630000e+02 |      6.400000e+01 |
|   3.000000e+00 |   1.000000e+00 |   1.000000e+00 |   1.640000e+02 |      6.500000e+01 |
+----------------+----------------+----------------+----------------+-------------------+

NEWTON:
timesteps = 30
+----------------+----------------+----------------+----------------+-------------------+
| time           | linear_its     | nonlinear_its  | tot_linear_its | tot_nonlinear_its |
+----------------+----------------+----------------+----------------+-------------------+
:                :                :                :                :                   :
|   1.600000e+00 |   1.000000e+00 |   1.000000e+00 |   3.900000e+01 |      3.900000e+01 |
|   1.700000e+00 |   1.000000e+00 |   1.000000e+00 |   4.000000e+01 |      4.000000e+01 |
|   1.800000e+00 |   8.000000e+00 |   8.000000e+00 |   4.800000e+01 |      4.800000e+01 |
|   1.900000e+00 |   6.000000e+00 |   6.000000e+00 |   5.400000e+01 |      5.400000e+01 |
|   2.000000e+00 |   7.000000e+00 |   7.000000e+00 |   6.100000e+01 |      6.100000e+01 |
|   2.100000e+00 |   5.000000e+00 |   5.000000e+00 |   6.600000e+01 |      6.600000e+01 |
|   2.200000e+00 |   3.000000e+00 |   3.000000e+00 |   6.900000e+01 |      6.900000e+01 |
|   2.300000e+00 |   1.000000e+00 |   1.000000e+00 |   7.000000e+01 |      7.000000e+01 |
|   2.400000e+00 |   1.000000e+00 |   1.000000e+00 |   7.100000e+01 |      7.100000e+01 |
|   2.500000e+00 |   1.000000e+00 |   1.000000e+00 |   7.200000e+01 |      7.200000e+01 |
|   2.600000e+00 |   1.000000e+00 |   1.000000e+00 |   7.300000e+01 |      7.300000e+01 |
|   2.700000e+00 |   1.000000e+00 |   1.000000e+00 |   7.400000e+01 |      7.400000e+01 |
|   2.800000e+00 |   1.000000e+00 |   1.000000e+00 |   7.500000e+01 |      7.500000e+01 |
|   2.900000e+00 |   1.000000e+00 |   1.000000e+00 |   7.600000e+01 |      7.600000e+01 |
|   3.000000e+00 |   1.000000e+00 |   1.000000e+00 |   7.700000e+01 |      7.700000e+01 |
+----------------+----------------+----------------+----------------+-------------------+


With normals and perturbations in pmat:
PJFNK
timesteps = 32
+----------------+----------------+----------------+----------------+-------------------+
| time           | linear_its     | nonlinear_its  | tot_linear_its | tot_nonlinear_its |
+----------------+----------------+----------------+----------------+-------------------+
:                :                :                :                :                   :
|   1.675000e+00 |   1.000000e+00 |   1.000000e+00 |   9.300000e+01 |      3.700000e+01 |
|   1.775000e+00 |   1.000000e+00 |   1.000000e+00 |   9.400000e+01 |      3.800000e+01 |
|   1.875000e+00 |   2.800000e+01 |   8.000000e+00 |   1.220000e+02 |      4.600000e+01 |
|   1.975000e+00 |   1.700000e+01 |   4.000000e+00 |   1.390000e+02 |      5.000000e+01 |
|   2.075000e+00 |   1.500000e+01 |   4.000000e+00 |   1.540000e+02 |      5.400000e+01 |
|   2.175000e+00 |   7.000000e+00 |   3.000000e+00 |   1.610000e+02 |      5.700000e+01 |
|   2.275000e+00 |   1.000000e+00 |   1.000000e+00 |   1.620000e+02 |      5.800000e+01 |
|   2.375000e+00 |   1.000000e+00 |   1.000000e+00 |   1.630000e+02 |      5.900000e+01 |
|   2.475000e+00 |   1.000000e+00 |   1.000000e+00 |   1.640000e+02 |      6.000000e+01 |
|   2.575000e+00 |   1.000000e+00 |   1.000000e+00 |   1.650000e+02 |      6.100000e+01 |
|   2.675000e+00 |   1.000000e+00 |   1.000000e+00 |   1.660000e+02 |      6.200000e+01 |
|   2.775000e+00 |   1.000000e+00 |   1.000000e+00 |   1.670000e+02 |      6.300000e+01 |
|   2.875000e+00 |   1.000000e+00 |   1.000000e+00 |   1.680000e+02 |      6.400000e+01 |
|   2.975000e+00 |   1.000000e+00 |   1.000000e+00 |   1.690000e+02 |      6.500000e+01 |
|   3.000000e+00 |   1.000000e+00 |   1.000000e+00 |   1.700000e+02 |      6.600000e+01 |
+----------------+----------------+----------------+----------------+-------------------+

NEWTON
timesteps = 30
+----------------+----------------+----------------+----------------+-------------------+
| time           | linear_its     | nonlinear_its  | tot_linear_its | tot_nonlinear_its |
+----------------+----------------+----------------+----------------+-------------------+
:                :                :                :                :                   :
|   1.600000e+00 |   1.000000e+00 |   1.000000e+00 |   4.000000e+01 |      4.000000e+01 |
|   1.700000e+00 |   1.000000e+00 |   1.000000e+00 |   4.100000e+01 |      4.100000e+01 |
|   1.800000e+00 |   8.000000e+00 |   8.000000e+00 |   4.900000e+01 |      4.900000e+01 |
|   1.900000e+00 |   7.000000e+00 |   7.000000e+00 |   5.600000e+01 |      5.600000e+01 |
|   2.000000e+00 |   7.000000e+00 |   7.000000e+00 |   6.300000e+01 |      6.300000e+01 |
|   2.100000e+00 |   4.000000e+00 |   4.000000e+00 |   6.700000e+01 |      6.700000e+01 |
|   2.200000e+00 |   3.000000e+00 |   3.000000e+00 |   7.000000e+01 |      7.000000e+01 |
|   2.300000e+00 |   1.000000e+00 |   1.000000e+00 |   7.100000e+01 |      7.100000e+01 |
|   2.400000e+00 |   1.000000e+00 |   1.000000e+00 |   7.200000e+01 |      7.200000e+01 |
|   2.500000e+00 |   1.000000e+00 |   1.000000e+00 |   7.300000e+01 |      7.300000e+01 |
|   2.600000e+00 |   1.000000e+00 |   1.000000e+00 |   7.400000e+01 |      7.400000e+01 |
|   2.700000e+00 |   1.000000e+00 |   1.000000e+00 |   7.500000e+01 |      7.500000e+01 |
|   2.800000e+00 |   1.000000e+00 |   1.000000e+00 |   7.600000e+01 |      7.600000e+01 |
|   2.900000e+00 |   1.000000e+00 |   1.000000e+00 |   7.700000e+01 |      7.700000e+01 |
|   3.000000e+00 |   1.000000e+00 |   1.000000e+00 |   7.800000e+01 |      7.800000e+01 |
+----------------+----------------+----------------+----------------+-------------------+

# 7/9/18

Test results

lindad@inl606513:~/projects2/moose/test((HEAD detached at 1b4ec50b65))$ ./run_tests -j24 --re lower_d -p 2
Final Test Results:
--------------------------------------------------------------------------------------------------------------
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK
mesh_modifiers/lower_d_block.second_order ................................................... FAILED (EXODIFF)
--------------------------------------------------------------------------------------------------------------
Ran 2 tests in 2.6 seconds.
1 passed, 0 skipped, 0 pending, 1 FAILED
lindad@inl606513:~/projects2/moose/test((HEAD detached at 1b4ec50b65))$ ./run_tests -j24 --re lower_d -p 4
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK
mesh_modifiers/lower_d_block.second_order ................................................................. OK
--------------------------------------------------------------------------------------------------------------
Ran 2 tests in 2.6 seconds.
2 passed, 0 skipped, 0 pending, 0 failed
lindad@inl606513:~/projects2/moose/test((HEAD detached at 1b4ec50b65))$ ./run_tests -j24 --re lower_d -p 3
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK
mesh_modifiers/lower_d_block.second_order ................................................................. OK
--------------------------------------------------------------------------------------------------------------
Ran 2 tests in 2.5 seconds.
2 passed, 0 skipped, 0 pending, 0 failed
lindad@inl606513:~/projects2/moose/test((HEAD detached at 1b4ec50b65))$ ./run_tests -j24 --re lower_d -p 2
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK
mesh_modifiers/lower_d_block.second_order ................................................................. OK
--------------------------------------------------------------------------------------------------------------
Ran 2 tests in 2.5 seconds.
2 passed, 0 skipped, 0 pending, 0 failed
lindad@inl606513:~/projects2/moose/test((HEAD detached at 1b4ec50b65))$ ./run_tests -j24 --re lower_d -p 2
mesh_modifiers/lower_d_block.second_order ................................................................. OK
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK
--------------------------------------------------------------------------------------------------------------
Ran 2 tests in 2.6 seconds.
2 passed, 0 skipped, 0 pending, 0 failed
lindad@inl606513:~/projects2/moose/test((HEAD detached at 1b4ec50b65))$ ./run_tests -j24 --re lower_d -p 2
mesh_modifiers/lower_d_block.second_order: Working Directory: /Users/lindad/projects2/moose/test/tests/mesh_modifiers/lower_d_block
mesh_modifiers/lower_d_block.second_order: Running command: mpiexec -n 2 /Users/lindad/projects2/moose/test/moose_test-opt --n-threads=1 -i lower_d.i Mesh/second_order=true Outputs/file_base=lower_d_second_order_out GlobalParams/order=SECOND Mesh/nx=5 Mesh/ny=5 --error --error-unused --error-override --no-gdb-backtrace
mesh_modifiers/lower_d_block.second_order:
mesh_modifiers/lower_d_block.second_order:
mesh_modifiers/lower_d_block.second_order: ===================================================================================
mesh_modifiers/lower_d_block.second_order: =   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
mesh_modifiers/lower_d_block.second_order: =   PID 51038 RUNNING AT inl606513.inl.gov
mesh_modifiers/lower_d_block.second_order: =   EXIT CODE: 11
mesh_modifiers/lower_d_block.second_order: =   CLEANING UP REMAINING PROCESSES
mesh_modifiers/lower_d_block.second_order: =   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
mesh_modifiers/lower_d_block.second_order: ===================================================================================
mesh_modifiers/lower_d_block.second_order: YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Segmentation fault: 11 (signal 11)
mesh_modifiers/lower_d_block.second_order: This typically refers to a problem with your application.
mesh_modifiers/lower_d_block.second_order: Please see the FAQ page for debugging suggestions
mesh_modifiers/lower_d_block.second_order:
mesh_modifiers/lower_d_block.second_order: ################################################################################
mesh_modifiers/lower_d_block.second_order: Tester failed, reason: CRASH
mesh_modifiers/lower_d_block.second_order:
mesh_modifiers/lower_d_block.second_order ..................................................... FAILED (CRASH)
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK


Final Test Results:
--------------------------------------------------------------------------------------------------------------
mesh_modifiers/lower_d_block.lower_d ...................................................................... OK
mesh_modifiers/lower_d_block.second_order ..................................................... FAILED (CRASH)
--------------------------------------------------------------------------------------------------------------
Ran 2 tests in 2.7 seconds.
1 passed, 0 skipped, 0 pending, 1 FAILED

So the problem is that p_my_neighbor is non-null and different between the two
processors. The elems are different

Still same error with distributed mesh:
  * frame #0: 0x000000010afb3f7a libmesh_dbg.0.dylib`libMesh::MacroFunctions::report_error(file="../src/mesh/distributed_mesh.C", line=344, date="nodate", time="notime") at libmesh_common.C:74
    frame #1: 0x000000010bb8d72b libmesh_dbg.0.dylib`libMesh::DistributedMesh::elem_ptr(this=0x00007f92bd5e7cc0, i=36) at distributed_mesh.C:344
    frame #2: 0x000000010ad71923 libmesh_dbg.0.dylib`libMesh::DofMap::elem_ptr(this=0x00007f92bec44060, mesh=0x00007f92bd5e7cc0, i=36) const at dof_map.C:307
...
    frame #8: 0x000000010ad96598 libmesh_dbg.0.dylib`void libMesh::DofMap::set_nonlocal_dof_objects<libMesh::MeshBase::element_iterator>(tis=0x00007f92bec44060, objects_begin=element_iterator @ 0x00007fff594eb800, objects_end=element_iterator @ 0x00007fff594eb7e0, mesh=0x0000f92bd5e7cc0, objects=f0 18 d7 0a 01 00 00 00 00 00 00 00 00 00 00 00)(libMesh::MeshBase&, unsigned int) const) at dof_map.C:438
    frame #9: 0x000000010ad84f90 libmesh_dbg.0.dylib`libMesh::DofMap::distribute_dofs(this=0x00007f92bec44060, mesh=0x00007f92bd5e7cc0) atdof_map.C:971
    frame #10: 0x000000010c97439f libmesh_dbg.0.dylib`libMesh::System::init_data(this=0x00007f92bec439b0) at system.C:270
    frame #11: 0x000000010c93d5ab libmesh_dbg.0.dylib`libMesh::ImplicitSystem::init_data(this=0x00007f92bec439b0) at implicit_system.C:93
    frame #12: 0x000000010c974234 libmesh_dbg.0.dylib`libMesh::System::init(this=0x00007f92bec439b0) at system.C:248
    frame #13: 0x000000010c8b0bbd libmesh_dbg.0.dylib`libMesh::EquationSystems::init(this=0x00007f92c0000fa0) at equation_systems.C:109
    frame #14: 0x0000000108094223 libmoose-dbg.0.dylib`FEProblemBase::init(this=0x00007f92c0000c18) at FEProblemBase.C:3786
    frame #15: 0x0000000107259adc libmoose-dbg.0.dylib`InitProblemAction::act(this=0x00007f92bd5b3308) at InitProblemAction.C:29
    frame #16: 0x0000000107218549 libmoose-dbg.0.dylib`ActionWarehouse::executeActionsWithAction(this=0x00007f92bd84cfb8, task="init_problm") at ActionWarehouse.C:367
    frame #17: 0x0000000107216e57 libmoose-dbg.0.dylib`ActionWarehouse::executeAllActions(this=0x00007f92bd84cfb8) at ActionWarehouse.C:33
    frame #18: 0x00000001088427f7 libmoose-dbg.0.dylib`MooseApp::runInputFile(this=0x00007f92bd84ca18) at MooseApp.C:757
    frame #19: 0x0000000108843ae1 libmoose-dbg.0.dylib`MooseApp::run(this=0x00007f92bd84ca18) at MooseApp.C:903
    frame #20: 0x0000000106711373 moose_test-dbg`main(argc=15, argv=0x00007fff594ef7c8) at main.C:31
    frame #21: 0x00007fffb5456235 libdyld.dylib`start + 1
    frame #22: 0x00007fffb5456235 libdyld.dylib`start + 1

# 7/10/18

I don't know wtf the problem with BISON is and why I always get the dman
variable group error:

*** ERROR ***
Error: could not map variable 0 to variable group.

It happens in parallel on BN3X15; on serial the problem just takes
forever. However, with sliding_block it runs great on the MeshModifier both with
distributed and replicated mesh. I remembered this problem...it was due to
setting `quadrature=true` for the `PenetrationLocator`

*** ERROR ***
Material property 'thermal_conductivity', requested by 'thermal_contact_conductivity_master_temp' is not defined on block 9000

[AuxKernels]
  [./gap_value_thermal_contact]
    type = GapValueAux
    boundary = 10
    execute_on = 'INITIAL LINEAR'
    normal_smoothing_distance = 0.1
    order = SECOND
    paired_boundary = 5
    paired_variable = temp
    variable = paired_temp
  [../]
  [./penetration_thermal_contact]
    type = PenetrationAux
    boundary = 10
    execute_on = 'INITIAL LINEAR'
    normal_smoothing_distance = 0.1
    order = SECOND
    paired_boundary = 5
    variable = penetration
  [../]
  [./paired_k_temp_slave_0]
    type = GapValueAux
    boundary = 10
    execute_on = TIMESTEP_BEGIN
    normal_smoothing_distance = 0.1
    order = SECOND
    paired_boundary = 5
    paired_variable = conductivity_temp
    variable = paired_k_temp
  [../]
  [./thermal_contact_conductivity_slave_temp]
    type = MaterialRealAux
    execute_on = TIMESTEP_END
    property = thermal_conductivity
    variable = conductivity_temp
  [../]
  [./thermal_contact_conductivity_master_temp]
    type = MaterialRealAux
    execute_on = TIMESTEP_END
    property = thermal_conductivity
    variable = conductivity_temp
  [../]
[]

[AuxVariables]
  [./penetration]
    order = SECOND
  [../]
  [./paired_temp]
    order = SECOND
  [../]
  [./paired_k_temp]
    order = SECOND
  [../]
  [./conductivity_temp]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./conductivity_temp]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[BCs]
  [./gap_bc_thermal_contact]
    type = GapHeatTransferLWR
    boundary = 10
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_y
    displacements = 'disp_x disp_y'
    gap_distance = penetration
    gap_temp = paired_temp
    order = SECOND
    variable = temp
  [../]
[]

[DiracKernels]
  [./GapHeatPointSourceMaster_thermal_contact]
    type = GapHeatPointSourceMaster
    boundary = 5
    normal_smoothing_distance = 0.1
    slave = 10
    variable = temp
  [../]
[]

[Materials]
  [./gap_value_0]
    type = GapConductanceLWR
    boundary = 10
    contact_pressure = contact_pressure
    gap_distance = penetration
    gap_k = paired_k_temp
    gap_temp = paired_temp
    gas_released = fis_gas_released
    initial_moles = initial_moles
    jump_distance_model = KENNARD
    order = SECOND
    plenum_pressure = plenum_pressure
    roughness_coef = 3.2
    roughness_fuel = 2e-06
    variable = temp
  [../]
[]

[AuxKernels]
  [./gap_value_left_to_right]
    type = GapValueAux
    boundary = leftright
    execute_on = 'INITIAL LINEAR'
    paired_boundary = rightleft
    paired_variable = temp
    variable = paired_temp
  [../]
  [./gap_value_master_left_to_right]
    type = GapValueAux
    boundary = rightleft
    execute_on = 'INITIAL LINEAR'
    paired_boundary = leftright
    paired_variable = temp
    variable = paired_temp
  [../]
  [./penetration_left_to_right]
    type = PenetrationAux
    boundary = leftright
    execute_on = 'INITIAL LINEAR'
    paired_boundary = rightleft
    variable = qpoint_penetration
  [../]
[]

[AuxVariables]
  [./qpoint_penetration]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./paired_temp]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[BCs]
  [./gap_bc_left_to_right]
    type = GapHeatTransfer
    boundary = leftright
    paired_boundary = rightleft
    quadrature = true
    use_displaced_mesh = true
    variable = temp
  [../]
  [./gap_bc_master_left_to_right]
    type = GapHeatTransfer
    boundary = rightleft
    paired_boundary = leftright
    quadrature = true
    use_displaced_mesh = true
    variable = temp
  [../]
[]

[Materials]
  [./left_to_right_gap_value]
    type = GapConductance
    boundary = leftright
    gap_conductivity = 3
    paired_boundary = rightleft
    quadrature = true
    variable = temp
  [../]
  [./left_to_right_gap_value_master]
    type = GapConductance
    boundary = rightleft
    gap_conductivity = 3
    paired_boundary = leftright
    quadrature = true
    variable = temp
  [../]
[]

Time Step 91, time = 2.63043e+07
          old time = 2.53043e+07
                dt = 1e+06
            old dt = 1e+06



Updating geometric search patches

    |residual|_2 of individual variables:
                     disp_x: 325524
                     disp_y: 5184.65
                     temp:   403.448
                     lm:     2.3792e-07

 0 Nonlinear |R| = 3.255659e+05

Time Step 1, time = 2.73043e+07
          old time = 2.63043e+07
                dt = 1e+06
            old dt = 1e+06



Updating geometric search patches

    |residual|_2 of individual variables:
                     disp_x: 3.11983e+06
                     disp_y: 11798.1
                     temp:   4770.2
                     lm:     2.3792e-07

 0 Nonlinear |R| = 3.119857e+06

superlu_dist: dlook_ahead_update.c:117

The lagrange multiplier method appears to be producing singular matrices for the
BISON problems. To wit with MUMPS:
 9 Nonlinear |R| = 1.201071e+03
    0 KSP unpreconditioned resid norm 1.201071323885e+03 true resid norm            nan ||r(i)||/||b||            nan
    0 KSP Residual norm 1.201071323885e+03 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
  Linear solve did not converge due to DIVERGED_PCSETUP_FAILED iterations 0
                 PCSETUP_FAILED due to FACTOR_NUMERIC_ZEROPIVOT

With superlu_dist, it just hangs at that nonlinear iteration.

this=0x0000000115006000
0x0000000115006800

0x0000000113010000

Bending contact PJFNK, mffd_err=1e-8, solve time = 30.5467, lagrange
+----------------+----------------+----------------+----------------+----------------+
| time           | linear         | linear_tot     | nl             | nl_tot         |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   3.500000e+01 |   3.500000e+01 |   1.000000e+01 |   1.000000e+01 |
|   2.000000e+01 |   2.800000e+01 |   6.300000e+01 |   7.000000e+00 |   1.700000e+01 |
|   3.000000e+01 |   3.400000e+01 |   9.700000e+01 |   9.000000e+00 |   2.600000e+01 |
|   4.000000e+01 |   3.500000e+01 |   1.320000e+02 |   8.000000e+00 |   3.400000e+01 |
|   5.000000e+01 |   3.700000e+01 |   1.690000e+02 |   1.000000e+01 |   4.400000e+01 |
|   6.000000e+01 |   4.400000e+01 |   2.130000e+02 |   1.000000e+01 |   5.400000e+01 |
|   7.000000e+01 |   5.500000e+01 |   2.680000e+02 |   1.300000e+01 |   6.700000e+01 |
|   8.000000e+01 |   3.700000e+01 |   3.050000e+02 |   7.000000e+00 |   7.400000e+01 |
|   9.000000e+01 |   1.180000e+02 |   4.230000e+02 |   2.700000e+01 |   1.010000e+02 |
|   1.000000e+02 |   4.000000e+01 |   4.630000e+02 |   7.000000e+00 |   1.080000e+02 |
+----------------+----------------+----------------+----------------+----------------+

Bending contact PJFNK, mffd_err=1e-7, lagrange
+----------------+----------------+----------------+----------------+----------------+
| time           | linear         | linear_tot     | nl             | nl_tot         |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   3.500000e+01 |   3.500000e+01 |   1.000000e+01 |   1.000000e+01 |
|   2.000000e+01 |   2.900000e+01 |   6.400000e+01 |   7.000000e+00 |   1.700000e+01 |
|   3.000000e+01 |   3.400000e+01 |   9.800000e+01 |   9.000000e+00 |   2.600000e+01 |
|   4.000000e+01 |   3.500000e+01 |   1.330000e+02 |   8.000000e+00 |   3.400000e+01 |
|   5.000000e+01 |   3.700000e+01 |   1.700000e+02 |   1.000000e+01 |   4.400000e+01 |
|   6.000000e+01 |   4.400000e+01 |   2.140000e+02 |   1.000000e+01 |   5.400000e+01 |
|   7.000000e+01 |   5.500000e+01 |   2.690000e+02 |   1.300000e+01 |   6.700000e+01 |
|   8.000000e+01 |   3.700000e+01 |   3.060000e+02 |   7.000000e+00 |   7.400000e+01 |
|   9.000000e+01 |   1.500000e+02 |   4.560000e+02 |   3.300000e+01 |   1.070000e+02 |
|   1.000000e+02 |   4.000000e+01 |   4.960000e+02 |   7.000000e+00 |   1.140000e+02 |
+----------------+----------------+----------------+----------------+----------------+

Bending contact NEWTON, solve time = 16.774, lagrange
+----------------+----------------+----------------+----------------+----------------+
| time           | linear         | linear_tot     | nl             | nl_tot         |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   1.400000e+01 |   1.400000e+01 |   1.400000e+01 |   1.400000e+01 |
|   2.000000e+01 |   1.300000e+01 |   2.700000e+01 |   1.300000e+01 |   2.700000e+01 |
|   3.000000e+01 |   1.400000e+01 |   4.100000e+01 |   1.300000e+01 |   4.000000e+01 |
|   4.000000e+01 |   1.300000e+01 |   5.400000e+01 |   1.300000e+01 |   5.300000e+01 |
|   5.000000e+01 |   1.300000e+01 |   6.700000e+01 |   1.300000e+01 |   6.600000e+01 |
|   6.000000e+01 |   3.600000e+01 |   1.030000e+02 |   3.600000e+01 |   1.020000e+02 |
|   7.000000e+01 |   1.600000e+01 |   1.190000e+02 |   1.500000e+01 |   1.170000e+02 |
|   8.000000e+01 |   1.500000e+01 |   1.340000e+02 |   1.400000e+01 |   1.310000e+02 |
|   9.000000e+01 |   2.100000e+01 |   1.550000e+02 |   1.900000e+01 |   1.500000e+02 |
|   1.000000e+02 |   1.500000e+01 |   1.700000e+02 |   1.500000e+01 |   1.650000e+02 |
+----------------+----------------+----------------+----------------+----------------+

Bending contact PJFNK, mffd_err=default, solve time = 12.948, penalty
+----------------+----------------+----------------+----------------+----------------+
| time           | linear         | linear_tot     | nl             | nl_tot         |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   1.200000e+01 |   1.200000e+01 |   4.000000e+00 |   4.000000e+00 |
|   2.000000e+01 |   8.000000e+00 |   2.000000e+01 |   4.000000e+00 |   8.000000e+00 |
|   3.000000e+01 |   4.000000e+00 |   2.400000e+01 |   3.000000e+00 |   1.100000e+01 |
|   4.000000e+01 |   2.700000e+01 |   5.100000e+01 |   7.000000e+00 |   1.800000e+01 |
|   5.000000e+01 |   1.600000e+01 |   6.700000e+01 |   5.000000e+00 |   2.300000e+01 |
|   6.000000e+01 |   9.000000e+00 |   7.600000e+01 |   3.000000e+00 |   2.600000e+01 |
|   7.000000e+01 |   7.600000e+01 |   1.520000e+02 |   1.500000e+01 |   4.100000e+01 |
|   8.000000e+01 |   5.400000e+01 |   2.060000e+02 |   6.000000e+00 |   4.700000e+01 |
|   9.000000e+01 |   9.000000e+00 |   2.150000e+02 |   3.000000e+00 |   5.000000e+01 |
|   1.000000e+02 |   1.000000e+01 |   2.250000e+02 |   4.000000e+00 |   5.400000e+01 |
+----------------+----------------+----------------+----------------+----------------+

Bending contact NEWTON, solve time = 4.89415, penalty
+----------------+----------------+----------------+----------------+----------------+
| time           | linear         | linear_tot     | nl             | nl_tot         |
+----------------+----------------+----------------+----------------+----------------+
|   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e+01 |   4.000000e+00 |   4.000000e+00 |   4.000000e+00 |   4.000000e+00 |
|   2.000000e+01 |   4.000000e+00 |   8.000000e+00 |   4.000000e+00 |   8.000000e+00 |
|   3.000000e+01 |   4.000000e+00 |   1.200000e+01 |   4.000000e+00 |   1.200000e+01 |
|   4.000000e+01 |   4.000000e+00 |   1.600000e+01 |   4.000000e+00 |   1.600000e+01 |
|   5.000000e+01 |   5.000000e+00 |   2.100000e+01 |   5.000000e+00 |   2.100000e+01 |
|   6.000000e+01 |   3.000000e+00 |   2.400000e+01 |   3.000000e+00 |   2.400000e+01 |
|   7.000000e+01 |   3.000000e+00 |   2.700000e+01 |   3.000000e+00 |   2.700000e+01 |
|   8.000000e+01 |   2.000000e+00 |   2.900000e+01 |   2.000000e+00 |   2.900000e+01 |
|   9.000000e+01 |   3.000000e+00 |   3.200000e+01 |   3.000000e+00 |   3.200000e+01 |
|   1.000000e+02 |   4.000000e+00 |   3.600000e+01 |   4.000000e+00 |   3.600000e+01 |
+----------------+----------------+----------------+----------------+----------------+

# 8/2/18

For default template paramters, there must be agreement between the LHS and
RHS. This can never work:

template<..., typename T = 0>

because the LHS is a type and the RHS is a non-type (e.g. it's integral). This
can work if the instantiation is consistent:

template<..., typename T::type = 0>

e.g. if typename T::type yields an integral type, e.g. not a type. So the
qualification, e.g. the nested :: specifier, makes the difference between an
immediate compilation error and an instantiation dependent compilation error.

Note that template<..., typename ...::... = void> will not work and will result
in this compilation error with clang: "expected '(' for function-style cast or
type construction". E.g. no qualification can occur in the default parameter
declaration when the RHS is a type (e.g. void).

# 8/8/18

Contact line search: 58 non-linears
Basic: No convergence
Bt: 53 non-linears

Contact with 1e-4, 3:

Postprocessor Values:
+----------------+------------------+----------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    | tot_nonlinears |
+----------------+------------------+----------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.442325e-02 |   4.000000e+00 |
|   2.000000e-01 |    -0.000000e+00 |   6.000000e+00 |  -1.305756e-03 |   1.000000e+01 |
|   3.000000e-01 |     8.068750e+03 |   7.000000e+00 |   1.613750e-04 |   1.700000e+01 |
|   4.000000e-01 |     9.417125e+03 |   5.000000e+00 |   1.890896e-04 |   2.200000e+01 |
|   5.000000e-01 |     6.944158e+03 |   8.000000e+00 |   1.397839e-04 |   3.000000e+01 |
|   6.000000e-01 |    -0.000000e+00 |   7.000000e+00 |  -2.981951e-03 |   3.700000e+01 |
|   7.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.660042e-02 |   4.100000e+01 |
|   8.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -3.233490e-02 |   4.500000e+01 |
|   9.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -4.770076e-02 |   4.900000e+01 |
|   1.000000e+00 |    -0.000000e+00 |   4.000000e+00 |  -6.027205e-02 |   5.300000e+01 |
+----------------+------------------+----------------+----------------+----------------+

# 8/9/18

Definitely observe oscillations with lagrange multipliers when mesh is made fine.

# 8/16/18

It's impossible to resolve a method pointer to method overloads that differ only
in their function arguments.

# 8/20/18

Here the correct perturbed residual vector:
1.49012e-08
0.
-3.72529e-09
0.
-3.72529e-09
-0.1
0.
-0.1

# 8/23/18

For SFINAE to work, the possible compilation errors **must** involve either
types coming from the template parameters, or it must involve values whose types
involve the template parameters. And for a templated method of a templated
class, the template parameters that matter are that of the method!

# 8/27/18

It doesn't matter that one template overload is a better match than another...if
any matches, it will have it's template arguments deduced and any side effects
will be incited. However, if an error occurs in the body, and the overload is
not as good a match as another, then things will be ok... e.g. the code in that
template will not be evaluated.

# 9/7/18

When compiling a file with c++14 enabled that includes <iostream>, this is error:
/usr/bin/../lib/gcc/x86_64-linux-gnu/4.8/../../../../include/c++/4.8/cstdio:120:11: error: no member named 'gets' in the global namespace
  using ::gets;

This is clearly a header error and happens because the stdc++ library (both lib
and headers) is incomplete. The c++ library however is c++14 complete, and so if
we include the c++ headers then we don't get any header error. But then we do
get a linking error because we don't provide the definitions for the
declarations in the c++ headers. Once we link to the c++ library, everything
works. Wait, let's qualify that last sentence. This line compiles:

clang++-3.5 -std=c++1y test.cpp -L$HOME/libcxx-3.5.0.src/installed/lib -lc++

ldd -r reveals:

	linux-vdso.so.1 =>  (0x00007ffd36f2c000)
	libc++.so.1 => /home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1 (0x00007fa0b3702000)
	libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fa0b33fe000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fa0b30f8000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fa0b2ee2000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fa0b2b19000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fa0b28fb000)
	librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007fa0b26f3000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fa0b3a9e000)

If we run ldd -r on the c++ library:

	linux-vdso.so.1 =>  (0x00007ffd51b28000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fbe8f162000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fbe8ed99000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fbe8ea93000)
	librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007fbe8e88b000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fbe8e675000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fbe8f71c000)
undefined symbol: __cxa_pure_virtual	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: _ZTVN10__cxxabiv117__class_type_infoE	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: _ZTVN10__cxxabiv121__vmi_class_type_infoE	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: _ZTVN10__cxxabiv120__si_class_type_infoE	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __gxx_personality_v0	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_end_catch	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_allocate_exception	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_guard_release	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_begin_catch	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_rethrow	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_throw	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_guard_abort	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_free_exception	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)
undefined symbol: __cxa_guard_acquire	(/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so.1)

We run into all these undefined references. These references either have to be
provided by libc++abi or libstdc++. So for example if we compile like this:

clang++-3.5 -std=c++1y -stdlib=libc++ test.cpp
-L$HOME/libcxx-3.5.0.src/installed/lib -lc++

then we get linker errors like the following:

/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `__cxa_guard_release'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `__cxa_pure_virtual'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `__cxa_rethrow'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `vtable for __cxxabiv1::__si_class_type_info'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `__cxa_guard_abort'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `vtable for __cxxabiv1::__class_type_info'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `vtable for __cxxabiv1::__vmi_class_type_info'
/home/lindad/libcxx-3.5.0.src/installed/lib/libc++.so: undefined reference to `__cxa_guard_acquire'

This happens because the linker no longer appends -lstdc++ which would have
provided those symbols. We can restore the compilation, however, if we run:

clang++-3.5 -std=c++1y -stdlib=libc++ test.cpp -L$HOME/libcxx-3.5.0.src/installed/lib -L$HOME/libcxxabi-3.5.0.src/installed/lib -lc++ -lc++abi

**and** make sure that libc++abi is in LD_LIBRARY_PATH come run time.

So, how does clang decide what c++ library to append to its link line?

Potentially important llvm flags:

CMAKE_CXX_COMPILER
CMAKE_CXX_FLAGS
CMAKE_C_COMPILER
CMAKE_C_FLAGS
CMAKE_EXE_LINKER_FLAGS
CMAKE_SHARED_LINKER_FLAGS
CMAKE_STATIC_LINKER_FLAGS
GCC_INSTALL_PREFIX
C_INCLUDE_DIRS

Potentiall important libc++ flags:
LIBCXX_CXX_ABI:STRING

Possible values for LIBCXX_CXX_ABI:libstdc++,libcxxabi,libcxxrt,none

# 9/25/18

300x300 ad-simple-diffusion (100): 6.444 s
300x300 ad-simple-diffusion (4): 3.774 s
300x300 simple-diffusion: 3.648 s

500x500 ad-simple-diffusion (100): 18.788
500x500 ad-simple-diffusion (4): 11.380
500x500 simple-diffusion: 11.135
500x500 ad-simple-diffusion (100) with MooseVariableFE opt: 13.361

# 10/3/18

  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 1.)  (2, 0.)
row 1: (0, 0.)  (1, 1.)  (2, 0.)  (3, 0.)
row 2: (0, -2.)  (2, 4.)  (4, -2.)
row 3: (0, 0.5)  (1, -1.)  (2, 0.)  (3, 2.)  (4, -0.5)  (5, -1.)
row 4: (2, 0.)  (4, 1.)
row 5: (2, 0.)  (3, 0.)  (4, 0.)  (5, 1.)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 1.)
row 1: (1, 1.)
row 2: (0, -2.)  (2, 4.)  (4, -2.)
row 3: (0, 0.5)  (1, -1.)  (2, -1.)  (3, 2.)  (4, 0.5)  (5, -1.)
row 4: (4, 1.)
row 5: (5, 1.)

Time Step 2, time = 2
                dt = 1

 0 Nonlinear |R| = 1.000000e+00
  ---------- Testing Jacobian -------------
  ||J - Jfd||_F/||J||_F = 0.411113, ||J - Jfd||_F = 2.44949
  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 1.)  (2, 0.)
row 1: (0, 0.)  (1, 1.)  (2, 0.)  (3, 0.)
row 2: (0, -2.)  (2, 4.)  (4, -2.)
row 3: (0, -0.5)  (1, -1.)  (2, 1.)  (3, 2.)  (4, -0.5)  (5, -1.)
row 4: (2, 0.)  (4, 1.)
row 5: (2, 0.)  (3, 0.)  (4, 0.)  (5, 1.)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 1.)
row 1: (1, 1.)
row 2: (0, -2.)  (2, 4.)  (4, -2.)
row 3: (0, 0.5)  (1, -1.)  (2, -1.)  (3, 2.)  (4, 0.5)  (5, -1.)
row 4: (4, 1.)
row 5: (5, 1.)
      0 Linear |R| = 1.000000e+00
      1 Linear |R| = 0.000000e+00
 1 Nonlinear |R| = 0.000000e+00
 Solve Converged!

 0 Nonlinear |R| = 1.000000e+00
  ---------- Testing Jacobian -------------
  ||J - Jfd||_F/||J||_F = 0.31334, ||J - Jfd||_F = 1.83712
  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 1.)  (2, 0.)
row 1: (0, 0.)  (1, 1.)  (2, 0.)  (3, 0.)
row 2: (0, -2.)  (2, 4.)  (4, -2.)
row 3: (0, -0.25)  (1, -1.)  (2, 0.5)  (3, 2.)  (4, -0.25)  (5, -1.)
row 4: (2, 0.)  (4, 1.)
row 5: (2, 0.)  (3, 0.)  (4, 0.)  (5, 1.)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 1.)
row 1: (1, 1.)
row 2: (0, -2.)  (2, 4.)  (4, -2.)
row 3: (0, 0.5)  (1, -1.)  (2, -1.)  (3, 2.)  (4, 0.5)  (5, -1.)
row 4: (4, 1.)
row 5: (5, 1.)
      0 Linear |R| = 1.000000e+00
      1 Linear |R| = 0.000000e+00
 1 Nonlinear |R| = 0.000000e+00
 Solve Converged!


# 11/21/18

  Hand-coded Jacobian ----------
Mat Object: () 1 MPI processes
  type: seqaij
row 0: (0, 1.)  (1, 0.)  (2, -1.)  (3, 0.)
row 1: (0, 0.25)  (1, 0.5)  (2, -0.25)  (3, -0.5)
row 2: (0, -1.)  (1, 0.)  (2, 1.)  (3, 0.)
row 3: (0, -0.25)  (1, -0.5)  (2, 0.25)  (3, 0.5)
  Finite difference Jacobian ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0: (0, 1.)  (2, -1.)
row 1: (0, -0.5)  (1, 0.5)  (2, 0.5)  (3, -0.5)
row 2: (0, -1.)  (2, 1.)
row 3: (0, 0.5)  (1, -0.5)  (2, -0.5)  (3, 0.5)
  Hand-coded minus finite-difference Jacobian with tolerance 1e-05 ----------
Mat Object: 1 MPI processes
  type: seqaij
row 0:
row 1: (0, 0.75)  (2, -0.75)
row 2:
row 3: (0, -0.75)  (2, 0.75)

In [1]: JxWp = 0.99999998174987925

In [2]: JxW = 1

In [3]: gradup = 1.0000000182501212

In [4]: gradu = 1

In [5]: gradtestp = -0.5000000091250606

In [6]: gradtest = -.5

In [7]: wscale = 2.7397079002971876e+07

In [8]: wscale * (JxWp - JxW)
Out[8]: -0.4999999999457628

In [9]: wscale * (gradup - gradu)
Out[9]: 0.5000000121125099

In [10]: wscale * (gradtestp - gradtest)
Out[10]: -0.25000000605625494

# 11/28/18

Regular time: 3.33 s
Updated AD: 7.96 s
Old AD: 22.05 s
Updated AD with right size: 5.56 s
With loop optimizations: 5.772s
With loop optimizations right size: 4.198s

# 11/30/18

PJFNK:
Num steps = 19
Total nonlinear = 85

NEWTON:
Num_steps = 19
Total nonlinear = 82

# 12/5/18

With PETSc scaling:

NEWTON without SUPG:
Num steps = 21
End time = 8.45688e-05
nonlinear = 84
linear = 1721

With PETSc scaling

NEWTON without SUPG:
Num steps = 15
End time = 8.675e-5
nonlinear = 63
linear = 1402

# 1/23/18

MOOSE 5f250a2486, libmesh after-libmesh-update: 154 seconds
MOOSE 5f250a2486, libmesh after-vector-data: 163 seconds
MOOSE 5f250a2486, libmesh before-libmesh-data: 136 seconds
MOOSE 5f250a2486, libmesh before-vector-data: 163 seconds
MOOSE 5f250a2486, libmesh 33da95a18: 154 seconds
MOOSE 5f250a2486, libmesh 710a816fd: 156 seconds
MOOSE 5f250a2486, libmesh 0dc5c33ba: 158.5 seconds
MOOSE 5f250a2486, libmesh 9bcf8a1ca: 157 seconds
MOOSE 0ea1a2d93a, libmesh ab2cf9725: 135 seconds
MOOSE 5f250a2486, libmesh ab2cf9725: 138 seconds

# 1/24/18

Failed golem history (random):

lindad@Alexanders-Pro:~/projects/golem/test/tests/crash_test(master)$ mpirun -np 11 ~/projects/golem/golem-dbg -i test.i

Framework Information:
MOOSE Version:           git commit bd8d51028f on 2019-01-21
LibMesh Version:         8ddbf2885ac645d18fb238a4c17cf65e5d8df9c6
PETSc Version:           3.10.3
Current Time:            Thu Jan 24 13:58:41 2019
Executable Timestamp:    Thu Jan 24 13:54:34 2019

Parallelism:
  Num Processors:          11
  Num Threads:             1

Mesh:
  Parallel Type:           replicated
  Mesh Dimension:          3
  Spatial Dimension:       3
  Nodes:
    Total:                 7547
    Local:                 808
  Elems:
    Total:                 39126
    Local:                 3615
  Num Subdomains:          4
  Num Partitions:          11
  Partitioner:             metis

Nonlinear System:
  Num DOFs:                37735
  Num Local DOFs:          4040
  Variables:               { "disp_x" "disp_y" "disp_z" "pore_pressure" "temperature" }
  Finite Element Types:    "LAGRANGE"
  Approximation Orders:    "FIRST"

Auxiliary System:
  Num DOFs:                234756
  Num Local DOFs:          21690
  Variables:               { "stress_xx" "stress_yy" "stress_zz" "stress_xy" "stress_xz" "stress_yz"
                             }
  Finite Element Types:    "MONOMIAL"
  Approximation Orders:    "CONSTANT"

Execution Information:
  Executioner:             Transient
  TimeStepper:             TimeSequenceStepper
  Solver Mode:             NEWTON


Time Step 0, time = 0
                dt = 0


Time Step 1, time = 1
                dt = 1

 0 Nonlinear |R| = 9.727970e+05
  0 SNES Function norm 9.727969590319e+05
      0 Linear |R| = 9.727970e+05
      1 Linear |R| = 8.534746e+01
 1 Nonlinear |R| = 8.534746e+01
  1 SNES Function norm 8.534746454278e+01
      0 Linear |R| = 8.534746e+01
      1 Linear |R| = 7.218560e-03
 2 Nonlinear |R| = 7.218560e-03
  2 SNES Function norm 7.218559876268e-03
Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 2
 Solve Converged!


 ----------------------------------------------------------------------------
| Reference count information                                                |
 ----------------------------------------------------------------------------
| N7libMesh10FEAbstractE reference count information:
|  Creations:    18190
|  Destructions: 18190
| N7libMesh10Parameters5ValueE reference count information:
|  Creations:    30363
|  Destructions: 30363
| N7libMesh12SparseMatrixIdEE reference count information:
|  Creations:    8
|  Destructions: 8
| N7libMesh13NumericVectorIdEE reference count information:
|  Creations:    73
|  Destructions: 73
| N7libMesh15EquationSystemsE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh15GhostingFunctorE reference count information:
|  Creations:    10
|  Destructions: 10
| N7libMesh15NonlinearSolverIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh4ElemE reference count information:
|  Creations:    159780
|  Destructions: 159780
| N7libMesh4NodeE reference count information:
|  Creations:    15094
|  Destructions: 15094
| N7libMesh5QBaseE reference count information:
|  Creations:    54261
|  Destructions: 54261
| N7libMesh6DofMapE reference count information:
|  Creations:    4
|  Destructions: 4
| N7libMesh6SystemE reference count information:
|  Creations:    4
|  Destructions: 4
| N7libMesh9DofObjectE reference count information:
|  Creations:    174874
|  Destructions: 174874
 ----------------------------------------------------------------------------
lindad@Alexanders-Pro:~/projects/golem/test/tests/crash_test(master)$ mpirun -np 11 ~/projects/golem/golem-dbg -i test.i

Framework Information:
MOOSE Version:           git commit bd8d51028f on 2019-01-21
LibMesh Version:         8ddbf2885ac645d18fb238a4c17cf65e5d8df9c6
PETSc Version:           3.10.3
Current Time:            Thu Jan 24 14:02:18 2019
Executable Timestamp:    Thu Jan 24 13:54:34 2019

Parallelism:
  Num Processors:          11
  Num Threads:             1

Mesh:
  Parallel Type:           replicated
  Mesh Dimension:          3
  Spatial Dimension:       3
  Nodes:
    Total:                 7547
    Local:                 808
  Elems:
    Total:                 39126
    Local:                 3615
  Num Subdomains:          4
  Num Partitions:          11
  Partitioner:             metis

Nonlinear System:
  Num DOFs:                37735
  Num Local DOFs:          4040
  Variables:               { "disp_x" "disp_y" "disp_z" "pore_pressure" "temperature" }
  Finite Element Types:    "LAGRANGE"
  Approximation Orders:    "FIRST"

Auxiliary System:
  Num DOFs:                234756
  Num Local DOFs:          21690
  Variables:               { "stress_xx" "stress_yy" "stress_zz" "stress_xy" "stress_xz" "stress_yz"
                             }
  Finite Element Types:    "MONOMIAL"
  Approximation Orders:    "CONSTANT"

Execution Information:
  Executioner:             Transient
  TimeStepper:             TimeSequenceStepper
  Solver Mode:             NEWTON


Time Step 0, time = 0
                dt = 0


Time Step 1, time = 1
                dt = 1

 0 Nonlinear |R| = 9.727970e+05
  0 SNES Function norm 9.727969590319e+05
      0 Linear |R| = 9.727970e+05
      1 Linear |R| = 8.534746e+01
 1 Nonlinear |R| = 8.534746e+01
  1 SNES Function norm 8.534746454279e+01
      0 Linear |R| = 8.534746e+01
      1 Linear |R| = 7.218560e-03
 2 Nonlinear |R| = 7.218560e-03
  2 SNES Function norm 7.218559875459e-03
Nonlinear solve converged due to CONVERGED_FNORM_ABS iterations 2
 Solve Converged!


 ----------------------------------------------------------------------------
| Reference count information                                                |
 ----------------------------------------------------------------------------
| N7libMesh10FEAbstractE reference count information:
|  Creations:    18190
|  Destructions: 18190
| N7libMesh10Parameters5ValueE reference count information:
|  Creations:    30363
|  Destructions: 30363
| N7libMesh12SparseMatrixIdEE reference count information:
|  Creations:    8
|  Destructions: 8
| N7libMesh13NumericVectorIdEE reference count information:
|  Creations:    73
|  Destructions: 73
| N7libMesh15EquationSystemsE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh15GhostingFunctorE reference count information:
|  Creations:    10
|  Destructions: 10
| N7libMesh15NonlinearSolverIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh4ElemE reference count information:
|  Creations:    159780
|  Destructions: 159780
| N7libMesh4NodeE reference count information:
|  Creations:    15094
|  Destructions: 15094
| N7libMesh5QBaseE reference count information:
|  Creations:    54261
|  Destructions: 54261
| N7libMesh6DofMapE reference count information:
|  Creations:    4
|  Destructions: 4
| N7libMesh6SystemE reference count information:
|  Creations:    4
|  Destructions: 4
| N7libMesh9DofObjectE reference count information:
|  Creations:    174874
|  Destructions: 174874
 ----------------------------------------------------------------------------
lindad@Alexanders-Pro:~/projects/golem/test/tests/crash_test(master)$
lindad@Alexanders-Pro:~/projects/golem/test/tests/crash_test(master)$ mpirun -np 11 ~/projects/golem/golem-dbg -i test.i

Framework Information:
MOOSE Version:           git commit bd8d51028f on 2019-01-21
LibMesh Version:         8ddbf2885ac645d18fb238a4c17cf65e5d8df9c6
PETSc Version:           3.10.3
Current Time:            Thu Jan 24 14:06:13 2019
Executable Timestamp:    Thu Jan 24 13:54:34 2019

Parallelism:
  Num Processors:          11
  Num Threads:             1

Mesh:
  Parallel Type:           replicated
  Mesh Dimension:          3
  Spatial Dimension:       3
  Nodes:
    Total:                 7547
    Local:                 844
  Elems:
    Total:                 39126
    Local:                 3629
  Num Subdomains:          4
  Num Partitions:          11
  Partitioner:             metis

Nonlinear System:
  Num DOFs:                37735
  Num Local DOFs:          4220
  Variables:               { "disp_x" "disp_y" "disp_z" "pore_pressure" "temperature" }
  Finite Element Types:    "LAGRANGE"
  Approximation Orders:    "FIRST"

Auxiliary System:
  Num DOFs:                234756
  Num Local DOFs:          21774
  Variables:               { "stress_xx" "stress_yy" "stress_zz" "stress_xy" "stress_xz" "stress_yz"
                             }
  Finite Element Types:    "MONOMIAL"
  Approximation Orders:    "CONSTANT"

Execution Information:
  Executioner:             Transient
  TimeStepper:             TimeSequenceStepper
  Solver Mode:             NEWTON

Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 4040
v.local_size() = 4220

Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3510
v.local_size() = 3365


Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3515
v.local_size() = 3650


Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3995
v.local_size() = 3775


[3] ../src/numerics/petsc_vector.CAssertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3305
v.local_size() = 3360


[4] ../src/numerics/petsc_vector.CAssertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3170
v.local_size() = 2885


[5] ../src/numerics/petsc_vector.C, line 567, compiled Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3360
v.local_size() = 3530


[6] ../src/numerics/petsc_vector.C, line 567, compiled nodate at Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3575
v.local_size() = 3535


[7] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime
Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3820
v.local_size() = 3595


[8] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime
Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 2380
v.local_size() = 3260


[9] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime
Assertion `this->local_size() == v.local_size()' failed.
this->local_size() = 3065
v.local_size() = 2560


[10] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime

[0] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime
[1] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime
[2] ../src/numerics/petsc_vector.C, line 567, compiled nodate at notime
, line 567, compiled nodate at notime
, line 567, compiled nodate at notime
nodate at notime
notime
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 6
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 1
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 2
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 3
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 4
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 5
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 7
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 8
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 9
application called MPI_Abort(MPI_COMM_WORLD, 1) - process 10
lindad@Alexanders-Pro:~/projects/golem/test/tests/crash_test(master)$ ls
