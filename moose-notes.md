Print position in vector in LLDB:
```
p elem_num_map.__begin_[8297]
```

How to access the multi app FEProblem's time from the master app's FEProblem:
```
p _multi_apps._active_objects.__begin_[0].__begin_[0].__ptr_->_apps.__begin_[0].__ptr_->_executioner.__ptr_->_fe_problem._time
```

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
