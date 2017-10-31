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

ComputeNodalStuffThreads:
  FEProblemBase::reinitNode
    *System*::reinitNode
	  if (var->isNodal())
      {
        var->reinitNode();
        var->computeNodalValues();
      }


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
