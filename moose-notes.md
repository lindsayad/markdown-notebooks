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