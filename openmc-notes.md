# 4/23/18

In tally.F90, we have these levels:

```
- TALLY_LOOP
  call filters(i_filt) % obj % get_all_bins(p, t % estimator, &
               filter_matches(i_filt))
  - FILTER_LOOP
    - SCORE_LOOP
    call expand_and_score(p, t, score_index, filter_index, score_bin, &
               score, k)
```

Backtrace for zernike tally:

```
* thread #1, queue = 'com.apple.main-thread', stop reason = breakpoint 1.1
  * frame #0: 0x000000010c1229cf libopenmc.dylib`__tally_filter_zernike_MOD_get_all_bins at tally_filter_zernike.F90:69
    frame #1: 0x000000010c0d8b2a libopenmc.dylib`__tally_MOD_score_collision_tally at tally.F90:2956
    frame #2: 0x000000010c0b42af libopenmc.dylib`__tracking_MOD_transport at tracking.F90:219
    frame #3: 0x000000010c090d25 libopenmc.dylib`openmc_next_batch._omp_fn.4 at simulation.F90:124
    frame #4: 0x000000010b73638e libomp.dylib`GOMP_parallel + 158
    frame #5: 0x000000010c08ef2e libopenmc.dylib`openmc_next_batch at simulation.F90:124
    frame #6: 0x000000010c08effe libopenmc.dylib`openmc_run at simulation.F90:69
    frame #7: 0x00000001000d6f16 libokapi-dbg.0.dylib`OpenMCTimeStepper::step(this=0x000000010f018358) at OpenMCTimeStepper.C:42
    frame #8: 0x00000001007e35ce libmoose-dbg.0.dylib`Transient::solveStep(this=0x000000010f036628, input_dt=-1) at Transient.C:552
    frame #9: 0x00000001007e2d0a libmoose-dbg.0.dylib`Transient::takeStep(this=0x000000010f036628, input_dt=-1) at Transient.C:465
    frame #10: 0x00000001007e2543 libmoose-dbg.0.dylib`Transient::execute(this=0x000000010f036628) at Transient.C:371
    frame #11: 0x00000001019151e3 libmoose-dbg.0.dylib`MooseApp::executeExecutioner(this=0x000000010f83d618) at MooseApp.C:710
    frame #12: 0x00000001019159db libmoose-dbg.0.dylib`MooseApp::run(this=0x000000010f83d618) at MooseApp.C:813
    frame #13: 0x00000001000017e3 okapi-dbg`main(argc=3, argv=0x00007fff5fbfea10) at main.C:24
    frame #14: 0x00007fff8f409235 libdyld.dylib`start + 1
```

Later in the `score_collision_tally` routine, there is a call to
`score_general`, which is a function pointer. For continuous energy simulations
it points to `score_general_ce`.

In April's implementation, the Zernike polynomials were analog tallies. In
Paul's, they are collision tallies.

In April's openmc:

First d_collision:

0.8587321037088177 -> this must be something else

Next simulation:

6.625716130346567

Next simulation:

6.625716130346567

**In master openmc**

After source sampling:

6.625716130346567

So we are perfectly matched up to this point!

# 4/25/18

A fixture that depends on another fixture will be torn down before the dependency.

# 5/24/18

Zernike moments are calculated in ZernikeFilter::get_all_bins and stored in `filter_weights`.
This is called from `score_analog_tally_ce`.
