Moltres
=======

# 6/14/17

Convergence is a word that often has two meanings in the context of finite
element modelling. There is:

- Convergence of the numerical solution to the analytic or true solution. This
  can be associated with order of accuracy testing etc.
- Convergence of the **solution method**, e.g. whether your non-linear or linear
  solver is able to solve the provided system of equations.

In my experience, mesh refinement has opposite effects on the two types of
convergence:

- It improves the numerical solution relative to the analytic solution (this is
  a fundamental principle of the FEM)
- It degrades convergence of the solution method, e.g. it makes it more
  difficult to solve the system of non-linear or linear equations.

These opposing characteristics are fascinating to me.

# 6/20/17

Observed something very interesting today. I've observed a fair number of times
over the years that solution of my governing equations often becomes quite slow
after I've resolved the main transients and am approaching steady-state. I saw
today that in that final phase of solution, not having a line search can
**greatly** increase solution speed.

2D axisymmetric simulation with moderator heating took 5209 seconds and 293 time
steps with PJFNK and line\_search on. On the other hand, turning line search off
resulted in **169** seconds and only 80 time steps! Wow, this could be a very
helpful discovery.

Neutronics only simulation with parametric variation of new control gain
parameter:

- gain = 1e-4; peak power density = 10.42 W/cm^3
- gain = 1e-5; peak power density = 14.20 W/cm^3

# 6/21/17

There are a couple of flow assumptions one can make:

- Plug flow: From [wikipedia](https://en.wikipedia.org/wiki/Plug_flow): For
  flows in pipes, if flow is turbulent then the laminar [or viscous] sublayer
  caused by the pipe wall is so thin that it is negligible. Plug flow will be
  achieved if the sublayer thickness is much less than the pipe diameter.
- Laminar pipe flow: laminar flow with a fully developed boundary layer. In a
  circular pipe, the velocity profile is parabolic.

# 6/22/17

There must be something wrong with the Navier Stokes Jacobians because the
difficult problem converges with FDP as the preconditioning matrix but not with
SMP. This is true for both `solve_type = NEWTON` and `solve_type = PJFNK`.

# 6/28/17

Single channel k:

- with vacuum bcs: .8269
- with refelctive bcs: 1.5622

# 6/30/17

Gavin ran the same `3d_auto_diff_rho.i` simulation that I ran back on April 28
and he saw temperatures get into the thousands of kelvin by 3 seconds. The only
difference was that his mesh was 1 cm **shorter**. When I ran the same input
file back in April, I never got above 938 Kelvin and I simulated out to 5 seconds.

# 7/3/17

Summary of control experiments on 2D axisymmetric simulation:

- With controller gain executed during the non-linear solve (gains up to 1e-4),
  have to use PJFNK with an iterative preconditioner to get good convergence
  during the initial time steps
  -  During the late time steps, convergence poor both with and without line
     search
- With controller gain executed just on `timestep_end`, the gain can't be more
  than 3e-7 or else the control response becomes oscillatory and unbounded
  - With the controller, the max power density is 46
  - Without the controller, the max power density is 51, so the controller at
    that low gain level doesn't have much effect
  - Can use all the best convergence features, however: direct preconditioning,
    no line search, and NEWTON

# 7/5/17

On a 2d-axisymmetric mesh scaled down by about .4:

FunctionMaterial: k = 1.006752
GenericMoltresMaterial: k = .239

The latter is much more reasonable!!! How the hell did we increase k by so
much???

Alright, setting fissxs and nsf equal to zero fixed this issue.

# 7/6/17

Eigenvalue for `3d_eigen_function_materials.i`: .997073603

Nice! So hopefully there's a decent chance that the long simulation works successfully.

# 7/7/17

For the 2D axisymmetric case, the eigenvalue calculation predicts a
super-critical reactor at a scale of .97. However, the transient case is
subcritical at this scale. The transient case is not super-critical until
.99. The corresponding k-eigenvalue from the eigen simulation for this scale is:
1.0051031918

Ok, with 3d case with scale factor of .98, k = 1.0026059425. Scale of .99, k =
1.0057546381. Let's use that value then.

I think we have a converged 3d case...

- 160 cpus
- active time = 8095.03 seconds = 2.53 hours

# 7/18/17

So far I haven't run across a single case where the following petsc options
cannot converge to steady state:
```
-pc_type lu -ksp_type gmres -snes_mf_operator -snes_linesearch_type basic
```
Sweet!

# 7/19/17

Trying to understand the overall design of the MSRE. The fuel salt is cooled by
a coolant salt in a heat exchanger. Then the coolant salt is cooled by an
air-cooled radiator. The whole experiment was designed to operate at 10 MW with
10 MW transferred between fuel salt and coolant salt and 10 MW transferred
between coolant salt and the air.

# 7/24/17

- Tf = 922, Tm = 922, k = 1.0052768822942
- Tf = 972, Tm = 922, k = 1.0009144326448
- Tf = 922, Tm = 972, k = 1.0048156828401
- Tf = 972, Tm = 972, k = 1.0004220217591

- alphaf = -8.671
- alpham = -0.913
- alpha total calculated = -9.584
- alpha total simulated = -9.655

# 8/1/17

### Channel flow cases

Correct:
- Laplace form, natural bcs for both components
- traction form, no bc for perpendicular, natural for parallel
- traction form, no bcs for both components
- laplace form, no bcs for both components
- laplace form, no bc for perpendicular, natural for parallel

Incorrect
- traction form, natural for both components

# 8/3/17

It looks to me like the leaks for gdb-8 are roughly the same between
`moltres-dbg` and some test program like `hello`. However, the gdb-8 leaks
appear to be quite different compared to gdb-7. For example the dominant leak in
gdb-8 appears to be in `libiberty/cp-demangle.c:d_growable_string_resize()`
whereas with gdb-7 it was in `record_thread`

Here's a headache of dealing with a large stack and pointers: if we set pointer
after pointer after pointer equal to each other and then I free the memory
somewhere far up the stack, I risk the fact that somewhere later I'm going to
get an invalid read or write.

`dgs.buf` is a pointer to `char`.

`cplus_demangle_print` returns `dgs.buf`. Its return value is a pointer to
`char`. A pointer points to an address in memory. Its caller,
`cp_comp_to_string` just returns the same return value, so it also returns a
pointer to char.

# 8/8/17

It's encouraging that the Steady simulation results and Transient
simulation marching to steady-state results match in the fluid flow simulations
that I'm doing.

# 8/11/17

Running cases with openfoam...zero gradient boundary conditions for k and
epsilon result in much different results than using the default tutorial bc
settings. :-(

# 11/27/17

```
- num_elems: 336, r: .8, z: 20, fuel: 922, group1: 4.25e-5
- num_elems: 1960, r: .4, z: 10, fuel: 927.6383, group1: 4.618033e6
- num_elems: 9338, r: .2, z: 5, fuel: 928.9675, group1: 5.747188e6
- num_elems: 40474, r: .1, z: 2.5, fuel: 929.3437, group1: 6.056901e6
- num_elems: 168000, r: .05, z: 1.25, fuel: 929.4508, group1: 6.144888e6
```
