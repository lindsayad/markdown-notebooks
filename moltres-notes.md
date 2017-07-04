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
file back in April, I never got about 938 Kelvin and I simulated out to 5 seconds.

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