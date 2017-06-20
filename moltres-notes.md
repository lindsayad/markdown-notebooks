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
