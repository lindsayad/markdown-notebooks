To run petsc tests:

make -j8 -f gmakefile test search=snes_tutorials-ex99_%

# 12/8/17

The preconditioning object `pc` has both a `mat` and a `pmat` member. `mat`
refers to the Jacobian operator (matrix if not using matrix-free methods).

As evidenced in PCApplyBAorAB, the Jacobian matrix multiplication with the
vector is explicitly done! And then the preconditioner is applied in the case of
left preconditioning.

# 1/10/18

Ok, the key petsc example to look at is indeed example 7 in
examples/tutorials. This demonstrates how to set a different function for
matrix-free differencing. Ugh, so much effort and it's in front of my nose.

# 1/11/18

Ok, the function I set will be called with these arguments:
(*ctx->func)(ctx->funcctx,U,F)

The function context is set with the last argument to `MatMFFDSetFunction`. For us this should be the PetscNonlinearSolver.

# 5/17/18

ksp: krylov-subspace projection

symmetric, positive, definite: CG (conjugate gradient)
otherwise: GMRES
