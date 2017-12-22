# 12/8/17

The preconditioning object `pc` has both a `mat` and a `pmat` member. `mat`
refers to the Jacobian operator (matrix if not using matrix-free methods).

As evidenced in PCApplyBAorAB, the Jacobian matrix multiplication with the
vector is explicitly done! And then the preconditioner is applied in the case of
left preconditioning.
