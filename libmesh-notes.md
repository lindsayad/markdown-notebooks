# 10/17/17

Ok, so `newton_iterate` in newton_solver.C comes from `_system.solution`. Ok for
just the vector FE problem with a `grid_size` of 6, the number of degrees of
freedom is 120. However, the size of the solution that is passed to exodus
is 338. This is consistent with our stupid nodal definition; we have 169 nodes
in the mesh and if we split our vector variable into two components, then we
have 169*2=338 values.

So in `newton_solver.C` `_system` is `(CurlCurlSystem *) $21 =
0x000000010d009670`. And then in `equation_systems.C`, `system` is
`(CurlCurlSystem *) $23 = 0x000000010d009670`. Ok, so they're the same object
which is good.

Ok, it's clear to me that operator()(int index) doesn't work in lldb. For some
reason the solution vector looks good when investigating it in
`libMesh::EquationSystems::build_parallel_solution_vector` but looks like
garbage from inside of `newton_solver.C`. Perhaps there's some intermediate
steps to that vector, I don't know.

# 10/19/17

Ok so FEMContext::interior_curl is only implemented for Gradient, which when
libmesh is compiled with complex enabled has the signature:
```
template void FEMContext::interior_curl<Gradient>(unsigned int, unsigned int, Gradient &) const;
template void FEMContext::interior_curl<Gradient>(unsigned int, unsigned int, std::complex<double> &) const;
```
thus if interior_curl is called with a RealGradient argument, the linker will fail.
