Air Force postdoc research fellowship deadline: May 1st

Here are my three threads I'm working concurrently on in the hope that one/some/all of them lead me to my goal of interfacing two physics domains (or as Yaqi I think aptly calls them InteriorNodalBC and InteriorIntegratedBC):

https://groups.google.com/forum/#!topic/moose-users/011fudvdHng (MeshModifiers)
https://groups.google.com/forum/#!topic/moose-users/u3qZrEXQvSg (Boundary ID in DGKernel)
https://groups.google.com/forum/#!topic/moose-users/iALxx76WhnU (Multiple-block, multiple-variable interface test)

More moose-user threads of interest:

Trying to get slope limiting into DG: https://groups.google.com/forum/#!search/yaqi$20convection/moose-users/EYaTNUsSkss/VZzZ-wZEdoQJ
Solving scalar pure advection problem, including DG and FV: https://groups.google.com/forum/#!searchin/moose-users/DG$20finite$20volume/moose-users/p6h7PqZNOrA/TZfG_cfUAwAJ

Matthew Hopkins: PIC & multiphysics plasma simulation guy at Sandia

Universal File Datasets Summaries:
http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-file-datasets-summary
https://docs.plm.automation.siemens.com/tdoc/nx/10/nx_help/#uid:index_advanced:xid602249:id625716:id625821

Really good articles on the moose-user list:

preconditioners (authored by Mariana Rodriguez)

Joule heating is working. It simply appears that currently the energy losses due to ionzation and elastic collisions outweigh the energy gains due to Joule heating.

Going to really start investigating what prevents taking long time steps. With no off diag jacobian elements for electrons and ions and off diag element for potential turned off, and with the following variable scalings:

EMACS COMMANDS:
just move the point to the end of any sexp and press C-xC-e to execute just that sexp. Usually it's not necessary to reload the whole file if you're just changing a line or two.

If you take care when choosing what keys to rebind, these will never conflict with keybindings of modes. When there is such a conflict, however, the more specific mode keymap will take precedence, and your global keybinding will be temporarily shadowed.

    (define-key text-mode-map (kbd "'") 'maybe-open-apostrophe)

No need to guess. If you want to find out the current active keymap in a given situation, like in the gnus group buffer, use this command:

    M-: (mapcar (lambda(x)(car(rassq x minor-mode-map-alist)))(current-minor-mode-maps))

To get the keybinding of something, do C-h k

Also, major and minor mode keybindings override global key bindings. Also, in order for changes in el library sources to take effect, they must be recompiled

em: 1e-11
Arp: 1e-6
potential: 1e4

the step size at 1e-7 seconds is 1.6e-9

em: 1e-6
Arp: 1e-6
potential: 1e4

the step size at 1e-7 seconds is 1.6e-9

Currently not able to even get to 1e-7 seconds with NEWTON. Perhaps my Jacobian is wrong. But it could also be right.

Best sim fail: somewhere in the microsecond range. Specs for this best simulation: nx = 4000, nl_rel_tol = 1e-2, hypre boomeramg, no scaling of any variables. Looks like fail time was: 3.0207e-6 seconds. I think a kink between .01 and .015 meters is what dooms the simulation. It's good to know that I can repeat this: after removing the stabilization which I implemented later in these notes, I reproduced this same later fail time, which is still the gold standard (assuming that that is a legitimate measure).

These simulation failures are fundamentally different from the problems I was having before. Before the solutions were all smooth and pretty but the solves were either diverging or just going too slowly. I believe that this was simply because my relative tolerance was too high. If making the tolerance higher doesn't introduce oscillations in the solution and it allows the time steps to be much longer, than I am absolutely going to increase the tolerance! My priority is getting a solution. I don't care if there's too much diffusion for example. I don't give a damn about that.

However, now I'm seeing oscillations. Fundamentally different problem. A better problem in my opinion.

Three individual things I'm going to try (only change one variable at a time! Be patient!!!)

Test suite #1

1) Scale up the potential residual. Fail time: 5.50e-7 seconds
   Worse!
2) Increase the mesh density (Increased from 4000 to 8000 elements). Petsc failed with message: "Computed Nan differencing parameter h." Fail tiem again at 1.27e-8 seconds. So all the dependent variable solutions look beautifully smooth. However, there is an oscillation in the electron temperature at the cathode. This could be relevant because the electron temperature is in the source term for the electrons and ions. There are two things that should be done here. 1) I am going to try doubling the mesh resolution again to see if I can completely remove the oscillations. 2) It's worth looking into whether there are suggestions in combating the "Nan differencing parameter" for Petsc. Perhaps the recent messages on the Moose message board about Petsc error on fine mesh?
a) From the fine mesh thread, which didn't make any reference to my particular error, a guy went from using ds for -mat_mffd_type to wp and his problem worked. -mat_mffd_type is for matrix free formulations (thus the mf). There is another petsc option, -mat_fd_type, which I believe is used if you're actually forming the Jacobian matrix (perhaps Newton?).
b) In one thread from MOOSE mentioned the "Nan differencing parameter", John Mangeri said he was getting this error when using hypre and boomeramg. He said that when he went away from that preconditioner, he no longer experienced that problem. In a separate thread, Derek suggested that this generally means that a Nan is occurring somewhere in the residual, e.g. a divide by zero or sqrt of negative number. I personally have a hard time believing this though, especially considering my logarithmic formulation.

3) Decrease the relative tolerance. (Decreased from 1e-2 to 1e-3). Petsc failed with message: "Computed Nan differencing parameter h." Fail time: 1.27e-8 seconds. I didn't check how the solutions looked when the problem failed. If there were no oscillations in the solution, then the problem wasn't with the solution.
Way worse!

Test suite #2

1) Increased from 8000 to 16000 elements to see whether I can completely eliminate oscillations or see the "Nan differencing" error again. Fail time = 2.7911e-6. Reason for fail: looks like oscillations. It's always hard to know where the oscillations originate from. Is it oscillations in the bulk? Is it the oscillations at the boundary? Is it the electron temperature? Is it the electron density? Is it the argon density? The residuals of em, Arp, and mean_en are all about the same value. em is the highest.

Oscillations are my problem. Paths to solution:

1) Stabilization -> Artificial diffusion
   Back to 4000 elements. With the stabilization, fail at 2.59839e-6 seconds. It looks like a kink forms around 5.5e-7 seconds. There aren't nearly as many oscillations in this simulation as in the 16000 element simulation.
2) Mesh refinement -> Salome tutorial

Table of fail times for energy formulation:
4000 elements. nl_rel_tol = 1e-2. No stabilization. Fail time = 3.0207e-6 (normal fail)
4000 elements. nl_rel_tol = 1e-2. With stabilization. Fail time = 2.59839e-6 (normal fail)
4000 elements. nl_rel_tol = 1e-3. With stabilization. Fail time = 9.9839e-07 (normal fail)
16000 elements. nl_rel_tol = 1e-2. No stabilization. Fail time = 2.7911e-6 (normal fail)
8000 elements. nl_rel_tol = 1e-2. No stabilization. Fail time = 1.27e-8 (Nan differencing parameter)
4000 elements. nl_rel_tol = 1e-3. No stabilization. Fail time = 1.27e-8 (Nan differencing parameter)
8000 elements. nl_rel_tol = 1e-4. With stabilization. Fail time = 9.278e-08 (normal fail)
4000 elements. nl_rel_tol = 1e-2. No stabilization. Townsend form. Fail time = 6.493e-8 seconds. At least there no kinks this time. The problems were all in the oscillations at the cathode.
4000 elements. nl_rel_tol = 1e-2. With stabilization. Townsend form. Fail time = 2.7639e-7 seconds. At least there no kinks this time. The problems were all in the oscillations at the cathode. And yes there are a lot of oscillations. Hahahahahahahah. Not really that funny. Guess what. I'm going to figure out how to mesh the hell out of that cathode region with salome. That's the next task! Hey at least stabilization helped that time!!!

Note that the current Jacobian for the electrons for the Townsend formulation is absolutely wrong because I don't really know how to take the derivative of the absolute value of the flux. If the simulation fails because of this (pending that I can determine that that was indeed the cause of the failure), then I will invest some time into figuring out how to do that derivative correctly.

Running some simple tests, it looks like my scheme does indeed help stabilize.

I went back to a local field formulation. And set a new world record for simulation time with (almost) all the physics included: Fail time = 1.92576e-5 seconds. Woooo! Hahaha. Going to add some stabilization for the ions because it looks like I'm getting some kinks and maybe some stabilization would help. I don't see any kinks in the electrons.

Table of fail times for LFA formulation:
4000 elements. nl_rel_tol = 1e-2. Electron stabilization. Fail time = 1.92576e-5 (normal fail)
4000 elements. nl_rel_tol = 1e-2. Electron and ion stabilization. Fail time = 8.15811e-5 seconds (normal fail). New world record. Lots of oscillations in the electrons and ions, particularly the ions. This is likely without having a functional source term. But do those oscillations at really low concentration numbers matter? Because MOOSE is trying to minimize the residual, and in all places we are using the exponential of the concentration, or we are multiplying times the exponential of the concentration.

I've added secondary electrons in and stuff is looking pretty damn good if I might say so myself. The problem is I got this error message after 3.8 microseconds (using ds for -mat_mffd_type):

[0]PETSC ERROR: Petsc has generated inconsistent data
[0]PETSC ERROR: Differencing parameter is not a number sum = nan dot = -nan norm = -nan
[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.

This error message occurred despite the fact that step sizes were still large, the previous time step's nonlinear solve had converged, and all the variable solutions were smooth. So if I can fix whathever problem caused the error message, I think the solver could keep going for quite a while.

With -mat_mffd_type = wp, the simulation ran until 4.0 microseconds (so a little better). Solutions still looking smooth. This is the error message:

[0]PETSC ERROR: Petsc has generated inconsistent data
[0]PETSC ERROR: Computed Nan differencing parameter h
[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.

So slightly different error message.

After scaling the potential by 1e13, the solve yielded the NaN differencing parameter at 3.071e-7 seconds (solve was proceeding well until that point.)

I changed the nl_rel_tol to 1e-3 (potential scaling still 1e13). With this setting, I incur a "normal" fail, e.g. no petsc errors. Solve failed at 3.3636e-6 seconds. Why did the solve fail? It's not because of the linear solve. In one of the later diverged non-linear solves, the linear residual is dropping just fine. However, the non-linear residual drops a couple orders of magnitude immediately but then refuses to drop any more. I don't know what the best strategy to approach this is.

After scaling both the electrons and the ions by 1e2, the solve failed at 3.6097e-6 seconds, so a little later. Now when the solve is failing at the end of the simulation, the linear residuals are not dropping. Maybe this is a suggestion that I need to fix my preconditioner; well not fix, but add more terms and examine the petsc options for the preconditioner.

I implemented all the Jacobian elements that I could think of and then switched to a Newton solve. The solve failed at 3.2548e-6 seconds. A normal fail with most of the crashes occurring because of diverged line search. Observe some oscillations where the electron temperature drops to zero in the mean_en plot. However, the mean_en is not the main contributor to the residual: the Arp kernel is. However, obviously em and Arp are both functions of mean_en and the electron temperature. However, the electron temperature actually looks smooth. The solve essentially fails as soon as the electron temperature touches zero. Is that a coincidence? I think not.

Removing the ionization source term makes everything beautiful. The newton solve works beautifully...I get these nice giant implicit time steps. So it's the ionization source term that is gumming things up.

At the point where the mean_en goes to its small value, the value of mean_en - em is approximately -20...this will result in exp(-Eiz/(2/3*exp(mean_en - em))) having a value of 0, e.g. it's too small foating point value to be stored (I think).

After going to 400 elements, the solve failed at 1.2209e-6 seconds.

With the perfect preconditioner (finite difference preconditioner), the solve still failed. It failed at 5.8375e-6 seconds, so that is a world record. And what's interesting is that the electron temperature didn't touch zero as quickly. Two things left to try: 1) scale the potential residual. 2) Tighten the non-linear tolerance.

Scaled the potential residual: very interestingly the electron temperature touched down at zero earlier and the solve failed at 3.9943e-6 seconds as opposed to the much later time of 5.8375e-6 seconds without the potential scaling.

Last attempt: tighten the non-linear tolerance. After tightening the tolerances, the solve failed at 3.6947e-6 seconds. Again occurred immediately after touchdown of the electron temperature. I think it's safe to say that the problem is with my physics or with my residual implementation of the physics. I'm using a perfect preconditioner and a Newton solve, with smart choices for my relative and absolute tolerances. The problem is with the physics.

I figured out why the source term wasn't having an effect on the solution: the pre-exponential factor simply wasn't large enough. If I increased it by two order of magnitude than the effects of the source were immediately apparent. Before making this change, I was using an alpha of 0.35 for air. How does that compare with number from bolos? It's incredibly small is what it is. For example, for an electric field of 1.5e7, alpha predicted by both bolos and bolsig is around 2e5. This is absurd compared to 0.35. 6 order of magnitude difference.

Yea, I solved this problem. Boom. Zapdos. Boom.

There's a lot of things to improve...for example Comsol took 10 minutes and 47 seconds to solve this problem on two processor. I believe that it took Zapdos 6 hours to solve on one processor. Sure, I could theoretically run Zapdos on many more processors, but then I need to fix my Jacobian. Comsol used a direct solver called MUMPS. I did not use a direct solver as far as I know. I would actually have to look at the default settings for Newton on one processor with a full preconditioning matrix formed from finite differencing of the residuals (I believe it's finite differencing of the residuals?).

When I try to run with pc_type = lu, I get an error around the simulation time of 2.7e-5 seconds. The error is a petsc error and it says:
[0]PETSC ERROR: Zero pivot in LU factorization: http://www.mcs.anl.gov/petsc/documentation/faq.html#ZeroPivot
[0]PETSC ERROR: Zero pivot row 1 value 6.9212e-27 tolerance 2.22045e-14

I got this error twice. First time was with mat_fd_type = wp. Second time was with mat_fd_type = ds and with sub_pc_factor_shift = NONZERO.

For p = 0, Discontinuous Galerkin is identical to first-order finite volume. That's really cool!

Trying the RF discharge. With default line search (I don't know what the default is, the moose list would lead me to believe that it's bt), the sim starts diverging at 3.63113e-8 seconds. With line_search = basic, the sim can't even take a single time step. Using line_search = cp, the sim starts diverging at t = 7.47679e-8 seconds. So I suppose that that's better.

Reading the petsc manual, it looks like lu is the only direct solver available. To use it, use the options: -ksp_type=preonly; -pc_type=lu

From petsc manual: "Sometimes one is required to solve linear systems that are singular. That is systems with the matrix has a null space. For example, the discretization of the Laplacian operator with Neumann boundary conditions as a null space of the constant functions."

If trying to solve singular systems with direct solvers (or an incomplete factorization), it may still detect a zero pivot. Can run with additional options: -pc_factor_shift_type NONZERO -pc_factor_shift_amount <dampingfactor> to prevent the zero pivot. A good choice for the dampingfactor is 1.e-
10.

Files I need to fix after constructor update: RFIon, DCIon, RFElectron, DCElectron (both C and h files)

So for the RF plasma, keep in mind that I only calculated transport parameters for an electric field ranging from 1e3 to 1e6 V/m. There are problems with convergence for alpha with high electric fields (for this low pressure case of 1 torr). So after simulating, check and see whether the electric field stayed within those bounds.

Here is Comsol's electron and electron energy stabilization terms:

N_A_const*exp(-Ne*ccp.zeta)*test(Ne)
N_A_const*exp(-En*ccp.zeta)*test(En)

Simple and elegant. Only becomes large when Ne or En is small. Remember that Ne and En are the log of the actual density and total electron energy. This can be easily implemented when I return to an energy equation formulation. (Hell it could be implemented right now if I wanted to). Zeta is a parameter set by the user. Please start using user parameters again.

Having a lot of problems with diverged line search when the solution approaches steady state.

We want to solve F<mat>(x<vec>) = 0 where F is a nonlinear system of equations. (n equations, n unknowns in x<vec>)

Newton's method: x_k+1<vec> = x_k<vec> - J_inverse<mat>(x_k<vec>) * F<mat>(x_k<vec>)

J = F_prime, e.g. the Jacobian is the derivative of F with respect to the system variables (x)

In order to use this method, the Jacobian must be invertible, e.g. it cannot be singular. Singular matrices cannot be inverted.

Analytic expressions for the Jacobian can be provided by the user, or they can be calculated somewhat automatically (and approximately) by finite differencing F (the residuals).

For a Newton solve in which the Jacobian elements are provided by the user, the -mat_fd_type makes no sense. The only time -mat_fd_type makes sense is when the Jacobian is calculated with finite differencing. One can also use a matrix free formulation, and then the option is -mat_mffd_type.

To my mind the Jacobian and preconditioning matrices are two different matrices. But I'm not sure.

Pre-conditioning is done in iterative linear solves. LU is a direct solver; it's a full inversion of the Jacobian matrix. Thus, if we're using LU, there is no preconditioning.

From Wikipedia for solving the system:

A*x = b, we instead solve A*P_inverse*y = b, where y = P*x (right preconditioning) or P_inverse*(A*x-b) = 0 (left preconditioning, apparently more common).

"Typically there is a trade-off in the choice of P. Since the operator P^{-1} must be applied at each step of the iterative linear solver, it should have a small cost (computing time) of applying the P^{-1} operation. The cheapest preconditioner would therefore be P=I since then P^{-1}=I. Clearly, this results in the original linear system and the preconditioner does nothing. At the other extreme, the choice P=A gives P^{-1}A = AP^{-1} = I, which has optimal condition number of 1, requiring a single iteration for convergence; however in this case P^{-1}=A^{-1}, and applying the preconditioner is as difficult as solving the original system. One therefore chooses P as somewhere between these two extremes, in an attempt to achieve a minimal number of linear iterations while keeping the operator P^{-1} as simple as possible. Some examples of typical preconditioning approaches are detailed below."

Thus to give the minimal number of linear iterations, P should equal A. (However, minimal number of linear iterations probably is not your only criteria for an efficient, accurate simulation). So A, wikipedia's notation, is actually the Jacobian in my problems. Thus the perfect preconditioner (in terms of least numer of linear iterations) for the linear solve, would actually be the Jacobian itself!; using the Jacobian itself would result in only one linear iteration.

-snes_mf_operator: Activates default matrix-free Jacobian-vector products, and a user-provided preconditioning matrix as set by SNESSetJacobian()

MOOSE delineates two things: the preconditioning matrix, and the preconditioning process.

Oh, very interesting. Looking at the doco for SNESSetJacobian, there are two input matrices: A, the Jacobian matrix, and B, the preconditioner matrix, which "is usually the same as the Jacobian."

Thus Jed Brown's suggestion of using: "-snes_mf_operator -pc_type lu to see if the Jacobian you are
using is wrong." still doesn't make total sense to me. Ok, I think all of it is mostly making sense now. -snes_mf_operator activates default matrix-free Jacobian-vector products, and a user-provided preconditioning matrix as set by SNESSetJacobian(). As soon as we move to a matrix-free formulation, I believe we have to be using iterative linear solves. And that means using a preconditioner. How do we create the preconditioning matrix? By default it just uses the diagonals, e.g. whatever diagonal jacobian elements we've coded. Using SMP, full=true, we can say to build the entire preconditioning matrix, setting it equal (hopefully) to the Jacobian. However, we can also create the preconditioning matrix using finite differencing of the residual statements in the hope of creating a good approximation of the jacobian.

My hypothesis if I turn off the FDP preconditioner and use Jed's -snes_mf_operator -pc_type lu is that my preconditioning matrix will be formed from my diagonal jacobians (which are not well coded) and then I will be doing a "direct" matrix-free solve using lu. So I think that's what I saw; I saw some linear solves that were quite terrible. I hypothesize that that's because my preconditioning matrix was terrible. However, something that was quite interesting, and that makes me kind of happy, is that with Jed's options, it only took one non-linear iteration to solve the first time step. I'm not sure that I've seen that before. The solve failed at 2.94834e-6 seconds. Ok, now I'm going to try with Jed's options but with my preconditiong matrix set equal to the finite differenced jacobian. I'm hoping that this eliminates any bad linear solves. That's my hypothesis at least.

When doing direct matrix solves, the Jacobian affects the nonlinear solution strategy. However, in matrix free solves, the Jacobian is not actually present. It's only loosely present in that you use some emulation of the Jacobian in constructing your preconditioner that's used in the linear iterative solve. So in a Newton matrix solve, the effect of the Jacobian will be pronounced in the non-linear steps. In a New matrix free solve, the effect of the Jacobian will be pronounced in the linear steps. So I could play around with -mat_mffd_type

I believe that in any matrix free formulation is going to involve finite differences and a differencing parameter. E.g. the action of the Jacobian on a vector is approximated using finite differences! So I believe I actually have two finite differencing operations coming into play: 1) I am constructing my Jacobian/Preconditionig matrix using finite differencing of my residuals. 2) Then I am approximating the action of my Jacoban/Preconditioning matrix on a vector using finite differences. #2 comes into play when I decide to go to a matrix free formulation. It's possible to have no finite differencing in your problem: 1) construct your jacobian using analytical expressions 2) don't use a matrix free formulation. Simulation using the FDP preconditioner seems quite stable. It would probably go all the way to completion if I let it run. First divergence occurred at 1.14239e-5 seconds. Going to go back to direct LU and see when first divergence occurs. With direct LU, first divergence occurs at 1.37097e-5 seconds, so a little better.

The default KSP type is GMRES with a restart of 30, using modified Gram-Schmidt orthogonalization.

With default differencing parameter (I believe wp), here's Jacobian results:

Norm of matrix ratio 4.08771e-07 difference 0.0207982 (user-defined state)
Norm of matrix ratio 6.25133e-09 difference 8.56263e+09 (constant state -1.0)
Norm of matrix ratio 5.08948e-09 difference 9.62509e+08 (constant state 1.0)

With mat_fd_type=ds, Jacobian results:

Norm of matrix ratio 0.00932033 difference 474.28 (user-defined state)
Norm of matrix ratio 6.4317e-29 difference 8.80968e-11 (constant state -1.0)
Norm of matrix ratio 3.21208e-27 difference 6.0746e-10 (constant state 1.0)

With mat_fd_type=ds, mat_fd_coloring_err=1e-6, Jacobian results:

Norm of matrix ratio 7.87438e-05 difference 4.00641 (user-defined state)
Norm of matrix ratio 4.66307e-07 difference 6.38715e+11 (constant state -1.0)
Norm of matrix ratio 5.01161e-07 difference 9.47782e+10 (constant state 1.0)

I'm curious what the default mat_fd_coloring_err is, or if there is a default So if low numbers are generally better, maybe I should try ds? The two constant states are amazing but the user-defined state doesn't seem to be as good as wp. Again, I don't know what that means. -snes_type = test, can only be run if solve_type is set to NEWTON. It cannot be used JFNK. (matrix free methods).

With ds, the solve first diverged at 1.37097e-5 seconds, so it appears to have been identical to wp. With mat_fd_coloring_err=1e-6, the solve crashed almost immediately. So that pretty well reflected the Jacobian results!

In the IntTD.C materials file, I took away eta, so all input files need to be free of an eta column.

So I've coded the Jacobians for the kernels; I forgot about the bloody boundary conditions!!!

A huge damn victory for Moose and Zapdos tonight!!!!!! Using my analytical jacobian, a Newton solve, a direct LU method, I solved the problem to steady state in 168 seconds :-) Take that comsol!!!!

Testing Jacobian after adding back in the Electron Energy equation. Hand-coded Jacobian:

Norm of matrix ratio 0.00756585 difference 384.97 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

FDP Jacobian:

Norm of matrix ratio 3.78445e-06 difference 0.192568 (user-defined state)
Norm of matrix ratio 1.21061e-07 difference 1.65821e+11 (constant state -1.0)
Norm of matrix ratio 1.22993e-07 difference 2.326e+10 (constant state 1.0)

Without BCs:

Norm of matrix ratio 0.00756585 difference 384.97 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

Without electron energy time derivative:

Norm of matrix ratio 0.00756608 difference 384.97 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

Without electron energy kernel:

Norm of matrix ratio 3.96561e-06 difference 0.201765 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

That last test kind of reveals where the error is, doesn't it?

After maybe modifying the source terms in the ElectronEnergyKernel:

Norm of matrix ratio 0.00756585 difference 384.97 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

After removing Joule heating:

Norm of matrix ratio 0.00756585 difference 384.97 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

After removing all but advection and diffusion:

Norm of matrix ratio 9.24413e-06 difference 0.470329 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

Same but using FDP:

Norm of matrix ratio 3.78303e-06 difference 0.192476 (user-defined state)
Norm of matrix ratio 1.21061e-07 difference 1.65821e+11 (constant state -1.0)
Norm of matrix ratio 1.22993e-07 difference 2.326e+10 (constant state 1.0)

With all physics in ElectronEnergyKernel removed:

Norm of matrix ratio 3.96561e-06 difference 0.201765 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

Two things to do:
1) Make a couple of results summary slides of what I've done so far in Zapdos
2) 5 modules (kernels, BCs) I would want to add to make Zapdos/Moose relevant for broader LTP community

The steady state electron density as determined by zapdos does not differ based on whether an iterative or direct solve is used. The peak electron density for zapdos (both direct and iterative solutions) and the peak electron density for comsol differ by 5.1\% (Which is actually fairly significant I suppose).

MKL PARDISO: Intel Math Kernel Library PARallel Direct SOlver

PARDISO in Intel MKL was originally developed by the Department of Computer Scienc at the University of Basel. It can be obtained at http://www.pardiso-project.org

After adding back in advection and diffusion and on-diagonal jacobians. I've reviewed these residuals and jacobians and they are exactly correct.

Norm of matrix ratio 7.10444e-06 difference 0.361464 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

After adding back advection off-diag jacobian. After review, this looks perfect.

Norm of matrix ratio 9.24413e-06 difference 0.470329 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

After adding in Joule-heating and associated potential off-diag jacobian:

Norm of matrix ratio 9.24413e-06 difference 0.470329 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

After adding in Joule-heating and associated em off-diag jacobian:

Norm of matrix ratio 9.24413e-06 difference 0.470329 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

After adding back reaction source and on-diag jacobian:

Norm of matrix ratio 0.00756587 difference 384.97 (user-defined state)
Norm of matrix ratio 1.57362e-07 difference 2.15543e+11 (constant state -1.0)
Norm of matrix ratio 1.56667e-07 difference 2.96284e+10 (constant state 1.0)

When trying to solve using super lu, I keep getting diverged_line_search after 0 iterations when I use the default line search. If I set line_search='none', then the initial divergence occurs at (at least approximately) the same time step, but the divergence message is: diverged_fnorm_nan.

When solving in parallel with moose's default parallel iterative options, the solve fails at 1.69359e-7 seconds. The final divergence reason is given by:

  Linear solve did not converge due to DIVERGED_NANORINF iterations 0
Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 2

Looking at all the divergences, the DIVERGED_LINEAR_SOLVE was about as frequent as DIVERGED_LINE_SEARCH.

Just confirming that in this test problem the Jacobian is good, when I run the fast_solve_parallel_iterative_attempt in test mode, with a mesh size of 1, here's what the results look like:

Norm of matrix ratio 4.44187e-08 difference 0.00226001 (user-defined state)
Norm of matrix ratio 2.33137e-08 difference 3.19335e+10 (constant state -1.0)
Norm of matrix ratio 2.20073e-08 difference 4.16197e+09 (constant state 1.0)

So excellent. When I change the mesh size to 100, here are the results (still good):

Norm of matrix ratio 1.58941e-07 difference 0.00808669 (user-defined state)
Norm of matrix ratio 3.48091e-08 difference 4.7679e+10 (constant state -1.0)
Norm of matrix ratio 3.26732e-08 difference 6.17907e+09 (constant state 1.0)

What am I trying to do right now? I'm trying to get my gold_serial_case to solve in parallel mode, using 1) iterative methods and 2) direct methods. Another thing that I've been working on is adding back an electron energy equation; I'm currently debugging the Jacobian for that.

I'm looking at the iterative solve right now. If I run with one mesh element, I get a floating point exception which appears to occur at line 52 in ElectronKernel.C. However, if I run with 100 mesh elements, than I no longer get any FPEs/program breaks when running ./zapdos-dbg (without gdb --args). So the FPE does indeed come from an exponential overflow. The last value of _u[_qp] before the floating point exception is signalled is: 7.88655e+07. Since I know that exp(1000) produces inf, then obviously 7.89e7 will also produce an inf. This exercise kind of just goes to show that you need an appropriate number of mesh elements or else your solution will have major problems.

I ran the direct solve with one mesh element (./zapdos-dbg) and got no FPEs.

Don't forget how to print values in the debugger:

p <variable>

or

print <variable>

Ok, using default MOOSE iterative with mpirun and an abs tolerance of 1e-6 (e.g. all the parameters the same from the fast serial lu gold solve except for the solve type), the solve fails at 2.49633e-07 seconds. All the divergence reasons are: DIVERGED_LINE_SEARCH. First divergence occurs at: 4.84966e-08 seconds. No line_search option specified for the settings in this passage (so default). Now the iterative solve has major problems with its solution variables, and this could easily explain why this solve fails. There are clear spikes at L/4, L/2, and 3L/4 in the dependent variable Arp and in the auxiliary variable EField. So again, parallel iterative solves are failing.

Using direct super LU with line_search='none', first divergence is at 1.28117e-07 seconds, and the message is DIVERGED_FNORM_NAN. All the divergence messages are DIVERGED_FNORM_NAN. Final failure occurs at 2.40013e-07 seconds. (So fairly comparable to the iterative attempt). With default line_serach, first divergence is identically at 1.28117e-07 seconds. All divergence reasons are DIVERGED_LINE_SEARCH. Final failure identically occurs at 2.40013e-07 seconds. So my initial intuitive reaction was correct: despite the change in line_search, the solve patterns are the exact same. Looking at the default line search results in paraview, every solution variable looks just fine at the premature end of the simulation: all variable profiles are smooth and continuous, so any reason for the solve fail can't really be found there.

The fun part now in both direct and iterative cases is to figure out: what's causing the diverged line search???

Alright, doing some cool stuff. Did a serial iterative solve. | Zapdos Performance: Alive time=150.776, Active time=149.795. pc_type=ilu, ksp_type=gmres, snes_type=newtonls. Results look perfect.

With pc_type=gamg, the problem wasn't even able to take a single time step.

hypre is meant for massively parallel computing (taken from its homepage).

With pc_type=asm, sub_pc_type=ilu, the problem solves. Solution summary: | Zapdos Performance: Alive time=159.017, Active time=158.009.

With:

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
  petsc_options_value = 'asm lu NONZERO 1.e-10'

Solution summary: | Zapdos Performance: Alive time=153.569, Active time=152.56

With pc_type=bjacobi, the problem solves. Solution summary: | Zapdos Performance: Alive time=149.462, Active time=148.471 (Record)

Next step: try suggestions from petsc email for direct solves.

This should be my compiler environment the next time I update_and_rebuild libmesh and recompile MOOSE and Zapdos:
export CC="ccache mpicc"
export CXX="ccache mpicxx"
export CCACHE_SLOPPINESS=time_macros

Oh heck yea: solved my problem in parallel!!!!!!!!! Wooooooooooooooooooooooo!!!!!!!!!!!!!!!!! Got mumps working correctly and it solved like a bloody damn champion. To summarize here were the petsc options:

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -ksp_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu NONZERO 1.e-10 preonly mumps'

And here is the execution summary: | Zapdos Performance: Alive time=133.099, Active time=131.212

So having the parallel solve didn't really speed up the solution, but it's a demonstration of principle gosh darn it!!!

Ok, this is really freaking weird. I went back to an old commit where there was no ZapdosRevision.h file, but presumably there should have been no new source code changes to compile, and when I typed make -j4, here's the first output:

MOOSE Generating Header /home/lindsayad/projects/zapdos/include/ZapdosRevision.h...

Then it proceeds to compile zapdos source files (presumably it was going to compile all of them before I manually stopped it.)

However, when I reset my pointer to point at the newest commit (where ZapdosRevision.h is included), it still recompiled all my source code, only this time it didn't have that first line:

MOOSE Generating Header /home/lindsayad/projects/zapdos/include/ZapdosRevision.h...

It's important to note that after recompilation, my gold mumps script executed without complaint, so that helps to assuage some of my fears that recompilation would prevent me from seeing petsc packages.

Using parallel Bjacobi with lu as the sub_pc_type, the solve failed at 2.35183e-07 seconds. And sure enough, I see the familiar kinks in the solution in the places where the processors come together. This was done on 4 processors.

However, when I went to 2 processors the problem solved! Performance summary: | Zapdos Performance: Alive time=111.615, Active time=110.712 (Record!)

In early time steps, I could still see some kinks in the Arp solution and occasionally even in the em solution. However, at longer time steps the kinks disappeared. I would hypothesize that these kinks only appear when the potential is almost perfectly linear and the electric field is incredibly uniform such that any tiny miscommunication between the processors will induce a very noted change in the electrostatics on an automatic scale. And these kinks in the electrostatics directly impact the particle drift motion. However, on longer time scales when the particles begin to affect the electric field, these affects are much stronger and on a much larger scale than the processor miscommunication so that the processor miscommunication no longer becomes a problem.

With that same run as the one immediately above on the desktop, the problem solved in 97 seconds with console output to the screen. When I tried running with 4 processors, I got the same crash that I've come to know. It also occurred at 2.35183e-07 seconds. So the success of running on 2 cores comes not from having the correct number of physical cores but from the formulation of the problem.

Reducing the amount of output to the log file doesn't seem to have helped the problem solve any faster. Scaling comparison:

2 processors: 97.0632 seconds
1 processor: 141.503 seconds

So computation time was reduced by 31\%. E.g. it took 69\% as long with 2 processors as it did with 1. With perfect parallelization that 69 number would be 50.

Here's the Jacobian test output for this ionization source residual and corresponding Jacobian:

    -_test[_i][_qp]*_rate_coeff_ion_en[_qp]*_Ar[_qp]*-_Eiz_en[_qp]*_em[_qp]

Norm of matrix ratio 7.86836e-08 difference 0.00400463 (user-defined state)
Norm of matrix ratio 1.76423e-08 difference 2.41652e+10 (constant state -1.0)
Norm of matrix ratio 1.9803e-08 difference 3.74509e+09 (constant state 1.0)

Here's the Jacobian test output for this ionization source residual and corresponding Jacobian:

    -_test[_i][_qp]*_rate_coeff_ion_en[_qp]*_Ar[_qp]*-_Eiz_en[_qp]*_em[_qp]*em[_qp]

Norm of matrix ratio 7.86836e-08 difference 0.00400463 (user-defined state)
Norm of matrix ratio 1.76423e-08 difference 2.41652e+10 (constant state -1.0)
Norm of matrix ratio 1.9803e-08 difference 3.74509e+09 (constant state 1.0)

Here's the Jacobian test output for this ionization source residual and corresponding Jacobian:

    -_test[_i][_qp]*_rate_coeff_ion_en[_qp]*_Ar[_qp]*-_Eiz_en[_qp]*std::sqrt(_em[_qp])

Norm of matrix ratio 7.86834e-08 difference 0.00400462 (user-defined state)
Norm of matrix ratio -nan difference -nan (constant state -1.0)
Norm of matrix ratio 1.9803e-08 difference 3.74509e+09 (constant state 1.0)

Here's the Jacobian test output for this ionization source residual and corresponding Jacobian:

    -_test[_i][_qp]*_rate_coeff_ion_en[_qp]*_Ar[_qp]*-_Eiz_en[_qp]*std::exp(_em[_qp])

Norm of matrix ratio 1.45574 difference 74090.6 (user-defined state)
Norm of matrix ratio 1.76423e-08 difference 2.41652e+10 (constant state -1.0)
Norm of matrix ratio 1.9803e-08 difference 3.74509e+09 (constant state 1.0)

Here's the Jacobian test output for this ionization source residual and corresponding Jacobian:

    -_test[_i][_qp]*_rate_coeff_ion_en[_qp]*_Ar[_qp]*-_Eiz_en[_qp]*std::exp(_u[_qp])

Norm of matrix ratio 4.08685e-08 difference 0.0123279 (user-defined state)
Norm of matrix ratio 1.76423e-08 difference 2.41652e+10 (constant state -1.0)
Norm of matrix ratio 1.9803e-08 difference 3.74509e+09 (constant state 1.0)

I think I'm figuring out what the constant states mean in the Jacobian test. I believe that constant state 1.0 is where all the dependent variables are set to 1.0 (an IC). User-defined state would be with my defined ICs. Yes my hypothesis is confirmed. So if I set all my variables to ICs of 1 or -1, my Jacobian looks really damn good!

It should also be noted that the variable scaling affects the "Norm of matrix ratios."

Trying with electron energy formulation again. Solve failed at 1.05637-07 seconds. All the variables look pretty nice and smooth, so I don't really know what the problem is. Divergence reasons were fairly evenly split between diverged_max_its and diverged_line_search. At the last convergence step, here were the variable residual norms:

                  potential: 4.042e-10
                  em:        2.6627e-08
                  Arp:       3.20928e-10
                  mean_en:   0.000475877

With an abs_tolerance=1e-2 and using mumps in parallel, the solve still failed at 1.05637e-07 seconds. With the ionization source terms elimination in the electron energy equation, the solve actually failed earlier at 7.59991e-08 seconds. This actually makes since because I still have physical production of electrons and ions in their respective kernels. Dumb me.

Tried again with electron energy formulation, but using PJFNK instead of NEWTON. Solve got a bit farther this time: failed at 3.39951e-07 seconds. Instead of any DIVERGED_LINE_SEARCH, all fails were due to DIVERGED_MAX_IT. i think diverged_line_search only really appears when doing Newton solves. The only reason it went farther is because I forgot to turn the ionization terms in the em and ion kernels back on again. With those turned back on, the solve fails at 1.08282e-07 seconds. Again, looking at the solution variables, I can see no reason why the solve should fail. However, if I look at the log, the linear residuals don't seem to be coming down nearly fast enough, suggesting to my continued bewilderment, that my Jacobian is somehow wrong. I'm going to try to solve using FDP.

Using FDP and a direct LU solve, the solve failes at 1.113e-07 seconds. Record, wooo! Haha...no. Why the hell does it fail? All the variables are gosh darn smooth! I don't see any oscillations anywhere!

With block-jacobi:

 9 Nonlinear |R| = 5.601308e-06
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 9
 Solve Converged!

With asm:

 9 Nonlinear |R| = 5.601308e-06
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 9
 Solve Converged!

So I have no bloody idea whether my command line petsc options are working. Actually it most certainly appears to not be working because when I pass -ksp_type=preonly on the command line, linear iterations still appear in the log file. Petsc options are still working from the input file, however. Another thing that appears to not be working is that when I pass the -help option on the command line, the values that appear in <> for the petsc variables don't represent the values being used in my program. Rather they appear to represent the defaults. This is strange because I feel like I have a memory of being able to see my values represented in <>.

I see the same behavior on my laptop. Perhaps it's a petsc version thing. Perhaps with version 3.5.x of petsc I could see the values of the variables I had changed instead of the defaults. It would be quite nice to confirm that I've changed the variables that I think I've changed.

With initial values of 1e17, the simulation fails at 8.81971e-08 seconds, so worse. With initial value of 1e7, the simulation fails at 1.26705e-07 seconds, so better. The reason for crashing has to be has to be boundary conditions. With my diffusion stabilization, the crash occurs at the same time: 1.26705e-07 seconds. I'm speculating that it's the electron energy boundary condition.

I changed the electron energy boundary condition in such a way that I would have assumed it would make a zero gradient boundary condition at both boundaries, but that didn't happen. In fact absolutely nothing appears to have changed. Solve still failed at 1.26705e-07 seconds. Changed to Neumann BCs for the electron energy. This did not change the fail time, but there were more oscillations at the right boundary...this doesn't surprise me since this BC change prevented energy from flowing out of the right boundary which it wants to be able to do. It's like a log jam with stuff piling up, leading to the oscillations. Remember, always try to make the problem as simple as possible!!!

Alright, so what I did is I left only the Joule heating and ionization energy loss terms in the Electron Energy Formulation, meaning that there's no energy loss from elastic collisions and there's no energy motion, negating the need for electron energy boundary conditions. The resulting simulation lasted a tiny bit longer: 1.29151e-07 seconds. Next flub: increase the ion mobility and diffusivity and see whether that removes the ion oscillations at the right boundary; i.e. also determine whether the ion oscillations might be responsbile for the simulation failures. I would think not, but can't know until we try.

Oh yes, in case I forgot to note this to myself: the fact that this problem won't solve even with the "theoretically perfect" Jacobian generated by FDP could mean that my hand coded Jacobians are not the reason for the failure of my hand coded runs.

With ion mobility and diffusivity increasd by 100 fold, the simulation lasted until 1.33596e-07 seconds. Improvement! Haha. The ion oscillations at the right boundary disappeared only to be replaced with electron oscillations. And interestingly the electric field magnitude increased in that anode region. Next thing I'm trying is decreasing the max voltage by a factor of 10, e.g. to 1e4. This seems to have done the trick. Currently the simulation is at 1e-4 seconds and still powering away. VERY IMPORTANT NOTE: the two unphysical things I think I currently have in this apparently successful simulation are: 1) the electron energy boundary condition is allowing energy to flow in (I think) 2) the ion mobility and diffusivity are a factor of 100 higher than what I've used in the past.

The problem did in fact solve. Boom goes the dyanamite. Dropping the mike :-)  Alive time: 26489.3 seconds, corresponding to 7.358 hours.

I tried going back to SMP formulation, and the solve failed at 8.28004e-07 seconds. Other changes I made other than just going back to SMP: I restored the ion transport coefficients to their default values, and I reduced the absolute tolerance to 1e-5.

Uhhhhh.....somehow my Jacobian is correct now? That would make no bloody sense. This is exceedingly mysterious to me but here are the results:

Norm of matrix ratio 4.58232e-08 difference 233.939 (user-defined state)
Norm of matrix ratio 1.32657e-08 difference 1.66194e+09 (constant state -1.0)
Norm of matrix ratio 1.54792e-08 difference 3.20123e+08 (constant state 1.0)

I don't think I changed a single thing to make this happen. Apparently reducing the potential_bc_func reduces the norm of the matrix-ratio of the user-defined state.

For potential = 1e3: Norm of matrix ratio 3.2554e-08 difference 23.3871 (user-defined state)
For pot = 1e4: Norm of matrix ratio 4.58232e-08 difference 233.939 (user-defined state)
For pot = 1e5: Norm of matrix ratio 4.64405e-08 difference 2359.25 (user-defined state)

With potential at 5e3, solve failed at 1.6731e-06 seconds, so better.
With potential at 1e3, the solve succeeded but only because the potential was below breakdown.

Maybe there could have been one thing different about my Jacobians: I had _a set to 1.0 for all boundaries in the Electron energy boundary condition source file. This is a change that should be in my simulation, but it's too far along right now for me to want to stop it.

1 mesh element, 4 processors: Norm of matrix ratio 2.05868e-08 difference 11.6923 (user-defined state)
4 mesh elements, 4 processors: Norm of matrix ratio 2.35982e-08 difference 13.4026 (user-defined state)
4 mesh elements, 1 processor, mpirun: Norm of matrix ratio 2.35982e-08 difference 13.4026 (user-defined state)
4 mesh elements, 1 processor, no mpirun: Norm of matrix ratio 2.35982e-08 difference 13.4026 (user-defined state)

Woooooooooooooooooooooooooooooooooooooooooooooooooooooo!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This truly represents the conquering of one of my bullet points in my MOOSE strategic plan. I got the electron energy formulation problem to solve in ~11 minutes :-) What was the trick? There's always a trick...the trick was to do as is oft repeated on the Moose discussion forum: don't allow outlier residual norms. Previously, my electron energy residual was 5 orders of magnitude larger than my em residual. After scaling it back by 5 orders of magnitude, I saw the orders of magnitude speed-up in my calculation. When I went back and solved in this manner, the active time was 864.56 seconds, corresponding to 14.4 minutes (this time I actually printed a performance log so 14.4 minutes should be taken as being the correct time compared to 11 minutes).

Perhaps the only thing that I can anticipate people complaining about is that I'm only applying 350 Volts over 1 cm: 3.5e2/1e-2 = 3.5e4 V/m. In Go's paper, he uses 210 V/1 mm = 2.1e2/1e-3 = 2.1e5 V/m. Reasons to explain this discrepancy: no possible attachment or recombination processes included in my model. Also, the transport/rate coefficient parameters were somewhat arbitrarily chosen.

So I see that my CoupledIntegratedBC is apparently having absolutely zero effect. Having the matched values BC causes v to match the value of u. u is 1 in it's whole domain (block 0). If the matched values BC is not active, then v is 0 everywhere in its domain.

From DGKernel.h: "The DGKernel class is responsible for calculating the residuals for physics on internal sides (edges/faces)."

Assembly::elemVolume(): Returns the reference to the current element volume
Assembly::elem(): Returns the current element
Assembly::side(): returns the current side
Assembly::neighborSide(): returns the current neighboring side
Assembly::sideElemVolume(): returns reference to the volume (or length [I believe that sideElemVolume should have units of Area if elemVolume has units of Volume, units of Length if elemVolume has units of Area, and be equal to one if elemVolume has units of Length. If that belief is true, then this calculation would have units of Length: h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.) which would make sense.]) of current side element
Assembly::neighborVolume(): returns reference to current neighbor volume

I think I'm starting to get an idea of what I need to do. I need to be able to pass a side (the only side in my interfacial case) that I want the DGKernel to act on. I need this kernel to have access both to u on the "inside" and to the coupled_var on the "outside." And in my particular case I want to know _grad_u at the "inside" side, and _grad_coupled_var at the "outside" side. Let's go explore IntegratedBC to see how one might access the values of those parameters (_u, _grad_u, _coupled_var, _grad_coupled_var) on sides (boundaries).

Here's constructor in KernelBase.C:

    MooseObject(parameters),
    BlockRestrictable(parameters),
    SetupInterface(parameters),
    CoupleableMooseVariableDependencyIntermediateInterface(parameters, false),
    FunctionInterface(parameters),
    UserObjectInterface(parameters),
    TransientInterface(parameters, "kernels"),
    PostprocessorInterface(parameters),
    MaterialPropertyInterface(parameters, blockIDs()),
    RandomInterface(parameters, *parameters.get<FEProblem *>("_fe_problem"), parameters.get<THREAD_ID>("_tid"), false),
    GeometricSearchInterface(parameters),
    Restartable(parameters, "Kernels"),
    ZeroInterface(parameters),
    MeshChangedInterface(parameters),
    _subproblem(*parameters.get<SubProblem *>("_subproblem")),
    _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
    _sys(*parameters.get<SystemBase *>("_sys")),
    _tid(parameters.get<THREAD_ID>("_tid")),
    _assembly(_subproblem.assembly(_tid)),
    _var(_sys.getVariable(_tid, parameters.get<NonlinearVariableName>("variable"))),
    _mesh(_subproblem.mesh()),
    _current_elem(_var.currentElem()),
    _current_elem_volume(_assembly.elemVolume()),
    _q_point(_assembly.qPoints()),
    _qrule(_assembly.qRule()),
    _JxW(_assembly.JxW()),
    _coord(_assembly.coordTransformation()),

    _test(_var.phi()),
    _grad_test(_var.gradPhi()),

    _phi(_assembly.phi()),
    _grad_phi(_assembly.gradPhi()),

Here's constructor in DGKernel.C:

    MooseObject(parameters),
    SetupInterface(parameters),
    TransientInterface(parameters, "dgkernels"),
    FunctionInterface(parameters),
    UserObjectInterface(parameters),
    NeighborCoupleableMooseVariableDependencyIntermediateInterface(parameters, false, false),
    TwoMaterialPropertyInterface(parameters),
    Restartable(parameters, "DGKernels"),
    MeshChangedInterface(parameters),
    _subproblem(*parameters.get<SubProblem *>("_subproblem")),
    _sys(*parameters.get<SystemBase *>("_sys")),
    _tid(parameters.get<THREAD_ID>("_tid")),
    _assembly(_subproblem.assembly(_tid)),
    _var(_sys.getVariable(_tid, parameters.get<NonlinearVariableName>("variable"))),
    _mesh(_subproblem.mesh()),

    _current_elem(_assembly.elem()),
    _current_elem_volume(_assembly.elemVolume()),

    _neighbor_elem(_assembly.neighbor()),
    _neighbor_elem_volume(_assembly.neighborVolume()),

    _current_side(_assembly.side()),
    _current_side_elem(_assembly.sideElem()),
    _current_side_volume(_assembly.sideElemVolume()),

    _coord_sys(_assembly.coordSystem()),
    _q_point(_assembly.qPointsFace()),
    _qrule(_assembly.qRuleFace()),
    _JxW(_assembly.JxWFace()),
    _coord(_assembly.coordTransformation()),

    _boundary_id(parameters.get<BoundaryID>("_boundary_id")),

    _u(_var.sln()),
    _grad_u(_var.gradSln()),

    _phi(_assembly.phiFace()),
    _grad_phi(_assembly.gradPhiFace()),

    _test(_var.phiFace()),
    _grad_test(_var.gradPhiFace()),

    _normals(_var.normals()),

    _phi_neighbor(_assembly.phiFaceNeighbor()),
    _grad_phi_neighbor(_assembly.gradPhiFaceNeighbor()),

    _test_neighbor(_var.phiFaceNeighbor()),
    _grad_test_neighbor(_var.gradPhiFaceNeighbor()),

    _u_neighbor(_var.slnNeighbor()),
    _grad_u_neighbor(_var.gradSlnNeighbor())

Here's constructor in IntegratedBC:

    BoundaryCondition(parameters),
    RandomInterface(parameters, _fe_problem, _tid, false),
    CoupleableMooseVariableDependencyIntermediateInterface(parameters, false),
    MaterialPropertyInterface(parameters),
    _current_elem(_assembly.elem()),
    _current_elem_volume(_assembly.elemVolume()),
    _current_side(_assembly.side()),
    _current_side_elem(_assembly.sideElem()),
    _current_side_volume(_assembly.sideElemVolume()),

    _normals(_var.normals()),

    _qrule(_assembly.qRuleFace()),
    _q_point(_assembly.qPointsFace()),
    _JxW(_assembly.JxWFace()),
    _coord(_assembly.coordTransformation()),

    _phi(_assembly.phiFace()),
    _grad_phi(_assembly.gradPhiFace()),

    _test(_var.phiFace()),
    _grad_test(_var.gradPhiFace()),

    _u(_is_implicit ? _var.sln() : _var.slnOld()),
    _grad_u(_is_implicit ? _var.gradSln() : _var.gradSlnOld()),

Conclusion: A buttload of similarities between DGKernel and IntegratedBC

I can do something exceedingly stupid. If I could set it up to determine whether a side is on the internal boundary, then the residual is non-zero. However, if the side is not on the boundary, then I set the residual to zero. Things I need to figure out:

1) How to tell whether a side is on a boundary
2) How to access a variable from another block

const std::vector<DGKernel *>& DGKernelWarehouse::active() ====> Get the list of all active kernels

void ComputeResidualThread::onInternalSide(const Elem * elem, unsigned int side) ===> elem: element we are on. Side: Local side number of the element 'elem'. Hmm local number may not be very helpful

ComputeResidualThread.C is where the magic happens. Here are three straight function definitions in the file:

ComputeResidualThread::onElement(const Elem *elem)
ComputeResidualThread::onBoundary(const Elem *elem, unsigned int side, BoundaryID bnd_id)
ComputeResidualThread::onInternalSide(const Elem *elem, unsigned int side)

onElement goes to the KernelWarehouse
onBoundary goes to the BCWarehouse and looks for active IntegratedBCs
onInternalSide goes to the DGKernelWarehouse

Here are the includes at the start of ComputeResidualThread:

#include "ComputeResidualThread.h"
#include "NonlinearSystem.h"
#include "Problem.h"
#include "FEProblem.h"
#include "KernelBase.h"
#include "IntegratedBC.h"
#include "DGKernel.h"
#include "Material.h"
// libmesh includes
#include "libmesh/threads.h"

Some member variables of ComputeResidualThread: _sys(sys), _kernel_type(type), _num_cached(0)

ComputeResidualThread also inherits from ThreadedElementLoop. ThreadedElementLoop inherits from ThreadedElementLoopBase.

THREAD_ID ThreadedElementLoopBase<RangeType>::_tid

To get active integrated boundary conditions, we call BCWarehouse::activeIntegrated(BoundaryID boundary_id, std::vector<IntegratedBC*>& active_integrated)

And here's the function content:
{
  active_integrated.clear();

  if (_bcs.find(boundary_id) != _bcs.end())
    for (std::vector<IntegratedBC *>::const_iterator it = _bcs.at(boundary_id).begin(); it != _bcs.at(boundary_id).end(); ++it)
      if ((*it)->isActive())
        active_integrated.push_back(*it);
}

std::map<BoundaryID, std::vector<IntegratedBC *> > BCWarehouse::_bcs

std::map has member functions end, find, begin, etc.
find: get iterator to element
begin: return iterator to beginning
end: return iterator to end
at: access element
All of the above are public member functions

 _sys.getBCWarehouse(_tid).activeIntegrated(bnd_id, bcs);
std::vector<DGKernel *> dgks = _sys.getDGKernelWarehouse(_tid).active();

The warehouses are important. In DGKernelWarehouse.C, the following member functions are defined: initialSetup, timestepSetup, residualSetup, jacobianSetup, addDGKernel(MooseSharedPointer<DGKernel> & dg_kernel), updateActiveDGKernels(Real /*t*/, Real /*dt*/)

const std::vector<DGKernel *> & active() const { return _active_dg_kernels; }

void
DGKernelWarehouse::updateActiveDGKernels(Real /*t*/, Real /*dt*/)
{
  _active_dg_kernels.clear();

  // add kernels that live everywhere
  for (std::vector<DGKernel *>::const_iterator it = _all_objects.begin(); it != _all_objects.end(); ++it)
  {
    DGKernel * dg_kernel = *it;
    // FIXME: add startTime/stopTime to DGKernel
//    if (dg_kernel->startTime() <= t + (1e-6 * dt) && dg_kernel->stopTime() >= t + (1e-6 * dt))
      _active_dg_kernels.push_back(dg_kernel);
  }
}

Following down the trail. Now pursuing _all_objects

_all_objects: All instances of objects (raw pointers)

onBoundary parameters:

elem : The element we are checking is on the boundary
side : The side of the element in question
bnd_id : ID of the boundary we are at

onBoundary gets called when doing boundary assembling

I might have been wrong about where the magic happens. The magic might happen in NonlinearSystem::computeResidual. I may return to this, but from Yaqi it sounds like Mortar could be what I need.

Sideset and nodeset naming are described in the MOOSE FAQ.

FaceFaceConstraint has available Params: "variable" ==> This is substituted with Lagrange multiplier; "master_variable" ==> This is the variable on the master side of the mortar interface, e.g. 'u'; "slave_variable" ==> This is the variable on the slave side of the mortar interface, e.g. 'v', however, if the "slave_variable" is not specified, then it is set to "master_variable".  Very cool

 NeighborCoupleable(parameters, nodal, neighbor_nodal)

nodal and neighbor_nodal are both of type bool.

NodeFaceConstraint actually has _grad in it

EqualValueConstraint inherits from FaceFaceConstraint
TiedValueConstraint inherits from NodeFaceConstraint
CoupledTiedValueConstraint also inherits from NodeFaceConstraint

getValue comes from the MooseVariable class

There are three fundamental constraints in the moose framework:

NodalConstraint
NodeFaceConstraint
FaceFaceConstraint

FaceFaceConstraint has to be the most natural fit for what I want to do right??????

NodeFaceConstraint has a grad definition. Nodal and FaceFace do not. The derivatives of NodeFaceConstraint are TiedValueConstraint and CoupledTiedValueConstraint. Neither of these have grad in them. None of the module source files that inherit from NodeFaceConstraint have grad in them either.

I think that I need to think about C++, arrays, and initialization. After declaration, access the values of the arrya using:

name[index]

MooseVariable::gradSln() returns _grad_u, e.g. it returns a VariableGradient

What does it mean when I do _grad_u_master[_qp] = _master_var.gradSln() ??

That just doesn't seem right to me.

It appears that my MeshModifier sidesets and nodesets don't get created until after my application is executed. E.g. they don't appear in peacock_run_tmp_mesh.e

Yes it appears that in the mesh file for conforming_two_var.i, there is a pre-existing subdomain called "Unnamed block ID: 1000 Type: EDGE2"

Recall that these mesh files are in the exodus format

After executing the conforming_two_var.i file, the mortar subdomain simply becomes "middle" which is its short name.

It's become apparent that SubdomainBoundingBox only works if you draw a box, e.g. you cannot try and pass a line and try and build a subdomain with it.

It looks like I either need to use an external package or write some new MooseModifier code.

I need type EDGE2

By default, when salome exports a mesh, it doesn't export any side or node sets. What should happen is that the side set and node set ids should be the same right? And there should be a set for every boundary/interface.

In the *.unv file, there are coordinates for each node. And perhaps the node ID is given as the first piece of information for each node.

2411 denotes nodes
2412 denotes elements
2467 denotes groups (of elements or nodes)

In the elements description, the description of 2D face elements differs from the description of 1D edge elements.

In my 3x3 mesh case example, there are 16 nodes. There are 9 face elements and 12 edge elements for a total of 21 elements. There appear to be 11 numbers that describe the face elements and 10 numbers that describe the edge elements.

A decent thread on the UNV file format: http://www.code-aster.org/forum2/viewtopic.php?id=14355

So much to learn here. Not quite sure why I seem to be struggling so much. Programmatic learning:

I'm looking at between.i. It loads a file mesh with element blocks 'left' and 'right.' Without conducting any mesh modification, no sidesets or nodesets appear when I load the _in.e file in paraview. However, when SideSetsBetweenSubdomains is called as a MeshModifier, a single sideset called in_between is created as well as a Nodeset called "Unnamed set ID: 1". So my question is: why the nodeset creation??

Full-dimensional elements have sides. But not all sides are what we consider "boundaries." Looking at the code in SideSetsBetweenSubdomains, it looks like a single side can be given multiple boundary IDs and correspondingly multiple boundary names. That seems a little peculiar to me but whatever. There are various versions of the function getBoundaryIDs, but the one I'm looking at is passed the arguments "boundary_names" and "true". True indicates that the getBoundaryIDs function needs to create ids for the boundary_names being passed because in our case we are creating a brand new boundary that doesn't have any ID yet.


add_side 	( 	const dof_id_type  	elem,
		const unsigned short int  	side,
		const boundary_id_type  	id
	)

Add side side of element number elem with boundary id id to the boundary information data structure.

std::string& libMesh::BoundaryInfo::sideset_name 	( 	boundary_id_type  	id	)

Returns a writable reference for setting an optional name for a sideset.


std::string & libMesh::BoundaryInfo::nodeset_name 	( 	boundary_id_type  	id	)

Returns a writable reference for setting an optional name for a nodeset.

The actions on boundary_info in SideSetsBetweenSubdomains are boundary.info.add_side and boundary_info.sideset_name.

Perhaps any time a boundary is created by adding sides, both a sideset and a nodeset are created.

Hell fffiinng yea, I think I'm starting to be able to read C++ and Moose better. I managed to successfully add a nodeset name!!

Alright, trying to add multiple sidesets at the same location does not work. What I imagine happening is that initially the sides are looped over and given the first name, then the same sides are looped over again and the first name is replaced with the second name. You have to think at the geomtric level! A set of sides is a sideset. A set of nodes is a nodeset. Two elements can share the same side. That happens at the interface between blocks left and right. Similarly, I envision that when I compound meshes in Salome, the meshes are physically joined such that we have the same case where elements on different blocks share the same side at the interface. Similarly we do not have an extra domain where the Lagrange multiplier variables live. The edge elements at the interface are the sides shared by the two blocks. Cannot create something out of nothing. Perhaps it's appropriate to think that a geomtric entity can only be used/named/given-an-id once. Perhaps.

Sideset and Nodeset names when viewed in paraview should absolutely match!!! This is because their name is taken from the boundary_name, which the sideset and nodeset should share.

In dump.txt:

      master                     = (required)                  # Master side ID
      slave                      = (required)                  # Slave side ID
      subdomain                  = (required)                  # Subdomain name that is the mortar interface

# Master side ID ==> More clearly, this is the boundary ID on the master side of the mortar interface
# Slave side ID ==> More clearly, this is the boundary ID on the slave side of the mortar interface

Left side slope: 3/0.5 = 6
Right side slope: -2/0.5 = -4

Difference in slopes = 6--4 = 10
Very weird hahaha

LHS slope = 5.5/0.5 = 11
RHS slope = -4.5/0.5 = -9
Difference in slopes = 11--9 = 20

So the computer's doing exactly what it's told, although it's apparently being told something quite different from the instructions that I thought I was giving. Apparently I was specifiying the change in the gradient...wouldn't that be kind of like the laplacian?? Whatever, this isn't what I'm trying to investigate.

An additional important note: it appears that only one BC can be specified for any given variable on a boundary. This would make sense since we are in the business of solving boundary value problems.

Either nodeset_name or sideset_name can be called in order to connect the boundary IDs with the boundary names in SideSetsBetweenSubdomains. And as soon as either of those naming functions has been called, the boundary ID-name connection has been made. Then regardless of whether nodeset_name or sideset_name was used, both IntegratedBC and NodalBC can access the required respective sideset or nodeset using the boundary name -> boundaryID -> sideset or nodeset connected to that boundaryID.

It's funny, sideset_name can be called even if you didn't add a sideset. Like with AddExtraNodeset, no sideset was added (you can't find one in the mesh file), but you can call sideset_name to make the connection between boundaryIDs and boundaryNames. That's the only purpose of sideset_name and nodeset_name that I can see, is making the connection between boundary_ids and boundary_names. Almost seems redundant to have both methods. Obviously sideset_name will only added a named sideset entity in the mesh file if a sideset was actually created.

As expected trying to implement a NeumannBC on a boundary that has no corresponding sideset does nothing. It's very important to remember then that if a boundary ID/name exists but there's not a corresponding sideset/nodeset, then MOOSE will not tell you that anything is wrong, but it will not enforce your Integrated/Nodal BC.

Very neat conversation on an Elmer forum originating from a guy wanting to do something very similar to me (includes mention of DG, mortar, mesh interfacing, etc.)

In the conforming_two_var_out.i test, u ranges from 0 to 1 and v ranges from 0 to 0.0388393. Do the lagrange multipliers actually play a role in the dependent variable solution? lm_u ranges from -4.68841 to 2.49956 and lm_v ranges from -0.353036 to 0.127077.

My test conclusively shows that my .unv file is not being read correctly. Specifically, I can change the mortar subdomain from 'Interface' (a group of edges on my mesh in Salome that should exist when read in by Moose) to 'Jabroni' (a totally fictional non-existent thing) and the problem computes the exact same results.

When creating libraries that may eventually be linked, The -fPIC flag enables the compilers to create position independent code, needed for shared libraries in Ubuntu on a 64 bit Intel processor. So in order to make exodus link with netcdf (I compiled netcdf before exodus), I had to make sure I configured netcdf to use CFLAGS=CXXFLAGS=CPPFLAGS=F90FLAGS=FFLAGS="-fPIC". Just using CFLAGS="-fPIC" did not do the trick. My guess would be that CPPFLAGS was the necessary addition, but I'm not going to spend the time to go check it out. The PIC in -fPIC stands for position independent code.

I believe first declaration of _all_objects (found by using my etags table) is in Warehouse.h. Here's the declaration line:

  /// All instances of objects (raw pointers)
  std::vector<T *> _all_objects;

Declaration vs definition: When you declare a variable, a function, or even a class all you are doing is saying: there is something with this name, and it has this type. Defining something means providing all of the necessary information to create that thing in its entirety. Defining a function means providing a function body; defining a class means giving all of the methods of the class and the fields. Once something is defined, that also counts as declaring it; so you can often both declare and define a function, class or variable at the same time. But you don't have to.

For example, having a declaration is often good enough for the compiler.

There are twelve tags of _dg_kernels. The declaration/definition is on line 508 of include/base/NonlinearSystem.h:

  std::vector<DGKernelWarehouse> _dg_kernels;

On line 655 of NonlinearSystem.C in the function NonlinearSystem::addDGKernel, we have _dg_kernels[tid].addDGKernel(dg_kernel)

So here's the thing: there are three definitions of addDGKernel. There's FEProblem::addDGKernel, NonlinearSystem::addDGKernel, and DGKernelWarehouse::addDGKernel. Which definition gets called?

All the operations on the _dg_kernels variable occur in the NonlinearSystem.C file (there are 11 operations). _dg_kernels is a protected member of the NonlinearSystem class.

Class variables and functions are both class "members." Members of a derived class can access protected members of the base class but not private members of the base class. Members who aren't in the derived or base class cannot access protected or private members of the base class. Public members can be accessed by members outside the class.

The keyword "this" represents a pointer to the object whose member function is being executed. It is used within a class's member function to refer to the object itself.

After the declarations of Rectangle and rect, any of the public members of object rect can be accessed (in Main or in other classes) as if they were normal functions or normal variables, by simply inserting a dot (.) between object name and member name. The only members of rect that cannot be accessed from outside the class are width and height, since they have private access and they can only be referred to from within other members of that same class.


    private members of a class are accessible only from within other members of the same class (or from their "friends").
    protected members are accessible from other members of the same class (or from their "friends"), but also from members of their derived classes.
    Finally, public members are accessible from anywhere where the object is visible.

Recall that _dg_kernels is a vector. This vector contains values of type DGKernelWarehouse. So _dg_kernels is a vector of objects of type (class) DGKernelWarehouse.

The class DGKernelWarehouse has a public member function updateActiveDGKernels(Real t, Real dt). It makes sense that this is public, otherwise we couldn't execute the command:

_dg_kernels[tid].updateActiveDGKernels(t, dt);

on line 2400 of NonlinearSystem.C (there is no friendship between DGKernelWarehouse and NonlinearSystem)

Note that _dg_kernels[tid] simply has the type DGKernelWarehouse because we've accessed an index of the vector.

void
DGKernelWarehouse::addDGKernel(MooseSharedPointer<DGKernel> & dg_kernel)
{
  _all_ptrs.push_back(dg_kernel);
  _all_objects.push_back(dg_kernel.get());
}

Where are the calls to addDGKernel?

void
FEProblem::addDGKernel(const std::string & dg_kernel_name, const std::string & name, InputParameters parameters)
{
  parameters.set<FEProblem *>("_fe_problem") = this;
  if (_displaced_problem != NULL && parameters.get<bool>("use_displaced_mesh"))
  {
    parameters.set<SubProblem *>("_subproblem") = _displaced_problem.get();
    parameters.set<SystemBase *>("_sys") = &_displaced_problem->nlSys();
    _reinit_displaced_face = true;
  }
  else
  {
    parameters.set<SubProblem *>("_subproblem") = this;
    parameters.set<SystemBase *>("_sys") = &_nl;
  }
  _nl.addDGKernel(dg_kernel_name, name, parameters);
}

So this is from the FEProblem class and it calls the addDGKernel function of the NonlinearSystem class. As expected, the addDGKernel member of the NonlinearSystem class is public. (It had to be.)

void
AddDGKernelAction::act()
{
  _problem->addDGKernel(_type, _name, _moose_object_pars);
}

I'm guessing that this calls the addDGKernel function in the FEProblem class.

void
NonlinearSystem::addDGKernel(std::string dg_kernel_name, const std::string & name, InputParameters parameters)
{
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    parameters.set<MaterialData *>("_material_data") = _fe_problem._bnd_material_data[tid];
    parameters.set<MaterialData *>("_neighbor_material_data") = _fe_problem._neighbor_material_data[tid];

    MooseSharedPointer<DGKernel> dg_kernel = MooseSharedNamespace::static_pointer_cast<DGKernel>(_factory.create(dg_kernel_name, name, parameters, tid));

    _dg_kernels[tid].addDGKernel(dg_kernel);
  }

  _doing_dg = true;
}

So recall that programming is logical and sequential. How do I think adding a DGKernel progresses at this point?:

AddDGKernelAction::act() ==> FEProblem::addDGKernel ==> NonlinearSystem::addDGKernel ==> DGKernelWarehouse::addDGKernel

There are of course 63 references to the function act(). There are also 57 definitions. All of the definitions are in the src/actions directory.

6 extra references (references that aren't simply declarations). 5 of these are in the ActionWarehouse.C. One is in the parent class Action.h

foo->bar is equivalent to (*foo).bar, i.e. it gets the member called bar from the struct that foo points to.

So foo points to an object of a class. foo->bar gets the member called bar from the class object that foo points to.

class Action : public ConsoleStreamInterface
{
public:
  Action(InputParameters parameters);
  virtual ~Action() {}                  // empty virtual destructor for proper memory release


act is a point that points to type/class Action

What are we doing when we do act->type(); ??

What this looks like to me is an attempt to access the member function type() of the class Action. The problem is that there is no member function type() in the class Action. However in the constructor of Action we have the line   Action(InputParameters parameters);

And the class InputParameters does have a member function type that returns std::string. However, it also accepts const std::string &name, and there doesn't appear to be any passing in act->type();

Man this is a heavy dose of c++.

      for (std::vector<Action *>::const_iterator k = _action_blocks.at(*j).begin(); k != _action_blocks.at(*j).end(); ++k)

k is a constant iterator over a vector of pointers that point to actions. That's my attempt at translating. Don't know whether it's right.

C++ STL is the C++ Standard Template Library.

std::vector is a sequence container that encapsulates dynamic size arrays

An interator is a pointer. You can think of an iterator as pointing to an item that is part of a larger container of items.  For instance, all containers support a function called begin, which will return an iterator pointing to the beginning of the container (the first element) and function, end, that returns an iterator corresponding to having reached the end of the container. In fact, you can access the element by "dereferencing" the iterator with a *, just as you would dereference a pointer.

A pointer is a variable which stores the address of another variable.

int firstvalue, secondvalue;
  int * mypointer;

  mypointer = &firstvalue;
  *mypointer = 10;


    & is the address-of operator, and can be read simply as "address of"
    * is the dereference operator, and can be read as "value pointed to by"

        Action * act = *k;

The above line means that the Action pointed to by act is equal to the Action pointed to by the const_iterator k

I'm conducting a sleuth hunt. There are four places where type() is defined: MaterialProperty.h, RestartableData.h, ArbitraryQuadrature.C, and InputParameters.C

global -v --result=grep --color=always --path-style=shorter --from-here=149:src/actions/ActionWarehouse.C type
include/materials/MaterialProperty.h:216:MaterialProperty<T>::type ()
include/restart/RestartableData.h:142:RestartableData<T>::type ()
src/utils/ArbitraryQuadrature.C:28:ArbitraryQuadrature::type() const
src/utils/InputParameters.C:444:InputParameters::type(const std::string &name)
4 objects located (using '/home/lindsayad/projects/moose/framework/GTAGS').

Global found 4 definitions at Mon Oct  5 16:59:56

Type Action does not have a "definition" of a member type() as determined by ggtags, but a quick search in the Action.h file reveals:

  const std::string & type() const { return _action_type; }

ggtags likely considers this to be a declaration as opposed to a definition, and unfortunately there are 783 references to type(). However, we know that the member type() either has to be defined (in the real sense, not the ggtags interpretation), in the Action class or in one of the parent classes of Action. We have to exercise our knowledge!

Action inherits from ConsoleStreamInterface.

act->type is executed in ActionWarehouse.C

ActionWarehouse inherits from Warehouse. If we assume that Action and ActionWarehouse don't share any inheritance, then we could have deduced that the member function type() of class Action had to be public (which it is).

All this knowledge gathering stems from me trying to determine when AddDGKernelAction::act() gets called. It has to get called in

void
ActionWarehouse::executeActionsWithAction(const std::string & task)
{
  // Set the current task name
  _current_task = task;

  for (ActionIterator act_iter = actionBlocksWithActionBegin(task);
       act_iter != actionBlocksWithActionEnd(task);
       ++act_iter)
  {
    if (_show_actions)
    {
      _console << "[DBG][ACT] "
               << "TASK (" << COLOR_YELLOW << std::setw (24) << task << COLOR_DEFAULT << ") "
               << "TYPE (" << COLOR_YELLOW << std::setw (32) << (*act_iter)->type() << COLOR_DEFAULT << ") "
               << "NAME (" << COLOR_YELLOW << std::setw (16) << (*act_iter)->getShortName() << COLOR_DEFAULT << ") ";

      (*act_iter)->act();
    }
    else
      (*act_iter)->act();
  }
}

I'm looking at actionBlocksWithActionBegin(task) ==>

ActionIterator
ActionWarehouse::actionBlocksWithActionBegin(const std::string & task)
{
  return _action_blocks[task].begin();
}

  /// Pointers to the actual parsed input file blocks
  std::map<std::string, std::vector<Action *> > _action_blocks;

So looks like we've gone all the way back to the parsing of the input file (which is one of the things I wanted to do!). I knew that the problem was translated from the input file into the Moose executable and problem solving using Actions. (Actions are input blocks right? I mean yes they bloody are. When I created new actions, I created new input blocks, with proper input block syntax specified. Yes input file syntax is done in ZapdosApp::associateSyntax, and one of the arguments we pass is action_factory.)

executeActionsWithAction is referenced in the following places:

include/actions/ActionWarehouse.h:119:  void executeActionsWithAction(const std::string & name);
src/actions/ActionWarehouse.C:308:    executeActionsWithAction(task);
src/base/MooseApp.C:457:  _action_warehouse.executeActionsWithAction("set_global_params");
src/base/MooseApp.C:458:  _action_warehouse.executeActionsWithAction("setup_mesh");
src/base/MooseApp.C:459:  _action_warehouse.executeActionsWithAction("prepare_mesh");
src/base/MooseApp.C:460:  _action_warehouse.executeActionsWithAction("add_mesh_modifier");
src/base/MooseApp.C:461:  _action_warehouse.executeActionsWithAction("execute_mesh_modifiers");
src/base/MooseApp.C:462:  _action_warehouse.executeActionsWithAction("uniform_refine_mesh");
src/base/MooseApp.C:463:  _action_warehouse.executeActionsWithAction("setup_mesh_complete");

Something that's interesting to note is that in:

ActionWarehouse::executeAllActions()
{
  if (_show_actions)
  {
    _console << "[DBG][ACT] Action Dependency Sets:\n";
    printActionDependencySets();

    _console << "\n[DBG][ACT] Executing actions:" << std::endl;
  }

  for (std::vector<std::string>::iterator it = _ordered_names.begin(); it != _ordered_names.end(); ++it)
  {
    std::string task = *it;
    executeActionsWithAction(task);
  }
}

the executeActionsWithAction member function of ActionWarehouse is called without any '.' or '->' or '::'. This is because there is no ambiguity about what namespace we are using because we are within a member function of ActionWarehouse calling a fellow member function of ActionWarehouse.

It looks like the actions of adding kernels is conducted during mesh preparation, and this process is initiated from lines 457-463 of MooseApp.C, which is within the meshOnly member function of MooseApp

meshOnly is executed only if user specifies --mesh-only when running his executable from the command line. Thus the important call of executeActionsWithAction happens in executeAllActions() which is called from MooseApp::runInputFile(), which is called from MooseApp::run(), which is called from these lines in the main.C file of Zapdos:

  // This creates dynamic memory that we're responsible for deleting
  MooseApp * app = AppFactory::createApp("ZapdosApp", argc, argv);

  // Execute the application
  app->run();

Wooooo!!!!!!!!!! We traced it all the way back to the beginning using the power of the ggtags!!!!!! Good day :-) The power of C++ is being channelled into my mind!!!!!!!!

There are three places where computeResidual gets called in ComputeResidualThread: in onElement (kernels), onBoundary (integrated BCs), and onInternalSide (DGKernels)

In ThreadedElementLoopBase we have the calls:

      onElement(elem);

      for (unsigned int side=0; side<elem->n_sides(); side++)
      {
        std::vector<BoundaryID> boundary_ids = _mesh.boundaryIDs(elem, side);

        if (boundary_ids.size() > 0)
          for (std::vector<BoundaryID>::iterator it = boundary_ids.begin(); it != boundary_ids.end(); ++it)
            onBoundary(elem, side, *it);

        if (elem->neighbor(side) != NULL)
          onInternalSide(elem, side);
      } // sides
      postElement(elem);

Need to investigate how ComputeResidualThread::onElement might get called here

Or need to say if there are other places where computeResidual gets called

There are 14 definitions of onElement as defined by ggtags:

include/base/ThreadedElementLoopBase.h:199:ThreadedElementLoopBase<RangeType>::onElement(const Elem * /*elem*/)
src/base/CacheChangedListsThread.C:34:CacheChangedListsThread::onElement(const Elem *elem)
src/base/ComputeDampingThread.C:44:ComputeDampingThread::onElement(const Elem *elem)
src/base/ComputeDiracThread.C:71:ComputeDiracThread::onElement(const Elem * elem)
src/base/ComputeElemAuxVarsThread.C:75:ComputeElemAuxVarsThread::onElement(const Elem * elem)
src/base/ComputeIndicatorThread.C:83:ComputeIndicatorThread::onElement(const Elem *elem)
src/base/ComputeJacobianThread.C:148:ComputeJacobianThread::onElement(const Elem *elem)
src/base/ComputeMarkerThread.C:76:ComputeMarkerThread::onElement(const Elem *elem)
src/base/ComputeMaterialsObjectThread.C:78:ComputeMaterialsObjectThread::onElement(const Elem *elem)
src/base/ComputeResidualThread.C:99:ComputeResidualThread::onElement(const Elem *elem)
src/base/ComputeUserObjectsThread.C:164:ComputeUserObjectsThread::onElement(const Elem * elem)
src/base/FlagElementsThread.C:57:FlagElementsThread::onElement(const Elem *elem)
src/base/ProjectMaterialProperties.C:75:ProjectMaterialProperties::onElement(const Elem *elem)
src/base/UpdateErrorVectorsThread.C:57:UpdateErrorVectorsThread::onElement(const Elem *elem)

kernel->computeResidual() in NonlinearSystem::computeResdiaulInternal(Moose::KernelType type) which is called from:

void
NonlinearSystem::computeResidual(NumericVector<Number> & residual, Moose::KernelType type)

computeResidual is called in several places as we have seen.

Let's understand polymorphism. Looking at our polymorphism example, ppoly1->area()

ppoly is defined as a pointer that points to an object of type Polygon, that in this case is assigned to the address of the the object Rectangle rect (this assignment is fine because rect is a Rectangle but it is also a Polygon), so ppoly1->area() accesses the area() member function of Rectangle rect. Because ppoly is defined as a pointer that points to an object of type Polygon, only members inherited from Polygon can be accessed using -> . Member functions defined in a base class but taggged with virtual can be re-defined in a derived class, and then those derived class function members can be accessed through a pointer to the base class. However, if the virtual tag was not given in the base class, then a member function with the same name in a derived class would be unaccessible through a pointer to the base class.

computeResidual calls:

src/base/ComputeDiracThread.C:113:          dirac_kernel->computeResidual(); ==> DiracKernel::computeResidual()
src/base/ComputeResidualThread.C:114:    (*it)->computeResidual(); ==> it points to type KernelBase. ==> KernelBase::computeResidual()
src/base/ComputeResidualThread.C:144:        bc->computeResidual(); ==> IntegratedBC::computeResidual
src/base/ComputeResidualThread.C:180:        dg->computeResidual(); ==> DGKernel::computeResidual
src/base/FEProblem.C:3281:  computeResidual(sys, u, residual); ==> FEProblem::computeResidual
src/base/FEProblem.C:3330:  _nl.computeResidual(residual, type); ==> _nl is of type NonlinearSystem, thus this calls ==> NonlinearSystem::computeResidual
src/base/FEProblem.C:3386:    computeResidual(sys, soln, *sys.rhs); ==> FEProblem::computeResidual
src/base/NonlinearSystem.C:86:    p->computeResidual(sys, soln, residual); ==> FEProblem::computeResidual
src/base/NonlinearSystem.C:233:    _fe_problem.computeResidual(_sys, *_current_solution, *_sys.rhs); ==> FEProblem::computeResidual
src/base/NonlinearSystem.C:866:        nc->computeResidual(residual); ==> nc points to a NodalConstraint ==> NodalConstratin::computeResidual
src/base/NonlinearSystem.C:1073:                nfc->computeResidual(); ==> NodeFaceConstraint::computeResidual
src/base/NonlinearSystem.C:1159:          ffc->computeResidual(); ==> FaceFaceConstraint:computeResidual
src/base/NonlinearSystem.C:1228:        kernel->computeResidual(); ==> ScalarKernel::computeResidual()
src/base/NonlinearSystem.C:1309:            bc->computeResidual(residual); ==> NodalBC::computeResidual()
src/bcs/NodalNormalBC.C:49:  NodalBC::computeResidual(residual);

computeResidual complete definitions (as defined by ggtags):

src/base/FEProblem.C:3263:FEProblem::computeResidual(NonlinearImplicitSystem &/*sys*/, const NumericVector<Number> & soln, NumericVector<Number> & residual)
src/base/NonlinearSystem.C:702:NonlinearSystem::computeResidual(NumericVector<Number> & residual, Moose::KernelType type)
src/bcs/IntegratedBC.C:102:IntegratedBC::computeResidual()
src/bcs/NodalBC.C:73:NodalBC::computeResidual(NumericVector<Number> & residual)
src/bcs/NodalNormalBC.C:45:NodalNormalBC::computeResidual(NumericVector<Number> & residual)
src/constraints/FaceFaceConstraint.C:122:FaceFaceConstraint::computeResidual()
src/constraints/NodalConstraint.C:48:NodalConstraint::computeResidual(NumericVector<Number> & residual)
src/constraints/NodeFaceConstraint.C:108:NodeFaceConstraint::computeResidual()
src/dgkernels/DGKernel.C:130:DGKernel::computeResidual()
src/dirackernels/DiracKernel.C:81:DiracKernel::computeResidual()
src/kernels/EigenKernel.C:71:EigenKernel::computeResidual()
src/kernels/Kernel.C:47:Kernel::computeResidual()
src/kernels/KernelGrad.C:37:KernelGrad::computeResidual()
src/kernels/KernelValue.C:37:KernelValue::computeResidual()
src/kernels/NodalEqualValueConstraint.C:46:NodalEqualValueConstraint::computeResidual()
src/kernels/ODEKernel.C:40:ODEKernel::computeResidual()
src/kernels/TimeKernel.C:34:TimeKernel::computeResidual()

Kernel is a direct child of KernelBase

I want to understand this function:

void
ComputeResidualThread::onElement(const Elem *elem)
{
  _fe_problem.prepare(elem, _tid);
  _fe_problem.reinitElem(elem, _tid);
  _fe_problem.reinitMaterials(_subdomain, _tid);

  const std::vector<KernelBase *> * kernels = NULL;
  switch (_kernel_type)
  // _kernel type is obviously a member variable of ComputeResidualThread and it is of type Moose::KernelType
  // _kernel_type is initialized in the constructor
  {
  case Moose::KT_ALL: kernels = & _sys.getKernelWarehouse(_tid).active(); break;
  case Moose::KT_TIME: kernels = & _sys.getKernelWarehouse(_tid).activeTime(); break;
  case Moose::KT_NONTIME: kernels = & _sys.getKernelWarehouse(_tid).activeNonTime(); break;
  }
  for (std::vector<KernelBase *>::const_iterator it = kernels->begin(); it != kernels->end(); ++it)
  {
    (*it)->computeResidual();
  }

  _fe_problem.swapBackMaterials(_tid);
}

The constructor tells us how to construct the class. (How about that ?!)

As an example:

ComputeResidualThread::ComputeResidualThread(FEProblem & fe_problem, NonlinearSystem & sys, Moose::KernelType type) :
    ThreadedElementLoop<ConstElemRange>(fe_problem, sys),
    _sys(sys),
    _kernel_type(type),
    _num_cached(0)
{
}

So when an object of type ComputeResidualThread is declared in another piece of code, it must be initialized with three arguments matching: (FEProblem & fe_problem, NonlinearSystem & sys, Moose::KernelType type)

Now ComputeResidualThread has three data members: _sys, _kernel_type, and _num_cached. It also publicly inherits the members of ThreadedElementLoop<ConstElemRange>. The constructor tells us that _sys should be set equal to sys passed by the user, _kernel_type should be set equal to tyoe passed by the user, and _num_cached should be set equal to zero (not passed by the user). Finally the members of ThreadedElementLoop can be created by passing the user supplied variables fe_problem and sys.

Now I recognize that in order to be able to use _kernel_type in ComputeResidualThread::onElement, an object of type ComputeResidualThread must have been created prior to that such that _kernel_type has been initialized.

And in fact, there is one place where an ComputeResdiaul object is created: in NonlinearSystem::computeResidualInternal(Moose::KernelType type) :

    ComputeResidualThread cr(_fe_problem, *this, type);

NonlinearSystem::computeResidualInternal is called from NonlinearSystem::computeResidual

_nl must be a data member of FEProblem. Sure enough on line 904 of FEProblem.h, I see the protected variable declaration:

NonlinearSystem & _nl;

Here's where _nl is initialized in the constructor of FEProblem:

     _nl(getParam<bool>("use_nonlinear") ? *(new NonlinearSystem(*this, name_sys("nl", _n))) : *(new EigenSystem(*this, name_sys("nl", _n)))),

The keyword this represents a pointer to the object whose member function is being executed. It is used within a class's member function to refer to the object itself. The only object that the FEProblem::FEProblem constructor reads in is const InputParameters & parameters.

I believe that * in the _nl case is dereferencing.

What does the initialization statement actually mean? It means this:

NonlinearSystem & _nl = getParam<bool>("use_nonlinear") ? *(new NonlinearSystem(*this, name_sys("nl", _n))) : *(new EigenSystem(*this, name_sys("nl", _n)))

Dynamic memory is allocated using operator new. new is followed by a data type specifier and, if a sequence of more than one element is required, the number of these within brackets []. It returns a pointer to the beginning of the new block of memory allocated. Its syntax is:

pointer = new type
pointer = new type [number_of_elements]

C++ references allow you to create a second name for the a variable that you can use to read or modify the original data stored in that variable. While this may not sound appealing at first, what this means is that when you declare a reference and assign it a variable, it will allow you to treat the reference exactly as though it were the original variable for the purpose of accessing and modifying the value of the original variable--even if the second name (the reference) is located within a different scope. This means, for instance, that if you make your function arguments references, and you will effectively have a way to change the original data passed into the function. This is quite different from how C++ normally works, where you have arguments to a function copied into new variables. It also allows you to dramatically reduce the amount of copying that takes place behind the scenes, both with functions and in other areas of C++, like catch clauses.

So _nl is a reference!

Ampersand on the LHS of an assignment ==> creating a reference
Ampersand on the RHS of an assignment ==> creating a pointer

Recall the NonlinearSystem constructor:

NonlinearSystem::NonlinearSystem(FEProblem & fe_problem, const std::string & name)

More assignments:

FEProblem & fe_problem = *this
const std::string & name = name_sys("nl", _n)

Alright so 'this' in this context is the pointer to the FEProblem object.

So when we construct a FEProblem object, we create a new NonlinearSystem. This NonlinearSystem and FEProblem are closely coupled because when we create the new NonlinearSystem, we pass it the FEProblem object that we're concurrently constructing. name_sys is a simple function that returns type std::string. I don't know whether it has to be included in the FEProblem.C file in order for the FEProblem contsructor to call it. It either has to be in that source file or in another included file. (It has to be somewhere!!)

The fe_problem might be created on line 338 of MooseApp.C within the runInputFile member function:

    MooseSharedPointer<FEProblem> fe_problem= _action_warehouse.problem();

_action_warehouse is a data member of MooseApp. It's initialized in the following way:

    _action_warehouse(*this, _syntax, _action_factory),

'this' is a pointer to a MooseApp object

ActionWarehouse::ActionWarehouse(MooseApp & app, Syntax & syntax, ActionFactory & factory)

ActionWarehouse (or one of its parents) must have a member function problem()

Here it is:

  MooseSharedPointer<FEProblem> & problem() { return _problem; }

This returns a reference to a MooseSharePointer<FEProblem> called _problem

I'm a little perplexed by the ActionWarehouse class. It has a protected data member:

MooseSharedPointer<FEProblem> _problem

Does that data member get initialized by the constructor?? Ah I see. It's not illegal to create data members that aren't initialized by constructors. It's just that if we try to do something with those declared variables before we define them we will get an undetermined result.

Forward declaration

A declaration of the following form
class-key attr identifier ;

Declares a class type which will be defined later in this scope. Until the definition appears, this class name has incomplete type. This allows classes that refer to each other:

class Vector; // forward declaration
class Matrix {
    // ...
    friend Vector operator*(const Matrix&, const Vector&);
};
class Vector {
    // ...
    friend Vector operator*(const Matrix&, const Vector&);
};

In ScalarKernel.h we have this:

  virtual void computeResidual() = 0;

There are children of the class ScalarKernel that define computeResidual and thus would override the virtual computeResidual() in ScalarKernel. On line 1225, I think it's conceivable that *it could dereference to children of ScalarKernels, and thus because of the use of the virtual computeResidual, those children's computeResidual function definitions could get called. The exact same reasoning is true for KernelBase. We have this line in KernelBase.h:

  virtual void computeResidual() = 0;

In fact of the all the computeResidual declarations in header files, only ScalarKernel and KernelBase have the = 0 assignments, such that the system is completely defined. This coding is beautiful.

Departing briefly from the beautiful moose framework, I am building my AVS talk.

I'm trying to compare energy and LFA formulations. In order to do that, I need to make sure I'm using corresponding transport and rate coefficients between the two cases. I'm going to use argon as my gas. I'm only going to consider elastic and ionization reactions. In order to make my energy formulation more closely aligned with physical reality, I need to examine the following parameters:

_Eiz_en ==> 15.76 eV (Lieberman)
_rate_coeff_ion_en ==> 5e-14 m^3/s (Lieberman)
_muem ==> TBD
_diffem ==> TBD
_muip ==> TBD
_diffip ==> TBD
_rate_coeff_elastic ==> Already correct (Lieberman)
Potential ==> -2.5 kV (Go experiment at atmosphere with argon)
Ballast resistor ==> 8.1 kOhm (Go experiment at atmosphere with argon)
Gap distance ==> 1-2 mm (Go experiment at atmosphere with argon)

Keep in mind that some of these params may have to be adjusted for the simulation (I have to be sure I can get convergence!). Whatever params I adjust though should be params that can also be adjusted for the LFA simulation. There may even be a physical reason to say reduce the applied potential or increase/decrease the area of the plasma such that it puts less stress on the numerical system. If a measure of stress on the system is the concentration of ions at the cathode, then there are several ways to decrease that stress:

Decrease applied potential
Increase ion mobility
Increase ballast resistance
Increase electrode/plasma area (==> this could derive a bit from what Go was saying, if a certain increase in curent is required and the potential is constant, then the diameter of the plasma must increase ==> Current = current density (should be constant?) * area (variable?))

I believe that in order for thunderbird to match a name to a email address, I must have received or sent an email from/to that person while I've had the thunderbird application installed on the particular computer. That's fair.

Again I note that bolos reproduces bolsig results. Figure 3.16 on pg. 80 of Lieberman looks kind of like a big sack of horse-poo compared to the results I predict with bolos or bolsig. Perhaps this is because of the 2-term approximation used in bolsig and currently in my bolos calculations. But it seems like Lieberman predicts much higher rates of ionization than bolos/bolsig does. What to do?

See here's a conundrum. The current functional form I've been using for my source terms is alpha s porportional to exp(-A/|Field|) or Kiz is proportional to exp(-B/|Te|). Neither of those functional forms are actually physically correct, e.g. both those terms continue to grow respectively with the field and the electron temperature. They most certainly do not plateau at some threshold value. I don't want to take the time to try and navigate any possible additional complexity introduced by lookup tables right now. There's also the pain in the ass realization that if I alter the functional forms of alpha and Kiz, I'm going to have to implement new Jacobians and new residuals.

For a p=5 polynomial fit of EField vs. alpha, without a penalty for negative model values (with initial guess of 1.e-14 for all coefficients): (a0,a1....,a5)
[  1.00000000e-14,  -7.33794821e-03,   1.02893016e-08,
        -1.50362990e-15,   1.09634198e-22,  -3.16558448e-30]

With a penalty for negative values:
[  1.00000000e-14,  -7.23069891e-03,   9.41410909e-09,
        -1.18445467e-15,   6.90018879e-23,  -1.43863034e-30]

Now I want to try fitting A*T^B*exp(-C/T)

Alright, I got a frieking awesome curve fit with the Lieberman function. Fit values for the electric field vs alpha:

[  1.43171672e-01,   9.05925536e-01,   3.04958892e+06]

Using the stupid function. Fit values for the electric field vs alpha:

array([  485001.06154102,  8073212.48740293])

Curve fit for mean_energy vs. kArIz, Lieberman function:
array([  1.43878529e-11,  -2.70610234e-01,   7.64727794e+01])

Curve fit for mean_energy vs. kArIz, stupid function:
array([  6.24642728e-12,   7.43060838e+01])

The mean_energy curve fit is not nearly as good as the EField curve fit. And the exponential factor of 76 Volts makes absolutely no sense. It could be because we're not assuming a Maxwellian distribution??

With a mesh size of 100, solve failed at 1.29159e-08 seconds because of oscillations in the ions. With 200 elements the result is even worse. Oscillations immediately. Solve fails after 3 steps at 3.64e-09 seconds.

With my sweet spot of 4000 elements, still see some spectacular oscillations. Solve fails at 4.54e-09 seconds.

I went back to 100 elements, and decreased the voltage by a factor of 100. This time the onset of oscillations in the ions was a little delayed but they did indeed appear, and the solve failed at 2.25989e-08 sconds. Going back to 4000 elements, still with the decreased voltage: oscillations witnessed and solve fails at 1.08299e-08 seconds.

The reason for the above failures is that I had a negative diffusion coefficient for the ions: dope!! With that fixed, the solve appears to start chugging along fine until it reaches around 1e-6 seconds (when the potential is rising at a rapid rate towards its max value), and then I get a terrible MUMPS error:

[3]PETSC ERROR: [0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Error reported by MUMPS in numerical factorization phase: INFO(1)=-1, INFO(2)=1

[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[0]PETSC ERROR: Petsc Release Version 3.6.0, Jun, 09, 2015
[0]PETSC ERROR: ./zapdos-opt on a arch-linux2-c-opt named lindsayad-OptiPlex-990 by lindsayad Mon Oct 12 16:05:03 2015
[0]PETSC ERROR: Configure options --prefix=/opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt --download-hypre=1 --with-ssl=0 --with-debugging=no --with-pic=1 --with-shared-libraries=1 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --download-fblaslapack=1 --download-metis=1 --download-parmetis=1 --download-superlu_dist=1 --download-scalapack=1 --download-mumps=1 CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 CFLAGS="-fPIC -fopenmp" CXXFLAGS="-fPIC -fopenmp" FFLAGS="-fPIC -fopenmp" FCFLAGS="-fPIC -fopenmp" F90FLAGS="-fPIC -fopenmp" F77FLAGS="-fPIC -fopenmp" PETSC_DIR=/home/lindsayad/projects/src/petsc-3.6.0


My translation: the problem started getting very hard.

With the potential decreased by a factor of 10, we appear to be below the ionization threshold.

Alright potential halved from experimental value to 1.25e3. Then we are above the ionization threshold, but the solve still fails with the MUMPS error, this time at 8.53805e-06 seconds (so it made it a bit farther this time.)

After increasing the ballast resistance by a factor of two, absolutely nothing changed. Well that's not true. Although the time at which the simulation failed didn't change, the residuals did change. An important thing to note: the very previous time step, the solve converged in just 6 nonlinear iterations and reached a final residual of 1e-10.

By moving just a little farther in the time domain (went to a uniprocessor solve without MUMPS), e.g. solve failed at 8.67973e-06 seconds, I think I can see why the MUMPS simulation was failing. Looks like we're starting to get steep gradients at the boundary that may be unresolved by the mesh. Good to see the potential at the left bound starting to decrease. I think that means that our applied potential is not too large if we're able to simulate out a point where the potential starts decreasing.

Yea, I would say I've come quite a ways in being able to diagnose, identify, and rectify problems in my simulation. Less than an hour, I solved a brand new implementation of my DC LFA problem. Woooo!!!!!!!!!!!!!!!!!!!!!!!!!! Alive time = 28.6948 seconds. Active time = 28.4932 seconds.

Not looking good. Mean energy simulation dies at 2.49269e-07 seconds with no obvious problems in the solutions. Everything looks nice and smooth. Thus I think the problem must be the transient I'm introducing in the boundary condition.

After changing 1/tau to 1e5, the simulation now fails at 2.02422e-06 seconds. The potential reaches -250 V. I believe it was around -350 when using 1/tau = 1e6. The solutions still look pretty and smooth.

With potential scaling changed to 1e0 (initial value 1e-5), the solve fails at 1.2981e-06 seconds.
With potential scaling changed to 1e-1 (initial value 1e-5), the solve fails at 1.2981e-06 seconds.
With potential scaling changed to 1e-2 (initial value 1e-5), the solve fails at 1.2981e-06 seconds.
With potential scaling changed to 1e-3 (initial value 1e-5), the solve fails at 1.2981e-06 seconds.
With potential scaling changed to 1e-4 (initial value 1e-5), the solve fails at 2.02298e-06 seconds. (just a tiny bit worse than when scaling set to 1e5)

With absolute tolerance set to 1e-2, problem failed at 2.02924e-06 seconds. A little better.

mean_en scaling at 1e-23, solve fails at 1.22402-07 seconds (a_tol = 1e-1, tau = 1e-6)
mean_en scaling at 1e-20, solve fails at 1.22402-07 seconds (a_tol = 1e-1, tau = 1e-6)

The ability to solve and the result of the solve is very much dependent on the voltage rise time. For the lfa formulation, I could not get a converged solution for tau = 1e-3 & 1e-4 seconds, but I can for 1e-6 seconds. For the mean_energy formulation, I could not get a converged solution for tau = 1e-6 or 1e-5 seconds, but I got a converged solution for tau = 1e-3 seconds. What appeared to help the energy form was giving the charged species enough time to form a significant enough current such that the voltage at the left boundary never went too high. With respect to the lfa form, if I wait too long such that all the electrons and ions disappear from the domain before I apply the potential, then the problem won't converge.

I was also able to get convergence with lfa for tau = 1e-5 seconds. There appears to be no difference in the results between the 1e-6 and 1e-5 cases.

I'm unable to get convergence using a ballast resistance of 1e6, and a rise time of 1e-5 seconds (voltage of 1.25e3, 600, 300, or 150 V) with the electron energy formulation.

It simply appears that it's way easier to create ionization using my energy functional formulation than it is using my lfa functional formulation. I don't quite know why that is.

There are a couple of things that are different in my ionization source term between lfa and energy: lfa uses a coefficient that is a function of the electric field and multiplies the electron flux. Energy uses a coefficient that is a function of the electron energy and multiplies the electron density.

Got quite a few options for trying to fit alpha vs mean_energy:

Letting all three parameters of Liberman's be variable:
array([  7.00485613e+07,  -3.54236666e-02,   5.31700097e+01])

Setting c = 15.76:
array([ 42.04011787,   4.588541  ])

Setting b = 1, c = 15.76 (looks terrible):
array([ 78868.24191509])

If instead of deleting the elements in alpha that are less than a certain value (1e-15 in this case), I set those elements equal to zero, the heres the three parameter fit for the Lieberman function:
array([  1.52165930e+08,  -2.87277596e-01,   5.51972192e+01])

Fit virtually doesn't change, and I still see an overshoot of the data in the 6-7 eV range.

While I was revamping my code to allow for a source flux energy formulation, I realized that my electron energy kernel was wrong: I didn't have energy being lost per ionization collision: in fact I was actually gaining energy from ionization collisions! So hopefully I fixed that.

With the revamped code, the solve fails at 7.93437e-07 seconds which is much better performance relative to tau (1e-06). Reducing the rel_tol actually hurt the solve: solve failed at 7.03832e-07 seconds. Returning rel_tol to default and changin abs_tol from 1e-1 to 1e1, the solve reached 8.85438e-07 seconds.

From the Comsol Plasma Module manual, recall thaat the Electron energy can change on scales less than a ns.

I got the energy formulation to solve with a rise time of 1e-6 seconds. The key was simply to decrease dtmin. (Use your knowledge of the physical time scales!!! That's where being a physicist/engineer is of greatest benefit! Knowing the physics! (Time & length scales. Dependent variable scales also)). The shape of the 1e-3 and 1e-6 rise time simulations are exactly the same. And quantitatively, they are also quite similar. This would lead me to believe that the dominant path for energy loss is through elastic collisions. This would also be supported by the fact that I was able to solve with the 1e-3 rise time a problem where I had actually configured my residuals such that ionziation collisions would actually create energy!

Zapdos active time = 1063.31 seconds ==> 17.7 minutes. Not bad I would say. (using all 4 processors)

Went to the flux formulation of the energy formulation ==> Active time = 593.726 seconds = less than 10 minutes. Took 603 time steps as opposed to 1272 time steps in the case of the rate coefficient formulation of the energy formulation. I guess it is true that townsend in DC discharges is more stable than rate coefficients!!

For elastic collision rate fit:
array([  1.60638169e-13,   3.17917979e-01,   4.66301096e+00])

Ran with variable elastic collisions and the peak density actually increased :-( Zapdos active time = 575.026 seconds < 10 minutes (Had to change the absolute tolerance because it was stalling (diverged line searches after great initial residual reduction) at a few percent over 1e-5) The reason it was stalling at a slightly higher value is likely because the peak densities were a few percent higher. With the absolute tolerance change, the number of steps went from many thousands to 505.

Better fit for mean_en vs. kArIz:     y = p[0]*x**p[1]*exp(p[2]*x**p[3])

[  9.61505373e-15   1.56694111e-01  -3.62336868e+03  -3.39520538e+00]

If I can't get a good model fit...why not just use the raw data?!

There are Standard Stream Objects for Console I/O: (cout, cin, cerr, clog, etc.) A 'stream' is internally nothing but a series of characters.

std::cout - Regular output (console output)
std::cerr - Error output (console error)
std::clog - Nobody cares (console log)

When writing error messages, use cerr rather than cout. In simple examples, the two appear to be the same. However, they behave differently when you use Unix redirection operators when you run your program. In particular, if you redirect output to a file, error messages printed to cerr will still appear on the user's terminal whereas error messages printed to cout will be mixed into the output file.

So running the interpolation problem with Newton...after about a half hour we were at time = 7.69011e-07 seconds after 1141 time steps.

Very peculiarly, when I tried running two different zapdos simulations, each one on two cores, my process scheduler scheduled the sapdos simulations to both run on the same two cores! So two cores were just sitting there doing nothing!

As might be expected from the above passage, solving the same case as argon_energy_mumps_tau_1e-6_townsend_variable_elastic_coeff took about four times as long as the gold run. Not surprising since I was running on half the processors and theoretically it was only being executed half the time (two processes running on the two processors). Number of time steps were essentially the same.

With jacobian_estimation (and residual given by linear interpolation, so jacob and residual not the same):

Time Step 100, time = 7.42153e-07
                dt = 4.92521e-11
Time Step 200, time = 7.4892e-07
                dt = 9.47565e-11
Time Step 300, time = 7.53893e-07
                dt = 6.07676e-11
Time Step 400, time = 7.57823e-07
                dt = 3.89705e-11
Time Step 500, time = 7.64239e-07
                dt = 5.45626e-11
Time Step 600, time = 7.67041e-07
                dt = 2.82938e-12
Time Step 875, time = 7.68072e-07
                dt = 4.98006e-12

Without jacobian_estimation:

Time Step 100, time = 2.94026e-07
                dt = 2.30869e-09
Time Step 200, time = 3.55389e-07
                dt = 5.92228e-10
Time Step 300, time = 3.82268e-07
                dt = 3.79798e-10
Time Step 600, time = 4.64467e-07
                dt = 3.00513e-10
Time Step 800, time = 5.15558e-07
                dt = 2.14569e-10

So yes it looks much worse without the jacobian_estimation

With jacob_estimation and residual the same
Time Step 100, time = 8.34043e-07
                dt = 7.2944e-11
Time Step 200, time = 8.38448e-07
                dt = 3.48495e-11
Time Step 300, time = 8.40077e-07
                dt = 6.38499e-12
Time Step 400, time = 8.59315e-07
                dt = 2.98505e-09
Time Step 500, time = 0.0457816
                dt = 0.00763011

So as expected the problem solves better when the residual and jacobian match. Who would have guessed?! The problem with that match is that my residuals don't define the right problem! And when I do define the right problem with my residuals, I can't define a matching Jacobian!!

I tried defining a matching Jacobian my doing second order differencing mean_energy vs alpha. I am using a linear interpolation for alpha in the residuals, and a linear interpolation of d_alpha_d_mean_en in the jacobians. However, Andy Wilkins suggested that in his experience linear interpolations in conjuction with implicit solves can lead to terrible convergence. I am in fact seeing terrible convergence. He gave me a reference for cubic splines. My only question then is what to do for the derivatives.

What are the coefficients that are variable? a_k, b_k, c_k, K

0,1,2 3,4,5 6,7,8 8,10,11
12,13,14

So when I run through the debugger, I get the error: Solving for an empty SplineInterpolation.

This could be the problem that's causing MPI to abort immediately. Nope. It was because of getEnv.

A get function/method is defined in the Parameters class in files included from libmesh. It's that inheritance.

When I tried using set in DGInterface with _boundary_id, I got this warning:

/home/lindsayad/moose/framework/src/dgkernels/DGInterface.C:26:42: warning: overflow in implicit constant conversion [-Woverflow]
   params.set<BoundaryID>("_boundary_id") = 78910;

  /**
   * Inserts a new Parameter into the object but does not return
   * a writable reference.  The value of the newly inserted
   * parameter may not be valid.
   */
  template <typename T>
  void insert (const std::string&);

The Elem class (from libmesh) has a function level:

  /**
   * @returns the refinement level of the current element.  If the
   * element's parent is \p NULL then by convention it is at
   * level 0, otherwise it is simply at one level greater than
   * its parent.
   */
  unsigned int level () const;

And then elem inherits members of dof_object. dof_object has a member function id ():

  /**
   * \returns the \p id for this \p DofObject
   */
  dof_id_type id () const;

Still just kind of working and navigating through the code here.

Should I think about everyone or just think about myself?

I feel like my changes should just be limited to DGKernel parts of the framework, and even then the changes should be as minimal as possible.

I think that my DGKernel should mimic an IntegratedBC as much as possible.

It's clear that I don't know what combining CG and DG in MOOSE does.

Current task: figure out what's going on with EArray vs. alpha.

Nothing wrong with EArray vs. alpha. Fit is still perfection incarnate.

A pure CG problem with MONOMIAL basis functions will not converge; at least not when I was twiddling around with the 2d_diffusion_dg_test case. DGFunctionDiffusionDirichletBC contains the necessary DG residuals from diffusion. This BC must contain those contributions because DGDiffusion is only executed on internal sides; thus to be complete, those residual contributions must be added to the boundary condition. Then the epsilon and sigma terms are added. These almost exactly analagous to the internal epsilon and sigma terms, except u_neighbor is modified to be equal to the desired boundary value (e.g. the value of our boundary condition). Similarly, since there is no test function on the other side of the boundary, {test} = test, [test]=test, {normals*grad_test} = normals*grad_test.

Also, in MOOSE, any time a DG kernel is used for transport, a corresponding CG kernel must be used because the CG kernel actually computes the necessary element integrals (e.g. it adds the residual contribution of the element volume integrals.) DG appears to converge with both monomials and lagrange functions.

Penalty terms are essentially artificial dissipative terms. Reference: http://scicomp.stackexchange.com/questions/10928/in-what-regime-do-the-continuous-and-discontinuous-galerkin-method-become-unstab.

DG does not require that the solution be continuous at nodes, but rather introduces a jump penalty term that adds an amount to the residual proportional to the solution jump discontinuity. CG is a limit of DG for some particular case of epsilon and sigma (I can't recall what the settings are. I know it's just taking one of the parameters to infinity, sigma I believe).

DG does not require solution continuity, however, the solution can be continuous. Thus Lagrange shape functions can be used. CG requires continuity, so Lagrange shape functions are desirable because they are continuous at the nodes.

The idea of using argon excitation to bring down the electron density is a bust. It does bring down the electron density...but way too much. This makes sense because for all energies in the domain, the excitation coefficient is larger than the ionization coefficient. Thus excitation processes will occur at a greater rate, bringing down the electron energy such that there is never appreciable ionization. This is not a problem in LFA where an excitation event wouldn't destroy the quantity that ionization depends on.

This morning's analysis validates Bolos IMHO. Now can proceed with introducing mean_energy dependent transport and rate parameters.

Hagelaar's paper is truly an excellent paper. I lot of excellent mathematical physics content. In our Boltzmann solver, the electron energy distribution is a function of E/N which is a function of E. E varies over our domain, thus our EEDF varies over the domain. D(EEDF(E(z))) -> thus in equation 54 of Hagelaar (the drift-diffusion coefficient), the diffusion coefficient shouldn't be extracted from the partial derivative. But to make my technical job easier, I might do just that. By doing this I'm essentially straying from the true physics. My residuals represent my physics.

In my GeneralUserObject, when I try and compile it, I get these errors:

/home/lindsayad/projects/zapdos/src/userobjects/ResidAndJacobTerms.C:27:22: error: 'coupledValue' was not declared in this scope
   _u(coupledValue("u")),
                      ^
/home/lindsayad/projects/zapdos/src/userobjects/ResidAndJacobTerms.C:28:30: error: 'coupledGradient' was not declared in this scope
   _grad_u(coupledGradient("u"))

I do not apparently have access to the functions coupledValue and coupledGraident. However, ElementIntegralVariableUserObject, obviously has access to it. How and why does it have access?

Where is the coupledValue function defined? It's defined in Coupleable.C. Thus I think that ElementIntegralVariableUserObject must somehow be a child of Coupleable. Looking at the inheritance diagram, this is exactly what we see:

Coupleable <- ElementUserObject <- ElementIntegralUserObject <- ElementIntegralVariableUserObject.

BlockAverageValue inherits from ElementIntegralVariablePostprocessor. ElementIntegralVariablePostprocessor -> ElementIntegralPostprocessor -> ElementPostprocessor -> ElementUserObject & Postprocessor

computeIntegral is defined in ElementIntegralPostprocessor:

  Real sum = 0;
   60
   61   for (_qp=0; _qp<_qrule->n_points(); _qp++)
   62     sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();
   63   return sum;

computeQpIntegral() is defined in ElementIntegralVariablePostprocessor:

38  return _u[_qp];

Combining, computeIntegral executes:

  Real sum = 0;

      for (_qp=0; _qp<_qrule->n_points(); _qp++)
        sum += _JxW[_qp]*_coord[_qp]*_u[_qp];
      return sum;

Children can override parent functions.

gatherSum() is a function defined in UserObject.C and looks like this:

  _communicator.sum(value);

In the header file for libmesh::ParallelObject, there are three different options for constructing an object of type ParallelObject because there are three different constructors. That's all I'm going to write down for now. Have to walk a fine line between learning important c++ knowledge and not falling down a rabbit hole away from what I need to be focusing on.

I want to be able to access _phi and _grad_phi. I need _assembly and _fe_type (this comes from looking at NodalNormalsPreprocessor).

_assembly is in UserObject -> Covered
_fe_type: first defined in Nodal... so I have to figure out how Nodal... got it

  /// Transformed Jacobian weights
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;

Active time for Zapdos with reduced ion diffusivity, variable electron energy coeff, and variable ionization coeff = 1911.58 seconds = 32 minutes.

I need to figure out why I'm getting a ramped density in the middle of the domain for my electron energy formulation. It makes no sense to me. From PIC we know that the density profile in the center should be flat.

Just one small typo in a Jacobian statement (using _u insteand of _em in the Ion kernel), can mess things up. And that mess up may not appear in the simple -snes_type=test configuration.

After correcting the Jacobian, the problem solved in a standard amount of time: 439 time steps, 24MB output file.

Adding variable elastic coefficient: 442 time steps, active time = 698.753 seconds, output file size = 24 MB

Reducing the plasma radius (still high mobility and diffusivity ions): 455 steps, active time = 755.409 seconds, output file size = 25 MB. Reducing the plasma radius has the effect of increasing the charged particle densities (same potential, same ballast resistor ==> in order to have the same current with a decreased cross sectional area, the current density must increase), which means that the magnitude of the residuals increases. Consequently, in order to get the simulation to converge like it should, I had to increase the nl_abs_tol from 1e-4 to 3.14e-4. One could see the problem by looking at the logs: excellent decrease in the residuals until it hit around 1.1e-4 and then the residuals stopped decreasing. Newton convergence should be initially linear and then quadratic, so if I notice non-Newton behavior, it's because I've set the absolute tolerance too low.

Going with low transport values for ions now: 868 time steps, active time = 1561.72 seconds, output file = 47 MB. Peak densities are again higher, so in order for the simulation to solve in a reasonable time, I upped the nl_abs_tol to 1e-3. I really should investigate whether there is an intelligent way to set nl_abs_tol as opposed to manually "guessing"; manual "guessing" in which I am informed mostly by the log files.

Looking at a convergence series:

1.35e8, 2.74e7, 1.22e6, 5.59e3, 7.24e1, 1.35e1, 2.12e0, 7.38e-1

x = x_k
y = x_(k+1)

Functions which reference _bcs:
activeBoundaries() ==> This function has no references, so I'm not going to worry about implementing an analogue for _intcs (no references either through MOOSE doxygen or ggtags)
activeIntegrated() ==> Referenced by ComputeFullJacobianThread::computeFaceJacobian(), ComputeJacobianThread::computeFaceJacobian(), ComputeJacobianThread::onBoundary(), ComputeResidualThread::onBoundary(), ComputeJacobianThread::subdomainChanged(), and ComputeResidualThread::subdomainChanged(), and NonlinearSystem::needMaterialOnSide() ==> Hmm, a lot of work to do here lol
addBC() ==> NonlinearSystem::addBoundaryCondition()
getBCs() ==> Referenced by ComputeFullJacobianThread::computeFaceJacobian().
initialSetup() ==> NonlinearSystem::initialSetupBCs()
jacobianSetup()==> NonlinearSystem::computeJacobianInternal()
residualSetup()==> NonlinearSystem::computeResidualInternal()
timestepSetup()==> NonlinearSystem::timestepSetup()

What I was going to write here before I got sidetracked into figuring out how to not exit navigation mode with RET, was that when I search for references to addBC() with ggtags, I find references. On the MOOSE doxygen page, there were no references given.

I think I am trying to do too much. I think that the implementation process doesn't have to be as large as I was preparing to make it.

DGKernelWarehouse::active() just returns _active_dg_kernels (type std::vector<DGKernel *>)

There are only two places that _active_dg_kernels is modified, and they're both in DGKernelWarehouse::updateActiveDGKernels:

  for (std::vector<DGKernel *>::const_iterator it = _all_objects.begin(); it != _all_objects.end(); ++it)
  {
    DGKernel * dg_kernel = *it;
    // FIXME: add startTime/stopTime to DGKernel
//    if (dg_kernel->startTime() <= t + (1e-6 * dt) && dg_kernel->stopTime() >= t + (1e-6 * dt))
      _active_dg_kernels.push_back(dg_kernel);
  }

Now what about _all_objects ? Comes from DGKernelWarehouse::addDGKernel:

{
  _all_ptrs.push_back(dg_kernel);
  _all_objects.push_back(dg_kernel.get());
}

DGKernelWarehouse::addDGKernel is called from NonlinearSystem::addDGKernel:

    MooseSharedPointer<DGKernel> dg_kernel = MooseSharedNamespace::static_pointer_cast<DGKernel>(_factory.create(dg_kernel_name, name, parameters, tid));

    _dg_kernels[tid].addDGKernel(dg_kernel);

I think that the static_pointer_cast is useful for casting children types to parent type. That's just a hypothesis though.

Going from _nl.addDGKernel to _problem->addDGKernel(_type, _name, _moose_object_pars) in AddDGKernelACtion::act(). _type refers to the DGKernel child type, e.g. it's communicating what type of DGKernel to build. Anyways, all child types of DGKernel are lumped together into the _all_objects class member in DGKernelWarehouse.

_sys.getDGKernelWarehouse(_tid).active() in ComputeResidualThread will summon all the DGKernel children from the input file (unless I have something to say about it!!!)

_boundary_id is an attribute that all children of DGKernel will have.

boundary is not a parameter that's introduced in the leaf children of BoundaryCondition (I use leaf to refere to children that themselves have no children). It's a parameter defined farther up the tree. So the boundary parameter could be and probably should be defined in the base DGKernel. For identifying whether to use BoundaryID or BoundaryName refer to the BoundaryCondition. The "boundary" parameter is not set in BoundaryCondition but rather in the ancestor of BoundaryCondition: BoundaryRestrictableRequired, and in the ancestor of BoundaryRestrictableRequired: BoundaryRestrictable. The difference between the two ancestors is that "boundary" is a RequiredParam in for types inheriting from BoundaryRestrictableREquired, whereas "boundary" is an optional parameter for types inheriting from BoundaryRestrictable. So DGKernel should inherit the parameters of BoundaryRestrictable and it should inherits the attributes and functions of BoundaryRestrictable as well.

Unrelated to Moose, I'm looking at my cpu usage by Thunderbird. If I look at the top command in the command line, the cpu usage referes to the percent of one cpu being used. If I press (shift + i) while top is running, then I can show the percent of total cpus running. This number will match much more closely with the number shown by the "System Monitor" application. Number of emails will be different from gmail to Thunderbird because, gmail displays it's number count as the "number of conversations." Thunderbird displays every message individually.

Each folder on the gmail web page has a matching folder on Thunderbird. This makes sense because all of these folders are subfolders of the locked adlinds3@ncsu.edu folder.

A member that I think will be very important for DGKernel: _boundary_restricted

We will inherit this member from BoundaryRestrictable. By default it is initialized to be false, unless the user does pass the parameter boundary. There is a problem, however. _boundary_restricted is a private attribute. But this is not a problem because there is a public member function of BoundaryRestrictable that returns _boundary_restricted! The function name is boundaryRestricted()

Potentially important code in ComputeResidualThread::onBoundary:

    // Set the active boundary id so that BoundaryRestrictable::_boundary_id is correct
    _fe_problem.setCurrentBoundaryID(bnd_id);

In class declaration of BoundaryRestrictable, under the protected attributes:

  /// Reference to active boundary id
  const BoundaryID & _current_boundary_id;

if _bnd_feproblem == NULL, then _current_boundary_id is a constant reference (constant alias) to _invalid_boundary_id (a private member of BoundaryRestrictable of type const BoundaryID). If _bnd_feproblem != NULL, the _current_boundary_id is a constant reference (constant alias) to _bnd_feproblem->getCurrentBoundaryID().

_bnd_feproblem is a pointer to type FEProblem.

_invalid_boundary_id is assigned to be equal to Moose::INVALID_BOUNDARY_ID in the standard BoundaryRestrictable constructor.

I believe that _bnd_ids is the data member in BoundaryRestrictable that holds the boundary ids where the object that is of my DGInterface type (and thus is also of the BoundaryRestrictable type) is active. I would probably name my DGInterface object something stupid like dg_interface. _bnd_ids are the the boundary ids where dg_interface lives.

_bnd_ids is a private member of BoundaryRestrictable, which means that even children of BoundaryRestrictable are unable to access _bnd_ids. But once again, there is a public member function called boundaryIDs, which does this:

const std::set<BoundaryID> &
BoundaryRestrictable::boundaryIDs() const
{
  return _bnd_ids;
}

This returns a constant alias for a set of BoundaryID's. It's the _bnd_ids where my dg_interface lives. Recall that const at the end of a type/class's method/function indicates that the method does not modify any of the class's data members.

_current_boundary_id: protected member of BoundaryRestrictable ==> passed onto DGInterface
_bnd_ids: private member of BoundaryRestrictable ==> not passed onto DGInterface

There's a method that returns the boundaryID's that dg_interface lives on (inherited from BoundaryRestrictable), but not a method that returns the _current_boundary_id. This makes sense because the current_boundary_id is not really owned by dg_interface; it's more owned by the finite element problem (rotating around looking at different piecies of the problem, different boundaries, etc.). However, the boundaries where dg_interface lives is entirely owned by dg_interface.

In FEProblem class declaration:

  void setCurrentBoundaryID(BoundaryID id){ _current_boundary_id = id; }

Protected member of FEProblem class:

  BoundaryID _current_boundary_id;

In constructor of BoundaryRestrictable:

    _current_boundary_id(_bnd_feproblem == NULL ? _invalid_boundary_id : _bnd_feproblem->getCurrentBoundaryID())

In FEProblem class declaration:

  const BoundaryID & getCurrentBoundaryID(){ return _current_boundary_id; }

As far as I can tell, the current_boundary_id doesn't really actually mean anything. Maybe it's just a way of monitoring whether things are happening at various bnd_id's during the program execution.

My question is whether as soon as setCurrentBoundaryID is run, does getCurrentBoundaryID() run for all the BoundaryRestrictable constructors? How would the computer know to do that? I'm so used to thinking of computers doing things sequentially, i.e. computers are dumb and they don't just guess that you need something to be done; you have to tell them that something needs to be done.

The answer is that Alex is dumb :-) All the _current_boundary_id's for all the BoundaryRestrictable objects created are references to the _current_boundary_id variable defined in the FEProblem class. Essentially, all the _current_boundary_id's belonging to the BC objects are slaves to the master FEProblem _current_boundary_id

MooseSharedNamespace:: = std:: or boost::
MooseSharedPointer:: = std::shared_ptr:: or boost::shared_ptr::

_factory.create(...) returns a MooseObjectPtr

from Factory.h:

typedef MooseSharedPtr<MooseObject> MooseObjectPtr;


Places where getDGKernelWarehouse.active() is called:

line 167 of ComputeFullJacobianThread -> Added analogous method to ComputeFullJacobianThread::computeInternalFaceJacobian. This means I need to track down all the original calls.
line 80 of ComputeJacobianThread
line 135 of ComputeJacobianThread
line 194 of ComputeJacobianThread
src/base/ComputeResidualThread.C:86 --> Check
src/base/ComputeResidualThread.C:198 --> This is the ComputeResidualThread::onInternalSide method for entire-domain DGKernels. I wrote an analogous method for restricted DGKernels. Because I wrote a new analogous method, that meant that I had to find where the original method was called and add the new method call. This is done. Check
src/base/ComputeResidualThread.C:211 --> This is in the same method as above. Check

Don't get carried away, doing more things than necessary!!!

My DGKernelWarehouse method names:

activeRestrictedDG(BoundaryID boundary_id, std::vector<DGKernel *> & active_restrictedDG) const
addRestrictedDG(const std::set<BoundaryID> & boundary_ids, MooseSharedPointer<DGKernel> & restricted_dg)

Unanswered question: since I'm putting my interface between subdomains, will the code I added in ComputeResidualThread::subdomainChanged() be effective?

ComputeFullJacobianThread is a child of ComputeJacobianThread. computeFaceJacobian is defined both in the parent ComputeJacobianThread and in the child ComputeFullJacobianThread. ComputeFullJacobianThread inherits onBoundary from ComputeJacobianThread. So if a ComputeFullJacThread object is created, a call to onBoundary will subsequently the ComputeFullJacThread computeFaceJacobian method.

From where is computeFaceJacobian called? It is only called in ComputeJacobianThread::onBoundary.

I'm not sure if this matters or not: in onBoundary in ResidThread, BoundaryID is set to invalid after the swapping of MaterialFace, whereas in JacThread, the boundary is set to invalid before the swapping of MaterialsFace. The fact that they're different between the two files suggests that the order of those stamements doesn't matter. But does it matter for dg where we have computeResidualNeighbor? Or some method like that.

Oh heavens, the only place where onInternalSide() and its analogous methods are actually called is in the base class ThreadedElementLoopBase. That means that there must be an onInterface() method defined in the ThreadedElementLoopBase class.

Public member functions of ComputeJacobianThread:

  	ComputeJacobianThread (FEProblem &fe_problem, NonlinearSystem &sys, SparseMatrix< Number > &jacobian)

 	ComputeJacobianThread (ComputeJacobianThread &x, Threads::split split)

virtual 	~ComputeJacobianThread ()

virtual void 	subdomainChanged ()
 	Called every time the current subdomain changes (i.e. More...

virtual void 	onElement (const Elem *elem)
 	Assembly of the element (not including surface assembly) More...

virtual void 	onBoundary (const Elem *elem, unsigned int side, BoundaryID bnd_id)
 	Called when doing boundary assembling. More...

virtual void 	onInternalSide (const Elem *elem, unsigned int side)
 	Called when doing internal edge assembling. More...

virtual void 	postElement (const Elem *)
 	Called after the element assembly is done (including surface assembling) More...

virtual void 	post ()
 	Called after the element range loop. More...

void 	join (const ComputeJacobianThread &)

virtual void 	caughtMooseException (MooseException &e)
 	Called if a MooseException is caught anywhere during the computation. More...

virtual bool 	keepGoing ()
 	Whether or not the loop should continue. More...

void 	operator() (const ConstElemRange &range, bool bypass_threading=false)

virtual void 	pre ()
 	Called before the element range loop. More...

Public member functions of ThreadedElementLoopBase:

 	ThreadedElementLoopBase (MooseMesh &mesh)

 	ThreadedElementLoopBase (ThreadedElementLoopBase &x, Threads::split split)

virtual 	~ThreadedElementLoopBase ()

void 	operator() (const RangeType &range, bool bypass_threading=false)

virtual void 	pre ()
 	Called before the element range loop. More...

virtual void 	post ()
 	Called after the element range loop. More...

virtual void 	onElement (const Elem *elem)
 	Assembly of the element (not including surface assembly) More...

virtual void 	postElement (const Elem *elem)
 	Called after the element assembly is done (including surface assembling) More...

virtual void 	onBoundary (const Elem *elem, unsigned int side, BoundaryID bnd_id)
 	Called when doing boundary assembling. More...

virtual void 	onInternalSide (const Elem *elem, unsigned int side)
 	Called when doing internal edge assembling. More...

virtual void 	subdomainChanged ()
 	Called every time the current subdomain changes (i.e. More...

virtual void 	caughtMooseException (MooseException &)
 	Called if a MooseException is caught anywhere during the computation. More...

virtual bool 	keepGoing ()
 	Whether or not the loop should continue. More...

number() returns _var_num from MooseVariableBase.h. _var_num is the variable number (from libMesh)

ElementElement or ElementNeighbor: i ==> master
NeighborElement or NeighborNeighbor: i ==> slave/neighbor

ElementElement or NeighborElement: j ==> master
ElementNeighbor or NeighborNeighbor: j ==> slave/neighbor

_sub_Kee, _sub_Ken, _sub_Kne, and _sub_Knn are all of type

  std::vector<std::vector<DenseMatrix<Number> > >

Reponse from JW Peterson about recompiling Moose code: "Your best bet is to settle on the interface early, thereby only changing Assembly.h once, and then as you develop the function in Assembly.C you will only have to recompile that one file."

Takeaway: limit modifications of .h files (e.g. decide on components of function declaration early so the .h file only has to be modified once), and then do all other function modifications (e.g. body modifications) in the .C files.

This is because the .h files are included all over the place, and wherever that .h file is included...that file will need to be recompiled, and then if that new file was a .h file, then any files that depend on the new .h file will have to be recompiled.

After successfully merging AlexDev into fork-devel (the local branch of the development branch of my Moose fork), zapdos-opt compiles as does moose-test-opt. However, I can't perform ./run_tests. Asked the list about this. Current state of things.

With 1e6 ballast resistor, low mobility and diffusivity ions, small plasma radius, no interpolation of transport parameters or the elastic rate coefficient, nl_abs_tol = 1e-4: 743 time steps, Active time = 656 seconds, output file size = 40 MB.

Material Property Checking occurs around line 325 of framework/src/base/SubProblem.C.

The "really cool" python script for comparing moose and sage was in the ~/gdrive/MooseOutput directory; well part of it was. All that was there was the matplotlib section. More importantly I recall that that work was done in the sage notebook after doing numpy and matplotlib imports.

New simulation with modular form (denoted by Coupling_plasma_liquid):

townsend_var_iz_const_el_old_ip_lg_plasma_radius_lg_b_resist_lg_abs_tol
# time steps = 911
Active time = 725 seconds (a little over 10 minutes)
Output file size = 49 MB

townsend_var_iz_const_el_old_ip_lg_plasma_radius_lg_b_resist_sm_abs_tol
# time steps = 1143

Solutions between the two look the exact same. Decreasing the tolerance did not speed the simulation up. It actually slowed the simulation down. The reason I tried speeding up this simulation is that 911 time steps seems like a lot. I thought that these might be very easy conditions that should solve quickly. However, it's possible that that is not the case. To illustrate, I just tried applying the modular form to an old test case:

Modular form
townsend_var_iz_const_el_old_ip_lg_plasma_radius_sm_b_resist
# time steps = 213

Old lumped form. All residuals lumped together for each variable.
townsend_var_iz_const_el_old_ip_lg_plasma_radius_sm_b_resist
# time steps = 439

Solutions to the two cases look the exact same. In other words it looks like I'm solving the exact same physics. In order words it looks like my residuals may be coded the exact same. That is great as long as I've coded the correct physics in both cases. However, since the # time steps is different, this suggests to me that the Jacobians are different, and that actually the jacobians in the modular form might be better. As a last test, I'm going to take this case:

townsend_var_iz_const_el_old_ip_lg_plasma_radius_lg_b_resist

and try running it with the old lumped form

The number of rows/columns in the Jacobian matrix is equal to: num_vars * num_nodes. In other words, it's equal to the number of test_functions i, times the number of variables j. Arp1, em1, pot1, en1, Arp2, em2, pot2, en2 ==> These are the values of each variable at the nodes for linear Lagrange and are also the coefficients of the basis functions, equal to the test functions for Continuous Galerkin.

Residuals: RArp1, Rem1, Rpot1, Ren1, RArp2, Rem2, Rpot2, Ren2

Can think of the residual RArp1 as being: the PDE (e.g. governing/conservation equaton) for Arp multipled by the test function 1, and then integrated over the domain

64 total elements obtained by differentiating each of the 8 residual terms, by each of the 8 coefficients.

# dofs = num_vars * num_nodes = num_sln_coefficients

List order is: pot1, em1, Arp1, en1, pot2, em2, Arp2, en2

Alright, using the finite difference jacobian just really embarrassed me. The problem solve in 265 time steps and the active time was 308.822 seconds. Wowzers, that's like five minutes. So where is the Jacobian wrong??

How did I get such awesome performance out of the FDP preconditioner when it's always been so slow in the past? My hypothesis is that in the past I always tried FDP with iterative solves where the barrier to convergence was the problem with the threading and really had nothing to do with the Jacobians.

Alright, so if I want to rapidly develop my physics (residuals), it looks like the best way to go is to just use FDP! But the fact that this problem solved perfectly, and matched my old lumped form's physics tells me that the 1) modular form residuals are the exact same as the old lumped form's residuals. This is also supported by the fact that I've been able to solve other problems today with the new modular form (in some cases better than the old lumped form, in some cases worse) and the final DC solutions have looked identical between the two cases. The second thing that this FDP exercise told me is that 2) whether the new modular form jacobians are better or worse than the old lumped form Jacobians, neither of them are perfect; e.g. the FDP performed way better than either of them.

Alright so the apparent plebian of the world is me because apparently I can't take derivatives to create jacobians. And the clear Jedi Master of the world is the FDP preconditioner combined with a direct LU solve. It doesn't even bloody matter that it's only one one processor. Power to the math and to efficiency! Down with the massive CPU clusters! Yep with the LFA formulation, the problem solved with the FDP preconditioner in 82.5272 seconds and only took 153 time steps. Jedi fucking master. With this unstoppable force, I'm tempted to try some other things: like Crank-Nicolson or higher order elements...however, I need to focus on getting results for now. Optimizing results can combe after obtaining results. For my thesis I need results...that's all I care about right now is results. I don't care whether I have beautiful analytic Jacobians which would then allow me to run on multiple cpu. For now I only care about results whether it's on one or many processors.

In the future I will return to my Jacobians. Remember that the world is computable!

Some things make sense a little more with respect to the jacobian debugging script. If the sum of all the jacobians is 0.0, then you will get the error message: a Jacobian needs to be implented.

Simulation failed at 2.07923e-06 seconds for no reason that I can see. All the solution variables are beautiful and smooth. The residual from the Newton solve is coming down beautifully. I really don't have much of a clue why MUMPS exited with a big old error message:

[0]PETSC ERROR: --------------------- Error Message --------------------------------------------------------------
[0]PETSC ERROR: Error in external library
[0]PETSC ERROR: Error reported by MUMPS in numerical factorization phase: INFO(1)=-9, INFO(2)=8

[0]PETSC ERROR: See http://www.mcs.anl.gov/petsc/documentation/faq.html for trouble shooting.
[0]PETSC ERROR: Petsc Release Version 3.6.0, Jun, 09, 2015
[0]PETSC ERROR: ./zapdos-opt on a arch-linux2-c-opt named lindsayad-OptiPlex-990 by lindsayad Mon Nov 16 15:32:49 2015
[0]PETSC ERROR: Configure options --prefix=/opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt --download-hypre=1 --with-ssl=0 --with-debugging=no --with-pic=1 --with-shared-libraries=1 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --download-fblaslapack=1 --download-metis=1 --download-parmetis=1 --download-superlu_dist=1 --download-scalapack=1 --download-mumps=1 CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 CFLAGS="-fPIC -fopenmp" CXXFLAGS="-fPIC -fopenmp" FFLAGS="-fPIC -fopenmp" FCFLAGS="-fPIC -fopenmp" F90FLAGS="-fPIC -fopenmp" F77FLAGS="-fPIC -fopenmp" PETSC_DIR=/home/lindsayad/projects/src/petsc-3.6.0
[0]PETSC ERROR: #1 MatFactorNumeric_MUMPS() line 1172 in /home/lindsayad/projects/src/petsc-3.6.0/src/mat/impls/aij/mpi/mumps/mumps.c
[0]PETSC ERROR: #2 MatLUFactorNumeric() line 2948 in /home/lindsayad/projects/src/petsc-3.6.0/src/mat/interface/matrix.c
[0]PETSC ERROR: #3 PCSetUp_LU() line 152 in /home/lindsayad/projects/src/petsc-3.6.0/src/ksp/pc/impls/factor/lu/lu.c
[0]PETSC ERROR: #4 PCSetUp() line 982 in /home/lindsayad/projects/src/petsc-3.6.0/src/ksp/pc/interface/precon.c
[0]PETSC ERROR: #5 KSPSetUp() line 332 in /home/lindsayad/projects/src/petsc-3.6.0/src/ksp/ksp/interface/itfunc.c
[0]PETSC ERROR: #6 KSPSolve() line 546 in /home/lindsayad/projects/src/petsc-3.6.0/src/ksp/ksp/interface/itfunc.c
[0]PETSC ERROR: #7 SNESSolve_NEWTONLS() line 233 in /home/lindsayad/projects/src/petsc-3.6.0/src/snes/impls/ls/ls.c
[0]PETSC ERROR: #8 SNESSolve() line 3894 in /home/lindsayad/projects/src/petsc-3.6.0/src/snes/interface/snes.c

Alright, I'm pretty confident that the reason the solve is failing is because of the mean_energy and specifically because of the mean_energy interface (right boundary) condition. I'm essentially abandoning the energy boundary condition while at the same not abandoning the electron boundary condition. That just doesn't seem to make much sense. Going to LFA.

Solve converges great with LFA with electron mobility and diffusivity the same between the two phases. Obviously that's not physically realistic.

Going back to more physically realistic transport coefficients in the liquid phase. Solve fails at 1.54524e-07 seconds and after 457 time steps. Solve fails because of diverged line search. I'm using an FDP preconditioner, so the Jacobian should be perfect. When the solve fails, here are the last set of variable residuals:

 4 Nonlinear |R| = 3.728681e+02
    |residual|_2 of individual variables:
                  potential: 7.1762e-07
                  em:        2.91941e-12
                  Arp:       4.63375e-13
                  OHm:       372.697
                  H3Op:      0.000331611

 5 Nonlinear |R| = 3.726967e+02

I think the solve could be failing for a couple of reasons: 1) the large disparity in variable residuals 2) There are some oscillations in the OHm and H3Op densities.

I'm going to try adding some stabilization to OHm and H3Op and see if that changes the solve fail time (for the better).

Let's assume that I don't have a _block_diagonal_matrix. But let's also know that I have a function for testing this now. It's an assembly method called Assembly::isBlockDiagonal. Most of the time we're calling Assembly::jacobianBlockNeighbor. To help, here is the call:

 DenseMatrix<Number> & Kxx = type == Moose::ElementElement ? _assembly.jacobianBlock(_var.number(), jvar) :
                              type == Moose::ElementNeighbor ? _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, _var.number(), jvar) :
                              type == Moose::NeighborElement ? _assembly.jacobianBlockNeighbor(Moose::NeighborElement, _neighbor_var.number(), jvar) :
                              _assembly.jacobianBlockNeighbor(Moose::NeighborNeighbor, _neighbor_var.number(), jvar);

What gets returned to the DGInterface is a DenseMatrix, something like this:

    case Moose::ElementNeighbor: return _sub_Ken[ivar][jvar];

Now what is _sub_Ken? It is defined to be:

  std::vector<std::vector<DenseMatrix<Number> > > _sub_Ken;

So it's a vector of a vector of DenseMatrix. Now how is the DenseMatrix formed? This must be something that's done in libMesh because there's no modification of _sub_Ken that occurs in Moose other than:

  _sub_Ken.resize(n_vars);

and:

      _sub_Ken[i].resize(n_vars);

What this means if the matrix is not block_diagonal: _sub_Ken is a vector of length n_vars of a vector of length n_vars of matrices of unknown sizes.

MooseMesh::elem returns a pointer to an Elem

MooseVariable::getDofIndices takes a const Elem * elem

Well that should be good because _current_elem is of type const Elem * !!!

Here's the full interface for MooseVariable::getDofIndices:

MooseVariable::getDofIndices(const Elem * elem, std::vector<dof_id_type> & dof_indices)

It looks like the .unv mesh I create with salome is not read perfectly into Moose because the values are not matched in a smooth way at the interface. I don't know quite what I'm going to do about that yet.

From DGMatDiffusionInt:

_D is a reference to memory location: @0xaaf280
_D data is at 0xb61e20

_D_neighbor is a reference to memory location: @0xae69b0
_D_neighbor data is at 0xb62420

From Water:

_diffpotentialliq is a reference to memory location: @0xace440
_diffpotentialliq data is at 0xb5ebf0

From ArgonConstTD:

_diffpotential is at: @0xa98900
_diffpotential data is at 0xb204b0

From CoeffDiffusionLin:

When _var.name() == "potential"
_diffusivity: @0xa98900
_diffusivity data: 0xb204b0

When _var.name() == "potentialliq"
_diffusivty: @0xace440
_diffusivity data: 0xb53bf0

Alright we have an exact match between the kernel object memory locations of the material property and the material object memory locations. However, the dgkernel object doesn't match either!!!

Well maybe there's some comfort in the fact that the memory locations of _D and _D_neighbor in the dgkernel aren't changing over time. So we need to figure out where the memory location gets set and how to make sure it's set to the right place, e.g. the correct material objects.

Currently at the interface (x=.001), em has a value of 11.8117 and emliq has a value of 31.348. That's not a very tied value!!! Of course if the matched value constraint isn't on then I supposed it makes sense that the values aren't tied. Dummy!

I think I've determined that the appearance that there is not a node at .001 is an artifact of paraview and not of the mesh that salome outputs. I would like to find a way to directly examine exodus files and see if there's a way that I can read node coordinates within the exodus file.

Although now that I think the artifact is from paraview, I can adroitly choose the resolution in Paraview to reveal the node at .001.

Something to keep in mind is that what we are trying to do is minimize a residual. NO conditions are strongly imposed! Everything is weak. Constraints are imposed weakly, e.g. as a component of the residual. I think it's beautiful how everything is tied together in one unified way. One global residual to rule them all!

But, this also leads to the fact that boundary/interface conditions will not be perfectly imposed. E.g. I'm looking at the interface between em and emliq and I see that em has a value of 34.3336 at the interface whereas emliq has a value of 34.3284. So not a perfect match.

Alright awesome, I can use paraview to write the exodus mesh data (after Moose converts the mesh from .unv to .e) to a .csv file and then examine node coordinates there. Interestingly, Paraview writes one file for each element block. I can look inside and beautifully see that the correct nodes are present. Schniceeeeee!!!!!

When paraview writes out a plot_over_line to a csv file it writes out a number of points equal to the setting of the resolution. This means that if the resolution is high enough, more points will be written than mesh nodes. I would like to know if paraview can just plot the points from the mesh. That would be very nice.

Actually I can look in the csv data file and see that em and emliq were perfectly matched at the interface! Both had values of 34.3336 at nodeID=401

Potential and potentialliq were also perfectly tied at the interface node with a value of -12.3762.

EField_gas and EField_liq are elemental variables, so when I write fields assocated with points, those two fields aren't written. Elemental vs. nodal!!

Alright, I've discovered that Plot Data is a better choice for my needs. And the plot shows expertly that I have some serious oscillations in the liquid electron concentration, almost certainly entirely due to the steep gradients there. I have a couple of choices then: 1) refine the mesh or 2) add some stabilization. This means that I need to return to the question of whether there's a minimum element size in salome.

With mesh refinement, beautifully able to eliminate oscillations :-)

Remember 4e11 for the TotalFlux of emliq at steady-state.

Now I'm looking at 6.3e12 for the TotalFlux of emliq at steady-state. Keep in mind, however, that the TotalFlux for emliq had a spatial, temporal maximum of 4.1e14, which is almost 100 times higher. We could just be a slave to the absolute tolerance here.

With the logarithmic stablization removed, I see a lot more of what I want to see: electrons flow out of the right boundary such that the electron density is reduced throughout the entire simulation until the concentration of electrons in the liquid is reduced below the absolute tolerance. What we're seeing here is a much more accurate reflection of reality. There are currently no physical sources in the system if all we have is charges and we're sweeping them with a potential. Consequently, there should never be a steady-state because there are no sources that can balance sinks!!! Thus I should be able to remove the absolute tolerance and the problem should solve all the way to 0.1 seconds without ever doing 0 non-linear iterations to solve the problem. Now it's possible that I might run into some machine-precision problems, but that's alright. Now the only thing I'm concerned about is enforcement of my boundary-conditions. I'm definitely not seeing a zero gradient in the electron liquid concentration at the right boundary. Could that be because the mesh is not fine enough? Or because my residual is wrong? I could also try scaling my boundary residual. That might be something interesting to try.

The boundary condition is definitely doing it's job in allowing electrons to flow out, that's for definites.

Absolutely beautiful: with a refined mesh, my zero gradient boundary condition is enforced beautifully :-) :-) :-) Rejoice when you encounter victories :-)

You can tell that the logarithmic stabilization doesn't start coming into effect until log_concentrations of about 19.

Clear oscillations in H3Op, OHm, and to a lesser extent Arp without advection stabilization. However, after adding the advection stabilization kernels, the oscillations disappeared :-) :-) Success!

Tomorrow I'm going to:

1) Look into Moose's finite difference preconditioner (nothing to learn about in there)
2) Look into the condition number of the Jacobian (follow PetSc's FAQ page)

Important to note that if I don't ramp up the boundary condition for the potential in a transient fashion, I can get the simulation to solve at least one time step (definitely more in this case) if I have these species: Arp, em, emliq, OHm, H3Op with certain values for scaling (see this input file: ~/zapdos/problems/Coupling_plasma_liquid_LFA_dg_with_ionization_sim_out_to_5.41e-07_seconds.i) and with the initial condition for the potential set to 0 in the whole domain. However, if I remove OHm and H3Op and keep the same scaling for the remaining variables, then the problem won't solve a single time step. My hypothesis is that in the former case the solver can reduce the residuals of H3Op and OHm sufficiently such that the global residual is reduced below the nl_rel_tol. However, with those variables removed it cannot. The potential residual is now too large to overcome. So, if I replace the constant 0 initial conditions for the potential with an initial condition consistent with the boundary conditions then it will solve at least one time step.

It's good to know that I can checkout my zapdos master branch and gold files still execute as gold. Thank you git!!! Restoring some confidence in myself :-)

Some very good notes on matrix condition number and some mentions of singular values:

If the condition number is not too much larger than one (but it can still be a multiple of one), the matrix is well conditioned which means its inverse can be computed with good accuracy. If the condition number is very large, then the matrix is said to be ill-conditioned. Practically, such a matrix is almost singular, and the computation of its inverse, or solution of a linear system of equations is prone to large numerical errors. A matrix that is not invertible has condition number equal to infinity.

I have a sample size of one so don't put a lot of stock in this statement, but it looks like the condition number of the A matrix (in my case this is the Jacobian matrix, formed by finite differencing the residuals) decreases with the # of non-linear iterations as long as the residual is decreasing.

It looks like -ksp_diag_scale_fix could be something very interesting and potentially something very useful. Applying those petsc options let to a reduction in my condition number of 8 orders of magnitude!!! And also it reported 0 of 851 singular values are (nearly) zero; before the diag_fix it was reporting 1 of 851. Unfortunately, the solve wasn't converging but maybe that's because I was using -pc_type=svd ? Does that actually apply a preconditioner or no?

Sure enough, the FDP preconditioner appears to reduce the nonlinear residual more than the SMP preconditioner (this also confirms that the preconditioner block is having an effect).

When using solve type = NEWTON, if the Jacobian is bad, the problem primarily shows up in the non-linear iterations. When solve type = PJFNK, then if the Jacobian is bad, the problem primarily shows in the linear iterations. At least that's my understanding.

SVD stands for singular value decomposition. It looks like svd is a real genuine preconditioner. However, the cases when it is useful and effective are unknown to me.

The number of "singular values" listed by the SVD process correspondson to the number of DOF in the nonlinear system.

Although both effective, the lu preconditioner appears to be more effective for my problem than the svd preconditioner (both quickly reduce the linear residual, however). However, both solves have their problems with the line search.

Changing the differencing parameter to ds had no effect on the condition number. Nor did adding a non-default coloring error.

Well this is fascinating. If I try and run my_coupled_tied_value.i with solve_type=NEWTON, then the problem doesn't solve. Things to note: my condition number for this test case is 1e17. Back to a huge number again! And might not this huge number explain why the NEWTON solve is failing?

This is going to be a very useful diagnostic test case I think. When I added an FDP, then the problem converged but it took 38 nonlinear iterations with pc_type=svd. With pc_type = lu, it took 34 non-linear iterations. I'm continuing to hypothesize that the slow convergence is because of the condition number. Interestingly, FDP yields norm of matrix ratios of 0.0356613. What's going on there? Without providing the preconditiong matrix, the norm of matrix ratios yields 0.0356617, which is almost the exact same.

Well this is perplexing. If I remove the DGKernel from my_coupled_tied_value.i, the condition number is still 1e17, but the problem solves in one nonlinear iteration. Remove matched value: still the same condition magnitude. Remove w, however, and the condition number drops to 1e4. Wowzers, that's a huge drop! Problem solves in one non-linear iteration. Add back DGKernel and MatchedValueBC: condition number = 6e3. But the problem's non-linear solve doesn't converged! E.g. we get a bad line search!!! Take away the dgkernel again: solves in one non-linear iteration; condition number = 2e4. What does this suggest? The condition number is not responsible for the bad line search. In fact it points to the DGKernel as being responsible for the bad line search.

Changing from FDP to SMP preconditioner didn't help the line search. Changing to line_search=none didn't make the solve work any better.

With PJFNK, the solve completes in 1 non-linear iteration (2 linear iterations) with DG and with SMP preconditioner. Same with FDP preconditioner. So I guess that PJFNK can handle DGKernels but Newton cannot?? Why?

PJFNK with pc_type=lu and ksp_type=preonly yields the same failed line_searches as NEWTON (I mean it basically is NEWTON right?) pc_type=lu, default ksp_type: 2 linear its (1.305e-11), 1 non-linear it (8.829e-9). Default pc_type, default ksp_type: 2 linear its (8.32e-12), 1 non-linear it (8.829e-9).

Alright, so apparently nothing happens if I use the -snes_linesearch_monitor flag during a PJFNK solve. Well how am I supposed to look at the line search if it can't be displayed? I'm being told that the snes solve did not converge because of DIVERGED_LINE_SEARCH. Wait if I run my_coupled_tied_value.i, then the line search gets printed. Wtf!?

If I specify line_search = none in my_coupled_tied_value.i with the dgkernels and solve_type=PJFNK, then it solves perfectly fine (and doesn't output a line search).

If I run PJFNK on my_coupled_tied_value  with -ksp_diagonal_scale and -ksp_diagonal_scale_fix, then the solve doesn't converge. So many options, and so much I don't understand!!!! What I'm coming to understand is how powerful math is, and how powerful having mathematical knowledge is.

Alright so trying a PJFNK solve on my Coupled plasma_liquid system, results in bad line searches, as well as non-monotonic changes in the true resid norm. What else could go wrong??

If I modify the DGKernel by getting rid of the 0.5 factors, then my test case, my_coupled_tied_value won't even perform a single linear iteration. It exits with the error: Linear solve did not converge due to DIVERGED_NANORINF iterations 0

For some reason, it just seems like NEWTON can't solve any of my problems with DGInterface kernels. It just can't!! Every time I tried to solve with NEWTON and I actually got a time step, the potential on the liquid side of the interface was positive!!! That is absoutely physically incorrect and is destined to never reach steady-state! Well I ended the day learning something important: for whatever reason I can't use NEWTON with my DG system. The moment I moved over to PJFNK, it got the potential sign on the interface side correct. Someday I would like to understand why this is the case.

Turning off the NeumannCircuitVoltage boundary condition and replacing it with a Dirichlet BC resulted in the problem converging perhaps to its steady state (no ionization, just motion and stabilization). Why?

With the dirichlet BC, my initial condition number appears to be around 2e29 with 105 of 1303 singular values (nearly) zero. That's not very good! Replacing the ConstantIC with a consistent FunctionIC didn't appear to make the condition number any better: 4e29 with 105 of 1303 singular values. Gross. With NeumannCircuitBC and constant ICs, the condition number is 5e29 with 105 of 1303 singular values. All terrible! Yet the dirichlet problem converges while the NeumannCircuitProblem does not!

Ok, so I found that changing the value of my boundary condition in simple.i from 1 to 1000 did not change the condition number at all. However, scaling the variable residual did change the condition number.

scaling = 1e-3  condition number = 5.10e2
no scaling 	condition number = 6.22e1
scaling = 1e3	condition number = 6.08e4
scaling = 1e1 	condition number = 6.085e2
scaling = 1e-1	condition number = 2.05e1

Ok so it appears that it is possible to improve the condition number with scaling, but it may not be quite what you think. Changing the IC also does not change the condition number.

When I can run a simulation with OHm and H3Op out a decent ways in time, and then all I do is eliminate OHm and H3Op and the solve no longer converges: this suggests to me that I have a scaling problem or that there is more that I need to learn here.

When I have the full variable set, the condition number is 1.56e29. When I remove H3Op and OHm, the condition number is 2.27e30. Doesn't seem to huge of a difference but what do I know?

On the Wikipedia page, we have this, which I know I've read before:

"Instead of solving the original linear system above, one may solve either the right preconditioned system (for my case, A is the Jacobian matrix formed by taking the derivatives of R with respect to u, and b is -R for the ith iterate, and P is the preconditioning matrix that we can choose):

    AP^{-1}Px = b

via solving

    AP^{-1}y=b

for y and

    Px=y

for x; or the left preconditioned system:

    P^{-1}(Ax-b)=0

both of which give the same solution as the original system so long as the preconditioner matrix P is nonsingular. The left preconditioning is more common.The goal of this preconditioned system is to reduce the condition number of the left or right preconditioned system matrix P^{-1}A or AP^{-1}, respectively. *The preconditioned matrix P^{-1}A [Left preconditioning] or AP^{-1} [Right preconditioning] is almost never explicitly formed. Only the action of applying the preconditioner solve operation P^{-1} to a given vector need to be computed in iterative methods.*

Typically there is a trade-off in the choice of P. Since the operator P^{-1} must be applied at each step of the iterative linear solver, it should have a small cost (computing time) of applying the P^{-1} operation. The cheapest preconditioner would therefore be P=I since then P^{-1}=I. Clearly, this results in the original linear system and the preconditioner does nothing. At the other extreme, the choice P=A gives P^{-1}A = AP^{-1} = I, which has optimal condition number of 1, requiring a single iteration for convergence; however in this case P^{-1}=A^{-1}, and applying the preconditioner is as difficult as solving the original system. One therefore chooses P as somewhere between these two extremes, in an attempt to achieve a minimal number of linear iterations while keeping the operator P^{-1} as simple as possible. Some examples of typical preconditioning approaches are detailed below."

I need to figure out what matrix free really means in practice for iterative systems.

It's clear that my SMP preconditioning matrix must suck because I can't get any convergence even with my "full" system.

I can get the reduced set to converge if I remove the electron dg advection and diffusion kernels. If I add back the diffusion dgkernel, no convergence.

Ok, I could get Newton to do some stuff if I turned off all the dgkernels and then properly scaled the potential.

Instant non-convergence when the potential dgkernel was added back :-(

Ok, by some miracle, I switched the the primary variable that the potential MatchedValueBC was acting on...and now Newton's method is convering like a bloody fucking champion!!! What the hell just happened?!?!?!?!

My current condition number is 2.73e16. That seems smaller than some I've had in the past. But if I switch my MatchedValueBC, which seems to be the key to convergence, then my condition number is about the same but my convergence is lost.

If I let my nl_rel_tol get too loose, then I see errors in my Arp solution, and the simulation fails (time step goes to dt min). I want to disect and understand reasons for convergence, slow convergence, lack of convergence.

So, I was having some problems...I kept noting the emergence of spikes in the Arp profile towards the right boundary. I speculated a couple of things: 1) bad residual scaling 2) bad boundary condition 3) loss of accuracy perhaps due to finite differencing

Scaling up the Arp residual didn't remove the spikes. Removing the ion boundary condition didn't remove the spikes. But changing from an FDP to SMP removed the spikes and the simulation ran out to a world record of 1.69e-5 seconds!!! There's so many things to think about... 1) residual scaling 2) boundary conditions 3) conditioning 4) quality of Jacobian 5) etc.....

Delving into solution context here. This is from Quarteroni's Nonlinear Systems and Numerical Optimization text:

"Theorem 7.1 thus confirms that Newtons method is quadratically conver-
gent only if x(0) is sufficiently close to the solution x and if the Jacobian
matrix is nonsingular. Moreover, it is worth noting that the computational
effort needed to solve the linear system (7.4) can be excessively high as n
gets large. Also, the Jacobian could be ill-conditioned, which makes it quite diffi-
cult to obtain an accurate solution. For these reasons, several modifications
to Newtons method have been proposed, which will be briefly considered
in the later sections, referring to the specialized literature for further details
(see [OR70], [DS83], [Erh97], [BS90] and the references therein)."

Root-finding and unconstrained minimization (optimization) are intricately linked. I love learning this stuff!

Alright, it sounds like the default bt line search method DOES NOT ALLOW increases in the residual. However, cp and l2 do. line_search=none -> no line search -> take full Newton step in direction indicated by the linear solve of J*delta_x = -R. cp tends to work well for problems derived from minimization. In our case, the physical problem is a system of PDEs, typically derived from conservation laws. Thus in reality they fall into category 1 of these two closely related categories 1) find roots of system of equations 2) unconstrained minimization. However, with Moose we convert from problem 1 to problem 2, attempting to minimize a residual in an unconstrained manner.

I can easily see now how a failed line search can be an indicator of a bad Jacobian. If we're using a bt line_search that doesn't allow the new residual iterate to increase from the old residual iterate, then we need to have solvev for a good delta_x in our linear solve. If our Jacobian is bad, then we might not get a good delta_x, and then we might not be able to take a step that decreases the residual. The Jacobian is so central and key!!!

Now there seems like a fair number of ways that we might get a bad Jacobian. If using an analytic Jacobian, then wrong analytic expressions or incomplete analytic expressions, would obviously lead us astray. However, what if we use a finite differenced Jacobian? What if our residuals or our variable values are all over the place? Then we might introduced some error into our Jacobian. Jed said that finite differencing leaves us with about 8 significant figures. I'm not quite sure what that means.

My variable values, my DOFs, should actually all be fairly small because I'm using a logarithmic formulation for all my species concentrations. However, my residuals as I've witnessed can be many orders of magnitude larger (we have exponentials and all that jazz).

There are many different ways to compute a matrix norm.

my_coupled_tied_value trials:

No scaling:
SMP Jacobian, MatchedValue on u: Converged, 19 its, 1.17e4 condition number
FDP Jacobian, MatchedValue on u: Not converged, 6.66e3 condition number
FDP Jacobian, MatchedValue on v: Not converged, 6.14e3 condition number
SMP Jacobian, MatchedValue on v: Not converged, 1.89e4 condition number
Block-diag Jacobian, MatchedValue on v: Not converged, 1.53e4 condition number
Block-diag Jacobian, MatchedValue on u: Converged, 36 its, 9.45e3 condition number

u scaling = 0.1:
SMP Jacobian, MatchedValue on v: Converged, 2 its, 1.61e5 condition number

So the condition number got worse and the Jacobian should have been no more correct, and the solve converged! There must be more reasons for why a solve does/doesn't converge than correctness of Jacobian and conditioning of Jacobian. The scaling of the residuals must play a role. Either that or SVD doesn't accurately compute the condition number. I want to understand all the fundamental mathematical reasons for why things may not converge.

Idea for testing analytic Jacobian: with different input files, test number of time steps with SMP vs FDP. If SMP < FDP, then that doesn't mean that the analytic Jacobian is exactly correct; however, if FDP < SMP, then I know that the anaytic Jacobian is not exactly correct.

For nxn matrices of real numbers, all norms are equivalent. Hmm, that's very interesting.

It's clear that different residual scaling leads to different solutions. Just compare 1.49e-5_fail_time_everything_smooth_emliq_resid_scaled_down to 1.67e-5_fail_time_everything_smooth_out.e. At the same point in time the scaled up emliq residual case has a peak density of 4.6e24 whereas in the other case the peak density is at 1.3e23. That's over an order of magnitude difference.

The left boundary potential doesn't start changing until about 1.4e-5 seconds.

From the beautiful book of Dennis: All our algorithms for nonlinear equations are implemented using a diagonal scaling matrix D_F on F(x) discussed in section 7.2. It causes:

f(x) = 1/2 * (L2Norm(D_F * F(x)))^2; grad_f(x) = transpose(J(x)) * D_F^2 * F(x); H = transpose(J(x)) * D_F^2 * J(x)

We use global algorithms for unconstrained minimization in solving systems of nonlinear equations (e.g. line search).

I think that residual scaling must have a role in conditioning of the jacobian matrix.

Ok, in the successful solve, the condition number is 360 (based off l1norms).

Alright it looks like when PetSc outputs the condition number using -pc_svd_monitor that it outputs the condition number based on l2norms. And it outputs it correctly! So indeed, I am getting a better solve with a worse conditioned Jacobian! So I can rule out Jacobian conditioning. Moreover, with no concern of precision loss, the Newton step should be the same regardless of the scaling (i.e. the problem is not with the linear solve). Thus the problem has to be with the line search. Interestingly if I use a cp line-search, then I can get convergence with the default scaling of 1 and 1. And I get worse convergence efficiency if I scale u. That seems more consistent with reality.

Alright so after I do the mesh conversion, there are 1008 nodes and 1007 elements in the .e file.

If I go into, the .msh file, it too shows 1008 nodes and 1007 elements. When gmsh is outputting info, it says 1008 vertices and 1023 elements. So it obviously has the number of nodes correct, but where is this 1023 number coming from?

As stated in Dennis's text: the residual scaling did not affect the Newton step. I verified that -ksp_view_solution does in fact show the solution for delta_x and not x. However, it actually appears to show -delta_x, as in x_(k+1) = x_k - delta_x. However, the difference is that in my scaled test case, the line search takes the full Newton step (i.e lambda = 1). In the unscaled case, however, lambda is set equal to 0.1. Why? Why does the line search fail here? It's got to be something to do with boundary vs. initial conditions right?

Side note: gmsh starts counting its nodes at 1. GeneratedMesh does the same.

It should be noted that I can get an options_table at the END of my simulation by specifying the petsc option: -options_left.

Tasks while waiting to figure out how to get my mesh into exodus II format:

1) Continue to improve understanding about how to evaluate the accuracy of the Jacobian matrix
2) Understand why a line search might fail.

A tolerance of 1e-3 is definitely too loose. The solution is inaccurate using nl_rel_tol=1e-3.

Moving the problem to a second order formulation actually made the solution less accurate for a given nl_rel_tol of 1e-4!! I'm going to look at reducing the rel_tol. That didn't work. Abandoning the second order formulation!!

From Wolfram MathWorld:

Condition number: the ratio C of the largest to smallest singular value in the singular value decomposition of a matrix. The base-b logarithm of C is an estimate of how many base-b digits are lost in solving a linear system with that matrix. In other words, it estimates worst-case loss of precision. A system is said to be singular if the condition number is infinite, and ill-conditioned if it is too large, where "too large" means roughly log(C)>~ the precision of matrix entries.

The biggest residual by far at my brick wall here is the em residual. I'm wondering whether I'm running into the situation where my largest residual contribution is coming from the interfacial flux of eletrons and I don't have off diagonal jacobians for the interface advection for example. I'm kind of shooting in the dark here, but that's a guess. To reduce the residual, you need a good direction, meaning you need a good Newton step. To get a good Newton step, you must have a good jacobian!!!

Sure enough, adding the binary reaction term for the electrons added a lot more OHm, so it's clear that the ProductAABBRxn kernel is having an effect.

It looks like the concentration profile of emliq didn't change at all, however, so if the ReactantAARxn kernel is having an effect, the sink term must be being balanced by an increased flux of electrons from the gas phase. I could test that if I had the flux auxkernels, but I didn't have them active. Yea, I think I need to look into this OffDiagJacobian business for shizzle.

Criticial to remember what error means:

Rleative error in y as an approximation to a non-zero x:

|x-y| / |x|

computeOffDiagJacobian for DGInterface has compiled. The thing to do tomorrow is see whether a test file can be run without generating any segmentation faults and then whether I can examine the Jacobian and see that it is what it should be.

Recompilation of a source file is done if the source file time stamp is more recent than the object file time stamp or if the object file does not exist. There you go! That should tell you whether files will have to be recompiled when using git or rsync between different machines! Before it was a mystery. Simple time stamping is how it's done! What do you think make is actually going to go into the file and investigate whether it's really changed to make a decision about whether to recompile? No that's just silly. Looking at the time stamps makes sense.

git stash reverts a file modified in the current working directory to it's previous state at HEAD. That changes the file time-stamp to whenever you press git stash. git stash apply also changes the file stamp, corresponding to whenever you entered the command. So pretty much any git command that modifies the state of a file will fast forward the time stamp of the modified file to the time that the git command is entered.

I think what I've been doing with my DGInterface kernels is wrong. Maybe.

It's important to think here! If jvar == _var.number() then OffDiag gets passed with NeighborElement. If jvar == _neighbor_var.number() then OffDiag gets passed with ElementNeighbor.

Ok the residual from a kernel is equal to: _JxW[_qp] * _coord[_qp] * computeQpResidual();

JxW is the current quadrature point weight value.

coord is the scaling factor to convert from cartesian to another coordiante system (e.g. rz, spherical, etc.)

CoeffParamDiffusion:

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _JxW[_qp] * _coord[_qp] * _D *  _grad_phi[_j][_qp] * _grad_test[_i][_qp];

In first order Lagrange (at least in 1-D), there are two non-zero shape functions within each element.

Alright, I know that n_point (number of quadrature points) = 2 (qp = 0,1)
phi.size() = 2 ==> j = 0,1
test.size() = 2 ==> i = 0,1

Here's a key: every variable must either be a degree of freedom, constant, or boil down to a function of the degrees of freedom. Those dof's are our unknowns!

I believe that block really is the apt term to describe the residuals and jacobians. I believe that the dimensions of these blocks are defined by the number of variables, NOT the number of DOF's (e.g. num_vars * num_nodes (in the case of linear Lagrange)).

This is confirmed by grepping for ResidualBlockNeighbor. The only place that function is called is in DGKernel, my beloved DGInterface, and NodeFaceConstraint. Thus it is only called in the places where there is a chance of having residual contributions from two different variables.

By looking at snes_test_display, it is quite clear that my hand-coded Jacobian for DGDiffusionInterface is wrong. Now I have to figure out why. To understand why, I should understand how the Jacobian is calculated and currently I don't have that!!!

Looking at the jacobian, it's like the Diffusion kernel doesn't even frieking exist!?!? What's going on here? That's not the case. Rather the Dirichlet BC looks weird!! How is the u1 residual a function of u2??

This makes some sense now I think: dof rows whose dof has a nodal boundary condition applied to it get zeroed. And then the diagonal point of the dof row is set to 1.

Alright I've got if figured out now. Wow this software is smart! It knows not to zero the u block when applying a nodal condition to v!

This also seems to indicate that the NodalBC's are being strongly enforced. This would make sense because I've never seen a NodalBC be in the least bit violated.

Zeroing can be confusing I think. Let's try eliminating the zeroing effects. Yes something's definitely wrong with my DGJacobian. Let's think about what it is.

The test of whether to compute the off-diagonal Jacobian should be:

  for (unsigned int i=0; i<boundary_ids.size(); ++i) {
    boundary_info.sideset_name(boundary_ids[i]) = boundary_names[i];
    boundary_info.nodeset_name(boundary_ids[i]) = boundary_names[i];
  }

With increased artifical diffusion, the solve actually fails sooner.

So WR is 5.70e-5 seconds. Without any artificial diff, solve failed at 5.65e-5 seconds.

With pot scaling at 1e4, solve failed at 9.06e-6 seconds, so worse.

With nl_rel_tol at 1e-8, solve failed at 4.25e-6 seconds, so worse.

With pot scaling at 1e0 and no artificial diff, solve failed at 5.16e-5 seconds, so worse.

Time Step 167: 4.96e-9
Time Step 166: 4.74e-9
Time Step 165: 4.98e-9
Time Step 165: 4.80e-9
Time Step 164: 3.85e-9

Very important note: at least in regular DGKernels, _normals[_qp] is the same whether your type is Element or Neighbor, i.e. the normal vector does not change when you switch over to the other side!!!

Alright, it looks like here in the case where I have one variable that's linear monomial, and another that's linear Lagrange, then all the variable dof's are grouped together in the jacobian. Linear monomial: dofs are associated with element, two dofs per element. Linear lagrange: dof associated with the node, one dof per node.

Problem is with the off-diagonal jacobian: derivatives wrt node 1 should be wrt node 0 (e.g. potential at node 0, not node 1)

Alright, at _qp=0 and _i=0, within the computeQpResidual, _grad_potential = -.33, _grad_potential_neighbor = .91, _velocity_average = -.32

Should only be one dgkernel executed right? Only one internal face?

With first monomial, there are still two separate functions in an element, the const a, and then b * x. Thus test.size() should still be 2.

Dopey Alex. The new computeQpResidual calls come from moving from the FD problem to the analytical problem. You're a dope face.

The u and v shape functions are not the same!!!

Alright, Clm, OHm, and Nap now have enormoush residuals after the switch to first monomial and DG formulation. E.g. their residuals are on the order 1e10, 1e11, whereas the other variables are on order 1e-10. Clearly something is wrong. Perhaps it's from sigma. That would be my guess actually. I'm just trying to think of thinks with divisors. After setting sigma = 0, the residuals are still enormous, so it's something else.

I believe the reason that it takes so long between the display of the initial residual, and any other output (e.g. beginning of nonlinear and linear iterations), is that we're creating the memory for the Jacobian. That's my belief. It makes sense that it is has to be an initialization process associated with the jacobian, because we've already seen the residual, and after the initial wait, there are no more similar waits after.

The enormous residuals are eliminated by removing DGCoeffDiffusion, so that appears to be the problem. So the problem is with _grad_test and _grad_test_neighbor. I believe this is becuause the problem is trying to calculate DGKernels at a non-subdomain-internal side. E.g. it's trying to do the calculation at a boundary, where the variable doesn't exist on the other side. This must lead to discontinuties and the huge test function gradients.

Alright, it's clear that for first order monomial basis functions, in order to achieve the correct solution, DGDiffusion kernels MUST be combined with Diffusion kernels. Cannot go without them.

Wrong!!!

Alright, with enough mesh points, numerical solution can match analytic solution exactly, using CG, DG, integrated BCs, and monomial basis functions.

Alright got oscillations with Crank-nic, but not with fully implicit. Now let's see if upwinding helps out the cause.

If we know nothing about Dirichlet boundaries, then the epsilon and sigma terms disappear from the boundary condition.

So I know exactly how the BC should look with just diffusion. Now how about when we add advection?

The Inflow and Outflow boundary conditions appear to be doing absolutely nothing! What the hell is going on here? E.g. when the inflow and outflow boundary conditions are imposed, the residuals are totally unchanged from when they're not imposed. Are we even contributing a residual from there? What's going on?

With CG+DG, oscillations don't grow in time, and the oscillations are local to the moving front.

These tests show that in order to take advantage of DG, you must use discontinuous shape functions! Using Lagrange when using DG+CG will essentially just make the method CG!

Yidong Xia and Austin Ladshaw are two folks very interested in advection dominated problems. Austin Ladshaw is a major contributer to the dgosprey application that I've downloaded.

_block_internal_side_user_objects is referenced by addUserObject, intneralSideUserObjects, and updateDependeObjects

InternalSideUserObject has a method called blockIDs()

InternalSideUserObject = DGKernel

So, I need DGKernel to have a blockIDs() method. Actually InternalSideUserObject inherits the blocks() method from the BlockRestrictable class, so after my inheritance, I should now have that method.

Modelling after KernelWarehouse, actions that happen on 3D dg_kernels that live in all volumes should happen in the same place that actions on 3D_kernels that live in just some volumes happen.

I'm looking to match the computeOffDiag threads between Kernels and DGKernels in ComputeFullJacobianThread since now DGKernels will be allowed to be block restricted

If I'm looking at changing the active dg kernel structure, I will have to make sure that I've addressed these:

src/base/ComputeFullJacobianThread.C:57:      const std::vector<KernelBase *> & kernels = kernel_warehouse.activeVar(ivar);
src/base/ComputeFullJacobianThread.C:81:        const std::vector<KernelBase *> & kernels = kernel_warehouse.activeVar(ivar.number());
src/base/ComputeFullJacobianThread.C:180:      std::vector<DGKernel *> dgks = dg_kernel_warehouse.activeVar(ivar);

Think about whether scaling might fix the liquid phase problem. Might be losing precision; that's possible. Make sure you have enough elements to resolve the gradients.

Jacobian could be ill-conditioned because of the dramatic difference in scales between the gas and the liquid, this could lead to a poor Newton step. Recall that the condition number is essentially the ratio of the largest singular value to the smallest singular value.

Think about possibility of iterating using multi-apps? E.g. splitting up your residuals and jacobians into two different problems with the idea of producing two different well conditioned Jacobians. Can you do that? Iterate within that one time step?

With nl_abs_tol = 1e-7, solve fails at 2.93e-7 seconds.

With nl_abs_tol = 1e-8, solve fails at 4.72e-7 seconds, so that's better.

Looking at condition numbers, it seems like stuff looks good!

The problem was the sign of the mobility of OHm. How ridiculous! It's always something with the physics!!!

Physics checks:

 - Gas electron and ion densities look fantastic.
 - Constantness of current in the gas phase looks fantastic.
 - Rholiq in the liquid looks fantastic.
 - Not current achieving continuity of electron concentration at the interface. Off by a factor of about 8
 - Achieving the correct change in electric field from gas to liquid. Changes by the ratio of the permittivities.
 - Not achieving continuity of electron current (directly proportional to total flux) across the interface. I'm pretty confident that that's because I've discretized the problem incorrectly for advection at the interface. I should be doing upwinding!

Moose active time = 420 seconds = 7 minutes ... pretty damn awesome :-)

Block 1 has elements [1,521] and nodes [1,522]
Block 2 has elements [522, 4706] and nodes [522,4707] (they share a node as they should)

-8.10658e-05 = electron current at interface on liquid side
-8.11127e-05 = electron current at interface on gas side

Those are pretty darn close! I would contend that that indicates current continuity.

May be time to move back to a single variable formulation for the elctrons and potential.

Without any electron stabilization, solve fails at 3.82592e-05 seconds, time step 110.

So far, my best run was with nl_abs_tol = 2e-9

Discontinuity in em between elements 725 and elements 726. Why???!!!

a: 8e-10
r: 1e-50
f: 3.26e-6

a: 1.6e-9
r: 1e-50
f: 3.23e-6

a: 3.2e-9
r: 1e-50
f: 3.38e-6

a: 6.4e-9
r: 1e-50
f: 3.84e-6

a: 1.3e-8
r: 1e-50
f: 4.45e-6

a: 2.6e-8
r: 1e-50
f: 5.99e-6
Without meeting abs tol, observed at least a 1e-6 decrease in the relative tolerance at the last diverged step

a: 5.2e-8
r: 1e-50
f: 5.93e-6

a: 1e-7
r: 1e-50
f: 6.02e-6
Observed oscillations

a: 2e-7
r: 1e-50
f: 5.77e-6

a: 2.6e-8
r: 1e-6
f: 6.02e-6

a: 2.6e-8
r: 1e-5
f: 9.60e-6
Looks like we may be hitting an absolute floor. Do I need to change my initial conditions?

Alright, I switched to a non-boundary condition transient.

a: 2.6e-8
r: 1e-5
f: 7.81e-6

a: 2e-6
r: 1e-6
Time step 312
time = 6.31e-5
dt = 1.01e-8

a: 4e-6
r: 1e-6
Solved to steady state
Yea, fuck the transient boundary condition. That pretty much artificially produces the need for two absolute tolerances.

Current at left boundary: -2391
Current at right boundary: -2398

I'd say that the above is a pretty good DC replica.

Active time was 1 minute on 1 processor. That's pretty nice.

In this simulation, the water width was only 100 nm. Screening specie concentrations were quite high.

I should double check and make sure that OHm reaction rates are correct and accurate. Presumably they should be since the OHm_lin concentration appears steady state.

Next thing to do: play around with some of the physics. Look again at extending out the liquid domain so that the em liquid concentration can effectively go to zero. See how that effects the total plasma densities.

Also, think about what could be done to make this simulation locally conservative, such that the electron current could be continuous across the interface.

Left potential = -.521882 kV = -521.882 Volts -> -522 V

Ion current = -1317 Amps/m^2
Electron current = -132.19 Amps/m^2
Total current = -1449.66 Amps/m^2

-> V_R = 727.73 Volts
-> V_p = 522.27 V -> 522 V

Perfect match in this case, so what's going on with the potential in the copuled plasma-liquid case?

With pot = 1e3 scaling, time steps = 169, voltage = ~2 Volts, abs_tol = 4e-6
1e6 scaling, tsteps = 198, voltage = 2.3 volts, abs_tol = 4e-6
1e7 scaling, tsteps = 296, voltage = 2.3 Volts, abs_tol = 4e-6
No scaling, tsteps = 170, voltage = 2.3 volts, abs_tol = 3e-6

So with scaling, the solve efficiency gets worse but the steady-state solution is the same.

Alright, after changing the electron mobility and diffusivity in the water to match their values in the gas, the left boundary potential condition is satisfied.

After turning off the advection_stabilization, I no longer got a discharge.

It's clear that 1e-8 is definitely too high of an absolute tolerance for early in the simulation.

With nl_abs_tol = 2e-9, nl_rel_tol = 1e-5, diverged line search at 1.78e-5 seconds with too high of residual

In LFA_lit_diff_mob_1e6_ballast_resist_no_H3Op_eps_80_gold_out.e, continuity of current and potential BC are both satisifed. Boom!!!

Emi does not include metastables in her PIC simulation. However, she most certainly does have excitation.

Unfortunately, with the current formulation for ElectronsFromIonization, the jacobian entries are so small that the kernel under investigations contribution to the matrix norms are too small to be registered. The finite difference jacobian doesn't even register any off-diagonal components even though they have to be present.

Good news!!! The finite differencing Jacobian will correctly handle material properties that are functions of dependent variables. Sweet!!!

Electron energy equation, unit checking check list:

ElectronTimeDerivative
EFieldAdvection
CoeffDiffusion
JouleHeating
ElectronEnergyLossFromIonization

Other kernel changes check list:

ElectronsFromIonization

Ok the trick variable appears to be em. That comes from studying IonsFromIonization. This is confirmed by studying ElectronsFromIonization.

So the question to answer is: should the on-diagonal of ElectronsFromIonization be non-zero even when _grad_em = 0 ? That will have to answered when I come back from my run.

When I get back, remove the interpolation of alpha_iz and its derivatives. Write down analytic expressions and see whether the weird behavior remains.

The problem here is an inescapable mathematical one. Trying to evaluate a function that contains a division by zero: in this case the electron flux magnitude. As soon as that value is non-zero the problems go away.

The problem in this case had nothing to do with interpolation.

Off diagonal jacobian for em is wrong.

Alright,
ElectronEnergyLossFromElastic no good
ElectronEnergyLossFromExcitation no good
ElectronEnergyLossFromIonization no good

With mean_en_log_stabilization at 27, fail time = 4.69e-6 seconds. nl_rel_tol = 1e-5
 nl_rel_tol = 1e-5
Stable = 22, fail = 1.09e-6, lowest_point = -34 nl_rel_tol = 1e-5
Stable = 17, fail = 2.69e-7, lowest_point = -17.5, no oscillations nl_rel_tol = 1e-5

Next going to try reducing nl_rel_tol

Stable = 27, nl_rel_tol = 1e-4, fail_time = 4.69e-6.  Same
Stable = 17, nl_rel_tol = 1e-4, nl_abs_tol = 5.6e-9, fail_time = 3.83e-7. Better
Stable = 17, nl_rel_tol = 1e-4, nl_abs_tol = 6e-7, fail_time = 4.75e-7. Better
Stable = 17, nl_rel_tol = 1e-4, nl_abs_tol = 8e-7, fail_time = 7.15e-7. Better

Electric field condition satisified at the interface.
Potential boundary condition satisfied.

We have a new warehouse structure in Moose, damn it. Alright, there are methods in NonlinearSystem for adding boundary conditions, dg kernels, etc. When we execute those methods' code we add objects to the attributes:

  MooseObjectWarehouse<DGKernel> _dg_kernels;
  MooseObjectWarehouse<KernelBase> _time_kernels;
  MooseObjectWarehouse<KernelBase> _non_time_kernels;

  ///@}

  ///@{
  /// BoundaryCondition Warhouses
  MooseObjectWarehouse<IntegratedBC> _integrated_bcs;
  MooseObjectWarehouse<NodalBC> _nodal_bcs;
  MooseObjectWarehouse<PresetNodalBC> _preset_nodal_bcs;

Alright, we're adding objects to the warehouses. So my question: do I want to make a new warehouse? And what should I inherit primarily from? DG or boundary condition?

What are the interface objects mathematically? They represent boundary conditions! I think they should be boundary condition objects. Now we need to investigate whether that is feasible. I think then it would have to inherit some properties of DGKernel, most importantly the access to variables on the other side of the interface.

In _factory.create(bc_name, name, ...), bc_name is the type of the object being constructed (e.g. from the input file, type = ____), and name is just the name of the object (e.g. in the input file [./right_side_bc]).

Factory::create returns a MooseObjectPtr

typedef MooseSharedPointer<MooseObject> MooseObjectPtr;

#ifdef LIBMESH_HAVE_CXX11_SHARED_PTR
#  include <memory>
#  define MooseSharedPointer std::shared_ptr
#  define MooseSharedNamespace std
#else

You see, here's where there might be a problem: DGInterface is a derived class of DGKernel so an _intc will be an _intc and a _dg_kernel :-(

That might be ok. The _dg_kernels in NonlinearSystem are not a reference to another variable; they are the root.

Alright, I need to go back and follow ComputeJacobianThread, ComputeFullJacobianThread, and DGInterface's jacobian calculation methods and clean things up and make them compile. That's the task at hand.

Good: the jacobians for the DGInterface test out clean.

So something is wrong with the coupling entries if I do block restriction with dg_kernels and apply SMP=true. If I don't apply SMP=true, then the problem solves perfectly fine. Additionally, if I do block restriction with just kernels and apply SMP=true, then the problem also solves perfectly fine.

Need to figure out how ce in ComputeFullJacobianThread is contstructed.

The assemblies get created in FEProblem here:

  _assembly.resize(n_threads);
  for (unsigned int i = 0; i < n_threads; ++i)
    _assembly[i] = new Assembly(_nl, couplingMatrix(), i);

So the coupling entries must come from couplingMatrix()

Alright, so SMP sets the coupling matrix. And for full=true, every element of the coupling matrix is set equal to 1.

Question, does offDiag get called for the kernel only case?

Alright all Moose tests pass, and my two jacobian check tests passed as well. All is well!!!!!!!

With mean_en offset at 17, solve fails at 2.06e-7 seconds, divet in mean_energy likely cause of solve failure.
Offset = 14, solve fails at 2.39e-7 seconds, no divet in mean_energy so that is not the reason for non-convergence. Rather it looks like abs_tol too high (2e-6)
abs_tol = 4e-6, fail = 3.61e-7

Changing the offset from 14 to 15 didn't affect the final results; that's good.

With r = 0.5, # time steps = 180.
r = 0.9, time steps = 192.

abs = 6e-8, fail = 4.20e-6
abs = 1.2e-7, rel = 1e-5, fail = 1.3e-7
abs = 1.2e-7, rel = 1e-4, fail = 1.3e-7

Assumptions that are better not to make:

1) Better not to assume that D = mu * Te
2) Better not to assume that mu_eps = 5/3 * mu, D_eps = 5/3 * D

Quantities that are defined in bolos:

mobility
diffusion
mean_energy
rate (the rate of any given process that x_sections were provided for)

solver.py contains the class BoltzmannSolver

To get a truly accurate boundary condition, a transport equation with vx, vy, and vz (or E, omx, omy, and omz) would have to be used.

From Duderstadt:

Unfortunately diffusion theory is capable of only appromixating even the integral condition, since it cannot yield the exact form for J_{+/-}.

Just accept and be happy.

Alright, final decision: I'm going to use Hagelaar's formulation for the electrons for my electron boundary conditions (and for the ion boundary conditions (but I don't have to worry about secondaries)) and finally for the electron energy. To obtain the electron energy I'm simply going to multiply the alpha electron boundary condition by 5/3 mean_energy_alpha and the gamma electron boundary condition by 5/3 mean_energy_gamma (e.g. I'm multiplying all electron speeds by 5/3;). Alright I'm done with thinking about this. My justification for these decisions comes from juxtaposing Hagelaar's BC paper and Yuki's electron energy boundary condition in the LFA vs mean_en comparison paper. There done!!!!! Damn it!


From FEProblem::checkProblemIntegrity, the checkBoundaryMatProps() function call is commented out. Let's imagine that it's uncommented. What happens: it calls SubProblem::checkBoundaryMatProps() and executes:

  checkMatProps(_map_boundary_material_props, _map_boundary_material_props_check, "boundary");


_map_boundary_material_props_check I believe is the map that stores what properties !should! exist on a given boundary ID.

Modified by: SubProblem::storeDelayedCheckMatProp

There are two different implementations of this method:

void
SubProblem::storeDelayedCheckMatProp(const std::string & requestor, SubdomainID block_id, const std::string & name)
{
  _map_block_material_props_check[block_id].insert(std::make_pair(requestor, name));
}

void
SubProblem::storeDelayedCheckMatProp(const std::string & requestor, BoundaryID boundary_id, const std::string & name)
{
  _map_boundary_material_props_check[boundary_id].insert(std::make_pair(requestor, name));
}

What calls storeDelayedCheckMatProp?

MaterialPropertyInterface::checkMaterialProperty calls storeDelayedCheckMatProp twice, once supposedly for blocks and once supposedly for boundaries. Yes it should call the correct storeDelayed method for both.

Now this is called only if _mi_boundary_ids has some stuff in it. How does that attribute get modified?

And what is _mi_name?

IntegratedBC inherits from MaterialPropertyInterface

u:0,4 em
v:1,5 pot
w:2,6 mean_en
p:3,7 ip

Recall that analytical derivatives (when correct) are fundamentally  more accurate than finite difference derivatives

Tests that failed when I tried implementing my boundary materials check:

misc/line_source.test........................................................................ FAILED (ERRMSG)
postprocessors/misc_pps.misc_pps_test........................................................ FAILED (ERRMSG)
materials/boundary_material.elem_aux_bc_on_bnd............................................... FAILED (ERRMSG)
materials/material.dg_test................................................................... FAILED (ERRMSG)
materials/material.test...................................................................... FAILED (ERRMSG)
materials/stateful_prop.spatial_test......................................................... FAILED (ERRMSG)
materials/stateful_prop.stateful_copy........................................................ FAILED (ERRMSG)
vectorpostprocessors/line_material_sampler.test.............................................. FAILED (ERRMSG)

Based off what Cody said, I think that my incredibly crude, hacked implementation should work fine for cases on which I'm only using external boundaries.

What I did today:

Implemented 5 new residual/jacobian classes. Successfully ran a simulation with all of Hagelaar's boundary conditions implemented. Discovered that they don't really change my solution!

WIP: implement energy dependence of muem, muel, diffem, and diffel

0 mesh refinement (522 vertices): fail time = 1.25e-6 (oscillations in ions)
1 mesh refinement (681 vertices): fail time = 1.26e-6
2 mesh refinments (765 vertices): fail time = 1.33e-6 (no more oscillations), nl_rel_tol = 1e-5, nl_abs_tol = 2e-9

I'm seeing no oscillations in the ion density now with a good mesh. No divets in the mean energy, so that says to me the problem is not an insufficient amount of stabilization. Also the last converged result had a residual of 1e-9

rel = 1e-4, abs = 2e-9, fail = 1.33e-6, ion offset = 17, mean_en offset = 27, em offset = 50

rel = 1e-4, abs = 2e-9, fail = 1.35e-10, ion offset = 27, mean_en offset = 27, em offset = 50

rel = 1e-4, abs = 2e-9, fail = 5.60e-11, ion offset = 17, mean_en offset = 17, em offset = 17

765 vertices, rel = 1e-4, abs = 2e-9, offsets = 40, ics = 28, fail = 4.26e-10 (oscillations in ions)

1677 vertices, rel = 1e-4, abs = 2e-9, offsets = 40, ics = 28, fail = no time step
1677 vertices, rel = 1e-4, abs = 2e-9, offsets = 50, ics = 28, fail = no time step
1677 vertices, rel = 1e-4, abs = 2e-9, offsets = 50, ics = 31, fail = 2.56e-11

Getting obscene oscillations in the ion density. In fact the more refined I make the mesh, the bigger the problem is.

Mean_en offset = 29, fail at 1.59e-9 seconds
Mean_en offset = 26, solve completed but not physical

Restored correct units, fail at 5.08e-7 seconds
Mean_en offset = 23, fail at 6.04e-7 seconds

Offsets = 50, fail at 1.23e-5 seconds, old mean_energy transport data
New mean_energy transport data, failt at 1.23e-5 seconds again

Guess what? At the end of the day, the boundary material implementation broke my app too. As soon as I removed the boundary material, my app worked perfectly fine again.

I lied, I didn't break it with the boundary material. So I guess what must have fixed everything was moving to offset = 50 for all objects?

It turns out that both the source and flux stabilization techniques are critical for a successful simulation. Particularly for the argon ions.

Simulation takes 209 time steps to solve with the interpolation of the transport coefficients.

I guess it's evidence enough that the simulation is doing something different since it fails when I try to run the exact same input file, only with interpolation of transport coefficients turned off.

Question is: how much ion current is diffusive and how much is from advection? Tomorrow's task is to correct the Neumann boundary condition (which also means doing a butt load of analytical derivatives.)

Yes currently the boundary conditions aren't matching: based off the total current through the plasma, the potential at the wall should be 475 volts, but instead it's 540 volts. Gotta fix that tomorrow.

On my own time, I would like to work on the expression simplification in sympy.

Ok, the numerator derivative is correct

Without ion stabilization, solve fails at 9.73e-6 seconds. Going to refine at left boundary
New fail at 9.00e-6 seconds.

After second level of refinement, solve completed successfully.

It's been my observation that at least with respect to log stabilization, as little as possible actually enhances convergence (although some non-zero amount is certainly required in the case of ions at the anode because there is no source for them like there are for electrons at the cathode). With respect to advection stabilization, under normal circumstances any oscillations arising from the domination of advection can be fixed with mesh refinement and without enlisting the aid of stabilization.

Alright, I'm using this as a scale length estimate for the change in the electric field: |E|/|dE/dx|. Using that I get about 40 to 55 microns in the cathode. This is actually indeed long on the scale of electron-neutral collisions.

pt999 and pt9999 had appropriate boundary conditions.

em dofs: 0, 1
emliq dofs: 2, 3
potential dofs: 4,6
mean_en dofs: 5,7

_mu_neighbor is a const reference to a Real MaterialProperty.

My prediction is that each object is associated with a geomtric element. And thus there will be many different kernel objects with a one to one correspondence to the number of elements there are. Well that doesn't appear to be true there...

Element1 JacMat: 0xa966b0
Left boundary jacmat: 0xaa1c50
Element2 JacMat: 0xa966b0

There was a second alternative: a single material object and a single kernel object that represent all similar element types, where the process is:

element 0:call Material (qp 0, qp1)
element 0:call Kernel (i0: qp 0, qp1; i1: qp0, qp1)
element 1:call material (qp0, qp1)
element1:call kernel (qp0, qp1)

so on and so forth. This is actually how the code works. The material properties are computed right before they are needed by the residual creating object

Ok this looks promising, there are two material objects that live on the geomtric object associated with my InerfaceKernel and they correspond to the correct sides. Yet that material object doesn't appear to have access to the variables.

#0  0x00007ffff78a77f0 in JacMat::JacMat(InputParameters const&)@plt ()
   from /home/lindsayad/projects/zapdos/lib/libzapdos-dbg.so.0
#1  0x00007ffff79646c5 in buildObject<JacMat> (parameters=...)
    at /home/lindsayad/moose/framework/include/base/Factory.h:130
#2  0x00007ffff67fdd56 in Factory::create (this=0x93c2a0, obj_name=..., name=..., parameters=...,
    tid=0) at /home/lindsayad/moose/framework/src/base/Factory.C:77
#3  0x00007ffff64b5821 in FEProblem::addMaterial (this=0xa51e30, mat_name=..., name=...,
    parameters=...) at /home/lindsayad/moose/framework/src/base/FEProblem.C:1626
#4  0x00007ffff6acbb37 in AddMaterialAction::act (this=0x9e3d50)
    at /home/lindsayad/moose/framework/src/actions/AddMaterialAction.C:32
#5  0x00007ffff6b00c3d in ActionWarehouse::executeActionsWithAction (this=0x93baa8, task=...)
    at /home/lindsayad/moose/framework/src/actions/ActionWarehouse.C:332
#6  0x00007ffff6b0081f in ActionWarehouse::executeAllActions (this=0x93baa8)
    at /home/lindsayad/moose/framework/src/actions/ActionWarehouse.C:308
#7  0x00007ffff66bb8d0 in MooseApp::runInputFile (this=0x93ae80)
    at /home/lindsayad/moose/framework/src/base/MooseApp.C:307
#8  0x00007ffff66bcea2 in MooseApp::run (this=0x93ae80)
    at /home/lindsayad/moose/framework/src/base/MooseApp.C:510
#9  0x0000000000420956 in main (argc=3, argv=0x7fffffffcf98)
    at /home/lindsayad/projects/zapdos/src/main.C:23

All objects appear to be built before any simulation starts. Thus it's possible that objects are built before initial values are read in.

Block 0 element Mat: 0xa8a4f0

Left boundary Mat: 0xa95a50

Left interface Mat: 0xa9d9d0


Block 1 element Mat: 0xaa3690

Right interface Mat: 0xaaaa40
Right boundary Mat: 0xaaaa40

This is interesting. Maybe the left interface should be the same object as the left boundary since the right interface is the same as the right boundary.

One can require coupled variables in a material class, and as long as those parameters are input in the input file, then moose won't complain. It doesn't even matter if the required variables don't even exist on the material block; moose won't complain.

Yes it looks like I can get the neighboring block material properties but am currently stuck with variable values in the material class corresponding to the element block.

Ok, one single block for JacMat:

Element Mats: 0xa85c90 (both left and right blocks)
boundary mats: 0xa922b0 (left boundary, right interface has acces to _emliq (master (right) side of InterfaceKernel))
neighbor mat: 0xa97530 (left interface has access to _emliq (master (right) side of InterfaceKernel))

So the problem appears to be that Coupleable is initialized before _neighbor gets set to be true.

I believe the cause of divergence is when the electron temperature starts to get too low (or too high!!!). I'm guessing that this produces problems for the interpolation?

Fail time with unaltered bolos data: 1.49e-5 seconds

After altering the property input file so that I effectively ensured positive and realistic order of magnitude estimates of all the coefficients, then the simulation converged on the next run! Look at what we're bloody learning about convergence! We're becoming convergence frigging masters!!!

Getting the thermodynamic simulation to converge is much more difficult than getting the kinetic version to converge. I am pretty confident that this is because with my boundary conditions, the electron temperature goes to zero which can create numerical issues.

To the collect method, I pass an expr. Let's say for example I am passing collect3. The type of collect3 is:

sympy.core.mul.Mul

This type appears to have the member function func. Let's try to learn about func.

mul has some members. These include Mul and prod. Mul also has some members. These include as_coeff_Mul, as_content_primitive, and others. However, I don't see func listed.

What do the multiple dots mean in python? In c++ there is only ever one dot, not a chain of dots. I'm in the process of researching this now.

['Add',
 'AssocOp',
 'Basic',
 'C',
 'Expr',
class 'Mul',
class 'NC_Marker',
 'Pow',
 'Rational',
 'S',
 '__builtins__',
 '__doc__',
 '__file__',
 '__name__',
 '__package__',
 '_addsort',
 '_args_sortkey',
def '_keep_coeff',
def '_mulsort',
 '_unevaluated_Add',
def '_unevaluated_Mul',
 'cacheit',
 'cmp_to_key',
 'defaultdict',
 'division',
def 'expand_2arg',
 'fuzzy_not',
 'operator',
 'print_function',
def 'prod',
 'reduce',
 'sympify',
 'xrange']

It seems that the mesh may play a major major role in whether something converges or not. Remember that with bad ProductAABBRxn Jacobian, the thermo simulation took 134 seconds to solve: that's not true. I wash still using the same executable, you big dumb dodo bird. With the correct Jacobian it takes the same amount of time. Thus the reason for my non-convergence with the kinetic simulation must be the mesh. It is not any other changes.

Ok with smallest mesh feature = 2e-9, fail time = 2.69e-7 seconds.

New kinetic mesh. Coarse liquid region. Fail time = 2.75e-7
Finer liquid region: 2.57e-7

Coarser mesh: fail time = 2.82e-7

Alright we have proof of things. I made the mesh coarser once more at the interface and the simulation solved beautifully. Thus in my opinion the lack of convergence before was coming from poor conditioning. I have seen time and again that bad scaling leads to bad convergence. I have scaled my dependent variables nicely so that they are not so far from 1. However, I have not scaled my position or time variables. Time might be ok...it starts at nanosecond values but then increases through the simulation. But I'm using nanometer type distances in my mesh and this doesn't change throughout the simulation. Just think about a gradient term where the dependent variable is not changing very much: then we have a small number divided by a big number: this leads to bad conditioning! Ok we're still learning things about convergence as we go here.

Dr. Graves was right: meshing can play a big old role in convergence!!!!

So apparently python will look over a whole method to determine whether variables are global or local and then execute the statements in the method. If at anyplace in a method, a global variable is reassigned without that variable being declared global in the method, then that variable will be treated as a local variable.

0a) Observation: total current changes by only 17\% between gamma extremes
0b) Consequently the total potential drop across the plasma only varies by that amount
1) From 0b): Total potential drop more or less constant
2) Sound backing: Build up in electrons at anode because of gamma
3) From 0a): Total current more or less constant. In order for this to be true, the cathode characteristics must be the same because we aren't changing the cathodic boundary conditions.
4) From 2) electric field becomes large and negative, leading to a large potential drop in the anode.
5) From 1) and 4) for smaller gammas, the potential drop in the bulk will be smaller because more potential is being lost in the anode.
6) From 5) the bulk electric field is smaller for smaller gammas
7) With a smaller bulk electric field, in order to maintain the same total current between cases, the *density* of electrons (and with quasineutrality consequently the ions) must be higher!

boom! Only thing that would make the argument more rock solid is if I could come up with a theoretical descripton for 0).

kinetic_r_0pt9_no*_gold.i ran without a hitch with the current commit version of liquidNew.geo
kinetic_r_0pt9999_no*_gold.i ran without a hitch with the current commit version of liquidNew.geo (so there was actually no need to run this with the new mesh, should have just added the Auxkernels to this input file)

Plots needed:

dens kinetic, gamma = 1, energy kinetic, gamma = 1, check
dens kinetic, gamma = 1e-4, energy kinetic, gamma = 1e-4, check
dens kinetic, gamma = 1, energy zero grad check (but poor liquid meshing)
dens kinetic, gamma = 1e-4, energy zero grad check (but poor liquid meshing)
dens thermo, H = 1, energy kinetic, gamma = 1, check (corresponds to high reflection density and low reflection energy -> inconsistent)
dens thermo, H = 1e5, energy kinetic, gamma = 1, check (corresponds to low reflection density and low reflection energy -> consistent)
dens thermo, H = 1e5, energy kinetic, gamma = 1e-4, check (corresponds to low reflection density and high reflection energy -> inconsistent)
dens thermo, H = 1, energy kinetic, gamma = 1e-4, check (corresponds to high reflection density and high reflection energy -> consistent)
dens thermo, H = 1, energy zero grad check
dens thermo, H = 1e5, energy zero grad check

Low gamma: high reflection
High H: low reflection

Next task: get the last case done, plot, and then add figures to paper. Figure out what figures should be kept for paper and start writing as Steve said.

MatPlotLib appears to only have 7 original colors. After 7 plots are made, then it starts to recycle colors from the top (blue, green, red, etc...)

Middle mesh: 100e-6, fail time = 1.09e-7
Middle mesh: 50e-6, right mesh = 5e-8, fail time = 1.09e-7, nl_abs_tol = 1e-12
Middle: 50e-6, right: 2e-8, fail: 1.10e-7
Middle: 50e-6, right: 1e-8, fail: 1.10093e-7
Middle: 50e-6, right: 5e-9, fail: 1.10151e-7, offset = 20
Middle: 50e-6, right: 5e-9, fail: 3.5833e-7, offset = 25
Middle: 50e-6, right: 5e-9, fail: 2.31938e-6, offset = 30
Middle: 20e-6, right: 5e-9, fail: 1.078e-7, offset = 20
Middle: 100e-6, right: 5e-9, fail: 1.08965e-7, offset = 20
Middle: 100e-6, right: 1e-9, fail: 1.09003e-7, offset = 20
Middle: 50e-6, right: 1e-9, fail: 1.10205e-7, offset = 20 (no oscillations or any problems in the solution...thus this has to be a conditioning problem in my opinion)
Middle: 50e-6, right: 1e-9, fail: 1.09347e-7, offset = 20, mean_en_offset = 15

I think lowering the offset just delays the onset of the discharge, it doesn't really improve convergence.

Could it be mesh conditioning or low electron temperatures or what?

So I proved that the gold H = 1, gamma = 1e-4 anomalous ion and electron temperature results arises solely from a glitch in initial conditions. If I only modify the initial conditions, then I can get rid of the ridiculous jump in electron temperature and ion density at the interface. Weird stuff.

Getting good convergence now after adding _test[_i][_qp] to the HagelaarEnergyBC. That could very well have been the reason that I've been struggling to get convergence with that boundary condition. Wowzers. Yea sure the Jacobian could have perfectly matched the residual, but if I'm supply a numerical method that's intrinsically unstable, then even with a perfect Jacobian, the problem won't converge. Wow, sometimes the weirdest things save you.

The bad dgkernel offender is: InterfaceLogDiffusionElectrons

Fix tomorrow. Fixed

em: 26
emliq: 26
Arp: 26
mean_en: 25
OHm: 15.6


0 Nonlinear |R| = 8.860315e+04
      Line search: gnorm after quadratic fit 1.885694427347e+05
      Line search: Cubically determined step, current gnorm 8.417309332120e+04 lambda=5.0000000000000003e-02
    |residual|_2 of individual variables:
                  potential: 84075
                  em:        0.0604987
                  emliq:     130.652
                  Arp:       0.0603611
                  mean_en:   3.29598
                  OHm:       4060.41

 1 Nonlinear |R| = 8.417309e+04
      Line search: Using full step: fnorm 8.417309332120e+04 gnorm 3.424562095658e+01
    |residual|_2 of individual variables:
                  potential: 0.14715
                  em:        0.0231411
                  emliq:     34.2238
                  Arp:       0.0202218
                  mean_en:   1.21432
                  OHm:       3.44686e-06

 2 Nonlinear |R| = 3.424562e+01

em: 26
emliq: 31
Arp: 26
mean_en: 25
OHm: 15.6

 0 Nonlinear |R| = 8.860315e+04

em: 26
emliq: 21
Arp: 26
mean_en: 25
OHm: 15.6


0 Nonlinear |R| = 8.860315e+04
      Line search: Using full step: fnorm 8.860314992199e+04 gnorm 3.278313816951e+01
    |residual|_2 of individual variables:
                  potential: 0.154749
                  em:        0.0240928
                  emliq:     32.7579
                  Arp:       0.021509
                  mean_en:   1.27748
                  OHm:       3.61508e-06

 1 Nonlinear |R| = 3.278314e+01

em: 21
emliq: 21
Arp: 21
mean_en: 20
OHm: 15.6

 0 Nonlinear |R| = 8.860872e+04
      Line search: gnorm after quadratic fit 7.974510602290e+04
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 79650
                  em:        0.569602
                  emliq:     601.907
                  Arp:       0.000610509
                  mean_en:   2.22051
                  OHm:       3846.7

 1 Nonlinear |R| = 7.974511e+04
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1
 Solve Did NOT Converge!


Got initial steps to go just fine after actually doing the correct meshing (dummy!).

Fail time with no scaling: 7.67636e-8
1e-3 scaling of liq phase: 7.60823e-8
1e-4 scaling of liq phase: 7.60823e-8
1e-4 scaling of liq phase, 1e-1 scaling of mean_en: 7.65912e-8

After increasing the maximum mean energy in my td_argon_mean_en.txt file: 7.86949e-8 (the electron temp had exceeded 20 volts at the anode)

After removing scaling: .0001102 (crashed because of wrong abs_tol)

So initially I identified the wrong problem: it wasn't a problem of scaling with the variables...in fact their relative scaling is probably just about in the perfectly optimable zone (there's a reason for the things we do and for the directions we go! Convergence leads us down the right path!). And the problem wasn't probably even with the scaling of the mesh. The problem was with the electron temperature exceeding my upper bound in my lookup table leading to bad numerics from that. And then the position scaling just made the solve go that much faster!

If Dr. Graves quibbles about factors, I will say that those nv/4 and 2kT flux quantities are derived for a Maxwellian distribution which assumes the absence of time variation, spatial gradients, and accelerations. In my DC discharge at the anode I have both spatial gradients and accelerations due to electric fields so those well known equations are not valid. For any further argument I would refer him to the Hagelaar paper which shows that many authors incorrectly leave off the -D/2 grad_n term of the boundary condition (this BC is familiar to me from neutran diffusion theory! The vacuum boundary condition!).

Current number of lines of code in .C and .h files in Zapdos: 23410

Yep, if things don't converge, it could be because:

1) Oscillations at large variable values
2) Bad material fitting (e.g. linear or spline). Check spline fits (see plot_mat_data.py)
3) Bad variable scaling
4) Problem is just amazingly non-linear. In this case, you can try ramping the non-linearity. I did this with great success for gamma_dens = 1e-4 and gamma_en = 1e-3. Without the ramping, I was getting incredibly inefficient convergence, e.g. overnight simulation without reaching steady state. However, adding in the nonlinear ramping, got to steady state in 93 seconds of simulation time. Nice work!

Unscaled solution time: 29.069 seconds
Scaled solution time: 13 seconds

JxW for a Kernel compute residual: 1.04e-6
JxW for IntegratedBC @ 0,0,0: 1

After removing the boundary condition and interface scaling factors, the solve no longer converged.

elem: 0xd4e5e0
side: 0
bnd_id: 1
interface kernel name: em_advection

elem: 0xd4e5e0
side: 0
bnd_id: 1
interface kernel name: em_diffusion

(this=0x7fffffffc810, elem=0xd4e690, side=0, bnd_id=1) passes the active interface kernels test.

(this=0x7fffffffc810, elem=0xd4e5e0, side=1, bnd_id=2) does not pass the active interface kernels test

(this=0x7fffffffc340, elem=0xd4e5e0, side=1, bnd_id=2) also does not pass

(this=0x7fffffffc340, elem=0xd4e690, side=0, bnd_id=1) does pass. elem centroid = 1.0000005489319999. elem id = 211

My question is: why does it get called twice??

There are four calls to NonlinearSystem::computeResidualInternal:

bt1:

#0  0x00007ffff6186690 in NonlinearSystem::computeResidualInternal(Moose::KernelType)@plt () from /home/lindsayad/moose/framework/libmoose-dbg.so.0
#1  0x00007ffff66de473 in NonlinearSystem::computeResidual (this=0xd62f80, residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:750
#2  0x00007ffff64a7493 in FEProblem::computeResidualType (this=0xd60110, soln=..., residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3334
#3  0x00007ffff64a6e33 in FEProblem::computeResidual (this=0xd60110, soln=..., residual=...) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3270
#4  0x00007ffff66da9cf in NonlinearSystem::solve (this=0xd62f80) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:242
#5  0x00007ffff64a5ed1 in FEProblem::solve (this=0xd60110) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3049
#6  0x00007ffff6c6eff3 in TimeStepper::step (this=0xd89e00) at /home/lindsayad/moose/framework/src/timesteppers/TimeStepper.C:188
#7  0x00007ffff6c6049a in Transient::solveStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:413
#8  0x00007ffff6c5fff2 in Transient::takeStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:340
#9  0x00007ffff6c5fba9 in Transient::execute (this=0xd932f0) at /home/lindsayad/moose/framework/src/executioners/Transient.C:253
#10 0x00007ffff66a10be in MooseApp::executeExecutioner (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:373
#11 0x00007ffff66a1ec2 in MooseApp::run (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:523
#12 0x0000000000420ac6 in main (argc=3, argv=0x7fffffffd078) at /home/lindsayad/projects/zapdos/src/main.C:23

bt2:

#0  NonlinearSystem::computeResidualInternal (this=0xd62f80, type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:1196
#1  0x00007ffff66de473 in NonlinearSystem::computeResidual (this=0xd62f80, residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:750
#2  0x00007ffff64a7493 in FEProblem::computeResidualType (this=0xd60110, soln=..., residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3334
#3  0x00007ffff64a6e33 in FEProblem::computeResidual (this=0xd60110, soln=..., residual=...) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3270
#4  0x00007ffff66da9cf in NonlinearSystem::solve (this=0xd62f80) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:242
#5  0x00007ffff64a5ed1 in FEProblem::solve (this=0xd60110) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3049
#6  0x00007ffff6c6eff3 in TimeStepper::step (this=0xd89e00) at /home/lindsayad/moose/framework/src/timesteppers/TimeStepper.C:188
#7  0x00007ffff6c6049a in Transient::solveStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:413
#8  0x00007ffff6c5fff2 in Transient::takeStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:340
#9  0x00007ffff6c5fba9 in Transient::execute (this=0xd932f0) at /home/lindsayad/moose/framework/src/executioners/Transient.C:253
#10 0x00007ffff66a10be in MooseApp::executeExecutioner (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:373
#11 0x00007ffff66a1ec2 in MooseApp::run (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:523
#12 0x0000000000420ac6 in main (argc=3, argv=0x7fffffffd078) at /home/lindsayad/projects/zapdos/src/main.C:23

bt4:

#0  NonlinearSystem::computeResidualInternal (this=0xd62f80, type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:1196
#1  0x00007ffff66de473 in NonlinearSystem::computeResidual (this=0xd62f80, residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:750
#2  0x00007ffff64a7493 in FEProblem::computeResidualType (this=0xd60110, soln=..., residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3334
#3  0x00007ffff64a6e33 in FEProblem::computeResidual (this=0xd60110, soln=..., residual=...) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3270
#4  0x00007ffff66d93d3 in Moose::compute_residual (soln=..., residual=..., sys=...) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:100
#5  0x00007ffff3a55c54 in libMesh::__libmesh_petsc_snes_residual (snes=0x102df70, x=0xd69560, r=0xe34470, ctx=0xd652f0) at ../src/solvers/petsc_nonlinear_solver.C:130
#6  0x00007ffff04db29b in SNESComputeFunction () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#7  0x00007ffff04bc1b6 in SNESSolve_NEWTONLS () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#8  0x00007ffff04e6719 in SNESSolve () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#9  0x00007ffff3a57d32 in libMesh::PetscNonlinearSolver<double>::solve (this=0xd652f0, jac_in=..., x_in=..., r_in=...) at ../src/solvers/petsc_nonlinear_solver.C:685
#10 0x00007ffff3acd23e in libMesh::NonlinearImplicitSystem::solve (this=0xd64740) at ../src/systems/nonlinear_implicit_system.C:183
#11 0x00007ffff62c0088 in TimeIntegrator::solve (this=0xd94600) at /home/lindsayad/moose/framework/src/timeintegrators/TimeIntegrator.C:54
#12 0x00007ffff66dab2a in NonlinearSystem::solve (this=0xd62f80) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:264
#13 0x00007ffff64a5ed1 in FEProblem::solve (this=0xd60110) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3049
#14 0x00007ffff6c6eff3 in TimeStepper::step (this=0xd89e00) at /home/lindsayad/moose/framework/src/timesteppers/TimeStepper.C:188
#15 0x00007ffff6c6049a in Transient::solveStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:413
#16 0x00007ffff6c5fff2 in Transient::takeStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:340
#17 0x00007ffff6c5fba9 in Transient::execute (this=0xd932f0) at /home/lindsayad/moose/framework/src/executioners/Transient.C:253
#18 0x00007ffff66a10be in MooseApp::executeExecutioner (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:373
#19 0x00007ffff66a1ec2 in MooseApp::run (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:523
#20 0x0000000000420ac6 in main (argc=3, argv=0x7fffffffd078) at /home/lindsayad/projects/zapdos/src/main.C:23

Then a residual is output to the terminal: 0 Nonlinear |R| = 8.860873e+04

How many bt's now?

bt2:

#0  NonlinearSystem::computeResidualInternal (this=0xd62f80, type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:1196
#1  0x00007ffff66de473 in NonlinearSystem::computeResidual (this=0xd62f80, residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:750
#2  0x00007ffff64a7493 in FEProblem::computeResidualType (this=0xd60110, soln=..., residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3334
#3  0x00007ffff64a6e33 in FEProblem::computeResidual (this=0xd60110, soln=..., residual=...) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3270
#4  0x00007ffff66d93d3 in Moose::compute_residual (soln=..., residual=..., sys=...) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:100
#5  0x00007ffff3a55c54 in libMesh::__libmesh_petsc_snes_residual (snes=0x102df70, x=0x1028530, r=0x1230b20, ctx=0xd652f0) at ../src/solvers/petsc_nonlinear_solver.C:130
#6  0x00007ffff04db29b in SNESComputeFunction () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#7  0x00007ffff04f8d01 in SNESLineSearchApply_BT () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#8  0x00007ffff04f520e in SNESLineSearchApply () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#9  0x00007ffff04bd0b6 in SNESSolve_NEWTONLS () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#10 0x00007ffff04e6719 in SNESSolve () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#11 0x00007ffff3a57d32 in libMesh::PetscNonlinearSolver<double>::solve (this=0xd652f0, jac_in=..., x_in=..., r_in=...) at ../src/solvers/petsc_nonlinear_solver.C:685
#12 0x00007ffff3acd23e in libMesh::NonlinearImplicitSystem::solve (this=0xd64740) at ../src/systems/nonlinear_implicit_system.C:183
#13 0x00007ffff62c0088 in TimeIntegrator::solve (this=0xd94600) at /home/lindsayad/moose/framework/src/timeintegrators/TimeIntegrator.C:54
#14 0x00007ffff66dab2a in NonlinearSystem::solve (this=0xd62f80) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:264
#15 0x00007ffff64a5ed1 in FEProblem::solve (this=0xd60110) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3049
#16 0x00007ffff6c6eff3 in TimeStepper::step (this=0xd89e00) at /home/lindsayad/moose/framework/src/timesteppers/TimeStepper.C:188
#17 0x00007ffff6c6049a in Transient::solveStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:413
#18 0x00007ffff6c5fff2 in Transient::takeStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:340
#19 0x00007ffff6c5fba9 in Transient::execute (this=0xd932f0) at /home/lindsayad/moose/framework/src/executioners/Transient.C:253
#20 0x00007ffff66a10be in MooseApp::executeExecutioner (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:373
#21 0x00007ffff66a1ec2 in MooseApp::run (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:523
#22 0x0000000000420ac6 in main (argc=3, argv=0x7fffffffd078) at /home/lindsayad/projects/zapdos/src/main.C:23

bt 4:

#0  NonlinearSystem::computeResidualInternal (this=0xd62f80, type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:1196
#1  0x00007ffff66de473 in NonlinearSystem::computeResidual (this=0xd62f80, residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:750
#2  0x00007ffff64a7493 in FEProblem::computeResidualType (this=0xd60110, soln=..., residual=..., type=Moose::KT_ALL) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3334
#3  0x00007ffff64a6e33 in FEProblem::computeResidual (this=0xd60110, soln=..., residual=...) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3270
#4  0x00007ffff66d93d3 in Moose::compute_residual (soln=..., residual=..., sys=...) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:100
#5  0x00007ffff3a55c54 in libMesh::__libmesh_petsc_snes_residual (snes=0x102df70, x=0x1028530, r=0x1230b20, ctx=0xd652f0) at ../src/solvers/petsc_nonlinear_solver.C:130
#6  0x00007ffff04db29b in SNESComputeFunction () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#7  0x00007ffff04f95a6 in SNESLineSearchApply_BT () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#8  0x00007ffff04f520e in SNESLineSearchApply () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#9  0x00007ffff04bd0b6 in SNESSolve_NEWTONLS () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#10 0x00007ffff04e6719 in SNESSolve () from /opt/moose/petsc/openmpi_petsc-3.6.0/gcc-opt/lib/libpetsc.so.3.6
#11 0x00007ffff3a57d32 in libMesh::PetscNonlinearSolver<double>::solve (this=0xd652f0, jac_in=..., x_in=..., r_in=...) at ../src/solvers/petsc_nonlinear_solver.C:685
#12 0x00007ffff3acd23e in libMesh::NonlinearImplicitSystem::solve (this=0xd64740) at ../src/systems/nonlinear_implicit_system.C:183
#13 0x00007ffff62c0088 in TimeIntegrator::solve (this=0xd94600) at /home/lindsayad/moose/framework/src/timeintegrators/TimeIntegrator.C:54
#14 0x00007ffff66dab2a in NonlinearSystem::solve (this=0xd62f80) at /home/lindsayad/moose/framework/src/base/NonlinearSystem.C:264
#15 0x00007ffff64a5ed1 in FEProblem::solve (this=0xd60110) at /home/lindsayad/moose/framework/src/base/FEProblem.C:3049
#16 0x00007ffff6c6eff3 in TimeStepper::step (this=0xd89e00) at /home/lindsayad/moose/framework/src/timesteppers/TimeStepper.C:188
#17 0x00007ffff6c6049a in Transient::solveStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:413
#18 0x00007ffff6c5fff2 in Transient::takeStep (this=0xd932f0, input_dt=-1) at /home/lindsayad/moose/framework/src/executioners/Transient.C:340
#19 0x00007ffff6c5fba9 in Transient::execute (this=0xd932f0) at /home/lindsayad/moose/framework/src/executioners/Transient.C:253
#20 0x00007ffff66a10be in MooseApp::executeExecutioner (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:373
#21 0x00007ffff66a1ec2 in MooseApp::run (this=0x93c4d0) at /home/lindsayad/moose/framework/src/base/MooseApp.C:523
#22 0x0000000000420ac6 in main (argc=3, argv=0x7fffffffd078) at /home/lindsayad/projects/zapdos/src/main.C:23

{<ThreadedElementLoop<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >> = {<ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >> = {_vptr.ThreadedElementLoopBase = 0x7ffff761db10 <vtable for ComputeResidualThread+16>, _mesh = @0xd37270, _tid = 0, _subdomain = 0, _old_subdomain = 0}, _system =
    @0xd62f80, _fe_problem = @0xd60110}, _sys = @0xd62f80, _kernel_type = Moose::KT_ALL, _num_cached = 0, _integrated_bcs = @0xd63a78, _dg_kernels = @0xd634f8, _interface_kernels = @0xd63658, _kernels =
    @0xd631f0, _time_kernels = @0xd637b8, _non_time_kernels = @0xd63918}

Breakpoint 3, ComputeResidualThread::onInterface (this=0x7fffffffc210, elem=0xd4e5d0, side=1, bnd_id=2) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:144
144	  if (_interface_kernels.hasActiveBoundaryObjects(bnd_id, _tid))

_interface_kernels address from inside ComputeResidualThread::onInterface =  @0xd63658

Matches that from cr within NonlinearSystem::computeResidualInternal

First time NonlinearSystem::computeResidualInternal is called, this is what cr =

{<ThreadedElementLoop<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >> = {<ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >> = {_vptr.ThreadedElementLoopBase = 0x7ffff761db10 <vtable for ComputeResidualThread+16>, _mesh = @0xd37270, _tid = 0, _subdomain = 0, _old_subdomain = 0},
    _system = @0xd62f80, _fe_problem = @0xd60110}, _sys = @0xd62f80, _kernel_type = Moose::KT_ALL, _num_cached = 0, _integrated_bcs = @0xd63a78, _dg_kernels = @0xd634f8, _interface_kernels = @0xd63658,
  _kernels = @0xd631f0, _time_kernels = @0xd637b8, _non_time_kernels = @0xd63918}

Second time, cr =
{<ThreadedElementLoop<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >> = {<ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >> = {_vptr.ThreadedElementLoopBase = 0x7ffff761db10 <vtable for ComputeResidualThread+16>, _mesh = @0xd37270, _tid = 0, _subdomain = 0, _old_subdomain = 0},
    _system = @0xd62f80, _fe_problem = @0xd60110}, _sys = @0xd62f80, _kernel_type = Moose::KT_ALL, _num_cached = 0, _integrated_bcs = @0xd63a78, _dg_kernels = @0xd634f8, _interface_kernels = @0xd63658,
  _kernels = @0xd631f0, _time_kernels = @0xd637b8, _non_time_kernels = @0xd63918}

Boundary tests:

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc810, elem=0xcedcb0, side=0, bnd_id=3) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check passes)

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc810, elem=0xcf6d10, side=1, bnd_id=1) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check fails)

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc810, elem=0xcf6dc0, side=0, bnd_id=2) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check fails)

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc810, elem=0xcf9860, side=1, bnd_id=4) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check fails)

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc340, elem=0xcedcb0, side=0, bnd_id=3) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check passes)

I can see that we're moving from left to right most likely: bnd_id = 3 -> 'left', bnd_id = 2 -> 'master0_interface', bnd_id = 1 -> 'master1_interface', bnd_id = 4 -> 'right'

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc340, elem=0xcf6d10, side=1, bnd_id=1) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check fails)

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc340, elem=0xcf6dc0, side=0, bnd_id=2) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check fails)

Breakpoint 6, ComputeResidualThread::onBoundary (this=0x7fffffffc340, elem=0xcf9860, side=1, bnd_id=4) at /home/lindsayad/moose/framework/src/base/ComputeResidualThread.C:116 (check fails)

A couple things to note:

1) At least in this problem, the convection appears to be that the left side of an element is side 0, and the right side of an element is side 1.
2) Two ComputeResidualThread objects get created in between time steps. Why is that??????
3) Each element looked at in the above debug session has a unique address.
4) onInterface will only get called from ThreadedElementLoopBase if elem(neighbor) != NULL; this is why we only see InterfaceKernels getting called twice, for boundary IDs = 1,2, as opposed to four times

So some things here make sense: we calculate an initial residual from the initial conditions of the variables. Then we solve Newton-Raphson, and get delta_u. We add delta_u to u0 and then calculate the new residual. If we don't like it then we back track on our step, and choose a new delta_u and then calculate a new residual. Thus we can get multiple calls to Moose::compute_residual in that way. If we like the full Newton step, then we only calculate one new residual.

Alright, after all that haranguing, I believe I've convinced myself that my onInterface method only gets called once as it should and that the InterfaceKernel is only associated with one boundary ID corresponding to the boundary ID of the master sideset.

Zapdos: 59-81

This is a test.

This is test 2.

Sure enough, removing the term with _grad_u in the vacuum boundary conditions removed the numerical oscillations at the boundaries, supporting Hagelaar's statement:  "In order to circumvent possible numerical difficulties in accurately evaluating the density gradient in this term, we will now rewrite Eq. 9."

Now my problem is that I'm not sure if I trust the solution that I'm looking at.

NEWTON and PJFNK with LU and SMP=full are the exact same solves. Taking out the full preconditioning, e.g. going to just block diagonal, is almost exactly the same.

Next step for Zapdos: add back in all the physics. I'll probably go all the way to start. I always do :-)

NEWTON incapable of solving 2d.i with the Postprocessor, both with timestep\_end and timestep\_begin. PJFNK can with both timestep\_end and timestep\_begin, however, not with nonlinear. NEWTON solves very very quickly once the Postprocessor is replaced with a Dirichlet anode.

PJFNK is the only method that can solve my problem. FD and NEWTON cannot. Now, I can run PJFNK with SMP with timestep\_end or timestep\_begin, but I cannot run it with nonlinear. However, with PJFNK and FDP, I can run with nonlinear.

Search your feelings Alex. Recognize that the FDP preconditioner did improve convergence of PJFNK, e.g. it was actually able to converge with postprocessor on nonlinear iterations. With SMP, it cannot converge on nonlinear. And recall that PJFNK is a form of Newton's method! These aren't some dramatically different class of solves! So the reality is that a better jacobian did help my convergence.

Ok, after fixing the Jacobian in the Circuit dirichlet boundary condition, NEWTON + SMP + timestep\_end converges beautifully. NEWTON + FDP does not converge. NEWTON + SMP + nonlinear does not converge and neither does NEWTON + FDP + nonlinear. So convergers:

PJFNK + SMP + timestep\_end
PJFNK + FDP + timestep\_end
PJFNK + FDP + nonlinear
NEWTON + SMP + timestep\_end

Divergers:

PJFNK + SMP + nonlinear
NEWTON + SMP + nonlinear
NEWTON + FDP + timestep\_end

Very cool...if I implement a penalty dirichlet version of my circuit potential, then I can see the Jacobian functional dependence coming from the postprocessor.

What turns a Newton-Krylov method into a Jacobian-free Newton-Krylov method? The approximation made in eq 10 of the Knoll JFNK paper. Moving to the Jac-free method requires the introduction of the parameter $\epsilon$, which is the finite differencing parameter which MOOSE users can pass to PetSc. The Wikipedia article on preconditioning gives an excellent passage on why P should be chosen somewhere between the extremes of P = I and P = A, when we are solving the linear system Ax = b, where the left-preconditioned problem is: $P^{-1}\left(Ax - b\right) = 0$. It also describes preconditioned iterative methods. \textbf{Preconditioning only makes sense in the context of using an iterative method to solve the linear system of equations. There's no need for preconditioning if performing a direct solve of Ax = b.} Preconditioning is important for iterative linear methods because the rate of convergence increases as the condition number of the matrix decreases. A preconditioner of a matrix A is a matrix such that $P^{-1}A$ has a smaller condition number than A. Preconditioned iterative solvers typically outperform direct solvers, e.g. Gaussian elimination, for large, and especiallly for sparse, matrices.

Wrap up for the night: The default solve in PetSc is a Preconditioned Newton-Krylov solve. We pass PetSc the Jacobian matrix (which we will now think of as A). We then create the preconditioning matrix P from ILU factorization of the Jacobian matrix A. We then proceed to iteratively solve the linear problem $P^{-1}\left(Ax - b\right) = 0$. We perform the preconditioning because it increases the rate of convergence of the iterative method because the condition number of $P^{-1}A$ is less than A. We can turn away from the PetSc default and move to the MOOSE default which is PJFNK by simply making the approximation of eq. 10 from Knoll and specifying the finite differencing parameter $\epsilon$. That is the only difference between the MOOSE and PetSc defaults! However, that seemingly small change can make it such that one can get away with much lazier Jacobian specification in MOOSE. When using PJFNK, the only purpose filling out the Jacobian in MOOSE has, is its effect on the preconditioning matrix P. It does not affect A, because we never explicitly form A. Instead we only estimate its action on a vector through use of the finite differencing parameter $\epsilon$. Boom. Done for the night. If one has the perfect jacobian (or at least a very good one), direct methods will always be more accurate than iterative methods, they just may be too computationally expensive. If one does not have the perfect or very good jacobian, than a jacobian-free method may be more likely converge because the ``true'' Jacobian is being felt through finite differencing of the residual. Cool!

Round-trip Raleigh-Calgary (June 18th-23rd): $550
Raleigh-Calgary (June 18th): $231
Seattle-Raleigh (June 30th): $231.60

So Zapdos will take time steps with all the kernels active and all the non-linearity in transport and townsend coefficients active and ion boundary conditions active; however, the moment I try to add electron or electron energy boundary conditions, then I can't even solve for one nonlinear iteration.

texlive-extra-utils
texlive-font-utils
texlive-fonts-extra
texlive-bibtex-extra
texlive-formats-extratexlive-generic-recommended
texlive-humanities

A lot of the above items are dependencies for texlive-full -> may just have to reinstall texlive-full

Solved the problem immediately above.

Current problem to fix with Zapdos: diffusion-only electron boundary condition leads to residual equal to 10$^{139}$. Fixed that. Alright convergence with HagelaarIonDiffusionBC for electrons. And now convergence with addition of HagelaarIonAdvectionBC for electrons. so why doesn't combination work? Is it the addition of the secondaries from ion fluxes? It does indeed appear that in the case of the electrons, the numerical problems are coming from the secondary electron terms. Alright removing the n\_gamma terms resulted in convergence. I am also able to restore the diffusive portion of the ion flux for seconday electron calculations.

Trying to add back energy boundary condition now. Won't work with n\_gamma; however, it does work without n\_gamma.

Also, somehow I'm now able to run my integral postprocessor on every nonlinear iteration and achieve convergence. Don't have a good reason for that one.

Explanation for different bulk profiles with varying $gamma_{en}$ reflection:



0a) Observation: total current changes by only 17\% between gamma extremes
0b) Consequently the total potential drop across the plasma only varies by that amount
1) From 0b): Total potential drop more or less constant
2) Sound backing: Build up in electrons at anode because of gamma
3) From 0a): Total current more or less constant. In order for this to be true, the cathode characteristics must be the same because we aren't changing the cathodic boundary conditions.
4) From 2) electric field becomes large and negative, leading to a large potential drop in the anode.
5) From 1) and 4) for smaller gammas, the potential drop in the bulk will be smaller because more potential is being lost in the anode.
6) From 5) the bulk electric field is smaller for smaller gammas
7) With a smaller bulk electric field, in order to maintain the same total current between cases, the *density* of electrons (and with quasineutrality consequently the ions) must be higher!

2d\_rate\_coeffs failed at 6.02826e-9 seconds. No advection stabilization. Arp\_lin and em\_lin results looked fine. No oscillations. Potential also looked fine. Some small oscillations observed in mean\_en. Going to try some stabilization.

Currently EfieldArtDiff does not have position scaling included in it. Also might need to check its Jacobian.

Conductivity of aluminum: 3.77x10$^7/(\Omega m)$.

MOOSE model has power deposition at interface of .1e10 W m$^{-3}$. -> 1e9 W $m^{-3}$
What about VHF source? 420 W. Let's say 3 cm diameter, 3 cm length -> 2e7 W $m^{-3}$. Looks like 100 times less if we're assuming an even power deposition. Of course power deposition in most of the DC discharge is zero. I wonder what the average power deposition for that system is.

Let's just go with the electron density. Electron density at the interface: 3e22 $m^{-3}$. Mobility = 1.74e-7 m2 / (s*V). Resulting conductivity = 8.3e-4. That's really low.

Total Concentration of salt in a .5 M NaCl solution = 6e26. Thats much higher than the free electron concentration I'm calculating in Zapdos.

Diffusivity = 4.5e-9 m2/s.

Ok, using 1 mm thick water with radius = 1.125 cm, eps\_r = 80 and a conductivity of 1e-3 S/m, I get an impedance for the water layer of 4.84e-3 - 3.49j

How does this compare with the plasma? Look at one of Brandon's early presentations or thesis. 420 W air plasma impedance: 37.7 - 263.5j.

Conduction current and displacement current. Conduction current: motion of free charges.

Bound current has two components! 1) Vacuum permittivity times the time variation of the electric field and 2) the time variation of the polarizability of the material. For linear isotropic dielectric:

\begin{equation}
  \vec{P} = \epsilon_0\left(\epsilon_r-1\right)\vec{E}
\end{equation}

Thus in a vacuum, obviously the polarizability is zero.

Able to reproduce the behavior suggested by Lieberman that when the collisional frequency is sufficiently low, the plasma behaves inductively, e.g. the reactance is positive. However, as the collision frequency increases (e.g. as we increase the gas pressure) the reactance becomes negative.

PJFNK is able to take much longer time steps with Preconditioning type = SMP than NEWTON is. I think this clearly shows the important role that the missing potential Jacobian is playing. PJFNK is able to feel out that missing Jacobian piece by using finite differencing of the residuals.

Tried doing a JFNK run and the benefit of preconditioning is quite obvious. Even with the first solve, we were at over 200 linear iterations without achieving convergence. Going to a preconditioner matrix with the block components filled resulted in convergence after 8 linear iterations. Going to a preconditioner matrix with all elements filled (full = true) resulted in convergence after 4 linear iterations.

Presentation notes from Caroline:

1) Abrupt transition between slides 26 and 27 x
2) Explain why we ever added baking soda and why we explored its effects on the solution chemistry x
3) Further bridge between modelling and experiments and applications. x

Presentation notes from Steve:

1) Flip interfacial concentration and hydrophobic uptake slides x
2) Highlight questions x
3) Slide 28: give full name for PFOS and PFOA x
4) Slide 36-37: delete character in title x
5) Slide 59: Get rid of separation between our group and community x

My notes:

\begin{enumerate}
\item Give description of nitrite vs nitrate before moving into description of water chemistry. Overall just make that transition into water chemistry better x
\item Fix quality of some of the figures in the water chemistry portion of the talk x
\item Delete some unimportant slides x
\item Add slide numbers x
\item Send out slides and print slides

\end{enumerate}

Ok, by running mean\_en\_postprocessor.i I've confirmed that the physics implemented in the combination of SideTotFluxIntegral and CircuitDirichletPotential are correct. E.g. the residuals implemented are correct. However, convergence is inefficient and slow. I want my correct Jacobians!!!!

Comparing analytic vs. finite differenced Jacobian with mean\_en.i test file:

Analytic: Active time = 17.286 seconds, number of time steps = 126
FDP: unable to converge

This is what existed in ComputeUserObjectsThread before Daniel started working on it:

ComputeUserObjectsThread::ComputeUserObjectsThread(FEProblem & problem, SystemBase & sys, const NumericVector<Number>& in\_soln, std::vector<UserObjectWarehouse> & user\_objects, UserObjectWarehouse::GROUP group) :
    ThreadedElementLoop<ConstElemRange>(problem, sys),
    \_soln(in\_soln),
    \_user\_objects(user\_objects),
    \_group(group)
{
}

Two separate flights
189.23
250.60

One round-trip flight
401.33

Alright, I believe my suspicions are correct. If I use mesh size of 1, then the Jacobian is perfect. However, if I use a mesh size of 2, then it is wrong. I believe that this is because we have a disagreement between -_phi.size and the size of the _jacobian_information which is assigned a length of n_dofs(). n_dofs() divided by the number of variables equals the number of nodes and the number of total shape functions if everyone is linear lagrange. But number of nodes only equals _phi.size() if the mesh size is one.

First call to computeOffDiagJacobian: _v_dofs[0] = 1, _v_dofs[1] = 4, _u_dofs[0] = 0, _u_dofs[1] = 3

_v_dofs.size() and _u_dofs.size() = 2.
_shp_jacobian.size() = 9

After a whole bunch of calls, I believe we move to a different element because then _v_dofs[0] = 4. (Note that the kernel object, denoted by the memory location *this*, is the same as for the other element). So it appears that the Jacobian handling process is perfectly correct in the kernel object. So now we have to move back and look at how the Jacobian is being computed in the UserObject.

Contents of the _shp_jacobian vector: {0.015625, 0.0078125, 0, 0.010416666666666668, 0.0052083333333333339, 0, 0.015625, 0.0078125, 0}

Ok, so we are passing the argument cj which is of type ComputeFullJacobianThread to a function whose second argument is supposed to be a reference to a variable of type Body. How does that work for type checking?

Here's output from action system:

[DBG][ACT] Action Dependency Sets:
[DBG][ACT] (deprecated_block, finish_input_file_output, meta_action, no_action, setup_oversampling)
[DBG][ACT] (dynamic_object_registration)
[DBG][ACT] (common_output)
[DBG][ACT]	CommonOutputAction
[DBG][ACT] (set_global_params)
[DBG][ACT] (setup_recover_file_base)
[DBG][ACT]	SetupRecoverFileBaseAction
[DBG][ACT] (check_copy_nodal_vars)
[DBG][ACT]	u
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	v
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	w
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT] (setup_mesh)
[DBG][ACT]	Mesh
[DBG][ACT]	SetupMeshAction (setup_mesh, init_mesh)
[DBG][ACT] (add_partitioner)
[DBG][ACT] (init_mesh)
[DBG][ACT]	Mesh
[DBG][ACT]	SetupMeshAction (init_mesh, setup_mesh)
[DBG][ACT] (prepare_mesh)
[DBG][ACT]	Mesh
[DBG][ACT]	SetupMeshCompleteAction (prepare_mesh, execute_mesh_modifiers, setup_mesh_complete, uniform_refine_mesh)
[DBG][ACT] (add_mesh_modifier)
[DBG][ACT] (execute_mesh_modifiers)
[DBG][ACT]	Mesh
[DBG][ACT]	SetupMeshCompleteAction (execute_mesh_modifiers, prepare_mesh, setup_mesh_complete, uniform_refine_mesh)
[DBG][ACT] (add_mortar_interface)
[DBG][ACT] (uniform_refine_mesh)
[DBG][ACT]	Mesh
[DBG][ACT]	SetupMeshCompleteAction (uniform_refine_mesh, execute_mesh_modifiers, prepare_mesh, setup_mesh_complete)
[DBG][ACT] (setup_mesh_complete)
[DBG][ACT]	Mesh
[DBG][ACT]	SetupMeshCompleteAction (setup_mesh_complete, execute_mesh_modifiers, prepare_mesh, uniform_refine_mesh)
[DBG][ACT] (determine_system_type)
[DBG][ACT]	Executioner
[DBG][ACT]	DetermineSystemType
[DBG][ACT] (create_problem)
[DBG][ACT]	CreateProblemAction
[DBG][ACT] (setup_time_integrator)
[DBG][ACT] (setup_executioner)
[DBG][ACT]	Executioner
[DBG][ACT]	CreateExecutionerAction
[DBG][ACT] (setup_time_stepper)
[DBG][ACT] (setup_predictor)
[DBG][ACT] (setup_postprocessor_data)
[DBG][ACT] (init_displaced_problem)
[DBG][ACT]	Mesh
[DBG][ACT]	CreateDisplacedProblemAction
[DBG][ACT] (add_aux_variable, add_elemental_field_variable, add_variable)
[DBG][ACT]	u
[DBG][ACT]	AddVariableAction
[DBG][ACT]	v
[DBG][ACT]	AddVariableAction
[DBG][ACT]	w
[DBG][ACT]	AddVariableAction
[DBG][ACT] (setup_variable_complete)
[DBG][ACT] (setup_quadrature)
[DBG][ACT]	SetupQuadratureAction
[DBG][ACT] (add_function)
[DBG][ACT] (add_periodic_bc)
[DBG][ACT] (add_user_object)
[DBG][ACT]	example_uo
[DBG][ACT]	AddUserObjectAction
[DBG][ACT] (setup_function_complete)
[DBG][ACT] (setup_adaptivity)
[DBG][ACT] (set_adaptivity_options)
[DBG][ACT] (add_ic)
[DBG][ACT]	InitialCondition
[DBG][ACT]	AddICAction
[DBG][ACT]	InitialCondition
[DBG][ACT]	AddICAction
[DBG][ACT] (add_constraint, add_preconditioning, add_split)
[DBG][ACT]	smp
[DBG][ACT]	SetupPreconditionerAction
[DBG][ACT] (ready_to_init)
[DBG][ACT]	EmptyAction
[DBG][ACT] (setup_dampers)
[DBG][ACT]	SetupDampersAction
[DBG][ACT] (setup_residual_debug)
[DBG][ACT] (add_bounds_vectors)
[DBG][ACT] (add_multi_app)
[DBG][ACT] (add_transfer)
[DBG][ACT] (copy_nodal_aux_vars, copy_nodal_vars)
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	u
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	v
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	w
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT] (add_material)
[DBG][ACT] (setup_material_output)
[DBG][ACT]	MaterialOutputAction
[DBG][ACT] (init_problem)
[DBG][ACT]	InitProblemAction
[DBG][ACT] (setup_debug)
[DBG][ACT] (add_output)
[DBG][ACT] (add_postprocessor)
[DBG][ACT] (add_vector_postprocessor)
[DBG][ACT] (add_aux_kernel, add_aux_scalar_kernel, add_bc, add_damper, add_dg_kernel, add_dirac_kernel, add_indicator, add_interface_kernel, add_kernel, add_marker, add_nodal_kernel, add_scalar_kernel)
[DBG][ACT]	diff_u
[DBG][ACT]	AddKernelAction
[DBG][ACT]	diff_v
[DBG][ACT]	AddKernelAction
[DBG][ACT]	shape_w
[DBG][ACT]	AddKernelAction
[DBG][ACT]	time_u
[DBG][ACT]	AddKernelAction
[DBG][ACT]	time_v
[DBG][ACT]	AddKernelAction
[DBG][ACT] (add_control)
[DBG][ACT] (check_output)
[DBG][ACT]	CheckOutputAction
[DBG][ACT] (check_integrity)
[DBG][ACT]	CheckIntegrityAction


[DBG][ACT] Executing actions:
[DBG][ACT] TASK (           common_output) TYPE (              CommonOutputAction) NAME (                )
[DBG][ACT] TASK ( setup_recover_file_base) TYPE (      SetupRecoverFileBaseAction) NAME (                )
[DBG][ACT] TASK (   check_copy_nodal_vars) TYPE (             CopyNodalVarsAction) NAME (               u)
[DBG][ACT] TASK (   check_copy_nodal_vars) TYPE (             CopyNodalVarsAction) NAME (               v)
[DBG][ACT] TASK (   check_copy_nodal_vars) TYPE (             CopyNodalVarsAction) NAME (               w)
[DBG][ACT] TASK (              setup_mesh) TYPE (                 SetupMeshAction) NAME (            Mesh)
[DBG][ACT] TASK (               init_mesh) TYPE (                 SetupMeshAction) NAME (            Mesh)
[DBG][ACT] TASK (            prepare_mesh) TYPE (         SetupMeshCompleteAction) NAME (            Mesh)
[DBG][ACT] TASK (  execute_mesh_modifiers) TYPE (         SetupMeshCompleteAction) NAME (            Mesh)
[DBG][ACT] TASK (     uniform_refine_mesh) TYPE (         SetupMeshCompleteAction) NAME (            Mesh)
[DBG][ACT] TASK (     setup_mesh_complete) TYPE (         SetupMeshCompleteAction) NAME (            Mesh)
[DBG][ACT] TASK (   determine_system_type) TYPE (             DetermineSystemType) NAME (     Executioner)
[DBG][ACT] TASK (          create_problem) TYPE (             CreateProblemAction) NAME (                )
[DBG][ACT] TASK (       setup_executioner) TYPE (         CreateExecutionerAction) NAME (     Executioner)
[DBG][ACT] TASK (  init_displaced_problem) TYPE (    CreateDisplacedProblemAction) NAME (            Mesh)
[DBG][ACT] TASK (            add_variable) TYPE (               AddVariableAction) NAME (               u)
[DBG][ACT] TASK (            add_variable) TYPE (               AddVariableAction) NAME (               v)
[DBG][ACT] TASK (            add_variable) TYPE (               AddVariableAction) NAME (               w)
[DBG][ACT] TASK (        setup_quadrature) TYPE (           SetupQuadratureAction) NAME (                )
[DBG][ACT] TASK (         add_user_object) TYPE (             AddUserObjectAction) NAME (      example_uo)
[DBG][ACT] TASK (                  add_ic) TYPE (                     AddICAction) NAME (InitialCondition)
[DBG][ACT] TASK (                  add_ic) TYPE (                     AddICAction) NAME (InitialCondition)
[DBG][ACT] TASK (     add_preconditioning) TYPE (       SetupPreconditionerAction) NAME (             smp)
[DBG][ACT] TASK (           ready_to_init) TYPE (                     EmptyAction) NAME (                )
[DBG][ACT] TASK (           setup_dampers) TYPE (              SetupDampersAction) NAME (                )
[DBG][ACT] TASK (     copy_nodal_aux_vars) TYPE (             CopyNodalVarsAction) NAME (                )
[DBG][ACT] TASK (         copy_nodal_vars) TYPE (             CopyNodalVarsAction) NAME (               u)
[DBG][ACT] TASK (         copy_nodal_vars) TYPE (             CopyNodalVarsAction) NAME (               v)
[DBG][ACT] TASK (         copy_nodal_vars) TYPE (             CopyNodalVarsAction) NAME (               w)
[DBG][ACT] TASK (   setup_material_output) TYPE (            MaterialOutputAction) NAME (                )
[DBG][ACT] TASK (            init_problem) TYPE (               InitProblemAction) NAME (                )
[DBG][ACT] TASK (              add_output) TYPE (                 AddOutputAction) NAME (         console)
[DBG][ACT] TASK (              add_kernel) TYPE (                 AddKernelAction) NAME (          diff_u)
[DBG][ACT] TASK (              add_kernel) TYPE (                 AddKernelAction) NAME (          diff_v)
[DBG][ACT] TASK (              add_kernel) TYPE (                 AddKernelAction) NAME (         shape_w)
[DBG][ACT] TASK (              add_kernel) TYPE (                 AddKernelAction) NAME (          time_u)
[DBG][ACT] TASK (              add_kernel) TYPE (                 AddKernelAction) NAME (          time_v)
[DBG][ACT] TASK (            check_output) TYPE (               CheckOutputAction) NAME (                )
[DBG][ACT] TASK (         check_integrity) TYPE (            CheckIntegrityAction) NAME (                )

Note the type of _sub_Kee:

  /// jacobian contributions
  std::vector<std::vector<DenseMatrix<Number> > > _sub_Kee;

  It's a vector of a vector of DenseMatrices; you access individual DenseMatrices by specifying the ivar and the jvar. So how are the dimensions of the DenseMatrices determined? And when?

  Exploring ability to ignore floating point exceptions. After adding ignore functions in NeumannCircuitVoltageMoles_KV.C, here is the backtrace:

    * frame #0: 0x00007ffff78818b6 libzapdos-dbg.so.0`_ZNK7libMesh10TypeVectorIdEmlIdEENS_9boostcopy11enable_if_cIXsr12ScalarTraitsIT_EE5valueENS0_INS_12CompareTypesIdS5_E9supertypeEEEE4typeES5_(this=0x00007fffffffc1d8, factor=+Inf) + 38 at type_vector.h:731
    frame #1: 0x00007ffff7873a06 libzapdos-dbg.so.0`_ZN7libMeshmlIddEENS_9boostcopy11enable_if_cIXsr12ScalarTraitsIT0_EE5valueENS_10TypeVectorINS_12CompareTypesIT_S3_E9supertypeEEEE4typeES3_RKNS4_IS6_EE(factor=+Inf, v=0x00007fffffffc1d8) + 38 at type_vector.h:746
    frame #2: 0x00007ffff79ceac0 libzapdos-dbg.so.0`NeumannCircuitVoltageMoles_KV::computeQpResidual(this=0x0000000000d79540) + 528 at NeumannCircuitVoltageMoles_KV.C:99
    frame #3: 0x00007ffff6bcbd05 libmoose-dbg.so.0`IntegratedBC::computeResidual(this=0x0000000000d79540) + 325 at IntegratedBC.C:117
    frame #4: 0x00007ffff673f333 libmoose-dbg.so.0`ComputeResidualThread::onBoundary(this=0x00007fffffffc718, elem=0x0000000000af4160, side=0, bnd_id=3) + 419 at ComputeResidualThread.C:132
    frame #5: 0x00007ffff660b164 libmoose-dbg.so.0`ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >::operator(this=0x00007fffffffc718, range=0x0000000000b4a530, bypass_threading=false)(libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> const&, bool) + 740 at ThreadedElementLoopBase.h:173
    frame #6: 0x00007ffff6863f4b libmoose-dbg.so.0`void libMesh::Threads::parallel_reduce<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*>, ComputeResidualThread>(range=0x0000000000b4a530, body=0x00007fffffffc718) + 155 at threads_tbb.h:160
    frame #7: 0x00007ffff684f1ad libmoose-dbg.so.0`NonlinearSystem::computeResidualInternal(this=0x0000000000b11cf0, type=KT_ALL) + 541 at NonlinearSystem.C:1271
    frame #8: 0x00007ffff684ed32 libmoose-dbg.so.0`NonlinearSystem::computeResidual(this=0x0000000000b11cf0, residual=0x0000000000be3370, type=KT_ALL) + 434 at NonlinearSystem.C:765
    frame #9: 0x00007ffff65fd20c libmoose-dbg.so.0`FEProblem::computeResidualType(this=0x0000000000b0f7d0, soln=0x0000000000b0d350, residual=0x0000000000be3370, type=KT_ALL) + 652 at FEProblem.C:3355
    frame #10: 0x00007ffff65fcb3d libmoose-dbg.so.0`FEProblem::computeResidual(this=0x0000000000b0f7d0, (null)=0x0000000000b0ce60, soln=0x0000000000b0d350, residual=0x0000000000be3370) + 93 at FEProblem.C:3295
    frame #11: 0x00007ffff684a4d1 libmoose-dbg.so.0`NonlinearSystem::solve(this=0x0000000000b11cf0) + 225 at NonlinearSystem.C:244
    frame #12: 0x00007ffff65fbc81 libmoose-dbg.so.0`FEProblem::solve(this=0x0000000000b0f7d0) + 129 at FEProblem.C:3074
    frame #13: 0x00007ffff6be24f4 libmoose-dbg.so.0`TimeStepper::step(this=0x0000000000b48f70) + 36 at TimeStepper.C:188
    frame #14: 0x00007ffff6bfba31 libmoose-dbg.so.0`Transient::solveStep(this=0x0000000000b518d0, input_dt=-1) + 625 at Transient.C:415
    frame #15: 0x00007ffff6bfb650 libmoose-dbg.so.0`Transient::takeStep(this=0x0000000000b518d0, input_dt=-1) + 480 at Transient.C:342
    frame #16: 0x00007ffff6bfb213 libmoose-dbg.so.0`Transient::execute(this=0x0000000000b518d0) + 211 at Transient.C:255
    frame #17: 0x00007ffff674f4e9 libmoose-dbg.so.0`MooseApp::executeExecutioner(this=0x00000000006c1f80) + 265 at MooseApp.C:428
    frame #18: 0x00007ffff675044c libmoose-dbg.so.0`MooseApp::run(this=0x00000000006c1f80) + 140 at MooseApp.C:578
    frame #19: 0x0000000000411473 zapdos-dbg`main(argc=3, argv=0x00007fffffffd1e8) + 211 at main.C:23
    frame #20: 0x00007fffee640830 libc.so.6`__libc_start_main(main=(zapdos-dbg`main at main.C:12), argc=3, argv=0x00007fffffffd1e8, init=<unavailable>, fini=<unavailable>, rtld_fini=<unavailable>, stack_end=0x00007fffffffd1d8) + 240 at libc-start.c:291
    frame #21: 0x00000000004112c9 zapdos-dbg`_start + 41

    Now, I'm trying to add a similar ignore function in libmesh type_vector.h.

      break libmesh_handleFPE
  run ...
  bt

  #0  0x00007ffff78818dd in libMesh::TypeVector<double>::operator*<double> (this=0x7fffffffc198, factor=inf)
    at /home/lindsayad/projects_devel/moose/scripts/../libmesh/installed/include/libmesh/type_vector.h:734
#1  0x00007ffff7873a06 in libMesh::operator*<double, double> (factor=inf, v=...)
    at /home/lindsayad/projects_devel/moose/scripts/../libmesh/installed/include/libmesh/type_vector.h:749
#2  0x00007ffff79ceaf0 in NeumannCircuitVoltageMoles_KV::computeQpResidual (this=0xd7a160)
    at /home/lindsayad/projects_devel/zapdos/src/bcs/NeumannCircuitVoltageMoles_KV.C:99
#3  0x00007ffff6bcbda5 in IntegratedBC::computeResidual (this=0xd7a160)
    at /home/lindsayad/projects_devel/moose/framework/src/bcs/IntegratedBC.C:117
#4  0x00007ffff673f3d3 in ComputeResidualThread::onBoundary (this=0x7fffffffc6d8, elem=0xaf4280, side=0, bnd_id=3)
    at /home/lindsayad/projects_devel/moose/framework/src/base/ComputeResidualThread.C:132
#5  0x00007ffff660b204 in ThreadedElementLoopBase<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*> >::operator() (this=0x7fffffffc6d8, range=..., bypass_threading=false)
    at /home/lindsayad/projects_devel/moose/framework/include/base/ThreadedElementLoopBase.h:173
#6  0x00007ffff6863feb in libMesh::Threads::parallel_reduce<libMesh::StoredRange<libMesh::MeshBase::const_element_iterator, libMesh::Elem const*>, ComputeResidualThread> (range=..., body=...)
    at /home/lindsayad/projects_devel/moose/scripts/../libmesh/installed/include/libmesh/threads_tbb.h:160
#7  0x00007ffff684f24d in NonlinearSystem::computeResidualInternal (this=0xb11ec0, type=Moose::KT_ALL)
    at /home/lindsayad/projects_devel/moose/framework/src/base/NonlinearSystem.C:1271
#8  0x00007ffff684edd2 in NonlinearSystem::computeResidual (this=0xb11ec0, residual=..., type=Moose::KT_ALL)
    at /home/lindsayad/projects_devel/moose/framework/src/base/NonlinearSystem.C:765
#9  0x00007ffff65fd2ac in FEProblem::computeResidualType (this=0xb0f9a0, soln=..., residual=..., type=Moose::KT_ALL)
    at /home/lindsayad/projects_devel/moose/framework/src/base/FEProblem.C:3355
#10 0x00007ffff65fcbdd in FEProblem::computeResidual (this=0xb0f9a0, soln=..., residual=...)
    at /home/lindsayad/projects_devel/moose/framework/src/base/FEProblem.C:3295
#11 0x00007ffff684a571 in NonlinearSystem::solve (this=0xb11ec0)
    at /home/lindsayad/projects_devel/moose/framework/src/base/NonlinearSystem.C:244
#12 0x00007ffff65fbd21 in FEProblem::solve (this=0xb0f9a0) at /home/lindsayad/projects_devel/moose/framework/src/base/FEProblem.C:3074
#13 0x00007ffff6be2594 in TimeStepper::step (this=0xb49100)
    at /home/lindsayad/projects_devel/moose/framework/src/timesteppers/TimeStepper.C:188
#14 0x00007ffff6bfbad1 in Transient::solveStep (this=0xb51b00, input_dt=-1)
    at /home/lindsayad/projects_devel/moose/framework/src/executioners/Transient.C:415
#15 0x00007ffff6bfb6f0 in Transient::takeStep (this=0xb51b00, input_dt=-1)
    at /home/lindsayad/projects_devel/moose/framework/src/executioners/Transient.C:342
#16 0x00007ffff6bfb2b3 in Transient::execute (this=0xb51b00)
    at /home/lindsayad/projects_devel/moose/framework/src/executioners/Transient.C:255
#17 0x00007ffff674f589 in MooseApp::executeExecutioner (this=0x6c2030)
    at /home/lindsayad/projects_devel/moose/framework/src/base/MooseApp.C:428
#18 0x00007ffff67504ec in MooseApp::run (this=0x6c2030) at /home/lindsayad/projects_devel/moose/framework/src/base/MooseApp.C:578
    #19 0x0000000000411473 in main (argc=3, argv=0x7fffffffd1a8) at /home/lindsayad/projects_devel/zapdos/src/main.C:23

    Through experimentation I discovered that making the time stamp of the Makefile newer than the files it's compiled does not cause a recompilation.

    agf93f

    Timing:
    \begin{itemize}
    \item direct LU solve time = 22 seconds (obviousy with solve\_type=Newton)
    \item solve\_type=PFJNK -pc\_type=bjacobi -sub\_pc\_type=lu solve time = 85 seconds; also much
      worse convergence without -sub\_pc\_type=lu
    \item solve\_type=PFJNK -pc\_type=asm -sub\_pc\_type=lu solve time = 85 seconds; also much
      worse convergence without -sub\_pc\_type=lu
      \item solve\_type=PFJNK -pc\_type=ilu very slow convergence (doesn't converge at all without
        -sub\_pc\_type=lu
    \item solve\_type=NEWTON, -pc\_type=asm -sub\_pc\_type=lu solve time = 22
      seconds; pretty close to the direct solve. I think this is proof of just
      how superior the NEWTON method is when the Jacobian is correct

Ok with whatever changes happened in MOOSE: -pc_type=asm, -sub_pc_type = lu, solve time = 238 seconds so way longer. What the hell changed? It has to be something with Jacobians right?

Gamg doesn't work well.

Alright dope, I can get NEWTON iterative (with pc_type = asm and sub_pc_type = lu) and NEWTON direct to solve at the same speed!!!! And thus I have a problem that can be naturally scaled in parallel!!!

### 11/21/16

PenaltyCircuitPotential uses the non-local jacobian user-object feature, and it
works!!! Important to note: changing the penalty from 1 to 1000 did not change
the solution at all.

# 12/7/17

First two time-steps in Zapdos with NEWTON, SMP, no line search:
```
Time Step  1, time = 1e-11
                dt = 1e-11
    |residual|_2 of individual variables:
                  potential: 88500
                  em:        0.676043
                  emliq:     860.265
                  Arp:       0.427591
                  mean_en:   40.3082
                  OHm:       4274.12

 0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 1.603806e-05
    |residual|_2 of individual variables:
                  potential: 0.0165603
                  em:        0.666791
                  emliq:     727.807
                  Arp:       0.638707
                  mean_en:   89.8353
                  OHm:       0.00068818

 1 Nonlinear |R| = 7.333310e+02
      0 Linear |R| = 7.333310e+02
      1 Linear |R| = 1.505504e-10
    |residual|_2 of individual variables:
                  potential: 0.00994833
                  em:        1.23299
                  emliq:     222208
                  Arp:       1.2318
                  mean_en:   30.774
                  OHm:       3.61972e-05

 2 Nonlinear |R| = 2.222078e+05
      0 Linear |R| = 2.222078e+05
      1 Linear |R| = 9.200163e-09
    |residual|_2 of individual variables:
                  potential: 0.0219212
                  em:        0.141632
                  emliq:     85817.1
                  Arp:       0.140243
                  mean_en:   8.35208
                  OHm:       0.000120854

 3 Nonlinear |R| = 8.581712e+04
      0 Linear |R| = 8.581712e+04
      1 Linear |R| = 2.035883e-09
    |residual|_2 of individual variables:
                  potential: 0.000988502
                  em:        0.00216072
                  emliq:     31529.1
                  Arp:       0.00136111
                  mean_en:   1.56318
                  OHm:       0.000110193

 4 Nonlinear |R| = 3.152906e+04
      0 Linear |R| = 3.152906e+04
      1 Linear |R| = 8.344147e-10
    |residual|_2 of individual variables:
                  potential: 9.59574e-06
                  em:        1.72857e-05
                  emliq:     11580.6
                  Arp:       1.94736e-07
                  mean_en:   0.00789034
                  OHm:       1.75771e-05

 5 Nonlinear |R| = 1.158060e+04
      0 Linear |R| = 1.158060e+04
      1 Linear |R| = 2.749675e-10
    |residual|_2 of individual variables:
                  potential: 5.51199e-08
                  em:        1.60183e-08
                  emliq:     4247.17
                  Arp:       6.46679e-11
                  mean_en:   7.45486e-07
                  OHm:       3.04818e-06

 6 Nonlinear |R| = 4.247168e+03
      0 Linear |R| = 4.247168e+03
      1 Linear |R| = 1.181364e-10
    |residual|_2 of individual variables:
                  potential: 2.15625e-08
                  em:        3.94818e-11
                  emliq:     1549.58
                  Arp:       5.06658e-14
                  mean_en:   2.88025e-11
                  OHm:       6.13661e-07

 7 Nonlinear |R| = 1.549583e+03
      0 Linear |R| = 1.549583e+03
      1 Linear |R| = 8.500159e-11
    |residual|_2 of individual variables:
                  potential: 8.44478e-09
                  em:        3.76462e-12
                  emliq:     555.328
                  Arp:       2.51319e-14
                  mean_en:   2.52167e-11
                  OHm:       1.51653e-07

 8 Nonlinear |R| = 5.553279e+02
      0 Linear |R| = 5.553279e+02
      1 Linear |R| = 2.197129e-11
    |residual|_2 of individual variables:
                  potential: 3.25608e-09
                  em:        3.43769e-12
                  emliq:     185.212
                  Arp:       2.77896e-14
                  mean_en:   1.98979e-11
                  OHm:       4.51212e-08

 9 Nonlinear |R| = 1.852125e+02
      0 Linear |R| = 1.852125e+02
      1 Linear |R| = 7.165906e-12
    |residual|_2 of individual variables:
                  potential: 1.19019e-09
                  em:        2.4609e-12
                  emliq:     46.6456
                  Arp:       2.15093e-14
                  mean_en:   2.02119e-11
                  OHm:       1.45765e-08

10 Nonlinear |R| = 4.664555e+01
      0 Linear |R| = 4.664555e+01
      1 Linear |R| = 1.462664e-12
    |residual|_2 of individual variables:
                  potential: 3.73756e-10
                  em:        3.31522e-12
                  emliq:     7.91837
                  Arp:       2.14571e-14
                  mean_en:   1.86068e-11
                  OHm:       4.33845e-09

11 Nonlinear |R| = 7.918367e+00
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 11
 Solve Converged!

Outlier Variable Residual Norms:
  emliq: 7.918367e+00

```
Without line search:
```
Time Step  1, time = 1e-11
                dt = 1e-11
    |residual|_2 of individual variables:
                  potential: 88500
                  em:        0.676043
                  emliq:     860.265
                  Arp:       0.427591
                  mean_en:   40.3082
                  OHm:       4274.12

 0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 1.603806e-05
      Line search: Using full step: fnorm 8.860733466940e+04 gnorm 7.333310287514e+02
    |residual|_2 of individual variables:
                  potential: 0.0165603
                  em:        0.666791
                  emliq:     727.807
                  Arp:       0.638707
                  mean_en:   89.8353
                  OHm:       0.00068818

 1 Nonlinear |R| = 7.333310e+02
      0 Linear |R| = 7.333310e+02
      1 Linear |R| = 1.505504e-10
      Line search: gnorm after quadratic fit 6.580493332020e+02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.0148226
                  em:        0.616675
                  emliq:     653.017
                  Arp:       0.591948
                  mean_en:   81.2205
                  OHm:       0.000619362

 2 Nonlinear |R| = 6.580493e+02
      0 Linear |R| = 6.580493e+02
      1 Linear |R| = 6.516030e-11
      Line search: Using full step: fnorm 6.580493332020e+02 gnorm 4.450273271312e+01
    |residual|_2 of individual variables:
                  potential: 0.00975053
                  em:        1.05859
                  emliq:     33.9477
                  Arp:       1.05742
                  mean_en:   28.7369
                  OHm:       4.92565e-08

 3 Nonlinear |R| = 4.450273e+01
      0 Linear |R| = 4.450273e+01
      1 Linear |R| = 2.436905e-10
      Line search: gnorm after quadratic fit 3.772670233073e+01
      Line search: Quadratically determined step, lambda=1.9710179627999808e-01
    |residual|_2 of individual variables:
                  potential: 0.00735178
                  em:        0.846267
                  emliq:     29.9376
                  Arp:       0.845299
                  mean_en:   22.9263
                  OHm:       4.02052e-08

 4 Nonlinear |R| = 3.772670e+01
      0 Linear |R| = 3.772670e+01
      1 Linear |R| = 4.003852e-11
      Line search: gnorm after quadratic fit 2.158673544182e+01
      Line search: Quadratically determined step, lambda=4.1724765619856824e-01
    |residual|_2 of individual variables:
                  potential: 0.00303388
                  em:        0.482791
                  emliq:     16.8887
                  Arp:       0.482103
                  mean_en:   13.4273
                  OHm:       2.56175e-08

 5 Nonlinear |R| = 2.158674e+01
      0 Linear |R| = 2.158674e+01
      1 Linear |R| = 3.748465e-11
      Line search: gnorm after quadratic fit 1.978660328162e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00271731
                  em:        0.434333
                  emliq:     15.6571
                  Arp:       0.433712
                  mean_en:   12.0825
                  OHm:       2.31348e-08

 6 Nonlinear |R| = 1.978660e+01
      0 Linear |R| = 1.978660e+01
      1 Linear |R| = 2.897005e-11
      Line search: gnorm after quadratic fit 1.898009983393e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00243508
                  em:        0.390754
                  emliq:     15.5475
                  Arp:       0.390196
                  mean_en:   10.8727
                  OHm:       2.08924e-08

 7 Nonlinear |R| = 1.898010e+01
      0 Linear |R| = 1.898010e+01
      1 Linear |R| = 5.518920e-11
      Line search: gnorm after quadratic fit 2.693537431437e+01
      Line search: Cubically determined step, current gnorm 1.873214841536e+01 lambda=6.3811928280707850e-02
    |residual|_2 of individual variables:
                  potential: 0.0022763
                  em:        0.365771
                  emliq:     15.7171
                  Arp:       0.36525
                  mean_en:   10.1784
                  OHm:       1.95839e-08

 8 Nonlinear |R| = 1.873215e+01
      0 Linear |R| = 1.873215e+01
      1 Linear |R| = 2.179061e-11
      Line search: gnorm after quadratic fit 1.801665678989e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00204139
                  em:        0.329092
                  emliq:     15.5076
                  Arp:       0.328623
                  mean_en:   9.15964
                  OHm:       1.76788e-08

 9 Nonlinear |R| = 1.801666e+01
      0 Linear |R| = 1.801666e+01
      1 Linear |R| = 1.851680e-11
      Line search: gnorm after quadratic fit 1.544497745388e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00183143
                  em:        0.296104
                  emliq:     13.0547
                  Arp:       0.295679
                  mean_en:   8.24294
                  OHm:       1.59613e-08

10 Nonlinear |R| = 1.544498e+01
      0 Linear |R| = 1.544498e+01
      1 Linear |R| = 1.386760e-11
      Line search: gnorm after quadratic fit 1.323920075399e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00164364
                  em:        0.26643
                  emliq:     10.9593
                  Arp:       0.266046
                  mean_en:   7.41805
                  OHm:       1.44053e-08

11 Nonlinear |R| = 1.323920e+01
      0 Linear |R| = 1.323920e+01
      1 Linear |R| = 3.503345e-11
      Line search: gnorm after quadratic fit 1.150146744120e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00147555
                  em:        0.239737
                  emliq:     9.35967
                  Arp:       0.239388
                  mean_en:   6.67575
                  OHm:       1.3005e-08

12 Nonlinear |R| = 1.150147e+01
      0 Linear |R| = 1.150147e+01
      1 Linear |R| = 2.987686e-11
      Line search: gnorm after quadratic fit 1.011708605815e+01
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00132501
                  em:        0.215722
                  emliq:     8.13444
                  Arp:       0.215407
                  mean_en:   6.00778
                  OHm:       1.17317e-08

13 Nonlinear |R| = 1.011709e+01
      0 Linear |R| = 1.011709e+01
      1 Linear |R| = 1.048194e-11
      Line search: gnorm after quadratic fit 8.976975383599e+00
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
    |residual|_2 of individual variables:
                  potential: 0.00119011
                  em:        0.194116
                  emliq:     7.16091
                  Arp:       0.193832
                  mean_en:   5.40668
                  OHm:       1.05847e-08

14 Nonlinear |R| = 8.976975e+00
      0 Linear |R| = 8.976975e+00
      1 Linear |R| = 4.367762e-11
      Line search: gnorm after quadratic fit 7.459271195424e+00
      Line search: Quadratically determined step, lambda=1.5343584416609685e-01
    |residual|_2 of individual variables:
                  potential: 0.00100297
                  em:        0.164269
                  emliq:     5.8858
                  Arp:       0.164026
                  mean_en:   4.57648
                  OHm:       9.01078e-09

15 Nonlinear |R| = 7.459271e+00
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 15
 Solve Converged!
```
The non-line-search performed better. If I change to FDP then the problem
diverges because of `DIVERGED_FNORM_NAN` and I ultimately hit the minimum time
step. Why? If I change to PJFNK, then the solve diverges for the first few
attempted time steps but eventually starts converging.

First converged time step:
```
Time Step  1, time = 1.6e-12
                dt = 1.6e-12
    |residual|_2 of individual variables:
                  potential: 88500
                  em:        0.676043
                  emliq:     860.265
                  Arp:       0.427591
                  mean_en:   40.3082
                  OHm:       4274.12

 0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 5.865874e+04
      2 Linear |R| = 3.389564e+04
      3 Linear |R| = 9.814945e+03
      4 Linear |R| = 5.744699e+03
      5 Linear |R| = 4.782629e+03
      6 Linear |R| = 3.458547e+03
      7 Linear |R| = 1.495322e+03
      8 Linear |R| = 1.031347e+03
      9 Linear |R| = 1.347923e+02
     10 Linear |R| = 1.043910e+02
     11 Linear |R| = 1.600047e+01
     12 Linear |R| = 1.598497e+01
     13 Linear |R| = 1.433009e+01
     14 Linear |R| = 5.494862e+00
     15 Linear |R| = 3.020107e+00
     16 Linear |R| = 2.121015e+00
     17 Linear |R| = 1.497558e+00
     18 Linear |R| = 1.327567e+00
     19 Linear |R| = 8.809741e-01
    |residual|_2 of individual variables:
                  potential: 95385.3
                  em:        0.261039
                  emliq:     255.853
                  Arp:       0.00278899
                  mean_en:   13.751
                  OHm:       4606.52

 1 Nonlinear |R| = 9.549684e+04
      0 Linear |R| = 9.549684e+04
      1 Linear |R| = 2.835309e+00
      2 Linear |R| = 7.803261e-03
    |residual|_2 of individual variables:
                  potential: 0.00962467
                  em:        0.0595431
                  emliq:     1.37893e+07
                  Arp:       6.1888e-06
                  mean_en:   0.726753
                  OHm:       0.0406746

 2 Nonlinear |R| = 1.378932e+07
      0 Linear |R| = 1.378932e+07
      1 Linear |R| = 2.673771e+00
    |residual|_2 of individual variables:
                  potential: 0.006444
                  em:        0.000753558
                  emliq:     5.15246e+06
                  Arp:       9.26303e-07
                  mean_en:   0.0524176
                  OHm:       0.00502752

 3 Nonlinear |R| = 5.152460e+06
      0 Linear |R| = 5.152460e+06
      1 Linear |R| = 1.267888e+00
    |residual|_2 of individual variables:
                  potential: 0.00586926
                  em:        7.10828e-05
                  emliq:     1.93305e+06
                  Arp:       7.03933e-07
                  mean_en:   0.000965679
                  OHm:       0.00539508

 4 Nonlinear |R| = 1.933045e+06
      0 Linear |R| = 1.933045e+06
      1 Linear |R| = 7.564543e-02
    |residual|_2 of individual variables:
                  potential: 0.00585886
                  em:        0.000332403
                  emliq:     710413
                  Arp:       7.67945e-07
                  mean_en:   0.000929986
                  OHm:       0.00147378

 5 Nonlinear |R| = 7.104130e+05
      0 Linear |R| = 7.104130e+05
      1 Linear |R| = 1.990726e-01
    |residual|_2 of individual variables:
                  potential: 0.00585891
                  em:        0.000335579
                  emliq:     260902
                  Arp:       7.78929e-07
                  mean_en:   0.000945384
                  OHm:       0.000226901

 6 Nonlinear |R| = 2.609015e+05
      0 Linear |R| = 2.609015e+05
      1 Linear |R| = 1.633487e-02
    |residual|_2 of individual variables:
                  potential: 0.00585897
                  em:        0.000339401
                  emliq:     95672
                  Arp:       7.9337e-07
                  mean_en:   0.000962707
                  OHm:       3.61212e-05

 7 Nonlinear |R| = 9.567197e+04
      0 Linear |R| = 9.567197e+04
      1 Linear |R| = 1.554469e-01
    |residual|_2 of individual variables:
                  potential: 0.00585904
                  em:        0.000343368
                  emliq:     34948
                  Arp:       8.08355e-07
                  mean_en:   0.000980682
                  OHm:       6.10755e-06

 8 Nonlinear |R| = 3.494800e+04
      0 Linear |R| = 3.494800e+04
      1 Linear |R| = 7.005794e-03
    |residual|_2 of individual variables:
                  potential: 0.0058591
                  em:        0.00034738
                  emliq:     12619
                  Arp:       8.23511e-07
                  mean_en:   0.000998862
                  OHm:       1.1582e-06

 9 Nonlinear |R| = 1.261899e+04
      0 Linear |R| = 1.261899e+04
      1 Linear |R| = 6.214382e-03
    |residual|_2 of individual variables:
                  potential: 0.00585917
                  em:        0.000351289
                  emliq:     4378.18
                  Arp:       8.38271e-07
                  mean_en:   0.00101657
                  OHm:       2.63668e-07

10 Nonlinear |R| = 4.378181e+03
      0 Linear |R| = 4.378181e+03
      1 Linear |R| = 5.976326e-03
    |residual|_2 of individual variables:
                  potential: 0.00585923
                  em:        0.000354855
                  emliq:     1322.82
                  Arp:       8.51736e-07
                  mean_en:   0.00103272
                  OHm:       7.43118e-08

11 Nonlinear |R| = 1.322820e+03
      0 Linear |R| = 1.322820e+03
      1 Linear |R| = 5.981777e-03
    |residual|_2 of individual variables:
                  potential: 0.00585927
                  em:        0.000357785
                  emliq:     281.993
                  Arp:       8.628e-07
                  mean_en:   0.00104599
                  OHm:       2.48914e-08

12 Nonlinear |R| = 2.819935e+02
      0 Linear |R| = 2.819935e+02
      1 Linear |R| = 5.980525e-03
      2 Linear |R| = 5.964502e-03
      3 Linear |R| = 5.963972e-03
      4 Linear |R| = 1.603342e-03
    |residual|_2 of individual variables:
                  potential: 0.000387588
                  em:        0.00664627
                  emliq:     97.3819
                  Arp:       7.06456e-06
                  mean_en:   12.3928
                  OHm:       0.0149275

13 Nonlinear |R| = 9.816733e+01
      0 Linear |R| = 9.816733e+01
      1 Linear |R| = 2.151425e-01
      2 Linear |R| = 2.150790e-01
      3 Linear |R| = 9.199626e-04
    |residual|_2 of individual variables:
                  potential: 2.31695e-05
                  em:        0.000497679
                  emliq:     21.3141
                  Arp:       9.7245e-07
                  mean_en:   0.13054
                  OHm:       0.0048875

14 Nonlinear |R| = 2.131448e+01
      0 Linear |R| = 2.131448e+01
      1 Linear |R| = 5.724278e-02
      2 Linear |R| = 3.228554e-02
      3 Linear |R| = 5.159937e-06
    |residual|_2 of individual variables:
                  potential: 8.8018e-08
                  em:        3.54163e-05
                  emliq:     7.05421
                  Arp:       4.25439e-07
                  mean_en:   0.000220283
                  OHm:       0.00137588

15 Nonlinear |R| = 7.054209e+00
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 15
 Solve Converged!

Outlier Variable Residual Norms:
  emliq: 7.054209e+00
```

Now going to SMP, PJFNK:
```
Time Step  1, time = 1e-11
                dt = 1e-11
    |residual|_2 of individual variables:
                  potential: 88500
                  em:        0.676043
                  emliq:     860.265
                  Arp:       0.427591
                  mean_en:   40.3082
                  OHm:       4274.12

 0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 3.294046e+03
      2 Linear |R| = 3.378421e+01
      3 Linear |R| = 8.642622e-02
    |residual|_2 of individual variables:
                  potential: 3679.28
                  em:        0.666326
                  emliq:     721.653
                  Arp:       0.638305
                  mean_en:   89.8144
                  OHm:       177.691

 1 Nonlinear |R| = 3.754668e+03
      0 Linear |R| = 3.754668e+03
      1 Linear |R| = 6.528383e-02
      2 Linear |R| = 8.526220e-06
    |residual|_2 of individual variables:
                  potential: 0.0100226
                  em:        1.23261
                  emliq:     625551
                  Arp:       1.23145
                  mean_en:   30.7695
                  OHm:       0.0190874

 2 Nonlinear |R| = 6.255510e+05
      0 Linear |R| = 6.255510e+05
      1 Linear |R| = 3.232792e+00
    |residual|_2 of individual variables:
                  potential: 0.0220638
                  em:        0.141517
                  emliq:     232001
                  Arp:       0.140158
                  mean_en:   8.38418
                  OHm:       0.000146046

 3 Nonlinear |R| = 2.320009e+05
      0 Linear |R| = 2.320009e+05
      1 Linear |R| = 9.071890e-01
    |residual|_2 of individual variables:
                  potential: 0.00104775
                  em:        0.0021738
                  emliq:     85297.1
                  Arp:       0.00135924
                  mean_en:   1.56687
                  OHm:       0.000136127

 4 Nonlinear |R| = 8.529715e+04
      0 Linear |R| = 8.529715e+04
      1 Linear |R| = 1.260183e-01
    |residual|_2 of individual variables:
                  potential: 9.58477e-06
                  em:        1.77893e-05
                  emliq:     31352.2
                  Arp:       1.90448e-07
                  mean_en:   0.00892927
                  OHm:       2.15712e-05

 5 Nonlinear |R| = 3.135218e+04
      0 Linear |R| = 3.135218e+04
      1 Linear |R| = 5.472914e-02
    |residual|_2 of individual variables:
                  potential: 6.12417e-08
                  em:        1.71207e-08
                  emliq:     11513.1
                  Arp:       6.17478e-11
                  mean_en:   9.17619e-07
                  OHm:       3.69371e-06

 6 Nonlinear |R| = 1.151308e+04
      0 Linear |R| = 1.151308e+04
      1 Linear |R| = 1.627385e-02
    |residual|_2 of individual variables:
                  potential: 2.39604e-08
                  em:        4.26009e-11
                  emliq:     4213.38
                  Arp:       5.32318e-14
                  mean_en:   2.56559e-11
                  OHm:       7.28331e-07

 7 Nonlinear |R| = 4.213376e+03
      0 Linear |R| = 4.213376e+03
      1 Linear |R| = 7.392047e-03
    |residual|_2 of individual variables:
                  potential: 9.38937e-09
                  em:        3.75505e-12
                  emliq:     1520.08
                  Arp:       4.39268e-14
                  mean_en:   2.15846e-11
                  OHm:       1.75672e-07

 8 Nonlinear |R| = 1.520081e+03
      0 Linear |R| = 1.520081e+03
      1 Linear |R| = 1.867816e-03
    |residual|_2 of individual variables:
                  potential: 3.63113e-09
                  em:        2.50022e-12
                  emliq:     511.967
                  Arp:       2.52437e-14
                  mean_en:   1.53428e-11
                  OHm:       5.1381e-08

 9 Nonlinear |R| = 5.119666e+02
      0 Linear |R| = 5.119666e+02
      1 Linear |R| = 1.082642e-03
    |residual|_2 of individual variables:
                  potential: 1.33979e-09
                  em:        2.7278e-12
                  emliq:     125.714
                  Arp:       3.08331e-14
                  mean_en:   2.64657e-11
                  OHm:       1.65598e-08

10 Nonlinear |R| = 1.257135e+02
      0 Linear |R| = 1.257135e+02
      1 Linear |R| = 7.161949e-04
    |residual|_2 of individual variables:
                  potential: 4.31829e-10
                  em:        2.93255e-12
                  emliq:     12.6985
                  Arp:       2.71921e-14
                  mean_en:   2.02786e-11
                  OHm:       5.03973e-09

11 Nonlinear |R| = 1.269852e+01
      0 Linear |R| = 1.269852e+01
      1 Linear |R| = 4.118541e-04
      2 Linear |R| = 4.107763e-07
    |residual|_2 of individual variables:
                  potential: 1.62114e-08
                  em:        5.59193e-10
                  emliq:     2.14082
                  Arp:       1.45729e-12
                  mean_en:   3.68463e-09
                  OHm:       2.38543e-09

12 Nonlinear |R| = 2.140816e+00
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 12
 Solve Converged!

Outlier Variable Residual Norms:
  emliq: 2.140816e+00
```

Using 1e-6:
```
0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 3.044054e-05
    |residual|_2 of individual variables:
                  potential: 9189.99
                  em:        1.1939e+11
                  emliq:     2.45558e+15
                  Arp:       0.716018
                  mean_en:   1.53715e+10
                  OHm:       1.87793e+15

 1 Nonlinear |R| = 3.091357e+15
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1
```
Using 1e-50
```
0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 3.044054e-05
    |residual|_2 of individual variables:
                  potential: 9189.99
                  em:        1.1939e+11
                  emliq:     2.45558e+15
                  Arp:       0.716018
                  mean_en:   1.53715e+10
                  OHm:       1.87793e+15

 1 Nonlinear |R| = 3.091357e+15
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1
 Solve Did NOT Converge!
```
Using default:
```
 0 Nonlinear |R| = 8.860733e+04
      0 Linear |R| = 8.860733e+04
      1 Linear |R| = 3.044054e-05
    |residual|_2 of individual variables:
                  potential: 9189.99
                  em:        1.1939e+11
                  emliq:     2.45558e+15
                  Arp:       0.716018
                  mean_en:   1.53715e+10
                  OHm:       1.87793e+15

 1 Nonlinear |R| = 3.091357e+15
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1
 Solve Did NOT Converge!
```

So changing the error at least in the way I'm specifying is doing absolutely
nothing. I need to try on a simpler problem where the NEWTON solve is definitely
going to converge with FDP and use `-options_left`.

# 3/5/18

Ok, IMO Zapdos clearly fails because of poor conditioning, witness:

Time Step 55, time = 1.13219e-11
                dt = 1.88707e-12
    |residual|_2 of individual variables:
                  potential: 4.65914e-11
                  em:        0.414222
                  emliq:     8.43238
                  Arp:       0.413861
                  mean_en:   20.3745
                  OHm:       0.00329452

 0 Nonlinear |R| = 2.205826e+01
    0 KSP unpreconditioned resid norm 2.205825692078e+01 true resid norm 2.205825692078e+01 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 2.205825692078e+01 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.384125056288e+00 true resid norm 2.384125056288e+00 ||r(i)||/||b|| 1.080831121358e-01
    1 KSP Residual norm 2.384125056288e+00 % max 1.005822073456e+00 min 1.005822073456e+00 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 8.446026833059e-01 true resid norm 8.693782758933e-01 ||r(i)||/||b|| 3.941282754189e-02
    2 KSP Residual norm 8.446026833059e-01 % max 3.410404257731e+00 min 1.001227270872e+00 max/min 3.406223898356e+00
    3 KSP unpreconditioned resid norm 8.414424096053e-01 true resid norm 8.435728879437e-01 ||r(i)||/||b|| 3.824295323848e-02
    3 KSP Residual norm 8.414424096053e-01 % max 1.880665401624e+01 min 9.315117418670e-01 max/min 2.018939018261e+01
    4 KSP unpreconditioned resid norm 4.302092907871e-01 true resid norm 2.551250542306e-01 ||r(i)||/||b|| 1.156596621151e-02
    4 KSP Residual norm 4.302092907871e-01 % max 2.693589616887e+01 min 8.816270674693e-01 max/min 3.055248320153e+01
    5 KSP unpreconditioned resid norm 2.494447433940e-02 true resid norm 3.267145574787e-01 ||r(i)||/||b|| 1.481144038951e-02
    5 KSP Residual norm 2.494447433940e-02 % max 2.727483225745e+01 min 8.140348617902e-01 max/min 3.350572996035e+01
    6 KSP unpreconditioned resid norm 3.505727283562e-03 true resid norm 3.772101451243e-01 ||r(i)||/||b|| 1.710063249689e-02
    6 KSP Residual norm 3.505727283562e-03 % max 3.968085580019e+01 min 5.216696675460e-01 max/min 7.606510071183e+01
    7 KSP unpreconditioned resid norm 1.426741614424e-04 true resid norm 4.490751218516e-01 ||r(i)||/||b|| 2.035859512673e-02
    7 KSP Residual norm 1.426741614424e-04 % max 4.038976013438e+01 min 5.070742457323e-01 max/min 7.965255674946e+01
      Line search: Using full step: fnorm 2.205825692078e+01 gnorm 5.409576886305e+00
    |residual|_2 of individual variables:
                  potential: 1.44398e-06
                  em:        0.00272412
                  emliq:     5.2894
                  Arp:       0.0025772
                  mean_en:   1.13389
                  OHm:       0.00478482

 1 Nonlinear |R| = 5.409577e+00
    0 KSP unpreconditioned resid norm 5.409576886305e+00 true resid norm 5.409576886305e+00 ||r(i)||/||b|| 1.000000000000e+00
    0 KSP Residual norm 5.409576886305e+00 % max 1.000000000000e+00 min 1.000000000000e+00 max/min 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.572102579936e-01 true resid norm 3.572102579936e-01 ||r(i)||/||b|| 6.603293852757e-02
    1 KSP Residual norm 3.572102579936e-01 % max 9.944054516878e-01 min 9.944054516878e-01 max/min 1.000000000000e+00
    2 KSP unpreconditioned resid norm 9.409709998170e-02 true resid norm 7.432426507433e-02 ||r(i)||/||b|| 1.373938602527e-02
    2 KSP Residual norm 9.409709998170e-02 % max 1.088529687336e+01 min 9.894736606428e-01 max/min 1.100109816596e+01
    3 KSP unpreconditioned resid norm 9.401884704269e-02 true resid norm 7.155770235289e-02 ||r(i)||/||b|| 1.322796659643e-02
    3 KSP Residual norm 9.401884704269e-02 % max 1.270318472876e+01 min 6.352144178087e-01 max/min 1.999826259074e+01
    4 KSP unpreconditioned resid norm 2.898141975363e-02 true resid norm 5.044378289867e-02 ||r(i)||/||b|| 9.324903584673e-03
    4 KSP Residual norm 2.898141975363e-02 % max 1.431702185898e+02 min 6.342005276083e-01 max/min 2.257491319499e+02
    5 KSP unpreconditioned resid norm 2.094076772138e-02 true resid norm 8.821712062322e-02 ||r(i)||/||b|| 1.630758236315e-02
    5 KSP Residual norm 2.094076772138e-02 % max 1.696704410372e+02 min 2.905480631076e-01 max/min 5.839668632530e+02
    6 KSP unpreconditioned resid norm 2.075085667714e-02 true resid norm 8.870406130169e-02 ||r(i)||/||b|| 1.639759692228e-02
    6 KSP Residual norm 2.075085667714e-02 % max 3.094260510923e+02 min 1.183064914499e-01 max/min 2.615461309859e+03
    7 KSP unpreconditioned resid norm 2.032077074790e-02 true resid norm 2.242633852539e-01 ||r(i)||/||b|| 4.145673311745e-02
    7 KSP Residual norm 2.032077074790e-02 % max 3.125720663423e+02 min 2.614556971125e-02 max/min 1.195506809736e+04
    8 KSP unpreconditioned resid norm 1.888566083704e-02 true resid norm 1.643064951885e+00 ||r(i)||/||b|| 3.037326183577e-01
    8 KSP Residual norm 1.888566083704e-02 % max 3.219988857588e+02 min 3.930828783756e-03 max/min 8.191628368284e+04
    9 KSP unpreconditioned resid norm 4.941432315092e-04 true resid norm 1.320663938108e+01 ||r(i)||/||b|| 2.441344241638e+00
    9 KSP Residual norm 4.941432315092e-04 % max 3.220011525510e+02 min 1.508550904326e-03 max/min 2.134506377130e+05
   10 KSP unpreconditioned resid norm 1.242091202663e-04 true resid norm 1.334597250807e+01 ||r(i)||/||b|| 2.467100993028e+00
   10 KSP Residual norm 1.242091202663e-04 % max 3.391809032844e+02 min 1.398430318278e-03 max/min 2.425440144218e+05
   11 KSP unpreconditioned resid norm 2.546092229864e-06 true resid norm 1.333058199913e+01 ||r(i)||/||b|| 2.464255944467e+00
   11 KSP Residual norm 2.546092229864e-06 % max 3.392316538481e+02 min 1.395205631358e-03 max/min 2.431409723583e+05
      Line search: gnorm after quadratic fit 5.760190942666e+00
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.460453562498e+00 lambda=1.9085662507565725e-02
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.419720425842e+00 lambda=4.0563605904697974e-03
      Line search: Cubic step no good, shrinking lambda, current gnorm 5.411743407567e+00 lambda=8.7876011150836105e-04
      Line search: unable to find good step length! After 3 tries
      Line search: fnorm=5.4095768863054809e+00, gnorm=5.4117434075673048e+00, ynorm=4.2119679961219519e-01, minlambda=1.0000000000000000e-03, lambda=8.7876011150836105e-04, initial slope=-1.3263729134982935e+01
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 1
 Solve Did NOT Converge!
