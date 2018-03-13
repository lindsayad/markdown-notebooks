### Non-linear residuals

Time Step  1, time = 0.5
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.666666668455e+01 true resid norm 5.666666668455e+01 ||r(i)||/||b|| 7.071067814097e-01
    2 KSP unpreconditioned resid norm 2.847481230244e-03 true resid norm 2.847448264429e-03 ||r(i)||/||b|| 3.553147017981e-05
    3 KSP unpreconditioned resid norm 2.697165536176e-11 true resid norm 5.364626739391e+10 ||r(i)||/||b|| 6.694171669337e+08
  Linear solve converged due to CONVERGED_RTOL iterations 3
 1 Nonlinear |R| = 3.766654e+01
    0 KSP unpreconditioned resid norm 3.766654077455e+01 true resid norm 3.766654077455e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.404259796105e-05 true resid norm 1.256713623141e+00 ||r(i)||/||b|| 3.336419000256e-02
  Linear solve converged due to CONVERGED_RTOL iterations 1
 2 Nonlinear |R| = 7.274320e-05
    0 KSP unpreconditioned resid norm 7.274320370973e-05 true resid norm 7.274320370973e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.178996051004e-05 true resid norm 2.178996051004e-05 ||r(i)||/||b|| 2.995463410848e-01
    2 KSP unpreconditioned resid norm 2.173400257865e-06 true resid norm 2.173398902384e-06 ||r(i)||/||b|| 2.987769016961e-02
    3 KSP unpreconditioned resid norm 1.441207565656e-09 true resid norm 1.440427799971e-09 ||r(i)||/||b|| 1.980154470125e-05
    4 KSP unpreconditioned resid norm 1.032848664134e-14 true resid norm 6.250904983865e-03 ||r(i)||/||b|| 8.593112022956e+01
  Linear solve converged due to CONVERGED_RTOL iterations 4
 3 Nonlinear |R| = 5.327259e-05
    0 KSP unpreconditioned resid norm 5.327258649792e-05 true resid norm 5.327258649792e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.766749845922e-05 true resid norm 3.766749845922e-05 ||r(i)||/||b|| 7.070709521622e-01
    2 KSP unpreconditioned resid norm 2.892464896149e-12 true resid norm 1.981023925990e-08 ||r(i)||/||b|| 3.718655421523e-04
  Linear solve converged due to CONVERGED_RTOL iterations 2
 4 Nonlinear |R| = 1.265893e-09
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 4
SNES Object: 1 MPI processes
  type: newtonls
  maximum iterations=15, maximum function evaluations=10000
  tolerances: relative=1e-10, absolute=1e-50, solution=1e-50
  total number of linear solver iterations=10
  total number of function evaluations=25
  norm schedule ALWAYS
  SNESLineSearch Object: 1 MPI processes
    type: basic
    maxstep=1.000000e+08, minlambda=1.000000e-12
    tolerances: relative=1.000000e-08, absolute=1.000000e-15, lambda=1.000000e-08
    maximum iterations=40
  KSP Object: 1 MPI processes
    type: gmres
      restart=30, using Classical (unmodified) Gram-Schmidt Orthogonalization with no iterative refinement
      happy breakdown tolerance 1e-30
    maximum iterations=200, initial guess is zero
    tolerances:  relative=1e-05, absolute=1e-50, divergence=10000.
    right preconditioning
    using UNPRECONDITIONED norm type for convergence test
  PC Object: 1 MPI processes
    type: none
    linear system matrix = precond matrix:
    Mat Object: 1 MPI processes
      type: mffd
      rows=8, cols=8
        Matrix-free approximation:
          err=1.49012e-08 (relative error in function evaluation)
          Using wp compute h routine
              Does not compute normU
 Solve Converged!

Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.652520137400e+01 true resid norm 5.652520137400e+01 ||r(i)||/||b|| 7.053415270498e-01
    2 KSP unpreconditioned resid norm 2.739718390433e+01 true resid norm 2.739712747531e+01 ||r(i)||/||b|| 3.418710815793e-01
    3 KSP unpreconditioned resid norm 5.424589545514e-02 true resid norm 5.426219993677e-02 ||r(i)||/||b|| 6.771029918363e-04
    4 KSP unpreconditioned resid norm 3.842311738285e-06 true resid norm 8.611287607979e-02 ||r(i)||/||b|| 1.074547034532e-03
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 1.043647e-04
    0 KSP unpreconditioned resid norm 1.043646505227e-04 true resid norm 1.043646505227e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.327805099082e-05 true resid norm 2.327805099082e-05 ||r(i)||/||b|| 2.230453594607e-01
    2 KSP unpreconditioned resid norm 2.731267022459e-06 true resid norm 2.731265413112e-06 ||r(i)||/||b|| 2.617040731160e-02
    3 KSP unpreconditioned resid norm 1.406938574106e-08 true resid norm 1.407334017569e-08 ||r(i)||/||b|| 1.348477679483e-04
    4 KSP unpreconditioned resid norm 1.765016369080e-13 true resid norm 2.364487015022e-07 ||r(i)||/||b|| 2.265601430351e-03
  Linear solve converged due to CONVERGED_RTOL iterations 4
 2 Nonlinear |R| = 2.576720e-10
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 2



### Every residual

Time Step  1, time = 0.5
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.666666668455e+01 true resid norm 5.666666668455e+01 ||r(i)||/||b|| 7.071067814097e-01
    2 KSP unpreconditioned resid norm 2.847481230244e-03 true resid norm 2.847448264429e-03 ||r(i)||/||b|| 3.553147017981e-05
    3 KSP unpreconditioned resid norm 2.697165536176e-11 true resid norm 5.364626739391e+10 ||r(i)||/||b|| 6.694171669337e+08
  Linear solve converged due to CONVERGED_RTOL iterations 3
 1 Nonlinear |R| = 3.766654e+01
    0 KSP unpreconditioned resid norm 3.766654077455e+01 true resid norm 3.766654077455e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.679986601811e-01 true resid norm 8.679618919429e-01 ||r(i)||/||b|| 2.304331308622e-02
    2 KSP unpreconditioned resid norm 2.217283365177e-06 true resid norm 4.067284126384e+02 ||r(i)||/||b|| 1.079813554085e+01
  Linear solve converged due to CONVERGED_RTOL iterations 2
 2 Nonlinear |R| = 2.239722e+00
    0 KSP unpreconditioned resid norm 2.239721657565e+00 true resid norm 2.239721657565e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.939342180510e+00 true resid norm 1.939343289026e+00 ||r(i)||/||b|| 8.658858490187e-01
    2 KSP unpreconditioned resid norm 1.279777191163e-05 true resid norm 1.304051874070e+00 ||r(i)||/||b|| 5.822383641578e-01
  Linear solve converged due to CONVERGED_RTOL iterations 2
 3 Nonlinear |R| = 2.914668e+00
    0 KSP unpreconditioned resid norm 2.914667707621e+00 true resid norm 2.914667707621e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.063758326266e+00 true resid norm 2.063760384187e+00 ||r(i)||/||b|| 7.080602631961e-01
    2 KSP unpreconditioned resid norm 1.092310519418e-05 true resid norm 2.309586139845e+00 ||r(i)||/||b|| 7.924011830940e-01
  Linear solve converged due to CONVERGED_RTOL iterations 2
 4 Nonlinear |R| = 1.499432e+00
    0 KSP unpreconditioned resid norm 1.499432455250e+00 true resid norm 1.499432455250e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.254139802509e+00 true resid norm 1.254140623904e+00 ||r(i)||/||b|| 8.364102160876e-01
    2 KSP unpreconditioned resid norm 1.332344028453e-05 true resid norm 1.654617301228e+00 ||r(i)||/||b|| 1.103495722955e+00
  Linear solve converged due to CONVERGED_RTOL iterations 2
 5 Nonlinear |R| = 4.099467e+00
    0 KSP unpreconditioned resid norm 4.099466712500e+00 true resid norm 4.099466712500e+00 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.163040290573e-01 true resid norm 7.162999864106e-01 ||r(i)||/||b|| 1.747300409164e-01
    2 KSP unpreconditioned resid norm 2.012124598479e-06 true resid norm 2.066972415816e+00 ||r(i)||/||b|| 5.042051956449e-01
  Linear solve converged due to CONVERGED_RTOL iterations 2
 6 Nonlinear |R| = 7.101350e-01
    0 KSP unpreconditioned resid norm 7.101350010808e-01 true resid norm 7.101350010808e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.147035853772e-01 true resid norm 5.147040751850e-01 ||r(i)||/||b|| 7.247975024491e-01
    2 KSP unpreconditioned resid norm 2.222905261982e-06 true resid norm 1.375223640778e+00 ||r(i)||/||b|| 1.936566481986e+00
  Linear solve converged due to CONVERGED_RTOL iterations 2
 7 Nonlinear |R| = 6.797180e-01
    0 KSP unpreconditioned resid norm 6.797180353188e-01 true resid norm 6.797180353188e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.808873149925e-01 true resid norm 4.808877953318e-01 ||r(i)||/||b|| 7.074812942198e-01
    2 KSP unpreconditioned resid norm 1.231327641547e-06 true resid norm 3.516607764423e-02 ||r(i)||/||b|| 5.173627271451e-02
  Linear solve converged due to CONVERGED_RTOL iterations 2
 8 Nonlinear |R| = 8.341716e-02
    0 KSP unpreconditioned resid norm 8.341715700058e-02 true resid norm 8.341715700058e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.194942898544e-02 true resid norm 8.194209384766e-02 ||r(i)||/||b|| 9.823170291825e-01
    2 KSP unpreconditioned resid norm 3.811962258143e-07 true resid norm 1.267748021567e-01 ||r(i)||/||b|| 1.519768914635e+00
  Linear solve converged due to CONVERGED_RTOL iterations 2
 9 Nonlinear |R| = 4.621135e-02
    0 KSP unpreconditioned resid norm 4.621135398979e-02 true resid norm 4.621135398979e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.333010951424e-02 true resid norm 3.333014152225e-02 ||r(i)||/||b|| 7.212543811119e-01
    2 KSP unpreconditioned resid norm 2.804464948244e-07 true resid norm 1.906697895000e-02 ||r(i)||/||b|| 4.126037716665e-01
  Linear solve converged due to CONVERGED_RTOL iterations 2
10 Nonlinear |R| = 1.001390e-02
    0 KSP unpreconditioned resid norm 1.001390165718e-02 true resid norm 1.001390165718e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.267786594352e-03 true resid norm 7.267793482895e-03 ||r(i)||/||b|| 7.257704071505e-01
    2 KSP unpreconditioned resid norm 1.744832581313e-07 true resid norm 2.823858096424e-03 ||r(i)||/||b|| 2.819937915408e-01
    3 KSP unpreconditioned resid norm 5.055232910182e-10 true resid norm 2.829966469682e-03 ||r(i)||/||b|| 2.826037808803e-01
  Linear solve converged due to CONVERGED_RTOL iterations 3
11 Nonlinear |R| = 1.257813e-03
    0 KSP unpreconditioned resid norm 1.257813196203e-03 true resid norm 1.257813196203e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.707695060217e-04 true resid norm 8.707704136443e-04 ||r(i)||/||b|| 6.922891382224e-01
    2 KSP unpreconditioned resid norm 5.143018585512e-10 true resid norm 1.211424671368e-04 ||r(i)||/||b|| 9.631197025319e-02
  Linear solve converged due to CONVERGED_RTOL iterations 2
12 Nonlinear |R| = 7.274340e-05
    0 KSP unpreconditioned resid norm 7.274340010252e-05 true resid norm 7.274340010252e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.125059067480e-05 true resid norm 5.125064229569e-05 ||r(i)||/||b|| 7.045400988058e-01
    2 KSP unpreconditioned resid norm 3.191021258140e-10 true resid norm 7.715363922994e-07 ||r(i)||/||b|| 1.060627343803e-02
  Linear solve converged due to CONVERGED_RTOL iterations 2
13 Nonlinear |R| = 1.106037e-06
    0 KSP unpreconditioned resid norm 1.106036587357e-06 true resid norm 1.106036587357e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.812593144716e-07 true resid norm 7.812600973448e-07 ||r(i)||/||b|| 7.063600845358e-01
    2 KSP unpreconditioned resid norm 2.149694114375e-10 true resid norm 2.117081769378e-08 ||r(i)||/||b|| 1.914115494532e-02
    3 KSP unpreconditioned resid norm 5.936432151859e-13 true resid norm 2.118498816276e-08 ||r(i)||/||b|| 1.915396688042e-02
  Linear solve converged due to CONVERGED_RTOL iterations 3
14 Nonlinear |R| = 2.156498e-08
    0 KSP unpreconditioned resid norm 2.156497620288e-08 true resid norm 2.156497620288e-08 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.524765410978e-08 true resid norm 1.524766935886e-08 ||r(i)||/||b|| 7.070570918052e-01
    2 KSP unpreconditioned resid norm 2.117805209279e-14 true resid norm 4.248472382985e-10 ||r(i)||/||b|| 1.970079791889e-02
  Linear solve converged due to CONVERGED_RTOL iterations 2
15 Nonlinear |R| = 4.245493e-10
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 15
SNES Object: 1 MPI processes
  type: newtonls
  maximum iterations=15, maximum function evaluations=10000
  tolerances: relative=1e-10, absolute=1e-50, solution=1e-50
  total number of linear solver iterations=33
  total number of function evaluations=82
  norm schedule ALWAYS
  SNESLineSearch Object: 1 MPI processes
    type: basic
    maxstep=1.000000e+08, minlambda=1.000000e-12
    tolerances: relative=1.000000e-08, absolute=1.000000e-15, lambda=1.000000e-08
    maximum iterations=40
  KSP Object: 1 MPI processes
    type: gmres
      restart=30, using Classical (unmodified) Gram-Schmidt Orthogonalization with no iterative refinement
      happy breakdown tolerance 1e-30
    maximum iterations=200, initial guess is zero
    tolerances:  relative=1e-05, absolute=1e-50, divergence=10000.
    right preconditioning
    using UNPRECONDITIONED norm type for convergence test
  PC Object: 1 MPI processes
    type: none
    linear system matrix = precond matrix:
    Mat Object: 1 MPI processes
      type: mffd
      rows=8, cols=8
        Matrix-free approximation:
          err=1.49012e-08 (relative error in function evaluation)
          Using wp compute h routine
              Does not compute normU
 Solve Converged!

Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.013876853447e+01 true resid norm 8.013876853447e+01 ||r(i)||/||b|| 9.999999999999e-01
    2 KSP unpreconditioned resid norm 8.013875528035e+01 true resid norm 1.365802293196e+02 ||r(i)||/||b|| 1.704296582257e+00
    3 KSP unpreconditioned resid norm 1.244041587808e+01 true resid norm 3.255468594125e+08 ||r(i)||/||b|| 4.062289268551e+06
    4 KSP unpreconditioned resid norm 1.134921504469e+00 true resid norm 3.335154947984e+08 ||r(i)||/||b|| 4.161724729460e+06
    5 KSP unpreconditioned resid norm 1.539553341317e-02 true resid norm 3.335824099331e+08 ||r(i)||/||b|| 4.162559720263e+06
    6 KSP unpreconditioned resid norm 1.090040178304e-02 true resid norm 3.335822396110e+08 ||r(i)||/||b|| 4.162557594924e+06
    7 KSP unpreconditioned resid norm 8.905094369648e-03 true resid norm 3.335825823457e+08 ||r(i)||/||b|| 4.162561871689e+06
    8 KSP unpreconditioned resid norm 7.713691258597e-03 true resid norm 3.335822426078e+08 ||r(i)||/||b|| 4.162557632319e+06
    9 KSP unpreconditioned resid norm 6.900480433749e-03 true resid norm 3.335825847501e+08 ||r(i)||/||b|| 4.162561901692e+06
   10 KSP unpreconditioned resid norm 6.299791644668e-03 true resid norm 3.335822427506e+08 ||r(i)||/||b|| 4.162557634101e+06
   11 KSP unpreconditioned resid norm 5.832933781807e-03 true resid norm 3.335825832316e+08 ||r(i)||/||b|| 4.162561882743e+06
   12 KSP unpreconditioned resid norm 5.456462806870e-03 true resid norm 3.335822416293e+08 ||r(i)||/||b|| 4.162557620108e+06
   13 KSP unpreconditioned resid norm 5.144638920350e-03 true resid norm 3.335825814849e+08 ||r(i)||/||b|| 4.162561860947e+06
   14 KSP unpreconditioned resid norm 4.880773927839e-03 true resid norm 3.335822440328e+08 ||r(i)||/||b|| 4.162557650100e+06
   15 KSP unpreconditioned resid norm 4.653775259025e-03 true resid norm 3.335825835867e+08 ||r(i)||/||b|| 4.162561887175e+06
   16 KSP unpreconditioned resid norm 4.455738359030e-03 true resid norm 3.335822422735e+08 ||r(i)||/||b|| 4.162557628148e+06
   17 KSP unpreconditioned resid norm 4.281025658881e-03 true resid norm 3.335825816598e+08 ||r(i)||/||b|| 4.162561863130e+06
   18 KSP unpreconditioned resid norm 4.125358515474e-03 true resid norm 3.335822470411e+08 ||r(i)||/||b|| 4.162557687639e+06
   19 KSP unpreconditioned resid norm 3.985537681872e-03 true resid norm 3.335825826041e+08 ||r(i)||/||b|| 4.162561874914e+06
   20 KSP unpreconditioned resid norm 3.859022111809e-03 true resid norm 3.335822451338e+08 ||r(i)||/||b|| 4.162557663838e+06
   21 KSP unpreconditioned resid norm 3.743846834544e-03 true resid norm 3.335825843492e+08 ||r(i)||/||b|| 4.162561896689e+06
   22 KSP unpreconditioned resid norm 3.638396083584e-03 true resid norm 3.335822449484e+08 ||r(i)||/||b|| 4.162557661525e+06
   23 KSP unpreconditioned resid norm 3.541388905011e-03 true resid norm 3.335825841170e+08 ||r(i)||/||b|| 4.162561893791e+06
   24 KSP unpreconditioned resid norm 3.451742585281e-03 true resid norm 3.335822433496e+08 ||r(i)||/||b|| 4.162557641575e+06
   25 KSP unpreconditioned resid norm 3.368582330542e-03 true resid norm 3.335825824841e+08 ||r(i)||/||b|| 4.162561873416e+06
   26 KSP unpreconditioned resid norm 3.291151853915e-03 true resid norm 3.335822444329e+08 ||r(i)||/||b|| 4.162557655093e+06
   27 KSP unpreconditioned resid norm 3.218831019678e-03 true resid norm 3.335825835424e+08 ||r(i)||/||b|| 4.162561886622e+06
   28 KSP unpreconditioned resid norm 3.151073432974e-03 true resid norm 3.335822447319e+08 ||r(i)||/||b|| 4.162557658824e+06
   29 KSP unpreconditioned resid norm 3.087425812572e-03 true resid norm 3.335825838221e+08 ||r(i)||/||b|| 4.162561890112e+06
   30 KSP unpreconditioned resid norm 3.335822444249e+08 true resid norm 3.335820804780e+08 ||r(i)||/||b|| 4.162555609206e+06
  Linear solve did not converge due to DIVERGED_DTOL iterations 30
 1 Nonlinear |R| = 3.344477e+08

### Every residual but at the end

Now with the contact state and contact force evaluated on every residual call **but** now at the end of the compute residual thread (although before nodal bcs); this is with JFNK:

Time Step  1, time = 0.5
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.666666668455e+01 true resid norm 5.666666668455e+01 ||r(i)||/||b|| 7.071067814097e-01
    2 KSP unpreconditioned resid norm 2.847481230244e-03 true resid norm 2.847448264429e-03 ||r(i)||/||b|| 3.553147017981e-05
    3 KSP unpreconditioned resid norm 2.833351389804e-03 true resid norm 2.847448262641e-03 ||r(i)||/||b|| 3.553147015749e-05
    4 KSP unpreconditioned resid norm 1.574461437767e-10 true resid norm 5.364626737586e+10 ||r(i)||/||b|| 6.694171667086e+08
  Linear solve converged due to CONVERGED_RTOL iterations 4
 1 Nonlinear |R| = 6.368477e+01
    0 KSP unpreconditioned resid norm 6.368477238213e+01 true resid norm 6.368477238213e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.613035793071e-05 true resid norm 1.613035793071e-05 ||r(i)||/||b|| 2.532843775263e-07
  Linear solve converged due to CONVERGED_RTOL iterations 1
 2 Nonlinear |R| = 2.592837e-05
    0 KSP unpreconditioned resid norm 2.592837440323e-05 true resid norm 2.592837440323e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.612941555715e-05 true resid norm 1.612941555715e-05 ||r(i)||/||b|| 6.220758504299e-01
    2 KSP unpreconditioned resid norm 2.164856724535e-06 true resid norm 2.164851377205e-06 ||r(i)||/||b|| 8.349352502929e-02
    3 KSP unpreconditioned resid norm 1.441082005553e-09 true resid norm 1.441500731926e-09 ||r(i)||/||b|| 5.559549200842e-05
    4 KSP unpreconditioned resid norm 1.382484428583e-14 true resid norm 2.572282442137e-11 ||r(i)||/||b|| 9.920723922502e-07
  Linear solve converged due to CONVERGED_RTOL iterations 4
 3 Nonlinear |R| = 2.197897e-10
 Solve Converged!

Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876862179e+01 true resid norm 8.013876862179e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.652519922723e+01 true resid norm 5.652519922723e+01 ||r(i)||/||b|| 7.053414994932e-01
    2 KSP unpreconditioned resid norm 5.425203302086e-02 true resid norm 5.426022278865e-02 ||r(i)||/||b|| 6.770783195426e-04
    3 KSP unpreconditioned resid norm 3.861931144240e-06 true resid norm 1.661458592719e-05 ||r(i)||/||b|| 2.073227005222e-07
  Linear solve converged due to CONVERGED_RTOL iterations 3
 1 Nonlinear |R| = 2.507082e+01
    0 KSP unpreconditioned resid norm 2.507081790712e+01 true resid norm 2.507081790712e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.083027563914e-05 true resid norm 2.083027563914e-05 ||r(i)||/||b|| 8.308574421589e-07
  Linear solve converged due to CONVERGED_RTOL iterations 1
 2 Nonlinear |R| = 2.106457e-05
    0 KSP unpreconditioned resid norm 2.106457323938e-05 true resid norm 2.106457323938e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.083003710327e-05 true resid norm 2.083003710327e-05 ||r(i)||/||b|| 9.888658491466e-01
    2 KSP unpreconditioned resid norm 2.766855053112e-06 true resid norm 2.766853985637e-06 ||r(i)||/||b|| 1.313510582054e-01
    3 KSP unpreconditioned resid norm 1.640631997317e-08 true resid norm 1.641307218646e-08 ||r(i)||/||b|| 7.791789560580e-04
    4 KSP unpreconditioned resid norm 1.963181052536e-13 true resid norm 1.781088245454e-10 ||r(i)||/||b|| 8.455373034214e-06
  Linear solve converged due to CONVERGED_RTOL iterations 4
 3 Nonlinear |R| = 2.306853e-10
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 3

Now with PJFNK, -pc_type none (consequently, ideally, the solve should be the exact same!!!):
Time Step  1, time = 0.5
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.013876853448e+01 true resid norm 8.013876853448e+01 ||r(i)||/||b|| 1.000000000000e+00
    2 KSP unpreconditioned resid norm 5.666666668339e+01 true resid norm 2.682244954493e+10 ||r(i)||/||b|| 3.347000463750e+08
    3 KSP unpreconditioned resid norm 4.598990406557e+01 true resid norm 3.614102095189e+10 ||r(i)||/||b|| 4.509804881310e+08
    4 KSP unpreconditioned resid norm 3.988651465795e+01 true resid norm 4.059883892809e+10 ||r(i)||/||b|| 5.066067231944e+08
    5 KSP unpreconditioned resid norm 3.570686665175e+01 true resid norm 4.328614865077e+10 ||r(i)||/||b|| 5.401399278072e+08
    6 KSP unpreconditioned resid norm 3.261484079094e+01 true resid norm 4.508296868728e+10 ||r(i)||/||b|| 5.625612860258e+08
    7 KSP unpreconditioned resid norm 3.020811193692e+01 true resid norm 4.636900545360e+10 ||r(i)||/||b|| 5.786089092902e+08
    8 KSP unpreconditioned resid norm 2.826597650826e+01 true resid norm 4.733495461134e+10 ||r(i)||/||b|| 5.906623657560e+08
    9 KSP unpreconditioned resid norm 2.665593183460e+01 true resid norm 4.808709219408e+10 ||r(i)||/||b|| 6.000478054937e+08
   10 KSP unpreconditioned resid norm 2.529298407419e+01 true resid norm 4.868933452889e+10 ||r(i)||/||b|| 6.075627991207e+08
   11 KSP unpreconditioned resid norm 2.411977513173e+01 true resid norm 4.918243053177e+10 ||r(i)||/||b|| 6.137158260751e+08
   12 KSP unpreconditioned resid norm 2.309600845970e+01 true resid norm 4.959358626034e+10 ||r(i)||/||b|| 6.188463732008e+08
   13 KSP unpreconditioned resid norm 2.219243356461e+01 true resid norm 4.994165949597e+10 ||r(i)||/||b|| 6.231897545878e+08
   14 KSP unpreconditioned resid norm 2.138723492170e+01 true resid norm 5.024013371969e+10 ||r(i)||/||b|| 6.269142218984e+08
   15 KSP unpreconditioned resid norm 2.066376626487e+01 true resid norm 5.049890531136e+10 ||r(i)||/||b|| 6.301432656734e+08
   16 KSP unpreconditioned resid norm 2.000907536300e+01 true resid norm 5.072540203507e+10 ||r(i)||/||b|| 6.329695721896e+08
   17 KSP unpreconditioned resid norm 1.941291274250e+01 true resid norm 5.092530758854e+10 ||r(i)||/||b|| 6.354640646447e+08
   18 KSP unpreconditioned resid norm 1.886704727044e+01 true resid norm 5.110304511196e+10 ||r(i)||/||b|| 6.376819365520e+08
   19 KSP unpreconditioned resid norm 1.836478230456e+01 true resid norm 5.126210829474e+10 ||r(i)||/||b|| 6.396667834082e+08
   20 KSP unpreconditioned resid norm 1.790060645485e+01 true resid norm 5.140529332644e+10 ||r(i)||/||b|| 6.414534970590e+08
   21 KSP unpreconditioned resid norm 1.746993682433e+01 true resid norm 5.153486468777e+10 ||r(i)||/||b|| 6.430703344986e+08
   22 KSP unpreconditioned resid norm 1.706892711482e+01 true resid norm 5.165267579888e+10 ||r(i)||/||b|| 6.445404233615e+08
   23 KSP unpreconditioned resid norm 1.669432208044e+01 true resid norm 5.176025824260e+10 ||r(i)||/||b|| 6.458828752819e+08
   24 KSP unpreconditioned resid norm 1.634334565486e+01 true resid norm 5.185888871485e+10 ||r(i)||/||b|| 6.471136213248e+08
   25 KSP unpreconditioned resid norm 1.601361391615e+01 true resid norm 5.194963992980e+10 ||r(i)||/||b|| 6.482460471981e+08
   26 KSP unpreconditioned resid norm 1.570306662467e+01 true resid norm 5.203341979356e+10 ||r(i)||/||b|| 6.492914820768e+08
   27 KSP unpreconditioned resid norm 1.540991282463e+01 true resid norm 5.211100188422e+10 ||r(i)||/||b|| 6.502595789427e+08
   28 KSP unpreconditioned resid norm 1.513258721787e+01 true resid norm 5.218304940855e+10 ||r(i)||/||b|| 6.511586135255e+08
   29 KSP unpreconditioned resid norm 1.486971487625e+01 true resid norm 5.225013420834e+10 ||r(i)||/||b|| 6.519957214699e+08
   30 KSP unpreconditioned resid norm 5.231275196988e+10 true resid norm 5.231275196988e+10 ||r(i)||/||b|| 6.527770881252e+08
  Linear solve did not converge due to DIVERGED_DTOL iterations 30
 1 Nonlinear |R| = 7.882883e+05
    0 KSP unpreconditioned resid norm 7.882882770519e+05 true resid norm 7.882882770519e+05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.481485980166e+01 true resid norm 5.481485980166e+01 ||r(i)||/||b|| 6.953656599672e-05
    2 KSP unpreconditioned resid norm 1.424810544089e+00 true resid norm 1.424823265114e+00 ||r(i)||/||b|| 1.807490110652e-06
  Linear solve converged due to CONVERGED_RTOL iterations 2
 2 Nonlinear |R| = 5.491415e+01
    0 KSP unpreconditioned resid norm 5.491414528393e+01 true resid norm 5.491414528393e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.482363716861e+01 true resid norm 5.482363716861e+01 ||r(i)||/||b|| 9.983518251109e-01
    2 KSP unpreconditioned resid norm 6.348263848927e-02 true resid norm 6.348218827720e-02 ||r(i)||/||b|| 1.156026155901e-03
    3 KSP unpreconditioned resid norm 5.738894637156e-06 true resid norm 4.343251050332e-04 ||r(i)||/||b|| 7.909166259213e-06
  Linear solve converged due to CONVERGED_RTOL iterations 3
 3 Nonlinear |R| = 1.741438e-01
    0 KSP unpreconditioned resid norm 1.741438159951e-01 true resid norm 1.741438159951e-01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.736821033825e-01 true resid norm 1.736821033825e-01 ||r(i)||/||b|| 9.973486706375e-01
    2 KSP unpreconditioned resid norm 1.682681670070e-06 true resid norm 1.723678229949e-06 ||r(i)||/||b|| 9.898015729698e-06
  Linear solve converged due to CONVERGED_RTOL iterations 2
 4 Nonlinear |R| = 1.902614e-06
    0 KSP unpreconditioned resid norm 1.902613555683e-06 true resid norm 1.902613555683e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.693336133901e-06 true resid norm 1.693336133901e-06 ||r(i)||/||b|| 8.900052923747e-01
    2 KSP unpreconditioned resid norm 4.388570424434e-09 true resid norm 4.388863055209e-09 ||r(i)||/||b|| 2.306754854185e-03
    3 KSP unpreconditioned resid norm 4.323046090618e-13 true resid norm 6.477780452255e-12 ||r(i)||/||b|| 3.404674813184e-06
  Linear solve converged due to CONVERGED_RTOL iterations 3
 5 Nonlinear |R| = 8.577249e-13
Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 5


Time Step  2, time = 1
                dt = 0.5
 0 Nonlinear |R| = 8.013877e+01
    0 KSP unpreconditioned resid norm 8.013876862179e+01 true resid norm 8.013876862179e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.652520357021e+01 true resid norm 5.652520357021e+01 ||r(i)||/||b|| 7.053415536865e-01
    2 KSP unpreconditioned resid norm 5.425812913652e-02 true resid norm 5.425162338339e-02 ||r(i)||/||b|| 6.769710131113e-04
    3 KSP unpreconditioned resid norm 3.862134058663e-06 true resid norm 1.921027115546e-05 ||r(i)||/||b|| 2.397125821352e-07
  Linear solve converged due to CONVERGED_RTOL iterations 3
 1 Nonlinear |R| = 2.507082e+01
    0 KSP unpreconditioned resid norm 2.507082225474e+01 true resid norm 2.507082225474e+01 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.093681131650e-05 true resid norm 2.093681131650e-05 ||r(i)||/||b|| 8.351066871187e-07
  Linear solve converged due to CONVERGED_RTOL iterations 1
 2 Nonlinear |R| = 2.116998e-05
    0 KSP unpreconditioned resid norm 2.116998389902e-05 true resid norm 2.116998389902e-05 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.093656388997e-05 true resid norm 2.093656388997e-05 ||r(i)||/||b|| 9.889740110262e-01
    2 KSP unpreconditioned resid norm 2.547729514728e-06 true resid norm 2.547729815384e-06 ||r(i)||/||b|| 1.203463274954e-01
    3 KSP unpreconditioned resid norm 1.711742308811e-08 true resid norm 1.711825480396e-08 ||r(i)||/||b|| 8.086097224076e-04
    4 KSP unpreconditioned resid norm 1.732897019926e-13 true resid norm 2.057200149901e-10 ||r(i)||/||b|| 9.717532898057e-06
  Linear solve converged due to CONVERGED_RTOL iterations 4
 3 Nonlinear |R| = 4.037723e-10
Nonlinear solve converged due to CONVERGED_FNORM_RELA
