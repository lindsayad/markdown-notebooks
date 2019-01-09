# 1/8/19

Initial condition number for FDP with manual unit scaling: 7.659023580798e+10
Initial condition number for FDP without any scaling: 1.349249207423e+26
Initial condition number for FDP with PETSc scaling: 3.782640170136e+13

`2d-one-spot-scaled-1st-order.i` 24 steps without surface tension:
SMP non-linear iterations: 133
SMP end time: 30.4495
final dt: 0.0593262

FDP non-linear iterations: 129
FDP end time: 31.1251
final dt: .133484

`2d-one-spot-scaled-1st-order.i` 24 steps with surface tension:
FDP non-linear: 143
FDP end time: 21.2793
final dt: .474609

I should run the `2d-one-spot-scaled-2nd-order.i` tests with surface tension and
compare SMP and FDP performance here. I know that it is worse. Combined with the
results for the 1st order case, it's definitely safe to say there is a bug in
the Jacobian somewhere. This is the number one priority for articuno currently.
