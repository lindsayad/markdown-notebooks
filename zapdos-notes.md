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
