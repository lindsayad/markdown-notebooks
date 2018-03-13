### Non-linear residuals

PJFNK, SMP, -pc_type lu, bt, no predictor, contact only evaluated on non-linear residuals (approximately):

Time Step  1, time = 0.1
                dt = 0.1
 0 Nonlinear |R| = 7.810169e-02
    0 KSP unpreconditioned resid norm 7.810169371831e-02 true resid norm 7.810169371831e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.310217323811e-02 true resid norm 1.310217323811e-02 ||r(i)||/||b|| 1.677578630416e-01
    2 KSP unpreconditioned resid norm 1.331047227087e-03 true resid norm 1.331050566366e-03 ||r(i)||/||b|| 1.704253138436e-02
    3 KSP unpreconditioned resid norm 1.972062470522e-04 true resid norm 1.972077436022e-04 ||r(i)||/||b|| 2.525012380826e-03
    4 KSP unpreconditioned resid norm 3.650798841283e-05 true resid norm 3.650827541262e-05 ||r(i)||/||b|| 4.674453737750e-04
    5 KSP unpreconditioned resid norm 5.576139303809e-06 true resid norm 5.576149903584e-06 ||r(i)||/||b|| 7.139601765483e-05
    6 KSP unpreconditioned resid norm 1.250213726376e-06 true resid norm 1.249983646812e-06 ||r(i)||/||b|| 1.600456516757e-05
    7 KSP unpreconditioned resid norm 2.313040436560e-07 true resid norm 2.307045731561e-07 ||r(i)||/||b|| 2.953899744968e-06
    8 KSP unpreconditioned resid norm 3.178315636695e-08 true resid norm 3.121668893311e-08 ||r(i)||/||b|| 3.996928548784e-07
      Line search: Using full step: fnorm 7.810169371831e-02 gnorm 1.369964436848e-02
 1 Nonlinear |R| = 1.369964e-02
    0 KSP unpreconditioned resid norm 1.369964436848e-02 true resid norm 1.369964436848e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.734334773015e-04 true resid norm 4.734334773015e-04 ||r(i)||/||b|| 3.455808520043e-02
    2 KSP unpreconditioned resid norm 1.681121848006e-05 true resid norm 1.681123072037e-05 ||r(i)||/||b|| 1.227128987308e-03
    3 KSP unpreconditioned resid norm 3.597620622885e-07 true resid norm 3.597619698927e-07 ||r(i)||/||b|| 2.626067949038e-05
    4 KSP unpreconditioned resid norm 6.848689139710e-09 true resid norm 6.848091120358e-09 ||r(i)||/||b|| 4.998736416921e-07
      Line search: Using full step: fnorm 1.369964436848e-02 gnorm 3.015701300514e-04
 2 Nonlinear |R| = 3.015701e-04
    0 KSP unpreconditioned resid norm 3.015701300514e-04 true resid norm 3.015701300514e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.658331124302e-07 true resid norm 4.658331124302e-07 ||r(i)||/||b|| 1.544692481151e-03
    2 KSP unpreconditioned resid norm 3.970593909978e-10 true resid norm 3.974983744669e-10 ||r(i)||/||b|| 1.318095974555e-06
    3 KSP unpreconditioned resid norm 1.848130802261e-13 true resid norm 2.430201742420e-12 ||r(i)||/||b|| 8.058496184638e-09
      Line search: Using full step: fnorm 3.015701300514e-04 gnorm 1.689796547995e-07
 3 Nonlinear |R| = 1.689797e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   3.000000e+00 |  -1.442323e-02 |
+----------------+------------------+----------------+----------------+


Time Step  2, time = 0.2
                dt = 0.1
 0 Nonlinear |R| = 8.727457e-02
    0 KSP unpreconditioned resid norm 8.727456965512e-02 true resid norm 8.727456965512e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.248222301654e-02 true resid norm 1.248222301654e-02 ||r(i)||/||b|| 1.430224527702e-01
    2 KSP unpreconditioned resid norm 2.039251443702e-03 true resid norm 2.039252216234e-03 ||r(i)||/||b|| 2.336593837463e-02
    3 KSP unpreconditioned resid norm 3.346550369001e-04 true resid norm 3.346537321439e-04 ||r(i)||/||b|| 3.834493065578e-03
    4 KSP unpreconditioned resid norm 7.431823068809e-05 true resid norm 7.431863240277e-05 ||r(i)||/||b|| 8.515496861967e-04
    5 KSP unpreconditioned resid norm 1.557943608725e-05 true resid norm 1.557949085153e-05 ||r(i)||/||b|| 1.785112308556e-04
    6 KSP unpreconditioned resid norm 2.500335558484e-06 true resid norm 2.500603348040e-06 ||r(i)||/||b|| 2.865214183148e-05
    7 KSP unpreconditioned resid norm 4.558708617648e-07 true resid norm 4.584492322590e-07 ||r(i)||/||b|| 5.252953226474e-06
    8 KSP unpreconditioned resid norm 9.215992031404e-08 true resid norm 9.297812505181e-08 ||r(i)||/||b|| 1.065351859301e-06
    9 KSP unpreconditioned resid norm 1.066604901502e-08 true resid norm 1.314190372894e-08 ||r(i)||/||b|| 1.505811346979e-07
      Line search: gnorm after quadratic fit 5.420139937811e-02
      Line search: Quadratically determined step, lambda=3.8748910962317201e-01
 1 Nonlinear |R| = 5.420140e-02
    0 KSP unpreconditioned resid norm 5.420139937811e-02 true resid norm 5.420139937811e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.746562685918e-03 true resid norm 4.746562685918e-03 ||r(i)||/||b|| 8.757269628419e-02
    2 KSP unpreconditioned resid norm 3.493697605954e-04 true resid norm 3.493700234406e-04 ||r(i)||/||b|| 6.445774969818e-03
    3 KSP unpreconditioned resid norm 3.867127216285e-05 true resid norm 3.867144833544e-05 ||r(i)||/||b|| 7.134769356353e-04
    4 KSP unpreconditioned resid norm 6.329866285776e-06 true resid norm 6.331562092039e-06 ||r(i)||/||b|| 1.168154727495e-04
    5 KSP unpreconditioned resid norm 6.318732588678e-07 true resid norm 6.326620553639e-07 ||r(i)||/||b|| 1.167243028082e-05
    6 KSP unpreconditioned resid norm 6.699231505182e-08 true resid norm 6.683552317048e-08 ||r(i)||/||b|| 1.233095896736e-06
    7 KSP unpreconditioned resid norm 4.300812205462e-09 true resid norm 4.974166983377e-09 ||r(i)||/||b|| 9.177192914664e-08
      Line search: Using full step: fnorm 5.420139937811e-02 gnorm 1.756961716177e-02
 2 Nonlinear |R| = 1.756962e-02
    0 KSP unpreconditioned resid norm 1.756961716177e-02 true resid norm 1.756961716177e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.729686326295e-02 true resid norm 1.729686326559e-02 ||r(i)||/||b|| 9.844758201801e-01
    2 KSP unpreconditioned resid norm 6.255437011293e-03 true resid norm 6.255436185470e-03 ||r(i)||/||b|| 3.560371366019e-01
    3 KSP unpreconditioned resid norm 4.102928530071e-03 true resid norm 4.102926979418e-03 ||r(i)||/||b|| 2.335239830009e-01
    4 KSP unpreconditioned resid norm 2.928213797981e-03 true resid norm 2.928211391472e-03 ||r(i)||/||b|| 1.666633578017e-01
    5 KSP unpreconditioned resid norm 1.591165721910e-03 true resid norm 1.591162380591e-03 ||r(i)||/||b|| 9.056329263979e-02
    6 KSP unpreconditioned resid norm 1.146135136986e-03 true resid norm 1.146132427857e-03 ||r(i)||/||b|| 6.523377358222e-02
    7 KSP unpreconditioned resid norm 1.891242678584e-05 true resid norm 1.891325891127e-05 ||r(i)||/||b|| 1.076475300351e-03
    8 KSP unpreconditioned resid norm 2.702957103082e-07 true resid norm 2.713516945487e-07 ||r(i)||/||b|| 1.544437149940e-05
    9 KSP unpreconditioned resid norm 2.334041398436e-09 true resid norm 9.493449413566e-01 ||r(i)||/||b|| 5.403333109740e+01
      Line search: gnorm after quadratic fit 8.710519797935e-03
      Line search: Quadratically determined step, lambda=4.9124132748063165e-01
 3 Nonlinear |R| = 8.710520e-03
    0 KSP unpreconditioned resid norm 8.710519797935e-03 true resid norm 8.710519797935e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.267396670620e-04 true resid norm 2.267396670620e-04 ||r(i)||/||b|| 2.603055527360e-02
    2 KSP unpreconditioned resid norm 1.215068430901e-05 true resid norm 1.215064889282e-05 ||r(i)||/||b|| 1.394939587383e-03
    3 KSP unpreconditioned resid norm 1.672185546120e-07 true resid norm 1.672704059166e-07 ||r(i)||/||b|| 1.920326338691e-05
    4 KSP unpreconditioned resid norm 2.518985030865e-09 true resid norm 2.542042223465e-09 ||r(i)||/||b|| 2.918358814898e-07
      Line search: Using full step: fnorm 8.710519797935e-03 gnorm 1.282591686427e-04
 4 Nonlinear |R| = 1.282592e-04
    0 KSP unpreconditioned resid norm 1.282591686427e-04 true resid norm 1.282591686427e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.581121939330e-07 true resid norm 7.581121939330e-07 ||r(i)||/||b|| 5.910783626276e-03
    2 KSP unpreconditioned resid norm 5.373716492204e-10 true resid norm 5.380994444727e-10 ||r(i)||/||b|| 4.195407238072e-06
    3 KSP unpreconditioned resid norm 4.376589167752e-13 true resid norm 3.526951292623e-12 ||r(i)||/||b|| 2.749862898651e-08
      Line search: Using full step: fnorm 1.282591686427e-04 gnorm 1.460116715192e-07
 5 Nonlinear |R| = 1.460117e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   3.000000e+00 |  -1.442323e-02 |
|   2.000000e-01 |    -0.000000e+00 |   5.000000e+00 |  -1.305816e-03 |
+----------------+------------------+----------------+----------------+


Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 8.270135e-02
    0 KSP unpreconditioned resid norm 8.270135001060e-02 true resid norm 8.270135001060e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.249176017799e-02 true resid norm 1.249176017799e-02 ||r(i)||/||b|| 1.510466295458e-01
    2 KSP unpreconditioned resid norm 1.760264073050e-03 true resid norm 1.760263098270e-03 ||r(i)||/||b|| 2.128457513746e-02
    3 KSP unpreconditioned resid norm 2.646258609147e-04 true resid norm 2.646242495083e-04 ||r(i)||/||b|| 3.199757313205e-03
    4 KSP unpreconditioned resid norm 6.526254014693e-05 true resid norm 6.526281371557e-05 ||r(i)||/||b|| 7.891384325311e-04
    5 KSP unpreconditioned resid norm 1.414131671245e-05 true resid norm 1.414395465571e-05 ||r(i)||/||b|| 1.710244712317e-04
    6 KSP unpreconditioned resid norm 2.396868399272e-06 true resid norm 2.396018455746e-06 ||r(i)||/||b|| 2.897193885516e-05
    7 KSP unpreconditioned resid norm 4.930504932335e-07 true resid norm 4.946972511073e-07 ||r(i)||/||b|| 5.981731266102e-06
    8 KSP unpreconditioned resid norm 9.718107554058e-08 true resid norm 9.629174525400e-08 ||r(i)||/||b|| 1.164330996310e-06
    9 KSP unpreconditioned resid norm 9.323535232231e-09 true resid norm 9.617038988544e-09 ||r(i)||/||b|| 1.162863603473e-07
      Line search: gnorm after quadratic fit 7.443321366477e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 1 Nonlinear |R| = 7.443321e-02
    0 KSP unpreconditioned resid norm 7.443321366477e-02 true resid norm 7.443321366477e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.044285765455e-02 true resid norm 1.044285763536e-02 ||r(i)||/||b|| 1.402983576981e-01
    2 KSP unpreconditioned resid norm 4.492864541457e-03 true resid norm 4.492863934508e-03 ||r(i)||/||b|| 6.036100973341e-02
    3 KSP unpreconditioned resid norm 4.163436835696e-03 true resid norm 4.163435938190e-03 ||r(i)||/||b|| 5.593518985948e-02
    4 KSP unpreconditioned resid norm 3.366835453085e-03 true resid norm 3.366836506516e-03 ||r(i)||/||b|| 4.523298593125e-02
    5 KSP unpreconditioned resid norm 1.106211903524e-03 true resid norm 1.106211757341e-03 ||r(i)||/||b|| 1.486180298922e-02
    6 KSP unpreconditioned resid norm 2.771727938598e-04 true resid norm 2.771718552922e-04 ||r(i)||/||b|| 3.723765798163e-03
    7 KSP unpreconditioned resid norm 2.132930377318e-04 true resid norm 2.132929800950e-04 ||r(i)||/||b|| 2.865561885527e-03
    8 KSP unpreconditioned resid norm 1.854286408804e-04 true resid norm 1.854287050029e-04 ||r(i)||/||b|| 2.491209177639e-03
    9 KSP unpreconditioned resid norm 5.392284450794e-05 true resid norm 5.392363603102e-05 ||r(i)||/||b|| 7.244566420830e-04
   10 KSP unpreconditioned resid norm 6.550155509635e-06 true resid norm 6.549278942901e-06 ||r(i)||/||b|| 8.798866286222e-05
   11 KSP unpreconditioned resid norm 5.943312951676e-07 true resid norm 5.937891190939e-07 ||r(i)||/||b|| 7.977475240666e-06
   12 KSP unpreconditioned resid norm 4.541500920333e-09 true resid norm 1.826932501109e-01 ||r(i)||/||b|| 2.454458716961e+00
      Line search: gnorm after quadratic fit 7.185519572874e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 2 Nonlinear |R| = 7.185520e-02
    0 KSP unpreconditioned resid norm 7.185519572874e-02 true resid norm 7.185519572874e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.383649271460e-02 true resid norm 2.383649227651e-02 ||r(i)||/||b|| 3.317295574073e-01
    2 KSP unpreconditioned resid norm 1.935666617096e-02 true resid norm 1.935666651760e-02 ||r(i)||/||b|| 2.693843683994e-01
    3 KSP unpreconditioned resid norm 1.742311930351e-02 true resid norm 1.742311941542e-02 ||r(i)||/||b|| 2.424754290726e-01
    4 KSP unpreconditioned resid norm 5.857415184732e-03 true resid norm 5.857414166156e-03 ||r(i)||/||b|| 8.151691894721e-02
    5 KSP unpreconditioned resid norm 4.434919839774e-03 true resid norm 4.434920122658e-03 ||r(i)||/||b|| 6.172024274208e-02
    6 KSP unpreconditioned resid norm 4.066700146117e-03 true resid norm 4.066703099266e-03 ||r(i)||/||b|| 5.659581131222e-02
    7 KSP unpreconditioned resid norm 3.211288184390e-03 true resid norm 3.211293669558e-03 ||r(i)||/||b|| 4.469118255110e-02
    8 KSP unpreconditioned resid norm 1.612600589292e-03 true resid norm 1.612604625708e-03 ||r(i)||/||b|| 2.244242200377e-02
    9 KSP unpreconditioned resid norm 1.438710945702e-03 true resid norm 1.438713600842e-03 ||r(i)||/||b|| 2.002240180757e-02
   10 KSP unpreconditioned resid norm 9.589729074772e-04 true resid norm 9.589743209349e-04 ||r(i)||/||b|| 1.334592872804e-02
   11 KSP unpreconditioned resid norm 1.964866903730e-04 true resid norm 1.964875880627e-04 ||r(i)||/||b|| 2.734493811755e-03
   12 KSP unpreconditioned resid norm 6.276513848915e-05 true resid norm 6.276874466884e-05 ||r(i)||/||b|| 8.735449682136e-04
   13 KSP unpreconditioned resid norm 5.392815333153e-05 true resid norm 5.392887143561e-05 ||r(i)||/||b|| 7.505215299837e-04
   14 KSP unpreconditioned resid norm 1.967866798059e-05 true resid norm 1.968140467744e-05 ||r(i)||/||b|| 2.739037097852e-04
   15 KSP unpreconditioned resid norm 7.207015292922e-07 true resid norm 7.202779446983e-07 ||r(i)||/||b|| 1.002402035640e-05
   16 KSP unpreconditioned resid norm 2.518354054349e-09 true resid norm 1.232467296549e+00 ||r(i)||/||b|| 1.715209713159e+01
      Line search: gnorm after quadratic fit 4.426463143824e-02
      Line search: Quadratically determined step, lambda=4.9718010874695862e-01
 3 Nonlinear |R| = 4.426463e-02
    0 KSP unpreconditioned resid norm 4.426463143824e-02 true resid norm 4.426463143824e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.548575679630e-02 true resid norm 2.548575612787e-02 ||r(i)||/||b|| 5.757589140537e-01
    2 KSP unpreconditioned resid norm 2.432869760798e-02 true resid norm 2.432869750315e-02 ||r(i)||/||b|| 5.496193397000e-01
    3 KSP unpreconditioned resid norm 1.019463036997e-02 true resid norm 1.019463685625e-02 ||r(i)||/||b|| 2.303111203009e-01
    4 KSP unpreconditioned resid norm 6.449405464008e-04 true resid norm 6.449597375703e-04 ||r(i)||/||b|| 1.457054349295e-02
    5 KSP unpreconditioned resid norm 4.162183278249e-05 true resid norm 4.162201815884e-05 ||r(i)||/||b|| 9.402996660417e-04
    6 KSP unpreconditioned resid norm 2.690672382372e-06 true resid norm 2.676856839125e-06 ||r(i)||/||b|| 6.047394391751e-05
    7 KSP unpreconditioned resid norm 1.766054635375e-07 true resid norm 1.878138009101e-07 ||r(i)||/||b|| 4.242976724478e-06
    8 KSP unpreconditioned resid norm 1.230508338462e-08 true resid norm 1.416871201501e+00 ||r(i)||/||b|| 3.200910423207e+01
      Line search: gnorm after quadratic fit 3.117032547383e-02
      Line search: Quadratically determined step, lambda=1.6075990840937815e-01
 4 Nonlinear |R| = 3.117033e-02
    0 KSP unpreconditioned resid norm 3.117032547383e-02 true resid norm 3.117032547383e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.704161770433e-03 true resid norm 1.704161770433e-03 ||r(i)||/||b|| 5.467256900680e-02
    2 KSP unpreconditioned resid norm 1.137309157638e-04 true resid norm 1.137305930093e-04 ||r(i)||/||b|| 3.648681599580e-03
    3 KSP unpreconditioned resid norm 7.334929069249e-06 true resid norm 7.334655101324e-06 ||r(i)||/||b|| 2.353089032542e-04
    4 KSP unpreconditioned resid norm 5.346442764444e-07 true resid norm 5.354392723334e-07 ||r(i)||/||b|| 1.717785310849e-05
    5 KSP unpreconditioned resid norm 3.181466224270e-08 true resid norm 3.264760072402e-08 ||r(i)||/||b|| 1.047393642117e-06
    6 KSP unpreconditioned resid norm 1.665637196894e-09 true resid norm 1.755593312815e-09 ||r(i)||/||b|| 5.632258521934e-08
      Line search: gnorm after quadratic fit 2.805931718847e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 5 Nonlinear |R| = 2.805932e-02
    0 KSP unpreconditioned resid norm 2.805931718847e-02 true resid norm 2.805931718847e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.374594148593e-03 true resid norm 1.374594148593e-03 ||r(i)||/||b|| 4.898886666986e-02
    2 KSP unpreconditioned resid norm 8.161474125246e-05 true resid norm 8.161465985113e-05 ||r(i)||/||b|| 2.908647395193e-03
    3 KSP unpreconditioned resid norm 4.855929140190e-06 true resid norm 4.855906332975e-06 ||r(i)||/||b|| 1.730586065355e-04
    4 KSP unpreconditioned resid norm 3.182197900292e-07 true resid norm 3.180562982979e-07 ||r(i)||/||b|| 1.133514034435e-05
    5 KSP unpreconditioned resid norm 1.677837961654e-08 true resid norm 1.634994463345e-08 ||r(i)||/||b|| 5.826921775621e-07
      Line search: gnorm after quadratic fit 2.525815336048e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 6 Nonlinear |R| = 2.525815e-02
    0 KSP unpreconditioned resid norm 2.525815336048e-02 true resid norm 2.525815336048e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.109505059285e-03 true resid norm 1.109505059285e-03 ||r(i)||/||b|| 4.392661028899e-02
    2 KSP unpreconditioned resid norm 5.870873954595e-05 true resid norm 5.870899292680e-05 ||r(i)||/||b|| 2.324358083068e-03
    3 KSP unpreconditioned resid norm 3.211047339602e-06 true resid norm 3.210750345013e-06 ||r(i)||/||b|| 1.271173826206e-04
    4 KSP unpreconditioned resid norm 1.889711383770e-07 true resid norm 1.889800639243e-07 ||r(i)||/||b|| 7.481943007759e-06
    5 KSP unpreconditioned resid norm 8.835318326413e-09 true resid norm 9.347945223484e-09 ||r(i)||/||b|| 3.700961463838e-07
      Line search: gnorm after quadratic fit 2.273611887956e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 7 Nonlinear |R| = 2.273612e-02
    0 KSP unpreconditioned resid norm 2.273611887956e-02 true resid norm 2.273611887956e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.960267451533e-04 true resid norm 8.960267451533e-04 ||r(i)||/||b|| 3.940983726818e-02
    2 KSP unpreconditioned resid norm 4.231458704770e-05 true resid norm 4.231447716501e-05 ||r(i)||/||b|| 1.861112593102e-03
    3 KSP unpreconditioned resid norm 2.120926184284e-06 true resid norm 2.120774713961e-06 ||r(i)||/||b|| 9.327778083829e-05
    4 KSP unpreconditioned resid norm 1.120275515118e-07 true resid norm 1.122862508409e-07 ||r(i)||/||b|| 4.938672753943e-06
    5 KSP unpreconditioned resid norm 4.649005256096e-09 true resid norm 4.913156805494e-09 ||r(i)||/||b|| 2.160947887158e-07
      Line search: gnorm after quadratic fit 2.046551066701e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 8 Nonlinear |R| = 2.046551e-02
    0 KSP unpreconditioned resid norm 2.046551066701e-02 true resid norm 2.046551066701e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 7.239490474353e-04 true resid norm 7.239490474353e-04 ||r(i)||/||b|| 3.537410129727e-02
    2 KSP unpreconditioned resid norm 3.054779058742e-05 true resid norm 3.054774805932e-05 ||r(i)||/||b|| 1.492645287789e-03
    3 KSP unpreconditioned resid norm 1.399416616160e-06 true resid norm 1.399394044520e-06 ||r(i)||/||b|| 6.837816398960e-05
    4 KSP unpreconditioned resid norm 6.632958972164e-08 true resid norm 6.609540434420e-08 ||r(i)||/||b|| 3.229599564830e-06
    5 KSP unpreconditioned resid norm 2.445556988820e-09 true resid norm 2.917237344635e-09 ||r(i)||/||b|| 1.425440777951e-07
      Line search: gnorm after quadratic fit 1.842134992577e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 9 Nonlinear |R| = 1.842135e-02
    0 KSP unpreconditioned resid norm 1.842134992577e-02 true resid norm 1.842134992577e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 5.851358754362e-04 true resid norm 5.851358754362e-04 ||r(i)||/||b|| 3.176400631842e-02
    2 KSP unpreconditioned resid norm 2.208276512913e-05 true resid norm 2.208270735789e-05 ||r(i)||/||b|| 1.198756195765e-03
    3 KSP unpreconditioned resid norm 9.224898779418e-07 true resid norm 9.224781208090e-07 ||r(i)||/||b|| 5.007657552383e-05
    4 KSP unpreconditioned resid norm 3.923606896161e-08 true resid norm 3.925158888171e-08 ||r(i)||/||b|| 2.130766151225e-06
    5 KSP unpreconditioned resid norm 1.286519308040e-09 true resid norm 1.428445535414e-09 ||r(i)||/||b|| 7.754293475611e-08
      Line search: gnorm after quadratic fit 1.658112019684e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
10 Nonlinear |R| = 1.658112e-02
    0 KSP unpreconditioned resid norm 1.658112019684e-02 true resid norm 1.658112019684e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.730863270135e-04 true resid norm 4.730863270135e-04 ||r(i)||/||b|| 2.853162641590e-02
    2 KSP unpreconditioned resid norm 1.598136702417e-05 true resid norm 1.598123067054e-05 ||r(i)||/||b|| 9.638209289133e-04
    3 KSP unpreconditioned resid norm 6.076068887155e-07 true resid norm 6.075486212191e-07 ||r(i)||/||b|| 3.664098770207e-05
    4 KSP unpreconditioned resid norm 2.319339152917e-08 true resid norm 2.307125787435e-08 ||r(i)||/||b|| 1.391417322862e-06
    5 KSP unpreconditioned resid norm 6.769667634337e-10 true resid norm 6.558154724067e-10 ||r(i)||/||b|| 3.955194007529e-08
      Line search: gnorm after quadratic fit 1.478867991420e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
11 Nonlinear |R| = 1.478868e-02
    0 KSP unpreconditioned resid norm 1.478867991420e-02 true resid norm 1.478867991420e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.109467030414e-03 true resid norm 1.109467043021e-03 ||r(i)||/||b|| 7.502137103909e-02
    2 KSP unpreconditioned resid norm 1.026924977810e-03 true resid norm 1.026924980639e-03 ||r(i)||/||b|| 6.943993558568e-02
    3 KSP unpreconditioned resid norm 3.801646422070e-04 true resid norm 3.801644920798e-04 ||r(i)||/||b|| 2.570645211645e-02
    4 KSP unpreconditioned resid norm 1.507562976634e-05 true resid norm 1.507585530510e-05 ||r(i)||/||b|| 1.019418595342e-03
    5 KSP unpreconditioned resid norm 4.576702204255e-07 true resid norm 4.577080845600e-07 ||r(i)||/||b|| 3.094989459611e-05
    6 KSP unpreconditioned resid norm 1.181970941189e-08 true resid norm 3.768461475530e-02 ||r(i)||/||b|| 2.548206802361e+00
      Line search: gnorm after quadratic fit 1.376819326973e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
12 Nonlinear |R| = 1.376819e-02
    0 KSP unpreconditioned resid norm 1.376819326973e-02 true resid norm 1.376819326973e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.838438478941e-03 true resid norm 2.838438381946e-03 ||r(i)||/||b|| 2.061591035468e-01
    2 KSP unpreconditioned resid norm 2.588193295185e-03 true resid norm 2.588193253250e-03 ||r(i)||/||b|| 1.879835068078e-01
    3 KSP unpreconditioned resid norm 3.578050501119e-04 true resid norm 3.578051863561e-04 ||r(i)||/||b|| 2.598780968182e-02
    4 KSP unpreconditioned resid norm 9.282728701153e-06 true resid norm 9.282631506202e-06 ||r(i)||/||b|| 6.742083964359e-04
    5 KSP unpreconditioned resid norm 2.748177047418e-07 true resid norm 2.744583957334e-07 ||r(i)||/||b|| 1.993423467819e-05
    6 KSP unpreconditioned resid norm 6.123656561037e-09 true resid norm 9.943267241448e-02 ||r(i)||/||b|| 7.221911435040e+00
      Line search: gnorm after quadratic fit 1.220687100766e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
13 Nonlinear |R| = 1.220687e-02
    0 KSP unpreconditioned resid norm 1.220687100766e-02 true resid norm 1.220687100766e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.418944862078e-04 true resid norm 2.418944862078e-04 ||r(i)||/||b|| 1.981625643919e-02
    2 KSP unpreconditioned resid norm 5.965846041420e-06 true resid norm 5.965900426812e-06 ||r(i)||/||b|| 4.887329785881e-04
    3 KSP unpreconditioned resid norm 1.783339426252e-07 true resid norm 1.784393780082e-07 ||r(i)||/||b|| 1.461794573697e-05
    4 KSP unpreconditioned resid norm 4.859181798409e-09 true resid norm 4.949197793944e-09 ||r(i)||/||b|| 4.054436055595e-07
      Line search: gnorm after quadratic fit 1.233977727332e-02
      Line search: Cubically determined step, current gnorm 1.124181750522e-02 lambda=5.0000000000000003e-02
14 Nonlinear |R| = 1.124182e-02
    0 KSP unpreconditioned resid norm 1.124181750522e-02 true resid norm 1.124181750522e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.062895212428e-03 true resid norm 2.062895208685e-03 ||r(i)||/||b|| 1.835019299794e-01
    2 KSP unpreconditioned resid norm 2.028400009105e-03 true resid norm 2.028399991315e-03 ||r(i)||/||b|| 1.804334566339e-01
    3 KSP unpreconditioned resid norm 4.037572984155e-04 true resid norm 4.037574673300e-04 ||r(i)||/||b|| 3.591567530273e-02
    4 KSP unpreconditioned resid norm 3.739140637180e-04 true resid norm 3.739140129542e-04 ||r(i)||/||b|| 3.326099296493e-02
    5 KSP unpreconditioned resid norm 1.344376208201e-04 true resid norm 1.344378271261e-04 ||r(i)||/||b|| 1.195872705314e-02
    6 KSP unpreconditioned resid norm 2.911919618513e-06 true resid norm 2.912013233610e-06 ||r(i)||/||b|| 2.590340247258e-04
    7 KSP unpreconditioned resid norm 6.319474202687e-08 true resid norm 6.326022355934e-08 ||r(i)||/||b|| 5.627223847919e-06
    8 KSP unpreconditioned resid norm 1.136756991998e-09 true resid norm 8.289793386169e-02 ||r(i)||/||b|| 7.374068634648e+00
      Line search: gnorm after quadratic fit 1.236656199818e-02
      Line search: Cubically determined step, current gnorm 1.084982075398e-02 lambda=3.9704611320898059e-02
15 Nonlinear |R| = 1.084982e-02
    0 KSP unpreconditioned resid norm 1.084982075398e-02 true resid norm 1.084982075398e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.370768857163e-03 true resid norm 1.370768857672e-03 ||r(i)||/||b|| 1.263402307516e-01
    2 KSP unpreconditioned resid norm 1.232582896930e-03 true resid norm 1.232582875448e-03 ||r(i)||/||b|| 1.136039851162e-01
    3 KSP unpreconditioned resid norm 2.121731188286e-04 true resid norm 2.121731171135e-04 ||r(i)||/||b|| 1.955544906451e-02
    4 KSP unpreconditioned resid norm 9.056497644294e-05 true resid norm 9.056504197218e-05 ||r(i)||/||b|| 8.347146374651e-03
    5 KSP unpreconditioned resid norm 6.549688826863e-05 true resid norm 6.549689154681e-05 ||r(i)||/||b|| 6.036679594250e-03
    6 KSP unpreconditioned resid norm 2.294674338373e-06 true resid norm 2.294903284546e-06 ||r(i)||/||b|| 2.115153177719e-04
    7 KSP unpreconditioned resid norm 4.673150727796e-08 true resid norm 4.651546052820e-08 ||r(i)||/||b|| 4.287210045488e-06
    8 KSP unpreconditioned resid norm 1.007814808541e-09 true resid norm 5.292598919613e-02 ||r(i)||/||b|| 4.878051941708e+00
      Line search: gnorm after quadratic fit 9.856304828714e-03
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
16 Nonlinear |R| = 9.856305e-03
    0 KSP unpreconditioned resid norm 9.856304828714e-03 true resid norm 9.856304828714e-03 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.906649099567e-03 true resid norm 1.906649106067e-03 ||r(i)||/||b|| 1.934446163345e-01
    2 KSP unpreconditioned resid norm 1.569127649138e-03 true resid norm 1.569127585886e-03 ||r(i)||/||b|| 1.592003913388e-01
    3 KSP unpreconditioned resid norm 2.873699305956e-04 true resid norm 2.873699067363e-04 ||r(i)||/||b|| 2.915594756152e-02
    4 KSP unpreconditioned resid norm 2.652413345526e-04 true resid norm 2.652413057476e-04 ||r(i)||/||b|| 2.691082615210e-02
    5 KSP unpreconditioned resid norm 4.931938161391e-05 true resid norm 4.931937387542e-05 ||r(i)||/||b|| 5.003840154349e-03
    6 KSP unpreconditioned resid norm 1.124945309014e-06 true resid norm 1.124950613669e-06 ||r(i)||/||b|| 1.141351280443e-04
    7 KSP unpreconditioned resid norm 2.231632674374e-08 true resid norm 2.202216183936e-08 ||r(i)||/||b|| 2.234322316738e-06
    8 KSP unpreconditioned resid norm 4.037062051928e-10 true resid norm 4.050031540588e-02 ||r(i)||/||b|| 4.109076992820e+00
      Line search: gnorm after quadratic fit 1.046012034126e-02
      Line search: Cubically determined step, current gnorm 9.039066354949e-03 lambda=4.3822192096195166e-02
17 Nonlinear |R| = 9.039066e-03

### Every residual

PJFNK, SMP, -pc_type lu, bt, no predictor, contact evaluated on every residual:

Time Step  1, time = 0.1
                dt = 0.1
 0 Nonlinear |R| = 7.810169e-02
    0 KSP unpreconditioned resid norm 7.810169371831e-02 true resid norm 7.810169371831e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.310217323811e-02 true resid norm 1.310217323811e-02 ||r(i)||/||b|| 1.677578630416e-01
    2 KSP unpreconditioned resid norm 1.331047227087e-03 true resid norm 1.331050566366e-03 ||r(i)||/||b|| 1.704253138436e-02
    3 KSP unpreconditioned resid norm 1.972062470522e-04 true resid norm 1.972077436022e-04 ||r(i)||/||b|| 2.525012380826e-03
    4 KSP unpreconditioned resid norm 3.650798841283e-05 true resid norm 3.650827541262e-05 ||r(i)||/||b|| 4.674453737750e-04
    5 KSP unpreconditioned resid norm 5.576139303809e-06 true resid norm 5.576149903584e-06 ||r(i)||/||b|| 7.139601765483e-05
    6 KSP unpreconditioned resid norm 1.250213726376e-06 true resid norm 1.249983646812e-06 ||r(i)||/||b|| 1.600456516757e-05
    7 KSP unpreconditioned resid norm 2.313040436560e-07 true resid norm 2.307045731561e-07 ||r(i)||/||b|| 2.953899744968e-06
    8 KSP unpreconditioned resid norm 3.178315636695e-08 true resid norm 3.121668893311e-08 ||r(i)||/||b|| 3.996928548784e-07
      Line search: Using full step: fnorm 7.810169371831e-02 gnorm 1.369964436848e-02
 1 Nonlinear |R| = 1.369964e-02
    0 KSP unpreconditioned resid norm 1.369964436848e-02 true resid norm 1.369964436848e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.734334773015e-04 true resid norm 4.734334773015e-04 ||r(i)||/||b|| 3.455808520043e-02
    2 KSP unpreconditioned resid norm 1.681121848006e-05 true resid norm 1.681123072037e-05 ||r(i)||/||b|| 1.227128987308e-03
    3 KSP unpreconditioned resid norm 3.597620622885e-07 true resid norm 3.597619698927e-07 ||r(i)||/||b|| 2.626067949038e-05
    4 KSP unpreconditioned resid norm 6.848689139710e-09 true resid norm 6.848091120358e-09 ||r(i)||/||b|| 4.998736416921e-07
      Line search: Using full step: fnorm 1.369964436848e-02 gnorm 3.015701300514e-04
 2 Nonlinear |R| = 3.015701e-04
    0 KSP unpreconditioned resid norm 3.015701300514e-04 true resid norm 3.015701300514e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.658331124302e-07 true resid norm 4.658331124302e-07 ||r(i)||/||b|| 1.544692481151e-03
    2 KSP unpreconditioned resid norm 3.970593909978e-10 true resid norm 3.974983744669e-10 ||r(i)||/||b|| 1.318095974555e-06
    3 KSP unpreconditioned resid norm 1.848130802261e-13 true resid norm 2.430201742420e-12 ||r(i)||/||b|| 8.058496184638e-09
      Line search: Using full step: fnorm 3.015701300514e-04 gnorm 1.689796547995e-07
 3 Nonlinear |R| = 1.689797e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   3.000000e+00 |  -1.442323e-02 |
+----------------+------------------+----------------+----------------+


Time Step  2, time = 0.2
                dt = 0.1
 0 Nonlinear |R| = 8.727457e-02
    0 KSP unpreconditioned resid norm 8.727456965512e-02 true resid norm 8.727456965512e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.248222301654e-02 true resid norm 1.248222301654e-02 ||r(i)||/||b|| 1.430224527702e-01
    2 KSP unpreconditioned resid norm 2.039251443702e-03 true resid norm 2.039252216234e-03 ||r(i)||/||b|| 2.336593837463e-02
    3 KSP unpreconditioned resid norm 3.346550369001e-04 true resid norm 3.346537321439e-04 ||r(i)||/||b|| 3.834493065578e-03
    4 KSP unpreconditioned resid norm 7.431823068809e-05 true resid norm 7.431863240277e-05 ||r(i)||/||b|| 8.515496861967e-04
    5 KSP unpreconditioned resid norm 1.557943608725e-05 true resid norm 1.557949085153e-05 ||r(i)||/||b|| 1.785112308556e-04
    6 KSP unpreconditioned resid norm 2.500335558484e-06 true resid norm 2.500603348040e-06 ||r(i)||/||b|| 2.865214183148e-05
    7 KSP unpreconditioned resid norm 4.558708617648e-07 true resid norm 4.584492322590e-07 ||r(i)||/||b|| 5.252953226474e-06
    8 KSP unpreconditioned resid norm 9.215992031404e-08 true resid norm 9.297812505181e-08 ||r(i)||/||b|| 1.065351859301e-06
    9 KSP unpreconditioned resid norm 1.066604901502e-08 true resid norm 1.314190372894e-08 ||r(i)||/||b|| 1.505811346979e-07
      Line search: gnorm after quadratic fit 5.420139937811e-02
      Line search: Quadratically determined step, lambda=3.8748910962317201e-01
 1 Nonlinear |R| = 5.420140e-02
    0 KSP unpreconditioned resid norm 5.420139937811e-02 true resid norm 5.420139937811e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.746562685918e-03 true resid norm 4.746562685918e-03 ||r(i)||/||b|| 8.757269628419e-02
    2 KSP unpreconditioned resid norm 3.493697605954e-04 true resid norm 3.493700234406e-04 ||r(i)||/||b|| 6.445774969818e-03
    3 KSP unpreconditioned resid norm 3.867127216285e-05 true resid norm 3.867144833544e-05 ||r(i)||/||b|| 7.134769356353e-04
    4 KSP unpreconditioned resid norm 6.329866285776e-06 true resid norm 6.331562092039e-06 ||r(i)||/||b|| 1.168154727495e-04
    5 KSP unpreconditioned resid norm 6.318732588678e-07 true resid norm 6.326620553639e-07 ||r(i)||/||b|| 1.167243028082e-05
    6 KSP unpreconditioned resid norm 6.699231505182e-08 true resid norm 6.683552317048e-08 ||r(i)||/||b|| 1.233095896736e-06
    7 KSP unpreconditioned resid norm 4.300812205462e-09 true resid norm 4.974166983377e-09 ||r(i)||/||b|| 9.177192914664e-08
      Line search: Using full step: fnorm 5.420139937811e-02 gnorm 1.756961716177e-02
 2 Nonlinear |R| = 1.756962e-02
    0 KSP unpreconditioned resid norm 1.756961716177e-02 true resid norm 1.756961716177e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.301487359162e-04 true resid norm 1.301487462819e-04 ||r(i)||/||b|| 7.407602856882e-03
    2 KSP unpreconditioned resid norm 3.087391513471e-06 true resid norm 3.087386048957e-06 ||r(i)||/||b|| 1.757230120913e-04
    3 KSP unpreconditioned resid norm 3.371705000643e-07 true resid norm 3.370977119948e-07 ||r(i)||/||b|| 1.918640052831e-05
    4 KSP unpreconditioned resid norm 1.247542885111e-08 true resid norm 1.251158428108e-08 ||r(i)||/||b|| 7.121147925920e-07
      Line search: Using full step: fnorm 1.756961716177e-02 gnorm 6.592553953246e-04
 3 Nonlinear |R| = 6.592554e-04
    0 KSP unpreconditioned resid norm 6.592553953246e-04 true resid norm 6.592553953246e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.526893236811e-06 true resid norm 1.526893236811e-06 ||r(i)||/||b|| 2.316087585539e-03
    2 KSP unpreconditioned resid norm 3.300525703301e-09 true resid norm 3.296842692530e-09 ||r(i)||/||b|| 5.000858113427e-06
    3 KSP unpreconditioned resid norm 6.696860251112e-12 true resid norm 1.009016859377e-11 ||r(i)||/||b|| 1.530540161723e-08
      Line search: Using full step: fnorm 6.592553953246e-04 gnorm 5.352464443926e-07
 4 Nonlinear |R| = 5.352464e-07
 Solve Converged!

Postprocessor Values:
+----------------+------------------+----------------+----------------+
| time           | contact_pressure | nonlinear_its  | penetration    |
+----------------+------------------+----------------+----------------+
|   0.000000e+00 |     0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
|   1.000000e-01 |    -0.000000e+00 |   3.000000e+00 |  -1.442323e-02 |
|   2.000000e-01 |    -0.000000e+00 |   4.000000e+00 |  -1.305742e-03 |
+----------------+------------------+----------------+----------------+


Time Step  3, time = 0.3
                dt = 0.1
 0 Nonlinear |R| = 8.270140e-02
    0 KSP unpreconditioned resid norm 8.270139583989e-02 true resid norm 8.270139583989e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.249162976070e-02 true resid norm 1.249162976070e-02 ||r(i)||/||b|| 1.510449688767e-01
    2 KSP unpreconditioned resid norm 1.994928746415e-03 true resid norm 1.994929906782e-03 ||r(i)||/||b|| 2.412208266283e-02
    3 KSP unpreconditioned resid norm 2.865500900947e-04 true resid norm 2.865523793586e-04 ||r(i)||/||b|| 3.464903783648e-03
    4 KSP unpreconditioned resid norm 6.561239134351e-05 true resid norm 6.561192725166e-05 ||r(i)||/||b|| 7.933593693955e-04
    5 KSP unpreconditioned resid norm 1.390852670147e-05 true resid norm 1.390873205341e-05 ||r(i)||/||b|| 1.681801366490e-04
    6 KSP unpreconditioned resid norm 2.520385248141e-06 true resid norm 2.520135094602e-06 ||r(i)||/||b|| 3.047270325982e-05
    7 KSP unpreconditioned resid norm 5.124692007057e-07 true resid norm 5.113062797394e-07 ||r(i)||/||b|| 6.182559248811e-06
    8 KSP unpreconditioned resid norm 9.477347074078e-08 true resid norm 9.396189342400e-08 ||r(i)||/||b|| 1.136158494905e-06
    9 KSP unpreconditioned resid norm 8.091661284678e-09 true resid norm 1.202285141644e-08 ||r(i)||/||b|| 1.453766444247e-07
      Line search: gnorm after quadratic fit 7.445772470565e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 1 Nonlinear |R| = 7.445772e-02
    0 KSP unpreconditioned resid norm 7.445772470565e-02 true resid norm 7.445772470565e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 1.021824982462e-02 true resid norm 1.021824985160e-02 ||r(i)||/||b|| 1.372355909611e-01
    2 KSP unpreconditioned resid norm 1.157885671391e-03 true resid norm 1.157887754391e-03 ||r(i)||/||b|| 1.555094194684e-02
    3 KSP unpreconditioned resid norm 1.430130714337e-04 true resid norm 1.430161577785e-04 ||r(i)||/||b|| 1.920769918016e-03
    4 KSP unpreconditioned resid norm 2.998894493294e-05 true resid norm 2.998824778235e-05 ||r(i)||/||b|| 4.027553608561e-04
    5 KSP unpreconditioned resid norm 4.649279695907e-06 true resid norm 4.648536454806e-06 ||r(i)||/||b|| 6.243188968213e-05
    6 KSP unpreconditioned resid norm 5.924328672653e-07 true resid norm 5.931708890293e-07 ||r(i)||/||b|| 7.966546001428e-06
    7 KSP unpreconditioned resid norm 8.543610682426e-08 true resid norm 8.697539896531e-08 ||r(i)||/||b|| 1.168117872379e-06
    8 KSP unpreconditioned resid norm 1.148676220866e-08 true resid norm 1.533822422448e-08 ||r(i)||/||b|| 2.059990992891e-07
      Line search: gnorm after quadratic fit 6.366652923580e-02
      Line search: Quadratically determined step, lambda=1.8856596447436538e-01
 2 Nonlinear |R| = 6.366653e-02
    0 KSP unpreconditioned resid norm 6.366652923580e-02 true resid norm 6.366652923580e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 8.132883645367e-03 true resid norm 8.132882573027e-03 ||r(i)||/||b|| 1.277418868383e-01
    2 KSP unpreconditioned resid norm 9.050257866029e-04 true resid norm 9.050257068383e-04 ||r(i)||/||b|| 1.421509414290e-02
    3 KSP unpreconditioned resid norm 1.947546817982e-04 true resid norm 1.947541659447e-04 ||r(i)||/||b|| 3.058972560345e-03
    4 KSP unpreconditioned resid norm 2.443905818014e-05 true resid norm 2.443865117433e-05 ||r(i)||/||b|| 3.838539883935e-04
    5 KSP unpreconditioned resid norm 3.078035636215e-06 true resid norm 3.076235340407e-06 ||r(i)||/||b|| 4.831793687094e-05
    6 KSP unpreconditioned resid norm 4.089002367993e-07 true resid norm 4.099035071354e-07 ||r(i)||/||b|| 6.438288878875e-06
    7 KSP unpreconditioned resid norm 4.375575845180e-08 true resid norm 5.030349232605e-08 ||r(i)||/||b|| 7.901089148389e-07
      Line search: Using full step: fnorm 6.366652923580e-02 gnorm 1.747232935912e-02
 3 Nonlinear |R| = 1.747233e-02
    0 KSP unpreconditioned resid norm 1.747232935912e-02 true resid norm 1.747232935912e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 4.126580102545e-04 true resid norm 4.126580151809e-04 ||r(i)||/||b|| 2.361780199419e-02
    2 KSP unpreconditioned resid norm 7.833804841745e-05 true resid norm 7.833800270617e-05 ||r(i)||/||b|| 4.483546589355e-03
    3 KSP unpreconditioned resid norm 7.041000286198e-06 true resid norm 7.041122450409e-06 ||r(i)||/||b|| 4.029870491615e-04
    4 KSP unpreconditioned resid norm 1.799934492270e-07 true resid norm 1.799839999998e-07 ||r(i)||/||b|| 1.030108786874e-05
    5 KSP unpreconditioned resid norm 4.747962252506e-09 true resid norm 4.792592460950e-09 ||r(i)||/||b|| 2.742961377642e-07
      Line search: gnorm after quadratic fit 1.575853046178e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 4 Nonlinear |R| = 1.575853e-02
    0 KSP unpreconditioned resid norm 1.575853046178e-02 true resid norm 1.575853046178e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.686254574996e-04 true resid norm 3.686253665087e-04 ||r(i)||/||b|| 2.339211561653e-02
    2 KSP unpreconditioned resid norm 7.058799746390e-05 true resid norm 7.058832668696e-05 ||r(i)||/||b|| 4.479372417254e-03
    3 KSP unpreconditioned resid norm 4.168887297993e-06 true resid norm 4.168998306779e-06 ||r(i)||/||b|| 2.645550177975e-04
    4 KSP unpreconditioned resid norm 1.835421781029e-07 true resid norm 1.838094501080e-07 ||r(i)||/||b|| 1.166412379338e-05
    5 KSP unpreconditioned resid norm 5.654125280681e-09 true resid norm 5.689908572853e-09 ||r(i)||/||b|| 3.610684756838e-07
      Line search: gnorm after quadratic fit 1.418047944959e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 5 Nonlinear |R| = 1.418048e-02
    0 KSP unpreconditioned resid norm 1.418047944959e-02 true resid norm 1.418047944959e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.069883562322e-04 true resid norm 3.069883295194e-04 ||r(i)||/||b|| 2.164865656417e-02
    2 KSP unpreconditioned resid norm 5.295387781758e-05 true resid norm 5.295382274553e-05 ||r(i)||/||b|| 3.734275906099e-03
    3 KSP unpreconditioned resid norm 2.850992985264e-06 true resid norm 2.850921032684e-06 ||r(i)||/||b|| 2.010454613202e-04
    4 KSP unpreconditioned resid norm 1.200519331392e-07 true resid norm 1.198826381452e-07 ||r(i)||/||b|| 8.454060990770e-06
    5 KSP unpreconditioned resid norm 3.383567642095e-09 true resid norm 3.581361206780e-09 ||r(i)||/||b|| 2.525557206660e-07
      Line search: gnorm after quadratic fit 1.276114941357e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 6 Nonlinear |R| = 1.276115e-02
    0 KSP unpreconditioned resid norm 1.276114941357e-02 true resid norm 1.276114941357e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.563057825506e-04 true resid norm 2.563058546908e-04 ||r(i)||/||b|| 2.008485649562e-02
    2 KSP unpreconditioned resid norm 3.994442443128e-05 true resid norm 3.994433312699e-05 ||r(i)||/||b|| 3.130151668352e-03
    3 KSP unpreconditioned resid norm 1.962686750006e-06 true resid norm 1.962712191215e-06 ||r(i)||/||b|| 1.538037152929e-04
    4 KSP unpreconditioned resid norm 7.935745397693e-08 true resid norm 7.938533180401e-08 ||r(i)||/||b|| 6.220860616174e-06
    5 KSP unpreconditioned resid norm 2.045989528038e-09 true resid norm 2.283380723088e-09 ||r(i)||/||b|| 1.789322144179e-07
      Line search: gnorm after quadratic fit 1.148434950579e-02
      Line search: Quadratically determined step, lambda=1.0000000000000001e-01
 7 Nonlinear |R| = 1.148435e-02
    0 KSP unpreconditioned resid norm 1.148434950579e-02 true resid norm 1.148434950579e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 2.145284109969e-04 true resid norm 2.145283878458e-04 ||r(i)||/||b|| 1.868006435520e-02
    2 KSP unpreconditioned resid norm 3.030013521147e-05 true resid norm 3.030015948939e-05 ||r(i)||/||b|| 2.638387091417e-03
    3 KSP unpreconditioned resid norm 1.360293150514e-06 true resid norm 1.360270262454e-06 ||r(i)||/||b|| 1.184455647025e-04
    4 KSP unpreconditioned resid norm 5.301317978538e-08 true resid norm 5.301854728179e-08 ||r(i)||/||b|| 4.616591236192e-06
    5 KSP unpreconditioned resid norm 1.250674455488e-09 true resid norm 1.346415141119e-09 ||r(i)||/||b|| 1.172391296904e-07
      Line search: gnorm after quadratic fit 1.221312038047e-02
      Line search: Cubically determined step, current gnorm 1.113904272651e-02 lambda=3.1282174556276493e-02
 8 Nonlinear |R| = 1.113904e-02
    0 KSP unpreconditioned resid norm 1.113904272651e-02 true resid norm 1.113904272651e-02 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 6.001129085887e-04 true resid norm 6.001131867060e-04 ||r(i)||/||b|| 5.387475400179e-02
    2 KSP unpreconditioned resid norm 1.591024896256e-05 true resid norm 1.591018669640e-05 ||r(i)||/||b|| 1.428326211420e-03
    3 KSP unpreconditioned resid norm 3.537121136724e-07 true resid norm 3.536442606569e-07 ||r(i)||/||b|| 3.174817345977e-05
    4 KSP unpreconditioned resid norm 2.481134919479e-08 true resid norm 2.495824502476e-08 ||r(i)||/||b|| 2.240609506360e-06
    5 KSP unpreconditioned resid norm 6.280939831947e-10 true resid norm 7.825230714614e-10 ||r(i)||/||b|| 7.025047759259e-08
      Line search: Using full step: fnorm 1.113904272651e-02 gnorm 3.585939794007e-04
 9 Nonlinear |R| = 3.585940e-04
    0 KSP unpreconditioned resid norm 3.585939794007e-04 true resid norm 3.585939794007e-04 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 9.330020496656e-06 true resid norm 9.330024784268e-06 ||r(i)||/||b|| 2.601835312422e-02
    2 KSP unpreconditioned resid norm 3.926801690618e-08 true resid norm 3.928161648344e-08 ||r(i)||/||b|| 1.095434355844e-04
    3 KSP unpreconditioned resid norm 3.065503041037e-10 true resid norm 3.098511110230e-10 ||r(i)||/||b|| 8.640722622862e-07
      Line search: Using full step: fnorm 3.585939794007e-04 gnorm 3.265151694809e-06
10 Nonlinear |R| = 3.265152e-06
    0 KSP unpreconditioned resid norm 3.265151694809e-06 true resid norm 3.265151694809e-06 ||r(i)||/||b|| 1.000000000000e+00
    1 KSP unpreconditioned resid norm 3.359141488324e-08 true resid norm 3.359147158166e-08 ||r(i)||/||b|| 1.028787472112e-02
    2 KSP unpreconditioned resid norm 1.124580382469e-10 true resid norm 1.124766996610e-10 ||r(i)||/||b|| 3.444761841841e-05
    3 KSP unpreconditioned resid norm 1.702558333457e-12 true resid norm 1.720775683673e-12 ||r(i)||/||b|| 5.270124773708e-07
      Line search: Using full step: fnorm 3.265151694809e-06 gnorm 7.965803178017e-11
11 Nonlinear |R| = 7.965803e-11
 Solve Converged!
