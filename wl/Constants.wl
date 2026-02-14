(* Constants core package (P,N,K) harness *)

BeginPackage["Constants`"];

BesselBlock::usage = "BesselBlock[p,k] gives J_p(pi k/2) p! (4/(pi k))^p with memoization.";
Heartbeat::usage = "Heartbeat[msg] prints a heartbeat status line at most once every $ConstantsHeartbeatPeriod seconds.";
ConstantsSetPrecision::usage = "ConstantsSetPrecision[prec] sets working precision for computations.";
ConstantsSetParameters::usage = "ConstantsSetParameters[{P,N,K}] sets default parameters.";
Sk::usage = "Sk[a,k] = sum_{p=0}^{P-1} a_p * mathcal J(p,k).";
SkSeriesCoefficients::usage = "SkSeriesCoefficients[a,r,K] gives coefficients s_t where Sk ~ k^(-1/2) Sum_{t=0}^K s_t k^-t for k congruent r mod 4.";
TailObjective::usage = "TailObjective[a,N,K] computes accelerated Kummer tail for objective.";
Objective::usage = "Objective[a,N,K] computes C(a)=1/2+Sum_{k>=1} Sk[a,k]^4 with mod-4 accelerated tail.";

Begin["`Private`"];

If[!ValueQ[$ConstantsHeartbeatPeriod], $ConstantsHeartbeatPeriod = 30;];
If[!ValueQ[$ConstantsLastHeartbeat], $ConstantsLastHeartbeat = -Infinity;];
If[!ValueQ[$ConstantsWorkingPrecision], $ConstantsWorkingPrecision = 120;];
If[!ValueQ[$ConstantsP], $ConstantsP = 8;];
If[!ValueQ[$ConstantsN], $ConstantsN = 20;];
If[!ValueQ[$ConstantsK], $ConstantsK = 4;];

ConstantsSetPrecision[prec_Integer?Positive] := ($ConstantsWorkingPrecision = prec);
ConstantsSetParameters[{p_Integer?Positive, n_Integer?Positive, k_Integer?NonNegative}] := ({ $ConstantsP,$ConstantsN,$ConstantsK} = {p,n,k});

ConstantsLog[msg_] := Print[DateString[{"TimeShort"}], " | ", msg];

Heartbeat[msg_String:"still running"] := Module[{now = AbsoluteTime[]},
  If[now - $ConstantsLastHeartbeat >= $ConstantsHeartbeatPeriod,
    ConstantsLog[msg];
    $ConstantsLastHeartbeat = now;
  ];
];

ClearAll[BesselBlock];
BesselBlock[p_Integer?NonNegative, k_Integer?Positive] := Module[{prec = $ConstantsWorkingPrecision, z},
  z = N[Pi*k/2, prec];
  BesselBlock[p,k] = N[BesselJ[p, z] * N[Factorial[p], prec] * (N[4/(Pi*k), prec])^p, prec]
];

ClearAll[AsymptCoeff];
AsymptCoeff[p_Integer?NonNegative, n_Integer?NonNegative] := AsymptCoeff[p,n] =
  If[n==0, 1, Product[4 p^2 - (2 j - 1)^2, {j,1,n}]/(8^n n!)];

CosTheta[r_Integer?NonNegative, p_Integer?NonNegative] := CosTheta[r,p] =
  Cos[Pi r/2] Cos[(2 p + 1) Pi/4] + Sin[Pi r/2] Sin[(2 p + 1) Pi/4];
SinTheta[r_Integer?NonNegative, p_Integer?NonNegative] := SinTheta[r,p] =
  Sin[Pi r/2] Cos[(2 p + 1) Pi/4] - Cos[Pi r/2] Sin[(2 p + 1) Pi/4];
PrefactorExact[p_Integer?NonNegative] := (2 Factorial[p] 4^p)/Pi^(p+1);

SkSeriesCoefficients[a_List, r_Integer?NonNegative, K_Integer?NonNegative] := Module[
  {prec = $ConstantsWorkingPrecision, twoOverPi, pow, s, p, ai},
  twoOverPi = N[2/Pi, prec];
  pow[n_Integer?NonNegative] := pow[n] = twoOverPi^n;
  s = ConstantArray[0, K+1];
  Do[
    p = i-1;
    ai = N[a[[i]], prec];
    If[p <= K,
      Do[
        t = p + 2 m;
        s[[t+1]] += ai * N[PrefactorExact[p]*CosTheta[r,p]*(-1)^m*AsymptCoeff[p,2m]*pow[2m], prec],
        {m, 0, Floor[(K - p)/2]}
      ];
      Do[
        t = p + (2 m + 1);
        s[[t+1]] += ai * N[PrefactorExact[p]*(-SinTheta[r,p])*(-1)^m*AsymptCoeff[p,2m+1]*pow[2m+1], prec],
        {m, 0, Floor[(K - p - 1)/2]}
      ];
    ],
    {i, 1, Length[a]}
  ];
  N[s, prec]
];

Sk[a_List, k_Integer?Positive] := Module[{prec = $ConstantsWorkingPrecision},
  Sum[N[a[[p+1]], prec] * BesselBlock[p, k], {p, 0, Length[a]-1}]
];

ObjectiveTailCoeff[s_List, K_Integer?NonNegative] := Module[{x, poly},
  poly = Normal@Series[Sum[s[[t+1]] x^t, {t,0,K}]^4, {x,0,K}];
  Table[Coefficient[poly, x, u], {u,0,K}]
];

ClearAll[TailObjective];
TailObjective[a_List, N_Integer?Positive, K_Integer?NonNegative] := Module[
  {prec = $ConstantsWorkingPrecision, sum = 0, r, s, c, u, sExp, first, n0},
  Do[
    s = SkSeriesCoefficients[a, r, K];
    c = ObjectiveTailCoeff[s, K];
    Do[
      sExp = 2 + u;
      first = N + 1;
      n0 = Ceiling[(first - r)/4];
      sum += N[c[[u+1]] * 4^(-sExp) * HurwitzZeta[sExp, n0 + r/4], prec],
      {u, 0, K}
    ],
    {r, 0, 3}
  ];
  N[sum, prec]
];

Objective[a_List, N_Integer:$ConstantsN, K_Integer:$ConstantsK] := Module[
  {prec = $ConstantsWorkingPrecision, Cfinite, Ctail},
  Cfinite = N[1/2 + Sum[Sk[a,k]^4, {k,1,N}], prec];
  Ctail = TailObjective[a, N, K];
  Cfinite + Ctail
];

End[];
EndPackage[];
