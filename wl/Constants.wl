(* ::Package:: *)

(* Constants.wl
   Core definitions for the Rechnitzer-style (P,N,K) ansatz.

   Public API:
     - ConstantsSetPrecision[prec]
     - ConstantsSetParameters[{P,N,K}]
     - ConstantsClearCaches[]
     - Heartbeat[msg]
     - BesselBlock[p,k]           (mathcal J(p,k))
     - Sk[a,k]
     - SkSeriesCoefficients[a,r,K]
     - TailObjective[a,N,K]
     - Objective[a,N,K]
*)

BeginPackage["Constants`"];

BesselBlock::usage =
  "BesselBlock[p,k] gives J_p(pi k/2) * p! * (4/(pi k))^p (the paper's \\[ScriptJ](p,k)). \
Values are cached by (p,k,workingPrecision).";

Heartbeat::usage =
  "Heartbeat[msg] prints a heartbeat line at most once every $ConstantsHeartbeatPeriod seconds.";

ConstantsSetPrecision::usage =
  "ConstantsSetPrecision[prec] sets the working precision used by this package.";

ConstantsSetParameters::usage =
  "ConstantsSetParameters[{P,N,K}] sets default parameters for Objective/series computations.";

ConstantsClearCaches::usage =
  "ConstantsClearCaches[] clears memoized caches (e.g. BesselBlock).";

Sk::usage =
  "Sk[a,k] = Sum_{p=0}^{Length[a]-1} a_p * BesselBlock[p,k].";

SkSeriesCoefficients::usage =
  "SkSeriesCoefficients[a,r,K] returns coefficients s_t (t=0..K) such that \
Sk[a,k] ~ k^(-1/2) * Sum_{t=0}^K s_t k^(-t) for k congruent to r (mod 4).";

TailObjective::usage =
  "TailObjective[a,N,K] computes the accelerated mod-4 Kummer tail for Sum_{k>N} Sk[a,k]^4.";

Objective::usage =
  "Objective[a,N,K] computes C(a) = 1/2 + Sum_{k>=1} Sk[a,k]^4 using \
a finite sum to N plus the mod-4 accelerated tail.";

Begin["`Private`"];

(* ------------------------------------------------------------------------- *)
(* Defaults / globals                                                         *)
(* ------------------------------------------------------------------------- *)

If[!ValueQ[$ConstantsHeartbeatPeriod], $ConstantsHeartbeatPeriod = 30];
If[!ValueQ[$ConstantsLastHeartbeat], $ConstantsLastHeartbeat = -Infinity];
If[!ValueQ[$ConstantsWorkingPrecision], $ConstantsWorkingPrecision = 120];

If[!ValueQ[$ConstantsP], $ConstantsP = 8];
If[!ValueQ[$ConstantsN], $ConstantsN = 20];
If[!ValueQ[$ConstantsK], $ConstantsK = 4];

If[!ValueQ[$ConstantsBesselBlockCache], $ConstantsBesselBlockCache = <||>];

ConstantsSetPrecision[prec_Integer?Positive] := Module[{},
  $ConstantsWorkingPrecision = prec;
  Null
];

ConstantsSetParameters[{p_Integer?Positive, n_Integer?Positive, k_Integer?NonNegative}] := Module[{},
  {$ConstantsP, $ConstantsN, $ConstantsK} = {p, n, k};
  Null
];

ConstantsClearCaches[] := Module[{},
  $ConstantsBesselBlockCache = <||>;
  ClearAll[AsymptCoeff, CosTheta, SinTheta];
  Null
];

ConstantsLog[msg_] := Print[DateString[{"TimeShort"}], " | ", msg];

Heartbeat[msg_String : "still running"] := Module[{now = AbsoluteTime[]},
  If[now - $ConstantsLastHeartbeat >= $ConstantsHeartbeatPeriod,
    ConstantsLog[msg];
    $ConstantsLastHeartbeat = now;
  ];
];

(* ------------------------------------------------------------------------- *)
(* BesselBlock = mathcal J(p,k)                                               *)
(* ------------------------------------------------------------------------- *)

ClearAll[BesselBlock];

BesselBlock[p_Integer?NonNegative, k_Integer?Positive] := Module[
  {prec = $ConstantsWorkingPrecision, key, cached, z, val},
  key = {p, k, prec};
  cached = Lookup[$ConstantsBesselBlockCache, key, Missing["NotFound"]];
  If[cached =!= Missing["NotFound"], Return[cached]];

  z = SetPrecision[Pi*k/2, prec];
  val = N[BesselJ[p, z] * Factorial[p] * (4/(Pi*k))^p, prec];

  $ConstantsBesselBlockCache[key] = val;
  val
];

(* ------------------------------------------------------------------------- *)
(* Asymptotic machinery for J(p,k) as k -> Infinity (mod 4 splitting)         *)
(* ------------------------------------------------------------------------- *)

ClearAll[AsymptCoeff];

(* a_n(p) = Product_{j=1..n} (4 p^2 - (2j-1)^2) / (8^n n!)  *)
AsymptCoeff[p_Integer?NonNegative, n_Integer?NonNegative] :=
  AsymptCoeff[p, n] =
    If[n == 0,
      1,
      Product[4 p^2 - (2 j - 1)^2, {j, 1, n}] / (8^n * n!)
    ];

(* For k = 4n + r, reduce theta = pi k/2 - p pi/2 - pi/4 mod 2pi. *)
ClearAll[CosTheta, SinTheta];

CosTheta[r_Integer?NonNegative, p_Integer?NonNegative] :=
  CosTheta[r, p] =
    Cos[Pi*r/2]*Cos[(2 p + 1) Pi/4] + Sin[Pi*r/2]*Sin[(2 p + 1) Pi/4];

SinTheta[r_Integer?NonNegative, p_Integer?NonNegative] :=
  SinTheta[r, p] =
    Sin[Pi*r/2]*Cos[(2 p + 1) Pi/4] - Cos[Pi*r/2]*Sin[(2 p + 1) Pi/4];

(* Prefactor for BesselBlock(p,k) after pulling out k^{-p-1/2} *)
ClearAll[PrefactorExact];
PrefactorExact[p_Integer?NonNegative] := (2 * Factorial[p] * 4^p) / Pi^(p + 1);

ClearAll[SkSeriesCoefficients];

SkSeriesCoefficients[a_List, r_Integer?NonNegative, K_Integer?NonNegative] := Module[
  {prec = $ConstantsWorkingPrecision, invPi, pow, s, i, p, ai, m, t},

  (* Important: for z = pi k/2, (2 z) = pi k, so powers are (pi k)^{-n}. *)
  invPi = N[1/Pi, prec];
  pow[n_Integer?NonNegative] := pow[n] = invPi^n;

  s = ConstantArray[0, K + 1];

  Do[
    p = i - 1;
    ai = N[a[[i]], prec];

    If[p <= K,
      (* even n = 2m terms: cos(theta) piece *)
      Do[
        t = p + 2*m;
        s[[t + 1]] += ai * N[
          PrefactorExact[p] * CosTheta[r, p] * (-1)^m *
            AsymptCoeff[p, 2*m] * pow[2*m],
          prec
        ];
      , {m, 0, Floor[(K - p)/2]}];

      (* odd n = 2m+1 terms: -sin(theta) piece *)
      Do[
        t = p + (2*m + 1);
        s[[t + 1]] += ai * N[
          PrefactorExact[p] * (-SinTheta[r, p]) * (-1)^m *
            AsymptCoeff[p, 2*m + 1] * pow[2*m + 1],
          prec
        ];
      , {m, 0, Floor[(K - p - 1)/2]}];
    ];
  , {i, 1, Length[a]}];

  N[s, prec]
];

(* ------------------------------------------------------------------------- *)
(* Objective: finite sum + mod-4 accelerated tail                              *)
(* ------------------------------------------------------------------------- *)

ClearAll[Sk];

Sk[a_List, k_Integer?Positive] := Module[{prec = $ConstantsWorkingPrecision},
  Sum[N[a[[p + 1]], prec] * BesselBlock[p, k], {p, 0, Length[a] - 1}]
];

ClearAll[ObjectiveTailCoeff];

ObjectiveTailCoeff[s_List, K_Integer?NonNegative] := Module[{x, poly},
  poly = Normal @ Series[(Sum[s[[t + 1]] * x^t, {t, 0, K}])^4, {x, 0, K}];
  Table[Coefficient[poly, x, u], {u, 0, K}]
];

ClearAll[TailObjective];

TailObjective[a_List, N_Integer?Positive, K_Integer?NonNegative] := Module[
  {prec = $ConstantsWorkingPrecision, sum = 0, r, s, c, u, sExp, first, n0},

  Do[
    s = SkSeriesCoefficients[a, r, K];
    c = ObjectiveTailCoeff[s, K];

    Do[
      sExp = 2 + u;                  (* because Sk^4 has leading k^{-2} *)
      first = N + 1;                 (* first index strictly greater than N *)
      n0 = Ceiling[(first - r)/4];   (* smallest n with k = 4n+r > N *)

      sum += N[
        c[[u + 1]] * 4^(-sExp) * HurwitzZeta[sExp, n0 + r/4],
        prec
      ];
    , {u, 0, K}];
  , {r, 0, 3}];

  N[sum, prec]
];

ClearAll[Objective];

Objective[a_List, N_Integer : $ConstantsN, K_Integer : $ConstantsK] := Module[
  {prec = $ConstantsWorkingPrecision, Cfinite, Ctail},
  Cfinite = N[1/2 + Sum[Sk[a, k]^4, {k, 1, N}], prec];
  Ctail = TailObjective[a, N, K];
  Cfinite + Ctail
];

End[];
EndPackage[];

