(* ::Package:: *)

(* Newton.wl
   Finite-sum (k <= N) objective/derivatives and one KKT Newton step.
   This module intentionally omits the Kummer tail. The Driver can add tail
   contributions via Constants`TailDerivatives`.
*)

BeginPackage["Constants`Newton`"];

NewtonObjective::usage =
  "NewtonObjective[a,P,N] computes the truncated objective 1/2 + Sum_{k=1..N} Sk[a,k]^4 \
(without the Kummer tail).";

NewtonDerivatives::usage =
  "NewtonDerivatives[a,P,N,lambda] returns an association with keys: \
\"objective\", \"gpart\" (gradient of truncated objective), \"hess\", \"res\" (= gpart - lambda), \
and \"constraint\" (= 1-Total[a]).";

NewtonStep::usage =
  "NewtonStep[a,lambda,P,N] performs one Newton step for the truncated KKT system \
(gradient + constraint Sum[a]==1).";

Begin["`Private`"];

Needs["Constants`"];

ClearAll[NewtonObjective];

NewtonObjective[a_List, P_Integer?Positive, N_Integer?Positive] := Module[{S, k},
  S[k_] := Sum[a[[j + 1]] * BesselBlock[j, k], {j, 0, P - 1}];
  1/2 + Sum[S[k]^4, {k, 1, N}]
];

ClearAll[NewtonDerivatives];

NewtonDerivatives[a_List, P_Integer?Positive, N_Integer?Positive, lambda_] := Module[
  {S, k, gpart, H, obj, constraint},

  S[k_] := Sum[a[[j + 1]] * BesselBlock[j, k], {j, 0, P - 1}];

  obj = 1/2 + Sum[S[k]^4, {k, 1, N}];

  gpart = Table[
    4 * Sum[BesselBlock[l, k] * S[k]^3, {k, 1, N}],
    {l, 0, P - 1}
  ];

  H = Table[
    12 * Sum[BesselBlock[i, k] * BesselBlock[l, k] * S[k]^2, {k, 1, N}],
    {i, 0, P - 1}, {l, 0, P - 1}
  ];

  constraint = 1 - Total[a];

  <|
    "objective" -> obj,
    "gpart" -> gpart,
    "hess" -> H,
    "res" -> (gpart - lambda),
    "constraint" -> constraint
  |>
];

ClearAll[NewtonStep];

NewtonStep[a_List, lambda_, P_Integer?Positive, N_Integer?Positive] := Module[
  {d, ones, M, rhs, sol, da, dl},

  d = NewtonDerivatives[a, P, N, lambda];
  ones = ConstantArray[1, P];

  (* Jacobian of F(a,lambda) = { g(a) - lambda*1 , 1-Total[a] } *)
  M = ArrayFlatten[{
    {d["hess"], -Transpose@{ones}},
    {-{ones}, {{0}}}
  }];

  rhs = Join[-d["res"], {-(d["constraint"])}];

  sol = LinearSolve[M, rhs];

  da = sol[[;; P]];
  dl = sol[[P + 1]];

  <|
    "deltaA" -> da,
    "deltaLambda" -> dl,
    "aNext" -> (a + da),
    "lambdaNext" -> (lambda + dl),
    "report" -> d
  |>
];

End[];
EndPackage[];

