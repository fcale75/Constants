BeginPackage["Constants`Newton`"]
NewtonObjective::usage="NewtonObjective[a,pmax,nmax] computes the truncated objective sum (finite k).";
NewtonDerivatives::usage="NewtonDerivatives[a,pmax,nmax,lambda] returns a report association with objective, gpart, Hessian and residual.";
NewtonStep::usage="NewtonStep[a,lambda,pmax,nmax] performs one Newton step on the KKT system (constraint Sum[a]==1).";

Begin["`Private`"]
Needs["Constants`"];

(* Truncated objective (no tail) *)
NewtonObjective[a_List,pmax_Integer,nmax_Integer] := Module[{S,k},
  S[k_] := Sum[a[[j+1]] BesselBlock[j,k], {j,0,pmax-1}];
  1/2 + Sum[(S[k])^4, {k,1,nmax}]
];

NewtonDerivatives[a_List,pmax_Integer,nmax_Integer,lambda_] := Module[{S,k,gpart,H},
  S[k_] := Sum[a[[j+1]] BesselBlock[j,k], {j,0,pmax-1}];
  gpart = Table[4*Sum[BesselBlock[l,k]*S[k]^3, {k,1,nmax}], {l,0,pmax-1}];
  H = Table[12*Sum[BesselBlock[i,k]*BesselBlock[l,k]*S[k]^2, {k,1,nmax}], {i,0,pmax-1}, {l,0,pmax-1}];
  <|
    "objective" -> (1/2 + Sum[(S[k])^4, {k,1,nmax}]),
    "gpart" -> gpart,
    "hess" -> H,
    "res" -> (gpart - lambda)
  |>
];

NewtonStep[a_List, lambda_, pmax_Integer, nmax_Integer] := Module[{d,ones,M,rhs,sol,da,dl},
  d = NewtonDerivatives[a,pmax,nmax,lambda];
  ones = ConstantArray[1,pmax];
  M = ArrayFlatten[{{d["hess"], -Transpose@{ones}}, {-{ones}, {{0}}}}];
  rhs = Join[-d["res"], {1-Total[a]}];
  sol = LinearSolve[M, rhs];
  da = sol[[;;pmax]]; dl = sol[[pmax+1]];
  <|"deltaA"->da,"deltaLambda"->dl,"aNext"->a+da,"lambdaNext"->lambda+dl|>
];

End[];
EndPackage[];
