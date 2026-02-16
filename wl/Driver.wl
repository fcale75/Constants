(* ::Package:: *)

(* Driver.wl
   Paper-faithful constrained Newton iteration for the Rechnitzer ansatz.

   Key fixes relative to earlier versions:
     - Uses the KKT residual (g(a) - lambda*1) as the stationarity test.
     - Optionally includes the mod-4 Kummer tail in objective/gradient/Hessian.
     - Includes a simple backtracking damping step for stability.
*)

BeginPackage["Constants`Driver`"];

NewtonOptimize::usage =
  "NewtonOptimize[P,Nmax,opts] runs a constrained Newton method for minimizing \
C(a)=1/2+Sum_{k>=1} Sk[a,k]^4 with the constraint Total[a]==1. \
By default it uses the mod-4 accelerated tail (UseTail->True).";

Begin["`Private`"];

Needs["Constants`"];
Needs["Constants`Newton`"];
Needs["Constants`TailDerivatives`"];

Options[NewtonOptimize] = {
  MaxIterations -> 25,
  Tolerance -> 10^-25,
  ConstraintTolerance -> 10^-25,
  Log -> True,

  WorkingPrecision -> Automatic,
  K -> Automatic,
  UseTail -> True,

  Damping -> True,
  StepMin -> 2^-20,

  InitialA -> Automatic,
  InitialLambda -> Automatic,
  ReturnHistory -> False
};

(* ------------------------------------------------------------------------- *)
(* Helpers                                                                    *)
(* ------------------------------------------------------------------------- *)

ClearAll[fullDerivatives];
fullDerivatives[a_List, P_Integer?Positive, Nmax_Integer?Positive, K_Integer?NonNegative, lambda_] := Module[
  {prec = Constants`Private`$ConstantsWorkingPrecision, rep, objTail, gTail, hTail},
  rep = Constants`Newton`NewtonDerivatives[a, P, Nmax, lambda];

  objTail = Constants`TailObjective[a, Nmax, K];
  gTail = Constants`TailDerivatives`TailGradient[a, Nmax, K];
  hTail = Constants`TailDerivatives`TailHessian[a, Nmax, K];

  rep["objective"] = N[rep["objective"] + objTail, prec];
  rep["gpart"] = N[rep["gpart"] + gTail, prec];
  rep["hess"] = N[rep["hess"] + hTail, prec];
  rep["res"] = rep["gpart"] - lambda; (* stationarity residual *)

  rep
];

ClearAll[kktStepFromReport];
kktStepFromReport[rep_Association, a_List] := Module[
  {P = Length[a], ones, M, rhs, sol, da, dl},
  ones = ConstantArray[1, P];

  M = ArrayFlatten[{
    {rep["hess"], -Transpose@{ones}},
    {-{ones}, {{0}}}
  }];

  rhs = Join[-rep["res"], {-(rep["constraint"])}];
  sol = LinearSolve[M, rhs];

  da = sol[[;; P]];
  dl = sol[[P + 1]];

  <|"deltaA" -> da, "deltaLambda" -> dl|>
];

ClearAll[currentObjective];
currentObjective[a_List, P_Integer?Positive, Nmax_Integer?Positive, K_Integer?NonNegative, useTail_] :=
  If[TrueQ[useTail],
    Constants`Objective[a, Nmax, K],
    Constants`Newton`NewtonObjective[a, P, Nmax]
  ];

(* ------------------------------------------------------------------------- *)
(* Main driver                                                                *)
(* ------------------------------------------------------------------------- *)

NewtonOptimize[P_Integer?Positive, Nmax_Integer?Positive, opts : OptionsPattern[]] := Module[
  {maxIt = OptionValue[MaxIterations],
   tol = OptionValue[Tolerance],
   cTol = OptionValue[ConstraintTolerance],
   logQ = TrueQ[OptionValue[Log]],
   precOpt = OptionValue[WorkingPrecision],
   Kopt = OptionValue[K],
   useTail = TrueQ[OptionValue[UseTail]],
   damping = TrueQ[OptionValue[Damping]],
   stepMin = OptionValue[StepMin],
   a0 = OptionValue[InitialA],
   lambda0 = OptionValue[InitialLambda],
   returnHistory = TrueQ[OptionValue[ReturnHistory]],

   Kval, a, lambda, it,
   rep, step, da, dl, t, obj0, objTry,
   resNorm, conNorm, history, aTry, lambdaTry},

  If[precOpt =!= Automatic, ConstantsSetPrecision[precOpt]];
  Kval = If[Kopt === Automatic, Constants`Private`$ConstantsK, Kopt];

  a = If[a0 === Automatic,
    Normalize[Table[(-1)^n/(8^n*(n + 1)), {n, 0, P - 1}], Total],
    a0
  ];

  (* Enforce the constraint at the start (and only at the start). *)
  a = a - (Total[a] - 1)/P * ConstantArray[1, P];

  (* Initialize lambda sensibly: mean gradient makes res have mean ~ 0. *)
  If[lambda0 === Automatic,
    rep = If[useTail,
      fullDerivatives[a, P, Nmax, Kval, 0],
      Constants`Newton`NewtonDerivatives[a, P, Nmax, 0]
    ];
    lambda = Mean[rep["gpart"]],
    lambda = lambda0
  ];

  history = {};

  For[it = 1, it <= maxIt, it++,
    Heartbeat["Newton iter " <> ToString[it]];

    rep = If[useTail,
      fullDerivatives[a, P, Nmax, Kval, lambda],
      Constants`Newton`NewtonDerivatives[a, P, Nmax, lambda]
    ];

    obj0 = rep["objective"];
    resNorm = Norm[rep["res"], Infinity];
    conNorm = Abs[rep["constraint"]];

    If[logQ,
      Print[
        "iter ", it,
        "  C=", NumberForm[obj0, {Infinity, 16}],
        "  ||g-lambda||_inf=", ScientificForm[resNorm],
        "  |1-Total[a]|=", ScientificForm[conNorm]
      ];
    ];

    If[returnHistory,
      AppendTo[history, <|
        "iter" -> it, "objective" -> obj0,
        "resNorm" -> resNorm, "constraintNorm" -> conNorm
      |>];
    ];

    If[resNorm < tol && conNorm < cTol,
      If[logQ, Print["converged"]];
      Return[<|
        "aFinal" -> a,
        "lambdaFinal" -> lambda,
        "objective" -> obj0,
        "report" -> rep,
        "iterations" -> it,
        "history" -> If[returnHistory, history, Missing["NotRequested"]],
        "params" -> <|"P" -> P, "N" -> Nmax, "K" -> Kval, "useTail" -> useTail|>
      |>];
    ];

    step = kktStepFromReport[rep, a];
    da = step["deltaA"];
    dl = step["deltaLambda"];

    t = 1;

    If[damping,
      While[t >= stepMin,
        aTry = a + t*da;
        lambdaTry = lambda + t*dl;
        objTry = currentObjective[aTry, P, Nmax, Kval, useTail];

        If[objTry < obj0, Break[]];
        t = t/2;
      ];

      If[t < stepMin, t = stepMin];

      If[logQ && t < 1,
        Print["  damping: step factor t=", t];
      ];
    ];

    a = a + t*da;
    lambda = lambda + t*dl;
  ];

  (* Return best available after max iterations. *)
  rep = If[useTail,
    fullDerivatives[a, P, Nmax, Kval, lambda],
    Constants`Newton`NewtonDerivatives[a, P, Nmax, lambda]
  ];

  <|
    "aFinal" -> a,
    "lambdaFinal" -> lambda,
    "objective" -> rep["objective"],
    "report" -> rep,
    "iterations" -> maxIt,
    "history" -> If[returnHistory, history, Missing["NotRequested"]],
    "params" -> <|"P" -> P, "N" -> Nmax, "K" -> Kval, "useTail" -> useTail|>
  |>
];

End[];
EndPackage[];
