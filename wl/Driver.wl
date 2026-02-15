(* Driver.wl - simple Newton iteration harness *)
BeginPackage["Constants`Driver`"]; NewtonOptimize;
Begin["`Private`"];
Needs["Constants`Heartbeat`"];
Needs["Constants`Newton`"];

Options[NewtonOptimize] = { MaxIterations -> 20, Tolerance -> 10^-25, Log -> True };

NewtonOptimize[pmax_Integer?Positive, nmax_Integer?Positive, opts : OptionsPattern[]] :=
 Module[{a = ConstantArray[1/pmax, pmax], lambda = 0, report, step,
   maxIt = OptionValue[MaxIterations],
   tol = OptionValue[Tolerance], it},
  Do[
   Constants`Heartbeat["Newton iter " <> ToString[it]];
   report = Constants`Newton`NewtonDerivatives[a, pmax, nmax, lambda];

   If[OptionValue[Log],
    Print["iter ", it, " C=", report["objective"], " ||g||=",
     Norm[report["gpart"]]]
    ];

   If[Norm[report["gpart"]] < tol,
    If[OptionValue[Log], Print["converged"]];
    Return[<|"aFinal" -> a, "lambdaFinal" -> lambda,
      "objective" -> report["objective"],
      "report" -> report|>]
    ];

   step = Constants`Newton`NewtonStep[a, lambda, pmax, nmax];
   a = a + step["deltaA"];
   lambda = lambda + step["deltaLambda"];
   , {it, 1, maxIt}
   ];

  <|"aFinal" -> a, "lambdaFinal" -> lambda,
    "objective" -> report["objective"],
    "report" -> report|>
  ];
End[];
EndPackage[];
