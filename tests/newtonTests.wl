(* newtonTests.wl
   Unit tests for Newton.wl (finite-sum objective/derivatives).
*)

root = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Newton.wl"}]];

ConstantsSetPrecision[80];

P = 3;
Nmax = 12;

a = ConstantArray[1/P, P];
lambda = 0;

rep = Constants`Newton`NewtonDerivatives[a, P, Nmax, lambda];

objective[a_] := Constants`Newton`NewtonObjective[a, P, Nmax];

eps = SetPrecision[10^-18, 50];
numGrad = Table[
  (objective[a + eps UnitVector[P, l]] - objective[a - eps UnitVector[P, l]])/(2 eps),
  {l, 1, P}
];

tests = {
  (* Hessian should be symmetric up to numerical noise. *)
  VerificationTest[Norm[rep["hess"] - Transpose[rep["hess"]], Infinity] < 10^-40, True],

  (* Gradient check (finite differences). *)
  VerificationTest[Norm[numGrad - rep["gpart"], Infinity] < 10^-10, True],

  (* KKT step keeps Total[a] close to 1 if we start on the constraint hyperplane. *)
  VerificationTest[
    Module[{step, aNext},
      step = Constants`Newton`NewtonStep[a, lambda, P, Nmax];
      aNext = step["aNext"];
      Abs[Total[aNext] - 1] < 10^-30
    ],
    True
  ]
};

TestReport[tests]

