(* newtonTests.wl - additional tests for Newton module *)
root = DirectoryName[$InputFileName];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Newton.wl"}]];

Constants`ConstantsSetPrecision[80];

Clear[objective];
objective[a_, p_, n_] := NewtonObjective[a,p,n];

p = 3; n = 8; a = ConstantArray[1/p, p]; lambda = 0;
report = NewtonDerivatives[a,p,n,lambda];

(* Hessian symmetry *)
hSym = VerificationTest[
  report["hess"] == Transpose[report["hess"]],
  True
];

(* Numerical gradient check for gpart (gradient of objective only) *)
eps = 10^-25;
numGrad = Table[
  (objective[a + eps UnitVector[p, l], p, n] - objective[a, p, n])/eps,
  {l, 1, p}
];
gCheck = VerificationTest[
  Max[Abs[numGrad - report["gpart"]]] < 10^-15,
  True
];

Print[TestReport[{hSym, gCheck}]];
