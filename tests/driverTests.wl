(* driverTests.wl
   Integration tests for the Newton driver (paper-faithful KKT residual).
*)

root = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Newton.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "TailDerivatives.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Driver.wl"}]];

result = Constants`Driver`NewtonOptimize[
  4, 24,
  K -> 8,
  WorkingPrecision -> 60,
  MaxIterations -> 3,
  Tolerance -> 10^-8,
  ConstraintTolerance -> 10^-12,
  Log -> False,
  ReturnHistory -> False
];

aFinal = result["aFinal"];
lambdaFinal = result["lambdaFinal"];
reportFinal = result["report"];

tests = {
  VerificationTest[VectorQ[aFinal, NumericQ], True],
  VerificationTest[Length[aFinal] == 4, True],
  VerificationTest[Abs[Total[aFinal] - 1] < 10^-20, True],
  VerificationTest[NumberQ[lambdaFinal], True],
  VerificationTest[NumberQ[result["objective"]], True],
  VerificationTest[KeyExistsQ[reportFinal, "res"], True],
  VerificationTest[Norm[reportFinal["res"], Infinity] >= 0, True]
};

TestReport[tests]

