(* driverTests.wl -- integration test for Newton driver *)
root = DirectoryName[$InputFileName];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Newton.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Driver.wl"}]];

Constants`ConstantsSetPrecision[80];

result = Constants`Driver`NewtonOptimize[8, 20, "MaxIterations" -> 3, "Log" -> False];
aFinal = result[[1]];
lambdaFinal = result[[2]];
reportFinal = Constants`Newton`NewtonDerivatives[aFinal, 8, 20, lambdaFinal];

tests = {
  VerificationTest[VectorQ[aFinal], True],
  VerificationTest[Abs[Total[aFinal] - 1] < 10^-10, True],
  VerificationTest[NumberQ[reportFinal["objective"]], True]
};

Print[TestReport[tests]];
