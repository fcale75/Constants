(* runHarness.wl
   Orchestrated test run + a small optimization smoke test.

   Designed to be executed via:
     wolframscript -file tests/runAll.wl
*)

root = DirectoryName[ExpandFileName[$InputFileName]];
repoRoot = DirectoryName[root];

(* Ensure package path includes wl/ *)
AppendTo[$Path, FileNameJoin[{repoRoot, "wl"}]];

(* Load core packages explicitly (tests also Get these, but this avoids path issues). *)
Get[FileNameJoin[{repoRoot, "wl", "Constants.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Newton.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "TailDerivatives.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Driver.wl"}]];

logsDir = FileNameJoin[{repoRoot, "logs"}];
If[!DirectoryQ[logsDir], CreateDirectory[logsDir]];

timestamp = DateString[{"Year", "Month", "Day", "Hour", "Minute", "Second"}];
logPath = FileNameJoin[{logsDir, "runHarness_" <> timestamp <> ".txt"}];
sumPath = FileNameJoin[{logsDir, "latest_summary.txt"}];

log = OpenWrite[logPath];
logPrint[args__] := WriteString[log, ToString[Row[{args}], OutputForm] <> "\n"];

reports = {};

$ConstantsTestRunner = True;

Block[{Print = logPrint},
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "runTests.wl"}]], $Failed]];
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "newtonTests.wl"}]], $Failed]];
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "driverTests.wl"}]], $Failed]];
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "tailTests.wl"}]], $Failed]];
];

Close[log];

failCounts = Table[
  If[Head[r] === TestReportObject,
    Quiet @ Check[r["TestsFailedCount"], Infinity],
    Infinity
  ],
  {r, reports}
];

failCount = Total[failCounts];

(* ------------------------------------------------------------------------- *)
(* Optimization smoke test                                                     *)
(* ------------------------------------------------------------------------- *)

ConstantsSetPrecision[80];

opt = Constants`Driver`NewtonOptimize[
  8, 32,
  K -> 8,
  WorkingPrecision -> 80,
  MaxIterations -> 10,
  Tolerance -> 10^-8,
  ConstraintTolerance -> 10^-12,
  Log -> False,
  ReturnHistory -> False
];

aFinal = Quiet @ Check[opt["aFinal"], {}];
lambdaFinal = Quiet @ Check[opt["lambdaFinal"], Indeterminate];
objFinal = Quiet @ Check[opt["objective"], Indeterminate];
resInf = Quiet @ Check[Norm[opt["report"]["res"], Infinity], Indeterminate];
conAbs = Quiet @ Check[Abs[opt["report"]["constraint"]], Indeterminate];

summary =
  "failCount=" <> ToString[failCount] <> "\n" <>
  "length[aFinal]=" <> ToString[Quiet @ Check[Length[aFinal], -1]] <> "\n" <>
  "Total[aFinal]=" <> ToString[Quiet @ Check[Total[aFinal], Indeterminate], InputForm] <> "\n" <>
  "lambdaFinal=" <> ToString[lambdaFinal, InputForm] <> "\n" <>
  "objectiveFinal=" <> ToString[objFinal, InputForm] <> "\n" <>
  "stationarityInfNorm=" <> ToString[resInf, InputForm] <> "\n" <>
  "constraintAbs=" <> ToString[conAbs, InputForm] <> "\n" <>
  "log=" <> logPath <> "\n";

Export[sumPath, summary, "String"];

Print["Summary:"];
Print[summary];

Exit[If[failCount === 0, 0, 1]];

