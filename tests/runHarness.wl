(* runHarness.wl
   Fast, deterministic harness: run unit tests, then do a cheap smoke check
   (no Newton iterations).
*)

root = DirectoryName[ExpandFileName[$InputFileName]];
repoRoot = DirectoryName[root];

AppendTo[$Path, FileNameJoin[{repoRoot, "wl"}]];

Get[FileNameJoin[{repoRoot, "wl", "Constants.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Newton.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "TailDerivatives.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Driver.wl"}]];

logsDir = FileNameJoin[{repoRoot, "logs"}];
If[!DirectoryQ[logsDir], CreateDirectory[logsDir]];

timestamp = DateString[{"Year", "Month", "Day", "Hour", "Minute", "Second"}];
logPath = FileNameJoin[{logsDir, "runHarness_" <> timestamp <> ".txt"}];
sumPath = FileNameJoin[{logsDir, "latest_summary.txt"}];
Export[sumPath, "STARTING RUN\n", "Text"];
log = OpenWrite[logPath];
logPrint[args__] := (WriteString[log, ToString[Row[{args}], OutputForm] <> "\n"]; Flush[log];);
stdoutLine[msg_] := (WriteString[$Output, msg <> "\n"]; Flush[$Output];);

reports = {};
$ConstantsTestRunner = True;

stdoutLine["[1/4] runTests.wl"];
Block[{Print = logPrint},
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "runTests.wl"}]], $Failed]];
];

stdoutLine["[2/4] newtonTests.wl"];
Block[{Print = logPrint},
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "newtonTests.wl"}]], $Failed]];
];

stdoutLine["[3/4] driverTests.wl"];
Export[sumPath, "REACHED STAGE 3/4 (driverTests)\n", "Text"];
Block[{Print = logPrint},
  AppendTo[reports, Quiet @ Check[Get[FileNameJoin[{repoRoot, "tests", "driverTests.wl"}]], $Failed]];
];

stdoutLine["[4/4] tailTests.wl"];
Block[{Print = logPrint},
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

(* If any test file failed to even return a TestReportObject, count it as a failure. *)
failCount = Total[Replace[failCounts, Infinity -> 1, {1}]];

(* ------------------------------------------------------------------------- *)
(* Smoke: single evaluation (NO Newton)                                       *)
(* ------------------------------------------------------------------------- *)

stdoutLine["[smoke] single-eval objective/derivatives"];

ConstantsClearCaches[];
ConstantsSetPrecision[60];

P = 8;
NMax = 32;
K = 8;

a0 = Normalize[Table[(-1)^n/(8^n*(n + 1)), {n, 0, P - 1}], Total];
a0 = a0 - (Total[a0] - 1)/P * ConstantArray[1, P];

obj = Quiet @ Check[Objective[a0, NMax, K], $Failed];
gFin = Quiet @ Check[Constants`Newton`NewtonDerivatives[a0, P, NMax, 0]["gpart"], $Failed];
gTail = Quiet @ Check[Constants`TailDerivatives`TailGradient[a0, NMax, K], $Failed];

(* Require fully numeric outputs (not symbolic leftovers). *)
smokeOK = And[
  NumericQ[obj],
  VectorQ[gFin, NumericQ] && Length[gFin] == P,
  VectorQ[gTail, NumericQ] && Length[gTail] == P
];

summary =
  "failCount=" <> ToString[failCount] <> "\n" <>
  "smokeOK=" <> ToString[smokeOK] <> "\n" <>
  "objective=" <> ToString[obj, InputForm] <> "\n" <>
  "Total[a0]=" <> ToString[Total[a0], InputForm] <> "\n" <>
  "log=" <> logPath <> "\n";

Export[sumPath, summary, "String"];

Print["Summary:"];
Print[summary];

Exit[If[failCount === 0 && smokeOK, 0, 1]];
