(* runHarness.wl - orchestrated test run and summary *)
repoRoot = DirectoryName[DirectoryName[$InputFileName]];
logsDir = FileNameJoin[{repoRoot, "logs"}];
If[Not@DirectoryQ[logsDir], CreateDirectory[logsDir]];

timestamp = DateString[{"Year","Month","Day","Hour","Minute","Second"}];
logPath = FileNameJoin[{logsDir, "runHarness_"<>timestamp<>".txt"}];
sumPath = FileNameJoin[{logsDir, "latest_summary.txt"}];
log = OpenWrite[logPath];
$ConstantsTestRunner = True;

logPrint[args__] := WriteString[log, ToString[Row[{args}], OutputForm]<>"\n"];

reports = {};
Block[{Print = logPrint, Exit = (Throw[{Exit, ##}] &), Quit = (Throw[{Quit, ##}] &)},
  Catch[
    AppendTo[reports, Get[FileNameJoin[{repoRoot,"tests","runTests.wl"}]]];
    AppendTo[reports, Get[FileNameJoin[{repoRoot,"tests","newtonTests.wl"}]]];
    AppendTo[reports, Get[FileNameJoin[{repoRoot,"tests","driverTests.wl"}]]];
    AppendTo[reports, Get[FileNameJoin[{repoRoot,"tests","tailTests.wl"}]]];
  ];
];
Close[log];

failCount = Total @ Map[
  If[Head[#] === TestReportObject,
    Quiet@Check[#"Statistics""FailureCount", Infinity],
    Infinity
  ]&,
  reports
];

(* Optimization check *)
ConstantsSetPrecision[80];
opt = Constants`Driver`NewtonOptimize[8, 8, MaxIterations->10, Tolerance->10^-8, Log->False];
aFinal = opt["aFinal"];
summary = "failCount="<>ToString[failCount]<>"\n"<>
  "length[aFinal]="<>ToString[Quiet@Check[Length[aFinal], -1]]<>"\n"<>
  "Total[aFinal]="<>ToString[Quiet@Check[Total[aFinal], Indeterminate], InputForm]<>"\n"<>
  "objectiveFinal="<>ToString[opt["objective"], InputForm]<>"\n"<>
  "log="<>logPath<>"\n";
Export[sumPath, summary, "String"];

Exit[If[failCount === 0, 0, 1]];
