If[Not@DirectoryQ["logs"], CreateDirectory["logs"]];

timestamp = DateString[{"Year","Month","Day","Hour","Minute","Second"}];
logPath = FileNameJoin[{"logs", "runAndLog_"<>timestamp<>".txt"}];
log = OpenWrite[logPath];
print0 = Print;
logPrint[args__] := (WriteString[log, ToString[Row[{args}], OutputForm]<>"\n"); print0[args];

$ConstantsTestRunner = True;

Block[{Print = logPrint, Exit = Function[{status}, Null]},
  print0["Running tests/runAndLog.wl; log=", logPath];
  report =.;

  Get["tests/runTests.wl"]; runTestsReport = report;
  Get["tests/newtonTests.wl"]; newtonTestsReport = report;
  Get["tests/driverTests.wl"]; driverTestsReport = report;
  Get["tests/tailTests.wl"]; tailTestsReport = report;

  WriteString[log, "\nReports: \n"]; 
  WriteString[log, ToString[{runTestsReport, newtonTestsReport, driverTestsReport, tailTestsReport}, InputForm]];
  Close[log];
];
