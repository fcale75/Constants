(* run_bootstrap_plan.wl
   Example persistent bootstrap run:
   - loads prior records from checkpoints/bootstrap_records.wl
   - runs staged continuation in P
   - appends fresh records
*)

repoRoot = "/Users/fcale/Dropbox/ChatGPT/Constants";
AppendTo[$Path, FileNameJoin[{repoRoot, "wl"}]];

Get[FileNameJoin[{repoRoot, "wl", "Constants.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Newton.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "TailDerivatives.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Driver.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Bootstrap.wl"}]];

recordsPath = FileNameJoin[{repoRoot, "checkpoints", "bootstrap_records.wl"}];

stages = {
  <|"P" -> 1, "Nmax" -> 256, "K" -> 32, "WorkingPrecision" -> 120, "MaxIterations" -> 25, "Log" -> False|>,
  <|"P" -> 2, "Nmax" -> 256, "K" -> 32, "WorkingPrecision" -> 120, "MaxIterations" -> 35, "Log" -> False|>,
  <|"P" -> 4, "Nmax" -> 256, "K" -> 32, "WorkingPrecision" -> 120, "MaxIterations" -> 35, "Log" -> False|>,
  <|"P" -> 8, "Nmax" -> 256, "K" -> 32, "WorkingPrecision" -> 120, "MaxIterations" -> 35, "Log" -> False|>
};

out = Constants`Bootstrap`RunBootstrapPlan[
  stages,
  RecordsPath -> recordsPath,
  SaveRecords -> True,
  Log -> True
];

If[FailureQ[out],
  Print["bootstrap run failed: ", out];
  Exit[1];
];

Print["bootstrap records file: ", recordsPath];
Do[
  rec = out[[i, "record"]];
  Print[
    "stage ", i,
    " P=", rec["params", "P"],
    " Nmax=", rec["params", "Nmax"],
    " wp=", rec["params", "WorkingPrecision"],
    " C=", ToString[rec["objective"], InputForm],
    " iter=", rec["iterations"],
    " init=", rec["initSource"]
  ];
, {i, 1, Length[out]}];

Exit[0];
