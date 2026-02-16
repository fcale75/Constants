(* ::Package:: *)

(* Bootstrap.wl
   Persistent bootstrap runner for coefficient/precision continuation.

   Records are stored as a Wolfram expression file (list of associations),
   preserving arbitrary-precision numbers exactly.
*)

BeginPackage["Constants`Bootstrap`"];

LoadBootstrapRecords::usage =
  "LoadBootstrapRecords[path] loads a list of saved bootstrap records (or {}).";

SaveBootstrapRecords::usage =
  "SaveBootstrapRecords[path, records] writes bootstrap records to path.";

AppendBootstrapRecord::usage =
  "AppendBootstrapRecord[path, record] appends one record to the persistent store.";

BestBootstrapRecord::usage =
  "BestBootstrapRecord[records, stage] returns the highest-precision compatible \
record for stage (same P/Nmax/K, precision <= stage precision), or Missing.";

BootstrapVector::usage =
  "BootstrapVector[aPrev, pTarget] zero-pads/truncates and renormalizes a coefficient vector.";

RunBootstrapStage::usage =
  "RunBootstrapStage[stage, opts] runs one NewtonOptimize stage with automatic \
initialization from previous result/records and returns <|\"result\", \"record\", \"initSource\"|>.";

RunBootstrapPlan::usage =
  "RunBootstrapPlan[stages, opts] runs stages sequentially, bootstrapping from \
previous stage and/or saved records, and returns a list of stage outputs.";

Options[RunBootstrapStage] = {
  PreviousResult -> Missing["None"],
  Records -> {},
  SavePath -> Automatic,
  Log -> True
};

Options[RunBootstrapPlan] = {
  RecordsPath -> Automatic,
  SaveRecords -> True,
  Log -> True
};

Begin["`Private`"];

Needs["Constants`Driver`"];

defaultRecordsPath[] := FileNameJoin[{Directory[], "checkpoints", "bootstrap_records.wl"}];

resolveNmax[stage_Association] := Module[{n1, n2},
  n1 = Lookup[stage, "Nmax", Missing["NotFound"]];
  n2 = Lookup[stage, "N", Missing["NotFound"]];
  If[n1 =!= Missing["NotFound"], n1, n2]
];

stageParamAssociation[stage_Association] := <|
  "P" -> Lookup[stage, "P", Missing["NotFound"]],
  "Nmax" -> resolveNmax[stage],
  "K" -> Lookup[stage, "K", Automatic],
  "WorkingPrecision" -> Lookup[stage, "WorkingPrecision", Automatic],
  "MaxIterations" -> Lookup[stage, "MaxIterations", 30],
  "Tolerance" -> Lookup[stage, "Tolerance", 10^-24],
  "ConstraintTolerance" -> Lookup[stage, "ConstraintTolerance", 10^-24],
  "UseTail" -> Lookup[stage, "UseTail", True]
|>;

LoadBootstrapRecords[path_String] := Module[{raw},
  If[!FileExistsQ[path], Return[{}]];
  raw = Quiet @ Check[Get[path], {}];
  If[ListQ[raw], raw, {}]
];

SaveBootstrapRecords[path_String, records_List] := Module[{dir},
  dir = DirectoryName[path];
  If[!DirectoryQ[dir], CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  Put[records, path];
  path
];

AppendBootstrapRecord[path_String, record_Association] := Module[{records},
  records = LoadBootstrapRecords[path];
  SaveBootstrapRecords[path, Append[records, record]]
];

compatibleRecordQ[rec_Association, stage_Association] := Module[
  {rp, sp, rk, sk, rn, sn, rwp, swp},
  rp = Lookup[Lookup[rec, "params", <||>], "P", Missing["NotFound"]];
  rn = Lookup[Lookup[rec, "params", <||>], "Nmax", Missing["NotFound"]];
  rk = Lookup[Lookup[rec, "params", <||>], "K", Missing["NotFound"]];
  rwp = Lookup[Lookup[rec, "params", <||>], "WorkingPrecision", Missing["NotFound"]];

  sp = Lookup[stage, "P", Missing["NotFound"]];
  sn = resolveNmax[stage];
  sk = Lookup[stage, "K", Missing["NotFound"]];
  swp = Lookup[stage, "WorkingPrecision", Missing["NotFound"]];

  rp === sp &&
  rn === sn &&
  rk === sk &&
  IntegerQ[rwp] && IntegerQ[swp] && rwp <= swp &&
  VectorQ[Lookup[rec, "aFinal", {}], NumericQ]
];

BestBootstrapRecord[records_List, stage_Association] := Module[{cands},
  cands = Select[records, compatibleRecordQ[#, stage] &];
  If[cands === {}, Return[Missing["NotFound"]]];
  First @ Reverse @ SortBy[cands, {
    Lookup[Lookup[#, "params", <||>], "WorkingPrecision", 0] &,
    Lookup[#, "timestamp", ""] &
  }]
];

BootstrapVector[aPrev_List, pTarget_Integer?Positive] := Module[
  {a},
  a = If[Length[aPrev] >= pTarget, Take[aPrev, pTarget], Join[aPrev, ConstantArray[0, pTarget - Length[aPrev]]]];
  a - (Total[a] - 1)/pTarget * ConstantArray[1, pTarget]
];

stageRunOptions[stage_Association] := {
  K -> Lookup[stage, "K", Automatic],
  WorkingPrecision -> Lookup[stage, "WorkingPrecision", Automatic],
  MaxIterations -> Lookup[stage, "MaxIterations", 30],
  Tolerance -> Lookup[stage, "Tolerance", 10^-24],
  ConstraintTolerance -> Lookup[stage, "ConstraintTolerance", 10^-24],
  Log -> Lookup[stage, "Log", False],
  UseTail -> Lookup[stage, "UseTail", True],
  ReturnHistory -> Lookup[stage, "ReturnHistory", False]
};

resolveStageInit[stage_Association, prevResult_, records_List] := Module[
  {initA = Lookup[stage, "InitialA", Automatic], initL = Lookup[stage, "InitialLambda", Automatic],
   p, best, source = "none"},
  p = Lookup[stage, "P", Missing["NotFound"]];

  If[initA =!= Automatic || initL =!= Automatic,
    Return[{initA, initL, "explicit"}]
  ];

  If[AssociationQ[prevResult] && VectorQ[Lookup[prevResult, "aFinal", {}], NumericQ] && IntegerQ[p],
    initA = BootstrapVector[prevResult["aFinal"], p];
    initL = Lookup[prevResult, "lambdaFinal", Automatic];
    Return[{initA, initL, "previous-stage"}]
  ];

  best = BestBootstrapRecord[records, stage];
  If[AssociationQ[best],
    initA = best["aFinal"];
    initL = Lookup[best, "lambdaFinal", Automatic];
    source = "saved-record";
    Return[{initA, initL, source}]
  ];

  {Automatic, Automatic, source}
];

buildRecord[stage_Association, result_Association, initSource_String] := Module[
  {obj, a, rep},
  obj = result["objective"];
  a = result["aFinal"];
  rep = result["report"];

  <|
    "timestamp" -> DateString[{"ISODate", " ", "Time"}],
    "initSource" -> initSource,
    "params" -> stageParamAssociation[stage],
    "iterations" -> result["iterations"],
    "objective" -> obj,
    "nu2Candidate" -> N[Sqrt[obj], 50],
    "kktInf" -> Norm[rep["res"], Infinity],
    "constraintAbs" -> Abs[1 - Total[a]],
    "aFinal" -> a,
    "lambdaFinal" -> result["lambdaFinal"]
  |>
];

RunBootstrapStage[stage_Association, opts : OptionsPattern[]] := Module[
  {prev = OptionValue[PreviousResult], records = OptionValue[Records], savePath = OptionValue[SavePath],
   logQ = TrueQ[OptionValue[Log]], p, nmax, initA, initL, initSource, runOpts, result, record},

  p = Lookup[stage, "P", Missing["NotFound"]];
  nmax = resolveNmax[stage];
  If[!IntegerQ[p] || !IntegerQ[nmax] || p <= 0 || nmax <= 0,
    Return[Failure["BadStage", <|"Message" -> "Stage must include positive integer P and Nmax/N."|>]]
  ];

  {initA, initL, initSource} = resolveStageInit[stage, prev, records];
  runOpts = stageRunOptions[stage];
  If[initA =!= Automatic, runOpts = Append[runOpts, InitialA -> initA]];
  If[initL =!= Automatic, runOpts = Append[runOpts, InitialLambda -> initL]];

  If[logQ,
    Print["Bootstrap stage P=", p, " Nmax=", nmax,
      " K=", Lookup[stage, "K", Automatic],
      " wp=", Lookup[stage, "WorkingPrecision", Automatic],
      " init=", initSource];
  ];

  result = Constants`Driver`NewtonOptimize[p, nmax, Sequence @@ runOpts];
  record = buildRecord[stage, result, initSource];

  If[StringQ[savePath], AppendBootstrapRecord[savePath, record]];

  <|"result" -> result, "record" -> record, "initSource" -> initSource|>
];

RunBootstrapPlan[stages_List, opts : OptionsPattern[]] := Module[
  {recordsPath = OptionValue[RecordsPath], saveQ = TrueQ[OptionValue[SaveRecords]],
   logQ = TrueQ[OptionValue[Log]], path, records, out = {}, prevResult = Missing["None"], stageOut, rec},

  path = If[recordsPath === Automatic, defaultRecordsPath[], recordsPath];
  records = If[StringQ[path], LoadBootstrapRecords[path], {}];

  Do[
    stageOut = RunBootstrapStage[
      stage,
      PreviousResult -> prevResult,
      Records -> records,
      SavePath -> If[saveQ && StringQ[path], path, Automatic],
      Log -> logQ
    ];

    If[FailureQ[stageOut], Return[stageOut]];

    AppendTo[out, stageOut];
    prevResult = stageOut["result"];
    rec = stageOut["record"];
    records = Append[records, rec];
  , {stage, stages}];

  out
];

End[];
EndPackage[];
