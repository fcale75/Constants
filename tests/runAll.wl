(* runAll.wl
   Entry point for wolframscript -file tests/runAll.wl
*)
root = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{root, "runHarness.wl"}]];

