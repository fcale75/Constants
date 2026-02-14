(* runTests.wl -- MUnit harness for Constants *)

root = DirectoryName[$InputFileName];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];

ConstantsSetPrecision[120];

baseline[p_, k_, prec_: 120] := Module[{z},
  z = N[Pi*k/2, prec];
  N[BesselJ[p, z] * N[Factorial[p], prec] * (N[4/(Pi*k), prec])^p, prec]
];

tests = {
  VerificationTest[
    BesselBlock[0, 1],
    baseline[0, 1, 120],
    SameTest -> (Abs[#1 - #2] < 10^-100 &)
  ],
  VerificationTest[
    BesselBlock[2, 7],
    baseline[2, 7, 120],
    SameTest -> (Abs[#1 - #2] < 10^-100 &)
  ],
  VerificationTest[
    NumberQ[BesselBlock[3, 4]],
    True
  ],
  VerificationTest[
    Block[{},
      v1 = BesselBlock[1, 5];
      v2 = BesselBlock[1, 5];
      SameQ[v1, v2]
    ],
    True
  ]
};

Print["Running tests...\n"];;
report = TestReport[tests];
Print[report];
If[FailureQ[report], Exit[1], Exit[0]];
