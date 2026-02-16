(* runTests.wl
   Unit tests for Constants.wl (core definitions).
*)

root = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];

baseline[p_, k_, prec_] := Module[{z},
  z = SetPrecision[Pi*k/2, prec];
  N[BesselJ[p, z] * Factorial[p] * (4/(Pi*k))^p, prec]
];

tests = Module[{v60, v120, b120, k = 7, p = 2},
  ConstantsSetPrecision[60];
  v60 = BesselBlock[p, k];

  ConstantsSetPrecision[120];
  v120 = BesselBlock[p, k];
  b120 = baseline[p, k, 120];

  {
    VerificationTest[NumberQ[v60], True],
    VerificationTest[NumberQ[v120], True],

    (* Cache must respect changes in working precision (no precision poisoning). *)
    VerificationTest[Precision[v60] >= 50, True],
    VerificationTest[Precision[v120] >= 110, True],
    VerificationTest[Precision[v120] > Precision[v60], True],

    (* Value check against a direct baseline at high precision. *)
    VerificationTest[
      Abs[v120 - b120] < 10^-100,
      True
    ],

    (* Memoization check: repeated call at same precision yields identical value. *)
    VerificationTest[
      Module[{x1, x2},
        ConstantsSetPrecision[80];
        x1 = BesselBlock[1, 5];
        x2 = BesselBlock[1, 5];
        SameQ[x1, x2]
      ],
      True
    ]
  }
];

TestReport[tests]

