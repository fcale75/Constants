(* Constants core package (P,N,K) harness *)

BeginPackage["Constants`"];

BesselBlock::usage = "BesselBlock[p,k] gives J_p(pi k/2) p! (4/(pi k))^p with memoization.";
Heartbeat::usage = "Heartbeat[msg] prints a heartbeat status line at most once every $ConstantsHeartbeatPeriod seconds.";
ConstantsSetPrecision::usage = "ConstantsSetPrecision[prec] sets working precision for computations.";

Begin["`Private`"];

If[!ValueQ[$ConstantsHeartbeatPeriod], $ConstantsHeartbeatPeriod = 30;];
If[!ValueQ[$ConstantsLastHeartbeat], $ConstantsLastHeartbeat = -Infinity;];
If[!ValueQ[$ConstantsWorkingPrecision], $ConstantsWorkingPrecision = 120;];

ConstantsSetPrecision[prec_Integer?Positive] := ($ConstantsWorkingPrecision = prec);

ConstantsLog[msg_] := Print[DateString[{"TimeShort"}], " | ", msg];

Heartbeat[msg_String:"still running"] := Module[{now = AbsoluteTime[]},
  If[now - $ConstantsLastHeartbeat >= $ConstantsHeartbeatPeriod,
    ConstantsLog[msg];
    $ConstantsLastHeartbeat = now;
  ];
];

ClearAll[BesselBlock];
BesselBlock[p_Integer?NonNegative, k_Integer?Positive] := Module[
  {prec = $ConstantsWorkingPrecision, z},
  z = N[Pi*k/2, prec];
  BesselBlock[p,k] = N[BesselJ[p, z] * N[Factorial[p], prec] * (N[4/(Pi*k), prec])^p, prec]
];

End[];
EndPackage[];
