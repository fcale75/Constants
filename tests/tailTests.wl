(* tailTests.wl
   Tests for mod-4 tail objective/derivatives consistency.

   These tests are designed to catch:
     - wrong asymptotic scaling in SkSeriesCoefficients (e.g., bad pi factors)
     - missing tail terms in gradient/Hessian
*)

root = DirectoryName[ExpandFileName[$InputFileName]];
Get[FileNameJoin[{root, "..", "wl", "Constants.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "Newton.wl"}]];
Get[FileNameJoin[{root, "..", "wl", "TailDerivatives.wl"}]];

ConstantsSetPrecision[80];

P = 4;
Nmax = 40;
K = 10;

a = ConstantArray[1/P, P];

obj[x_] := Objective[x, Nmax, K];

eps = SetPrecision[10^-10, 50];
numGrad = Table[
  (obj[a + eps UnitVector[P, j]] - obj[a - eps UnitVector[P, j]])/(2 eps),
  {j, 1, P}
];

rep = Constants`Newton`NewtonDerivatives[a, P, Nmax, 0];
g = rep["gpart"] + Constants`TailDerivatives`TailGradient[a, Nmax, K];
h = rep["hess"] + Constants`TailDerivatives`TailHessian[a, Nmax, K];

(* Series approximation sanity check at a reasonably large k. *)
kLarge = 201;
r = Mod[kLarge, 4];
s = SkSeriesCoefficients[a, r, K];
approx = kLarge^(-1/2) * Sum[s[[t + 1]] * kLarge^(-t), {t, 0, K}];
direct = Sk[a, kLarge];

tests = {
  VerificationTest[VectorQ[g, NumericQ] && Length[g] == P, True],
  VerificationTest[MatrixQ[h, NumericQ] && Dimensions[h] == {P, P}, True],

  (* Gradient must match finite differences to modest tolerance. *)
  VerificationTest[Norm[g - numGrad, Infinity] < 10^-6, True],

  (* Hessian symmetry. *)
  VerificationTest[Norm[h - Transpose[h], Infinity] < 10^-30, True],

  (* Asymptotic series should approximate Sk at large k to a modest relative error. *)
  VerificationTest[Abs[direct - approx]/Max[Abs[direct], 10^-30] < 5*10^-4, True]
};

TestReport[tests]
