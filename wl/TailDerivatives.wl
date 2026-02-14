(* Tail derivatives for mod-4 accelerated series *)
BeginPackage["Constants`TailDerivatives`"];
TailGradient::usage = "TailGradient[a,N,K] gives tail contribution to \ngradient for objective with k>N via mod-4 accelerated series.";
TailHessian::usage = "TailHessian[a,N,K] gives tail contribution to \nHessian for objective with k>N via mod-4 accelerated series.";
Begin["`Private`"];
Needs["Constants`"];
ClearAll[powSeries, convolve];
powSeries[s_List, p_Integer?NonNegative, K_Integer?NonNegative] := Module[{x, poly},
  poly = Normal@Series[Sum[s[[t+1]] x^t, {t, 0, K}]^p, {x, 0, K}];
  Table[Coefficient[poly, x, u], {u, 0, K}]
];
convolve[a_List, b_List, K_Integer?NonNegative] :=
  Table[Sum[a[[i+1]] b[[t-i+1]], {i, 0, t}], {t, 0, K}];
TailGradient[a_List, N_Integer?Positive, K_Integer?NonNegative] := Module[
  {prec = Constants`Private`$ConstantsWorkingPrecision,
   P = Length[a], grad, r, coeffs, sk, sk3, n0},
  grad = ConstantArray[0, P];
  Do[
    sk = Constants`SkSeriesCoefficients[a, r, K];
    sk3 = powSeries[sk, 3, K];
    n0 = Ceiling[(N + 1 - r)/4];
    Do[
      coeffs = 4 convolve[Constants`SkSeriesCoefficients[UnitVector[P, l+1], r, K], sk3, K];
      Do[grad[[l+1]] +=
        N[coeffs[[u+1]] 4^(-(2+u)) HurwitzZeta[2+u, n0 + r/4], prec],
        {u, 0, K}],
      {l, 0, P-1}],
    {r, 0, 3}];
  N[grad, prec]
];
TailHessian[a_List, N_Integer?Positive, K_Integer?NonNegative] := Module[
  {prec = Constants`Private`$ConstantsWorkingPrecision,
   P = Length[a], hess, r, sk, sk2, n0, Ji, Jl, coeffs, prod},
  hess = ConstantArray[0, {P, P}];
  Do[
    sk = Constants`SkSeriesCoefficients[a, r, K];
    sk2 = powSeries[sk, 2, K];
    n0 = Ceiling[(N + 1 - r)/4];
    Do[
      Ji = Constants`SkSeriesCoefficients[UnitVector[P, i+1], r, K];
      Do[
        Jl = Constants`SkSeriesCoefficients[UnitVector[P, l+1], r, K];
        prod = convolve[convolve[Ji, Jl, K], sk2, K];
        coeffs = 12 prod;
        Do[hess[[i+1, l+1]] +=
          N[coeffs[[u+1]] 4^(-(2+u)) HurwitzZeta[2+u, n0 + r/4], prec],
          {u, 0, K}],
        {l, 0, P-1}],
      {i, 0, P-1}],
    {r, 0, 3}];
  (* enforce symmetry explicitly *)
  hess = (hess + Transpose[hess])/2;
  N[hess, prec]
];
End[];
EndPackage[];
