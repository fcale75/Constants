(* ::Package:: *)

(* TailDerivatives.wl
   Mod-4 Kummer tail contributions for the gradient and Hessian of
     Sum_{k>N} Sk[a,k]^4
   using the same residue-class asymptotic series used by TailObjective.
*)

BeginPackage["Constants`TailDerivatives`"];

TailGradient::usage =
  "TailGradient[a,Nmax,K] returns the tail contribution to the gradient of \
C(a)=1/2+Sum_{k>=1} Sk[a,k]^4 for k>Nmax, using the mod-4 asymptotic series.";

TailHessian::usage =
  "TailHessian[a,Nmax,K] returns the tail contribution to the Hessian of \
C(a)=1/2+Sum_{k>=1} Sk[a,k]^4 for k>Nmax, using the mod-4 asymptotic series.";

Begin["`Private`"];

Needs["Constants`"];

(* Power series helper: given s_t for Sum_{t=0..K} s_t x^t, return coefficients of that^p. *)
ClearAll[powSeries];
powSeries[s_List, p_Integer?NonNegative, K_Integer?NonNegative] := Module[{x, poly},
  poly = Normal @ Series[(Sum[s[[t + 1]] * x^t, {t, 0, K}])^p, {x, 0, K}];
  Table[Coefficient[poly, x, u], {u, 0, K}]
];

(* Convolution of two coefficient lists truncated to order K. *)
ClearAll[convolve];
convolve[a_List, b_List, K_Integer?NonNegative] := Table[
  Sum[
    a[[i + 1]] * b[[t - i + 1]],
    {i, 0, t}
  ],
  {t, 0, K}
];

ClearAll[TailGradient];

TailGradient[a_List, Nmax_Integer?Positive, K_Integer?NonNegative] := Module[
  {prec = Constants`Private`$ConstantsWorkingPrecision, P = Length[a], grad,
   r, sk, sk3, n0, l, Jl, coeffs, u},

  grad = ConstantArray[0, P];

  Do[
    sk = Constants`SkSeriesCoefficients[a, r, K];
    sk3 = powSeries[sk, 3, K];
    n0 = Ceiling[(Nmax + 1 - r)/4];

    Do[
      Jl = Constants`SkSeriesCoefficients[UnitVector[P, l + 1], r, K];
      coeffs = 4 * convolve[Jl, sk3, K];

      Do[
        grad[[l + 1]] += N[
          coeffs[[u + 1]] * 4^(-(2 + u)) * HurwitzZeta[2 + u, n0 + r/4],
          prec
        ];
      , {u, 0, K}];
    , {l, 0, P - 1}];
  , {r, 0, 3}];

  N[grad, prec]
];

ClearAll[TailHessian];

TailHessian[a_List, Nmax_Integer?Positive, K_Integer?NonNegative] := Module[
  {prec = Constants`Private`$ConstantsWorkingPrecision, P = Length[a], hess,
   r, sk, sk2, n0, i, l, Ji, Jl, prod, u},

  hess = ConstantArray[0, {P, P}];

  Do[
    sk = Constants`SkSeriesCoefficients[a, r, K];
    sk2 = powSeries[sk, 2, K];
    n0 = Ceiling[(Nmax + 1 - r)/4];

    Do[
      Ji = Constants`SkSeriesCoefficients[UnitVector[P, i + 1], r, K];
      Do[
        Jl = Constants`SkSeriesCoefficients[UnitVector[P, l + 1], r, K];
        prod = convolve[convolve[Ji, Jl, K], sk2, K];

        Do[
          hess[[i + 1, l + 1]] += N[
            12 * prod[[u + 1]] * 4^(-(2 + u)) * HurwitzZeta[2 + u, n0 + r/4],
            prec
          ];
        , {u, 0, K}];
      , {l, 0, P - 1}];
    , {i, 0, P - 1}];
  , {r, 0, 3}];

  (* enforce symmetry explicitly *)
  hess = (hess + Transpose[hess]) / 2;
  N[hess, prec]
];

End[];
EndPackage[];
