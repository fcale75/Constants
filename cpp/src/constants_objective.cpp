#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#include <flint/flint.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cctype>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

class ArbNum {
public:
  ArbNum() {
    arb_init(v_);
    arb_zero(v_);
  }

  ArbNum(const ArbNum& other) {
    arb_init(v_);
    arb_set(v_, other.v_);
  }

  ArbNum& operator=(const ArbNum& other) {
    if (this != &other) {
      arb_set(v_, other.v_);
    }
    return *this;
  }

  ArbNum(ArbNum&& other) noexcept {
    arb_init(v_);
    arb_swap(v_, other.v_);
  }

  ArbNum& operator=(ArbNum&& other) noexcept {
    if (this != &other) {
      arb_swap(v_, other.v_);
    }
    return *this;
  }

  ~ArbNum() {
    arb_clear(v_);
  }

  arb_t& raw() { return v_; }
  const arb_t& raw() const { return v_; }

private:
  arb_t v_;
};

using ArbVec = std::vector<ArbNum>;
using ArbMat = std::vector<ArbVec>;

struct Options {
  std::vector<std::string> coeffStrings;
  slong nmax = 256;
  slong kTerms = 32;
  slong precDps = 120;
  slong printDps = 60;
  bool useTail = true;
  bool derivatives = false;
  bool computeHessian = true;
};

struct FiniteEval {
  ArbNum objective;
  ArbVec gradient;
  ArbMat hessian;
};

struct TailPrecompute {
  // basis[r][p][t]
  std::vector<std::vector<ArbVec>> basis;
  // weights[r][u] = 4^{-(2+u)} * HurwitzZeta(2+u, n0+r/4)
  std::vector<ArbVec> weights;
};

struct TailEval {
  ArbNum objective;
  ArbVec gradient;
  ArbMat hessian;
};

[[noreturn]] void UsageAndExit(const char* argv0) {
  std::cerr
      << "Usage:\n"
      << "  " << argv0
      << " --coeffs a0,a1,... [--nmax 256] [--k 32] [--prec-dps 120] [--print-dps 60]"
      << " [--tail true|false] [--derivatives true|false] [--hessian true|false]\n";
  std::exit(1);
}

std::vector<std::string> SplitCsv(const std::string& csv) {
  std::vector<std::string> out;
  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    std::string cleaned;
    for (char c : token) {
      if (!std::isspace(static_cast<unsigned char>(c))) {
        cleaned.push_back(c);
      }
    }
    if (!cleaned.empty()) {
      out.push_back(cleaned);
    }
  }
  return out;
}

bool ParseBool(const std::string& s) {
  if (s == "1" || s == "true" || s == "True" || s == "TRUE") {
    return true;
  }
  if (s == "0" || s == "false" || s == "False" || s == "FALSE") {
    return false;
  }
  throw std::runtime_error("Expected boolean argument");
}

Options ParseArgs(int argc, char** argv) {
  Options opt;

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--coeffs") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.coeffStrings = SplitCsv(argv[++i]);
      continue;
    }
    if (arg == "--nmax") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.nmax = std::stol(argv[++i]);
      continue;
    }
    if (arg == "--k") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.kTerms = std::stol(argv[++i]);
      continue;
    }
    if (arg == "--prec-dps") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.precDps = std::stol(argv[++i]);
      continue;
    }
    if (arg == "--print-dps") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.printDps = std::stol(argv[++i]);
      continue;
    }
    if (arg == "--tail") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.useTail = ParseBool(argv[++i]);
      continue;
    }
    if (arg == "--derivatives") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.derivatives = ParseBool(argv[++i]);
      continue;
    }
    if (arg == "--hessian") {
      if (i + 1 >= argc) UsageAndExit(argv[0]);
      opt.computeHessian = ParseBool(argv[++i]);
      continue;
    }
    if (arg == "--help" || arg == "-h") {
      UsageAndExit(argv[0]);
    }

    std::cerr << "Unknown argument: " << arg << "\n";
    UsageAndExit(argv[0]);
  }

  if (opt.coeffStrings.empty()) {
    throw std::runtime_error("Missing --coeffs");
  }
  if (opt.nmax <= 0) {
    throw std::runtime_error("--nmax must be positive");
  }
  if (opt.kTerms < 0) {
    throw std::runtime_error("--k must be non-negative");
  }
  if (opt.precDps <= 0 || opt.printDps <= 0) {
    throw std::runtime_error("Precision settings must be positive");
  }

  return opt;
}

slong DigitsToBits(slong dps) {
  constexpr double kLog2_10 = 3.32192809488736234787;
  return static_cast<slong>(std::ceil(static_cast<double>(dps) * kLog2_10)) + 24;
}

std::string ArbToString(const ArbNum& x, slong digits) {
  char* raw = arb_get_str(x.raw(), digits, ARB_STR_NO_RADIUS);
  std::string out(raw);
  flint_free(raw);
  return out;
}

std::string VecToCsv(const ArbVec& v, slong digits) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (i > 0) oss << ",";
    oss << ArbToString(v[i], digits);
  }
  return oss.str();
}

std::string MatToRowCsv(const ArbMat& m, slong digits) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < m.size(); ++i) {
    if (i > 0) oss << ";";
    for (std::size_t j = 0; j < m[i].size(); ++j) {
      if (j > 0) oss << ",";
      oss << ArbToString(m[i][j], digits);
    }
  }
  return oss.str();
}

ArbVec ParseCoefficients(const std::vector<std::string>& coeffs, slong precBits) {
  ArbVec out(coeffs.size());
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (arb_set_str(out[i].raw(), coeffs[i].c_str(), precBits) != 0) {
      throw std::runtime_error("Failed to parse coefficient: " + coeffs[i]);
    }
  }
  return out;
}

ArbMat ZeroMatrix(std::size_t n) {
  return ArbMat(n, ArbVec(n));
}

std::vector<std::vector<ArbNum>> BuildBesselBlockTable(std::size_t pCount, slong nmax, slong precBits) {
  std::vector<std::vector<ArbNum>> table(pCount, std::vector<ArbNum>(static_cast<std::size_t>(nmax + 1)));

  ArbNum pi, z, nu, besselJ, fac, ratio, ratioPow, denom, tmp;
  arb_const_pi(pi.raw(), precBits);

  for (std::size_t p = 0; p < pCount; ++p) {
    arb_set_ui(nu.raw(), static_cast<ulong>(p));
    arb_fac_ui(fac.raw(), static_cast<ulong>(p), precBits);

    for (slong k = 1; k <= nmax; ++k) {
      arb_mul_si(z.raw(), pi.raw(), k, precBits);
      arb_div_ui(z.raw(), z.raw(), 2, precBits);  // z = pi*k/2

      arb_hypgeom_bessel_j(besselJ.raw(), nu.raw(), z.raw(), precBits);

      arb_mul_si(denom.raw(), pi.raw(), k, precBits);   // pi*k
      arb_set_ui(ratio.raw(), 4);
      arb_div(ratio.raw(), ratio.raw(), denom.raw(), precBits);  // 4/(pi*k)
      arb_pow_ui(ratioPow.raw(), ratio.raw(), static_cast<ulong>(p), precBits);

      arb_mul(tmp.raw(), besselJ.raw(), fac.raw(), precBits);
      arb_mul(table[p][static_cast<std::size_t>(k)].raw(), tmp.raw(), ratioPow.raw(), precBits);
    }
  }

  return table;
}

ArbVec ComputeSkValues(const ArbVec& a, const std::vector<std::vector<ArbNum>>& bessel, slong nmax, slong precBits) {
  ArbVec skVals(static_cast<std::size_t>(nmax + 1));
  ArbNum tmp;

  for (slong k = 1; k <= nmax; ++k) {
    arb_zero(skVals[static_cast<std::size_t>(k)].raw());
    for (std::size_t p = 0; p < a.size(); ++p) {
      arb_mul(tmp.raw(), a[p].raw(), bessel[p][static_cast<std::size_t>(k)].raw(), precBits);
      arb_add(skVals[static_cast<std::size_t>(k)].raw(), skVals[static_cast<std::size_t>(k)].raw(), tmp.raw(), precBits);
    }
  }

  return skVals;
}

FiniteEval EvaluateFinite(
    const ArbVec& a,
    const std::vector<std::vector<ArbNum>>& bessel,
    slong nmax,
    bool wantDerivatives,
    bool wantHessian,
    slong precBits) {
  FiniteEval out;
  if (wantDerivatives) {
    out.gradient = ArbVec(a.size());
    if (wantHessian) {
      out.hessian = ZeroMatrix(a.size());
    }
  }

  ArbNum sk2, sk3, sk4, tmp, term;
  arb_set_ui(out.objective.raw(), 1);
  arb_div_ui(out.objective.raw(), out.objective.raw(), 2, precBits);

  const ArbVec skVals = ComputeSkValues(a, bessel, nmax, precBits);

  for (slong k = 1; k <= nmax; ++k) {
    const ArbNum& sk = skVals[static_cast<std::size_t>(k)];
    arb_pow_ui(sk2.raw(), sk.raw(), 2, precBits);
    arb_pow_ui(sk3.raw(), sk.raw(), 3, precBits);
    arb_pow_ui(sk4.raw(), sk.raw(), 4, precBits);
    arb_add(out.objective.raw(), out.objective.raw(), sk4.raw(), precBits);

    if (wantDerivatives) {
      for (std::size_t l = 0; l < a.size(); ++l) {
        arb_mul(tmp.raw(), bessel[l][static_cast<std::size_t>(k)].raw(), sk3.raw(), precBits);
        arb_mul_ui(term.raw(), tmp.raw(), 4, precBits);
        arb_add(out.gradient[l].raw(), out.gradient[l].raw(), term.raw(), precBits);
      }

      if (wantHessian) {
        for (std::size_t i = 0; i < a.size(); ++i) {
          for (std::size_t l = i; l < a.size(); ++l) {
            arb_mul(tmp.raw(), bessel[i][static_cast<std::size_t>(k)].raw(), bessel[l][static_cast<std::size_t>(k)].raw(), precBits);
            arb_mul(tmp.raw(), tmp.raw(), sk2.raw(), precBits);
            arb_mul_ui(term.raw(), tmp.raw(), 12, precBits);
            arb_add(out.hessian[i][l].raw(), out.hessian[i][l].raw(), term.raw(), precBits);
            if (l != i) {
              arb_add(out.hessian[l][i].raw(), out.hessian[l][i].raw(), term.raw(), precBits);
            }
          }
        }
      }
    }
  }

  return out;
}

std::vector<std::vector<ArbNum>> BuildAsymCoeff(std::size_t pCount, slong kTerms, slong precBits) {
  std::vector<std::vector<ArbNum>> asym(pCount, std::vector<ArbNum>(static_cast<std::size_t>(kTerms + 1)));
  ArbNum tmp;

  for (std::size_t p = 0; p < pCount; ++p) {
    arb_set_ui(asym[p][0].raw(), 1);
    for (slong n = 1; n <= kTerms; ++n) {
      const slong termInt =
          4 * static_cast<slong>(p) * static_cast<slong>(p) - (2 * n - 1) * (2 * n - 1);
      const slong denom = 8 * n;

      arb_mul_si(tmp.raw(), asym[p][static_cast<std::size_t>(n - 1)].raw(), termInt, precBits);
      arb_div_si(asym[p][static_cast<std::size_t>(n)].raw(), tmp.raw(), denom, precBits);
    }
  }

  return asym;
}

std::vector<ArbNum> BuildPrefactors(std::size_t pCount, slong precBits) {
  std::vector<ArbNum> pref(pCount);
  ArbNum fac, pow4, num, pi, piPow;
  arb_const_pi(pi.raw(), precBits);

  for (std::size_t p = 0; p < pCount; ++p) {
    arb_fac_ui(fac.raw(), static_cast<ulong>(p), precBits);
    arb_ui_pow_ui(pow4.raw(), 4, static_cast<ulong>(p), precBits);

    arb_mul(num.raw(), fac.raw(), pow4.raw(), precBits);
    arb_mul_ui(num.raw(), num.raw(), 2, precBits);

    arb_pow_ui(piPow.raw(), pi.raw(), static_cast<ulong>(p + 1), precBits);
    arb_div(pref[p].raw(), num.raw(), piPow.raw(), precBits);
  }

  return pref;
}

void BuildCosSinTheta(
    std::vector<std::vector<ArbNum>>* cosTheta,
    std::vector<std::vector<ArbNum>>* sinTheta,
    std::size_t pCount,
    slong precBits) {
  cosTheta->assign(4, std::vector<ArbNum>(pCount));
  sinTheta->assign(4, std::vector<ArbNum>(pCount));

  ArbNum x;
  for (slong r = 0; r < 4; ++r) {
    for (std::size_t p = 0; p < pCount; ++p) {
      const slong numer = 2 * r - 2 * static_cast<slong>(p) - 1;
      arb_set_si(x.raw(), numer);
      arb_div_ui(x.raw(), x.raw(), 4, precBits);  // x = (2r - 2p - 1) / 4
      arb_cos_pi((*cosTheta)[static_cast<std::size_t>(r)][p].raw(), x.raw(), precBits);
      arb_sin_pi((*sinTheta)[static_cast<std::size_t>(r)][p].raw(), x.raw(), precBits);
    }
  }
}

ArbVec BuildInvPiPowers(slong kTerms, slong precBits) {
  ArbVec out(static_cast<std::size_t>(kTerms + 1));
  ArbNum pi, invPi;
  arb_const_pi(pi.raw(), precBits);
  arb_inv(invPi.raw(), pi.raw(), precBits);

  arb_set_ui(out[0].raw(), 1);
  for (slong n = 1; n <= kTerms; ++n) {
    arb_mul(out[static_cast<std::size_t>(n)].raw(),
            out[static_cast<std::size_t>(n - 1)].raw(),
            invPi.raw(),
            precBits);
  }
  return out;
}

TailPrecompute BuildTailPrecompute(std::size_t pCount, slong nmax, slong kTerms, slong precBits) {
  TailPrecompute pre;
  pre.basis.assign(4, std::vector<ArbVec>(pCount, ArbVec(static_cast<std::size_t>(kTerms + 1))));
  pre.weights.assign(4, ArbVec(static_cast<std::size_t>(kTerms + 1)));

  const auto asym = BuildAsymCoeff(pCount, kTerms, precBits);
  const auto pref = BuildPrefactors(pCount, precBits);
  std::vector<std::vector<ArbNum>> cosTheta;
  std::vector<std::vector<ArbNum>> sinTheta;
  BuildCosSinTheta(&cosTheta, &sinTheta, pCount, precBits);
  const auto invPiPowers = BuildInvPiPowers(kTerms, precBits);

  ArbNum term, sVal, aVal, zetaVal, fourPow, fourNeg;
  const slong first = nmax + 1;

  for (slong r = 0; r < 4; ++r) {
    const slong n0 = (first - r + 3) / 4;

    for (std::size_t p = 0; p < pCount; ++p) {
      if (static_cast<slong>(p) > kTerms) {
        continue;
      }

      for (slong m = 0; static_cast<slong>(p) + 2 * m <= kTerms; ++m) {
        const slong t = static_cast<slong>(p) + 2 * m;
        arb_mul(term.raw(), pref[p].raw(), cosTheta[static_cast<std::size_t>(r)][p].raw(), precBits);
        arb_mul(term.raw(), term.raw(), asym[p][static_cast<std::size_t>(2 * m)].raw(), precBits);
        arb_mul(term.raw(), term.raw(), invPiPowers[static_cast<std::size_t>(2 * m)].raw(), precBits);
        if (m % 2 == 1) {
          arb_neg(term.raw(), term.raw());
        }
        arb_set(pre.basis[static_cast<std::size_t>(r)][p][static_cast<std::size_t>(t)].raw(), term.raw());
      }

      for (slong m = 0; static_cast<slong>(p) + (2 * m + 1) <= kTerms; ++m) {
        const slong t = static_cast<slong>(p) + (2 * m + 1);
        arb_mul(term.raw(), pref[p].raw(), sinTheta[static_cast<std::size_t>(r)][p].raw(), precBits);
        arb_mul(term.raw(), term.raw(), asym[p][static_cast<std::size_t>(2 * m + 1)].raw(), precBits);
        arb_mul(term.raw(), term.raw(), invPiPowers[static_cast<std::size_t>(2 * m + 1)].raw(), precBits);
        if (m % 2 == 0) {
          arb_neg(term.raw(), term.raw());
        }
        arb_set(pre.basis[static_cast<std::size_t>(r)][p][static_cast<std::size_t>(t)].raw(), term.raw());
      }
    }

    for (slong u = 0; u <= kTerms; ++u) {
      const slong sExp = 2 + u;
      arb_set_si(sVal.raw(), sExp);
      arb_set_si(aVal.raw(), 4 * n0 + r);
      arb_div_ui(aVal.raw(), aVal.raw(), 4, precBits);
      arb_hurwitz_zeta(zetaVal.raw(), sVal.raw(), aVal.raw(), precBits);

      arb_ui_pow_ui(fourPow.raw(), 4, static_cast<ulong>(sExp), precBits);
      arb_inv(fourNeg.raw(), fourPow.raw(), precBits);

      arb_mul(pre.weights[static_cast<std::size_t>(r)][static_cast<std::size_t>(u)].raw(),
              fourNeg.raw(),
              zetaVal.raw(),
              precBits);
    }
  }

  return pre;
}

ArbVec SkSeriesFromBasis(const ArbVec& a, const std::vector<ArbVec>& basisR, slong kTerms, slong precBits) {
  ArbVec s(static_cast<std::size_t>(kTerms + 1));
  ArbNum tmp;

  for (std::size_t p = 0; p < a.size(); ++p) {
    for (slong t = 0; t <= kTerms; ++t) {
      arb_mul(tmp.raw(), a[p].raw(), basisR[p][static_cast<std::size_t>(t)].raw(), precBits);
      arb_add(s[static_cast<std::size_t>(t)].raw(), s[static_cast<std::size_t>(t)].raw(), tmp.raw(), precBits);
    }
  }

  return s;
}

ArbVec Convolve(const ArbVec& a, const ArbVec& b, slong kTerms, slong precBits) {
  ArbVec out(static_cast<std::size_t>(kTerms + 1));
  ArbNum tmp;
  for (slong i = 0; i <= kTerms; ++i) {
    for (slong j = 0; j + i <= kTerms; ++j) {
      arb_mul(tmp.raw(), a[static_cast<std::size_t>(i)].raw(), b[static_cast<std::size_t>(j)].raw(), precBits);
      arb_add(out[static_cast<std::size_t>(i + j)].raw(),
              out[static_cast<std::size_t>(i + j)].raw(),
              tmp.raw(),
              precBits);
    }
  }
  return out;
}

ArbVec ObjectiveTailCoefficients(const ArbVec& s, slong kTerms, slong precBits) {
  const ArbVec sq = Convolve(s, s, kTerms, precBits);
  return Convolve(sq, sq, kTerms, precBits);
}

ArbNum WeightedDot(const ArbVec& coeff, const ArbVec& weights, slong kTerms, slong precBits) {
  ArbNum out, tmp;
  arb_zero(out.raw());
  for (slong u = 0; u <= kTerms; ++u) {
    arb_mul(tmp.raw(), coeff[static_cast<std::size_t>(u)].raw(), weights[static_cast<std::size_t>(u)].raw(), precBits);
    arb_add(out.raw(), out.raw(), tmp.raw(), precBits);
  }
  return out;
}

void AddScaledInPlace(ArbNum* target, const ArbNum& addend, slong scale, slong precBits) {
  ArbNum tmp;
  arb_mul_si(tmp.raw(), addend.raw(), scale, precBits);
  arb_add(target->raw(), target->raw(), tmp.raw(), precBits);
}

TailEval EvaluateTail(
    const ArbVec& a,
    const TailPrecompute& pre,
    slong kTerms,
    bool wantDerivatives,
    bool wantHessian,
    slong precBits) {
  TailEval out;
  const std::size_t pCount = a.size();

  if (wantDerivatives) {
    out.gradient = ArbVec(pCount);
    if (wantHessian) {
      out.hessian = ZeroMatrix(pCount);
    }
  }

  for (slong r = 0; r < 4; ++r) {
    const auto& basisR = pre.basis[static_cast<std::size_t>(r)];
    const auto& weightsR = pre.weights[static_cast<std::size_t>(r)];

    const ArbVec s = SkSeriesFromBasis(a, basisR, kTerms, precBits);
    const ArbVec c4 = ObjectiveTailCoefficients(s, kTerms, precBits);
    const ArbNum rObj = WeightedDot(c4, weightsR, kTerms, precBits);
    arb_add(out.objective.raw(), out.objective.raw(), rObj.raw(), precBits);

    if (wantDerivatives) {
      const ArbVec s2 = Convolve(s, s, kTerms, precBits);
      const ArbVec s3 = Convolve(s2, s, kTerms, precBits);

      for (std::size_t l = 0; l < pCount; ++l) {
        const ArbVec coeff = Convolve(basisR[l], s3, kTerms, precBits);
        const ArbNum d = WeightedDot(coeff, weightsR, kTerms, precBits);
        AddScaledInPlace(&out.gradient[l], d, 4, precBits);
      }

      if (wantHessian) {
        for (std::size_t i = 0; i < pCount; ++i) {
          for (std::size_t l = i; l < pCount; ++l) {
            const ArbVec basisProd = Convolve(basisR[i], basisR[l], kTerms, precBits);
            const ArbVec prod = Convolve(basisProd, s2, kTerms, precBits);
            const ArbNum h = WeightedDot(prod, weightsR, kTerms, precBits);

            ArbNum scaled;
            arb_mul_ui(scaled.raw(), h.raw(), 12, precBits);
            arb_add(out.hessian[i][l].raw(), out.hessian[i][l].raw(), scaled.raw(), precBits);
            if (i != l) {
              arb_add(out.hessian[l][i].raw(), out.hessian[l][i].raw(), scaled.raw(), precBits);
            }
          }
        }
      }
    }
  }

  return out;
}

ArbNum SqrtArb(const ArbNum& x, slong precBits) {
  ArbNum out;
  arb_sqrt(out.raw(), x.raw(), precBits);
  return out;
}

void AddVectorInPlace(ArbVec* base, const ArbVec& add, slong precBits) {
  for (std::size_t i = 0; i < base->size(); ++i) {
    arb_add((*base)[i].raw(), (*base)[i].raw(), add[i].raw(), precBits);
  }
}

void AddMatrixInPlace(ArbMat* base, const ArbMat& add, slong precBits) {
  for (std::size_t i = 0; i < base->size(); ++i) {
    for (std::size_t j = 0; j < (*base)[i].size(); ++j) {
      arb_add((*base)[i][j].raw(), (*base)[i][j].raw(), add[i][j].raw(), precBits);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options opt = ParseArgs(argc, argv);
    const slong precBits = DigitsToBits(opt.precDps);

    const ArbVec a = ParseCoefficients(opt.coeffStrings, precBits);
    const std::size_t pCount = a.size();

    const auto t0 = std::chrono::steady_clock::now();

    const auto tBesselStart = std::chrono::steady_clock::now();
    const auto bessel = BuildBesselBlockTable(pCount, opt.nmax, precBits);
    const auto tBesselEnd = std::chrono::steady_clock::now();

    const auto tFiniteStart = std::chrono::steady_clock::now();
    FiniteEval finite = EvaluateFinite(a, bessel, opt.nmax, opt.derivatives, opt.derivatives && opt.computeHessian, precBits);
    const auto tFiniteEnd = std::chrono::steady_clock::now();

    TailEval tail;
    std::chrono::steady_clock::time_point tTailPreStart, tTailPreEnd, tTailEvalStart, tTailEvalEnd;
    if (opt.useTail) {
      tTailPreStart = std::chrono::steady_clock::now();
      const TailPrecompute pre = BuildTailPrecompute(pCount, opt.nmax, opt.kTerms, precBits);
      tTailPreEnd = std::chrono::steady_clock::now();

      tTailEvalStart = std::chrono::steady_clock::now();
      tail = EvaluateTail(a, pre, opt.kTerms, opt.derivatives, opt.derivatives && opt.computeHessian, precBits);
      tTailEvalEnd = std::chrono::steady_clock::now();
    }

    ArbNum objective;
    arb_add(objective.raw(), finite.objective.raw(), tail.objective.raw(), precBits);
    ArbNum nu2 = SqrtArb(objective, precBits);

    ArbVec gradient;
    ArbMat hessian;
    if (opt.derivatives) {
      gradient = finite.gradient;
      if (opt.useTail) {
        AddVectorInPlace(&gradient, tail.gradient, precBits);
      }

      if (opt.computeHessian) {
        hessian = finite.hessian;
        if (opt.useTail) {
          AddMatrixInPlace(&hessian, tail.hessian, precBits);
        }
      }
    }

    const auto t1 = std::chrono::steady_clock::now();

    const double elapsed =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    const double tBessel =
        std::chrono::duration_cast<std::chrono::duration<double>>(tBesselEnd - tBesselStart).count();
    const double tFinite =
        std::chrono::duration_cast<std::chrono::duration<double>>(tFiniteEnd - tFiniteStart).count();

    double tTailPre = 0.0;
    double tTailEval = 0.0;
    if (opt.useTail) {
      tTailPre = std::chrono::duration_cast<std::chrono::duration<double>>(tTailPreEnd - tTailPreStart).count();
      tTailEval = std::chrono::duration_cast<std::chrono::duration<double>>(tTailEvalEnd - tTailEvalStart).count();
    }

    std::cout << "P=" << pCount << "\n";
    std::cout << "Nmax=" << opt.nmax << "\n";
    std::cout << "K=" << opt.kTerms << "\n";
    std::cout << "precDps=" << opt.precDps << "\n";
    std::cout << "useTail=" << (opt.useTail ? "true" : "false") << "\n";
    std::cout << "derivatives=" << (opt.derivatives ? "true" : "false") << "\n";
    std::cout << "hessian=" << (opt.derivatives && opt.computeHessian ? "true" : "false") << "\n";

    std::cout << "finiteObjective=" << ArbToString(finite.objective, opt.printDps) << "\n";
    std::cout << "tailObjective=" << ArbToString(tail.objective, opt.printDps) << "\n";
    std::cout << "objective=" << ArbToString(objective, opt.printDps) << "\n";
    std::cout << "nu2Candidate=" << ArbToString(nu2, opt.printDps) << "\n";

    if (opt.derivatives) {
      std::cout << "gradientCsv=" << VecToCsv(gradient, opt.printDps) << "\n";
      if (opt.computeHessian) {
        std::cout << "hessianRows=" << MatToRowCsv(hessian, opt.printDps) << "\n";
      }
    }

    std::cout << "timeBessel=" << tBessel << "\n";
    std::cout << "timeFinite=" << tFinite << "\n";
    std::cout << "timeTailPrecompute=" << tTailPre << "\n";
    std::cout << "timeTailEval=" << tTailEval << "\n";
    std::cout << "timeSec=" << elapsed << "\n";

    flint_cleanup();
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "error: " << ex.what() << "\n";
    return 1;
  }
}
