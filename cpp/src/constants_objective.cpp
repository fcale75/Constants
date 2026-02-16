#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#include <flint/flint.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
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

struct Options {
  std::vector<std::string> coeffStrings;
  slong nmax = 256;
  slong kTerms = 32;
  slong precDps = 120;
  slong printDps = 60;
  bool useTail = true;
};

[[noreturn]] void UsageAndExit(const char* argv0) {
  std::cerr
      << "Usage:\n"
      << "  " << argv0
      << " --coeffs a0,a1,... [--nmax 256] [--k 32] [--prec-dps 120] [--print-dps 60] [--tail true|false]\n";
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
  throw std::runtime_error("Expected boolean value for --tail");
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
    throw std::runtime_error("Precision arguments must be positive");
  }

  return opt;
}

slong DigitsToBits(slong dps) {
  constexpr double kLog2_10 = 3.32192809488736234787;
  return static_cast<slong>(std::ceil(static_cast<double>(dps) * kLog2_10)) + 24;
}

std::vector<ArbNum> ParseCoefficients(const std::vector<std::string>& coeffs, slong precBits) {
  std::vector<ArbNum> out(coeffs.size());
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (arb_set_str(out[i].raw(), coeffs[i].c_str(), precBits) != 0) {
      throw std::runtime_error("Failed to parse coefficient: " + coeffs[i]);
    }
  }
  return out;
}

std::string ArbToString(const ArbNum& x, slong digits) {
  char* raw = arb_get_str(x.raw(), digits, ARB_STR_NO_RADIUS);
  std::string out(raw);
  flint_free(raw);
  return out;
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

ArbNum FiniteObjective(
    const std::vector<ArbNum>& a,
    const std::vector<std::vector<ArbNum>>& bessel,
    slong nmax,
    slong precBits) {
  ArbNum objective, sk, term, tmp;
  arb_set_ui(objective.raw(), 1);
  arb_div_ui(objective.raw(), objective.raw(), 2, precBits);

  const std::size_t pCount = a.size();
  for (slong k = 1; k <= nmax; ++k) {
    arb_zero(sk.raw());
    for (std::size_t p = 0; p < pCount; ++p) {
      arb_mul(tmp.raw(), a[p].raw(), bessel[p][static_cast<std::size_t>(k)].raw(), precBits);
      arb_add(sk.raw(), sk.raw(), tmp.raw(), precBits);
    }
    arb_pow_ui(term.raw(), sk.raw(), 4, precBits);
    arb_add(objective.raw(), objective.raw(), term.raw(), precBits);
  }

  return objective;
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

std::vector<ArbNum> BuildInvPiPowers(slong kTerms, slong precBits) {
  std::vector<ArbNum> out(static_cast<std::size_t>(kTerms + 1));
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

std::vector<ArbNum> SkSeriesCoefficients(
    const std::vector<ArbNum>& a,
    slong residue,
    slong kTerms,
    const std::vector<std::vector<ArbNum>>& asymCoeff,
    const std::vector<ArbNum>& prefactors,
    const std::vector<std::vector<ArbNum>>& cosTheta,
    const std::vector<std::vector<ArbNum>>& sinTheta,
    const std::vector<ArbNum>& invPiPowers,
    slong precBits) {
  std::vector<ArbNum> s(static_cast<std::size_t>(kTerms + 1));
  ArbNum term, tmp;

  for (std::size_t p = 0; p < a.size(); ++p) {
    if (static_cast<slong>(p) > kTerms) {
      continue;
    }

    for (slong m = 0; static_cast<slong>(p) + 2 * m <= kTerms; ++m) {
      const slong t = static_cast<slong>(p) + 2 * m;
      arb_mul(term.raw(), a[p].raw(), prefactors[p].raw(), precBits);
      arb_mul(term.raw(), term.raw(), cosTheta[static_cast<std::size_t>(residue)][p].raw(), precBits);
      arb_mul(term.raw(), term.raw(), asymCoeff[p][static_cast<std::size_t>(2 * m)].raw(), precBits);
      arb_mul(term.raw(), term.raw(), invPiPowers[static_cast<std::size_t>(2 * m)].raw(), precBits);
      if (m % 2 == 1) {
        arb_neg(term.raw(), term.raw());
      }
      arb_add(s[static_cast<std::size_t>(t)].raw(), s[static_cast<std::size_t>(t)].raw(), term.raw(), precBits);
    }

    for (slong m = 0; static_cast<slong>(p) + (2 * m + 1) <= kTerms; ++m) {
      const slong t = static_cast<slong>(p) + (2 * m + 1);
      arb_mul(term.raw(), a[p].raw(), prefactors[p].raw(), precBits);
      arb_mul(term.raw(), term.raw(), sinTheta[static_cast<std::size_t>(residue)][p].raw(), precBits);
      arb_mul(term.raw(), term.raw(), asymCoeff[p][static_cast<std::size_t>(2 * m + 1)].raw(), precBits);
      arb_mul(term.raw(), term.raw(), invPiPowers[static_cast<std::size_t>(2 * m + 1)].raw(), precBits);

      // odd terms use (-SinTheta) * (-1)^m
      if (m % 2 == 0) {
        arb_neg(term.raw(), term.raw());
      }
      arb_add(s[static_cast<std::size_t>(t)].raw(), s[static_cast<std::size_t>(t)].raw(), term.raw(), precBits);
    }
  }

  return s;
}

std::vector<ArbNum> Convolve(const std::vector<ArbNum>& a, const std::vector<ArbNum>& b, slong kTerms, slong precBits) {
  std::vector<ArbNum> out(static_cast<std::size_t>(kTerms + 1));
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

std::vector<ArbNum> ObjectiveTailCoefficients(const std::vector<ArbNum>& s, slong kTerms, slong precBits) {
  const auto sq = Convolve(s, s, kTerms, precBits);
  return Convolve(sq, sq, kTerms, precBits);
}

ArbNum TailObjective(
    const std::vector<ArbNum>& a,
    slong nmax,
    slong kTerms,
    slong precBits) {
  const std::size_t pCount = a.size();
  const auto asym = BuildAsymCoeff(pCount, kTerms, precBits);
  const auto pref = BuildPrefactors(pCount, precBits);
  std::vector<std::vector<ArbNum>> cosTheta;
  std::vector<std::vector<ArbNum>> sinTheta;
  BuildCosSinTheta(&cosTheta, &sinTheta, pCount, precBits);
  const auto invPiPowers = BuildInvPiPowers(kTerms, precBits);

  ArbNum sum, sVal, aVal, zetaVal, fourPow, fourNeg, term, tmp;
  arb_zero(sum.raw());

  const slong first = nmax + 1;
  for (slong r = 0; r < 4; ++r) {
    const auto s = SkSeriesCoefficients(a, r, kTerms, asym, pref, cosTheta, sinTheta, invPiPowers, precBits);
    const auto c = ObjectiveTailCoefficients(s, kTerms, precBits);

    const slong n0 = (first - r + 3) / 4;  // smallest n with (4n+r) >= first
    for (slong u = 0; u <= kTerms; ++u) {
      const slong sExp = 2 + u;

      arb_set_si(sVal.raw(), sExp);
      arb_set_si(aVal.raw(), 4 * n0 + r);
      arb_div_ui(aVal.raw(), aVal.raw(), 4, precBits);  // n0 + r/4
      arb_hurwitz_zeta(zetaVal.raw(), sVal.raw(), aVal.raw(), precBits);

      arb_ui_pow_ui(fourPow.raw(), 4, static_cast<ulong>(sExp), precBits);
      arb_inv(fourNeg.raw(), fourPow.raw(), precBits);

      arb_mul(term.raw(), c[static_cast<std::size_t>(u)].raw(), fourNeg.raw(), precBits);
      arb_mul(tmp.raw(), term.raw(), zetaVal.raw(), precBits);
      arb_add(sum.raw(), sum.raw(), tmp.raw(), precBits);
    }
  }

  return sum;
}

ArbNum SqrtArb(const ArbNum& x, slong precBits) {
  ArbNum out;
  arb_sqrt(out.raw(), x.raw(), precBits);
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options opt = ParseArgs(argc, argv);
    const slong precBits = DigitsToBits(opt.precDps);

    const auto a = ParseCoefficients(opt.coeffStrings, precBits);
    const std::size_t pCount = a.size();

    const auto t0 = std::chrono::steady_clock::now();
    const auto bessel = BuildBesselBlockTable(pCount, opt.nmax, precBits);
    ArbNum finite = FiniteObjective(a, bessel, opt.nmax, precBits);
    ArbNum tail;
    if (opt.useTail) {
      tail = TailObjective(a, opt.nmax, opt.kTerms, precBits);
    }

    ArbNum objective;
    arb_add(objective.raw(), finite.raw(), tail.raw(), precBits);
    ArbNum nu2 = SqrtArb(objective, precBits);

    const auto t1 = std::chrono::steady_clock::now();
    const double elapsed =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();

    std::cout << "P=" << pCount << "\n";
    std::cout << "Nmax=" << opt.nmax << "\n";
    std::cout << "K=" << opt.kTerms << "\n";
    std::cout << "precDps=" << opt.precDps << "\n";
    std::cout << "useTail=" << (opt.useTail ? "true" : "false") << "\n";
    std::cout << "finiteObjective=" << ArbToString(finite, opt.printDps) << "\n";
    std::cout << "tailObjective=" << ArbToString(tail, opt.printDps) << "\n";
    std::cout << "objective=" << ArbToString(objective, opt.printDps) << "\n";
    std::cout << "nu2Candidate=" << ArbToString(nu2, opt.printDps) << "\n";
    std::cout << "timeSec=" << elapsed << "\n";

    flint_cleanup();
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "error: " << ex.what() << "\n";
    return 1;
  }
}
