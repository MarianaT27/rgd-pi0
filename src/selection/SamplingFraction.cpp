#include "selection/SamplingFraction.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace pi0::selection {
namespace {

constexpr int kNumSectors = 6;

// ---------------------------------------------------------------------------
// CALIBRATION DATA -- TRANSCRIBED VERBATIM.
//
// Source: clas-analysis-1/clas12/include/clas12/services/SamplingFractionService.hpp,
//         SamplingFractionService::initialize_data(), the m_data["OB"] and
//         m_data["IB"] blocks. Six sectors each, eight coefficients each
//         (mu: a b c d, then sigma: a b c d), in sector order 1..6.
//         SINCE v1.0.0 the IB block comes from the common aspects analysis note
//         instead (see kInbending below); the OB block is unchanged.
//
// These are MEASUREMENTS, not code. Every digit is significant and none has
// been rounded, reformatted or "cleaned up": the E-notation and the digit
// count below match the source character-for-character so that a diff against
// it is meaningful. Do not reflow this table. Do not fold a common factor out
// of it. If it must change, it must change because the calorimeter was
// recalibrated, and then the whole table is replaced from the new fit output.
//
// Layout note: the source keyed these on the strings "OB"/"IB" in an
// unordered_map. Here they are two flat arrays indexed by (sector - 1), which
// is the same data with the string lookup removed.
// ---------------------------------------------------------------------------

// Target-averaged Outbending parameters (6 sectors). Equal, digit for digit, to the
// 'Average' rows of the common aspects analysis note (v3.1), Table 22 (OB).
constexpr std::array<SectorSf, kNumSectors> kOutbending = {{
    {{2.300896E-01, 4.862074E-03, -6.441645E-04, 1.247456E-05},
     {2.226749E-02, -3.186742E-03, 3.087387E-04, -1.247691E-05}},
    {{2.286452E-01, 4.985634E-03, -3.668113E-04, -1.110160E-05},
     {2.150934E-02, -2.549354E-03, 1.855833E-04, -4.394923E-06}},
    {{2.258053E-01, 9.073840E-03, -1.240745E-03, 3.579321E-05},
     {2.399538E-02, -4.208668E-03, 4.503779E-04, -1.783519E-05}},
    {{2.281016E-01, 6.073295E-03, -5.916479E-04, 2.569914E-06},
     {2.217527E-02, -2.788589E-03, 2.266746E-04, -7.717066E-06}},
    {{2.302059E-01, 3.020237E-03, -4.307581E-05, -2.216133E-05},
     {2.277357E-02, -3.941980E-03, 4.455085E-04, -1.912118E-05}},
    {{2.262583E-01, 6.849767E-03, -8.027577E-04, 1.288497E-05},
     {2.253050E-02, -3.334553E-03, 3.374127E-04, -1.413483E-05}},
}};

// Target-averaged Inbending parameters (6 sectors)
//
// SINCE v1.0.0: the 'Average' rows of the RG-D common aspects analysis note (v3.1),
// Table 23 (IB), digit for digit. Selected by Cuts::load(..., "inbending").
// Up to v0.2.8 this array held the legacy SamplingFractionService m_data["IB"]
// block, which was never used (the project was outbending-only) and does NOT
// match the note (S1 mu a: 2.311088E-01 there, 2.172761E-01 in the note).
// The outbending array above already equals the note's Table 22 'Average' rows.
constexpr std::array<SectorSf, kNumSectors> kInbending = {{
    {{2.172761E-01, 9.874863E-03, -1.460141E-03, 4.856991E-05},
     {2.970359E-02, -6.146500E-03, 8.519948E-04, -4.106624E-05}},
    {{2.295669E-01, 5.180274E-03, -9.473454E-06, -8.477999E-05},
     {2.903474E-02, -5.861235E-03, 7.084581E-04, -1.055444E-05}},
    {{2.287942E-01, 6.720455E-03, -6.166277E-04, -2.783523E-05},
     {2.597826E-02, -4.096651E-03, 3.504637E-04, 4.854673E-07}},
    {{2.279671E-01, 7.922469E-03, -8.822416E-04, 6.664547E-07},
     {2.932152E-02, -5.793094E-03, 7.024527E-04, -2.655757E-05}},
    {{2.293399E-01, 6.470559E-03, -7.222006E-04, 5.215601E-06},
     {2.609526E-02, -4.366779E-03, 4.176203E-04, -9.356880E-06}},
    {{2.307398E-01, 5.277840E-03, -3.128345E-04, -5.064927E-05},
     {2.886612E-02, -5.828894E-03, 6.541773E-04, -1.334127E-05}},
}};

[[nodiscard]] bool is_valid_sector(int sector) { return sector >= 1 && sector <= kNumSectors; }

}  // namespace

Polarity polarity_from_string(const std::string& polarity) {
    if (polarity == "outbending") return Polarity::Outbending;
    if (polarity == "inbending") return Polarity::Inbending;
    throw std::runtime_error("pi0::selection::polarity_from_string: unknown polarity '" + polarity +
                             "'; expected 'inbending' or 'outbending' (stageA_skim --polarity)");
}

const char* to_string(Polarity polarity) {
    return polarity == Polarity::Outbending ? "outbending" : "inbending";
}

const SectorSf& sf_params(int sector, Polarity polarity) {
    if (!is_valid_sector(sector)) {
        throw std::out_of_range("pi0::selection::sf_params: sector " + std::to_string(sector) +
                                " is outside [1, 6]");
    }
    const auto& table = (polarity == Polarity::Outbending) ? kOutbending : kInbending;
    return table[static_cast<std::size_t>(sector - 1)];
}

double sf_mu(double p_gev, int sector, Polarity polarity) {
    return sf_params(sector, polarity).mu.eval(p_gev);
}

double sf_sigma(double p_gev, int sector, Polarity polarity) {
    return sf_params(sector, polarity).sigma.eval(p_gev);
}

bool pass(double sf, double p_gev, int sector, Polarity polarity, double n_sigma) {
    // A track with no valid PCAL sector has no band to be tested against. The
    // old service returned false here too; keep that, and keep it BEFORE
    // sf_params() so this function never throws.
    if (!is_valid_sector(sector)) return false;

    const SectorSf& params = sf_params(sector, polarity);
    const double mu = params.mu.eval(p_gev);
    const double half_width = n_sigma * params.sigma.eval(p_gev);

    // Strict at both ends, matching every other comparison in this analysis.
    return (mu - half_width) < sf && sf < (mu + half_width);
}

}  // namespace pi0::selection
