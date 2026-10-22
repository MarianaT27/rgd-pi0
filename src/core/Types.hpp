#pragma once

/// \file Types.hpp
/// \brief Plain-old-data types shared by the whole analysis.
///
/// These are deliberately dumb aggregates: no ROOT, no HIPO, no invariants
/// enforced in constructors. They exist so that the physics functions in
/// core/ can be exercised from a unit test without a data file.
///
/// Units: momenta and energies in GeV, angles in DEGREES.

#include <cmath>
#include <cstddef>
#include <limits>

namespace pi0 {

/// A photon candidate: momentum in GeV, energy taken massless (E = |p|).
struct Photon {
    double px{}, py{}, pz{};
    double e_gamma_deg{};  ///< angle to the scattered electron, DEGREES

    /// \return |p| in GeV.
    [[nodiscard]] double p() const { return std::sqrt(px * px + py * py + pz * pz); }

    /// \return the photon energy in GeV. Photons are massless here, so E = |p|.
    [[nodiscard]] double e() const { return p(); }
};

/// A reconstructed gamma-gamma pair. NOT necessarily a pi0 -- the mass window is
/// deliberately wide (see cuts.json pairing.mass_window_gev), because the signal
/// is extracted by fitting the m_gg spectrum downstream. Named honestly.
struct GGPair {
    std::size_t i{}, j{};  ///< indices into the photon list
    double mgg{};          ///< GeV
    double px{}, py{}, pz{}, e{};
};

/// DIS kinematics of one event.
struct DisKin {
    double q2{};  ///< GeV^2
    double nu{};  ///< GeV
    double xb{};  ///< dimensionless
    double w{};   ///< GeV
    double y{};   ///< dimensionless
};

/// SIDIS kinematics of one hadron w.r.t. the virtual photon.
struct SidisKin {
    double z{};           ///< E_h / nu, dimensionless
    double pt2{};         ///< GeV^2, transverse to the virtual photon
    double phi_h_deg{};   ///< Trento convention, DEGREES, range [-180, 180]
    /// Feynman x = p_L* / p_L,max* in the gamma*-nucleon CM frame (nucleon = proton
    /// mass, as for x_B and W). Positive: hadron forward along the virtual photon
    /// (current region); negative: backward (target region). p_L,max* is the
    /// exclusive limit e N -> e' pi0 N. NaN when undefined (W <= M_p + m_pi0, or no
    /// virtual photon), so a NaN can never land in a histogram bin by accident.
    double xf{std::numeric_limits<double>::quiet_NaN()};
};

}  // namespace pi0
