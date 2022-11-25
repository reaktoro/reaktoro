// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class Phase;
class Species;
class SpeciesList;

/// The chemical properties of an aqueous phase.
class AqueousProps
{
public:
    /// Construct an uninitialized AqueousProps object with given chemical system.
    explicit AqueousProps(ChemicalSystem const& system);

    /// Construct an AqueousProps object with given chemical state of the system.
    explicit AqueousProps(ChemicalState const& state);

    /// Construct an AqueousProps object with given chemical properties of the system.
    explicit AqueousProps(ChemicalProps const& props);

    /// Construct a copy of a AqueousProps object.
    AqueousProps(AqueousProps const& other);

    /// Destroy this AqueousProps object.
    virtual ~AqueousProps();

    /// Assign a AqueousProps object to this object.
    auto operator=(AqueousProps other) -> AqueousProps&;

    /// Compute an AqueousProps object with given ChemicalProps object.
    /// Use this method if your model (e.g., reaction rate model) depends on an
    /// AqueousProps object and you are repeatedly constructing an AqueousProps
    /// object and updating it with the same ChemicalProps object, with same
    /// internal state. For example, when calculating mineral saturation indices
    /// on mineral reaction rates, each individual rate model requires building
    /// an AqueousProps object and updating it with the system's current
    /// ChemicalProps. By using this method, you take advantage of memoization
    /// optimization. If you supply a ChemicalProps object twice to this method
    /// with unchanged internal state, on the second call you will get a cached
    /// AqueousProps object.
    static auto compute(ChemicalProps const& props) -> AqueousProps const&;

    /// Set an activity model for a non-aqueous species that will be used in the calculation of its saturation index.
    /// @param species The name or index of the non-aqueous species in the list of species returned by @ref saturationSpecies.
    /// @param generator The activity model generator to be assigned for the non-aqueous species.
    auto setActivityModel(StringOrIndex const& species, ActivityModelGenerator const& generator) -> void;

    /// Update the aqueous properties with given chemical state of the system.
    auto update(ChemicalState const& state) -> void;

    /// Update the aqueous properties with given chemical properties of the system.
    auto update(ChemicalProps const& props) -> void;

    /// Return the temperature of the aqueous phase (in K).
    auto temperature() const -> real;

    /// Return the pressure of the aqueous phase (in Pa).
    auto pressure() const -> real;

    /// Return the amount of solvent water in the aqueous phase (in mol).
    auto waterAmount() const -> real;

    /// Return the mass of solvent water in the aqueous phase (in kg).
    auto waterMass() const -> real;

    /// Return the electric charge in the aqueous phase (in mol).
    auto charge() const -> real;

    /// Return the molality concentration in the aqueous phase (in molal).
    auto chargeMolality() const -> real;

    /// Return the molality of an element (in molal).
    auto elementMolality(StringOrIndex const& symbol) const -> real;

    /// Return the molality concentrations of the elements (in molal).
    auto elementMolalities() const -> ArrayXr;

    /// Return the molality of an aqueous solute species (in molal).
    auto speciesMolality(StringOrIndex const& name) const -> real;

    /// Return the molality concentrations of the species (in molal).
    auto speciesMolalities() const -> ArrayXr;

    /// Return the effective ionic strength of the aqueous phase (in molal). Equivalent to @ref ionicStrengthEffective.
    auto ionicStrength() const -> real;

    /// Return the effective ionic strength of the aqueous phase (in molal).
    auto ionicStrengthEffective() const -> real;

    /// Return the stoichiometric ionic strength of the aqueous phase (in molal).
    auto ionicStrengthStoichiometric() const -> real;

    /// Return the pH of the aqueous phase.
    auto pH() const -> real;

    /// Return the pE of the aqueous phase.
    auto pE() const -> real;

    /// Return the reduction potential of the aqueous phase (in V).
    auto Eh() const -> real;

    /// Return the total alkalinity of the aqueous phase (in eq/L).
    /// The total alkalinity (Alk) of the aqueous phase is by default
    /// calculated as the *acid neutralizing capacity* (ANC) of the solution
    /// using the formula:
    /// @f[
    /// \mathrm{Alk=[Na^{+}]+[K^{+}]+2[Ca^{2+}]+2[Mg^{2+}]-[Cl^{-}]-2[SO_{4}^{2-}]},
    /// @f] where @f$[\mathrm{species}]@f$ is the free molar concentration
    /// (mol/L) of the species in the solution. This formula is simpler,
    /// derived from the charge balance condition, and equivalent to the
    /// standard formula of total alkalinity.
    auto alkalinity() const -> real;

    /// Return the non-aqueous species that could be formed from the aqueous solution.
    /// This method returns, for example, gaseous and mineral species that
    /// could form when the aqueous solution is saturated with respect to them.
    auto saturationSpecies() const -> SpeciesList;

    /// Return the saturation index (@eq{\mathrm{SI} \equiv \log \Omega = \log \mathrm{IAP}/K}) of a non-aqueous species.
    /// @param species The name or index of the non-aqueous species in the list of species returned by @ref saturationSpecies.
    auto saturationIndex(StringOrIndex const& species) const -> real;

    /// Return the saturation indices of all non-aqueous species.
    /// These non-aqueous species can be obtained with @ref saturationSpecies.
    auto saturationIndices() const -> ArrayXr;

    /// Return the saturation ratio (@eq{\mathrm{SR} \equiv \Omega = \mathrm{IAP}/K}) of a non-aqueous species.
    /// @param species The name or index of the non-aqueous species in the list of species returned by @ref saturationSpecies.
    auto saturationRatio(StringOrIndex const& species) const -> real;

    /// Return the saturation ratios of all non-aqueous species.
    /// These non-aqueous species can be obtained with @ref saturationSpecies.
    auto saturationRatios() const -> ArrayXr;

    /// Return the saturation ratios of all non-aqueous species (in natural log).
    /// These non-aqueous species can be obtained with @ref saturationSpecies.
    auto saturationRatiosLn() const -> ArrayXr;

    /// Return the underlying ChemicalProps object.
    auto props() const -> ChemicalProps const&;

    /// Return the underlying ChemicalSystem object.
    auto system() const -> ChemicalSystem const&;

    /// Return the underlying Phase object corresponding to the aqueous phase.
    auto phase() const -> Phase const&;

    /// Output the properties of the aqueous phase to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the properties of the aqueous phase to a file.
    auto output(String const& filename) const -> void;

    // DEPRECATED METHODS : TO BE REMOVED IN THE NEAR FUTURE

    /// Return the saturation index of a given species (in natural log).
    /// @param species The name or index of the non-aqueous species in the list of species returned by @ref saturationSpecies.
    [[deprecated("Rely on the use of saturationIndex(species) instead.")]]
    auto saturationIndexLn(StringOrIndex const& species) const -> real;

    /// Return the saturation index of a given species (in log base 10).
    /// @param species The name or index of the non-aqueous species in the list of species returned by @ref saturationSpecies.
    [[deprecated("Rely on the use of saturationIndex(species) instead.")]]
    auto saturationIndexLg(StringOrIndex const& species) const -> real;

    /// Return the saturation indices of all non-aqueous species (in natural log).
    /// These non-aqueous species can be obtained with @ref saturationSpecies.
    [[deprecated("Rely on the use of saturationIndices() instead.")]]
    auto saturationIndicesLn() const -> ArrayXr;

    /// Return the saturation indices of all non-aqueous species (in log base 10).
    /// These non-aqueous species can be obtained with @ref saturationSpecies.
    [[deprecated("Rely on the use of saturationIndices() instead.")]]
    auto saturationIndicesLg() const -> ArrayXr;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output an AqueousProps object to an output stream.
auto operator<<(std::ostream& out, AqueousProps const& state) -> std::ostream&;

} // namespace Reaktoro
