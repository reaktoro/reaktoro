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

#include "DoubleLayer.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Singletons/DissociationReactions.hpp>

namespace Reaktoro {

struct DoubleLayer::Impl
{
    /// All species on the diffusive double layer (DDL).
    SpeciesList species;

    /// The neutral DDL species in the mixture.
    SpeciesList neutral;

    /// The charged DDL species in the mixture.
    SpeciesList charged;

    /// The electric charges of the DDL species.
    ArrayXd z;

    /// The index of the water species.
    Index idx_water;

    /// The indices of the neutral DDL species.
    Indices idx_neutral_species;

    /// The indices of the charged DDL species.
    Indices idx_charged_species;

    /// The matrix that represents the dissociation of the DDL complexes into ions.
    MatrixXd dissociation_matrix;

    /// Construct a default DoubleLayer::Impl instance.
    Impl()
    {}

    /// Construct an DoubleLayer::Impl instance with given species.
    Impl(const SpeciesList& species)
    : species(species)
    {
        // Initialize the index related data
        initializeIndices();

        /// Initialize the array of electric charges of the species
        initializeCharges();

        // Initialize the dissociation matrix of the neutral species w.r.t. the charged species
        initializeDissociationMatrix();
    }

    /// Initialize the index related data of the species.
    auto initializeIndices() -> void
    {
        // Initialize the index of the water species
        idx_water = species.indexWithFormula("H2O!");
        std::cout << "index H2O! = " << idx_water << std::endl;

        // Initialize the indices of the charged and neutral species
        for(auto i = 0; i < species.size(); ++i)
        {
            if(i == idx_water)
                continue;
            if(species[i].charge() == 0.0)
            {
                idx_neutral_species.push_back(i);
                neutral.push_back(species[i]);
            }
            else
            {
                idx_charged_species.push_back(i);
                charged.push_back(species[i]);
            }
        }
    }

    /// Initialize the array of electric charges of the species
    auto initializeCharges() -> void
    {
        const auto charges = vectorize(species, RKT_LAMBDA(x, x.charge()));
        z = ArrayXd::Map(charges.data(), charges.size());
    }

    /// Initialize the dissociation matrix of the neutral species w.r.t. the charged species.
    auto initializeDissociationMatrix() -> void
    {
        // Return the stoichiometry of the i-th charged species in the j-th neutral species
        auto stoichiometry = [&](Index i, Index j)
        {
            const auto ineutral = idx_neutral_species[i];
            const auto icharged = idx_charged_species[j];
            const auto neutral = species[ineutral].formula();
            const auto charged = species[icharged].formula();
            return DissociationReactions::coefficient(neutral, charged);
        };

        // Assemble the dissociation matrix of the neutral species with respect to the charged species
        const auto num_charged_species = idx_charged_species.size();
        const auto num_neutral_species = idx_neutral_species.size();
        dissociation_matrix.resize(num_neutral_species, num_charged_species);
        for(auto i = 0; i < num_neutral_species; ++i)
            for(auto j = 0; j < num_charged_species; ++j)
                dissociation_matrix(i, j) = stoichiometry(i, j);
    }

    /// Return the molalities of the DDL species with given mole fractions.
    auto molalities(ArrayXrConstRef x) const -> ArrayXr
    {
        const auto xw = x[idx_water];
        if(xw == 0.0)
            return ArrayXr::Zero(x.size());
        return x/(waterMolarMass * xw);
    }

    /// Return the stoichiometric molalities of the charged species with given molalities.
    auto stoichiometricMolalities(ArrayXrConstRef m) const -> ArrayXr
    {
        // The molalities of the charged species
        const auto mc = m(idx_charged_species).matrix(); // convert from array to matrix expression

        // The molalities of the neutral species
        const auto mn = m(idx_neutral_species).matrix(); // convert from array to matrix expression

        // The stoichiometric molalities of the charged species
        const auto ms = mc + dissociation_matrix.transpose() * mn;

        return ms;
    }

    /// Return the effective ionic strength of the DDL mixture with given molalities.
    auto effectiveIonicStrength(ArrayXrConstRef m) const -> real
    {
        return 0.5 * (z * z * m).sum();
    }

    /// Return the stoichiometric ionic strength of the DDL mixture with given stoichiometric molalities of the charged species.
    auto stoichiometricIonicStrength(ArrayXrConstRef ms) const -> real
    {
        const auto zc = z(idx_charged_species);
        return 0.5 * (zc * zc * ms).sum();
    }

    /// Return the state of the diffusive double layer.
    auto state(real T, real P, ArrayXrConstRef x) -> DoubleLayerState
    {
        DoubleLayerState state;
        state.T = T;
        state.P = P;
        state.m  = molalities(x);
        state.ms = stoichiometricMolalities(state.m);
        state.Ie = effectiveIonicStrength(state.m);
        state.Is = stoichiometricIonicStrength(state.ms);
        return state;
    }

};

DoubleLayer::DoubleLayer()
: pimpl(new Impl())
{}

DoubleLayer::DoubleLayer(const SpeciesList& species)
: pimpl(new Impl(species))
{}

auto DoubleLayer::clone() const -> DoubleLayer
{
    DoubleLayer copy;
    *copy.pimpl = *pimpl;
    return copy;
}

auto DoubleLayer::species(Index idx) const -> const Species&
{
    return pimpl->species[idx];
}

auto DoubleLayer::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto DoubleLayer::indexWater() const -> Index
{
    return pimpl->idx_water;
}

auto DoubleLayer::charges() const -> ArrayXdConstRef
{
    return pimpl->z;
}

auto DoubleLayer::state(real T, real P, ArrayXrConstRef x) -> DoubleLayerState
{
    return pimpl->state(T, P, x);
}

} // namespace Reaktoro
