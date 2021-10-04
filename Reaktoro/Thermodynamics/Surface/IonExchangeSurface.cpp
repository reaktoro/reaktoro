// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "IonExchangeSurface.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct IonExchangeSurface::Impl
{
    /// All species on the ion exchange surface.
    SpeciesList species;

    /// The ion exchange species on the surface.
    SpeciesList exchange_species;

    /// The exchanger on the surface.
    SpeciesList exchanger_species;

    /// The array of exchanger's equivalents numbers for exchange species only (ze.size() = exchange_species.size()).
    ArrayXd ze;

    /// The index of the exchanger.
    Index idx_exchanger = 0;

    /// The indices of ion exchange species.
    Indices idx_exchange_species;

    /// Construct a default IonExchangeSurface::Impl instance.
    Impl()
    {}

    /// Construct an IonExchangeSurface::Impl instance with given species.
    Impl(const SpeciesList& species)
    : species(species)
    {
        /// Initialize the indices related data of the species
        initializeIndices();

        /// Initialize the array of exchanger's equivalents numbers.
        initializeExchangerEquivalentsNumbers();
    }

    // Return the number of exchanger's equivalents (the charge of cations) for the ion exchange species.
    auto exchangerEquivalentsNumber(const Species& species, const String& exchanger_symbol) -> real
    {
        // Run through the elements of the current species and return the coefficient of the exchanger
        for(auto [element, coeff] : species.elements())
            if(element.symbol() == exchanger_symbol)
                return coeff;

        // If all the elements are part of the periodic table then the exchanger is missing
        errorif(true, "Could not get information about the exchanger equivalents number. "
                      "Ensure the ion exchange phase contains correct species")
    }

    /// Initialize the array of exchanger's equivalents numbers (or cation charges) in all species.
    /// Note: exchanger's equivalents of the exchanger is assumed zero
    auto initializeExchangerEquivalentsNumbers() -> void
    {
        // The number of species in the ion exchange phase only
        const auto num_species = species.size();

        // The numbers of exchanger's equivalents for exchange species
        ze = ArrayXr::Zero(num_species);

        // Define the element symbol presenting the exchanger
        auto exchanger_symbol = species[idx_exchanger].elements().symbols()[0];

        // Initialize exchanger's equivalents by parsing the elements of the ion exchange species
        for(Index i : idx_exchange_species)
            ze[i] = exchangerEquivalentsNumber(species[i], exchanger_symbol);
    }

    /// Initialize the indices related data of the species.
    auto initializeIndices() -> void
    {
        // Initialize the array of indices with charged species
        Indices idx_charged_species;

        // Initialize the indices of the ion exchange species
        for(auto i = 0; i < species.size(); ++i)
            if(species[i].charge() == 0.0)
                idx_exchange_species.push_back(i);
            else
                idx_charged_species.push_back(i);

        // Initialize the index of the exchanger (assuming that it is the only charged species)
        idx_exchanger = idx_charged_species[0];
    }

    /// Initialize the lists of species separated into exchanger and exchange species on the surface
    auto initializeSpecies() -> void
    {
        // Separate exchanger (e.g., X-) from exchange species (NaX, KX, NaY, NaY, ect)
        exchanger_species = filter(species, [](const Species& s){ return s.charge() != 0.0;});
        exchange_species = filter(species, [](const Species& s){ return s.charge() == 0.0;});
    }

};

IonExchangeSurface::IonExchangeSurface()
: pimpl(new Impl())
{}

IonExchangeSurface::IonExchangeSurface(const SpeciesList& species)
: pimpl(new Impl(species))
{}

auto IonExchangeSurface::clone() const -> IonExchangeSurface
{
    IonExchangeSurface copy;
    *copy.pimpl = *pimpl;
    return copy;
}

auto IonExchangeSurface::species(Index idx) const -> const Species&
{
    return pimpl->species[idx];
}

auto IonExchangeSurface::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto IonExchangeSurface::ze() const -> ArrayXdConstRef
{
    return pimpl->ze;
}

auto IonExchangeSurface::indexExchanger() const -> Index
{
    return pimpl->idx_exchanger;
}

auto IonExchangeSurface::indicesExchange() const -> const Indices&
{
    return pimpl->idx_exchange_species;
}

} // namespace Reaktoro
