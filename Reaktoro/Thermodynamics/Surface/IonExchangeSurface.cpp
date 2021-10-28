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

    /// The element symbol representing the exchanger
    String exchanger_symbol;

    /// The array of exchanger's equivalents numbers for exchange species only (ze.size() = exchange_species.size()).
    ArrayXd ze;

    /// Construct a default IonExchangeSurface::Impl instance.
    Impl()
    {}

    /// Construct an IonExchangeSurface::Impl instance with given species.
    Impl(const SpeciesList& species)
    : species(species)
    {
        /// Initialize the symbol of the exchanger
        intializeExchanger();

        /// Initialize the array of exchanger's equivalents numbers.
        initializeExchangerEquivalentsNumbers();
    }

    /// Initialize the symbol representing the exchanger.
    /// The method parses the ion exchange species list and identifies a common element that will be regarded as exchanger
    /// For example, for the list of species NaX, CaX2, MgX2, KX, the `exchanger_symbol` is `X`
    auto intializeExchanger() -> void
    {
        errorif(species.size() == 0, "There is no species in the IonExchangePhase")

        // Fetch elements' symbols from the first species
        auto esymbols = species[0].elements().symbols();
        auto num_esymbols = esymbols.size();

        // Create the auxiliary vector to count the matches of the symbol
        std::vector<Index> count_elements(esymbols.size());

        // Loop over symbols
        for(auto i = 0; i < num_esymbols; i++)
        {
            // Loop over all species to check for the existence of the current symbol in them
            for(auto s : species)
            {
                if(contains(s.elements().symbols(), esymbols[i]))
                    count_elements[i]++;
            }
            // Step out if the common symbol has been found
            if(count_elements[i] == species.size())
            {
                exchanger_symbol = esymbols[i];
                break;
            }
        }
        // Raise an error if the exchanger_symbol wasn't found
        errorif(exchanger_symbol.empty(), "There is no common element amount species to represent an exchanger.")
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

        // Initialize exchanger's equivalents by parsing the species on the ion exchange surface
        for(int i = 0; i < num_species; ++i)
            ze[i] = exchangerEquivalentsNumber(species[i], exchanger_symbol);
    }

    /// Return the equivalences fractions of the species on ion exchange surface if molar fractions are provided.
    auto equivalencesFractions(ArrayXrConstRef x) const -> ArrayXr
    {
        // beta_i = xi * zi / sum_c (xc * zc)
        return x*ze/(x*ze).sum();
    }

    /// Return the state of the ion exchange surface.
    auto state(real T, real P, ArrayXrConstRef x) -> IonExchangeSurfaceState
    {
        IonExchangeSurfaceState exchange_state;
        exchange_state.beta = equivalencesFractions(x);
        return exchange_state;
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

auto IonExchangeSurface::state(real T, real P, ArrayXrConstRef x) -> IonExchangeSurfaceState
{
    return pimpl->state(T, P, x);
}

} // namespace Reaktoro
