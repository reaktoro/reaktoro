// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "ChemicalQuantity.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

struct ChemicalQuantity::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The reactions in the chemical system
    ReactionSystem reactions;

    /// The chemical state of the system
    ChemicalState state;

    /// The thermodynamic properties of the chemical system at (*T*, *P*, **n**)
    ChemicalProperties properties;

    /// The progress variable at which the chemical state is referred (in units of s)
    double t;

    /// The temperature of the chemical system (in units of K).
    double T;

    /// The pressure of the chemical system (in units of Pa).
    double P;

    /// The molar amounts of the species in the chemical system (in units of mol).
    Vector n;

    /// The rates of the reactions in the chemical system (in units of mol/s).
    ChemicalVector r;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance with given ChemicalSystem object
    Impl(const ChemicalSystem& system)
    : system(system)
    {}

    /// Construct a custom Impl instance with given ReactionSystem object
    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions)
    {}

    /// Update the state of the chemical quantity instance
    auto update(const ChemicalState& state) -> void
    {
        update(state, 0.0);
    }

    /// Update the state of the chemical quantity instance
    auto update(const ChemicalState& state_, double t_) -> void
    {
        // Update the chemical state of the system
        state = state_;

        // Update the progress variable
        t = t_;

        // Update the temperature, pressure and molar composition of the system
        T = state.temperature();
        P = state.pressure();
        n = state.speciesAmounts();

        // Update the thermodynamic properties of the system
        properties = system.properties(T, P, n);

        // Update the rates of the reactions
        if(!reactions.reactions().empty())
            r = reactions.rates(properties);
    }

    auto value(std::string str) const -> double
    {
        auto words = split(str, ":");

        std::string quantity = words[0];
        std::string units = words.size() > 1 ? words[1] : "";

        if(quantity[0] == 't')
        {
            units = units.empty() ? "s" : units;
            return units::convert(t, "s", units);
        }
        if(quantity[0] == 'n')
        {
            units = units.empty() ? "mol" : units;
            auto name = split(quantity, "[]").back();
            const double ni = state.speciesAmount(name);
            return units::convert(ni, "mol", units);
        }
        if(quantity[0] == 'b')
        {
            units = units.empty() ? "mol" : units;
            auto names = split(quantity, "[]");
            std::string element = names[1];
            std::string phase = names.size() > 2 ? names[2] : "";
            const double bi = phase.empty() ?
                state.elementAmount(element, units) :
                state.elementAmountInPhase(element, phase, units);
            return bi;
        }
        if(quantity[0] == 'x')
        {
            auto name = split(quantity, "[]").back();
            auto ispecies = system.indexSpeciesWithError(name);
            auto iphase = system.indexPhaseWithSpecies(ispecies);
            auto ifirst = system.indexFirstSpeciesInPhase(iphase);
            auto size = system.numSpeciesInPhase(iphase);
            auto nt = sum(rows(n, ifirst, size));
            auto ni = n[ispecies];
            auto xi = nt ? ni/nt : 0.0;
            return xi;
        }
        if(quantity[0] == 'm')
        {
            units = units.empty() ? "molal" : units;
            std::string name = split(quantity, "[]").back();
            auto ispecies = system.indexSpecies(name);
            auto ielement = system.indexElement(name);
            Assert(ielement < system.numElements() || ispecies < system.numSpecies(),
                "Cannot calculate the molality of `" + name + "`.",
                "There is no species or element in the chemical system with such name.");
            auto amount = ispecies < system.numSpecies()  ?
                state.speciesAmount(name) :
                state.elementAmountInPhase(name, "Aqueous");
            const double nH2O = state.speciesAmount("H2O(l)");
            const double mi = nH2O ? amount/(nH2O * waterMolarMass) : 0.0;
            return units::convert(mi, "molal", units);
        }
        if(quantity[0] == 'a')
        {
            std::string name = split(quantity, "[]").back();
            Index index = system.indexSpecies(name);
            const double ln_ai = properties.lnActivities().val[index];
            return std::exp(ln_ai);
        }
        if(quantity[0] == 'g')
        {
            std::string name = split(quantity, "[]").back();
            Index index = system.indexSpecies(name);
            const double ln_gi = properties.lnActivityCoefficients().val[index];
            return std::exp(ln_gi);
        }
        if(quantity == "pH")
        {
            const Index iH = system.indexSpecies("H+");
            const double ln_aH = properties.lnActivities().val[iH];
            return -ln_aH/std::log(10);
        }
        if(quantity[0] == 'r')
        {
            // Return zero if there are no reactions
            if(r.val.rows() == 0) return 0.0;
            units = units.empty() ? "mol/s" : units;
            std::string reaction = split(quantity, "[]").back();
            Index index = reactions.indexReactionWithError(reaction);
            const double ri = r.val[index];
            return units::convert(ri, "mol/s", units);
        }

        RuntimeError("Cannot calculate the chemical quantity `" + str + "`.",
            "The formatted string `" + str + "` does not contain a valid quantity.");

        return 0.0;
    }
};

ChemicalQuantity::ChemicalQuantity()
: pimpl(new Impl())
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalQuantity& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalQuantity::ChemicalQuantity(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

ChemicalQuantity::~ChemicalQuantity()
{}

auto ChemicalQuantity::operator=(ChemicalQuantity other) -> ChemicalQuantity&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalQuantity::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto ChemicalQuantity::update(const ChemicalState& state, double t) -> void
{
    pimpl->update(state, t);
}

auto ChemicalQuantity::value(std::string str) const -> double
{
    return pimpl->value(str);
}

auto ChemicalQuantity::operator[](std::string quantity) const -> double
{
    return value(quantity);
}

} // namespace Reaktoro
