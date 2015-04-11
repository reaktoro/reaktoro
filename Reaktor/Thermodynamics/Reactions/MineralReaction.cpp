// Reaktor is a C++ library for computational reaction modelling.
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

#include "MineralReaction.hpp"

// C++ includes
#include <math.h>

// Reaktor includes
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/ReactionEquation.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Common/StringUtils.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Reaction.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Utils.hpp>

namespace Reaktor {
namespace internal {

using MineralCatalystFunction = std::function<ChemicalScalar(double, double, const Vector&, const ChemicalVector&)>;
using MineralMechanismFunction = std::function<ChemicalScalar(double, double, const Vector&, const ChemicalVector&)>;

auto mineralCatalystFunctionActivity(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    const double power = catalyst.power;
    const std::string species = catalyst.species;
    const Index ispecies = system.indexSpeciesWithError(species);
    ChemicalScalar ai, res;
    ChemicalVector ln_a;

    MineralCatalystFunction fn = [=](double T, double P, const Vector& n, const ChemicalVector& a) mutable
    {
        ai = a.row(ispecies);
        return pow(ai, power);
    };

    return fn;
}

auto mineralCatalystFunctionPartialPressure(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    const auto gas         = catalyst.species;                         // the species of the catalyst
    const auto power       = catalyst.power;                           // the power of the catalyst
    const auto idx_phase   = system.indexPhase("Gaseous");             // the index of the gaseous phase
    const auto gases       = names(system.phase(idx_phase).species()); // the names of the gaseous species
    const auto igases      = system.indicesSpecies(gases);             // the indices of the gaseous species
    const auto igas        = index(gas, gases);                        // the index of the gaseous species
    const auto num_species = system.numSpecies();                      // the number of species
    const auto num_gases   = gases.size();                             // the number of gases

    Vector ng;
    Vector dxidng;
    Vector dxidn;
    ChemicalScalar res;

    MineralCatalystFunction fn = [=](double T, double P, const Vector& n, const ChemicalVector& a) mutable
    {
        // The molar composition of the gaseous species
        ng = rows(n, igases);

        // The total number of moles in the gaseous phase
        const double ngsum = ng.sum();

        // The molar fraction of the gas
        const double xi = ng[igas]/ngsum;

        // The derivative of the molar fraction of the gas with respect to the all gaseous species
        dxidng.noalias() = (Vector::Unit(num_gases, igas) - Vector::Constant(num_gases, xi))/ngsum;

        // The derivative of the molar fraction of the gas with respect to the all species
        dxidn = zeros(num_species);
        rows(dxidn, igases) = dxidng;

        // The pressure in units of bar
        P = convert<Pa,bar>(P);

        // Evaluate the mineral catalyst function
        res.val = std::pow(xi * P, power);
        res.ddn = res.val * power/xi * dxidn;

        return res;
    };

    return fn;
}

auto mineralCatalystFunction(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    if (catalyst.quantity == "a" || catalyst.quantity == "activity")
        return mineralCatalystFunctionActivity(catalyst, system);
    else
        return mineralCatalystFunctionPartialPressure(catalyst, system);
}

auto mineralMechanismFunction(const MineralMechanism& mechanism, const Reaction& reaction, const ChemicalSystem& system) -> ReactionRateFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();

    // The universal gas constant (in units of kJ/(mol*K))
    const double R = 8.3144621e-3;

    // Create the mineral catalyst functions
    std::vector<MineralCatalystFunction> catalysts;
    for(const MineralCatalyst& catalyst : mechanism.catalysts)
        catalysts.push_back(mineralCatalystFunction(catalyst, system));

    // Auxiliary variables
    ChemicalScalar aux, f, g;

    // Define the mineral mechanism function
    ReactionRateFunction fn = [=](double T, double P, const Vector& n, const ChemicalVector& a) mutable
    {
        // The result of this function evaluation
        ChemicalScalar res(num_species);

        // Calculate the saturation index of the mineral
        ChemicalScalar lnOmega = reaction.lnReactionQuotient(a) - reaction.lnEquilibriumConstant(T, P);

        // Calculate the rate constant for the current mechanism
        const double kappa = mechanism.kappa * std::exp(-mechanism.Ea/R * (1.0/T - 1.0/298.15));

        // Calculate the saturation index
        const double Omega = std::exp(lnOmega.val);

        // Calculate the p and q powers of the saturation index Omega
        const double pOmega = std::pow(Omega, mechanism.p);
        const double qOmega = std::pow(1 - pOmega, mechanism.q);

        // Calculate the function f
        f.val = kappa * qOmega;
        f.ddn = -mechanism.p * mechanism.q * pOmega/(1.0 - pOmega) * f.val * lnOmega.ddn;

        // Calculate the function g
        g.val = 1.0;
        g.ddn = zeros(num_species);

        for(const MineralCatalystFunction& catalyst : catalysts)
        {
            aux = catalyst(T, P, n, a);
            g.val *= aux.val;
            g.ddn += aux.ddn/aux.val;
        }

        g.ddn *= g.val;

        // Calculate the resulting mechanism function
        res.val = f.val * g.val;
        res.ddn = f.val * g.ddn + f.ddn * g.val;

        return res;
    };

    return fn;
}

inline auto surfaceAreaUnitError(const std::string& unit) -> void
{
    Exception exception;
    exception.error << "Cannot set the specific surface area of the mineral reaction";
    exception.reason << "The provided specific surface area unit " << unit << " cannot be converted to m2/g or m2/m3";
    RaiseError(exception);
}

inline auto errroZeroSurfaceArea(const MineralReaction& reaction) -> void
{
    Exception exception;
    exception.error << "Cannot calculate the molar surface area of the mineral " << reaction.mineral() << ".";
    exception.reason << "The specific surface area of the mineral was not set in reaction " << reaction.equation() << ".";
    RaiseError(exception);
}

} /* namespace internal */

using namespace internal;

class MineralReaction::Impl
{
private:
    /// The name of the mineral species
    std::string mineral$;

    /// The equation of the mineral reaction
    ReactionEquation equation$;

    /// The equilibrium constant of the mineral reaction
    ThermoScalarFunction lnk;

    /// The volumetric surface area of the mineral
    double volumetric_surface_area$;

    /// The specific surface area of the mineral
    double specific_surface_area$;

    /// The mineral rate mechanisms of the mineral dismixture/precipitation equation
    std::vector<MineralMechanism> mechanisms$;

public:
    Impl()
    {}

    Impl(const std::string& mineral)
    : mineral$(mineral)
    {}

    auto setMineral(const std::string& mineral) -> void
    {
        mineral$ = mineral;
    }

    auto setEquation(const ReactionEquation& equation) -> void
    {
        equation$ = equation;
    }

    auto setEquation(const std::string& equation) -> void
    {
        equation$ = ReactionEquation(equation);
    }

    auto setEquilibriumConstant(const ThermoScalarFunction& lnk) -> void
    {
        this->lnk = lnk;
    }

    auto setSpecificSurfaceArea(double value, const std::string& unit) -> void
    {
        // Reset both specific and volumetric surface area instances
        specific_surface_area$ = 0.0;
        volumetric_surface_area$ = 0.0;

        // Check the appropriate unit of the surface area and set the corresponding variable
        if(units::convertible(unit, "m2/kg"))
            specific_surface_area$ = units::convert(value, unit, "m2/kg");
        else if(units::convertible(unit, "m2/m3"))
            volumetric_surface_area$ = units::convert(value, unit, "m2/kg");
        else surfaceAreaUnitError(unit);
    }

    auto addMechanism(const std::string& mechanism) -> void
    {
        addMechanism(MineralMechanism(mechanism));
    }

    auto addMechanism(const MineralMechanism& mechanism) -> void
    {
        mechanisms$.push_back(mechanism);
    }

    auto setMechanisms(const std::vector<MineralMechanism>& mechanisms) -> void
    {
        mechanisms$ = mechanisms;
    }

    auto mineral() const -> const std::string&
    {
        return mineral$;
    }

    auto equation() const -> const ReactionEquation&
    {
        return equation$;
    }

    auto equilibriumConstant() const -> const ThermoScalarFunction&
    {
        return lnk;
    }

    auto specificSurfaceArea() const -> double
    {
        return specific_surface_area$;
    }

    auto volumetricSurfaceArea() const -> double
    {
        return volumetric_surface_area$;
    }

    auto mechanisms() const -> const std::vector<MineralMechanism>&
    {
        return mechanisms$;
    }
};

MineralReaction::MineralReaction()
: pimpl(new Impl())
{}

MineralReaction::MineralReaction(const std::string& mineral)
: pimpl(new Impl(mineral))
{}

MineralReaction::MineralReaction(const MineralReaction& other)
: pimpl(new Impl(*other.pimpl))
{}

MineralReaction::~MineralReaction()
{}

auto MineralReaction::operator=(MineralReaction other) -> MineralReaction&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto MineralReaction::setMineral(const std::string& mineral) -> MineralReaction&
{
    pimpl->setMineral(mineral);
    return *this;
}

auto MineralReaction::setEquation(const ReactionEquation& equation) -> MineralReaction&
{
    pimpl->setEquation(equation);
    return *this;
}

auto MineralReaction::setEquation(const std::string& equation) -> MineralReaction&
{
    pimpl->setEquation(equation);
    return *this;
}

auto MineralReaction::setEquilibriumConstant(const ThermoScalarFunction& lnk) -> MineralReaction&
{
    pimpl->setEquilibriumConstant(lnk);
    return *this;
}

auto MineralReaction::setSpecificSurfaceArea(double value, const std::string& unit) -> MineralReaction&
{
    pimpl->setSpecificSurfaceArea(value, unit);
    return *this;
}

auto MineralReaction::addMechanism(const std::string& mechanism) -> MineralReaction&
{
    pimpl->addMechanism(mechanism);
    return *this;
}

auto MineralReaction::addMechanism(const MineralMechanism& mechanism) -> MineralReaction&
{
    pimpl->addMechanism(mechanism);
    return *this;
}

auto MineralReaction::setMechanisms(const std::vector<MineralMechanism>& mechanisms) -> MineralReaction&
{
    pimpl->setMechanisms(mechanisms);
    return *this;
}

auto MineralReaction::mineral() const -> const std::string&
{
    return pimpl->mineral();
}

auto MineralReaction::equation() const -> const ReactionEquation&
{
    return pimpl->equation();
}

auto MineralReaction::equilibriumConstant() const -> const ThermoScalarFunction&
{
    return pimpl->equilibriumConstant();
}

auto MineralReaction::specificSurfaceArea() const -> double
{
    return pimpl->specificSurfaceArea();
}

auto MineralReaction::volumetricSurfaceArea() const -> double
{
    return pimpl->volumetricSurfaceArea();
}

auto MineralReaction::mechanisms() const -> const std::vector<MineralMechanism>&
{
    return pimpl->mechanisms();
}

/// Calculate the molar surface area of a mineral (in units of m2/mol)
/// @param reaction The mineral reaction instance
/// @param system The chemical system instance
/// @return The molar surface area of the mineral (in units of m2/mol)
auto molarSurfaceArea(const MineralReaction& reaction, const ChemicalSystem& system) -> double
{
    // The temperature and pressure for the calculation of the mineral density
    // Note: These values do not matter much, since the density of the minerals is a constant function
    const double T = 298.15; // in units of kelvin
    const double P = 1.0e5;  // in units of pascal

    // The specific surface area of the mineral (in units of m2/kg)
    const double specific_surface_area = reaction.specificSurfaceArea();

    // The molar mass of the mineral (in units of kg/mol)
    const double molar_mass = system.species(reaction.mineral()).molarMass();

    // Check if the specific surface area of the mineral was set
    if(specific_surface_area) return specific_surface_area * molar_mass;

    // The volumetric surface area of the mineral (in units of m2/m3)
    const double volumetric_surface_area = reaction.volumetricSurfaceArea();

    // The molar volume of the mineral species (in units of m3/mol)
    const double molar_volume = system.species(reaction.mineral()).standardVolume(T, P).val;

    // Check if the volumetric surface area of the mineral was set
    if(volumetric_surface_area) return volumetric_surface_area * molar_volume;

    errroZeroSurfaceArea(reaction);

    return 0.0;
}

auto createReaction(const MineralReaction& mineralrxn, const ChemicalSystem& system) -> Reaction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();

    // The index of the mineral
    const Index imineral = system.indexSpeciesWithError(mineralrxn.mineral());

    // The molar surface area of the mineral
    const double molar_surface_area = molarSurfaceArea(mineralrxn, system);

    // Create a reaction instance
    Reaction reaction(mineralrxn.equation(), system);

    // Set the name of the reaction
    reaction.setName(mineralrxn.mineral());

    // Check if an equilibrium constant was provided to the mineral reaction
    if(mineralrxn.equilibriumConstant())
        reaction.setEquilibriumConstantFunction(mineralrxn.equilibriumConstant());

    // Create the mineral mechanism functions
    std::vector<ReactionRateFunction> mechanisms;
    for(const MineralMechanism& mechanism : mineralrxn.mechanisms())
        mechanisms.push_back(mineralMechanismFunction(mechanism, reaction, system));

    // Create the mineral rate function
    ReactionRateFunction rate = [=](double T, double P, const Vector& n, const ChemicalVector& a)
    {
        // The number of moles of the mineral
        const double nm = n[imineral];

        // Calculate the reactive surface area of the mineral
        const double sa = molar_surface_area * nm;

        // The reaction rate and its partial molar derivatives
        ChemicalScalar rate(num_species);

        // Iterate over all mechanism functions
        for(const ReactionRateFunction& mechanism : mechanisms)
            rate += sa * mechanism(T, P, n, a);

        rate.ddn[imineral] += rate.val/nm;

        return rate;
    };

    // Set the rate of the reaction
    reaction.setRate(rate);

    return reaction;
}

} // namespace Reaktor
