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
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Reaction.hpp>
#include <Reaktor/Utils/ConvertUtils.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>
#include <Reaktor/Utils/SetUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {
namespace internal {

using MineralCatalystFunction = std::function<PartialScalar(double T, double P, const Vector& n, const PartialVector& a)>;

using MineralMechanismFunction = std::function<PartialScalar(double T, double P, const Vector& n, const PartialVector& a, const PartialScalar& omega)>;

auto createMineralCatalystFunction(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction;

auto createMineralMechanismFunction(const MineralMechanism& mechanism, const ChemicalSystem& system) -> MineralMechanismFunction;

inline auto surfaceAreaUnitError(const std::string& unit) -> void
{
    Exception exception;
    exception.error << "Cannot set the specific surface area of the mineral reaction";
    exception.reason << "The provided specific surface area unit " << unit << " cannot be converted to m2/g or m2/m3";
    raise(exception);
}

inline auto zeroSurfaceAreaError(const MineralReaction& reaction) -> void
{
    Exception exception;
    exception.error << "Cannot calculate the molar surface area of the mineral " << reaction.mineral() << ".";
    exception.reason << "The specific surface area of the mineral was not set in reaction " << reaction.equation() << ".";
    raise(exception);
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
    ThermoScalarFunction equilibrium_constant$;

    /// The volumetric surface area of the mineral
    units::VolumetricSurfaceArea volumetric_surface_area$;

    /// The specific surface area of the mineral
    units::SpecificSurfaceArea specific_surface_area$;

    /// The mineral rate mechanisms of the mineral dissolution/precipitation equation
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

    auto setEquilibriumConstant(const ThermoScalarFunction& equilibrium_constant) -> void
    {
        equilibrium_constant$ = equilibrium_constant;
    }

    auto setSpecificSurfaceArea(double value, const std::string& unit) -> void
    {
        // Reset both specific and volumetric surface area instances
        specific_surface_area$ = 0.0;
        volumetric_surface_area$ = 0.0;

        // Check the appropriate unit of the surface area and set the corresponding variable
        if(units::convertible(unit, "m2/g"))
            specific_surface_area$ = units::SpecificSurfaceArea(value, unit);
        else if(units::convertible(unit, "m2/m3"))
            volumetric_surface_area$ = units::VolumetricSurfaceArea(value, unit);
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

    auto equilibriumConstant() const -> const EquilibriumConstant&
    {
        return equilibrium_constant$;
    }

    auto specificSurfaceArea() const -> units::SpecificSurfaceArea
    {
        return specific_surface_area$;
    }

    auto volumetricSurfaceArea() const -> units::VolumetricSurfaceArea
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

auto MineralReaction::setEquilibriumConstant(const ThermoScalarFunction& equilibrium_constant) -> MineralReaction&
{
    pimpl->setEquilibriumConstant(equilibrium_constant);
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

auto MineralReaction::equilibriumConstant() const -> const EquilibriumConstant&
{
    return pimpl->equilibriumConstant();
}

auto MineralReaction::specificSurfaceArea() const -> units::SpecificSurfaceArea
{
    return pimpl->specificSurfaceArea();
}

auto MineralReaction::volumetricSurfaceArea() const -> units::VolumetricSurfaceArea
{
    return pimpl->volumetricSurfaceArea();
}

auto MineralReaction::mechanisms() const -> const std::vector<MineralMechanism>&
{
    return pimpl->mechanisms();
}

/**
 * Calculates the molar surface area of a mineral (in units of m2/mol)
 * @param reaction The mineral reaction instance
 * @param system The chemical system instance
 * @return The molar surface area of the mineral (in units of m2/mol)
 */
auto molarSurfaceArea(const MineralReaction& reaction, const ChemicalSystem& system) -> double
{
    // The temperature and pressure for the calculation of the mineral density
    // Note: These values do not matter much, since the density of the minerals is a constant function
    const double T = (25.0 * unit(celsius))();
    const double P = (1.0 * unit(bar))();

    // The specific surface area of the mineral (in units of m2/kg)
    const double specific_surface_area = reaction.specificSurfaceArea().in(unit(m2)/unit(kg));

    // The molar mass of the mineral (in units of kg/mol)
    const double molar_mass = system.species(reaction.mineral()).molarMass().in(unit(kg)/unit(mol));

    // Check if the specific surface area of the mineral was set
    if(specific_surface_area) return specific_surface_area * molar_mass;

    // The volumetric surface area of the mineral (in units of m2/m3)
    const double volumetric_surface_area = reaction.volumetricSurfaceArea().in(unit(m2)/unit(m3));

    // The molar density of the mineral species (in units of mol/m3)
    const double molar_density = system.species(reaction.mineral()).molarDensity(T, P);

    // Check if the volumetric surface area of the mineral was set
    if(volumetric_surface_area) return volumetric_surface_area / molar_density;

    zeroSurfaceAreaError(reaction);

    return 0.0;
}

auto createReaction(const MineralReaction& reaction, const ChemicalSystem& system) -> Reaction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();

    // The index of the mineral
    const Index idx_mineral = system.idxSpeciesWithError(reaction.mineral());

    // The molar surface area of the mineral
    const double molar_surface_area = molarSurfaceArea(reaction, system);

    // Create a reaction instance
    Reaction converted(reaction.equation(), system);

    // Check if an equilibrium constant was provided to the mineral reaction
    if(reaction.equilibriumConstant())
        converted.setEquilibriumConstant(reaction.equilibriumConstant());

    // Create the mineral mechanism functions
    std::vector<MineralMechanismFunction> mechanisms;
    for(const MineralMechanism& mechanism : reaction.mechanisms())
        mechanisms.push_back(createMineralMechanismFunction(mechanism, system));

    // Create the mineral rate function
    ReactionRate rate = [=](double T, double P, const Vector& n, const PartialVector& a)
    {
        // Calculate the equilibrium constant of the mineral reaction
        const double K = converted.equilibriumConstant(T, P);

        // Calculate the saturation index of the mineral
        PartialScalar omega = converted.reactionQuotient(a);
        func(omega) /= K;
        grad(omega) /= K;

        // The number of moles of the mineral
        const double nm = n[idx_mineral];

        // Calculate the reactive surface area of the mineral
        const double sa = molar_surface_area * nm;

        // The reaction rate and its partial molar derivatives
        PartialScalar rate = partialScalar(0.0, zeros(num_species));

        // Iterate over all mechanism functions
        for(const MineralMechanismFunction& mechanism : mechanisms)
        {
            PartialScalar aux = mechanism(T, P, n, a, omega);
            func(rate) += sa * func(aux);
            grad(rate) += sa * grad(aux);
        }

        grad(rate)[idx_mineral] += func(rate)/nm;

        return rate;
    };

    // Set the rate of the reaction
    converted.setRate(rate);

    return converted;
}

namespace internal {

auto mineralCatalystActivity(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    const auto& species = catalyst.species;
    const auto& power = catalyst.power;
    const Index idx_species = system.idxSpeciesWithError(species);
    PartialScalar res;

    MineralCatalystFunction fn = [=](double T, double P, const Vector& n, const PartialVector& a) mutable
    {
        // Evaluate the mineral catalyst function
        func(res) = std::pow(func(a)[idx_species], power);
        grad(res) = func(res) * power/func(a)[idx_species] * grad(a).row(idx_species);

        return res;
    };

    return fn;
}

auto mineralCatalystPartialPressure(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    const auto gas         = catalyst.species;                       // the species of the catalyst
    const auto power       = catalyst.power;                         // the power of the catalyst
    const auto idx_phase   = system.idxPhase("Gaseous");             // the index of the gaseous phase
    const auto gases       = system.phase(idx_phase).speciesNames(); // the names of the gaseous species
    const auto idx_gases   = system.idxSpecies(gases);               // the indices of the gaseous species
    const auto idx_gas     = find(gas, gases);                       // the index of the gaseous species
    const auto num_species = system.numSpecies();                    // the number of species
    const auto num_gases   = gases.size();                           // the number of gases

    Vector ng;
    Vector dxidng;
    Vector dxidn;
    PartialScalar res;

    MineralCatalystFunction fn = [=](double T, double P, const Vector& n, const PartialVector& a) mutable
    {
        // The molar composition of the gaseous species
        ng = rows(idx_gases, n);

        // The total number of moles in the gaseous phase
        const double ngsum = ng.sum();

        // The molar fraction of the gas
        const double xi = ng[idx_gas]/ngsum;

        // The derivative of the molar fraction of the gas with respect to the all gaseous species
        dxidng.noalias() = (Vector::Unit(num_gases, idx_gas) - Vector::Constant(num_gases, xi))/ngsum;

        // The derivative of the molar fraction of the gas with respect to the all species
        dxidn = zeros(num_species);
        setRows(idx_gases, dxidng, dxidn);

        // The pressure in units of bar
        P = convert<Pa,bar>(P);

        // Evaluate the mineral catalyst function
        func(res) = std::pow(xi * P, power);
        grad(res) = func(res) * power/xi * dxidn;

        return res;
    };

    return fn;
}

auto createMineralCatalystFunction(const MineralCatalyst& catalyst, const ChemicalSystem& system) -> MineralCatalystFunction
{
    if(catalyst.quantity == "a" or catalyst.quantity == "activity")
        return mineralCatalystActivity(catalyst, system);
    else
        return mineralCatalystPartialPressure(catalyst, system);
}

auto createMineralMechanismFunction(const MineralMechanism& mechanism, const ChemicalSystem& system) -> MineralMechanismFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();

    // The universal gas constant (in units of kJ/(mol*K))
    const double R = 8.3144621e-3;

    // Create the mineral catalyst functions
    std::vector<MineralCatalystFunction> catalysts;
    for(const MineralCatalyst& catalyst : mechanism.catalysts)
        catalysts.push_back(createMineralCatalystFunction(catalyst, system));

    // Auxiliary variables
    PartialScalar aux, f, g;

    // Define the mineral mechanism function
    MineralMechanismFunction fn = [=](double T, double P, const Vector& n, const PartialVector& a, const PartialScalar& omega) mutable
    {
        // The result of this function evaluation
        PartialScalar res = partialScalar(0.0, zeros(num_species));

        // Calculate the rate constant for the current mechanism
        const double kappa = mechanism.kappa * std::exp(-mechanism.Ea/R * (1.0/T - 1.0/298.15));

        // Calculate the p and q powers of the saturation index Omega
        const double pOmega = std::pow(func(omega), mechanism.p);
        const double qOmega = std::pow(1 - pOmega, mechanism.q);

        // Calculate the function f
        func(f) = kappa * qOmega;
        grad(f) = -mechanism.p * mechanism.q * pOmega/(1.0 - pOmega) * func(f)/func(omega) * grad(omega);

        // Calculate the function g
        func(g) = 1.0;
        grad(g).setZero(num_species);

        for(const MineralCatalystFunction& catalyst : catalysts)
        {
            aux = catalyst(T, P, n, a);
            func(g) *= func(aux);
            grad(g) += grad(aux)/func(aux);
        }

        grad(g) *= func(g);

        // Calculate the resulting mechanism function
        func(res) = func(f) * func(g);
        grad(res) = func(f) * grad(g) + grad(f) * func(g);

        return res;
    };

    return fn;
}

} // namespace internal
} // namespace Reaktor
