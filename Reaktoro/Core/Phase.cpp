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

#include "Phase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {
namespace {

auto defaultMolarVolumeFunction(const std::vector<Species>& species) -> ChemicalScalarFunction
{
    ChemicalScalarFunction fn = [=](double T, double P, const Vector& n) -> ChemicalScalar
    {
        ChemicalScalar res;
        const double nt = sum(n);
        if(nt == 0.0) return res;
        for(unsigned i = 0; i < n.size(); ++i)
            res += n[i]/nt * species[i].standardVolume(T, P);
        return res;
    };
    return fn;
}

} // namespace

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The list of Species instances defining the phase
    std::vector<Species> species;

    /// The list of Element instances in the phase
    std::vector<Element> elements;

    /// The function for the concentrations of the species (no uniform units).
    ChemicalVectorFunction concentration_fn;

    /// The function for the natural log of the activity coefficients of the species.
    ChemicalVectorFunction activity_coefficient_fn;

    /// The function for the natural log of the activities of the species.
    ChemicalVectorFunction activity_fn;

    /// The function for the molar volume of the phase (in units of m3/mol).
    ChemicalScalarFunction molar_volume_fn;
};

Phase::Phase()
: pimpl(new Impl())
{}

Phase::Phase(const Phase& other)
: pimpl(new Impl(*other.pimpl))
{}

Phase::~Phase()
{}

auto Phase::operator=(Phase other) -> Phase&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Phase::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> void
{
    pimpl->species = species;

    if(not pimpl->molar_volume_fn)
        setMolarVolumeFunction(defaultMolarVolumeFunction(species));
}

auto Phase::setConcentrationFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->concentration_fn = function;
}

auto Phase::setActivityCoefficientFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->activity_coefficient_fn = function;
}

auto Phase::setActivityFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->activity_fn = function;
}

auto Phase::setMolarVolumeFunction(const ChemicalScalarFunction& function) -> void
{
    pimpl->molar_volume_fn = function;
}

auto Phase::numElements() const -> unsigned
{
    return elements().size();
}

auto Phase::numSpecies() const -> unsigned
{
    return species().size();
}

auto Phase::name() const -> std::string
{
    return pimpl->name;
}

auto Phase::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::species(Index index) const -> const Species&
{
    return pimpl->species[index];
}

auto Phase::concentrationFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->concentration_fn;
}

auto Phase::activityCoefficientFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->activity_coefficient_fn;
}

auto Phase::activityFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->activity_fn;
}

auto Phase::molarVolumeFunction() const -> const ChemicalScalarFunction&
{
    return pimpl->molar_volume_fn;
}

auto Phase::standardGibbsEnergies(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardGibbsEnergy(T, P);
    return res;
}

auto Phase::standardEnthalpies(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardEnthalpy(T, P);
    return res;
}

auto Phase::standardHelmholtzEnergies(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardHelmholtzEnergy(T, P);
    return res;
}

auto Phase::standardEntropies(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardEntropy(T, P);
    return res;
}

auto Phase::standardVolumes(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardVolume(T, P);
    return res;
}

auto Phase::standardInternalEnergies(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardInternalEnergy(T, P);
    return res;
}

auto Phase::standardHeatCapacities(double T, double P) const -> ThermoVector
{
    const unsigned num_species = numSpecies();
    ThermoVector res(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        res.row(i) = species(i).standardHeatCapacity(T, P);
    return res;
}

auto Phase::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->concentration_fn(T, P, n);
}

auto Phase::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->activity_coefficient_fn(T, P, n);
}

auto Phase::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    return pimpl->activity_fn(T, P, n);
}

auto Phase::chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector
{
    const double R = universalGasConstant;
    ThermoScalar RT = R*ThermoScalar(T, 1.0, 0.0);
    ThermoVector u0 = standardGibbsEnergies(T, P);
    ChemicalVector ln_a = log(activities(T, P, n));
    ChemicalVector u = u0 + RT*ln_a;
    return u;
}

auto Phase::molarVolume(double T, double P, const Vector& n) const -> ChemicalScalar
{
    return pimpl->molar_volume_fn(T, P, n);
}

auto operator<(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
