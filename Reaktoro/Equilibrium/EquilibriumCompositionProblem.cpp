// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "EquilibriumCompositionProblem.hpp"

// C++ includes
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/PhaseThermoProperties.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

struct EquilibriumCompositionProblem::Impl
{
    /// The chemical system.
    ChemicalSystem system;

    /// The partition of the chemical system.
    Partition partition;

    /// The titrants used to fix the saturations of fluid phases
    std::map<std::string, std::string> titrants_fluid_phases;

    /// The titrant used to fix the porosity of the solid matrix
    std::string titrant_solid;

    /// The volume fractions of fluid phases
    std::map<std::string, double> saturation_fluid_phases;

    /// The volume fractions of solid phases
    std::map<std::string, double> volume_fractions_solid_phases;

    /// The temperature for the equilibrium calculations (in units of K)
    double T = 298.15;

    /// The pressure for the equilibrium calculations (in units of Pa)
    double P = 1e5;

    /// The porosity of the solid matrix
    double porosity = 1.0;

    /// Construct a default EquilibriumCompositionProblem instance.
    Impl()
    {}

    /// Construct a custom EquilibriumCompositionProblem instance.
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;
    }

    /// Set the composition of the aqueous phase using molalities of compounds.
    auto setAqueousFluid(std::string molalities) -> void
    {
        // Add 1 kg by to the titrant used to set the aqueous phase
        std::string aqueous_titrant = "1 kg H2O;";

        // Split molalities in triplets (value, molal_units, compound)
        for(auto compound : split(molalities, ";"))
        {
            // Get the current triplet (value, molal_units, compound)
            auto words = split(compound);

            // Convert the given molality value to molal
            auto amount = units::convert(tofloat(words[0]), words[1], "molal");

            // Update aqueous_titrant with current compound
            aqueous_titrant += " " + std::to_string(amount) + " mol " + words[2] + ";";
        }

        // Set the titrant of the aqueous phase
        titrants_fluid_phases["Aqueous"] = aqueous_titrant;
    }

    /// Set the composition of the gaseous phase using molar fractions of compounds.
    auto setGaseousFluid(std::string molarfractions) -> void
    {
        // The gaseous titrant
        std::string gaseous_titrant;

        // The sum of the molar fractions
        double sum = 0.0;

        // Split molarfractions in pairs (molarfraction, gasname)
        for(auto compound : split(molarfractions, ";"))
        {
            // Get the current pair (molarfraction, gasname)
            auto words = split(compound);

            // Check if only gas name is provided
            if(words.size() == 1)
            {
                gaseous_titrant += " 1 mol " + words[0] + ";";
                sum += 1.0;
            }
            // Check if both molar fraction and gas name are provided
            else if(words.size() == 2)
            {
                gaseous_titrant += " " + words[0] + " mol " + words[1] + ";";
                sum += tofloat(words[0]);
            }
            // None of above - throw runtime error
            else RuntimeError("Could not set the composition of the gaseous phase.",
                "Expecting either a list of compounds like `0.80 CO2; 0.20 H2O` "
                "or a single compound like `CO2`");
        }

        // Assert the provided molar fractions sum to one
        Assert(std::abs(sum-1.0) < 1e-14, "Could not set the composition of the gaseous phase.",
            "Expecting a list of compounds whose molar fractions sum to one.");

        // Set the titrant of the gaseous phase
        titrants_fluid_phases["Gaseous"] = gaseous_titrant;
    }

    /// Set the volume fractions of the solid phases.
    auto setSolid(std::string volumefractions) -> void
    {
        // The indices of the solid phases
        const Indices& isp = partition.indicesSolidPhases();

        // The molar volumes of the solid phases at (T, P), defined as
        // the average of the standard molar volumes of their end-members.
        Vector molar_volumes_solid_phases(isp.size());
        for(Index i = 0; i < isp.size(); ++i)
            molar_volumes_solid_phases[i] = sum(system.phase(isp[i]).properties(T, P).
                standardPartialMolarVolumes().val)/isp.size();

        // The sum of the volume fractions
        double sum = 0.0;

        // Clear the current volume fractions of solid phases
        volume_fractions_solid_phases.clear();

        // Clear the current titrant for the solid volume constraint
        titrant_solid.clear();

        // Split volumefractions in pairs (fraction phasename)
        for(auto pair : split(volumefractions, ";"))
        {
            // Split current (fraction phasename) pair
            auto words = split(pair);

            // Get the index of the phase and its solid phase index
            const Index iphase = system.indexPhase(words[1]);
            const Index isolidphase = index(iphase, isp);

            // Check if there is any solid phase with given name
            Assert(isolidphase < isp.size(), "Could not set the solid composition.",
                "There is no solid phase named `" << words[1] << "`.");

            // The volume fraction of the current phase
            const double vphase = tofloat(words[0]);

            // The molar volume of the current phase
            const double molar_volume = molar_volumes_solid_phases[isolidphase];

            // The number of species in current phase
            const Index size = system.numSpeciesInPhase(iphase);

            // Update the sum of volume fractions
            sum += vphase;

            // Set the volume fraction of the solid phase
            volume_fractions_solid_phases[words[1]] = vphase;

            // The amount of current solid phase
            const double amount = vphase/(molar_volume*size);

            // Loop over all species in the current solid phase
            for(const Species& species : system.phase(iphase).species())
                titrant_solid += " " + std::to_string(amount) + " mol " + species.formula() + ";";
        }

        // Assert the provided volume fractions sum to one
        Assert(std::abs(sum - 1.0) < 1e-14, "Could not set the volume fractions of solid phases.",
            "Expecting a list of solid phases whose volume fractions sum to one.");
    }

    /// Return the equilibrium problem defining the equilibrium of fluid and solid phases.
    auto equilibriumInverseProblem() -> EquilibriumInverseProblem
    {
        // Assert the sum of fluid phase saturations is one.
        assertSaturations();

        // The indices of the solid phases
        const Indices& isp = partition.indicesSolidPhases();

        // The names of the solid phases
        auto names_solid_phases = extract(names(system.phases()), isp);

        // Create an equilibrium problem
        EquilibriumInverseProblem problem(system);
        problem.setPartition(partition);
        problem.setTemperature(T);
        problem.setPressure(P);

        // Set the volume fractions of fluid phases
        for(auto pair : saturation_fluid_phases)
        {
            // Assert the composition of current fluid phase was given
            Assert(titrants_fluid_phases.count(pair.first),
                "Could not equilibrate fluid and solid phases.",
                "No composition was given for fluid phase `" + pair.first + "`.");

            // Fix the fluid phase volume in the equilibrium calculation
            problem.fixPhaseVolume(pair.first, porosity * pair.second, "m3",
                titrants_fluid_phases.at(pair.first));
        }

        // Add the constraint on the sum of solid phase volumes
        problem.fixPhaseSetVolume(names_solid_phases, 1 - porosity, "m3", titrant_solid);

        return problem;
    }

    // Assert the sum of fluid phase saturations is one.
    void assertSaturations()
    {
        double saturation_sum = 0.0;
        for(auto iter : saturation_fluid_phases)
            saturation_sum += iter.second;
        Assert(std::abs(saturation_sum - 1.0) < 1e-14,
            "Could not equilibrate fluid and solid phases.",
            "Expecting the sum of fluid phase saturations to be one, "
            "but got " + std::to_string(saturation_sum) + " instead.");
    }
};

EquilibriumCompositionProblem::EquilibriumCompositionProblem()
: pimpl(new Impl())
{}

EquilibriumCompositionProblem::EquilibriumCompositionProblem(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumCompositionProblem::EquilibriumCompositionProblem(const EquilibriumCompositionProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumCompositionProblem::~EquilibriumCompositionProblem()
{}

auto EquilibriumCompositionProblem::operator=(EquilibriumCompositionProblem other) -> EquilibriumCompositionProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumCompositionProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumCompositionProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

auto EquilibriumCompositionProblem::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumCompositionProblem::setTemperature(double value, std::string units) -> void
{
    pimpl->T = units::convert(value, units, "K");
}

auto EquilibriumCompositionProblem::setPressure(double value, std::string units) -> void
{
    pimpl->P = units::convert(value, units, "Pa");
}

auto EquilibriumCompositionProblem::setAqueousComposition(std::string molalities) -> void
{
    pimpl->setAqueousFluid(molalities);
}

auto EquilibriumCompositionProblem::setGaseousComposition(std::string molarfractions) -> void
{
    pimpl->setGaseousFluid(molarfractions);
}

auto EquilibriumCompositionProblem::setAqueousSaturation(double value) -> void
{
    pimpl->saturation_fluid_phases["Aqueous"] = value;
}

auto EquilibriumCompositionProblem::setGaseousSaturation(double value) -> void
{
    pimpl->saturation_fluid_phases["Gaseous"] = value;
}

auto EquilibriumCompositionProblem::setSolidComposition(std::string volumefractions) -> void
{
    pimpl->setSolid(volumefractions);
}

auto EquilibriumCompositionProblem::setPorosity(double value) -> void
{
    pimpl->porosity = value;
}

EquilibriumCompositionProblem::operator EquilibriumInverseProblem()
{
    return pimpl->equilibriumInverseProblem();
}

}  // namespace Reaktoro
