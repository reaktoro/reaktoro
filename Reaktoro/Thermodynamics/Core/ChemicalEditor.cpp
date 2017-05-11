// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "ChemicalEditor.hpp"

// C++ includes
#include <set>

// Reaktoro includes
#include <Reaktoro/Common/ElementUtils.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Core/Thermo.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {
namespace {

auto collectElementsInCompounds(std::vector<std::string> compounds) -> std::vector<std::string>
{
    std::set<std::string> elemset;
    for(auto compound : compounds)
        for(auto pair : elements(compound))
            elemset.insert(pair.first);
    return {elemset.begin(), elemset.end()};
}

} // namespace

struct ChemicalEditor::Impl
{
private:
    /// The database instance
    Database database;

    /// The definition of the aqueous phase
    AqueousPhase aqueous_phase;

    /// The definition of the gaseous phase
    GaseousPhase gaseous_phase;

    /// The definition of the mineral phases
    std::vector<MineralPhase> mineral_phases;

    /// The mineral reactions of the chemical system
    std::vector<MineralReaction> mineral_reactions;

    /// The temperatures for constructing interpolation tables of thermodynamic properties (in units of K).
    std::vector<double> temperatures;

    /// The pressures for constructing interpolation tables of thermodynamic properties (in units of Pa).
    std::vector<double> pressures;

public:
    Impl()
    : Impl(Database("supcrt98"))
    {
    }

    explicit Impl(const Database& database)
    : database(database)
    {
        // The default temperatures for the interpolation of the thermodynamic properties (in units of celsius)
        temperatures = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300 };

        // The default pressures for the interpolation of the thermodynamic properties (in units of bar)
        pressures = { 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600 };

        // Convert the temperatures and pressures to units of kelvin and pascal respectively
        for(auto& x : temperatures) x = x + 273.15;
        for(auto& x : pressures)    x = x * 1.0e+5;
    }

    auto setTemperatures(std::vector<double> values, std::string units) -> void
    {
        temperatures = values;
        for(auto& x : temperatures)
            x = units::convert(x, units, "kelvin");
    }

    auto setPressures(std::vector<double> values, std::string units) -> void
    {
        pressures = values;
        for(auto& x : pressures)
            x = units::convert(x, units, "pascal");
    }

    auto addPhase(const AqueousPhase& phase) -> AqueousPhase&
    {
        aqueous_phase = phase;
        aqueous_phase.setInterpolationPoints(temperatures, pressures);
        return aqueous_phase;
    }

    auto addPhase(const GaseousPhase& phase) -> GaseousPhase&
    {
        gaseous_phase = phase;
        return gaseous_phase;
    }

    auto addPhase(const MineralPhase& phase) -> MineralPhase&
    {
        for(MineralPhase& x : mineral_phases)
            if(x.name() == phase.name())
                return x = phase;
        mineral_phases.push_back(phase);
        return mineral_phases.back();
    }

    auto addAqueousPhaseHelper(std::vector<AqueousSpecies> species) -> AqueousPhase&
    {
        AqueousMixture mixture(species);
        return addPhase(AqueousPhase(mixture));
    }

    auto addAqueousPhaseWithSpecies(std::vector<std::string> species) -> AqueousPhase&
    {
        Assert(species.size(), "Could not create the AqueousPhase object.",
            "Expecting at least one species name.");
        std::vector<AqueousSpecies> aqueous_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            aqueous_species[i] = database.aqueousSpecies(species[i]);
        return addAqueousPhaseHelper(aqueous_species);
    }

    auto addAqueousPhaseWithElements(std::vector<std::string> elements) -> AqueousPhase&
    {
        Assert(elements.size(), "Could not create the AqueousPhase object.",
            "Expecting at least one chemical element or compound name.");
        return addAqueousPhaseHelper(database.aqueousSpeciesWithElements(elements));
    }

    auto addAqueousPhaseWithCompounds(std::vector<std::string> compounds) -> AqueousPhase&
    {
        return addAqueousPhaseWithElements(collectElementsInCompounds(compounds));
    }

    auto addAqueousPhase(std::vector<std::string> species) -> AqueousPhase&
    {
        return addAqueousPhaseWithSpecies(species);
    }

    auto addAqueousPhase(std::string compounds) -> AqueousPhase&
    {
        auto words = split(compounds, " ");
        return addAqueousPhaseWithCompounds(words);
    }

    auto addGaseousPhaseHelper(const std::vector<GaseousSpecies>& species) -> GaseousPhase&
    {
        GaseousMixture mixture(species);
        gaseous_phase = GaseousPhase(mixture);
        return gaseous_phase;
    }

    auto addGaseousPhaseWithSpecies(std::vector<std::string> species) -> GaseousPhase&
    {
        Assert(species.size(), "Could not create the GaseousPhase object.",
            "Expecting at least one species name.");
        std::vector<GaseousSpecies> gaseous_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            gaseous_species[i] = database.gaseousSpecies(species[i]);
        return addGaseousPhaseHelper(gaseous_species);
    }

    auto addGaseousPhaseWithElements(std::vector<std::string> elements) -> GaseousPhase&
    {
        Assert(elements.size(), "Could not create the GaseousPhase object.",
            "Expecting at least one chemical element or compound name.");
        return addGaseousPhaseHelper(database.gaseousSpeciesWithElements(elements));
    }

    auto addGaseousPhaseWithCompounds(std::vector<std::string> compounds) -> GaseousPhase&
    {
        return addGaseousPhaseWithElements(collectElementsInCompounds(compounds));
    }

    auto addGaseousPhase(std::vector<std::string> species) -> GaseousPhase&
    {
        return addGaseousPhaseWithSpecies(species);
    }

    auto addGaseousPhase(std::string compounds) -> GaseousPhase&
    {
        auto words = split(compounds, " ");
        return addGaseousPhaseWithCompounds(words);
    }

    auto addMineralPhaseHelper(const std::vector<MineralSpecies>& species) -> MineralPhase&
    {
        MineralMixture mixture(species);
        return addPhase(MineralPhase(mixture));
    }

    auto addMineralPhaseWithSpecies(std::vector<std::string> species) -> MineralPhase&
    {
        Assert(species.size(), "Could not create the MineralPhase object.",
            "Expecting at least one species name.");
        std::vector<MineralSpecies> mineral_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            mineral_species[i] = database.mineralSpecies(species[i]);
        return addMineralPhaseHelper(mineral_species);
    }

    auto addMineralPhaseWithElements(std::vector<std::string> elements) -> MineralPhase&
    {
        Assert(elements.size(), "Could not create the MineralPhase object.",
            "Expecting at least one chemical element or compound name.");
        return addMineralPhaseHelper(database.mineralSpeciesWithElements(elements));
    }

    auto addMineralPhaseWithCompounds(std::vector<std::string> compounds) -> MineralPhase&
    {
        return addMineralPhaseWithElements(collectElementsInCompounds(compounds));
    }

    auto addMineralPhase(std::vector<std::string> species) -> MineralPhase&
    {
        return addMineralPhaseWithSpecies(species);
    }

    auto addMineralPhase(std::string compounds) -> MineralPhase&
    {
        auto words = split(compounds, " ");
        if(words.size() == 1 && database.containsMineralSpecies(words[0]))
            return addMineralPhaseWithSpecies({words[0]});
        return addMineralPhaseWithCompounds(words);
    }

    auto addReaction(const MineralReaction& reaction) -> MineralReaction&
    {
        mineral_reactions.push_back(reaction);
        return mineral_reactions.back();
    }

    auto addMineralReaction(const MineralReaction& reaction) -> MineralReaction&
    {
        mineral_reactions.push_back(reaction);
        return mineral_reactions.back();
    }

    auto addMineralReaction(std::string mineral) -> MineralReaction&
    {
        return addMineralReaction(MineralReaction(mineral));
    }

    auto aqueousPhase() const -> const AqueousPhase&
    {
        return aqueous_phase;
    }

    auto aqueousPhase() -> AqueousPhase&
    {
        return aqueous_phase;
    }

    auto gaseousPhase() const -> const GaseousPhase&
    {
        return gaseous_phase;
    }

    auto gaseousPhase() -> GaseousPhase&
    {
        return gaseous_phase;
    }

    auto mineralPhases() const -> const std::vector<MineralPhase>&
    {
        return mineral_phases;
    }

    auto mineralPhases() -> std::vector<MineralPhase>&
    {
        return mineral_phases;
    }

    template<typename SpeciesType>
    auto convertSpecies(const SpeciesType& species) const -> Species
    {
        // Create the Species instance
        Species converted;
        converted.setName(species.name());
        converted.setFormula(species.formula());
        converted.setElements(species.elements());

        return converted;
    }

    template<typename PhaseType>
    auto convertPhase(const PhaseType& phase) const -> Phase
    {
        // The number of species in the phase
        const unsigned nspecies = phase.numSpecies();

        // Define the lambda functions for the calculation of the essential thermodynamic properties
        Thermo thermo(database);

        std::vector<ThermoScalarFunction> standard_gibbs_energy_fns(nspecies);
        std::vector<ThermoScalarFunction> standard_enthalpy_fns(nspecies);
        std::vector<ThermoScalarFunction> standard_volume_fns(nspecies);
        std::vector<ThermoScalarFunction> standard_heat_capacity_cp_fns(nspecies);
        std::vector<ThermoScalarFunction> standard_heat_capacity_cv_fns(nspecies);

        // Create the ThermoScalarFunction instances for each thermodynamic properties of each species
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const std::string name = phase.species(i).name();

            standard_gibbs_energy_fns[i]     = [=](double T, double P) { return thermo.standardPartialMolarGibbsEnergy(T, P, name); };
            standard_enthalpy_fns[i]         = [=](double T, double P) { return thermo.standardPartialMolarEnthalpy(T, P, name); };
            standard_volume_fns[i]           = [=](double T, double P) { return thermo.standardPartialMolarVolume(T, P, name); };
            standard_heat_capacity_cp_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarHeatCapacityConstP(T, P, name); };
            standard_heat_capacity_cv_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarHeatCapacityConstV(T, P, name); };
        }

        // Create the interpolation functions for thermodynamic properties of the species
        ThermoVectorFunction standard_gibbs_energies_interp     = interpolate(temperatures, pressures, standard_gibbs_energy_fns);
        ThermoVectorFunction standard_enthalpies_interp         = interpolate(temperatures, pressures, standard_enthalpy_fns);
        ThermoVectorFunction standard_volumes_interp            = interpolate(temperatures, pressures, standard_volume_fns);
        ThermoVectorFunction standard_heat_capacities_cp_interp = interpolate(temperatures, pressures, standard_heat_capacity_cp_fns);
        ThermoVectorFunction standard_heat_capacities_cv_interp = interpolate(temperatures, pressures, standard_heat_capacity_cv_fns);

        // Define the thermodynamic model function of the species
        PhaseThermoModel thermo_model = [=](double T, double P)
        {
            // Calculate the standard thermodynamic properties of each species
            PhaseThermoModelResult res;
            res.standard_partial_molar_gibbs_energies     = standard_gibbs_energies_interp(T, P);
            res.standard_partial_molar_enthalpies         = standard_enthalpies_interp(T, P);
            res.standard_partial_molar_volumes            = standard_volumes_interp(T, P);
            res.standard_partial_molar_heat_capacities_cp = standard_heat_capacities_cp_interp(T, P);
            res.standard_partial_molar_heat_capacities_cv = standard_heat_capacities_cv_interp(T, P);

            return res;
        };

        // Create the Phase instance
        Phase converted = phase;
        converted.setThermoModel(thermo_model);

        return converted;
    }

    auto createChemicalSystem() const -> ChemicalSystem
    {
        std::vector<Phase> phases;
        phases.reserve(2 + mineral_phases.size());

        if(aqueous_phase.numSpecies())
            phases.push_back(convertPhase(aqueous_phase));

        if(gaseous_phase.numSpecies())
            phases.push_back(convertPhase(gaseous_phase));

        for(const MineralPhase& mineral_phase : mineral_phases)
            phases.push_back(convertPhase(mineral_phase));

        return ChemicalSystem(phases);
    }

    auto createReactionSystem() const -> ReactionSystem
    {
        ChemicalSystem system = createChemicalSystem();

        std::vector<Reaction> reactions;
        for(const MineralReaction& rxn : mineral_reactions)
            reactions.push_back(createReaction(rxn, system));

        return ReactionSystem(system, reactions);
    }
};

ChemicalEditor::ChemicalEditor()
: pimpl(new Impl())
{}

ChemicalEditor::ChemicalEditor(const Database& database)
: pimpl(new Impl(database))
{}

ChemicalEditor::ChemicalEditor(const ChemicalEditor& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalEditor::~ChemicalEditor()
{}

auto ChemicalEditor::operator=(const ChemicalEditor& other) -> ChemicalEditor&
{
    pimpl.reset(new Impl(*other.pimpl));
    return *this;
}

auto ChemicalEditor::setTemperatures(std::vector<double> values, std::string units) -> void
{
    pimpl->setTemperatures(values, units);
}

auto ChemicalEditor::setPressures(std::vector<double> values, std::string units) -> void
{
    pimpl->setPressures(values, units);
}

auto ChemicalEditor::addPhase(const AqueousPhase& phase) -> AqueousPhase&
{
    return pimpl->addPhase(phase);
}

auto ChemicalEditor::addPhase(const GaseousPhase& phase) -> GaseousPhase&
{
    return pimpl->addPhase(phase);
}

auto ChemicalEditor::addPhase(const MineralPhase& phase) -> MineralPhase&
{
    return pimpl->addPhase(phase);
}

auto ChemicalEditor::addReaction(const MineralReaction& reaction) -> MineralReaction&
{
    return pimpl->addReaction(reaction);
}

auto ChemicalEditor::addAqueousPhase(std::vector<std::string> species) -> AqueousPhase&
{
    return pimpl->addAqueousPhase(species);
}

auto ChemicalEditor::addAqueousPhase(std::initializer_list<std::string> species) -> AqueousPhase&
{
    return pimpl->addAqueousPhase(std::vector<std::string>(species.begin(), species.end()));
}

auto ChemicalEditor::addAqueousPhase(std::string compounds) -> AqueousPhase&
{
    return pimpl->addAqueousPhase(compounds);
}

auto ChemicalEditor::addGaseousPhase(std::vector<std::string> species) -> GaseousPhase&
{
    return pimpl->addGaseousPhase(species);
}

auto ChemicalEditor::addGaseousPhase(std::initializer_list<std::string> species) -> GaseousPhase&
{
    return pimpl->addGaseousPhase(std::vector<std::string>(species.begin(), species.end()));
}

auto ChemicalEditor::addGaseousPhase(std::string compounds) -> GaseousPhase&
{
    return pimpl->addGaseousPhase(compounds);
}

auto ChemicalEditor::addMineralPhase(std::vector<std::string> species) -> MineralPhase&
{
    return pimpl->addMineralPhase(species);
}

auto ChemicalEditor::addMineralPhase(std::initializer_list<std::string> species) -> MineralPhase&
{
    return pimpl->addMineralPhase(std::vector<std::string>(species.begin(), species.end()));
}

auto ChemicalEditor::addMineralPhase(std::string compounds) -> MineralPhase&
{
    return pimpl->addMineralPhase(compounds);
}

auto ChemicalEditor::addMineralReaction(const MineralReaction& reaction) -> MineralReaction&
{
    return pimpl->addMineralReaction(reaction);
}

auto ChemicalEditor::addMineralReaction(std::string mineral) -> MineralReaction&
{
    return pimpl->addMineralReaction(mineral);
}

auto ChemicalEditor::aqueousPhase() const -> const AqueousPhase&
{
    return pimpl->aqueousPhase();
}

auto ChemicalEditor::aqueousPhase() -> AqueousPhase&
{
    return pimpl->aqueousPhase();
}

auto ChemicalEditor::gaseousPhase() const -> const GaseousPhase&
{
    return pimpl->gaseousPhase();
}

auto ChemicalEditor::gaseousPhase() -> GaseousPhase&
{
    return pimpl->gaseousPhase();
}

auto ChemicalEditor::mineralPhases() const -> const std::vector<MineralPhase>&
{
    return pimpl->mineralPhases();
}

auto ChemicalEditor::mineralPhases() -> std::vector<MineralPhase>&
{
    return pimpl->mineralPhases();
}

auto ChemicalEditor::createChemicalSystem() const -> ChemicalSystem
{
    return pimpl->createChemicalSystem();
}

auto ChemicalEditor::createReactionSystem() const -> ReactionSystem
{
    return pimpl->createReactionSystem();
}

ChemicalEditor::operator ChemicalSystem() const
{
    return createChemicalSystem();
}

ChemicalEditor::operator ReactionSystem() const
{
    return createReactionSystem();
}

} // namespace Reaktoro
