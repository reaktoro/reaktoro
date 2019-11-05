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

#include "ChemicalEditor.hpp"

// C++ includes
#include <set>

// Reaktoro includes
#include <Reaktoro/Common/ElementUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/StringList.hpp>
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
#include <Reaktoro/Thermodynamics/Mixtures/LiquidMixture.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/LiquidPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/LiquidSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

auto collectElementsInCompounds(const std::vector<std::string>& compounds) -> std::vector<std::string>
{
    std::set<std::string> elemset;
    for(const auto& compound : compounds)
        for(auto pair : elements(compound))
            elemset.insert(pair.first);
    return {elemset.begin(), elemset.end()};
}

auto lnActivityConstants(const AqueousPhase& phase) -> ThermoVectorFunction
{
    // The ln activity constants of the aqueous species
    ThermoVector ln_c(phase.numSpecies());

    // The index of solvent water species
    const Index iH2O = phase.indexSpeciesAnyWithError(alternativeWaterNames());

    // Set the ln activity constants of aqueous species to ln(55.508472)
    ln_c = std::log(1.0 / waterMolarMass);

    // Set the ln activity constant of water to zero
    ln_c[iH2O] = 0.0;

    ThermoVectorFunction f = [=](Temperature T, Pressure P) mutable {
        return ln_c;
    };

    return f;
}

auto lnActivityConstants(const FluidPhase& phase) -> ThermoVectorFunction
{
    // The ln activity constants of the generic species
    ThermoVector ln_c(phase.numSpecies());

    ThermoVectorFunction f = [=](Temperature T, Pressure P) mutable {
        ln_c = log(P * 1e-5); // ln(Pbar)
        return ln_c;
    };

    return f;
}

auto lnActivityConstants(const MineralPhase& phase) -> ThermoVectorFunction
{
    // The ln activity constants of the mineral species
    ThermoVector ln_c(phase.numSpecies());

    ThermoVectorFunction f = [=](Temperature T, Pressure P) mutable {
        return ln_c;
    };

    return f;
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

    /// The definition of the liquid phase
    LiquidPhase liquid_phase;

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
        temperatures = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300};

        // The default pressures for the interpolation of the thermodynamic properties (in units of bar)
        pressures = {1, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};

        // Convert the temperatures and pressures to units of kelvin and pascal respectively
        for(auto& x : temperatures)
            x = x + 273.15;
        for(auto& x : pressures)
            x = x * 1.0e+5;
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

    auto initializePhasesWithElements(const std::vector<std::string>& elements) -> void
    {
        aqueous_phase = {};
        gaseous_phase = {};
        liquid_phase = {};
        mineral_phases.clear();

        auto aqueous_species = database.aqueousSpeciesWithElements(elements);
        auto gaseous_species = database.gaseousSpeciesWithElements(elements);
        auto liquid_species = database.liquidSpeciesWithElements(elements);
        auto mineral_species = database.mineralSpeciesWithElements(elements);

        if(aqueous_species.size())
            addPhase(AqueousPhase(AqueousMixture(aqueous_species)));

        if(gaseous_species.size())
            addPhase(GaseousPhase(GaseousMixture(gaseous_species)));

        if(liquid_species.size())
            addPhase(LiquidPhase(LiquidMixture(liquid_species)));

        for(auto mineral : mineral_species)
            addPhase(MineralPhase(MineralMixture(mineral)));
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

    auto addPhase(const LiquidPhase& phase) -> LiquidPhase&
    {
        liquid_phase = phase;
        return liquid_phase;
    }

    auto addPhase(const MineralPhase& phase) -> MineralPhase&
    {
        for(MineralPhase& x : mineral_phases)
            if(x.name() == phase.name())
                return x = phase;
        mineral_phases.push_back(phase);
        return mineral_phases.back();
    }

    auto addAqueousPhaseHelper(const std::vector<AqueousSpecies>& species) -> AqueousPhase&
    {
        AqueousMixture mixture(species);
        return addPhase(AqueousPhase(mixture));
    }

    auto addAqueousPhaseWithSpecies(const std::vector<std::string>& species) -> AqueousPhase&
    {
        Assert(species.size(), "Could not create the AqueousPhase object.",
               "Expecting at least one species name.");
        std::vector<AqueousSpecies> aqueous_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            aqueous_species[i] = database.aqueousSpecies(species[i]);
        return addAqueousPhaseHelper(aqueous_species);
    }

    auto addAqueousPhaseWithElements(const std::vector<std::string>& elements) -> AqueousPhase&
    {
        Assert(elements.size(), "Could not create the AqueousPhase object.",
               "Expecting at least one chemical element or compound name.");
        return addAqueousPhaseHelper(database.aqueousSpeciesWithElements(elements));
    }

    auto addAqueousPhaseWithCompounds(const std::vector<std::string>& compounds) -> AqueousPhase&
    {
        return addAqueousPhaseWithElements(collectElementsInCompounds(compounds));
    }

    auto addGaseousPhaseHelper(const std::vector<GaseousSpecies>& species) -> GaseousPhase&
    {
        GaseousMixture mixture(species);
        gaseous_phase = GaseousPhase(mixture);
        return gaseous_phase;
    }

    auto addGaseousPhaseWithSpecies(const std::vector<std::string>& species) -> GaseousPhase&
    {
        Assert(species.size(), "Could not create the GaseousPhase object that represents a gas.",
               "Expecting at least one species name.");
        std::vector<GaseousSpecies> gaseous_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            gaseous_species[i] = database.gaseousSpecies(species[i]);
        return addGaseousPhaseHelper(gaseous_species);
    }

    auto addGaseousPhaseWithElements(const std::vector<std::string>& elements) -> GaseousPhase&
    {
        Assert(elements.size(), "Could not create the GaseousPhase object that represents a gas.",
               "Expecting at least one chemical element or compound name.");
        return addGaseousPhaseHelper(database.gaseousSpeciesWithElements(elements));
    }

    auto addGaseousPhaseWithCompounds(const std::vector<std::string>& compounds) -> GaseousPhase&
    {
        return addGaseousPhaseWithElements(collectElementsInCompounds(compounds));
    }

    auto addLiquidPhaseHelper(const std::vector<LiquidSpecies>& species) -> LiquidPhase&
    {
        LiquidMixture mixture(species);
        liquid_phase = LiquidPhase(mixture);
        return liquid_phase;
    }

    auto addLiquidPhaseWithSpecies(const std::vector<std::string>& species) -> LiquidPhase&
    {
        Assert(species.size(), "Could not create the LiquidPhase object that represents a liquid.",
               "Expecting at least one species name.");
        std::vector<LiquidSpecies> liquid_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            liquid_species[i] = database.liquidSpecies(species[i]);
        return addLiquidPhaseHelper(liquid_species);
    }

    auto addLiquidPhaseWithElements(const std::vector<std::string>& elements) -> LiquidPhase&
    {
        Assert(elements.size(), "Could not create the FluidPhase object that represents a liquid.",
               "Expecting at least one chemical element or compound name.");
        return addLiquidPhaseHelper(database.liquidSpeciesWithElements(elements));
    }

    auto addLiquidPhaseWithCompounds(const std::vector<std::string>& compounds) -> LiquidPhase&
    {
        return addLiquidPhaseWithElements(collectElementsInCompounds(compounds));
    }

    auto addMineralPhaseHelper(const std::vector<MineralSpecies>& species) -> MineralPhase&
    {
        MineralMixture mixture(species);
        return addPhase(MineralPhase(mixture));
    }

    auto addMineralPhaseWithSpecies(const std::vector<std::string>& species) -> MineralPhase&
    {
        Assert(species.size(), "Could not create the MineralPhase object.",
               "Expecting at least one species name.");
        std::vector<MineralSpecies> mineral_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            mineral_species[i] = database.mineralSpecies(species[i]);
        return addMineralPhaseHelper(mineral_species);
    }

    auto addMineralPhaseWithElements(const std::vector<std::string>& elements) -> MineralPhase&
    {
        Assert(elements.size(), "Could not create the MineralPhase object.",
               "Expecting at least one chemical element or compound name.");
        return addMineralPhaseHelper(database.mineralSpeciesWithElements(elements));
    }

    auto addMineralPhaseWithCompounds(const std::vector<std::string>& compounds) -> MineralPhase&
    {
        return addMineralPhaseWithElements(collectElementsInCompounds(compounds));
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

    auto liquidPhase() const -> const LiquidPhase&
    {
        return liquid_phase;
    }

    auto liquidPhase() -> LiquidPhase&
    {
        return liquid_phase;
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

    template<typename Phase_>
    auto convertPhase(const Phase_& phase) const -> Phase
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
        for(unsigned i = 0; i < nspecies; ++i) {
            const std::string name = phase.species(i).name();

            standard_gibbs_energy_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarGibbsEnergy(T, P, name); };
            standard_enthalpy_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarEnthalpy(T, P, name); };
            standard_volume_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarVolume(T, P, name); };
            standard_heat_capacity_cp_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarHeatCapacityConstP(T, P, name); };
            standard_heat_capacity_cv_fns[i] = [=](double T, double P) { return thermo.standardPartialMolarHeatCapacityConstV(T, P, name); };
        }

        // Create the interpolation functions for thermodynamic properties of the species
        ThermoVectorFunction standard_gibbs_energies_interp = interpolate(temperatures, pressures, standard_gibbs_energy_fns);
        ThermoVectorFunction standard_enthalpies_interp = interpolate(temperatures, pressures, standard_enthalpy_fns);
        ThermoVectorFunction standard_volumes_interp = interpolate(temperatures, pressures, standard_volume_fns);
        ThermoVectorFunction standard_heat_capacities_cp_interp = interpolate(temperatures, pressures, standard_heat_capacity_cp_fns);
        ThermoVectorFunction standard_heat_capacities_cv_interp = interpolate(temperatures, pressures, standard_heat_capacity_cv_fns);
        ThermoVectorFunction ln_activity_constants_func = lnActivityConstants(phase);

        // Define the thermodynamic model function of the species
        PhaseThermoModel thermo_model = [=](PhaseThermoModelResult& res, Temperature T, Pressure P) {
            // Calculate the standard thermodynamic properties of each species
            res.standard_partial_molar_gibbs_energies = standard_gibbs_energies_interp(T, P);
            res.standard_partial_molar_enthalpies = standard_enthalpies_interp(T, P);
            res.standard_partial_molar_volumes = standard_volumes_interp(T, P);
            res.standard_partial_molar_heat_capacities_cp = standard_heat_capacities_cp_interp(T, P);
            res.standard_partial_molar_heat_capacities_cv = standard_heat_capacities_cv_interp(T, P);
            res.ln_activity_constants = ln_activity_constants_func(T, P);

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
        const auto number_of_fluid_phases = 3;
        phases.reserve(number_of_fluid_phases + mineral_phases.size());

        if(aqueous_phase.numSpecies())
            phases.push_back(convertPhase(aqueous_phase));

        if(gaseous_phase.numSpecies())
            phases.push_back(convertPhase(gaseous_phase));

        if(liquid_phase.numSpecies())
            phases.push_back(convertPhase(liquid_phase));

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

auto ChemicalEditor::initializePhasesWithElements(const StringList& elements) -> void
{
    pimpl->initializePhasesWithElements(elements);
}

auto ChemicalEditor::addPhase(const AqueousPhase& phase) -> AqueousPhase&
{
    return pimpl->addPhase(phase);
}

auto ChemicalEditor::addPhase(const GaseousPhase& phase) -> GaseousPhase&
{
    return pimpl->addPhase(phase);
}

auto ChemicalEditor::addPhase(const LiquidPhase& phase) -> LiquidPhase&
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

auto ChemicalEditor::addAqueousPhase(const StringList& species) -> AqueousPhase&
{
    return pimpl->addAqueousPhaseWithSpecies(species);
}

auto ChemicalEditor::addAqueousPhaseWithElements(const StringList& elements) -> AqueousPhase&
{
    return pimpl->addAqueousPhaseWithElements(elements);
}

auto ChemicalEditor::addAqueousPhaseWithElementsOf(const StringList& compounds) -> AqueousPhase&
{
    return pimpl->addAqueousPhaseWithCompounds(compounds);
}

auto ChemicalEditor::addGaseousPhase(const StringList& species) -> GaseousPhase&
{
    return pimpl->addGaseousPhaseWithSpecies(species);
}

auto ChemicalEditor::addGaseousPhaseWithElements(const StringList& elements) -> GaseousPhase&
{
    return pimpl->addGaseousPhaseWithElements(elements);
}

auto ChemicalEditor::addGaseousPhaseWithElementsOf(const StringList& compounds) -> GaseousPhase&
{
    return pimpl->addGaseousPhaseWithCompounds(compounds);
}

auto ChemicalEditor::addLiquidPhase(const StringList& species) -> LiquidPhase&
{
    return pimpl->addLiquidPhaseWithSpecies(species);
}

auto ChemicalEditor::addLiquidPhaseWithElements(const StringList& elements) -> LiquidPhase&
{
    return pimpl->addLiquidPhaseWithElements(elements);
}

auto ChemicalEditor::addLiquidPhaseWithElementsOf(const StringList& compounds) -> LiquidPhase&
{
    return pimpl->addLiquidPhaseWithCompounds(compounds);
}

auto ChemicalEditor::addMineralPhase(const StringList& species) -> MineralPhase&
{
    return pimpl->addMineralPhaseWithSpecies(species);
}

auto ChemicalEditor::addMineralPhaseWithElements(const StringList& elements) -> MineralPhase&
{
    return pimpl->addMineralPhaseWithElements(elements);
}

auto ChemicalEditor::addMineralPhaseWithElementsOf(const StringList& compounds) -> MineralPhase&
{
    return pimpl->addMineralPhaseWithCompounds(compounds);
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

auto ChemicalEditor::liquidPhase() const -> const LiquidPhase&
{
    return pimpl->liquidPhase();
}

auto ChemicalEditor::liquidPhase() -> LiquidPhase&
{
    return pimpl->liquidPhase();
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
