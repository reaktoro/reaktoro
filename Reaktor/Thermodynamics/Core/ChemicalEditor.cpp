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

#include "ChemicalEditor.hpp"

// Reaktor includes
#include <Reaktor/Common/StringUtils.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Thermodynamics/Core/Database.hpp>
#include <Reaktor/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktor/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktor/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktor/Thermodynamics/Reactions/MineralReaction.hpp>

namespace Reaktor {
namespace {

auto convertSpecies(const AqueousSpecies& species, const Database& database) -> Species
{
//    species.setStandardGibbsEnergy();
    return species;
}

} // namespace
struct ChemicalEditor::Impl
{
private:
    /// The database instance
    Database database;

    /// The current state of the aqueous phase instance
    AqueousPhase aqueous_phase;

    /// The current state of the gaseous phase instance
    GaseousPhase gaseous_phase;

    /// The current state of the mineral phase instance
    std::vector<MineralPhase> mineral_phases;

    /// The mineral reactions of the chemical system
    std::vector<MineralReaction> mineral_reactions;

public:
    explicit Impl(const Database& database)
    : database(database)
    {}

    auto addPhase(const AqueousPhase& phase) -> AqueousPhase&
    {
        aqueous_phase = phase;
        return aqueous_phase;
    }

    auto addPhase(const GaseousPhase& phase) -> GaseousPhase&
    {
        gaseous_phase = phase;
        return gaseous_phase;
    }

    auto addPhase(const MineralPhase& phase) -> MineralPhase&
    {
        mineral_phases.push_back(phase);
        return mineral_phases.back();
    }

    auto addAqueousPhase(const std::vector<std::string>& species) -> AqueousPhase&
    {
        // Collect the aqueous species instances from the database
        std::vector<AqueousSpecies> aqueous_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            aqueous_species[i] = database.aqueousSpecies(species[i]);

        // Create the aqueous phase
        aqueous_phase = AqueousPhase(aqueous_species);

        // Set the default activity models
        aqueous_phase.setActivityModelHKFWater();
        aqueous_phase.setActivityModelHKFChargedSpecies();
        aqueous_phase.setActivityModelDuanSunCO2();

        return aqueous_phase;
    }

    auto addAqueousPhase(const std::string& species) -> AqueousPhase&
    {
        std::vector<std::string> words = split(species, " ");

        return addAqueousPhase(words);
    }

    auto addGaseousPhase(const std::vector<std::string>& species) -> GaseousPhase&
    {
        std::vector<GaseousSpecies> gaseous_species(species.size());

        for(unsigned i = 0; i < species.size(); ++i)
            gaseous_species[i] = database.gaseousSpecies(species[i]);

        gaseous_phase = GaseousPhase(gaseous_species);

        gaseous_phase.setActivityModelDuanSunCO2();
        gaseous_phase.setActivityModelIdeal("H2O(g)");

        return gaseous_phase;
    }

    auto addGaseousPhase(const std::string& species) -> GaseousPhase&
    {
        std::vector<std::string> words = split(species, " ");

        return addGaseousPhase(words);
    }

    auto addMineralPhase(const std::vector<std::string>& species) -> MineralPhase&
    {
        std::vector<MineralSpecies> mineral_species(species.size());
        for(unsigned i = 0; i < species.size(); ++i)
            mineral_species[i] = database.mineralSpecies(species[i]);

        mineral_phases.push_back(MineralPhase(mineral_species));

        return mineral_phases.back();
    }

    auto addMineralPhase(const std::string& species) -> MineralPhase&
    {
        std::vector<std::string> words = split(species, " ");

        return addMineralPhase(words);
    }

    auto addReaction(const MineralReaction& reaction) -> MineralReaction&
    {
        mineral_reactions.push_back(reaction);
        return mineral_reactions.back();
    }

    auto addMineralReaction(const std::string& mineral) -> MineralReaction&
    {
        mineral_reactions.push_back(MineralReaction(mineral));
        return mineral_reactions.back();
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

    template<typename PhaseType>
    auto convert(const PhaseType& phase) const -> Phase
    {
        Phase converted = createPhase(phase);
        return converted;
    }

    auto createChemicalSystem() const -> ChemicalSystem
    {
//        std::vector<Phase> phases;
//        phases.reserve(2 + mineral_phases.size());
//
//        if(aqueous_phase.numSpecies())
//            phases.push_back(convert(aqueous_phase));
//
//        if(gaseous_phase.numSpecies())
//            phases.push_back(convert(gaseous_phase));
//
//        for(const MineralPhase& mineral_phase : mineral_phases)
//            phases.push_back(convert(mineral_phase));
//
//        return ChemicalSystem(phases);
    }
};

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

auto ChemicalEditor::addAqueousPhase(const std::vector<std::string>& species) -> AqueousPhase&
{
    return pimpl->addAqueousPhase(species);
}

auto ChemicalEditor::addAqueousPhase(const std::string& species) -> AqueousPhase&
{
    return pimpl->addAqueousPhase(species);
}

auto ChemicalEditor::addGaseousPhase(const std::vector<std::string>& species) -> GaseousPhase&
{
    return pimpl->addGaseousPhase(species);
}

auto ChemicalEditor::addGaseousPhase(const std::string& species) -> GaseousPhase&
{
    return pimpl->addGaseousPhase(species);
}

auto ChemicalEditor::addMineralPhase(const std::vector<std::string>& species) -> MineralPhase&
{
    return pimpl->addMineralPhase(species);
}

auto ChemicalEditor::addMineralPhase(const std::string& species) -> MineralPhase&
{
    return pimpl->addMineralPhase(species);
}

auto ChemicalEditor::addMineralReaction(const std::string& mineral) -> MineralReaction&
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

ChemicalEditor::operator ChemicalSystem() const
{
    return createChemicalSystem();
}

} // namespace Reaktor
