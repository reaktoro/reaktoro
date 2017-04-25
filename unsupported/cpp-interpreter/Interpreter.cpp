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

#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// Reaktoro includes
#include <unsupported/cpp-interpreter/Interpreter.hpp>
#include <unsupported/cpp-interpreter/Processors.hpp>
#include <unsupported/cpp-interpreter/Utils.hpp>

namespace Reaktoro {

struct Interpreter::Impl
{
    /// The state of the interpreter
    InterpreterState istate;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Initialize essential members of the interpreter state.
    /// This method is used to initialize a ChemicalEditor instance by firstly
    /// identifying all compound and species names in the input file. After this,
    /// an aqueous phase is created with all possible species found in the database
    /// that contains the elements composing the list of found compounds. A gaseous
    /// phase is created only if names of gaseous species present in the database are
    /// found in the input script. Pure mineral phases are created for each mineral name
    /// found in the input script file that is also present in the database.
    auto init(const Node& root) -> void
    {
        // Initialize the database
        istate.database = Database("supcrt98");

        // Initialize the list of compounds found in the script file
        istate.compounds = collectCompounds(root);

        // Initialize the list of elements that compose the compounds
        istate.elements = identifyElements(istate.compounds, istate.database);

        // Initialize the automatic chemical system
        initChemicalEditor(root);
    }

    // Initialize the default state of chemical editor based on the compounds found in the script
    auto initChemicalEditor(const Node& root) -> void
    {
        // Auxiliary references
        const auto& database = istate.database;
        const auto& elements = istate.elements;

        // Collect all aqueous species that can be formed out of the elements found in the script
        auto aqueous_species = names(database.aqueousSpeciesWithElements(elements));

        // Determine if there are gaseous and mineral species among the found compounds in the script
        auto mineral_species = filterMineralSpecies(istate.compounds, istate.database);
        auto gaseous_species = filterGaseousSpecies(istate.compounds, istate.database);

        // Check if there is any Speciation, in which case all possible gases and minerals are considered
        if(hasSpeciation(root))
        {
            gaseous_species = names(database.gaseousSpeciesWithElements(elements));
            mineral_species = names(database.mineralSpeciesWithElements(elements));
        }

        // Add the aqueous phase using the compounds as the initializer
        istate.editor.addAqueousPhaseWithSpecies(aqueous_species);

        // Add a gaseous phase if there are gaseous species
        if(gaseous_species.size())
            istate.editor.addGaseousPhaseWithSpecies(gaseous_species);

        // Add a mineral phase for each mineral species
        for(auto x : mineral_species)
            istate.editor.addMineralPhaseWithSpecies({x});
    }

    /// Execute a Reaktoro input script as string.
    auto execute(std::string str) -> void
    {
        // Preprocess the input script so that it conforms with YAML rules.
        str = preprocess(str);

        // The root node of the yaml script
        Node root = YAML::Load(str);

        // Initialize some essential members of the interpreter state
        init(root);

        // Auxiliary type alias to a yaml node process function
        using ftype = std::function<void(InterpreterState&, const Node&)>;

        // The map of process functions (from keyword to respective process function)
        std::map<std::string, ftype> fmap = {
            {"database"           , processDatabaseNode},
            {"aqueousphase"       , processAqueousPhaseNode},
            {"gaseousphase"       , processGaseousPhaseNode},
            {"mineralphase"       , processMineralPhaseNode},
            {"minerals"           , processMineralsNode},
            {"mineralreaction"    , processMineralReactionNode},
            {"equilibrium"        , processEquilibriumNode},
            {"equilibriumproblem" , processEquilibriumNode},
            {"equilibriumpath"    , processEquilibriumPathNode},
            {"kinetics"           , processKineticPathNode},
            {"kineticpath"        , processKineticPathNode},
            {"phreeqc"            , processPhreeqcNode},
        };

        // Keyword triggers that start the initialization of chemical system
        std::set<std::string> triggers = {
            "equilibrium",
            "equilibriumproblem",
            "speciation",
            "speciationproblem",
        };

        // For every child node in the root node...
        for(auto child : root)
        {
            // The keyword of the current yaml node
            std::string key = lowercase(keyword(child));

            // Check if current key is a chemical system initializer
            if(triggers.count(key) && istate.system.numSpecies() == 0)
                istate.system = istate.editor;

            // Find an entry in the process function map with that key
            auto it = fmap.find(key);

            // Assert that this entry was found
            Assert(it != fmap.end(), "Could not parse `" + child + "`.",
                "Expecting a valid keyword. Did you misspelled it?");

            // Process the current child node and update the interpreter state
            it->second(istate, child);
        }
    }
};

Interpreter::Interpreter()
: pimpl(new Impl())
{}

Interpreter::Interpreter(const Interpreter& other)
: pimpl(new Impl(*other.pimpl))
{}

Interpreter::~Interpreter()
{}

auto Interpreter::operator=(Interpreter other) -> Interpreter&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Interpreter::execute(std::string str) -> void
{
    pimpl->execute(str);
}

auto Interpreter::execute(std::istream& stream) -> void
{
    std::stringstream ss; ss << stream.rdbuf();
    execute(ss.str());
}

} // namespace Reaktoro
