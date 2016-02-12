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

#include "Interpreter.hpp"

// ReaktoroInterpreter includes
#include <ReaktoroInterpreterCpp/ParserUtils.hpp>
#include <ReaktoroInterpreterCpp/ProcessUtils.hpp>
#include <ReaktoroInterpreterCpp/InterpreterState.hpp>

namespace Reaktoro {

struct Interpreter::Impl
{
    /// The state of the interpreter
    InterpreterState istate;

    /// Construct a default Impl instance
    Impl()
    {
        // Initialize the database with a built-in database
        istate.database = Database("supcrt98");
    }

    /// Execute a Reaktoro input script as string.
    auto execute(std::string str) -> void
    {
        // Preprocess the input script so that it conforms with YAML rules.
        str = preprocess(str);

        // The root node of the yaml script
        Node root = YAML::Load(str);

        // Auxiliary type alias to a yaml node process function
        using ftype = std::function<void(InterpreterState&, const Node&)>;

        // The map of process functions (from keyword to respective process function)
        std::map<std::string, ftype> fmap = {
            {"mineralreaction"    , processMineralReactionNode},
            {"equilibrium"        , processEquilibriumNode},
            {"equilibriumproblem" , processEquilibriumNode},
            {"kinetics"           , processKineticsNode},
            {"kineticpath"        , processKineticsNode},
        };

        // For every child node in the root node...
        for(auto child : root)
        {
            // The keyword of the current yaml node
            std::string key = lowercase(keyword(child));

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
