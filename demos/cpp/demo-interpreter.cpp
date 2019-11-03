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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

/// The json input for the interpreter as a string
const std::string input = R"(
	{
	  "system": {
		"database": "supcrt98.xml",
		"elements": ["H", "O", "C", "Na", "Cl"]
	  },

	  "calculations": [
		{
		  "equilibrium": {
			"stateReference": "state-h2o-nacl-co2",
			"temperature": { "value": 100.0, "units": "celsius" },
			"pressure": { "value": 300.0, "units": "bar" },
			"substances": [
			  { "substance": "H2O", "quantity": 1.0, "units": "kg" },
			  { "substance": "CO2", "quantity": 2.0, "units": "mol" },
			  { "substance": "NaCl", "quantity": 1.0, "units": "mol" }
			]
		  }
		}
	  ]
	}
)";

int main()
{
	// Create the interpreter to execute an input script in json format
    Interpreter interpreter;

	// Execute the json input script (as a string)
    interpreter.executeJsonString(input);

    // Output to a file the saved chemical state computed during the execution
    interpreter.state("state-h2o-nacl-co2").output("demo-interpreter.txt");
}
