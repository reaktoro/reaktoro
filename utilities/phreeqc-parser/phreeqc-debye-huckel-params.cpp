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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

// Phreeqc includes
#define Phreeqc PHREEQC
#define protected public
#include <phreeqc/Phreeqc.h>
#include <phreeqc/GasPhase.h>
#undef Phreeqc

int main(int argc, char **argv)
{
	Phreeqc phreeqc("databases/phreeqc/phreeqc.dat");

	PHREEQC& obj = phreeqc.phreeqc();

	std::cout << std::left << std::setw(20) << "Species";
	std::cout << std::left << std::setw(20) << "DebyeHuckel(a)";
	std::cout << std::left << std::setw(20) << "DebyeHuckel(b)";
	std::cout << std::endl;
	for(int i = 0; i < obj.count_s; i++)
	{
		if(obj.s[i]->dha == 0.0 && obj.s[i]->dhb == 0.0)
			continue;

		if(obj.s[i]->dha == 0.0 && obj.s[i]->dhb == 0.1)
			continue;

		std::cout << std::left << std::setw(20) << obj.s[i]->name;
		std::cout << std::left << std::setw(20) << obj.s[i]->dha;
		std::cout << std::left << std::setw(20) << obj.s[i]->dhb;
		std::cout << std::endl;
	}
}
