// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// Phreeqc includes
#define Phreeqc PHREEQC
#define protected public
#include <phreeqc/Phreeqc.h>
#include <phreeqc/GasPhase.h>
#undef Phreeqc

const unsigned ncols = 2;

void markdown(PHREEQC& obj);
void cppmap(PHREEQC& obj);

int main(int argc, char **argv)
{
	Phreeqc phreeqc("databases/phreeqc/phreeqc.dat");
	PHREEQC& obj = phreeqc.phreeqc();

	markdown(obj);
	cppmap(obj);
}

void markdown(PHREEQC& obj)
{
	std::vector<std::string> lines;
	for(int i = 0; i < obj.count_s; i++)
	{
		// Skip species with both a and b equal to 0.0
		if(obj.s[i]->dha == 0.0 && obj.s[i]->dhb == 0.0)
			continue;

		// Skip neutral species with default b = 0.1 (which are many)
		if(obj.s[i]->dha == 0.0 && obj.s[i]->dhb == 0.1)
			continue;

		// Skip exchange species
		if(obj.s[i]->exch_gflag)
			continue;

		std::string name = (obj.s[i]->z) ?
			conventionalChargedSpeciesName(obj.s[i]->name) :
			conventionalNeutralSpeciesName(obj.s[i]->name);
		name = "`" + name + "`";

		std::stringstream ss;
		ss << "| " << std::left << std::setw(15) << name;
		ss << "| " << std::left << std::setw(8) << obj.s[i]->dha;
		ss << "| " << std::left << std::setw(6) << obj.s[i]->dhb;
		lines.push_back(ss.str());
	}

	for(unsigned i = 0; i < ncols; ++i)
	{
		std::cout << "| " << std::left << std::setw(15) << "Ion";
		std::cout << "| " << std::left << std::setw(8) << "*å* (Å)";
		std::cout << "| " << std::left << std::setw(6) << " *b*";
	}

	std::cout << std::endl;

	for(unsigned i = 0; i < ncols; ++i)
	{
		std::cout << "| " << std::left << std::setw(15) << "-";
		std::cout << "| " << std::left << std::setw(8) << "-";
		std::cout << "| " << std::left << std::setw(6) << "-";
	}

	std::cout << std::endl;

	for(unsigned j = 0; j < lines.size(); j += ncols)
	{
		for(unsigned i = 0; i < ncols; ++i)
		{
			if(j + i < lines.size())
				std::cout << lines[j + i];
			else
			{
                std::cout << "| " << std::left << std::setw(15) << "---";
                std::cout << "| " << std::left << std::setw(8) << "---";
                std::cout << "| " << std::left << std::setw(6) << "---";
			}
		}

		std::cout << std::endl;
	}

	std::cout << std::endl;
}

void cppmap(PHREEQC& obj)
{
	std::vector<std::string> alines, blines;
	for(int i = 0; i < obj.count_s; i++)
	{
		// Skip species with both a and b equal to 0.0
		if(obj.s[i]->dha == 0.0 && obj.s[i]->dhb == 0.0)
			continue;

		// Skip neutral species with default b = 0.1 (which are many)
		if(obj.s[i]->dha == 0.0 && obj.s[i]->dhb == 0.1)
			continue;

		// Skip exchange species
		if(obj.s[i]->exch_gflag)
			continue;

		std::string name = (obj.s[i]->z) ?
			conventionalChargedSpeciesName(obj.s[i]->name) :
			conventionalNeutralSpeciesName(obj.s[i]->name);
		name = "\"" + name + "\"";

		std::stringstream as, bs;
		as << "{" << name << ", " << obj.s[i]->dha << "}";
		bs << "{" << name << ", " << obj.s[i]->dhb << "}";

		if(obj.s[i]->dha)
			alines.push_back(as.str());
		if(obj.s[i]->dhb)
			blines.push_back(bs.str());
	}

	std::cout << "const std::map<std::string, double> aion_phreeqc = {";
	for(unsigned i = 0; i < alines.size(); ++i)
		std::cout << (i == 0 ? "" : ", ") << alines[i];
	std::cout << "};" << std::endl;

	std::cout << std::endl;

	std::cout << "const std::map<std::string, double> bion_phreeqc = {";
	for(unsigned i = 0; i < blines.size(); ++i)
		std::cout << (i == 0 ? "" : ", ") << blines[i];
	std::cout << "};" << std::endl;

	std::cout << std::endl;
}
