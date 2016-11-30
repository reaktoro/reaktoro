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

// C++ includes
namespace Reaktoro {

template<typename NamedValues>
auto names(const NamedValues& values) -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(values.size());
    for(const auto& entry : values)
        names.push_back(entry.name());
    return names;
}

template<typename SpeciesValues>
auto charges(const SpeciesValues& species) -> Vector
{
    Vector charges(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        charges[i] = species[i].charge();
    return charges;
}

template<typename SpeciesValues>
auto molarMasses(const SpeciesValues& species) -> Vector
{
    Vector molar_masses(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        molar_masses[i] = species[i].molarMass();
    return molar_masses;
}

template<typename Derived>
auto molarFractions(const Eigen::MatrixBase<Derived>& n) -> ChemicalVector
{
	const auto nc = composition(n);
    const unsigned nspecies = n.size();
    if(nspecies == 1)
        return ChemicalVector(nspecies, nspecies, 1.0);
    const ChemicalScalar nt = sum(nc);
    return (nt != 0.0) ? nc/nt : ChemicalVector(nspecies);
}

} // namespace Reaktoro
