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

template<typename NamedValues>
auto index(const std::string& name, const NamedValues& values) -> Index
{
    Index idx = 0;
    for(const auto& value : values)
        if(value.name() == name) return idx; else ++idx;
    return idx;
}

template<typename NamedValue, typename NamedValues>
auto index(const NamedValue& value, const NamedValues& values) -> Index
{
    return index(value.name(), values);
}

template<typename Names, typename NamedValues>
auto indexAny(const Names& names, const NamedValues& values) -> Index
{
    for(auto& name : names)
    {
        const Index i = index(name, values);
        if(i < values.size())
            return i;
    }
    return values.size();
}

template<typename NamedValues>
auto indices(const std::vector<std::string>& names, const NamedValues& values) -> Indices
{
    Indices idxs;
    idxs.reserve(names.size());
    for(const auto& name : names)
        idxs.push_back(index(name, values));
    return idxs;
}

template<typename NamedValues>
auto indices(const NamedValues& subvalues, const NamedValues& values) -> Indices
{
    Indices idxs;
    idxs.reserve(subvalues.size());
    for(const auto& value : subvalues)
        idxs.push_back(index(value, values));
    return idxs;
}

template<typename NamedValues>
auto contained(const std::string& name, const NamedValues& values) -> bool
{
    return index(name, values) < values.size();
}

} // namespace Reaktoro
