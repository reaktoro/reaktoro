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

#include "Phase.hpp"

// C++ includes
#include <algorithm>

// Reaktor includes
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {

class Phase::Impl
{
private:
    /// The name of the phase
    std::string name$;

    /// The chemical species that compose the phase
    std::vector<Species> species$;

    /// The concentration function of the phase
    ConcentrationFn concentration$;

    /// The activity function of the phase
    ActivityFn activity$;

public:
    Impl()
    {}

    auto setName(const std::string& name) -> void
    {
        name$ = name;
    }

    auto setSpecies(const std::vector<Species>& species) -> void
    {
        species$ = species;
    }

    auto setConcentration(const ConcentrationFn& concentration) -> void
    {
        concentration$ = concentration;
    }

    auto setActivity(const ActivityFn& activity) -> void
    {
        activity$ = activity;
    }

    auto name() const -> const std::string&
    {
        return name$;
    }

    auto numSpecies() const -> unsigned
    {
        return species$.size();
    }

    auto idxSpecies(const std::string& name) const -> Index
    {
        const auto compare = [&](const Species& s) { return s.name() == name; };
        return std::find_if(species$.begin(), species$.end(), compare) - species$.begin();
    }

    auto species(const Index& idxSpecies) const -> const Species&
    {
        return species$[idxSpecies];
    }

    auto species(const Index& idxSpecies) -> Species&
    {
        return species$[idxSpecies];
    }

    auto species(const std::string& name) const -> const Species&
    {
        return species(idxSpecies(name));
    }

    auto species(const std::string& name) -> Species&
    {
        return species(idxSpecies(name));
    }

    auto species() const -> const std::vector<Species>&
    {
        return species$;
    }

    auto species() -> std::vector<Species>&
    {
        return species$;
    }

    auto speciesNames() const -> std::vector<std::string>
    {
        return names(species$);
    }

    auto chemicalPotentials(double T, double P) const -> Vector
    {
        Vector mu0(numSpecies());

        for(unsigned i = 0; i < numSpecies(); ++i)
            mu0[i] = species$[i].chemicalPotential(T, P);

        return mu0;
    }

    auto molarFractions(const Vector& n) const -> Vector
    {
        const double ntotal = n.sum();

        if(ntotal == 0.0) return zeros(n.rows());

        return n/ntotal;
    }

    auto concentrations(const Vector& n) const -> Vector
    {
        return concentration$(n);
    }

    auto activities(double T, double P, const Vector& n) const -> VectorResult
    {
        return activity$(T, P, n);
    }
};

Phase::Phase()
: pimpl(new Impl())
{}

Phase::Phase(const Phase& other)
: pimpl(new Impl(*other.pimpl))
{}

Phase::~Phase()
{}

auto Phase::operator=(Phase other) -> Phase&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Phase::setName(const std::string& name) -> Phase&
{
    pimpl->setName(name);
    return *this;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> Phase&
{
    pimpl->setSpecies(species);
    return *this;
}

auto Phase::setConcentrationFn(const ConcentrationFn& concentration) -> Phase&
{
    pimpl->setConcentration(concentration);
    return *this;
}

auto Phase::setActivityFn(const ActivityFn& activity) -> Phase&
{
    pimpl->setActivity(activity);
    return *this;
}

auto Phase::name() const -> const std::string&
{
    return pimpl->name();
}

auto Phase::numSpecies() const -> unsigned
{
    return pimpl->numSpecies();
}

auto Phase::idxSpecies(const std::string& name) const -> Index
{
    return pimpl->idxSpecies(name);
}

auto Phase::species(const Index& idxSpecies) const -> const Species&
{
    return pimpl->species(idxSpecies);
}

auto Phase::species(const Index& idxSpecies) -> Species&
{
    return pimpl->species(idxSpecies);
}

auto Phase::species(const std::string& name) const -> const Species&
{
    return pimpl->species(name);
}

auto Phase::species(const std::string& name) -> Species&
{
    return pimpl->species(name);
}

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->species();
}

auto Phase::species() -> std::vector<Species>&
{
    return pimpl->species();
}

auto Phase::speciesNames() const -> std::vector<std::string>
{
    return pimpl->speciesNames();
}

auto Phase::chemicalPotentials(double T, double P) const -> Vector
{
    return pimpl->chemicalPotentials(T, P);
}

auto Phase::molarFractions(const Vector& n) const -> Vector
{
    return pimpl->molarFractions(n);
}

auto Phase::concentrations(const Vector& n) const -> Vector
{
    return pimpl->concentrations(n);
}

auto Phase::activities(double T, double P, const Vector& n) const -> VectorResult
{
    return pimpl->activities(T, P, n);
}

auto Phase::operator==(const Phase& phase) const -> bool
{
    return name() == phase.name();
}

auto operator<<(std::ostream& out, const Phase& phase) -> std::ostream&
{
    out << phase.name() << std::endl;
    out << "  ";
    const auto& species = phase.species();
    for(unsigned i = 0; i < species.size(); ++i)
        out << (i > 0 ? ", " : "") << species[i].name();
    return out;
}

} // namespace Reaktor
