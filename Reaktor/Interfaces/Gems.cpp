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

#include "Gems.hpp"

// Gems includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#include <Reaktor/gems3k/node.h>

namespace Reaktor {

struct Gems::Impl
{
    /// The TNode instance from Gems
    TNode node;
};

Gems::Gems()
: pimpl(new Impl())
{}

Gems::Gems(std::string filename)
: pimpl(new Impl())
{
    if(node().GEM_init(filename.c_str()))
        throw std::runtime_error("Error reading the Gems chemical system specification file.");
}

auto Gems::setTemperature(double val) -> void
{
    node().setTemperature(val);
}

auto Gems::setPressure(double val) -> void
{
    node().setPressure(val);
}

auto Gems::setSpeciesAmounts(const Vector& n) -> void
{
    node().setSpeciation(n.memptr());
}

auto Gems::setComponentAmounts(const Vector& b) -> void
{
    for(unsigned i = 0; i < numComponents(); ++i)
        node().Set_IC_b(b[i], i);
}

auto Gems::numComponents() const -> unsigned
{
    return node().pCSD()->nIC;
}

auto Gems::numSpecies() const -> unsigned
{
    return node().pCSD()->nDC;
}

auto Gems::numPhases() const -> unsigned
{
    return node().pCSD()->nPH;
}

auto Gems::numSpeciesInPhase(unsigned index) const -> unsigned
{
    return node().pCSD()->nDCinPH[index];
}

auto Gems::componentName(unsigned index) const -> std::string
{
    return node().pCSD()->ICNL[index];
}

auto Gems::speciesName(unsigned index) const -> std::string
{
    return node().pCSD()->DCNL[index];
}

auto Gems::phaseName(unsigned index) const -> std::string
{
    return node().pCSD()->PHNL[index];
}

auto Gems::componentIndex(std::string name) const -> unsigned
{
    unsigned index = 0;
    const unsigned size = numComponents();
    for(; index < size; ++index)
        if(componentName(index) == name)
            return index;
    return size;
}

auto Gems::speciesIndex(std::string name) const -> unsigned
{
    unsigned index = 0;
    const unsigned size = numSpecies();
    for(; index < size; ++index)
        if(speciesName(index) == name)
            return index;
    return size;
}

auto Gems::phaseIndex(std::string name) const -> unsigned
{
    unsigned index = 0;
    const unsigned size = numPhases();
    for(; index < size; ++index)
        if(phaseName(index) == name)
            return index;
    return size;
}

auto Gems::componentMolarMass(unsigned index) const -> double
{
    return node().ICmm(index);
}

auto Gems::speciesMolarMass(unsigned index) const -> double
{
    return node().DCmm(index);
}

auto Gems::temperature() const -> double
{
    return node().Get_TK();
}

auto Gems::pressure() const -> double
{
    return node().Get_P();
}

auto Gems::componentAmounts() const -> Vector
{
    Vector b(numComponents());
    for(unsigned i = 0; i < b.size(); ++i)
        b[i] = node().Get_bIC(i);
    return b;
}

auto Gems::speciesAmounts() const -> Vector
{
    Vector n(numSpecies());
    for(unsigned i = 0; i < n.size(); ++i)
        n[i] = node().Get_nDC(i);
    return n;
}

auto Gems::balanceMatrix() const -> Matrix
{
    const unsigned N = numSpecies();
    const unsigned C = numComponents();
    Matrix A(C, N);
    for(unsigned i = 0; i < C; ++i)
        for(unsigned j = 0; j < N; ++j)
            A(i, j) = node().DCaJI(j, i);
    return A;
}

auto Gems::gibbsEnergies() -> Vector
{
    const unsigned num_species = numSpecies();
    Vector u0(num_species);
    node().updateStandardGibbsEnergies();
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        u0[i] = ap->tpp_G[i];
    return u0;
}

auto Gems::chemicalPotentials() -> Vector
{
    const unsigned num_species = numSpecies();
    Vector u(num_species);
    node().updateStandardGibbsEnergies();
    node().initActivityCoefficients();
    node().updateConcentrations();
    node().updateActivityCoefficients();
    node().updateChemicalPotentials();
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        u[i] = ap->F[i];
    return u;
}

auto Gems::equilibrate() -> void
{
    node().pCNode()->NodeStatusCH = NEED_GEM_SIA;
    node().GEM_run(false);
}

auto Gems::converged() const -> bool
{
    const auto status = node().pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA or status == OK_GEM_SIA;
}

auto Gems::numIterations() const -> unsigned
{
    return node().pCNode()->IterDone;
}

auto Gems::wallTime() const -> double
{
    return node().GEM_CalcTime();
}

auto Gems::node() -> TNode&
{
    return pimpl->node;
}

auto Gems::node() const -> const TNode&
{
    return pimpl->node;
}

} // namespace Reaktor
