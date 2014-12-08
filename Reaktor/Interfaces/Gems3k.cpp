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

#include "Gems3k.hpp"

namespace Reaktor {
namespace internal {

double temperature(const GemsHandle& node)
{
    return node->Get_TK();
}

double pressure(const GemsHandle& node)
{
    return node->Get_P();
}

unsigned numComponents(const GemsHandle& node)
{
    return node->pCSD()->nIC;
}

unsigned numSpecies(const GemsHandle& node)
{
    return node->pCSD()->nDC;
}

unsigned numPhases(const GemsHandle& node)
{
    return node->pCSD()->nPH;
}

unsigned numSpeciesInPhase(const GemsHandle& node, unsigned iphase)
{
    return node->pCSD()->nDCinPH[iphase];
}

std::string componentName(const GemsHandle& node, unsigned icomponent)
{
    return node->pCSD()->ICNL[icomponent];
}

std::string speciesName(const GemsHandle& node, unsigned ispecies)
{
    return node->pCSD()->DCNL[ispecies];
}

std::string phaseName(const GemsHandle& node, unsigned iphase)
{
    return node->pCSD()->PHNL[iphase];
}

unsigned componentIndex(const GemsHandle& node, std::string component)
{
    unsigned index = 0;
    const unsigned size = numComponents(node);
    for(; index < size; ++index)
        if(componentName(node, index) == component)
            return index;
    return size;
}

unsigned speciesIndex(const GemsHandle& node, std::string species)
{
    unsigned index = 0;
    const unsigned size = numSpecies(node);
    for(; index < size; ++index)
        if(speciesName(node, index) == species)
            return index;
    return size;
}

unsigned phaseIndex(const GemsHandle& node, std::string phase)
{
    unsigned index = 0;
    const unsigned size = numPhases(node);
    for(; index < size; ++index)
        if(phaseName(node, index) == phase)
            return index;
    return size;
}

double componentMolarMass(const GemsHandle& node, unsigned icomponent)
{
    return node->ICmm(icomponent);
}

double speciesMolarMass(const GemsHandle& node, unsigned ispecies)
{
    return node->DCmm(ispecies);
}

Vector speciesAmounts(const GemsHandle& node)
{
    Vector n(numSpecies(node));
    for(unsigned i = 0; i < n.size(); ++i)
        n[i] = node->Get_nDC(i);
    return n;
}

Vector componentAmounts(const GemsHandle& node)
{
    Vector b(numComponents(node));
    for(unsigned i = 0; i < b.size(); ++i)
        b[i] = node->Get_bIC(i);
    return b;
}

Matrix formulaMatrix(const GemsHandle& node)
{
    const unsigned N = numSpecies(node);
    const unsigned E = numComponents(node);
    Matrix A(E, N);
    for(unsigned i = 0; i < E; ++i)
        for(unsigned j = 0; j < N; ++j)
            A(i, j) = node->DCaJI(j, i);
    return A;
}

Vector standardChemicalPotentials(GemsHandle& node, double T, double P)
{
    const unsigned num_species = numSpecies(node);
    Vector u0(num_species);
    node->setTemperature(T);
    node->setPressure(P);
    node->updateStandardGibbsEnergies();
    ACTIVITY* ap = node->pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        u0[i] = ap->tpp_G[i];
    return u0;
}

Vector chemicalPotentials(GemsHandle& node, double T, double P, const Vector& n)
{
    const unsigned num_species = numSpecies(node);
    Vector u(num_species);

    if(T != internal::temperature(node) or P != internal::pressure(node))
    {
        node->setTemperature(T);
        node->setPressure(P);
        node->updateStandardGibbsEnergies();
    }

    node->setTemperature(T);
    node->setPressure(P);
    node->updateStandardGibbsEnergies();

    node->initActivityCoefficients();
    node->setSpeciation(n.memptr());
    node->updateConcentrations();
    node->updateActivityCoefficients();
    node->updateChemicalPotentials();
    ACTIVITY* ap = node->pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        u[i] = ap->F[i];
    return u;
}

bool converged(const GemsHandle& node)
{
    const auto status = node->pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA or status == OK_GEM_SIA;
}

} // namespace internal
} // namespace Reaktor
