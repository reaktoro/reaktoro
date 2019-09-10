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


int main() 
{
    Mesh mymesh{};

    TransportSolver transpsolver{};

    transpsolver.setMesh(mymesh);

    transpsolver.setVelocity(1);

    transpsolver.setDiffusionCoeff(0.0);

    transpsolver.setBoundaryValue(1);

    transpsolver.setTimeStep(0.001);
    transpsolver.initialize();

    VectorXd u(mymesh.numCells());
    u = zeros(u.size());

    transpsolver.step(u);

    std::cout << u << std::endl;

    for (int i = 0; i < 10000; i++) {
        transpsolver.step(u);
    }

    transpsolver.step(u);

    std::cout << u << std::endl;
}