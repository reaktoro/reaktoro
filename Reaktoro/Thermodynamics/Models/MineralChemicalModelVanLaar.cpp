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
//
//#include "MineralChemicalModelIdeal.hpp"
//
//// Reaktoro includes
//#include <Reaktoro/Common/Constants.hpp>
//#include <Reaktoro/Common/TableUtils.hpp>
//#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>
//
//namespace Reaktoro {
//
//auto mineralChemicalModelVanLaar(const MineralMixture& mixture, VectorConstRef a, MatrixConstRef W) -> PhaseChemicalModel
//{
//    MineralMixtureState state;
//    PhaseChemicalModelResult res(2);
//
//    VectorXd phi;
//
//    ChemicalScalar avg;
//
//    const Index nspecies = mixture.numSpecies();
//
//    Table2D<ChemicalScalar> B = table2D<ChemicalScalar>(nspecies, nspecies);
//
//    PhaseChemicalModel f = [=](double T, double P, VectorConstRef n) mutable
//    {
//        state = mixture.state(T, P, n);
//
//        const auto RT = universalGasConstant * state.T;
//        const auto& x = state.x;
//
//        phi = a % x;
//        avg = sum(phi);
//        phi /= avg;
//
//        for(Index i = 0; i < nspecies - 1; ++i)
//            for(Index j = i + 1; j < nspecies; ++j)
//                B[i][j] = 2*avg/(a[i] + a[j]) * W(i, j);
//
//        for(Index k = 0; k < nspecies; ++k)
//        {
//
//        }
//
////        const auto x1 = state.x[0];
////        const auto x2 = state.x[1];
////
////        res.ln_activity_coefficients[0] = x2*x2*(a0 + a1*(3*x1 - x2) + a2*(x1 - x2)*(5*x1 - x2));
////        res.ln_activity_coefficients[1] = x1*x1*(a0 - a1*(3*x2 - x1) + a2*(x2 - x1)*(5*x2 - x1));
////
////        res.ln_activities = res.ln_activity_coefficients + log(state.x);
////
////        res.residual_molar_gibbs_energy = (x1*x2*(a0 + a1*(x1 - x2) + a2*pow((x1 - x2), 2))) * RT;
////        res.residual_molar_enthalpy = res.residual_molar_gibbs_energy;
//
//        return res;
//    };
//
//    return f;
//}
//
//} // namespace Reaktoro
//
//
//
