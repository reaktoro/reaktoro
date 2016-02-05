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
int main()
{
//    double n1 = 3.0;
//    double n2 = 6.0;
//
//    ChemicalVector u(2, 2), v(2, 2);
//    u.val << n1*n2, n1+n2;
//    u.ddn << n2, n1, 1, 1;
//    v.val << n1*n1, n2*n2;
//    v.ddn << 2*n1, 0, 0, 2*n2;
//
//    ChemicalVector t = u % v;
//    ChemicalScalar s = sum(u % v);
//
//    ChemicalVector texp(2, 2);
//    texp.val << n1*n1*n1*n2, (n1+n2)*n2*n2;
//    texp.ddn << 3*n1*n1*n2, n1*n1*n1, n2*n2, 2*n2*(n1+n2)+n2*n2;
//
//    ChemicalScalar sexp(2);
//    sexp.val = n1*n1*n1*n2 + n2*n2*(n1+n2);
//    sexp.ddn << 3*n1*n1*n2 + n2*n2, n1*n1*n1 + 2*n2*(n1+n2)+n2*n2;
//
//    std::cout << "t.val(n) = \n" << t.val << std::endl;
//    std::cout << "t.val(e) = \n" << texp.val << std::endl;
//
//    std::cout << "t.ddn(n) = \n" << t.ddn << std::endl;
//    std::cout << "t.ddn(e) = \n" << texp.ddn << std::endl;
//
//    std::cout << "s.val(n) = \n" << s.val << std::endl;
//    std::cout << "s.val(e) = \n" << sexp.val << std::endl;
//
//    std::cout << "s.ddn(n) = \n" << s.ddn << std::endl;
//    std::cout << "s.ddn(e) = \n" << sexp.ddn << std::endl;
//
//    return 0;


    Database database("databases/supcrt/supcrt98.xml");

    ChemicalEditor editor(database);
    editor.addAqueousPhase("H2O NaCl CaCO3");
    editor.addGaseousPhase("H2O(g) CO2(g)");
    editor.addMineralPhase("Calcite");

    ChemicalSystem system(editor);

    std::cout << system << std::endl;

    EquilibriumProblem problem(system);
    problem.add("H2O", 1, "kg");
    problem.add("NaCl", 0.1, "mol");
//    problem.add("CaCO3", 10.0, "mol");
//    problem.add("HCl", 2.0, "mol");
//    problem.add("CO2", 1, "umol");
//    problem.setSpeciesAmount("Calcite", 10, "moles");
//    problem.setSpeciesAmount("CO2(g)", 10, "moles");
//    problem.setSpeciesAmount("H2O(l)", 55, "moles");
//    problem.setPhaseAmount("Aqueous", 60, "moles", "H2O");
//    problem.setPhaseAmount("Aqueous", 60, "moles")
//        .titrateWith("1 kg H2O; 1 mol NaCl");
//    problem.setPhaseAmount("Calcite", 10, "moles");
//    problem.setPhaseAmount("Gaseous", 10, "moles");
//    problem.pH(4.0).titrateWith("CO2");
//    problem.pH(4.0).titrateWithEither("HCl", "NaOH");
//    problem.setPhaseVolume("Aqueous", 0.5. "m3").
//        titrateWith("1 kg H2O; 1 mol NaCl");

//    problem.pH(3.0, "HCl");
//    problem.pH(3.0, "HCl", "NaOH");
    problem.setPhaseVolume("Gaseous", 0.5, "m3", "(1:mol:CO2)");
    problem.setPhaseVolume("Aqueous", 0.5, "m3", "(1:kg:H2O)(0.1:mol:NaCl)");
//    problem.setPhaseVolume("Aqueous", 0.5, "m3", "H2O");
//    problem.setPhaseAmount("Aqueous", 100, "mol", "H2O");
//    problem.setPhaseAmount("Aqueous", 100, "mol", "(1:kg:H2O)(0.1:mol:NaCl)");
//    problem.setPhaseVolume("Calcite", 0.5, "m3", "CaCO3");

    EquilibriumOptions options;
    options.optimum.output.active = true;
    options.hessian = GibbsHessian::Exact;

    ChemicalState state = equilibrate(problem, options);

    std::cout << state << std::endl;
}
