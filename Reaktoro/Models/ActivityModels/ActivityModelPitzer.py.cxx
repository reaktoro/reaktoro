// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzer.hpp>
using namespace Reaktoro;

void exportActivityModelPitzer(py::module& m)
{
    auto parent = py::class_<ActivityModelParamsPitzer>(m, "ActivityModelParamsPitzer");

    py::enum_<ActivityModelParamsPitzer::CorrectionModel>(parent, "CorrectionModel")
        .value("Constant", ActivityModelParamsPitzer::CorrectionModel::Constant, "Correction model based on constant expression.")
        .value("Phreeqc", ActivityModelParamsPitzer::CorrectionModel::Phreeqc, "Correction model based on formula from PHREEQC v3 where")
        .value("HeMorse1993", ActivityModelParamsPitzer::CorrectionModel::HeMorse1993, "Correction model based on formula from He and Morse (1993) (doi: 10.1016/0016-7037(93)90137-L)")
        .value("Dai2013", ActivityModelParamsPitzer::CorrectionModel::Dai2013, "Correction model based on formula from Dai et al. (2013) (doi: 10.2118/164045-ms)")
        .value("Dai2014", ActivityModelParamsPitzer::CorrectionModel::Dai2014, "Correction model based on formula from Dai et al. (2014) (doi: 10.2118/169786-ms)")
        .value("ChristovMoller2004", ActivityModelParamsPitzer::CorrectionModel::ChristovMoller2004, "Correction model based on formula from Christov and Møller (2004) (see Table 5 in 10.1016/j.chemgeo.2007.07.023)")
        .value("Holmes1987", ActivityModelParamsPitzer::CorrectionModel::Holmes1987, "Correction model based on formula from Holmes et al. (1987) (see Table 6 in 10.1016/j.chemgeo.2007.07.023)")
        .value("Pitzer1984", ActivityModelParamsPitzer::CorrectionModel::Pitzer1984, "Correction model based on formula from Pitzer et al. (1984) (see Table 7 in 10.1016/j.chemgeo.2007.07.023)")
        .value("PalabanPitzer1987", ActivityModelParamsPitzer::CorrectionModel::PalabanPitzer1987, "Correction model based on formula from Palaban and Pitzer (1987) (see Table 8 in 10.1016/j.chemgeo.2007.07.023)")
        .value("Polya2001", ActivityModelParamsPitzer::CorrectionModel::Polya2001, "Correction model based on formula from Polya et al. (2001) (see Table 10 in 10.1016/j.chemgeo.2007.07.023)")
        .value("LiDuan2007", ActivityModelParamsPitzer::CorrectionModel::LiDuan2007, "Correction model based on formula from Li and Duan (2007) (see Table 12 in 10.1016/j.chemgeo.2007.07.023)")
        ;

    py::class_<ActivityModelParamsPitzer::InteractionParamAttribs>(parent, "InteractionParamAttribs")
        .def(py::init<>())
        .def_readwrite("formulas", &ActivityModelParamsPitzer::InteractionParamAttribs::formulas, "The chemical formulas of the species associated to this species interaction parameter sorted in descending order of charge.")
        .def_readwrite("model", &ActivityModelParamsPitzer::InteractionParamAttribs::model, "The model used for temperature-pressure correction of this species interaction parameter (options).")
        .def_readwrite("parameters", &ActivityModelParamsPitzer::InteractionParamAttribs::parameters, "The parameters for the temperature-pressure correction model of this species interaction parameter.")
        ;

    parent
        .def(py::init<>())
        .def_readwrite("beta0", &ActivityModelParamsPitzer::beta0, "The parameters `beta^{(0)}_{ij}(T, P)` in the Pitzer model for cation-anion interactions.")
        .def_readwrite("beta1", &ActivityModelParamsPitzer::beta1, "The parameters `beta^{(1)}_{ij}(T, P)` in the Pitzer model for cation-anion interactions.")
        .def_readwrite("beta2", &ActivityModelParamsPitzer::beta2, "The parameters `beta^{(2)}_{ij}(T, P)` in the Pitzer model for cation-anion interactions.")
        .def_readwrite("Cphi", &ActivityModelParamsPitzer::Cphi, "The parameters `C^{phi}_{ij}(T, P)` in the Pitzer model for cation-anion interactions.")
        .def_readwrite("theta", &ActivityModelParamsPitzer::theta, "The parameters `theta_{ij}(T, P)` in the Pitzer model for cation-cation and anion-anion interactions.")
        .def_readwrite("psi", &ActivityModelParamsPitzer::psi, "The parameters `psi_{ijk}(T, P)` in the Pitzer model for cation-cation-anion and anion-anion-cation interactions.")
        .def_readwrite("lambda", &ActivityModelParamsPitzer::lambda, "The parameters `lambda_{ij}(T, P)` in the Pitzer model for neutral-cation and neutral-anion interactions.")
        .def_readwrite("zeta", &ActivityModelParamsPitzer::zeta, "The parameters `zeta_{ijk}(T, P)` in the Pitzer model for neutral-cation-anion interactions.")
        .def_readwrite("mu", &ActivityModelParamsPitzer::mu, "The parameters `mu_{ijk}(T, P)` in the Pitzer model for neutral-neutral-neutral, neutral-neutral-cation, and neutral-neutral-anion interactions.")
        .def_readwrite("eta", &ActivityModelParamsPitzer::eta, "The parameters `eta_{ijk}(T, P)` in the Pitzer model for neutral-cation-cation and neutral-anion-anion interactions.")
        .def_readwrite("alpha1", &ActivityModelParamsPitzer::alpha1, "The parameters `alpha_1_{ij}` associated to the parameters `beta^{(1)}_{ij}`.")
        .def_readwrite("alpha2", &ActivityModelParamsPitzer::alpha2, "The parameters `alpha_2_{ij}` associated to the parameters `beta^{(2)}_{ij}`.")
        ;

    m.def("ActivityModelPitzer", py::overload_cast<>(ActivityModelPitzer));
    m.def("ActivityModelPitzer", py::overload_cast<ActivityModelParamsPitzer const&>(ActivityModelPitzer));
    m.def("ActivityModelPitzer", py::overload_cast<Params const&>(ActivityModelPitzer));
}
