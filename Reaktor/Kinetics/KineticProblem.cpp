//// Reaktor is a C++ library for computational reaction modelling.
////
//// Copyright (C) 2014 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#include "KineticProblem.hpp"
//
//// Reaktor includes
//#include <Reaktor/Core/ChemicalSystem.hpp>
//#include <Reaktor/Core/Partition.hpp>
//#include <Reaktor/Core/Reaction.hpp>
//
//namespace Reaktor {
//
//struct KineticProblem::Impl
//{
//    const ChemicalSystem& system;
//
//    const Reactions& reactions;
//
//    const Partition& partition;
//
//    double temperature;
//
//    double pressure;
//
//    Vector n;
//
//    double initial_time;
//
//    double final_time;
//
//    Impl(const ChemicalSystem& system, const Reactions& reactions, const Partition& partition)
//    : system(system), reactions(reactions), partition(partition),
//      temperature(INFINITY), pressure(INFINITY),
//      initial_time(0), final_time(INFINITY)
//    {}
//};
//
//KineticProblem::KineticProblem(const ChemicalSystem& system, const Reactions& reactions)
//: KineticProblem(system, reactions, Partition(system))
//{}
//
//KineticProblem::KineticProblem(const ChemicalSystem& system, const Reactions& reactions, const Partition& partition)
//: pimpl(new Impl(system, reactions, partition))
//{}
//
//KineticProblem::KineticProblem(const KineticProblem& other)
//: pimpl(new Impl(*other.pimpl))
//{}
//
//KineticProblem::~KineticProblem()
//{}
//
//auto KineticProblem::operator=(KineticProblem other) -> KineticProblem&
//{
//    pimpl = std::move(other.pimpl);
//    return *this;
//}
//
//auto KineticProblem::setTemperature(double val) -> KineticProblem&
//{
//    pimpl->temperature = val;
//    return *this;
//}
//
//auto KineticProblem::setPressure(double val) -> KineticProblem&
//{
//    pimpl->pressure = val;
//    return *this;
//}
//
//auto KineticProblem::setInitialState(const Vector& n) -> KineticProblem&
//{
//    pimpl->n = n;
//    return *this;
//}
//
//auto KineticProblem::setInitialTime(double val) -> KineticProblem&
//{
//    pimpl->initial_time = val;
//    return *this;
//}
//
//auto KineticProblem::setFinalTime(double val) -> KineticProblem&
//{
//    pimpl->final_time = val;
//     return *this;
//}
//
//auto KineticProblem::temperature() const -> double
//{
//    return pimpl->temperature;
//}
//
//auto KineticProblem::pressure() const -> double
//{
//    return pimpl->pressure;
//}
//
//auto KineticProblem::initialState() const -> const Vector&
//{
//    return pimpl->n;
//}
//
//auto KineticProblem::initialTime() const -> double
//{
//    return pimpl->initial_time;
//}
//
//auto KineticProblem::finalTime() const -> double
//{
//    return pimpl->final_time;
//}
//
//} // namespace Reaktor
//
