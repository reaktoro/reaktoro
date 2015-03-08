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
//#pragma once
//
//// C++ includes
//#include <vector>
//#include <string>
//
//namespace Reaktor {
//
//// Forward declarations
//class Database;
//class AqueousSpecies;
//class GaseousSpecies;
//class MineralSpecies;
//class SpeciesThermoProperties;
//class BilinearInterpolator;
//
//class Thermo
//{
//public:
//    explicit Thermo(const Database& database);
//
//    auto setTemperatures(const std::vector<double>& temperatures) -> void;
//
//    auto setPressures(const std::vector<double>& pressures) -> void;
//
//    auto database() const -> const Database&;
//
//    auto temperatures() const -> const std::vector<double>&;
//
//    auto pressures() const -> const std::vector<double>&;
//
//    auto standardGibbsEnergy(const AqueousSpecies& species) const -> ThermoScalarFunction;
//
//    auto standardHelmholtzEnergy(const AqueousSpecies& species) const -> ThermoScalarFunction;
//
//    auto standardInternalEnergy(const AqueousSpecies& species) const -> ThermoScalarFunction;
//
//    auto standardEnthalpy(const AqueousSpecies& species) const -> ThermoScalarFunction;
//
//    auto standardEntropy(const AqueousSpecies& species) const -> ThermoScalarFunction;
//
//    auto standardVolume(const AqueousSpecies& species) const -> ThermoScalarFunction;
//
//    auto thermoDataHKF(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataHKF(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataHKF(const MineralSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataReaction(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataReaction(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataReaction(const MineralSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataInterpolated(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataInterpolated(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoDataInterpolated(const MineralSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoData(const std::string& species) const -> SpeciesThermoProperties;
//
//    auto thermoData(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoData(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//
//    auto thermoData(const MineralSpecies& species) const -> SpeciesThermoProperties;
//
//private:
//    /// The database instance
//    const Database& database$;
//
//    /// The temperatures for the interpolation (in units of K)
//    std::vector<double> temperatures$;
//
//    /// The pressures for the interpolation (in units of Pa)
//    std::vector<double> pressures$;
//};
//
// namespace Reaktor
