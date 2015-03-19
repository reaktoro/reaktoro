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
//#include <Reaktor/Common/Exception.hpp>
//#include <Reaktor/Common/ReactionEquation.hpp>
//#include <Reaktor/Common/ReactionEquation.hpp>
//#include <Reaktor/Math/BilinearInterpolator.hpp>
//#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
//#include <Reaktor/Thermodynamics/Species/GaseousSpecies.hpp>
//#include <Reaktor/Thermodynamics/Species/GeneralSpecies.hpp>
//#include <Reaktor/Thermodynamics/Species/MineralSpecies.hpp>
//#include <Reaktor/Thermodynamics/Species/SpeciesThermoProperties.hpp>
//#include <Reaktor/Thermodynamics/Core/Database.hpp>
//#include <Reaktor/Thermodynamics/Core/Thermo.hpp>
//#include <Reaktor/Thermodynamics/Models/SpeciesThermoState.hpp>
//#include <Reaktor/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>
//
//namespace Reaktor {
//namespace internal {
//
//template<typename SpeciesType>
//auto thermoDataHKF(const Thermo& editor, const SpeciesType& species) -> SpeciesThermoProperties
//{
//    const auto& temperatures = editor.temperatures();
//    const auto& pressures = editor.pressures();
//
//    std::vector<double> gibbs;
//    std::vector<double> enthalpy;
//    std::vector<double> entropy;
//    std::vector<double> volume;
//    std::vector<double> cp;
//
//    for(auto P : pressures) for(auto T : temperatures)
//    {
//        SpeciesThermoState state = SpeciesThermoStateHKF(T, P, species);
//
//        gibbs.push_back(state.gibbs_energy);
//        enthalpy.push_back(state.enthalpy);
//        entropy.push_back(state.entropy);
//        volume.push_back(state.volume);
//        cp.push_back(state.heat_capacity);
//    }
//
//    SpeciesThermoProperties data;
//    data.gibbs_energy = BilinearInterpolator(temperatures, pressures, gibbs);
//    data.enthalpy     = BilinearInterpolator(temperatures, pressures, enthalpy);
//    data.entropy      = BilinearInterpolator(temperatures, pressures, entropy);
//    data.volume       = BilinearInterpolator(temperatures, pressures, volume);
//    data.cp           = BilinearInterpolator(temperatures, pressures, cp);
//
//    return data;
//}
//
//template<typename SpeciesType>
//auto thermoDataReaction(const Thermo& editor, const SpeciesType& species) -> SpeciesThermoProperties
//{
//    // Get a reference to the thermodynamic data of the species
//    const ThermoDataReaction& thermo_data = species.thermoData().reaction.get();
//
//    // The thermodynamic data of each species in the reaction
//    std::vector<std::pair<double, SpeciesThermoProperties>> thermo_data_each_species;
//
//    // The stoichiometry of the species in the reaction for which thermodynamic properties are being calculated
//    double stoich = 0.0;
//
//    // Iterate over all species in the reaction skip the dependent species
//    for(const auto& name_stoichiometry : thermo_data.equation)
//    {
//        const auto& name = name_stoichiometry.first;
//        const auto& stoichiometry = name_stoichiometry.second;
//
//        if(name != species.name())
//            thermo_data_each_species.push_back(
//                {stoichiometry, editor.thermoData(name)});
//        else
//            stoich = stoichiometry;
//    }
//
//    auto gibbs = [&](double T, double P)
//    {
//        double res = thermo_data.gibbs(T, P);
//        for(const auto& pair : thermo_data_each_species)
//            res -= pair.first * pair.second.gibbs(T, P);
//        return res/stoich;
//    };
//
//    auto enthalpy = [&](double T, double P)
//    {
//        double res = thermo_data.enthalpy(T, P);
//        for(const auto& pair : thermo_data_each_species)
//            res -= pair.first * pair.second.enthalpy(T, P);
//        return res/stoich;
//    };
//
//    auto entropy = [&](double T, double P)
//    {
//        double res = thermo_data.entropy(T, P);
//        for(const auto& pair : thermo_data_each_species)
//            res -= pair.first * pair.second.entropy(T, P);
//        return res/stoich;
//    };
//
//    auto volume = [&](double T, double P)
//    {
//        double res = thermo_data.volume(T, P);
//        for(const auto& pair : thermo_data_each_species)
//            res -= pair.first * pair.second.volume(T, P);
//        return res/stoich;
//    };
//
//    const auto& temperatures = editor.temperatures();
//    const auto& pressures = editor.pressures();
//
//    SpeciesThermoProperties data;
//
//    if(thermo_data.gibbs.initialised())
//        data.gibbs = BilinearInterpolator(temperatures, pressures, gibbs);
//
//    if(thermo_data.enthalpy.initialised())
//        data.enthalpy = BilinearInterpolator(temperatures, pressures, enthalpy);
//
//    if(thermo_data.entropy.initialised())
//        data.entropy = BilinearInterpolator(temperatures, pressures, entropy);
//
//    if(thermo_data.volume.initialised())
//        data.volume = BilinearInterpolator(temperatures, pressures, volume);
//
//    return data;
//}
//
//template<typename SpeciesType>
//auto thermoData(const Thermo& editor, const SpeciesType& species) -> SpeciesThermoProperties
//{
//    const auto& thermo_data = species.thermoData();
//
//    if(thermo_data.interpolated)
//        return editor.thermoDataInterpolated(species);
//
//    if(thermo_data.reaction)
//        return editor.thermoDataReaction(species);
//
//    if(thermo_data.hkf)
//        return editor.thermoDataHKF(species);
//
//    if(species.name() == "H2O(l)")
//        return editor.thermoDataHKF(species);
//
//    Exception exception;
//    exception.error  << "Unable to generate standard chemical potential function for species " << species.name() << ".";
//    exception.reason << "This species does not have any thermodynamic data available.";
//    RaiseError(exception);
//
//    return SpeciesThermoProperties();
//}
//
//} /* namespace internal */
//
//Thermo::Thermo(const Database& database)
//: database$(database)
//{
//    // The default temperature points for the interpolation of the thermodynamic quantities (in units of celsius)
//    temperatures$ = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300 };
//
//    // The default pressure points for the interpolation of the thermodynamic quantities (in units of bar)
//    pressures$ = { 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600 };
//
//    // Convert the temperatures and pressures to units of kelvin and pascal respectively
//    for(auto& x : temperatures$) x = x + 273.15;
//    for(auto& x : pressures$)    x = x * 1.0e+5;
//}
//
//auto Thermo::setTemperatures(const std::vector<double>& temperatures) -> void;
//{
//    temperatures$ = temperatures;
//}
//
//auto Thermo::setPressures(const std::vector<double>& pressures) -> void;
//{
//    pressures$ = pressures;
//}
//
//auto Thermo::database() const -> const Database&;
//{
//    return database$;
//}
//
//auto Thermo::temperatures() const -> const std::vector<double>&;
//{
//    return temperatures$;
//}
//
//auto Thermo::pressures() const -> const std::vector<double>&;
//{
//    return pressures$;
//}
//
//auto Thermo::thermoDataHKF(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoDataHKF(*this, species);
//}
//
//auto Thermo::thermoDataHKF(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoDataHKF(*this, species);
//}
//
//auto Thermo::thermoDataHKF(const MineralSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoDataHKF(*this, species);
//}
//
//auto Thermo::thermoDataReaction(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoDataReaction(*this, species);
//}
//
//auto Thermo::thermoDataReaction(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoDataReaction(*this, species);
//}
//
//auto Thermo::thermoDataReaction(const MineralSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoDataReaction(*this, species);
//}
//
//auto Thermo::thermoDataInterpolated(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return species.thermoData().interpolated.get();
//}
//
//auto Thermo::thermoDataInterpolated(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return species.thermoData().interpolated.get();
//}
//
//auto Thermo::thermoDataInterpolated(const MineralSpecies& species) const -> SpeciesThermoProperties;
//{
//    return species.thermoData().interpolated.get();
//}
//
//auto Thermo::thermoData(const std::string& species) const -> SpeciesThermoProperties;
//{
//    if(database$.containsAqueousSpecies(species))
//        return thermoData(database$.aqueousSpecies(species));
//
//    if(database$.containsGaseousSpecies(species))
//        return thermoData(database$.gaseousSpecies(species));
//
//    if(database$.containsMineralSpecies(species))
//        return thermoData(database$.mineralSpecies(species));
//
//    Exception exception;
//    exception.error  << "Unable to generate standard chemical potential function for species " << species << ".";
//    exception.reason << "This species is not contained in the database.";
//    RaiseError(exception);
//
//    return SpeciesThermoProperties();
//}
//
//auto Thermo::thermoData(const AqueousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoData(*this, species);
//}
//
//auto Thermo::thermoData(const GaseousSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoData(*this, species);
//}
//
//auto Thermo::thermoData(const MineralSpecies& species) const -> SpeciesThermoProperties;
//{
//    return internal::thermoData(*this, species);
//}
//
// namespace Reaktor
