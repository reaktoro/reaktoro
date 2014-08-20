/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <functional>
#include <string>
#include <vector>

namespace Reaktor {

// Reaktor forward declarations
class Database;
class ReactionEquation;
class BilinearInterpolator;
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;

/**
 * Generates a interpolation function for the chemical potential of a species in the database
 * The default temperature and pressure points are:
 *  T(celsius) = 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300
 *  P(bar) = 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600
 * @param database The database instance
 * @param species The name of the species
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see Database, BilinearInterpolator
 */
auto chemicalPotentialFn(const Database& database, std::string species) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of a species in the database
 * @param database The database instance
 * @param species The name of the species
 * @param tPoints The temperature points for the interpolation
 * @param pPoints The pressure points for the interpolation (in units of Pa)
 * @param tUnit The string representing the unit of the temperature points (must be convertible to K)
 * @param pUnit The string representing the unit of the pressure points (must be convertible to Pa)
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see Database, BilinearInterpolator
 */
auto chemicalPotentialFn(const Database& database, std::string species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of an aqueous species
 * The default temperature and pressure points are:
 *  T(celsius) = 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300
 *  P(bar) = 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600
 * @param species The instance of the aqueous species
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see AqueousSpecies, BilinearInterpolator
 */
auto chemicalPotentialFn(const AqueousSpecies& species) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of an aqueous species
 * @param species The instance of the aqueous species
 * @param tPoints The temperature points for the interpolation
 * @param pPoints The pressure points for the interpolation (in units of Pa)
 * @param tUnit The string representing the unit of the temperature points (must be convertible to K)
 * @param pUnit The string representing the unit of the pressure points (must be convertible to Pa)
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see AqueousSpecies, BilinearInterpolator
 */
auto chemicalPotentialFn(const AqueousSpecies& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of a gaseous species
 * The default temperature and pressure points are:
 *  T(celsius) = 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300
 *  P(bar) = 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600
 * @param species The instance of the gaseous species
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see GaseousSpecies, BilinearInterpolator
 */
auto chemicalPotentialFn(const GaseousSpecies& species) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of a gaseous species
 * @param species The instance of the gaseous species
 * @param tPoints The temperature points for the interpolation
 * @param pPoints The pressure points for the interpolation (in units of Pa)
 * @param tUnit The string representing the unit of the temperature points (must be convertible to K)
 * @param pUnit The string representing the unit of the pressure points (must be convertible to Pa)
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see GaseousSpecies, BilinearInterpolator
 */
auto chemicalPotentialFn(const GaseousSpecies& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of a mineral species
 * The default temperature and pressure points are:
 *  T(celsius) = 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300
 *  P(bar) = 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600
 * @param species The instance of the mineral species
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see MineralSpecies, BilinearInterpolator
 */
auto chemicalPotentialFn(const MineralSpecies& species) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the chemical potential of a mineral species
 * @param species The instance of the mineral species
 * @param tPoints The temperature points for the interpolation
 * @param pPoints The pressure points for the interpolation (in units of Pa)
 * @param tUnit The string representing the unit of the temperature points (must be convertible to K)
 * @param pUnit The string representing the unit of the pressure points (must be convertible to Pa)
 * @return A bilinear interpolation function for the chemical potential of the given species (in units of J/mol)
 * @see MineralSpecies, BilinearInterpolator
 */
auto chemicalPotentialFn(const MineralSpecies& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the equilibrium constant of the reaction
 * The default temperature and pressure points are:
 *  T(celsius) = 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300
 *  P(bar) = 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600
 * @param database The database containing the thermodynamic parameters of the species
 * @param equation The equation of the reaction
 * @return A bilinear interpolation function for the equilibrium constant of the reaction equation
 * @see Database, ReactionEquation, BilinearInterpolator
 */
auto createEquilibriumConstant(const Database& database, const ReactionEquation& equation) -> BilinearInterpolator;

/**
 * Generates a interpolation function for the equilibrium constant of the reaction
 * @param database The database containing the thermodynamic parameters of the species
 * @param equation The equation of the reaction
 * @param tPoints The temperature points for the interpolation
 * @param pPoints The pressure points for the interpolation (in units of Pa)
 * @param tUnit The string representing the unit of the temperature points (must be convertible to K)
 * @param pUnit The string representing the unit of the pressure points (must be convertible to Pa)
 * @return A bilinear interpolation function for the equilibrium constant of the reaction equation
 * @see Database, ReactionEquation, BilinearInterpolator
 */
auto createEquilibriumConstant(const Database& database, const ReactionEquation& equation,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator;

} /* namespace Reaktor */
