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
#include <vector>

// Boost includes
#include <boost/optional.hpp>

// Reaktor includes
#include <Reaktor/Core/ReactionEquation.hpp>
#include <Reaktor/Math/BilinearInterpolator.hpp>

namespace Reaktor {

/**
 * Stores discrete thermodynamic properties of a species over a range of temperatures and pressures
 */
struct ThermoDataSpecies
{
    /**
     * The interpolator of the apparent standard molar Gibbs free energy of formation @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator gibbs;

    /**
     * The interpolator of the apparent standard molar enthalpy of formation @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator enthalpy;

    /**
     * The interpolator of the standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator entropy;

    /**
     * The interpolator of the standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator volume;

    /**
     * The interpolator of the standard molar isobaric heat capacity @f$ C_{p}^{\circ}@f$ of the species (in units of J/(mol*K))
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator cp;
};

/**
 * Stores discrete thermodynamic properties of a reaction over a range of temperatures and pressures
 */
struct ThermoDataReaction
{
    /**
     * The reaction equation for which the thermodynamic properties are known
     *
     * For example, if the apparent standard molar Gibbs free energy of formation of calcite
     * is to be calculated in terms of the log(K) of its reaction, then the reaction equation
     * could be defined as:
     *
     * ~~~
     * ReactionEquation equation = "-1:Calcite, -1:H+, 1:Ca++, 1:HCO3-";
     * ~~~
     */
    ReactionEquation equation;

    /**
     * The interpolator of the equilibrium constant @f$\log K_{\mathrm{eq}}@f$ of the reaction
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator logK;

    /**
     * The interpolator of the equilibrium constant @f$\ln K_{\mathrm{eq}}@f$ of the reaction
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator lnK;

    /**
     * The interpolator of the equilibrium constant @f$\mathrm{p}K_{\mathrm{eq}}@f$ of the reaction
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator pK;

    /**
     * The interpolator of the standard molar Gibbs free energy @f$\Delta G_{r}^{\circ}@f$ of the reaction (in units of J/mol)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator gibbs;

    /**
     * The interpolator of the standard molar enthalpy @f$\Delta H_{r}^{\circ}@f$ of the reaction (in units of J/mol)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator enthalpy;

    /**
     * The interpolator of the standard molar entropy @f$\Delta S_{r}^{\circ}@f$ of the reaction (in units of J/K)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator entropy;

    /**
     * The interpolator of the standard molar volume @f$\Delta V_{r}^{\circ}@f$ of the reaction (in units of m3/mol)
     *
     * The temperature values (in units of K) and pressure values (in units of Pa) are
     * along the x and y coordinates of the bilinear interpolator.
     */
    BilinearInterpolator volume;
};

struct AqueousThermoDataHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    double Sr;

    /// The coefficient a1 of the HKF equation of state of the aqueous species (in units of cal/(mol*bar))
    double a1;

    /// The coefficient a2 of the HKF equation of state of the aqueous species (in units of cal/mol)
    double a2;

    /// The coefficient a3 of the HKF equation of state of the aqueous species (in units of (cal*K)/(mol*bar))
    double a3;

    /// The coefficient a4 of the HKF equation of state of the aqueous species (in units of (cal*K)/mol)
    double a4;

    /// The coefficient c1 of the HKF equation of state of the aqueous species (in units of cal/(mol*K))
    double c1;

    /// The coefficient c2 of the HKF equation of state of the aqueous species (in units of (cal*K)/mol)
    double c2;

    /// The conventional Born coefficient of the aqueous species at reference temperature 298.15 K and pressure 1 bar (in units of cal/mol)
    double wref;
};

struct GaseousThermoDataHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol�K))
    double Sr;

    /// The coefficient a of the HKF equation of state of the gaseous species (in units of cal/(mol�K))
    double a;

    /// The coefficient b of the HKF equation of state of the gaseous species (in units of cal/(mol�K^2))
    double b;

    /// The coefficient c of the HKF equation of state of the gaseous species (in units of (cal�K)/mol)
    double c;

    /// The maximum temperature at which the HKF equation of state can be applied for the gaseous species (in units of K)
    double Tmax;
};

struct MineralThermoDataHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol�K))
    double Sr;

    /// The standard molal volume of the mineral species at reference temperature and pressure (in units of cm^3/mol)
    double Vr;

    /// The number of phase transitions that the mineral may undergo
    int nptrans;

    /// The coefficients ai of the HKF equation of state of the mineral species for each phase region (in units of cal/(mol�K))
    std::vector<double> a;

    /// The coefficients bi of the HKF equation of state of the mineral species for each phase region (in units of cal/(mol�K^2))
    std::vector<double> b;

    /// The coefficients ci of the HKF equation of state of the mineral species for each phase region (in units of (cal�K)/mol)
    std::vector<double> c;

    /// The temperatures at which the mineral experiences phase transition along the line of reference pressure (in units of K)
    std::vector<double> Ttr;

    /// The change in the standard molal enthalpy of each mineral phase transition (in units of cal/mol)
    std::vector<double> Htr;

    /// The change in the standard molal volume of each mineral phase transition (in units of cm^3/mol)
    std::vector<double> Vtr;

    /// The Clapeyron slote at each mineral phase transition (in units of bar/K)
    std::vector<double> dPdTtr;

    /// The maximum temperature at which the HKF equation of state can be applied for the mineral (in units of K)
    double Tmax;
};

struct ThermoData
{
    /// The thermodynamic data of the species in terms of interpolated standard chemical potentials
    boost::optional<ThermoDataSpecies> interpolated;

    /// The thermodynamic data of the species in terms of interpolated equilibrium constants
    boost::optional<ThermoDataReaction> reaction;
};

struct AqueousThermoData : public ThermoData
{
    /// The thermodynamic data of the aqueous species from the HKF model
    boost::optional<AqueousThermoDataHKF> hkf;
};

struct GaseousThermoData : public ThermoData
{
    /// The thermodynamic data of the gaseous species from the HKF model
    boost::optional<GaseousThermoDataHKF> hkf;
};

struct MineralThermoData : public ThermoData
{
    /// The thermodynamic data of the aqueous species from the HKF model
    boost::optional<MineralThermoDataHKF> hkf;
};

} /* namespace Reaktor */
