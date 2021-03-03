// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

using autodiff::jacobian;
using autodiff::wrt;
using autodiff::at;

TEST_CASE("Testing EquilibriumProps", "[EquilibriumProps]")
{
    // Create the Param objects for the G0 values of the species.
    // This is needed so that we can reference these parameters
    // later for derivative calculations using autodiff.
    Map<String, Param> G0;
    G0["H2O"]           = Param("G0[H2O]"          ,  -237181.72);
    G0["H+"]            = Param("G0[H+]"           ,        0.00);
    G0["OH-"]           = Param("G0[OH-]"          ,  -157297.48);
    G0["H2"]            = Param("G0[H2]"           ,    17723.42);
    G0["O2"]            = Param("G0[O2]"           ,    16543.54);
    G0["Na+"]           = Param("G0[Na+]"          ,  -261880.74);
    G0["Cl-"]           = Param("G0[Cl-]"          ,  -131289.74);
    G0["NaCl"]          = Param("G0[NaCl]"         ,  -388735.44);
    G0["HCl"]           = Param("G0[HCl]"          ,  -127235.44);
    G0["NaOH"]          = Param("G0[NaOH]"         ,  -417981.60);
    G0["Ca++"]          = Param("G0[Ca++]"         ,  -552790.08);
    G0["Mg++"]          = Param("G0[Mg++]"         ,  -453984.92);
    G0["CH4"]           = Param("G0[CH4]"          ,   -34451.06);
    G0["CO2"]           = Param("G0[CO2]"          ,  -385974.00);
    G0["HCO3-"]         = Param("G0[HCO3-]"        ,  -586939.89);
    G0["CO3--"]         = Param("G0[CO3--]"        ,  -527983.14);
    G0["CaCl2"]         = Param("G0[CaCl2]"        ,  -811696.00);
    G0["CaCO3"]         = Param("G0[CaCO3]"        , -1099764.40);
    G0["MgCO3"]         = Param("G0[MgCO3]"        ,  -998971.84);
    G0["SiO2"]          = Param("G0[SiO2]"         ,  -833410.96);
    G0["CO2(g)"]        = Param("G0[CO2(g)]"       ,  -394358.74);
    G0["O2(g)"]         = Param("G0[O2(g)]"        ,        0.00);
    G0["H2(g)"]         = Param("G0[H2(g)]"        ,        0.00);
    G0["H2O(g)"]        = Param("G0[H2O(g)]"       ,  -228131.76);
    G0["CH4(g)"]        = Param("G0[CH4(g)]"       ,   -50720.12);
    G0["CO(g)"]         = Param("G0[CO(g)]"        ,  -137168.26);
    G0["NaCl(s)"]       = Param("G0[NaCl(s)]"      ,  -384120.49);
    G0["CaCO3(s)"]      = Param("G0[CaCO3(s)]"     , -1129177.92);
    G0["MgCO3(s)"]      = Param("G0[MgCO3(s)]"     , -1027833.07);
    G0["CaMg(CO3)2(s)"] = Param("G0[CaMg(CO3)2(s)]", -2166307.84);
    G0["SiO2(s)"]       = Param("G0[SiO2(s)]"      ,  -856238.86);

    // Create the Database object from the data above.
    Database db;
    for(auto [key, value] : G0)
        db.addSpecies( Species(key).withStandardGibbsEnergy(value) );

    SECTION("there is only pure water")
    {
        Phases phases(db);
        // phases.add( AqueousPhase(speciate("H O C")) );
        phases.add( AqueousPhase("H2O H+ OH- O2 H2 HCO3- CO2 CO3--" ) );
        phases.add( GaseousPhase("CO2(g) H2O(g) CH4(g)") );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(50.0, "celsius");
        state.setPressure(100.0, "bar");
        state.setSpeciesAmounts(1e-16); // don't let zeros for amounts
        state.setSpeciesAmount("H2O",    55.0,    "mol");
        state.setSpeciesAmount("H+",     1.0e-7,  "mol");
        state.setSpeciesAmount("OH-",    1.0e-7,  "mol");
        state.setSpeciesAmount("O2",     1.0e-31, "mol");
        state.setSpeciesAmount("H2",     2.0e-31, "mol");
        state.setSpeciesAmount("HCO3-",  1.0e-3,  "mol");
        state.setSpeciesAmount("CO2",    1.0e-1,  "mol");
        state.setSpeciesAmount("CO3--",  1.0e-5,  "mol");
        state.setSpeciesAmount("CO2(g)", 1.0,     "mol");
        state.setSpeciesAmount("H2O(g)", 1.0e-3,  "mol");
        state.setSpeciesAmount("CH4(g)", 1.0e-6,  "mol");

        EquilibriumSpecs specs(system);
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.pH();
        specs.openTo("O2"); // this introduces a *p* control variable n[O2] with the amount of O2 in/out
        specs.addParameter(G0["H2O"]); // do this to enable sensitivity derivatives with respect to G0 param of H2O
        specs.addParameter(G0["HCO3-"]); // do this to enable sensitivity derivatives with respect to G0 param of HCO3-

        EquilibriumDims dims(specs);

        const auto Nn = dims.Nn; // number of species in {H2O, H+, OH-, H2, O2, ...}
        const auto Np = dims.Np; // number of control variables in p = (n[O2])
        const auto Nw = dims.Nw; // number of input parameters in w = (T, P, V, pH, G0[H2O], G0[HCO3-])

        VectorXr p = zeros(Np);
        VectorXr n = state.speciesAmounts();
        Params w = specs.params();

        w["T"] = 300.0; // in K
        w["P"] = 1.0e5; // in Pa

        EquilibriumProps eprops(specs);

        // Start recording derivatives of the chemical properties
        // with respect to seeded variables in n, p, w
        eprops.assembleFullJacobianBegin();

        // Start a series of updates call with seeded n[i] variables
        for(auto i = 0; i < Nn; ++i)
        {
            autodiff::seed(n[i]);
            eprops.update(n, p, w);
            autodiff::unseed(n[i]);
        }

        // Start a series of updates call with seeded p[i] variables
        for(auto i = 0; i < Np; ++i)
        {
            autodiff::seed(p[i]);
            eprops.update(n, p, w);
            autodiff::unseed(p[i]);
        }

        // Start a series of updates call with seeded w[i] variables
        for(auto i = 0; i < Nw; ++i)
        {
            autodiff::seed(w[i]);
            eprops.update(n, p, w);
            autodiff::unseed(w[i]);
        }

        // End the recording of derivatives in further calls to EquilibriumProps::update.
        eprops.assembleFullJacobianEnd();

        // Construct the u = u(n, p, w) function that returns all chemical properties in a vector
        auto u = [&](VectorXrConstRef n, VectorXrConstRef p, const Params& w) -> VectorXr
        {
            const auto T = w.get("T");
            const auto P = w.get("P");
            ChemicalProps props(system);
            props.update(T, P, n);
            ArrayStream<real> stream;
            props.serialize(stream);
            return stream.data();
        };

        // Calculate the expected Jacobian matrices of chemical properties u wrt n, p, w
        const auto dudn = jacobian(u, wrt(n), at(n, p, w));
        const auto dudp = jacobian(u, wrt(p), at(n, p, w));
        const auto dudw = jacobian(u, wrt(w), at(n, p, w));

        // Compare to those recorded inside the EquilibriumProps object
        INFO("eprops.dudn()\n" << eprops.dudn());
        INFO("dudn\n" << dudn);
        INFO("eprops.dudp()\n" << eprops.dudp());
        INFO("dudp\n" << dudp);
        INFO("eprops.dudw()\n" << eprops.dudw());
        INFO("dudw\n" << dudw);
        CHECK( eprops.dudn().isApprox(dudn) );
        CHECK( eprops.dudp().isApprox(dudp) );
        CHECK( eprops.dudw().isApprox(dudw) );
    }
}
