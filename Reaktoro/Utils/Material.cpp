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

#include "Material.hpp"

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>

namespace Reaktoro {

Material::Material(ChemicalSystem const& system)
: m_system(system)
{
}

auto Material::addSpeciesAmount(StringOrIndex const& species, double amount) -> void
{
    const auto ispecies = detail::resolveSpeciesIndex(m_system, species);
    errorif(ispecies >= m_system.species().size(),
        "Could not add species with name or index `", stringfy(species), "` "
        "to your custom material because it does not exist "
        "in the underlying chemical system of the material.");
    const auto idx = indexfn(m_species, RKT_LAMBDA(x, x.first == ispecies)); // check if m_species already contains ispecies
    if(idx < m_species.size())
        m_species[idx].second += amount; // accumulate on existing entry for ispecies
    else m_species.emplace_back(ispecies, amount); // create a new entry for ispecies
}

auto Material::addSpeciesAmount(StringOrIndex const& species, double amount, Chars unit) -> void
{
    amount = units::convert(amount, unit, "mol");
    addSpeciesAmount(species, amount);
}

auto Material::addSpeciesMass(StringOrIndex const& species, double mass, Chars unit) -> void
{
    const auto ispecies = detail::resolveSpeciesIndex(m_system, species);
    errorif(ispecies >= m_system.species().size(),
        "Could not add species with name or index `", stringfy(species), "` "
        "to your custom material because it does not exist "
        "in the underlying chemical system of the material.");
    mass = units::convert(mass, unit, "kg");
    const auto amount = mass / m_system.species(ispecies).molarMass(); // kg / (kg/mol) = mol
    addSpeciesAmount(ispecies, amount);
}

auto Material::addSubstanceAmount(ChemicalFormula const& substance, double amount) -> void
{
    auto const& system_elements = m_system.elements();
    for(auto const& [element, coefficient] : substance.elements())
        errorifnot(system_elements.find(element) < system_elements.size(),
            "While adding substance `", substance, "` to this material, I found that "
            "element `", element, "` does not exist in the underlying chemical system. "
            "Please consider this element in the definition of your chemical system.");
    const auto idx = indexfn(m_substances, RKT_LAMBDA(x, x.first == substance)); // check if m_substances already contains substance
    if(idx < m_substances.size())
        m_substances[idx].second += amount; // accumulate on existing entry for substance
    else m_substances.emplace_back(substance, amount); // create a new entry for substance
}

auto Material::addSubstanceAmount(ChemicalFormula const& substance, double amount, Chars unit) -> void
{
    amount = units::convert(amount, unit, "mol");
    addSubstanceAmount(substance, amount);
}

auto Material::addSubstanceMass(ChemicalFormula const& substance, double mass, Chars unit) -> void
{
    mass = units::convert(mass, unit, "kg");
    const auto amount = mass / substance.molarMass(); // kg / (kg/mol) = mol
    addSubstanceAmount(substance, amount);
}

auto Material::addMaterialAmount(Material const& material, double amount) -> void
{
    errorif(!identical(m_system.elements(), material.m_system.elements()),
        "While adding another material into your custom material, "
        "I found that they don't have the same underlying chemical system.")

    const auto scale_factor = amount / material.amount();
    for(const auto& [substance, substance_amount] : material.m_substances)
        addSubstanceAmount(substance, substance_amount * scale_factor);
    for(const auto& [ispecies, species_amount] : material.m_species)
        addSpeciesAmount(ispecies, species_amount * scale_factor);
}

auto Material::addMaterialAmount(Material const& material, double amount, Chars unit) -> void
{
    amount = units::convert(amount, unit, "mol");
    addMaterialAmount(material, amount);
}

auto Material::addMaterialMass(Material const& material, double mass, Chars unit) -> void
{
    mass = units::convert(mass, unit, "kg");
    const auto amount = mass / material.molarMass(); // kg / (kg/mol) = mol
    addMaterialAmount(material, amount);
}

auto Material::add(String const& substance, double value, Chars unit) -> void
{
    const auto ispecies = m_system.species().findWithName(substance);
    if(ispecies < m_system.species().size())
    {
        if(units::convertible(unit, "mol"))
            addSpeciesAmount(ispecies, value, unit);
        else if(units::convertible(unit, "kg"))
            addSpeciesMass(ispecies, value, unit);
        else errorif(true, "While adding species `", substance, "` to your custom material, "
            "I encountered your given unit `", unit, "` that is neither convertible to mol nor kg.");
    }
    else
    {
        if(units::convertible(unit, "mol"))
            addSubstanceAmount(substance, value, unit);
        else if(units::convertible(unit, "kg"))
            addSubstanceMass(substance, value, unit);
        else errorif(true, "While adding substance `", substance, "` to your custom material, "
            "I encountered your given unit `", unit, "` that is neither convertible to mol nor kg.");
    }
}

auto Material::add(Material const& material, double value, Chars unit) -> void
{
    if(units::convertible(unit, "mol"))
        addMaterialAmount(material, value, unit);
    else if(units::convertible(unit, "kg"))
        addMaterialMass(material, value, unit);
    else errorif(true, "While adding a material to your custom material, "
        "I encountered your given unit `", unit, "` that is neither convertible to mol nor kg.")
}

auto Material::scaleAmount(double value, Chars unit) -> void
{
    value = units::convert(value, unit, "mol");
    const auto scale_factor = value / amount();
    for(auto& [substance, substance_amount] : m_substances)
        substance_amount *= scale_factor;
    for(auto& [ispecies, species_amount] : m_species)
        species_amount *= scale_factor;
}

auto Material::scaleMass(double value, Chars unit) -> void
{
    value = units::convert(value, unit, "kg");
    const auto amount = value / molarMass(); // kg / (kg/mol) = mol
    scaleAmount(amount, "mol");
}

auto Material::scale(double value, Chars unit) -> void
{
    if(units::convertible(unit, "mol"))
        scaleAmount(value, unit);
    else if(units::convertible(unit, "kg"))
        scaleMass(value, unit);
    else errorif(true, "While scaling the amount/mass of your custom material, "
        "I encountered your given unit `", unit, "` that is neither convertible to mol nor kg.")
}

auto Material::with(double value, Chars unit) const -> Material
{
    Material copy(*this);
    copy.scale(value, unit);
    return copy;
}

auto Material::system() const -> const ChemicalSystem&
{
    return m_system;
}

auto Material::substances() const -> const Pairs<ChemicalFormula, double>&
{
    return m_substances;
}

auto Material::species() const -> const Pairs<Index, double>&
{
    return m_species;
}

auto Material::componentAmounts() const -> ArrayXd
{
    const auto E = m_system.elements().size();
    ArrayXd res = ArrayXd::Zero(E + 1);
    res << elementAmounts(), charge();
    return res;
}

auto Material::elementAmounts() const -> ArrayXd
{
    const auto E = m_system.elements().size();

    const auto idx = [&](auto const& symbol)
    {
        auto ielement = m_system.elements().find(symbol);
        errorifnot(ielement < E, "The Material object contains substances with element `", symbol, "`, which does not exist in the chemical system.");
        return ielement;
    };

    ArrayXd res = ArrayXd::Zero(E);
    for(const auto& [substance, substance_amount] : m_substances) // iterate over all substances in the material
        for(const auto& [symbol, element_coeff] : substance.elements()) // iterate over all elements in current substance
            res[idx(symbol)] += substance_amount * element_coeff;
    for(const auto& [ispecies, species_amount] : m_species) // iterate over all species in the material
        for(const auto& [element, element_coeff] : m_system.species(ispecies).elements()) // iterate over all elements in current species
            res[idx(element.symbol())] += species_amount * element_coeff;

    return res;
}

auto Material::charge() const -> double
{
    double res = 0.0;
    for(const auto& [substance, substance_amount] : m_substances)
        res += substance_amount * substance.charge();
    for(const auto& [ispecies, species_amount] : m_species)
        res += species_amount * m_system.species(ispecies).charge();
    return res;
}

auto Material::amount() const -> double
{
    double sum = 0.0;
    for(const auto& [substance, substance_amount] : m_substances)
        sum += substance_amount;
    for(const auto& [ispecies, species_amount] : m_species)
        sum += species_amount;
    return sum;
}

auto Material::mass() const -> double
{
    double sum = 0.0;
    for(const auto& [substance, substance_amount] : m_substances)
        sum += substance_amount * substance.molarMass(); // mol * kg/mol = kg
    for(const auto& [ispecies, species_amount] : m_species)
        sum += species_amount * m_system.species(ispecies).molarMass(); // mol * kg/mol = kg
    return sum;
}

auto Material::molarMass() const -> double
{
    return mass() / amount();
}

auto Material::equilibrate() -> ChemicalState
{
    return equilibrate(EquilibriumOptions());
}

auto Material::equilibrate(const EquilibriumOptions& options) -> ChemicalState
{
    return equilibrate(EquilibriumRestrictions(m_system), options);
}

auto Material::equilibrate(const EquilibriumRestrictions& restrictions) -> ChemicalState
{
    return equilibrate(restrictions, EquilibriumOptions());
}

auto Material::equilibrate(const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> ChemicalState
{
    return equilibrate(298.15, "K", 1.0e5, "Pa", restrictions, options);
}

auto Material::equilibrate(double T, Chars unitT, double P, Chars unitP) -> ChemicalState
{
    return equilibrate(T, unitT, P, unitP, EquilibriumOptions());
}

auto Material::equilibrate(double T, Chars unitT, double P, Chars unitP, const EquilibriumOptions& options) -> ChemicalState
{
    return equilibrate(T, unitT, P, unitP, EquilibriumRestrictions(m_system), options);
}

auto Material::equilibrate(double T, Chars unitT, double P, Chars unitP, const EquilibriumRestrictions& restrictions) -> ChemicalState
{
    return equilibrate(T, unitT, P, unitP, restrictions, EquilibriumOptions());
}

auto Material::equilibrate(double T, Chars unitT, double P, Chars unitP, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> ChemicalState
{
    // Convert temperature and pressure to SI units
    T = units::convert(T, unitT, "K");
    P = units::convert(P, unitP, "Pa");

    // Construct a suitable initial chemical state for the equilibrium calculation
    ChemicalState state = initialState(T, P);

    // Get the amounts of elements and charge in the material
    const auto b0 = componentAmounts();

    // Finally, perform the equilibrium calculation enforcing element/charge amounts in b0
    m_result = Reaktoro::equilibrate(state, restrictions, options, b0);

    return state;
}

auto Material::result() const -> const EquilibriumResult&
{
    return m_result;
}

auto Material::initialState(double T, double P) const -> ChemicalState
{
    ChemicalState state(m_system);
    state.setTemperature(T);
    state.setPressure(P);

    // Skip the rest if no substances or species were added into the material!
    if(m_substances.empty() && m_species.empty())
        return state;

    // Use the standard chemical potentials of the species at (T,P) to
    // determine which substance amounts are used as initial amounts of
    // chemical species in the initial chemical state. For example, assume 55
    // mols of H2O has been given as substance in the material. If the system
    // has species H2O(aq) and H2O(g), we must decide which one gets
    // initialized with 55 mols. We should give priority to the chemical
    // species with least standard Gibbs energy at (T,P) (most stable).
    for(const auto& [substance, substance_amount] : m_substances)
    {
        // Find the species with same formula as this substance. The species
        // with least standard Gibbs energy at (T, P) will be chosen to receive
        // the substance amount as its initial amount in the chemical state.
        Indices ispecies_canditates;
        for(auto const& [i, species] : enumerate(m_system.species()))
            if(species.formula().equivalent(substance))
                ispecies_canditates.push_back(i);

        // Skip to the next substance if this substance has no species candidate (with same formula)
        if(ispecies_canditates.empty())
            continue;

        // Compute the standard chemical potentials of the candidate species
        const auto G0s = vectorize(ispecies_canditates,
            RKT_LAMBDA(i, m_system.species(i).standardThermoProps(T, P).G0));

        // Find the index corresponding to smallest G0
        const auto iminG0 = std::min_element(G0s.begin(), G0s.end()) - G0s.begin();

        // We've now found the species among the candidates with least standard Gibbs energy.
        const auto ichosen_species = ispecies_canditates[iminG0];

        // Initialize the amount of this species with the given substance amount
        state.setSpeciesAmount(ichosen_species, substance_amount, "mol");
    }

    // Add the constribution of the species in the initial chemical state
    for(const auto& [ispecies, species_amount] : m_species)
        state.add(ispecies, species_amount, "mol"); // use add to accumulate in case there were substances above with same formula as this species

    return state;
}

auto Material::operator()(double value, Chars unit) const -> Material
{
    return with(value, unit);
}

auto operator+(const Material& l, const Material& r) -> Material
{
    Material res(l);
    res.addMaterialAmount(r, r.amount());
    return res;
}

auto operator<<(std::ostream& out, const Material& material) -> std::ostream&
{
    // const auto species = props.system().species();
    // const auto elements = props.system().elements();
    // const auto b   = props.elementAmounts();
    // const auto n   = props.speciesAmounts();
    // const auto x   = props.speciesMoleFractions();
    // const auto lng = props.speciesActivityCoefficientsLn();
    // const auto lna = props.speciesActivitiesLn();
    // const auto mu  = props.speciesChemicalPotentials();
    // const auto G0  = props.speciesStandardGibbsEnergies();
    // const auto H0  = props.speciesStandardEnthalpies();
    // const auto V0  = props.speciesStandardVolumes();
    // const auto S0  = props.speciesStandardEntropies();
    // const auto U0  = props.speciesStandardInternalEnergies();
    // const auto A0  = props.speciesStandardHelmholtzEnergies();
    // const auto Cp0 = props.speciesStandardHeatCapacitiesConstP();
    // const auto Cv0 = props.speciesStandardHeatCapacitiesConstV();

    // Table table;
    // table.add_row({ "Property", "Value", "Unit" });
    // table.add_row({ "Temperature", str(props.temperature()), "K" });
    // table.add_row({ "Pressure", str(props.pressure()), "Pa" });
    // table.add_row({ "Volume", str(props.volume()), "m3" });
    // table.add_row({ "Gibbs Energy", str(props.gibbsEnergy()), "J" });
    // table.add_row({ "Enthalpy", str(props.enthalpy()), "J" });
    // table.add_row({ "Entropy", str(props.entropy()), "J/K" });
    // table.add_row({ "Internal Energy", str(props.internalEnergy()), "J" });
    // table.add_row({ "Helmholtz Energy", str(props.helmholtzEnergy()), "J" });
    // table.add_row({ "Charge", str(props.charge()), "mol" });

    // table.add_row({ "Element Amount:" }); for(auto i = 0; i < b.size(); ++i) table.add_row({ ":: " + elements[i].symbol(), str(b[i]), "mol" });
    // table.add_row({ "Species Amount:" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(n[i]), "mol" });
    // table.add_row({ "Mole Fraction:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(x[i]), "mol/mol" });
    // table.add_row({ "Activity Coefficient:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(exp(lng[i])), "-" });
    // table.add_row({ "Activity:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(exp(lna[i])), "-" });
    // table.add_row({ "lg(Activity):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(lna[i]/ln10), "-" });
    // table.add_row({ "ln(Activity):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(lna[i]), "-" });
    // table.add_row({ "Chemical Potential:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(mu[i]), "J/mol" });
    // table.add_row({ "Standard Volume:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(V0[i]), "m3/mol" });
    // table.add_row({ "Standard Gibbs Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(G0[i]), "J/mol" });
    // table.add_row({ "Standard Enthalpy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(H0[i]), "J/mol" });
    // table.add_row({ "Standard Entropy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(S0[i]), "J/(mol*K)" });
    // table.add_row({ "Standard Internal Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(U0[i]), "J/mol" });
    // table.add_row({ "Standard Helmholtz Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(A0[i]), "J/mol" });
    // table.add_row({ "Standard Heat Capacity (constant P):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(Cp0[i]), "J/(mol*K)" });
    // table.add_row({ "Standard Heat Capacity (constant V):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].name(), str(Cv0[i]), "J/(mol*K)" });

    // auto i = 0;
    // for(auto& row : table)
    // {
    //     if(i >= 2)  // apply from the third row
    //         table[i]
    //             .format()
    //             .border_top("")
    //             .column_separator("")
    //             .corner_top_left("")
    //             .corner_top_right("");
    //     i += 1;
    // }

    // table.row(0).format().font_style({FontStyle::bold});  // Bold face for header
    // table.column(1).format().font_align(FontAlign::right); // Value column with right alignment
    // table.column(2).format().font_align(FontAlign::right); // Unit column with right alignment

    // out << table;
    return out;
}

} // namespace Reaktoro
