import numpy as np

from reaktoro import ChemicalProperty


def convert_reaktoro_state_to_dict(state, exclude):

    T = state.temperature()
    P = state.pressure()
    n = state.speciesAmounts()
    y = state.elementDualPotentials()
    z = state.speciesDualPotentials()
    b = state.elementAmounts()
    R = 8.314462618  # universal gas constant in J/(mol*K)

    properties = state.properties()

    system = state.system()

    molar_fractions = properties.moleFractions().val

    lnactivity_coeffs = properties.lnActivityCoefficients().val

    lnactivities = properties.lnActivities().val

    chemical_potentials = properties.chemicalPotentials().val

    phase_moles = properties.phaseAmounts().val

    phase_masses = properties.phaseMasses().val

    phase_molar_volumes = properties.phaseMolarVolumes().val

    phase_volumes = properties.phaseVolumes().val

    phase_densities = properties.phaseDensities().val

    phase_stability_indices = state.phaseStabilityIndices()

    pH = ChemicalProperty.pH(system)(properties).val

    elementAmountsInPhase = []
    for i in range(0, system.numPhases()):
        elementAmountsInPhase.append(state.elementAmountsInPhase(i))

    output = {}
    output["Temperature [K]"] = np.asarray([T])
    output["Pressure [Pa]"] = np.asarray([P])
    output["Element amounts [mol]"] = np.asarray(b)
    output["Gibbs energy [-]"] = np.asarray([n.dot(chemical_potentials) / (R*T)])  # normalized by RT
    # output["Gibbs energy (dual) [-]"] = np.asarray([b.dot(y) / (R*T)])  # normalized by RT           # the current optimization algorithm produces slightly different y values across operating systems
    for i in range(0, system.numPhases()):
        output["Element amounts in " + system.phase(i).name() + " [mol]"] = \
            np.asarray(elementAmountsInPhase[i])
    # output["Element Dual Potential [kJ/mol]"] = y / 1000
    output["Species amounts [mol]"] = n
    output["Mole fractions [mol/mol]"] = molar_fractions
    output["Activity coefficients [-]"] = np.exp(lnactivity_coeffs)
    output["Activities [-]"] = np.exp(lnactivities)
    # output["Chemical potentials [kJ/mol]"] = chemical_potentials / 1000
    output["Phase amounts [mol]"] = phase_moles
    output["Stability indices [-]"] = phase_stability_indices
    output["Phase masses [kg]"] = phase_masses
    output["Phase volumes [m³]"] = phase_volumes
    output["Phase densities [kg/m³]"] = phase_densities
    output["Phase molar volumes [m³/mol]"] = phase_molar_volumes
    output["pH [-]"] = np.asarray([pH])

    if exclude is not None:
        output = { key: val for key, val in output.items() if key not in exclude }

    return output


def convert_table_to_dict(table):

    output = {}
    for elem in table.head():
        output[elem] = table[elem].astype(float).values

    return output
