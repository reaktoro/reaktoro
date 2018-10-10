import numpy as np

from reaktoro import ChemicalProperty

def convert_reaktoro_state_to_dict(state):
    
    T = state.temperature()
    P = state.pressure() 
    n = state.speciesAmounts()
    y = state.elementDualPotentials()
    z = state.speciesDualPotentials()

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
    
    phase_volume_fractions = phase_volumes/np.sum(phase_volumes)
    
    phase_densities = phase_masses/phase_volumes
    
    phase_stability_indices = state.phaseStabilityIndices()
    
    total_element_amounts = state.elementAmounts()
    
    pH = ChemicalProperty.pH(system)(properties).val
    
    elementAmountInPhase = []
    for i in range(0,system.numPhases()):
        elementAmountInPhase.append(state.elementAmountsInPhase(i))
    
    output ={}
    output['Temperature [K]'] = np.asarray([T]) 
    output['Pressure [Pa]'] = np.asarray([P])
    output['Total element amount [mol]'] = np.asarray(total_element_amounts)
    for i in range(0,system.numPhases()):
        output['Elements amount in '+ system.phase(i).name() + ' [mol]'] = np.asarray(elementAmountInPhase[i]) 
    output['Element Dual Potential [kJ/mol]'] = y/1000
    output['Species amount [mol]'] = n
    output['Mole Fraction [mol/mol]'] = molar_fractions
    output['Activity coefficient [-]'] = np.exp(lnactivity_coeffs)
    output['Activity [-]'] = np.exp(lnactivities)
    output['Potential [kJ/mol]'] = chemical_potentials/1000
    output['Phase Amount [mol]'] = phase_moles
    output['Stability Index [-]'] = phase_stability_indices
    output['Phase Mass [kg]'] = phase_masses
    output['Phase Volume [m³]'] = phase_volumes
    output['Density [kg/m³]'] = phase_densities
    output['Molar Volume [m³/mol]'] = phase_molar_volumes
    output['Volume Fraction [m³/m³]'] = phase_volume_fractions
    output['pH [-]'] = np.asarray([pH])
    
    return output

def convert_dataframe_to_dict(table):

    output = {}
    for elem in table.head():
        output[elem] = table[elem].astype(float).values
    
    return output