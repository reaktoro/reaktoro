import numpy as np

def stateDictionary(state):
    T = state.temperature()
    P = state.pressure()
    R = 8.3144621 #universalGasConstant
    F = 96485.3329 #faradayConstant 
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
    
    elementAmountInPhase = []
    for i in range(0,system.numPhases()):
        elementAmountInPhase.append(state.elementAmountsInPhase(i))
    
    outputDic ={}
    outputDic['Temperature [K]'] = np.asarray([T]) 
    outputDic['Pressure [Pa]'] = np.asarray([P])
    outputDic['Total element amount [mol]'] = np.asarray(total_element_amounts)
    for i in range(0,system.numPhases()):
        outputDic['Elements amount in '+ system.phase(i).name() + ' [mol]'] = np.asarray(elementAmountInPhase[i]) 
    outputDic['Dual Potential [kJ/mol]'] = y/1000
    outputDic['Species amount [mol]'] = n
    outputDic['Mole Fraction [mol/mol]'] = molar_fractions
    outputDic['Activity coefficient [-]'] = np.exp(lnactivity_coeffs)
    outputDic['Activity [-]'] = np.exp(lnactivities)
    outputDic['Potential [kJ/mol]'] = chemical_potentials/1000
    outputDic['Phase Amount [mol]'] = phase_moles
    outputDic['Stability Index [-]'] = phase_stability_indices
    outputDic['Phase Mass [kg]'] = phase_masses
    outputDic['Phase Volume [m³]'] = phase_volumes
    outputDic['Density [kg/m³]'] = phase_densities
    outputDic['Molar Volume [m³/mol]'] = phase_molar_volumes
    outputDic['Volume Fraction [m³/m³]'] = phase_volume_fractions
    
    return outputDic