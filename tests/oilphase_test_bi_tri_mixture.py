import pytest
from reaktoro import *

"""
All the tests tries to ensure the capability of Reaktoro on solving system that
has liquidlike phases and single mixture
"""
from _pytest.python_api import approx

@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (380.35, 75.6  , [5.56315e-05,
                          0.000521173,
                          0.0107749,
                          0.48522,
                          0.0499444,
                          0.0494788,
                          0.389225,
                          0.0147805,
                          0.000000,
                          0.000000,
                          0.000000,
                          0.000000] ),
        (380.35, 122.7  , [0.000111446,
                          0.000756515,
                          0.0122444,
                          0.477391,
                          0.0498886,
                          0.0492435,
                          0.387756,
                          0.0226085,
                          0.000000,
                          0.000000,
                          0.000000,
                          0.000000] ),
        (310.95, 130.0  , [0.000270695,
                          0.00118197,
                          0.0114215,
                          0.479987,
                          0.000000,
                          0.000000,
                          0.000000,
                          0.000000,
                          0.0497293,
                          0.048818,
                          0.388579,
                          0.0200132] ),
        (310.95, 62.6  , [0.000219954,
                          0.0011231,
                          0.0114764,
                          0.481836,
                          0.00801333,
                          0.00252232,
                          0.00988838,
                          4.60857e-05,
                          0.0417667,
                          0.0463546,
                          0.378635,
                          0.0181178] ),
    ],
    ids=[
        "temperature equal 380.35 K and 75.6 bar",
        "temperature equal 380.35 K and 122.7 bar",
        "temperature equal 310.95 K and 130.0 bar",
        "temperature equal 310.95 K and 62.6 bar",
    ],
)
def test_equilibrium_CH4_H2S_CO2_H2O(temperature, pressure, result):
    """
    This test checks the capability of solving a system that has CH4, H2S,
    CO2, H2O with
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CH4(aq)
        result[0] = amount of H2S(aq)
        result[1] = amount of CO2(aq)
        result[2] = amount of H2O(aq)
        result[3] = amount of CH4(g)
        result[4] = amount of H2S(g)
        result[5] = amount of CO2(g)
        result[6] = amount of H2O(g)
        result[7] = amount of CH4(oil)
        result[8] = amount of H2S(oil)
        result[9] = amount of CO2(oil)
        result[10] = amount of H2O(oil)
    """
    
    db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    """
    Current version of Reaktoro can only work with system that has H2O, it was
    added with amount equal 0
    """ 
    editor.addAqueousPhase(["Methane(aq)", "CO2(aq)",  "H2S(aq)", "H2O(l)" ])
    
    phases = []
    
    phases.append(editor.convertAqueousPhase(editor.aqueousPhase()))
    
    #Add gas phase
    gas_species = []
    gas_species.append(db.gaseousSpecies("CH4(g)"))
    gas_species.append(db.gaseousSpecies("CO2(g)"))
    gas_species.append(db.gaseousSpecies("H2S(g)"))
    gas_species.append(db.gaseousSpecies("H2O(g)"))
    
    
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelPengRobinson()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add oil phase
    oil_species = [] #HydrocarbonSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("H2O(oil)")))    
    
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = LiquidPhase(oil_mixture)
    oil_phase.setChemicalModelPengRobinson()
    
    phases.append(editor.convertLiquidPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("H2O(g)", 0.50, "mol")
    problem.add("CO2(g)", 0.05, "mol")
    problem.add("H2S(g)", 0.40, "mol")
    problem.add("CH4(g)", 0.05, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions();
    options.hessian = GibbsHessian.Exact;
    options.nonlinear.max_iterations = 100;
    options.optimum.max_iterations = 200;
    options.optimum.ipnewton.step = StepMode.Conservative;
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    res_final = solver.solve(state, 
                             problem, 
                             phaseIdentificationMethod.Gibbs_residual_based, 
                             0)
    

    assert (state.speciesAmount("Methane(aq)"),
            state.speciesAmount("CO2(aq)"),
            state.speciesAmount("H2S(aq)"),
            state.speciesAmount("H2O(l)"),
            state.speciesAmount("CH4(g)"),
            state.speciesAmount("CO2(g)"),
            state.speciesAmount("H2S(g)"),
            state.speciesAmount("H2O(g)"),
            state.speciesAmount("CH4(oil)"),
            state.speciesAmount("CO2(oil)"),
            state.speciesAmount("H2S(oil)"),
            state.speciesAmount("H2O(oil)")) == approx(result, abs=1e-6)

