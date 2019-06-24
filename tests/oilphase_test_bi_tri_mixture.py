import pytest
from reaktoro import *

"""
All the tests tries to ensure the capability of Reaktoro on solving system that
has liquidlike phases and single mixture
"""
from _pytest.python_api import approx

@pytest.fixture()
def result1():
    return [0.00011203134060024244, 0.0010495435393212702, 0.02169859291970615, 0.9771398322003723, 0.0992084064592162, 0.0982836648316285, 0.773148274101816, 0.029359654607339174]

@pytest.fixture()
def result2():
    return [0.0002272061857777754, 0.0015423224743338398, 0.024962982938293878, 0.9732674884015945, 0.09791743893448875, 0.09665134556058955, 0.7610569547352207, 0.044374260769701024]

@pytest.fixture()
def result3():
    return [0.0003673993311776902, 0.0018805983412758444, 0.024677927357265494, 0.973074074970281, 0.09441698407292762, 0.09306279881112409, 0.7358815435448768, 0.07663867357107167]

@pytest.fixture()
def result4():
    return [0.00022060924973954155, 0.0012639209346541965, 0.02334365876241143, 0.9751718110531948, 0.08858766001037123, 0.0877789125352502, 0.6919739799416967, 0.13165944751268185]


@pytest.fixture()
def result5():
    return [0.0003903551197939781, 0.0018934267624614442, 0.029578541078506923, 0.9681376770392376, 0.08806134279215526, 0.08690816128539626, 0.684193490189778, 0.14083700573267047]

@pytest.fixture()
def result6():
    return [0.0005492314119343933, 0.0023981722832135378, 0.0231738569146855, 0.9738787393901666, 0.09805851780867789, 0.09626163253618183, 0.7662168743847, 0.03946297527044041]

@pytest.fixture()
def result7():
    return [0.000575137181710875, 0.002453852363103085, 0.023287167041705362, 0.9736838434134807, 0.09787280894172562, 0.09605308971934498, 0.7648831873227171, 0.041190914016212214]

@pytest.fixture()
def result8():
    return [0.000575137181710875, 0.002453852363103085, 0.023287167041705362, 0.9736838434134807, 0.09787280894172562, 0.09605308971934498, 0.7648831873227171, 0.041190914016212214]

@pytest.fixture()
def result9():
    return [0.0004446602888606334, 0.0022704656294320175, 0.023200732362584835, 0.9740841417191225, 0.39146506496612105, 0.12321944982374416, 0.4830641207935819, 0.002251364416552807, 0.08613925643978738, 0.09560124167018687, 0.7808936234358999, 0.037365878454125936]

@pytest.fixture()
def result10():
    return [0.0004446602888606334, 0.0022704656294320175, 0.023200732362584835, 0.9740841417191225, 0.39146506496612105, 0.12321944982374416, 0.4830641207935819, 0.002251364416552807, 0.08613925643978738, 0.09560124167018687, 0.7808936234358999, 0.037365878454125936]

@pytest.fixture()
def result11():
    return [0.00031837419639163517, 0.001863801347719783, 0.02290086207222947, 0.9749169623836591, 0.2473321292759114, 0.12504373633776059, 0.6193144040855454, 0.0083097303007828, 0.08183205790091432, 0.09310435483761828, 0.7765921452677544, 0.04847144199371302]




def test_liquid_like_with_CH4(result4):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was CH4 which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CH4(g)
        result[1] = amount of CH4(oil)
    """
    
    db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    """
    Current version of Reaktoro can only work with system that has H2O, it was
    added with amount equal 0
    """ 
    editor.addAqueousPhase(["Methane(aq)", "CO2(aq)",  "H2S(aq)", "H2O(l)" ])
    
    phases = []
    
    phases.append(convertPhase(editor.aqueousPhase(), db))
    
    #Add gas phase
    gas_species = []
    gas_species.append(db.gaseousSpecies("CH4(g)"))
    gas_species.append(db.gaseousSpecies("CO2(g)"))
    gas_species.append(db.gaseousSpecies("H2S(g)"))
    gas_species.append(db.gaseousSpecies("H2O(g)"))
    
    
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelPengRobinson()
    
    phases.append(convertPhase(gas_phase, db))
    
    #Add oil phase
    oil_species = [] #OilSpecies_vector()
    oil_species.append(OilSpecies(db.gaseousSpecies("CH4(oil)")))
    oil_species.append(OilSpecies(db.gaseousSpecies("CO2(oil)")))
    oil_species.append(OilSpecies(db.gaseousSpecies("H2S(oil)")))
    oil_species.append(OilSpecies(db.gaseousSpecies("H2O(oil)")))    
    
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = OilPhase(oil_mixture)
    oil_phase.setChemicalModelPengRobinson()
    
    phases.append(convertPhase(oil_phase, db))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(338.75 , "K")
    problem.setPressure(8.43 , "MPa")
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
    
    molar_fractions = []
    
    amounts = state.speciesAmounts()
    properties = state.properties();
    
    
    for i in range(system.numSpecies()):
        if (amounts[i] > 1.e-10):
             molar_fractions.append(properties.moleFractions().val[i]) 
    
    print(state)
    print(molar_fractions)
    
    assert 1 == 1 ## molar_fractions == approx(result4)

