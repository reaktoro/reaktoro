import sys
sys.path.append('/home/leal_a/git/Reaktoro/reaktoro/')
from reaktoro.interpreter.ireaktoro import interpret

input = """
ChemicalSystem:
    Database: supcrt98.xml
    AqueousPhase:
        Species: H2O(l) H+ OH- HCO3- CO2(aq) Ca+2 Mg+2 Na+ Cl- O2(aq) H2(aq)
    GaeousPhase:
        Species: H2O(g) CO2(g)
    MineralPhases: Calcite

ReactionSystem:
    MineralReaction Calcite:
        Equation: -1:Calcite -1:H+ 1:Ca+2 1:HCO3-
        SpecificSurfaceArea: 9.8 cm2/g
        Mechanism Acid:
            RateConstant: 10**(-0.30) mol/(m2*s)
            ActivationEnergy: 14.4 kJ/mol
            ActivityPower H+: 1.0
            PressurePower CO2(g): 1.0
        Mechanism Neutral:
            RateConstant: 10**(-5.81) mol/(m2*s)
            ActivationEnergy: 23.5 kJ/mol

Equilibrium State1:
    Temperature: 60 celsius
    Pressure: 150 bar
    Mixture:
        H2O: 1 kg
        NaCl: 1 umol
        HCl: 1 umol
        O2: 1 umol
        CaCO3: 10 mol

Equilibrium State2:
    Temperature: 60 celsius
    Pressure: 150 bar
    Mixture:
        H2O: 1 kg
        NaCl: 1 umol
        HCl: 1 mol
        O2: 1 umol
        CaCO3: 10 mol

EquilibriumPath:
    From: State1
    To: State2
    Plot 1:
        x: t
        y: m[Ca] pH
"""

interpret(input)
