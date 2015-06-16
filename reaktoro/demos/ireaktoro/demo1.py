from reaktoro.interpreter.ireaktoro import interpret

input = """
ChemicalSystem:
    Database: supcrt98.xml
    AqueousPhase:
        Species: H2O(l) H+ OH- Na+ Cl- Ca+2 Mg+2 CO2(aq) HCO3- O2(aq) H2(aq)
    MineralPhases: Calcite Dolomite Quartz

Equilibrium StateIC:
    Temperature: 100 celsius
    Pressure: 300 bar
    Mix: |
        H2O 1 kg
        NaCl 1 mol
        CaCO3 200 g
        MgCO3 2 g
        SiO2 10 g
        O2 1 umol
"""

interpret(input)
