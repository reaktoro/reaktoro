from reaktoro import *

def test_GaseousSpecies_FluidSpecies_LiquidSpecies_import():
    from reaktoro import FluidSpecies
    from reaktoro import GaseousSpecies
    from reaktoro import LiquidSpecies

def test_FluidSpecies_use():
    Tc = 190.6
    Pc = 45.99
    wc = 0.012
    fluidSpecies = FluidSpecies()
    fluidSpecies.setCriticalTemperature(Tc)
    fluidSpecies.setCriticalPressure(Pc)
    fluidSpecies.setAcentricFactor(wc)
    assert fluidSpecies.criticalTemperature() == Tc
    assert fluidSpecies.criticalPressure() == Pc
    assert fluidSpecies.acentricFactor() == wc

def test_GaseousSpecies_use():
    Tc = 190.6
    Pc = 45.99
    wc = 0.012
    gaseousSpecies = GaseousSpecies()
    gaseousSpecies.setCriticalTemperature(Tc)
    gaseousSpecies.setCriticalPressure(Pc)
    gaseousSpecies.setAcentricFactor(wc)
    assert gaseousSpecies.criticalTemperature() == Tc
    assert gaseousSpecies.criticalPressure() == Pc
    assert gaseousSpecies.acentricFactor() == wc

def test_LiquidSpecies_use():
    Tc = 190.6
    Pc = 45.99
    wc = 0.012
    liquidSpecies = LiquidSpecies()
    liquidSpecies.setCriticalTemperature(Tc)
    liquidSpecies.setCriticalPressure(Pc)
    liquidSpecies.setAcentricFactor(wc)
    assert liquidSpecies.criticalTemperature() == Tc
    assert liquidSpecies.criticalPressure() == Pc
    assert liquidSpecies.acentricFactor() == wc
    