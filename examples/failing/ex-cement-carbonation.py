from reaktoro import *

db = ThermoFunDatabase("cemdata18")

params = ActivityModelDebyeHuckelParams()
params.aiondefault = 3.67
params.biondefault = 0.123
params.bneutraldefault = 0.123

solution = AqueousPhase(speciate("H O K Na S Si Ca Mg Al C Cl"))
solution.setActivityModel(ActivityModelDebyeHuckel(params))

minerals = MineralPhases("Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag ")

ss_C3AFS084H  = SolidPhase("C3FS0.84H4.32 C3AFS0.84H4.32")
ss_C3AFS084H.setName("ss_C3AFS084H")

ss_ettringite = SolidPhase("ettringite ettringite30")
ss_ettringite.setName("ss_Ettrignite")

ss_OH_SO4_AFm = SolidPhase("C4AH13 monosulphate12")
ss_OH_SO4_AFm.setName("ss_Monosulfate")

ss_CSHQ = SolidPhase("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH")
ss_CSHQ.setName("ss_CSHQ")

system = ChemicalSystem(db, solution, minerals, ss_C3AFS084H, ss_ettringite, ss_OH_SO4_AFm, ss_CSHQ)

cement_clinker = Material(system)
cement_clinker.add("SiO2" , 20.47, "g")
cement_clinker.add("CaO"  , 65.70, "g")
cement_clinker.add("Al2O3",  4.90, "g")
cement_clinker.add("Fe2O3",  3.20, "g")
cement_clinker.add("K2O"  ,  0.79, "g")
cement_clinker.add("Na2O" ,  0.42, "g")
cement_clinker.add("MgO"  ,  1.80, "g")
cement_clinker.add("SO3"  ,  2.29, "g")
cement_clinker.add("CO2"  ,  0.26, "g")
cement_clinker.add("O2"   ,  0.15, "g")

water = Material(system)
water.add("H2O", 1.0, "g")

calcite = Material(system)
calcite.add("CaCO3", 1, "g")

cement_mix = cement_clinker(84.0, "g") + calcite(16.0, "g") + water(50.0, "g")

state = cement_mix.equilibrate(25.0, "celsius", 1.0, "bar")

assert cement_mix.result().succeeded(), "Equilibrium calculation failed"
