# API Reference {#mainpage}

Welcome to %Reaktoro's C++ API reference. You can also use this API guide for
your Python applications.

Find below a list of important methods and classes in %Reaktoro that you might
want to read about.

### Thermochemical databases

- @ref Reaktoro::Database
- @ref Reaktoro::PhreeqcDatabase
- @ref Reaktoro::ThermoFunDatabase
- @ref Reaktoro::SupcrtDatabase (*use Database instead*)
- @ref Reaktoro::NasaDatabase (*in development*)

### Chemical system constituents

- @ref Reaktoro::Element
- @ref Reaktoro::Species
- @ref Reaktoro::Phase

### Chemical system definition

- @ref Reaktoro::AqueousPhase
- @ref Reaktoro::GaseousPhase
- @ref Reaktoro::LiquidPhase
- @ref Reaktoro::MineralPhase
- @ref Reaktoro::MineralPhases
- @ref Reaktoro::SolidPhase
- @ref Reaktoro::Phases
- @ref Reaktoro::ChemicalSystem

### Chemical state and thermochemical properties of a chemical system

- @ref Reaktoro::ChemicalState
- @ref Reaktoro::ChemicalPropsPhase
- @ref Reaktoro::ChemicalProps
- @ref Reaktoro::ThermoPropsPhase
- @ref Reaktoro::ThermoProps
- @ref Reaktoro::AqueousProps

### Chemical equilibrium problem definition and calculations

- @ref Reaktoro::EquilibriumSpecs
- @ref Reaktoro::EquilibriumRestrictions
- @ref Reaktoro::EquilibriumConditions
- @ref Reaktoro::EquilibriumOptions
- @ref Reaktoro::EquilibriumSolver
- @ref Reaktoro::EquilibriumResult
- @ref Reaktoro::EquilibriumSensitivity
- @ref Reaktoro::EquilibriumPredictor

### Activity models for aqueous phases

- @ref Reaktoro::ActivityModelIdealAqueous
- @ref Reaktoro::ActivityModelDebyeHuckel
- @ref Reaktoro::ActivityModelPitzerHMW
- @ref Reaktoro::ActivityModelHKF
- @ref Reaktoro::ActivityModelDrummond
- @ref Reaktoro::ActivityModelDuanSun
- @ref Reaktoro::ActivityModelRumpf
- @ref Reaktoro::ActivityModelSetschenow

### Activity models for fluid phases (gaseous or liquid)

- @ref Reaktoro::ActivityModelIdealGas
- @ref Reaktoro::ActivityModelIdealSolution
- @ref Reaktoro::ActivityModelVanDerWaals
- @ref Reaktoro::ActivityModelRedlichKwong
- @ref Reaktoro::ActivityModelSoaveRedlichKwong
- @ref Reaktoro::ActivityModelPengRobinson
- @ref Reaktoro::ActivityModelSpycherPruessEnnis
- @ref Reaktoro::ActivityModelSpycherReed

### Activity models for solid solutions

- @ref Reaktoro::ActivityModelIdealSolution
- @ref Reaktoro::ActivityModelRedlichKister
- @ref Reaktoro::ActivityModelVanLaar (*needs adaptation to v2.0*)

### Standard thermodynamic properties and models

- @ref Reaktoro::StandardThermoProps
- @ref Reaktoro::StandardThermoModelConstant
- @ref Reaktoro::StandardThermoModelHKF
- @ref Reaktoro::StandardThermoModelHollandPowell
- @ref Reaktoro::StandardThermoModelMaierKelley
- @ref Reaktoro::StandardThermoModelMineralHKF
- @ref Reaktoro::StandardThermoModelWaterHKF
- @ref Reaktoro::StandardVolumeModel
- @ref Reaktoro::StandardVolumeModelConstant

### Reaction thermodynamic properties and models

- @ref Reaktoro::ReactionThermoProps
- @ref Reaktoro::ReactionThermoModelConstLgK
- @ref Reaktoro::ReactionThermoModelGemsLgK
- @ref Reaktoro::ReactionThermoModelPhreeqcLgK
- @ref Reaktoro::ReactionThermoModelPressureCorrection
- @ref Reaktoro::ReactionThermoModelVantHoff

### Water thermodynamic and electrostatic properties

- @ref Reaktoro::WaterElectroProps
- @ref Reaktoro::WaterThermoProps
- @ref Reaktoro::WaterHelmholtzState (to be renamed `WaterHelmholtzProps` in a future release)
- @ref Reaktoro::waterDensityHGK
- @ref Reaktoro::waterDensityWagnerPruss
- @ref Reaktoro::waterLiquidDensityHGK
- @ref Reaktoro::waterLiquidDensityWagnerPruss
- @ref Reaktoro::waterVaporDensityHGK
- @ref Reaktoro::waterVaporDensityWagnerPruss
- @ref Reaktoro::waterPressureHGK
- @ref Reaktoro::waterPressureWagnerPruss
- @ref Reaktoro::waterSaturatedPressureWagnerPruss
- @ref Reaktoro::waterSaturatedLiquidDensityWagnerPruss
- @ref Reaktoro::waterSaturatedVapourDensityWagnerPruss
- @ref Reaktoro::waterElectroPropsJohnsonNorton
- @ref Reaktoro::waterThermoPropsHGK
- @ref Reaktoro::waterThermoPropsWagnerPruss
- @ref Reaktoro::waterHelmholtzStateHGK (to be renamed `waterHelmholtzPropsHGK` in a future release)
- @ref Reaktoro::waterHelmholtzStateWagnerPruss (to be renamed `waterHelmholtzPropsWagnerPruss` in a future release)

### Worth checking classes and methods

- @ref Reaktoro::ElementList
- @ref Reaktoro::PhaseList
- @ref Reaktoro::SpeciesList
- @ref Reaktoro::FormationReaction
- @ref Reaktoro::ChemicalFormula
- @ref Reaktoro::Elements
- @ref Reaktoro::CriticalProps
- @ref Reaktoro::DissociationReactions
- @ref Reaktoro::ReactionEquation
- @ref Reaktoro::AggregateState
- @ref Reaktoro::ElementalComposition
- @ref Reaktoro::Reaction
- @ref Reaktoro::StateOfMatter

### For all other needs

Check the namespace Reaktoro for all available classes and methods.
