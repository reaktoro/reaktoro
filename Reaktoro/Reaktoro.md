# API Reference {#index}

Welcome to %Reaktoro's C++ API reference. You can also use this API guide for
your Python applications.

Find below a list of important methods and classes in %Reaktoro that you might
want to read about.

## Thermochemical databases

- @ref Reaktoro::Database
- @ref Reaktoro::PhreeqcDatabase
- @ref Reaktoro::ThermoFunDatabase
- @ref Reaktoro::SupcrtDatabase
- @ref Reaktoro::NasaDatabase

## Chemical system constituents

- @ref Reaktoro::Element
- @ref Reaktoro::Species
- @ref Reaktoro::Phase
- @ref Reaktoro::Reaction
- @ref Reaktoro::Surface

## Chemical system definition

- @ref Reaktoro::ChemicalSystem
- @ref Reaktoro::Phases
- @ref Reaktoro::Reactions
- @ref Reaktoro::Surfaces
- @ref Reaktoro::AqueousPhase
- @ref Reaktoro::CondensedPhase
- @ref Reaktoro::CondensedPhases
- @ref Reaktoro::GaseousPhase
- @ref Reaktoro::IonExchangePhase
- @ref Reaktoro::LiquidPhase
- @ref Reaktoro::MineralPhase
- @ref Reaktoro::MineralPhases
- @ref Reaktoro::SolidPhase
- @ref Reaktoro::GeneralReaction
- @ref Reaktoro::MineralReaction
- @ref Reaktoro::GeneralSurface
- @ref Reaktoro::MineralSurface

## Chemical state and thermochemical properties of a chemical system

- @ref Reaktoro::ChemicalState
- @ref Reaktoro::ChemicalPropsPhase
- @ref Reaktoro::ChemicalProps
- @ref Reaktoro::ThermoPropsPhase
- @ref Reaktoro::ThermoProps
- @ref Reaktoro::AqueousProps

## Chemical equilibrium problem definition and calculations

- @ref Reaktoro::EquilibriumSpecs
- @ref Reaktoro::EquilibriumRestrictions
- @ref Reaktoro::EquilibriumConditions
- @ref Reaktoro::EquilibriumOptions
- @ref Reaktoro::EquilibriumSolver
- @ref Reaktoro::EquilibriumResult
- @ref Reaktoro::EquilibriumSensitivity
- @ref Reaktoro::EquilibriumPredictor

## Chemical kinetics problem definition and calculations

- @ref Reaktoro::KineticsSolver
- @ref Reaktoro::KineticsOptions
- @ref Reaktoro::KineticsResult
- @ref Reaktoro::KineticsSensitivity

## Machine learning accelerated chemical equilibrium calculations

- @ref Reaktoro::SmartEquilibriumSolver
- @ref Reaktoro::SmartEquilibriumOptions
- @ref Reaktoro::SmartEquilibriumResult

## Machine learning accelerated chemical kinetics calculations

- @ref Reaktoro::SmartKineticsSolver
- @ref Reaktoro::SmartKineticsOptions
- @ref Reaktoro::SmartKineticsResult

## Activity models for aqueous phases

- @ref Reaktoro::ActivityModelIdealAqueous
- @ref Reaktoro::ActivityModelDavies
- @ref Reaktoro::ActivityModelDebyeHuckel
- @ref Reaktoro::ActivityModelPitzerHMW
- @ref Reaktoro::ActivityModelHKF
- @ref Reaktoro::ActivityModelDrummond
- @ref Reaktoro::ActivityModelDuanSun
- @ref Reaktoro::ActivityModelRumpf
- @ref Reaktoro::ActivityModelSetschenow

## Activity models for fluid phases (gaseous or liquid)

- @ref Reaktoro::ActivityModelIdealGas
- @ref Reaktoro::ActivityModelIdealSolution
- @ref Reaktoro::ActivityModelVanDerWaals
- @ref Reaktoro::ActivityModelRedlichKwong
- @ref Reaktoro::ActivityModelSoaveRedlichKwong
- @ref Reaktoro::ActivityModelPengRobinson
- @ref Reaktoro::ActivityModelSpycherPruessEnnis
- @ref Reaktoro::ActivityModelSpycherReed

## Activity models for solid solutions

- @ref Reaktoro::ActivityModelIdealSolution
- @ref Reaktoro::ActivityModelRedlichKister
- @ref Reaktoro::ActivityModelVanLaar (*needs adaptation to v2.0*)

## Standard thermodynamic properties and models

- @ref Reaktoro::StandardThermoProps
- @ref Reaktoro::StandardThermoModelConstant
- @ref Reaktoro::StandardThermoModelHKF
- @ref Reaktoro::StandardThermoModelHollandPowell
- @ref Reaktoro::StandardThermoModelMaierKelley
- @ref Reaktoro::StandardThermoModelMineralHKF
- @ref Reaktoro::StandardThermoModelWaterHKF
- @ref Reaktoro::StandardVolumeModel
- @ref Reaktoro::StandardVolumeModelConstant

## Reaction thermodynamic properties and models

- @ref Reaktoro::ReactionStandardThermoProps
- @ref Reaktoro::ReactionStandardThermoModelConstLgK
- @ref Reaktoro::ReactionStandardThermoModelGemsLgK
- @ref Reaktoro::ReactionStandardThermoModelPhreeqcLgK
- @ref Reaktoro::ReactionStandardThermoModelPressureCorrection
- @ref Reaktoro::ReactionStandardThermoModelVantHoff

## Water thermodynamic and electrostatic properties

- @ref Reaktoro::WaterElectroProps
- @ref Reaktoro::WaterThermoProps
- @ref Reaktoro::WaterHelmholtzProps
- @ref Reaktoro::waterDensityHGK
- @ref Reaktoro::waterDensityWagnerPruss
- @ref Reaktoro::waterLiquidDensityHGK
- @ref Reaktoro::waterLiquidDensityWagnerPruss
- @ref Reaktoro::waterVaporDensityHGK
- @ref Reaktoro::waterVaporDensityWagnerPruss
- @ref Reaktoro::waterPressureHGK
- @ref Reaktoro::waterPressureWagnerPruss
- @ref Reaktoro::waterSaturationPressureWagnerPruss
- @ref Reaktoro::waterSaturationLiquidDensityWagnerPruss
- @ref Reaktoro::waterSaturationVapourDensityWagnerPruss
- @ref Reaktoro::waterElectroPropsJohnsonNorton
- @ref Reaktoro::waterThermoPropsHGK
- @ref Reaktoro::waterThermoPropsWagnerPruss
- @ref Reaktoro::waterHelmholtzPropsHGK
- @ref Reaktoro::waterHelmholtzPropsWagnerPruss

## Worth checking classes and methods

- @ref Reaktoro::Material
- @ref Reaktoro::Param
- @ref Reaktoro::Params
- @ref Reaktoro::Table
- @ref Reaktoro::Data
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
- @ref Reaktoro::StateOfMatter

## For all other needs

Check the namespace Reaktoro for all available classes and methods.
