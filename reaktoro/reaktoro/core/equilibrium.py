from reaktoro.core import *
import numpy as npy

class EquilibriumResult:
    def __init__(self):
        # The flag that indicates if the equilibrium calculations succeeded
        self.succeeded = False

        # The index of the degree of freedom where the equilibrium calculation failed
        self.failed_where = None

        # The tuple `(T, P, b)` containing the temperature `T`, pressure `P`
        # and elemental molar amounts `b` that caused a failure
        self.failed_with = None

        # The list of number of iterations performed for every degree of freedom
        self.iterations = []

        # The maximum number of iterations
        self.iterations_max = 0

        # The minimum number of iterations
        self.iterations_min = float('inf')

        # The average number of iterations
        self.iterations_average = 0

        # The total time of the equilibrium calculations
        self.time = []

        # The total time spent on objective evaluations
        self.time_objective_evals = 0

        # The total time spent on constraint evaluations
        self.time_constraint_evals = 0

        # The total time spent on linear systems
        self.time_linear_systems = 0

class Equilibrium:
    def __init__(self, system, partition, **kwargs):
        # Set the data members
        self.system = system
        self.partition = partition
        self.solver = Reaktor.EquilibriumSolver()
        self.problem = Reaktor.EquilibriumProblem(system, partition)
        self.gems = kwargs.get("gems", None)

        # Set some options for the equilibrium solver
        options = Reaktor.EquilibriumOptions()
        options.optimum.max_iterations = 1000
        options.epsilon = min(options.epsilon, 1e-20)
        self.solver.setOptions(options)

        self.result = EquilibriumResult()

    def setOptions(self, options):
        self.solver.setOptions(options)

    def solve(self, be, states):
        # Check if Gems has been set, if so, use it as equilibrator
        if self.gems is not None:
            return self._solveWithGems(be, states)

        # Create a result instance for access to info about the calculation
        self.result = EquilibriumResult()
        self.result.iterations = npy.full(len(states), 9999)
        self.result.time = npy.empty(len(states))

        # Iterate over all degrees of freedom and calculate their equilibrium
        for k in range(len(states)):
            # Define auxiliary variables
            T = states[k].temperature()
            P = states[k].pressure()
            b = Reaktor.Vector(be[k])

            # Set temperature, pressure and element molar amounts from the current degree of freedom
            self.problem.setTemperature(T)
            self.problem.setPressure(P)
            self.problem.setElementAmounts(b)

            # Solve the equilibrium problem
            res = self.solver.solve(self.problem, states[k])

            # Extract result information of the last equilibrium calculation
            self.result.succeeded = res.optimum.succeeded
            self.result.time[k] = res.optimum.time
            self.result.time_objective_evals += res.optimum.time_objective_evals
            self.result.time_constraint_evals += res.optimum.time_constraint_evals
            self.result.time_linear_systems += res.optimum.time_linear_systems
            self.result.iterations[k] = res.optimum.iterations

            # Check if the calculation failed
            if not self.result.succeeded:
                self.result.failed_where = k
                self.result.failed_with = (T, P, b)
                return self.result

        self.result.iterations_min = min(self.result.iterations)
        self.result.iterations_max = max(self.result.iterations)
        self.result.iterations_average = sum(self.result.iterations)/float(len(self.result.iterations))

        return self.result

    def _solveWithGems(self, be, states):
        self.result = EquilibriumResult()
        self.result.iterations = npy.full(len(states), 9999)
        self.result.time = npy.empty(len(states))

        # Iterate over all degrees of freedom and calculate their equilibrium
        for k in range(len(states)):
            # Define auxiliary variables
            T = states[k].temperature()
            P = states[k].pressure()
            b = Reaktor.Vector(be[k])

            # Set temperature, pressure and element molar amounts from the current degree of freedom
            self.gems.setTemperature(T)
            self.gems.setPressure(P)
            self.gems.setElementAmounts(b)

            # Solve the equilibrium problem
            self.gems.equilibrate()

            # Set the equilibrium molar composition at the current degree of freedom
            states[k].setSpeciesAmounts(self.gems.speciesAmounts())

            # Extract result information of the last equilibrium calculation
            self.result.succeeded = self.gems.converged()
            self.result.time[k] = self.gems.elapsedTime()
            self.result.iterations[k] = self.gems.numIterations()

            # Check if the calculation failed
            if not self.result.succeeded:
                self.result.failed_where = k
                self.result.failed_with = (T, P, b)
                return self.result

        self.result.iterations_min = min(self.result.iterations)
        self.result.iterations_max = max(self.result.iterations)
        self.result.iterations_average = sum(self.result.iterations)/float(len(self.result.iterations))

        return self.result
