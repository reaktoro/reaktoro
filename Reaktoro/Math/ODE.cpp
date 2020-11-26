// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ODE.hpp"

// Sundials includes
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Math.hpp>

namespace Reaktoro {

#define VecEntry(v, i)    NV_Ith_S(v, i)
#define MatEntry(A, i, j) DENSE_ELEM(A, i, j)

#define CheckInitialize(r) \
    Assert(r == CV_SUCCESS, \
        "Cannot proceed with ODESolver::initialize to initialize the solver.", \
        "There was a failure in the CVODE call `" + std::string(#r) + "`."); \

#define CheckIntegration(r) \
    Assert(r == CV_SUCCESS, \
        "Cannot proceed with ODESolver::integrate to integrate the differential equations.", \
        "There was a failure in the CVODE call `" + std::string(#r) + "`."); \

using N_Vector = struct _generic_N_Vector*;

inline int CVODEStep(const ODEStepMode& step);
inline int CVODEIteration(const ODEIterationMode& iteration);
inline int CVODEMaxStepOrder(const ODEOptions& options);

int CVODEFunction(realtype t, N_Vector y, N_Vector ydot, void* user_data);
int CVODEJacobian(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

struct ODEData
{
    ODEData(const ODEProblem& problem, VectorRef y, VectorRef f, MatrixRef J)
    : problem(problem), y(y), f(f), J(J), num_equations(problem.numEquations())
    {}

    const ODEProblem& problem;
    VectorRef y;
    VectorRef f;
    MatrixRef J;
    Index num_equations;
};

struct ODEProblem::Impl
{
    /// The number of ordinary differential equations
    Index num_equations = 0;

    /// The right-hand side function of the system of ordinary differential equations
    ODEFunction ode_function;

    /// The Jacobian of the right-hand side function of the system of ordinary differential equations
    ODEJacobian ode_jacobian;
};

struct ODESolver::Impl
{
    /// The ODE problem
    ODEProblem problem;

    /// The options for the ODE integration
    ODEOptions options;

    /// The CVODE context pointer
    void* cvode_mem;

    /// The CVODE vector y
    N_Vector cvode_y;

    /// The auxiliary vector f for the function evaluation
    Vector f;

    /// The auxiliary matrix J for the Jacobian evaluation
    Matrix J;

    /// Construct a default ODESolver::Impl instance
    Impl()
    : cvode_mem(nullptr), cvode_y(nullptr)
    {}

    ~Impl()
    {
        // Free the dynamic memory allocated for cvode context
        if(cvode_mem) CVodeFree(&cvode_mem);

        // Free dynamic memory allocated for y
        if(cvode_y) N_VDestroy_Serial(cvode_y);
    }

    /// Initializes the ODE solver
    auto initialize(double tstart, VectorConstRef y) -> void
    {
        // Check if the ordinary differential problem has been initialized
        Assert(problem.initialized(),
            "Cannot proceed with ODESolver::initialize to initialize the solver.",
            "The provided ODEProblem instance was not properly initialized.");

        // Check if the dimension of 'y' matches the number of equations
        Assert(y.size() == problem.numEquations(),
            "Cannot proceed with ODESolver::initialize to initialize the solver.",
            "The dimension of the vector parameter `y` does not match the number of equations.");

        // The number of differential equations
        const unsigned int num_equations = problem.numEquations();

        // Allocate memory for f and J
        f.resize(num_equations);
        J.resize(num_equations, num_equations);

        // Free any dynamic memory allocated for cvode_mem context
        if(cvode_mem) CVodeFree(&cvode_mem);

        // Free dynamic memory allocated for y
        if(cvode_y) N_VDestroy_Serial(cvode_y);

        // Initialize a new vector y
        cvode_y = N_VNew_Serial(num_equations);
        for(Index i = 0; i < num_equations; ++i)
            VecEntry(cvode_y, i) = y[i];

        // Initialize a new cvode context
        cvode_mem = CVodeCreate(CVODEStep(options.step), CVODEIteration(options.iteration));

        // Set the maximum order of the Adams or BDF methods
        CheckInitialize(CVodeSetMaxOrd(cvode_mem, CVODEMaxStepOrder(options)));

        // Check if the cvode creation succeeded
        Assert(cvode_mem != nullptr,
            "Cannot proceed with ODESolver::initialize to initialize the solver.",
            "There was an error creating the CVODE context.")

        // Initialize the cvode context
        CheckInitialize(CVodeInit(cvode_mem, CVODEFunction, tstart, cvode_y));

        // Initialize the vector of absolute tolerances
        N_Vector abstols = N_VNew_Serial(num_equations);

        if(options.abstols.size() == num_equations)
            for(Index i = 0; i < num_equations; ++i)
                VecEntry(abstols, i) = options.abstols[i];
        else
            for(Index i = 0; i < num_equations; ++i)
                VecEntry(abstols, i) = options.abstol;

        // Set the parameters for the calculation
        CheckInitialize(CVodeSetStabLimDet(cvode_mem, options.stability_limit_detection));
        CheckInitialize(CVodeSetInitStep(cvode_mem, options.initial_step));
        if(options.stop_time) CheckInitialize(CVodeSetStopTime(cvode_mem, options.stop_time));
        CheckInitialize(CVodeSetMinStep(cvode_mem, options.min_step));
        CheckInitialize(CVodeSetMaxStep(cvode_mem, options.max_step));
        CheckInitialize(CVodeSetMaxNumSteps(cvode_mem, int(options.max_num_steps)));
        CheckInitialize(CVodeSetMaxHnilWarns(cvode_mem, int(options.max_hnil_warnings)));
        CheckInitialize(CVodeSetMaxErrTestFails(cvode_mem, int(options.max_num_error_test_failures)));
        CheckInitialize(CVodeSetMaxNonlinIters(cvode_mem, int(options.max_num_nonlinear_iterations)));
        CheckInitialize(CVodeSetMaxConvFails(cvode_mem, int(options.max_num_convergence_failures)));
        CheckInitialize(CVodeSetNonlinConvCoef(cvode_mem, options.nonlinear_convergence_coefficient));
        CheckInitialize(CVodeSVtolerances(cvode_mem, options.reltol, abstols));

        // Call CVDense to specify the CVDENSE dense linear solver
        CheckInitialize(CVDense(cvode_mem, num_equations));

        // Set the Jacobian function
        CheckInitialize(CVDlsSetDenseJacFn(cvode_mem, CVODEJacobian));

        // Free dynamic memory allocated for `yc`
        N_VDestroy_Serial(abstols);
    }

    auto initializeWithoutCVODE(VectorConstRef y) -> void
    {
        // Check if the ordinary differential problem has been initialized
        Assert(problem.initialized(),
               "Cannot proceed with ODESolver::initialize to initialize the solver.",
               "The provided ODEProblem instance was not properly initialized.");

        // Check if the dimension of 'y' matches the number of equations
        Assert(y.size() == problem.numEquations(),
               "Cannot proceed with ODESolver::initialize to initialize the solver.",
               "The dimension of the vector parameter `y` does not match the number of equations.");

        // The number of differential equations
        const unsigned int num_equations = problem.numEquations();

        // Allocate memory for f and J
        f.resize(num_equations);
        J.resize(num_equations, num_equations);

    }

    /// Integrate the ODE performing a single step.
    auto integrate(double& t, VectorRef y) -> void
    {
        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // Define an infinite time.
        double tfinal = 10*(t + 1);

        // Set the user-defined data to cvode_mem
        CheckIntegration(CVodeSetUserData(cvode_mem, &data));

        // Solve the ode problem from `tstart` to `tfinal`
        CheckIntegration(CVode(cvode_mem, tfinal, cvode_y, &t, CV_ONE_STEP));

        // Transfer the result from cvode_y to y
        for(Index i = 0; i < data.num_equations; ++i)
            y[i] = VecEntry(cvode_y, i);
    }

    /// Integrate the ODE performing a single step not going over a given time.
    auto integrate(double& t, VectorRef y, double tfinal) -> void
    {
        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // Set the user-defined data to cvode_mem
        CheckIntegration(CVodeSetUserData(cvode_mem, &data));

        // Solve the ode problem from `tstart` to `tfinal`
        CheckIntegration(CVode(cvode_mem, tfinal, cvode_y, &t, CV_ONE_STEP));

        // Check if the current time is now greater than the final time
        if(t > tfinal)
        {
            // Interpolate y at t using its old state and the new one
            CheckIntegration(CVodeGetDky(cvode_mem, tfinal, 0, cvode_y));

            // Adjust the current time
            t = tfinal;
        }

        // Transfer the result from cvode_y to y
        for(unsigned int i = 0; i < data.num_equations; ++i)
            y[i] = VecEntry(this->cvode_y, i);
    }

    /// Integrate the ODE performing a single step not going over a given time `tfinal`.
    auto integrate(double& t, VectorRef y, double tfinal, MatrixRef S) -> void
    {
        // Initialize the cvode context
        initialize(t, y);

        while(t < tfinal)
        {
            // Initialize the ODE data
            ODEData data(problem, y, f, J);

            // Set the user-defined data to cvode_mem
            CheckIntegration(CVodeSetUserData(cvode_mem, &data));

            double t0 = t;

            // Solve the ode problem from `tstart` to `tfinal`
            CheckIntegration(CVode(cvode_mem, tfinal, cvode_y, &t, CV_ONE_STEP));

            // Check if the current time is now greater than the final time
            if(t > tfinal)
            {
                // Interpolate y at t using its old state and the new one
                CheckIntegration(CVodeGetDky(cvode_mem, tfinal, 0, cvode_y));

                // Adjust the current time
                t = tfinal;
            }
            // Transfer the result from cvode_y to y
            for(unsigned int i = 0; i < data.num_equations; ++i)
                y[i] = VecEntry(this->cvode_y, i);

            // Define dt
            double dt = t - t0;

            // Identity matrix
            Matrix I = Matrix::Identity(data.num_equations, data.num_equations);

            // Set the value of the current jacobian
            problem.function(t, y, f);
            problem.jacobian(t, y, J);

            // Perform LU decomposition for matrix A = I - dt * J^{k+1}
            LU lu(I - dt * J);

            // Solve system of equation (I - dt * J^{k+1}) * S^{k+1} = S^k
            S = lu.solve(S);
        }
    }

    /*
    /// Solve the ODE equations from a given start time to a final one.
    auto solve(double& t, double dt, VectorRef y) -> void
    {
        // Initialize the cvode context
        initialize(t, y);

        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // Set the user-defined data to cvode_mem
        CheckIntegration(CVodeSetUserData(cvode_mem, &data));

        // Solve the ode problem from `tstart` to `tfinal`
        // CV_NORMAL - indicates that the step is controled by `tstart` and `tfinal`
        // CVODE will trigger integration from t until t + dt
        CheckIntegration(CVode(cvode_mem, t + dt, cvode_y, &t, CV_NORMAL));

        // Transfer the result from cvode_y to y
        for(Index i = 0; i < data.num_equations; ++i)
            y[i] = VecEntry(this->cvode_y, i);
    }
    */
    /// Solve the ODE equations from a given start time to a final one.
    auto solve(double& t, double dt, VectorRef y) -> void
    {
        // Initialize the cvode context
        initialize(t, y);

        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // Set the user-defined data to cvode_mem
        CheckIntegration(CVodeSetUserData(cvode_mem, &data));

        double t_start = t;
        // Solve the ode problem from `tstart` to `tfinal`
        // CV_NORMAL - indicates that the step is controlled by `tstart` and `tfinal`
        // CVODE will trigger integration from t until t + dt
        auto res = CVode(cvode_mem, t + dt, cvode_y, &t, CV_NORMAL);

        // If CVODE was not successful, run the 1st order implicit scheme
        if(res != CV_SUCCESS){
            // Apply auxiliary 1st order implicit scheme
            solveWithImplicitRKI(t_start, dt, y);
        }else{
            // Transfer the result from cvode_y to y
            for(unsigned int i = 0; i < data.num_equations; ++i)
                y[i] = VecEntry(this->cvode_y, i);
        }
    }

    /// Solve the ODE equations from a given start time to a final one.
    auto solve(double& t, double dt, VectorRef y, MatrixRef S) -> void
    {
        // Initialize the cvode context
        initialize(t, y);

        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // Set the user-defined data to cvode_mem
        CheckIntegration(CVodeSetUserData(cvode_mem, &data));

        double t_start = t;
        // Solve the ode problem from `tstart` to `tfinal`
        // CV_NORMAL - indicates that the step is controlled by `tstart` and `tfinal`
        // CVODE will trigger integration from t until t + dt
        auto res = CVode(cvode_mem, t + dt, cvode_y, &t, CV_NORMAL);


        // If CVODE was not successful, run the 1st order implicit scheme
        if(res != CV_SUCCESS){
            // Apply 1st order implicit scheme
            solveWithImplicitRKI(t_start, dt, y);
        }else{
            // Transfer the result from cvode_y to y
            for(unsigned int i = 0; i < data.num_equations; ++i)
                y[i] = VecEntry(this->cvode_y, i);
        }

        // -------------------------------------------------------------------
        // Reconstruct sensitivities using implicit scheme
        // -------------------------------------------------------------------

        // Initialize identity matrix and Jacobian on the t = t^{k+1}
        Matrix I = Matrix::Identity(data.num_equations, data.num_equations);

        // Evaluate current value of Jacobian and RHS
        problem.function(t, y, f); // provides evaluation of properties and rates at t = t_{k+1}
        problem.jacobian(t, y, J); // provides evaluation of sensitivities at t = t_{k+1}

        // Perform LU decomposition for matrix A = I - dt
        LU lu(I - dt * J);

        // Solve system of equation (I - dt * J^{k+1}) * S^{k+1} = S^k(=I)
        //S = lu.solve(I);
        S = lu.solve(S);
    }

    /// Solve the ODE equations from a given start time to a final one.
    auto solveWithImplicitRKI(double& t, double dt, VectorRef y) -> void
    {
        // Initialize y
        initializeWithoutCVODE(y);

        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // -------------------------------------------------------------------
        // Reconstruct unknown only y using implicit 1st order RK scheme
        // -------------------------------------------------------------------

        // Initialize identity matrix and Jacobian on the t = t^{k+1}
        Matrix I = Matrix::Identity(data.num_equations, data.num_equations);

        // Initial vector
        Vector yk1 = y, yk = y;
        // Residual and increment for the Newton method
        Vector r, dy;
        // Final time of integration
        double tk1 = t + dt;

        double error_rel = 1e8;
        double alpha = 1.0;
        Index i = 0;

        while(error_rel > 1e-4 * options.reltol) { // reltol = 1e-10

            // Evaluate current value of Jacobian and right-hand side
            problem.function(tk1, yk, f); // evaluate right-hand side at t_{k+1}
            problem.jacobian(tk1, yk, J); // evaluate Jacobian at t_{k+1}

            // Calculate residual
            r = yk1 - yk - dt * f;
            // Find dy solving J * dy = - r
            LU lu(I - dt * J);
            dy = lu.solve(-r);

            // Update yk1 and yk
            y = yk1; // remember old value of yk1
            yk1 = yk + alpha * dy; // update yk1
            yk = y; // update yk

            // Calculate the error
            error_rel = dy.squaredNorm() / y.squaredNorm();
            i++;

        }

        // Update final value
        y = yk1;
        // Update time
        t = t + dt;

    }

    /// Solve the ODE equations from a given start time to a final one.
    auto solveWithImplicitRKI(double& t, double dt, VectorRef y,  MatrixRef S) -> void
    {
        // Initialize y
        initializeWithoutCVODE(y);

        // Initialize the ODE data
        ODEData data(problem, y, f, J);

        // -------------------------------------------------------------------
        // Reconstruct unknown y using implicit 1st order RK scheme
        // -------------------------------------------------------------------

        // Initialize identity matrix and Jacobian on the t = t^{k+1}
        Matrix I = Matrix::Identity(data.num_equations, data.num_equations);

        // Initial vector
        Vector yk1 = y, yk = y;
        // Residual and increment for the Newton method
        Vector r, dy;
        // Final time of integration
        double tk1 = t + dt;

        double error_rel = 1e8;
        Index i = 0;

        while(error_rel > 1e-4 * options.reltol) {

            // Evalute current value of Jacobian and right-hand side
            problem.function(tk1, yk, f); // evaluate right-hand side at t_{k+1}
            problem.jacobian(tk1, yk, J); // evaluate jacobian at t_{k+1}

            // Calculate residual
            r = yk1 - yk - dt * f;
            // Find dy solving J * dy = - r
            LU lu(I - dt * J);
            dy = lu.solve(-r);

            // Update yk1 and yk
            y = yk1; // remember old value of yk1
            yk1 = yk + dy; // update yk1
            yk = y; // update yk

            // Calculate the error between two
            error_rel = dy.squaredNorm() / y.squaredNorm();
            i++;

        }

        // Update final value
        y = yk1;
        // Update time
        t = t + dt;

        // -------------------------------------------------------------------
        // Reconstruct sensitivities using implicit scheme
        // -------------------------------------------------------------------

        // Evaluate current value of Jacobian and RHS
        problem.function(t, y, f); // provides evaluation of properties and rates at t = t_{k+1}
        problem.jacobian(t, y, J); // provides evaluation of sensitivities at t = t_{k+1}

        // Perform LU decomposition for matrix  I - dt * J^{k+1}
        LU lu(I - dt * J);

        // Solve system of equation (I - dt * J^{k+1}) * S^{k+1} = S^k(=I)
        S = lu.solve(S);
    }
};

inline int CVODEStep(const ODEStepMode& step)
{
    switch(step)
    {
        case ODEStepMode::Adams: return CV_ADAMS;
        default: return CV_BDF;
    }
}

inline int CVODEIteration(const ODEIterationMode& iteration)
{
    switch(iteration)
    {
        case ODEIterationMode::Functional: return CV_FUNCTIONAL;
        default: return CV_NEWTON;
    }
}

inline int CVODEMaxStepOrder(const ODEOptions& options)
{
    switch(options.step)
    {
        case ODEStepMode::Adams: return options.max_order_adams;
        default: return options.max_order_bdf;
    }
}

int CVODEFunction(realtype t, N_Vector y, N_Vector f, void* user_data)
{
    ODEData& data = *static_cast<ODEData*>(user_data);

    for(Index i = 0; i < data.num_equations; ++i)
        data.y[i] = VecEntry(y, i);

    int result = data.problem.function(t, data.y, data.f);

    for(Index i = 0; i < data.num_equations; ++i)
        VecEntry(f, i) = data.f[i];

    return result;
}

int CVODEJacobian(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    ODEData& data = *static_cast<ODEData*>(user_data);

    for(Index i = 0; i < data.num_equations; ++i)
        data.y[i] = VecEntry(y, i);

    int result = data.problem.jacobian(t, data.y, data.J);

    for(Index i = 0; i < data.num_equations; ++i)
        for(int j = 0; j < data.num_equations; ++j)
            MatEntry(J, i, j) = data.J(i, j);

    return result;
}

ODEProblem::ODEProblem()
: pimpl(new Impl())
{}

ODEProblem::ODEProblem(const ODEProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

ODEProblem::~ODEProblem()
{}

auto ODEProblem::operator=(ODEProblem other) -> ODEProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ODEProblem::setNumEquations(unsigned num) -> void
{
    pimpl->num_equations = num;
}

auto ODEProblem::setFunction(const ODEFunction& f) -> void
{
    pimpl->ode_function = f;
}

auto ODEProblem::setJacobian(const ODEJacobian& J) -> void
{
    pimpl->ode_jacobian = J;
}

auto ODEProblem::initialized() const -> bool
{
    return numEquations() && function();
}

auto ODEProblem::numEquations() const -> unsigned
{
    return pimpl->num_equations;
}

auto ODEProblem::function() const -> const ODEFunction&
{
    return pimpl->ode_function;
}

auto ODEProblem::jacobian() const -> const ODEJacobian&
{
    return pimpl->ode_jacobian;
}

auto ODEProblem::function(double t, VectorConstRef y, VectorRef f) const -> int
{
    return function()(t, y, f);
}

auto ODEProblem::jacobian(double t, VectorConstRef y, MatrixRef J) const -> int
{
    return jacobian()(t, y, J);
}

ODESolver::ODESolver()
: pimpl(new Impl())
{}

ODESolver::ODESolver(const ODESolver& other)
: pimpl(new Impl(*other.pimpl))
{}

ODESolver::~ODESolver()
{}

auto ODESolver::operator=(ODESolver other) -> ODESolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ODESolver::setOptions(const ODEOptions& options) -> void
{
    pimpl->options = options;
}

auto ODESolver::setProblem(const ODEProblem& problem) -> void
{
    pimpl->problem = problem;
}

auto ODESolver::initialize(double tstart, VectorConstRef y) -> void
{
    pimpl->initialize(tstart, y);
}

auto ODESolver::integrate(double& t, VectorRef y) -> void
{
    pimpl->integrate(t, y);
}

auto ODESolver::integrate(double& t, VectorRef y, double tfinal) -> void
{
    pimpl->integrate(t, y, tfinal);
}

auto ODESolver::integrate(double& t, VectorRef y, double tfinal, MatrixRef S) -> void
{
    pimpl->integrate(t, y, tfinal, S);
}

auto ODESolver::solve(double& t, double dt, VectorRef y) -> void
{
    pimpl->solve(t, dt, y);
}

auto ODESolver::solve(double& t, double dt, VectorRef y, MatrixRef S) -> void
{
    pimpl->solve(t, dt, y, S);
}

// Methods independent of CVODE

auto ODESolver::solveWithImplicitRKI(double& t, double dt, VectorRef y) -> void
{
    pimpl->solveWithImplicitRKI(t, dt, y);
}

auto ODESolver::solveWithImplicitRKI(double& t, double dt, VectorRef y,  MatrixRef S) -> void
{
    pimpl->solveWithImplicitRKI(t, dt, y, S);
}

} // namespace Reaktoro
