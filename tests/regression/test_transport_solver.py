import numpy as np
import pytest
import reaktoro as rkt

from collections import namedtuple

# a and b values of a linear source: q = a*x+b
source_parameters = namedtuple("source_parameters", ["a", "b"])


@pytest.mark.parametrize(
    "source_parameters",
    [
        pytest.param(source_parameters(0, 0), id="diffusion-source rate q=0"),
        pytest.param(source_parameters(0, 1), id="diffusion-source rate q=1"),
        pytest.param(source_parameters(1, 1), id="diffusion-source rate q=x+1"),
    ],
)
def test_transport_solver_diffusion(source_parameters, num_regression):
    """
    A test to check the solution of a advection-diffusion equation with v = 0 
    and du/dt = 0
    Eq: 
        du/dt + v du/dx = D d²u/dx² + q
    
    The result were compared with the following analytic solution and got
    near results when we increase the number of cells. We decided to compare with
    numerical solve to save some computational time and analytic result, but with
    relative error of 1e-1.
    
    analytic_u = -((a*x**3)/(6*D)) - ((b*x**2)/(2*D)) + ((a*x*xr**2)/(2*D)) + ((b*x*xr)/(D)) + ul
    
    @param source_parameters
        a tuple that has values of a and b coefficient of a source term
        that behaves as q a*x+b
    """
    a = source_parameters.a
    b = source_parameters.b
    D = 0.0002
    v = 0
    ul = 1
    dt = 100
    num_steps = 1000
    num_cells = 10
    xl = 0
    xr = 1.0

    mesh = rkt.Mesh(num_cells, xl, xr)

    x = mesh.xcells()
    q = a * x + b

    transp_solver = rkt.TransportSolver()

    transp_solver.setMesh(mesh)
    transp_solver.setVelocity(v)
    transp_solver.setDiffusionCoeff(D)
    transp_solver.setBoundaryValue(ul)
    transp_solver.setTimeStep(dt)

    transp_solver.initialize()

    numerical_u = np.zeros(num_cells)
    transp_solver.step(numerical_u, q)

    for i in range(num_steps):
        transp_solver.step(numerical_u, q)

    num_regression.check({"u": numerical_u})
    
    analytic_u = -(a*x**3)/(6*D) - (b*x**2)/(2*D) + (a*x*xr**2)/(2*D) + (b*x*xr)/D + ul
    
    assert numerical_u == pytest.approx(np.array(analytic_u), rel=0.1)

@pytest.mark.parametrize(
    "source_parameters",
    [
        pytest.param(source_parameters(0, 0), id="advection-with source rate q=0"),
        pytest.param(source_parameters(0, 1), id="advection-with source rate q=1"),
        pytest.param(source_parameters(1, 1), id="advection-with source rate q=x+1"),
    ],
)
def test_transport_solver_advection(source_parameters, num_regression):
    """
    A test to check the solution of a advection-diffusion equation with D = 0 
    and du/dt = 0
    Eq: 
        du/dt + v du/dx = D d²u/dx² + q
    
    The result was compared with the following analytic solution and got
    near results when we increase the number of cells. We decided to compare with
    numerical solve to save some computational time and analytic result, but with
    relative error of 1e-1
    
    analytic_u = (a*x**2)/(2*v) + (b*x)/(v) + ul
    
    @param source_parameters
        a tuple that has values of "a" and "b" coefficients of a source term
        that behaves as q=a*x+b
    """
    a = source_parameters.a
    b = source_parameters.b
    D = 0
    v = 1
    ul = 1
    dt = 0.01
    num_steps = 5000
    num_cells = 10
    xl = 0
    xr = 1.0

    mesh = rkt.Mesh(num_cells, xl, xr)

    x = mesh.xcells()
    q = a * x + b

    transp_solver = rkt.TransportSolver()

    transp_solver.setMesh(mesh)
    transp_solver.setVelocity(v)
    transp_solver.setDiffusionCoeff(D)
    transp_solver.setBoundaryValue(ul)
    transp_solver.setTimeStep(dt)

    transp_solver.initialize()

    numerical_u = np.ones(num_cells)
    transp_solver.step(numerical_u, q)

    for i in range(num_steps):
        transp_solver.step(numerical_u, q)

    num_regression.check({"u": numerical_u})
    
    analytic_u = (a*x**2)/(2*v) + (b*x)/v + ul
    
    assert numerical_u == pytest.approx(np.array(analytic_u), rel=0.1)
