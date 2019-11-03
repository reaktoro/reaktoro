# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import *

T, P = 298.15, 1e5

wts = waterThermoStateWagnerPruss(T, P, StateOfMatter.Liquid)

print("temperature.val = ", wts.temperature.val)
print("temperature.ddT = ", wts.temperature.ddT)
print("temperature.ddP = ", wts.temperature.ddP)
print("volume.val = ", wts.volume.val)
print("volume.ddT = ", wts.volume.ddT)
print("volume.ddP = ", wts.volume.ddP)
print("entropy.val = ", wts.entropy.val)
print("entropy.ddT = ", wts.entropy.ddT)
print("entropy.ddP = ", wts.entropy.ddP)
print("helmholtz.val = ", wts.helmholtz.val)
print("helmholtz.ddT = ", wts.helmholtz.ddT)
print("helmholtz.ddP = ", wts.helmholtz.ddP)
print("internal_energy.val = ", wts.internal_energy.val)
print("internal_energy.ddT = ", wts.internal_energy.ddT)
print("internal_energy.ddP = ", wts.internal_energy.ddP)
print("enthalpy.val = ", wts.enthalpy.val)
print("enthalpy.ddT = ", wts.enthalpy.ddT)
print("enthalpy.ddP = ", wts.enthalpy.ddP)
print("gibbs.val = ", wts.gibbs.val)
print("gibbs.ddT = ", wts.gibbs.ddT)
print("gibbs.ddP = ", wts.gibbs.ddP)
print("cv.val = ", wts.cv.val)
print("cv.ddT = ", wts.cv.ddT)
print("cv.ddP = ", wts.cv.ddP)
print("cp.val = ", wts.cp.val)
print("cp.ddT = ", wts.cp.ddT)
print("cp.ddP = ", wts.cp.ddP)
print("density.val = ", wts.density.val)
print("density.ddT = ", wts.density.ddT)
print("density.ddP = ", wts.density.ddP)
print("densityT.val = ", wts.densityT.val)
print("densityT.ddT = ", wts.densityT.ddT)
print("densityT.ddP = ", wts.densityT.ddP)
print("densityP.val = ", wts.densityP.val)
print("densityP.ddT = ", wts.densityP.ddT)
print("densityP.ddP = ", wts.densityP.ddP)
print("densityTT.val = ", wts.densityTT.val)
print("densityTT.ddT = ", wts.densityTT.ddT)
print("densityTT.ddP = ", wts.densityTT.ddP)
print("densityTP.val = ", wts.densityTP.val)
print("densityTP.ddT = ", wts.densityTP.ddT)
print("densityTP.ddP = ", wts.densityTP.ddP)
print("densityPP.val = ", wts.densityPP.val)
print("densityPP.ddT = ", wts.densityPP.ddT)
print("densityPP.ddP = ", wts.densityPP.ddP)
print("pressure.val = ", wts.pressure.val)
print("pressure.ddT = ", wts.pressure.ddT)
print("pressure.ddP = ", wts.pressure.ddP)
print("pressureT.val = ", wts.pressureT.val)
print("pressureT.ddT = ", wts.pressureT.ddT)
print("pressureT.ddP = ", wts.pressureT.ddP)
print("pressureD.val = ", wts.pressureD.val)
print("pressureD.ddT = ", wts.pressureD.ddT)
print("pressureD.ddP = ", wts.pressureD.ddP)
print("pressureTT.val = ", wts.pressureTT.val)
print("pressureTT.ddT = ", wts.pressureTT.ddT)
print("pressureTT.ddP = ", wts.pressureTT.ddP)
print("pressureTD.val = ", wts.pressureTD.val)
print("pressureTD.ddT = ", wts.pressureTD.ddT)
print("pressureTD.ddP = ", wts.pressureTD.ddP)
print("pressureDD.val = ", wts.pressureDD.val)
print("pressureDD.ddT = ", wts.pressureDD.ddT)
print("pressureDD.ddP = ", wts.pressureDD.ddP)
