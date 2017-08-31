# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2015 Allan Leal
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import sys; sys.path.append("/home/allan/git/Reaktoro/build/debug/lib/python2.7/site-packages")

from reaktoro import *

def printstate(T, P):
    D = waterDensityWagnerPruss(T, P)
    D = waterDensityWagnerPruss(T, P)

    print "D.val", D.val
    print "D.ddT", D.ddT
    print "D.ddP", D.ddP, "\n"

    whs = waterHelmholtzStateWagnerPruss(T, D)
    whsT = waterHelmholtzStateWagnerPruss(T*(1+1e-8), D)
    # D = waterDensityHGK(T, P)
    # whs = waterHelmholtzStateHGK(T, D)
    print "helmholtz.val", whs.helmholtz.val
    print "helmholtz.ddT", whs.helmholtz.ddT
    print "helmholtz.ddTnum", (whsT.helmholtz.val-whs.helmholtz.val)/(T*1e-8)
    print "helmholtz.ddP", whs.helmholtz.ddP, "\n"
    print "helmholtzT.val", whs.helmholtzT.val
    print "helmholtzT.ddT", whs.helmholtzT.ddT
    print "helmholtzT.ddTnum", (whsT.helmholtzT.val-whs.helmholtzT.val)/(T*1e-8)
    print "helmholtzT.ddP", whs.helmholtzT.ddP, "\n"
    print "helmholtzD.val", whs.helmholtzD.val
    print "helmholtzD.ddT", whs.helmholtzD.ddT
    print "helmholtzD.ddTnum", (whsT.helmholtzD.val-whs.helmholtzD.val)/(T*1e-8)
    print "helmholtzD.ddP", whs.helmholtzD.ddP, "\n"
    print "helmholtzTT.val", whs.helmholtzTT.val
    print "helmholtzTT.ddT", whs.helmholtzTT.ddT
    print "helmholtzTT.ddTnum", (whsT.helmholtzTT.val-whs.helmholtzTT.val)/(T*1e-8)
    print "helmholtzTT.ddP", whs.helmholtzTT.ddP, "\n"
    print "helmholtzTD.val", whs.helmholtzTD.val
    print "helmholtzTD.ddT", whs.helmholtzTD.ddT
    print "helmholtzTD.ddTnum", (whsT.helmholtzTD.val-whs.helmholtzTD.val)/(T*1e-8)
    print "helmholtzTD.ddP", whs.helmholtzTD.ddP, "\n"
    print "helmholtzDD.val", whs.helmholtzDD.val
    print "helmholtzDD.ddT", whs.helmholtzDD.ddT
    print "helmholtzDD.ddTnum", (whsT.helmholtzDD.val-whs.helmholtzDD.val)/(T*1e-8)
    print "helmholtzDD.ddP", whs.helmholtzDD.ddP, "\n"
    print "helmholtzTTT.val", whs.helmholtzTTT.val
    print "helmholtzTTT.ddT", whs.helmholtzTTT.ddT
    print "helmholtzTTT.ddTnum", (whsT.helmholtzTTT.val-whs.helmholtzTTT.val)/(T*1e-8)
    print "helmholtzTTT.ddP", whs.helmholtzTTT.ddP, "\n"
    print "helmholtzTTD.val", whs.helmholtzTTD.val
    print "helmholtzTTD.ddT", whs.helmholtzTTD.ddT
    print "helmholtzTTD.ddTnum", (whsT.helmholtzTTD.val-whs.helmholtzTTD.val)/(T*1e-8)
    print "helmholtzTTD.ddP", whs.helmholtzTTD.ddP, "\n"
    print "helmholtzTDD.val", whs.helmholtzTDD.val
    print "helmholtzTDD.ddT", whs.helmholtzTDD.ddT
    print "helmholtzTDD.ddTnum", (whsT.helmholtzTDD.val-whs.helmholtzTDD.val)/(T*1e-8)
    print "helmholtzTDD.ddP", whs.helmholtzTDD.ddP, "\n"
    print "helmholtzDDD.val", whs.helmholtzDDD.val
    print "helmholtzDDD.ddT", whs.helmholtzDDD.ddT
    print "helmholtzDDD.ddTnum", (whsT.helmholtzDDD.val-whs.helmholtzDDD.val)/(T*1e-8)
    print "helmholtzDDD.ddP", whs.helmholtzDDD.ddP, "\n"

    wts  = waterThermoStateWagnerPruss(T, P)
    wtsT = waterThermoStateWagnerPruss(T+T*1e-8, P)
    wtsP = waterThermoStateWagnerPruss(T, P+P*1e-8)
    # wts  = waterThermoStateHGK(T, P)
    # wtsT = waterThermoStateHGK(T+T*1e-8, P)
    # wtsP = waterThermoStateHGK(T, P+P*1e-8)
    print "temperature.val", wts.temperature.val
    print "temperature.ddT", wts.temperature.ddT
    print "temperature.ddP", wts.temperature.ddP
    print "temperature.ddTe", abs(wts.temperature.ddT - (wtsT.temperature.val - wts.temperature.val)/(T*1e-8))/abs(wts.temperature.val)
    print "temperature.ddPe", abs(wts.temperature.ddP - (wtsP.temperature.val - wts.temperature.val)/(P*1e-8))/abs(wts.temperature.val), "\n"
    print "volume.val", wts.volume.val
    print "volume.ddT", wts.volume.ddT
    print "volume.ddP", wts.volume.ddP
    print "volume.ddTe", abs(wts.volume.ddT - (wtsT.volume.val - wts.volume.val)/(T*1e-8))/abs(wts.volume.val)
    print "volume.ddPe", abs(wts.volume.ddP - (wtsP.volume.val - wts.volume.val)/(P*1e-8))/abs(wts.volume.val), "\n"
    print "entropy.val", wts.entropy.val
    print "entropy.ddT", wts.entropy.ddT
    print "entropy.ddP", wts.entropy.ddP
    print "entropy.ddTe", abs(wts.entropy.ddT - (wtsT.entropy.val - wts.entropy.val)/(T*1e-8))/abs(wts.entropy.val)
    print "entropy.ddPe", abs(wts.entropy.ddP - (wtsP.entropy.val - wts.entropy.val)/(P*1e-8))/abs(wts.entropy.val), "\n"
    print "helmholtz.val", wts.helmholtz.val
    print "helmholtz.ddT", wts.helmholtz.ddT
    print "helmholtz.ddP", wts.helmholtz.ddP
    print "helmholtz.ddTe", abs(wts.helmholtz.ddT - (wtsT.helmholtz.val - wts.helmholtz.val)/(T*1e-8))/abs(wts.helmholtz.val)
    print "helmholtz.ddPe", abs(wts.helmholtz.ddP - (wtsP.helmholtz.val - wts.helmholtz.val)/(P*1e-8))/abs(wts.helmholtz.val), "\n"
    print "internal_energy.val", wts.internal_energy.val
    print "internal_energy.ddT", wts.internal_energy.ddT
    print "internal_energy.ddP", wts.internal_energy.ddP
    print "internal_energy.ddTe", abs(wts.internal_energy.ddT - (wtsT.internal_energy.val - wts.internal_energy.val)/(T*1e-8))/abs(wts.internal_energy.val)
    print "internal_energy.ddPe", abs(wts.internal_energy.ddP - (wtsP.internal_energy.val - wts.internal_energy.val)/(P*1e-8))/abs(wts.internal_energy.val), "\n"
    print "enthalpy.val", wts.enthalpy.val
    print "enthalpy.ddT", wts.enthalpy.ddT
    print "enthalpy.ddP", wts.enthalpy.ddP
    print "enthalpy.ddTe", abs(wts.enthalpy.ddT - (wtsT.enthalpy.val - wts.enthalpy.val)/(T*1e-8))/abs(wts.enthalpy.val)
    print "enthalpy.ddPe", abs(wts.enthalpy.ddP - (wtsP.enthalpy.val - wts.enthalpy.val)/(P*1e-8))/abs(wts.enthalpy.val), "\n"
    print "gibbs.val", wts.gibbs.val
    print "gibbs.ddT", wts.gibbs.ddT
    print "gibbs.ddP", wts.gibbs.ddP
    print "gibbs.ddTe", abs(wts.gibbs.ddT - (wtsT.gibbs.val - wts.gibbs.val)/(T*1e-8))/abs(wts.gibbs.val)
    print "gibbs.ddPe", abs(wts.gibbs.ddP - (wtsP.gibbs.val - wts.gibbs.val)/(P*1e-8))/abs(wts.gibbs.val), "\n"
    print "cv.val", wts.cv.val
    print "cv.ddT", wts.cv.ddT
    print "cv.ddP", wts.cv.ddP
    print "cv.ddTe", abs(wts.cv.ddT - (wtsT.cv.val - wts.cv.val)/(T*1e-8))/abs(wts.cv.val)
    print "cv.ddPe", abs(wts.cv.ddP - (wtsP.cv.val - wts.cv.val)/(P*1e-8))/abs(wts.cv.val), "\n"
    print "cp.val", wts.cp.val
    print "cp.ddT", wts.cp.ddT
    print "cp.ddP", wts.cp.ddP
    print "cp.ddTe", abs(wts.cp.ddT - (wtsT.cp.val - wts.cp.val)/(T*1e-8))/abs(wts.cp.val)
    print "cp.ddPe", abs(wts.cp.ddP - (wtsP.cp.val - wts.cp.val)/(P*1e-8))/abs(wts.cp.val), "\n"
    print "density.val", wts.density.val
    print "density.ddT", wts.density.ddT
    print "density.ddP", wts.density.ddP
    print "density.ddTe", abs(wts.density.ddT - (wtsT.density.val - wts.density.val)/(T*1e-8))/abs(wts.density.val)
    print "density.ddPe", abs(wts.density.ddP - (wtsP.density.val - wts.density.val)/(P*1e-8))/abs(wts.density.val), "\n"
    print "densityT.val", wts.densityT.val
    print "densityT.ddT", wts.densityT.ddT
    print "densityT.ddP", wts.densityT.ddP
    print "densityT.ddTe", abs(wts.densityT.ddT - (wtsT.densityT.val - wts.densityT.val)/(T*1e-8))/abs(wts.densityT.val)
    print "densityT.ddPe", abs(wts.densityT.ddP - (wtsP.densityT.val - wts.densityT.val)/(P*1e-8))/abs(wts.densityT.val), "\n"
    print "densityP.val", wts.densityP.val
    print "densityP.ddT", wts.densityP.ddT
    print "densityP.ddP", wts.densityP.ddP
    print "densityP.ddTe", abs(wts.densityP.ddT - (wtsT.densityP.val - wts.densityP.val)/(T*1e-8))/abs(wts.densityP.val)
    print "densityP.ddPe", abs(wts.densityP.ddP - (wtsP.densityP.val - wts.densityP.val)/(P*1e-8))/abs(wts.densityP.val), "\n"
    print "densityTT.val", wts.densityTT.val
    print "densityTT.ddT", wts.densityTT.ddT
    print "densityTT.ddP", wts.densityTT.ddP
    print "densityTT.ddTe", abs(wts.densityTT.ddT - (wtsT.densityTT.val - wts.densityTT.val)/(T*1e-8))/abs(wts.densityTT.val)
    print "densityTT.ddPe", abs(wts.densityTT.ddP - (wtsP.densityTT.val - wts.densityTT.val)/(P*1e-8))/abs(wts.densityTT.val), "\n"
    print "densityTP.val", wts.densityTP.val
    print "densityTP.ddT", wts.densityTP.ddT
    print "densityTP.ddP", wts.densityTP.ddP
    print "densityTP.ddTe", abs(wts.densityTP.ddT - (wtsT.densityTP.val - wts.densityTP.val)/(T*1e-8))/abs(wts.densityTP.val)
    print "densityTP.ddPe", abs(wts.densityTP.ddP - (wtsP.densityTP.val - wts.densityTP.val)/(P*1e-8))/abs(wts.densityTP.val), "\n"
    print "densityPP.val", wts.densityPP.val
    print "densityPP.ddT", wts.densityPP.ddT
    print "densityPP.ddP", wts.densityPP.ddP
    print "densityPP.ddTe", abs(wts.densityPP.ddT - (wtsT.densityPP.val - wts.densityPP.val)/(T*1e-8))/abs(wts.densityPP.val)
    print "densityPP.ddPe", abs(wts.densityPP.ddP - (wtsP.densityPP.val - wts.densityPP.val)/(P*1e-8))/abs(wts.densityPP.val), "\n"
    print "pressure.val", wts.pressure.val
    print "pressure.ddT", wts.pressure.ddT
    print "pressure.ddP", wts.pressure.ddP
    print "pressure.ddTe", abs(wts.pressure.ddT - (wtsT.pressure.val - wts.pressure.val)/(T*1e-8))/abs(wts.pressure.val)
    print "pressure.ddPe", abs(wts.pressure.ddP - (wtsP.pressure.val - wts.pressure.val)/(P*1e-8))/abs(wts.pressure.val), "\n"
    print "pressureT.val", wts.pressureT.val
    print "pressureT.ddT", wts.pressureT.ddT
    print "pressureT.ddP", wts.pressureT.ddP
    print "pressureT.ddTe", abs(wts.pressureT.ddT - (wtsT.pressureT.val - wts.pressureT.val)/(T*1e-8))/abs(wts.pressureT.val)
    print "pressureT.ddPe", abs(wts.pressureT.ddP - (wtsP.pressureT.val - wts.pressureT.val)/(P*1e-8))/abs(wts.pressureT.val), "\n"
    print "pressureD.val", wts.pressureD.val
    print "pressureD.ddT", wts.pressureD.ddT
    print "pressureD.ddP", wts.pressureD.ddP
    print "pressureD.ddTe", abs(wts.pressureD.ddT - (wtsT.pressureD.val - wts.pressureD.val)/(T*1e-8))/abs(wts.pressureD.val)
    print "pressureD.ddPe", abs(wts.pressureD.ddP - (wtsP.pressureD.val - wts.pressureD.val)/(P*1e-8))/abs(wts.pressureD.val), "\n"
    print "pressureTT.val", wts.pressureTT.val
    print "pressureTT.ddT", wts.pressureTT.ddT
    print "pressureTT.ddP", wts.pressureTT.ddP
    print "pressureTT.ddTe", abs(wts.pressureTT.ddT - (wtsT.pressureTT.val - wts.pressureTT.val)/(T*1e-8))/abs(wts.pressureTT.val)
    print "pressureTT.ddPe", abs(wts.pressureTT.ddP - (wtsP.pressureTT.val - wts.pressureTT.val)/(P*1e-8))/abs(wts.pressureTT.val), "\n"
    print "pressureTD.val", wts.pressureTD.val
    print "pressureTD.ddT", wts.pressureTD.ddT
    print "pressureTD.ddP", wts.pressureTD.ddP
    print "pressureTD.ddTe", abs(wts.pressureTD.ddT - (wtsT.pressureTD.val - wts.pressureTD.val)/(T*1e-8))/abs(wts.pressureTD.val)
    print "pressureTD.ddPe", abs(wts.pressureTD.ddP - (wtsP.pressureTD.val - wts.pressureTD.val)/(P*1e-8))/abs(wts.pressureTD.val), "\n"
    print "pressureDD.val", wts.pressureDD.val
    print "pressureDD.ddT", wts.pressureDD.ddT
    print "pressureDD.ddP", wts.pressureDD.ddP
    print "pressureDD.ddTe", abs(wts.pressureDD.ddT - (wtsT.pressureDD.val - wts.pressureDD.val)/(T*1e-8))/abs(wts.pressureDD.val)
    print "pressureDD.ddPe", abs(wts.pressureDD.ddP - (wtsP.pressureDD.val - wts.pressureDD.val)/(P*1e-8))/abs(wts.pressureDD.val), "\n"

printstate(298.15, 1e5)