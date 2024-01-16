# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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


def testSurface():

    props = ChemicalProps()

    def areafn(props: ChemicalProps):
        return 1.23

    # When using constructor Surface()
    surface = Surface()

    surface = surface.withName("AqueousPhase:GaseousPhase")
    assert surface.name() == "AqueousPhase:GaseousPhase"

    surface = surface.withAreaModel(areafn)
    assert surface.area(props) == 1.23

    # When using constructor Surface(name)
    surface = Surface("Quartz")

    assert surface.name() == "Quartz"

    # When using constructor Surface(name, model)
    surface = Surface("Calcite", areafn)

    assert surface.name() == "Calcite"
    assert surface.area(props) == 1.23
