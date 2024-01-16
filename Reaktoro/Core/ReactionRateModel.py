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


def testReactionRateModel():

    db = SupcrtDatabase("supcrtbl")

    phases = Phases(db)
    phases.add( AqueousPhase("H2O(aq)") )

    system = ChemicalSystem(phases)
    props = ChemicalProps(system)

    # Define a rate function that returns a ReactionRate object
    def ratefn1(props: ChemicalProps) -> ReactionRate:
        return ReactionRate(1.23)

    # Define a rate function that returns a float value
    def ratefn2(props: ChemicalProps) -> float:
        return 2.34

    model1 = ReactionRateModel(ratefn1)
    model2 = ReactionRateModel(ratefn2)

    assert model1(props).value() == 1.23
    assert model2(props).value() == 2.34
