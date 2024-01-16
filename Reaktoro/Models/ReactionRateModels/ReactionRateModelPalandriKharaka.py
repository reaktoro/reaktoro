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


def testReactionRateModelPalandriKharaka():
    params = Params.embedded("PalandriKharaka.yaml")

    # GeneralReaction("Calcite").setRateModel(ReactionRateModelPalandriKharaka(params))

    db = SupcrtDatabase("supcrtbl")

    def areafn(props: ChemicalProps):
        return 1.0

    system = ChemicalSystem(db,
        AqueousPhase("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)").setActivityModel(ActivityModelDavies()),
        MineralPhase("Calcite"),
        GeneralReaction("Calcite").setRateModel(ReactionRateModelPalandriKharaka(params)),
        Surface("Calcite").withAreaModel(areafn)
    )
