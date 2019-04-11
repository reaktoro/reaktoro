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

from reaktoro import Database
import pytest

import locale

def test_elements_molar_mass():
    database = Database("supcrt07.xml")
    molar_mass_map = {element.name(): element.molarMass() for element in database.elements()}
    assert pytest.approx(0.00100794) == molar_mass_map['H']
    assert pytest.approx(0.0120107) == molar_mass_map['C']
    assert pytest.approx(0.040078) == molar_mass_map['Ca']
    assert pytest.approx(0.0159994) == molar_mass_map['O']
    assert pytest.approx(0.03206499999999) == molar_mass_map['S']
    assert pytest.approx(0.02298977) == molar_mass_map['Na']
    assert pytest.approx(0.035453) == molar_mass_map['Cl']
    assert pytest.approx(0.08762) == molar_mass_map['Sr']

def test_locale_problem_with_pugixml():
    old_locale = locale.setlocale(locale.LC_NUMERIC)
    locale.setlocale(locale.LC_NUMERIC, 'pt_BR.utf8')

    database = Database("supcrt07.xml")
    molar_mass_map = {element.name(): element.molarMass() for element in database.elements()}
    assert 0.001007 > molar_mass_map['H']
    assert 0.01201 > molar_mass_map['C']
    assert 0.04007 > molar_mass_map['Ca']
    assert 0.0159 > molar_mass_map['O']
    assert 0.03206 > molar_mass_map['S']
    assert 0.0229 > molar_mass_map['Na']
    assert 0.0354 > molar_mass_map['Cl']
    assert 0.0876 > molar_mass_map['Sr']

    locale.setlocale(locale.LC_NUMERIC, old_locale)

def test_invariant_database():
    database = Database("supcrt07.xml")

    elements = {element.name(): element for element in database.elements()}
    elements['H'].setMolarMass(10)
    assert 10 == elements['H'].molarMass()

    elements = {element.name(): element for element in database.elements()}
    assert pytest.approx(0.00100794) == elements['H'].molarMass()
