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
import subprocess

def check_molar_mass(molar_mass_map):
    assert pytest.approx(0.00100794) == molar_mass_map['H']
    assert pytest.approx(0.0120107) == molar_mass_map['C']
    assert pytest.approx(0.040078) == molar_mass_map['Ca']
    assert pytest.approx(0.0159994) == molar_mass_map['O']
    assert pytest.approx(0.03206499999999) == molar_mass_map['S']
    assert pytest.approx(0.02298977) == molar_mass_map['Na']
    assert pytest.approx(0.035453) == molar_mass_map['Cl']
    assert pytest.approx(0.08762) == molar_mass_map['Sr']

def find_locales():
    out = subprocess.run(['locale', '-a'], stdout=subprocess.PIPE).stdout
    try:
        res = out.decode('utf-8')
    except:
        res = out.decode('latin-1')
    return res.rstrip('\n').splitlines()

def set_locale():
    available_locales = find_locales()
    for loc in available_locales:
        locale.setlocale(locale.LC_NUMERIC, loc)
        convention = locale.localeconv()
        if convention['decimal_point'] == ',':
            return True
    return False

def test_elements_molar_mass():
    database = Database("supcrt07.xml")
    check_molar_mass({element.name(): element.molarMass() for element in database.elements()})

def test_locale_problem_with_pugixml():
    old_locale = locale.setlocale(locale.LC_NUMERIC)
    if not set_locale():
        locale.setlocale(locale.LC_NUMERIC, old_locale)
        pytest.skip("Couldn't find any valid locale to test")

    database = Database("supcrt07.xml")
    locale.setlocale(locale.LC_NUMERIC, old_locale)
    check_molar_mass({element.name(): element.molarMass() for element in database.elements()})

def test_invariant_database():
    database = Database("supcrt07.xml")

    elements = {element.name(): element for element in database.elements()}
    elements['H'].setMolarMass(10)
    assert 10 == elements['H'].molarMass()

    elements = {element.name(): element for element in database.elements()}
    assert pytest.approx(0.00100794) == elements['H'].molarMass()
