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
import sys

LANGUAGES = {
    'bg_BG': 'Bulgarian',
    'cs_CZ': 'Czech',
    'da_DK': 'Danish',
    'de_DE': 'German',
    'el_GR': 'Greek',
    'en_US': 'English',
    'es_ES': 'Spanish',
    'et_EE': 'Estonian',
    'fi_FI': 'Finnish',
    'fr_FR': 'French',
    'hr_HR': 'Croatian',
    'hu_HU': 'Hungarian',
    'it_IT': 'Italian',
    'lt_LT': 'Lithuanian',
    'lv_LV': 'Latvian',
    'nl_NL': 'Dutch',
    'no_NO': 'Norwegian',
    'pl_PL': 'Polish',
    'pt_PT': 'Portuguese',
    'ro_RO': 'Romanian',
    'ru_RU': 'Russian',
    'sk_SK': 'Slovak',
    'sl_SI': 'Slovenian',
    'sv_SE': 'Swedish',
    'tr_TR': 'Turkish',
    'zh_CN': 'Chinese',
}


def try_set_locale(loc):
    try:
        locale.setlocale(locale.LC_NUMERIC, loc)
        return True
    except:
        return False


def get_linux_locales():
    locales = list(LANGUAGES.keys())
    locales += [loc + '.utf8' for loc in locales]
    return locales


def get_locales():
    system = sys.platform
    if 'win32' in system:
        return list(LANGUAGES.values())

    return get_linux_locales()


def locale_has_comma_for_decimal_point():
    return locale.localeconv()['decimal_point'] == ','


def check_molar_mass(molar_mass_map):
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

    locales = get_locales()
    at_least_one_locale_has_comma_for_decimal_separator = False
    at_least_one_valid_locale = False
    for loc in locales:
        if try_set_locale(loc):
            at_least_one_valid_locale = True
            locale.setlocale(locale.LC_NUMERIC, loc)
            if not at_least_one_locale_has_comma_for_decimal_separator:
                at_least_one_locale_has_comma_for_decimal_separator = locale_has_comma_for_decimal_point()
            database = Database("supcrt07.xml")
            locale.setlocale(locale.LC_NUMERIC, old_locale)
            check_molar_mass({element.name(): element.molarMass() for element in database.elements()})

    assert at_least_one_valid_locale, "Couldn't find any valid locale to test"
    assert at_least_one_locale_has_comma_for_decimal_separator, "Couldn't find any locale with comma fot decimal separator"


def test_elements_molar_mass():
    old_locale = locale.setlocale(locale.LC_NUMERIC)
    database = Database("supcrt07.xml")
    assert old_locale == locale.setlocale(locale.LC_NUMERIC)
    check_molar_mass({element.name(): element.molarMass() for element in database.elements()})


def test_invariant_database():
    database = Database("supcrt07.xml")

    elements = {element.name(): element for element in database.elements()}
    elements['H'].setMolarMass(10)
    assert 10 == elements['H'].molarMass()

    elements = {element.name(): element for element in database.elements()}
    assert pytest.approx(0.00100794) == elements['H'].molarMass()
