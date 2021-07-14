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

from reaktoro import (
    Database,
    Element,
)

from pathlib import Path
import locale
import os
import pytest
import sys


def get_test_data_dir():
    return Path(os.path.abspath(__file__)).parents[0] / "data"


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


@pytest.fixture
def guard_locale():
    old_locale = locale.setlocale(locale.LC_NUMERIC)
    yield old_locale
    locale.setlocale(locale.LC_NUMERIC, old_locale)


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


def locale_has_comma_for_decimal_separator():
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


@pytest.mark.usefixtures("guard_locale")
def test_locale_problem_with_pugixml():
    locales = get_locales()
    at_least_one_locale_has_comma_for_decimal_separator = False
    at_least_one_valid_locale = False
    for loc in locales:
        if try_set_locale(loc):
            at_least_one_valid_locale = True
            if not at_least_one_locale_has_comma_for_decimal_separator:
                at_least_one_locale_has_comma_for_decimal_separator = locale_has_comma_for_decimal_separator()
            database = Database("supcrt07.xml")
            check_molar_mass({element.name(): element.molarMass() for element in database.elements()})

    assert at_least_one_valid_locale, "Couldn't find any valid locale to test"
    # assert at_least_one_locale_has_comma_for_decimal_separator, "Couldn't find any locale with comma for decimal separator"


def test_elements_molar_mass(guard_locale):
    old_locale = guard_locale
    database = Database("supcrt07.xml")
    assert old_locale == locale.setlocale(locale.LC_NUMERIC)
    check_molar_mass({element.name(): element.molarMass() for element in database.elements()})


def test_invariant_database():
    database = Database("supcrt07.xml")

    elements = {element.name(): element for element in database.elements()}
    assert len(elements) == 117

    elements['H'].setMolarMass(10)
    assert 10 == elements['H'].molarMass()

    elements = {element.name(): element for element in database.elements()}
    assert pytest.approx(0.00100794) == elements['H'].molarMass()


def test_database_instantiation_with_wrong_filename():
    with pytest.raises(RuntimeError):
        database = Database("wrong_name.xml")


def test_adding_and_getting_database_elements():
    database = Database(str(get_test_data_dir() / "supcrt98_simplified.xml"))

    new_element = Element()
    new_element.setName("He")
    new_element.setMolarMass(4.002602e-3)

    database.addElement(new_element)

    elements = database.elements()

    assert elements[0].name() == "Fe"
    assert elements[0].molarMass() == pytest.approx(55.845e-3)
    assert elements[1].name() == "H"
    assert elements[1].molarMass() == pytest.approx(1.00794e-3)
    assert elements[2].name() == "He"
    assert elements[2].molarMass() == pytest.approx(4.002602e-3)
    assert elements[3].name() == "S"
    assert elements[3].molarMass() == pytest.approx(32.065e-3)


def test_database_parse():
    """
    Test the fact that species should be added as
    liquid species even if the Type is Gaseous
    expected result:
    - liquid_species[0] = "H2S(liq)" -- added as liquid
    - gaseous_species[1] = "H2S(g)" -- added as gas
    """
    database = Database(str(get_test_data_dir() / "supcrt98_simplified.xml"))

    gaseous_species = database.gaseousSpecies()
    liquid_species = database.liquidSpecies()

    assert gaseous_species[0].name() == "H2S(g)"
    assert liquid_species[0].name() == "H2S(liq)"



def test_database_species_adding_and_getting():
    database = Database(str(get_test_data_dir() / "supcrt98_simplified.xml"))
    no_species_database = Database(str(get_test_data_dir() / "supcrt98_no_species.xml"))

    aqueous_species = database.aqueousSpecies()
    gaseous_species = database.gaseousSpecies()
    liquid_species = database.liquidSpecies()
    mineral_species = database.mineralSpecies()

    for aqueous_specie, gaseous_specie, mineral_specie, liquid_specie in  zip(aqueous_species, gaseous_species, mineral_species, liquid_species):
        no_species_database.addAqueousSpecies(aqueous_specie)
        no_species_database.addGaseousSpecies(gaseous_specie)
        no_species_database.addLiquidSpecies(liquid_specie)
        no_species_database.addMineralSpecies(mineral_specie)


    for aqueous_specie in aqueous_species:
        assert no_species_database.aqueousSpecies(aqueous_specie.name()).name() == aqueous_specie.name()

    for gaseous_specie in gaseous_species:
        assert no_species_database.gaseousSpecies(gaseous_specie.name()).name() == gaseous_specie.name()

    for liquid_specie in liquid_species:
        assert no_species_database.liquidSpecies(liquid_specie.name()).name() == liquid_specie.name()

    for mineral_specie in mineral_species:
        assert no_species_database.mineralSpecies(mineral_specie.name()).name() == mineral_specie.name()


def test_database_contains():
    database = Database(str(get_test_data_dir() / "supcrt98_simplified.xml"))

    aqueous_species = database.aqueousSpecies()
    gaseous_species = database.gaseousSpecies()
    liquid_species = database.liquidSpecies()
    mineral_species = database.mineralSpecies()

    assert database.containsAqueousSpecies(aqueous_species[0].name())
    assert database.containsGaseousSpecies(gaseous_species[0].name())
    assert database.containsLiquidSpecies(liquid_species[0].name())
    assert database.containsMineralSpecies(mineral_species[0].name())


def test_database_looking_for_species_with_element():
    database = Database(str(get_test_data_dir() / "supcrt98_simplified.xml"))

    aqueous_species_with_H_or_Fe = database.aqueousSpeciesWithElements(["H", "S"])
    gaseous_species_with_H_or_Fe = database.gaseousSpeciesWithElements(["H", "S"])
    liquid_species_with_H_or_Fe = database.liquidSpeciesWithElements(["H", "S"])
    mineral_species_with_H_or_Fe = database.mineralSpeciesWithElements(["Fe", "S"])

    assert aqueous_species_with_H_or_Fe[0].name() == "H2S(aq)"
    assert gaseous_species_with_H_or_Fe[0].name() == "H2S(g)"
    assert liquid_species_with_H_or_Fe[0].name() == "H2S(liq)"
    assert mineral_species_with_H_or_Fe[0].name() == "Pyrrhotite"

