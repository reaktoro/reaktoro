// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "CriticalProps.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/AggregateState.hpp>

namespace Reaktoro {
namespace detail {

const std::deque<SubstanceCriticalProps> preset_critical_props_data =
{
//     Tc/K     Pc/Pa     omega    unique identifiers (case insensitive)
    { {190.60,  45.99e5,  0.0120}, {"METHANE", "CH4"}           },
    { {305.30,  48.72e5,  0.1000}, {"ETHANE", "C2H6"}           },
    { {369.80,  42.48e5,  0.1520}, {"PROPANE", "C3H8"}          },
    { {425.10,  37.96e5,  0.2000}, {"N-BUTANE"}                 },
    { {469.70,  33.70e5,  0.2520}, {"N-PENTANE", "C5H12"}       },
    { {507.60,  30.25e5,  0.3010}, {"N-HEXANE", "C6H14"}        },
    { {540.20,  27.40e5,  0.3500}, {"N-HEPTANE"}                },
    { {568.70,  24.90e5,  0.4000}, {"N-OCTANE"}                 },
    { {594.60,  22.90e5,  0.4440}, {"N-NONANE"}                 },
    { {617.70,  21.10e5,  0.4920}, {"N-DECANE"}                 },
    { {408.10,  36.48e5,  0.1810}, {"ISOBUTANE"}                },
    { {544.00,  25.68e5,  0.3020}, {"ISOOCTANE"}                },
    { {511.80,  45.02e5,  0.1960}, {"CYCLOPENTANE"}             },
    { {553.60,  40.73e5,  0.2100}, {"CYCLOHEXANE", "C6H12"}     },
    { {532.80,  37.85e5,  0.2300}, {"METHYLCYCLOPENTANE"}       },
    { {572.20,  34.71e5,  0.2350}, {"METHYLCYCLOHEXANE"}        },
    { {282.30,  50.40e5,  0.0870}, {"ETHYLENE", "C2H4"}         },
    { {365.60,  46.65e5,  0.1400}, {"PROPYLENE", "C3H6"}        },
    { {420.00,  40.43e5,  0.1910}, {"1-BUTENE", "C4H8"}         },
    { {435.60,  42.43e5,  0.2050}, {"CIS-2-BUTENE"}             },
    { {428.60,  41.00e5,  0.2180}, {"TRANS-2-BUTENE"}           },
    { {504.00,  31.40e5,  0.2800}, {"1-HEXENE"}                 },
    { {417.90,  40.00e5,  0.1940}, {"ISOBUTYLENE"}              },
    { {425.20,  42.77e5,  0.1900}, {"1,3-BUTADIENE"}            },
    { {560.40,  43.50e5,  0.2120}, {"CYCLOHEXENE"}              },
    { {308.30,  61.39e5,  0.1870}, {"ACETYLENE", "C2H2"}        },
    { {562.20,  48.98e5,  0.2100}, {"BENZENE", "C6H6"}          },
    { {591.80,  41.06e5,  0.2620}, {"TOLUENE", "C7H8"}          },
    { {617.20,  36.06e5,  0.3030}, {"ETHYLBENZENE"}             },
    { {631.10,  32.09e5,  0.3260}, {"CUMENE"}                   },
    { {630.30,  37.34e5,  0.3100}, {"O-XYLENE"}                 },
    { {617.10,  35.36e5,  0.3260}, {"M-XYLENE"}                 },
    { {616.20,  35.11e5,  0.3220}, {"P-XYLENE"}                 },
    { {636.00,  38.40e5,  0.2970}, {"STYRENE"}                  },
    { {748.40,  40.51e5,  0.3020}, {"NAPHTHALENE"}              },
    { {789.30,  38.50e5,  0.3650}, {"BIPHENYL"}                 },
    { {408.00,  65.90e5,  0.2820}, {"FORMALDEHYDE"}             },
    { {466.00,  55.50e5,  0.2910}, {"ACETALDEHYDE"}             },
    { {506.60,  47.50e5,  0.3310}, {"METHYL-ACETATE"}           },
    { {523.30,  38.80e5,  0.3660}, {"ETHYL-ACETATE"}            },
    { {508.20,  47.01e5,  0.3070}, {"ACETONE"}                  },
    { {535.50,  41.50e5,  0.3230}, {"METHYL-ETHYL-KETONE"}      },
    { {466.70,  36.40e5,  0.2810}, {"DIETHYL-ETHER"}            },
    { {497.10,  34.30e5,  0.2660}, {"METHYL-T-BUTYL-ETHER"}     },
    { {512.60,  80.97e5,  0.5640}, {"METHANOL", "CH4O"}         },
    { {513.90,  61.48e5,  0.6450}, {"ETHANOL", "C2H6O"}         },
    { {536.80,  51.75e5,  0.6220}, {"1-PROPANOL"}               },
    { {563.10,  44.23e5,  0.5940}, {"1-BUTANOL"}                },
    { {611.40,  35.10e5,  0.5790}, {"1-HEXANOL"}                },
    { {508.30,  47.62e5,  0.6680}, {"2-PROPANOL"}               },
    { {694.30,  61.30e5,  0.4440}, {"PHENOL"}                   },
    { {719.70,  77.00e5,  0.4870}, {"ETHYLENE-GLYCOL"}          },
    { {592.00,  57.86e5,  0.4670}, {"ACETIC-ACID"}              },
    { {615.70,  40.64e5,  0.6810}, {"N-BUTYRIC-ACID"}           },
    { {751.00,  44.70e5,  0.6030}, {"BENZOIC-ACID"}             },
    { {545.50,  48.30e5,  0.3380}, {"ACETONITRILE"}             },
    { {430.10,  74.60e5,  0.2810}, {"METHYLAMINE"}              },
    { {456.20,  56.20e5,  0.2850}, {"ETHYLAMINE"}               },
    { {588.20,  63.10e5,  0.3480}, {"NITROMETHANE"}             },
    { {556.40,  45.60e5,  0.1930}, {"CARBON-TETRACHLORIDE"}     },
    { {536.40,  54.72e5,  0.2220}, {"CHLOROFORM"}               },
    { {510.00,  60.80e5,  0.1990}, {"DICHLOROMETHANE"}          },
    { {416.30,  66.80e5,  0.1530}, {"METHYL-CHLORIDE"}          },
    { {460.40,  52.70e5,  0.1900}, {"ETHYL-CHLORIDE"}           },
    { {632.40,  45.20e5,  0.2500}, {"CHLOROBENZENE"}            },
    { {374.20,  40.60e5,  0.3270}, {"TETRAFLUOROETHANE"}        },
    { {150.90,  48.98e5,  0.0000}, {"ARGON", "Ar"}              },
    { {209.40,  55.02e5,  0.0000}, {"KRYPTON", "Kr"}            },
    { {289.70,  58.40e5,  0.0000}, {"XENON", "Xe"}              },
    { {  5.20,   2.28e5, -0.3900}, {"HELIUM", "He"}             },
    { { 33.19,  13.13e5, -0.2160}, {"HYDROGEN", "H2"}           },
    { {154.60,  50.43e5,  0.0220}, {"OXYGEN", "O2"}             },
    { {126.20,  34.00e5,  0.0380}, {"NITROGEN", "N2"}           },
    { {132.20,  37.45e5,  0.0350}, {"AIR"}                      },
    { {417.20,  77.10e5,  0.0690}, {"CHLORINE", "Cl2"}          },
    { {132.90,  34.99e5,  0.0480}, {"CARBON-MONOXIDE", "CO"}    },
    { {304.20,  73.83e5,  0.2240}, {"CARBON-DIOXIDE", "CO2"}    },
    { {552.00,  79.00e5,  0.1110}, {"CARBON-DISULFIDE", "CS2"}  },
    { {373.50,  89.63e5,  0.0940}, {"HYDROGEN-SULFIDE", "H2S"}  },
    { {430.80,  78.84e5,  0.2450}, {"SULFUR-DIOXIDE", "SO2"}    },
    { {490.90,  82.10e5,  0.4240}, {"SULFUR-TRIOXIDE", "SO3"}   },
    { {180.20,  64.80e5,  0.5830}, {"NITRIC-OXIDE", "NO"}       },
    { {309.60,  72.45e5,  0.1410}, {"NITROUS-OXIDE", "N2O"}     },
    { {324.70,  83.10e5,  0.1320}, {"HYDROGEN-CHLORIDE", "HCl"} },
    { {456.70,  53.90e5,  0.4100}, {"HYDROGEN-CYANIDE", "HCN"}  },
    { {647.10, 220.55e5,  0.3450}, {"WATER", "H2O"}             },
    { {405.70, 112.80e5,  0.2530}, {"AMMONIA", "NH3"}           },
    { {520.00,  68.90e5,  0.7140}, {"NITRIC-ACID", "HNO3"}      },
    { {924.00,  64.00e5,  0.0000}, {"SULFURIC-ACID", "H2SO4"}   },
    { {377.00,  62.80e5,  0.0000}, {"RADON", "Ra"}              },
    { { 44.40 , 27.60e5,  0.0000}, {"NEON", "Ne"}               },
};

/// Return a substance name corrected for checking in the ChemicalProps database.
/// This function removes suffix from given substance name. For
/// example, `H2O(aq)` becomes `H2O`. It also replaces spaces by dashes and
/// transforms all characters to upper case, converting, for example, `carbon
/// dioxide` into `CARBON-DIOXIDE`.
auto correctName(String name) -> String
{
    const auto [name0, _] = splitSpeciesNameSuffix(name);
    name = name0;
    name = replace(name, " ", "-");
    name = uppercase(name);
    return name;
}

/// Return a vector of substance names corrected for checking in the ChemicalProps database.
auto correctNames(Strings names) -> Strings
{
    for(auto& name : names)
        name = correctName(name);
    return names;
}

} // namespace detail

SubstanceCriticalProps::SubstanceCriticalProps(const StringList& names)
: SubstanceCriticalProps({}, names)
{}

SubstanceCriticalProps::SubstanceCriticalProps(const SubstanceCriticalPropsData& data, const StringList& names)
: m_data(std::move(data)), m_names(detail::correctNames(names))
{}

auto SubstanceCriticalProps::setTemperature(real value, String unit) -> void
{
    value = units::convert(value, unit, "K");
    errorif(value <= 0.0, "Cannot set non-positive critical temperature value (", value, " K) to substance with name identifiers ", names(), ".");
    m_data.Tcr = value;
}

auto SubstanceCriticalProps::setPressure(real value, String unit) -> void
{
    errorif(value <= 0.0, "Cannot set non-positive critical pressure value (", value, " ", unit, ") to substance with name identifiers ", names(), ".");
    m_data.Pcr = units::convert(value, unit, "Pa");
}

auto SubstanceCriticalProps::setAcentricFactor(real value) -> void
{
    m_data.omega = value;
}

auto SubstanceCriticalProps::names() const -> const Strings&
{
    return m_names;
}

auto SubstanceCriticalProps::temperature() const -> real
{
    return m_data.Tcr;
}

auto SubstanceCriticalProps::pressure() const -> real
{
    return m_data.Pcr;
}

auto SubstanceCriticalProps::acentricFactor() const -> real
{
    return m_data.omega;
}

auto SubstanceCriticalProps::data() const -> const SubstanceCriticalPropsData&
{
    return m_data;
}

SubstanceCriticalProps::operator SubstanceCriticalPropsData() const
{
    return data();
}

CriticalProps::CriticalProps()
: m_data(detail::preset_critical_props_data)
{}

CriticalProps::~CriticalProps()
{}

auto CriticalProps::instance() -> CriticalProps&
{
    static CriticalProps obj;
    return obj;
}

auto CriticalProps::data() -> const std::deque<SubstanceCriticalProps>&
{
    return instance().m_data;
}

auto CriticalProps::defaultCriticalProps() -> const Optional<SubstanceCriticalProps>&
{
    return instance().m_default_crprops;
}

auto CriticalProps::append(SubstanceCriticalProps substance) -> void
{
    // Ensure there are no equivalent substances in the database (same name or same aliases).
    auto& substances = instance().m_data;
    for(auto& current : substances)
    {
        const auto has_common_name = !disjoint(substance.names(), current.names());

        errorif(has_common_name,
            "Appending critical property data for substance with names {", substance.names(), "}.\n"
            "However, one of these names conflic with one or more names of another already\n"
            "stored substance with names {", current.names(), "}.\n"
            "Use CriticalProps::overwrite instead of CriticalProps::append if you want to\n"
            "force append and overwrite.");
    }
    substances.push_back(substance);
}

auto CriticalProps::overwrite(SubstanceCriticalProps substance) -> void
{
    // Ensure there are no equivalent substances in the database (same name or same aliases).
    auto& substances = instance().m_data;
    for(auto& current : substances)
    {
        const auto has_common_name = !disjoint(substance.names(), current.names());
        if(has_common_name)
        {
            current = substance;
            return;
        }
    }
    substances.push_back(substance);
}

auto CriticalProps::setMissingAs(const String& substance) -> void
{
    const auto idx = find(substance);
    errorif(idx >= size(),
        "CriticalProps::setMissingAs requires an existing substance "
        "in the CriticalProps database, and ", substance, " is not present in it.");
    instance().m_default_crprops = data()[idx];
}

auto CriticalProps::size() -> Index
{
    return data().size();
}

auto CriticalProps::find(const String& substance) -> Index
{
    const auto name = detail::correctName(substance);
    return indexfn(data(), [&](auto&& s) { return contains(s.names(), name); });
}

auto CriticalProps::get(const String& substance) -> Optional<SubstanceCriticalProps>
{
    const auto idx = find(substance);
    if(idx < size()) return data()[idx];
    return defaultCriticalProps();
}

auto CriticalProps::get(const StringList& substances) -> Optional<SubstanceCriticalProps>
{
    for(auto&& substance : substances)
        if(const auto subs = get(substance); subs)
            return subs;
    return defaultCriticalProps();
}

} // namespace Reaktoro
