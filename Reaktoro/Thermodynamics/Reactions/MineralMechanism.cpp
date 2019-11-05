// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "MineralMechanism.hpp"

// C++ includes
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {
namespace internal {

inline auto unknownOptionError(std::string option) -> void
{
    Exception exception;
    exception.error << "Cannot set the option " << option << " in the mineral mechanism.";
    exception.reason << "This option has incorrect format or is not supported.";
    RaiseError(exception);
}

inline auto missingUnitError(std::string quantity) -> void
{
    Exception exception;
    exception.error << "Cannot set the quantity " << quantity << " in the mineral mechanism.";
    exception.reason << "The units of quantity " << quantity << " have not been specified.";
    RaiseError(exception);
}

inline auto checkRateConstantUnit(std::string unit) -> void
{
    if(!units::convertible(unit, "mol/(m2*s)"))
        RuntimeError("Cannot set the kinetic rate constant of the mineral reaction",
                     "The provided unit cannot be converted to mol/(m2*s)");
}

inline auto checkActivationEnergyUnit(std::string unit) -> void
{
    if(!units::convertible(unit, "kJ/mol"))
        RuntimeError("Cannot set the Arrhenius activation energy of the mineral reaction",
                     "The provided unit cannot be converted to kJ/mol");
}

} /* namespace internal */

using namespace internal;

MineralMechanism::MineralMechanism()
{}

MineralMechanism::MineralMechanism(std::string mechanism)
{
    const auto options = split(mechanism, ";");

    for(const auto& option : options) {
        if(option.find("a[") != std::string::npos || option.find("activity[") != std::string::npos ||
           option.find("p[") != std::string::npos || option.find("pressure[") != std::string::npos) {
            catalysts.push_back(MineralCatalyst(option));
            continue;
        }

        const auto words = split(option, "= ");
        if(words.empty() || words.size() > 3)
            unknownOptionError(option);

        const auto& quantity = words[0];
        const auto& value = tofloat(words[1]);

        if(quantity == "logk" || quantity == "Ea")
            if(words.size() < 3)
                missingUnitError(quantity);

        const auto& unit = (words.size() == 3) ? words[2] : "";

        if(quantity == "logk")
            setRateConstant(std::pow(10.0, value), unit);
        else if(quantity == "Ea")
            setActivationEnergy(value, unit);
        else if(quantity == "p")
            p = value;
        else if(quantity == "q")
            q = value;
        else
            unknownOptionError(option);
    }
}

auto MineralMechanism::setRateConstant(double value, std::string unit) -> MineralMechanism&
{
    checkRateConstantUnit(unit);
    kappa = units::convert(value, unit, "mol/(m2*s)");
    return *this;
}

auto MineralMechanism::setActivationEnergy(double value, std::string unit) -> MineralMechanism&
{
    checkActivationEnergyUnit(unit);
    Ea = units::convert(value, unit, "kJ/mol");
    return *this;
}

auto MineralMechanism::setPowerP(double value) -> MineralMechanism&
{
    p = value;
    return *this;
}

auto MineralMechanism::setPowerQ(double value) -> MineralMechanism&
{
    q = value;
    return *this;
}

auto MineralMechanism::setCatalysts(std::string strcatalysts) -> MineralMechanism&
{
    catalysts.clear();
    catalysts.push_back(MineralCatalyst(strcatalysts));
    return *this;
}

auto MineralMechanism::setCatalysts(const MineralCatalyst& catalyst) -> MineralMechanism&
{
    catalysts.clear();
    catalysts.push_back(catalyst);
    return *this;
}

auto MineralMechanism::setCatalysts(const std::vector<MineralCatalyst>& veccatalysts) -> MineralMechanism&
{
    catalysts = veccatalysts;
    return *this;
}

auto operator<(const MineralMechanism& lhs, const MineralMechanism& rhs) -> bool
{
    return lhs.kappa < rhs.kappa;
}

auto operator==(const MineralMechanism& lhs, const MineralMechanism& rhs) -> bool
{
    return lhs.kappa == rhs.kappa;
}

} // namespace Reaktoro
