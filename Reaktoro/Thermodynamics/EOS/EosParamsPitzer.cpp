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

#include "EosParamsPitzer.hpp"

// yaml-cpp includes
#include <yaml-cpp/yaml.h>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace {

auto parseBlockBinaryMultiParams(std::string block) -> std::vector<std::tuple<std::string, std::string, std::vector<double>>>
{
    std::vector<std::tuple<std::string, std::string, std::vector<double>>> table;

    auto lines = split(block, "\n");

    for(const auto& line : lines)
    {
        auto words = split(line, " ");

        Assert(words.size() > 2, "Cannot parse the Pitzer database.",
            "Expecting two species names plus at least one parameter in the line `" + line + "`.");

        std::string species1 = words[0];
        std::string species2 = words[1];

        std::vector<double> values;
        for(unsigned i = 2; i < words.size(); ++i)
            values.push_back(tofloat(words[i]));

        table.push_back(std::make_tuple(species1, species2, values));
    }

    return table;
}

auto parseBlockBinarySingleParam(std::string block) -> std::vector<std::tuple<std::string, std::string, double>>
{
    std::vector<std::tuple<std::string, std::string, double>> table;

    auto lines = split(block, "\n");

    for(const auto& line : lines)
    {
        auto words = split(line, " ");

        Assert(words.size() == 3, "Cannot parse the Pitzer database.",
            "Expecting two species names plus one parameter in the line `" + line + "`.");

        std::string species1 = words[0];
        std::string species2 = words[1];
        double value = tofloat(words[2]);

        table.push_back(std::make_tuple(species1, species2, value));
    }

    return table;
}

auto parseBlockTernarySingleParam(std::string block) -> std::vector<std::tuple<std::string, std::string, std::string, double>>
{
    std::vector<std::tuple<std::string, std::string, std::string, double>> table;

    auto lines = split(block, "\n");

    for(const auto& line : lines)
    {
        auto words = split(line, " ");

        Assert(words.size() == 4, "Cannot parse the Pitzer database.",
            "Expecting three species names plust one parameter in the line `" + line + "`.");

        std::string species1 = words[0];
        std::string species2 = words[1];
        std::string species3 = words[2];
        double value = tofloat(words[3]);

        table.push_back(std::make_tuple(species1, species2, species3, value));
    }

    return table;
}

//auto theta(std::string ion1, std::string ion2, std::string theta_block) -> double
//{
//    std::set<std::string> ions = {ion1, ion2};
//
//    auto lines = split(theta_block, "\n");
//
//    for(const auto& line : lines)
//    {
//        std::vector<std::string> words = split(line, " ");
//        std::set<std::string> names(words.begin(), words.begin() + 2);
//
//        if(ions == names)
//            return tofloat(words[2]);
//    }
//
//    return 0.0;
//}
//
//auto psi(std::string ion1, std::string ion2, std::string ion3, std::string psi_block) -> double
//{
//    std::set<std::string> ions = {ion1, ion2, ion3};
//
//    auto lines = split(psi_block, "\n");
//
//    for(const auto& line : lines)
//    {
//        std::vector<std::string> words = split(line, " ");
//        std::set<std::string> names(words.begin(), words.begin() + 3);
//
//        if(ions == names)
//            return tofloat(words[3]);
//    }
//
//    return 0.0;
//}
//
//auto lambda(std::string neutral, std::string ion, std::string lambda_block) -> double
//{
//    std::set<std::string> species = {neutral, ion};
//
//    auto lines = split(lambda_block, "\n");
//
//    for(const auto& line : lines)
//    {
//        std::vector<std::string> words = split(line, " ");
//        std::set<std::string> names(words.begin(), words.begin() + 2);
//
//        if(species == names)
//            return tofloat(words[2]);
//    }
//
//    return 0.0;
//}
//
//auto zeta(std::string neutral, std::string cation, std::string anion, std::string zeta_block) -> double
//{
//    std::set<std::string> species = {neutral, cation, anion};
//
//    auto lines = split(zeta_block, "\n");
//
//    for(const auto& line : lines)
//    {
//        std::vector<std::string> words = split(line, " ");
//        std::set<std::string> names(words.begin(), words.begin() + 3);
//
//        if(species == names)
//            return tofloat(words[3]);
//    }
//
//    return 0.0;
//}
//
//auto createSingleSaltParamFunction(std::string cation, std::string anion, std::string block) -> std::function<double(double)>
//{
//    // Split the given block of interaction parameters in lines of text to be parsed below
//    auto lines = split(block, "\n");
//
//    // Iterate over all lines of data and find the one with the pair cation and anion
//    for(const auto& line : lines)
//    {
//        auto words = split(line, " ");
//
//        if(cation == words[0] && anion == words[1])
//        {
//            const double Tr = 298.15;
//
//            std::vector<double> c(words.size() - 2);
//
//            for(unsigned i = 0; i < c.size(); ++i)
//                c[i] = tofloat(words[i + 2]);
//
//            if(c.size() == 1)
//                return [=](double T) { return c[0]; };
//
//            if(c.size() == 2)
//                return [=](double T) { return c[0] + c[1]*(T - Tr); };
//
//            if(c.size() == 5)
//                return [=](double T) { return c[0] + c[1]*(1/T - 1/Tr) + c[2]*std::log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr); };
//
//            RuntimeError("Cannot create the single salt parameter function of Pitzer model.",
//                "The number of coefficients for the equation is not supported");
//        }
//    }
//
//    // Return a zero function in case the pair cation and anion does not have Pitzer data
//    return [=](double T) { return 0.0; };
//}
//
//auto createSingleSaltParamTable(const std::vector<std::string>& cations, const std::vector<std::string>& anions, std::string block) -> Table2D<std::function<double(double)>>
//{
//    Table2D<std::function<double(double)>> table = table2D<std::function<double(double)>>(cations.size(), anions.size());
//
//    for(unsigned i = 0; i < cations.size(); ++i)
//        for(unsigned j = 0; j < anions.size(); ++j)
//            table[i][j] = createSingleSaltParamFunction(cations[i], anions[j], block);
//
//    return table;
//}
//
//auto createBeta0Table(const std::vector<std::string>& cations, const std::vector<std::string>& anions, std::string beta0_block) -> Table2D<std::function<double(double)>>
//{
//    return createSingleSaltParamTable(cations, anions, beta0_block);
//}
//
//auto createBeta1Table(const std::vector<std::string>& cations, const std::vector<std::string>& anions, std::string beta1_block) -> Table2D<std::function<double(double)>>
//{
//    return createSingleSaltParamTable(cations, anions, beta1_block);
//}
//
//auto createBeta2Table(const std::vector<std::string>& cations, const std::vector<std::string>& anions, std::string beta2_block) -> Table2D<std::function<double(double)>>
//{
//    return createSingleSaltParamTable(cations, anions, beta2_block);
//}
//
//auto createCphiTable(const std::vector<std::string>& cations, const std::vector<std::string>& anions, std::string cphi_block) -> Table2D<std::function<double(double)>>
//{
//    return createSingleSaltParamTable(cations, anions, cphi_block);
//}
//
//auto createThetaTable(const std::vector<std::string>& ions1, const std::vector<std::string>& ions2, std::string theta_block) -> Table2D<double>
//{
//    Table2D<double> table = table2D<double>(ions1.size(), ions2.size());
//
//    for(unsigned i = 0; i < ions1.size(); ++i)
//        for(unsigned j = 0; j < ions2.size(); ++j)
//            table[i][j] = theta(ions1[i], ions2[j], theta_block);
//
//    return table;
//}
//
//auto createPsiTable(const std::vector<std::string>& ions1, const std::vector<std::string>& ions2, const std::vector<std::string>& ions3, std::string psi_block) -> Table3D<double>
//{
//    Table3D<double> table = table3D<double>(ions1.size(), ions2.size(), ions3.size());
//
//    for(unsigned i = 0; i < ions1.size(); ++i)
//        for(unsigned j = 0; j < ions2.size(); ++j)
//            for(unsigned k = 0; k < ions3.size(); ++k)
//                table[i][j][k] = psi(ions1[i], ions2[j], ions3[k], psi_block);
//
//    return table;
//}
//
//auto createLambdaTable(const std::vector<std::string>& neutrals, const std::vector<std::string>& ions, std::string lambda_block) -> Table2D<double>
//{
//    Table2D<double> table = table2D<double>(neutrals.size(), ions.size());
//
//    for(unsigned i = 0; i < neutrals.size(); ++i)
//        for(unsigned j = 0; j < ions.size(); ++j)
//            table[i][j] = lambda(neutrals[i], ions[j], lambda_block);
//
//    return table;
//}
//
//auto createZetaTable(const std::vector<std::string>& neutrals, const std::vector<std::string>& cations, const std::vector<std::string>& anions, std::string zeta_block) -> Table3D<double>
//{
//    Table3D<double> table = table3D<double>(neutrals.size(), cations.size(), anions.size());
//
//    for(unsigned i = 0; i < neutrals.size(); ++i)
//        for(unsigned j = 0; j < cations.size(); ++j)
//            for(unsigned k = 0; k < anions.size(); ++k)
//                table[i][j][k] = zeta(neutrals[i], cations[j], anions[k], zeta_block);
//
//    return table;
//}

} // namespace

EosParamsPitzer::EosParamsPitzer()
{

}

EosParamsPitzer::EosParamsPitzer(std::string filename)
{
    YAML::Node doc = YAML::LoadFile(filename);

    auto beta0_node  = doc["Beta0"];
    auto beta1_node  = doc["Beta1"];
    auto beta2_node  = doc["Beta2"];
    auto cphi_node   = doc["Cphi"];
    auto theta_node  = doc["Theta"];
    auto lambda_node = doc["Lambda"];
    auto psi_node    = doc["Psi"];
    auto zeta_node   = doc["Zeta"];

    Assert(beta0_node,  "Cannot parse the Pitzer database.", "Expecting a `Beta0` data block in the Pitzer database.");
    Assert(beta1_node,  "Cannot parse the Pitzer database.", "Expecting a `Beta1` data block in the Pitzer database.");
    Assert(beta2_node,  "Cannot parse the Pitzer database.", "Expecting a `Beta2` data block in the Pitzer database.");
    Assert(cphi_node,   "Cannot parse the Pitzer database.", "Expecting a `Cphi` data block in the Pitzer database.");
    Assert(theta_node,  "Cannot parse the Pitzer database.", "Expecting a `Theta` data block in the Pitzer database.");
    Assert(lambda_node, "Cannot parse the Pitzer database.", "Expecting a `Lambda` data block in the Pitzer database.");
    Assert(psi_node,    "Cannot parse the Pitzer database.", "Expecting a `Psi` data block in the Pitzer database.");
    Assert(zeta_node,   "Cannot parse the Pitzer database.", "Expecting a `Zeta` data block in the Pitzer database.");

    beta0  = parseBlockBinaryMultiParams(beta0_node.as<std::string>());
    beta1  = parseBlockBinaryMultiParams(beta1_node.as<std::string>());
    beta2  = parseBlockBinaryMultiParams(beta2_node.as<std::string>());
    cphi   = parseBlockBinaryMultiParams(cphi_node.as<std::string>());
    theta  = parseBlockBinarySingleParam(theta_node.as<std::string>());
    lambda = parseBlockBinarySingleParam(lambda_node.as<std::string>());
    psi    = parseBlockTernarySingleParam(psi_node.as<std::string>());
    zeta   = parseBlockTernarySingleParam(zeta_node.as<std::string>());
}

} // namespace Reaktoro
