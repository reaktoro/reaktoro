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

#include "Gnuplot.hpp"

// C++ includes
#include <stdio.h>

// Ensure appropriate popen or pclose calls when compiling with MSVC
#ifdef _MSC_VER
#define popen _popen
#define pclose _pclose
#endif

namespace Reaktoro {
namespace {

const std::string style = R"(
#--------------------------------------
# Macros
#--------------------------------------
set macros

# Line width
LW = "2"

#--------------------------------------
# Line Styles
#--------------------------------------
set style line 1 lc rgb '#1B9E77' # dark teal
set style line 2 lc rgb '#D95F02' # dark orange
set style line 3 lc rgb '#7570B3' # dark lilac
set style line 4 lc rgb '#E7298A' # dark magenta
set style line 5 lc rgb '#66A61E' # dark lime green
set style line 6 lc rgb '#E6AB02' # dark banana
set style line 7 lc rgb '#A6761D' # dark tan
set style line 8 lc rgb '#666666' # dark gray

#--------------------------------------
# Pallete
#--------------------------------------
set palette maxcolors 8
set palette defined ( 0 '#1B9E77',\
                      1 '#D95F02',\
                      2 '#7570B3',\
                      3 '#E7298A',\
                      4 '#66A61E',\
                      5 '#E6AB02',\
                      6 '#A6761D',\
                      7 '#666666' )

#--------------------------------------
# Grid Style
#--------------------------------------
set style line 102 lc rgb '#d6d7d9' lt 0 lw @LW
set grid back ls 102

#--------------------------------------
# Border Style
#--------------------------------------
set style line 101 lc rgb '#808080' lt 1 lw @LW
set border front ls 101
set tics nomirror out scale 0.75
set format '%g'

)";

} // namespace

Gnuplot::Gnuplot()
    : pipe(popen("gnuplot -persist", "w"))
{
    *this << style;
}

Gnuplot::~Gnuplot()
{
    pclose(pipe);
}

auto Gnuplot::operator<<(std::string str) -> Gnuplot&
{
    str.append("\n");
    fputs(str.c_str(), pipe);
    return *this;
}

auto Gnuplot::operator<<(const std::stringstream& ss) -> Gnuplot&
{
    fputs(ss.str().c_str(), pipe);
    return *this;
}

} // namespace Reaktoro
