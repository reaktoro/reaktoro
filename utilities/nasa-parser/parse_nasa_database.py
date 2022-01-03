# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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


class NasaPolynomial:
    """The coefficients and integration constants for the NASA polynomial model
    for standard thermodynamic properties calculation. These coefficients and
    constants are documented in [1] (Section 4.2, page 19) and in [2] (Appendix
    A, page 74).
       1. Mcbride, B.J., Gordon, S., Mcbride, B.J. (1994). Computer program for
          calculation of complex chemical equilibrium compositions and
          applications. I: Analysis. In NASA Reference Publication 1311.
          https://doi.org/NASA RP-1311
       2. McBride, B.J., Gordon, S. (1996). Computer Program for Calculation of
          Complex Chemical Equilibrium Compositions and Applications: II-User
          Manual and Program Description. In NASA Reference Publication 1311.
          https://doi.org/NASA RP-1311
    """
    Tmin:float # The minimum temperature (in K) for which the coefficients below are valid.
    Tmax:float # The maximum temperature (in K) for which the coefficients below are valid.
    qN:int     # The number of exponent coefficients qi, which is always 7.
    q1:float   # The exponent q1 in the regression model for Cp0.
    q2:float   # The exponent q2 in the regression model for Cp0.
    q3:float   # The exponent q3 in the regression model for Cp0.
    q4:float   # The exponent q4 in the regression model for Cp0.
    q5:float   # The exponent q5 in the regression model for Cp0.
    q6:float   # The exponent q6 in the regression model for Cp0.
    q7:float   # The exponent q7 in the regression model for Cp0.
    a1:float   # The least-square coefficient a1 in the regression model for Cp0.
    a2:float   # The least-square coefficient a2 in the regression model for Cp0.
    a3:float   # The least-square coefficient a3 in the regression model for Cp0.
    a4:float   # The least-square coefficient a4 in the regression model for Cp0.
    a5:float   # The least-square coefficient a5 in the regression model for Cp0.
    a6:float   # The least-square coefficient a6 in the regression model for Cp0.
    a7:float   # The least-square coefficient a7 in the regression model for Cp0.
    b1:float   # The integration constant b1 used to compute H0.
    b2:float   # The integration constant b2 used to compute S0.


    def __repr__(self):
        """Useful for comparing data during pytest failure"""
        d = {}
        d["Tmin"] = self.Tmin
        d["Tmax"] = self.Tmax
        d["qN"] = self.qN
        d["q1"] = self.q1
        d["q2"] = self.q2
        d["q3"] = self.q3
        d["q4"] = self.q4
        d["q5"] = self.q5
        d["q6"] = self.q6
        d["q7"] = self.q7
        d["a1"] = self.a1
        d["a2"] = self.a2
        d["a3"] = self.a3
        d["a4"] = self.a4
        d["a5"] = self.a5
        d["a6"] = self.a6
        d["a7"] = self.a7
        d["b1"] = self.b1
        d["b2"] = self.b2
        return repr(d)


class NasaSpecies:
    """Used to represent a species in a NASA thermodynamic database."""
    name:str                               # The name of the species.
    comment:str                            # The comment and data source of the species.
    idcode:str                             # The identification code of the species.
    elements:list[tuple]                   # The element symbols and coefficients in the formula of the species.
    aggregatecode:int                      # The aggregate state code of the species (zero for gas, nonzero for condensed phases).
    aggregatestate:str                     # The aggregate state of the species (liquid, gas, or solid) determined from name and aggregate state code.
    speciestype:str                        # The type of the species (product, reactant, product-reactant).
    molarmass:float                        # The molar mass of the species (in g/mol).
    dHf:float                              # The heat of formation \eq{\Delta H_{f}^{\circ}} at 298.15 K (in J/mol).
    dH0:float                              # The value of \eq{\Delta H_{0}^{\circ}=H^{\circ}(298.15)-H^{\circ}(0)} (in J/mol).
    H0:float                               # The assigned enthalpy (in J/mol) of the species when there are no temperature intervals.
    T0:float                               # The temperature (in K) corresponding to the assigned enthalpy when there are no temperature intervals.
    polynomials:list[NasaPolynomial] = []  # The polynomials used to compute standard thermodynamic properties of the species at different temperature ranges.};


def formatNumber(num:float):
    """Return int if possible, to avoid decimals in the constructed database."""
    return int(num) if int(num) == num else num


def correctNasaSpeciesName(name:str):
    """Replace CL and AL by Cl and Al in the name of a species.
    Replace also (L) to (l) as liquid designator."""
    return name.replace("CL", "Cl").replace("AL", "Al").replace("(L)", "(l)")


def identifyCommonSpeciesName(names:list[str]):
    """Given names Na2SO4(V), Na2SO4(IV), Na2SO4(I), Na2SO4(L), return Na2SO4"""
    if names == []:
        return ""
    if len(names) == 1:
        return names[0]
    i = 0
    commonname = ""

    def _cleanCommonNameEnding(name:str):
        iopen = name.rfind('(')
        iclose = name.rfind(')')
        return name if iopen <= iclose else name[:iopen]

    n = min([len(name) for name in names])
    for i in range(n):
        commonchar = names[0][i]
        for name in names:
            if name[i] != commonchar:
                return _cleanCommonNameEnding(commonname)
        commonname += commonchar
    return _cleanCommonNameEnding(commonname)


def determineSpeciesNameSuffix(specieslist:list[NasaSpecies]) -> str:
    """Determine if suffix (cd), (s) or (l) must be returned for a species formed by combining one or more others."""
    if len(specieslist) == 1:
        return ""  # there is just one species, so use whatever suffix it already has!
    names = set([x.name for x in specieslist])
    if len(names) == 1:
        return ""  # all species have the same name, so use whatever suffix it already has!
    aggregatestates = set([x.aggregatestate for x in specieslist])
    if aggregatestates == set(["Liquid", "Solid"]):
        return "(cd)"  # the species have different aggregate state, so use cd for condensed
    if aggregatestates == set(["Liquid"]):
        return "(l)"  # the species have common liquid aggregate state, so return (l) to denote liquid
    if aggregatestates == set(["Solid"]):
        return "(s)"  # the species have common solid aggregate state, so return (s) to denote solid
    raise RuntimeError("Could not identify a suffix for name of the common species formed out of {names}")


def determineSpeciesAggregateState(specieslist:list[NasaSpecies]) -> str:
    """Determine the aggregate state for a species formed by combining one or more others."""
    if len(specieslist) == 1:
        return specieslist[0].aggregatestate  # there is just one species, so use whatever agregate state it already has!
    aggregatestates = set([x.aggregatestate for x in specieslist])
    if aggregatestates == set(["Liquid", "Solid"]):
        return "CondensedPhase"  # the species have different aggregate state, liquid/solid, so assign condensed phase to it
    if aggregatestates == set(["Liquid", "Gas"]):
        return "Fluid"  # the species have different aggregate state, liquid/gas, so assign fluid phase to it
    if aggregatestates == set(["Liquid"]):
        return "Liquid"  # the species have common liquid aggregate state
    if aggregatestates == set(["Solid"]):
        return "Solid"  # the species have common solid aggregate state
    if aggregatestates == set(["Gas"]):
        return "Gas"  # the species have common gas aggregate state
    raise RuntimeError("Could not identify the aggregate state for the common species formed out of {names}")


def combineSpeciesBlocksWhenPossible(specieslist:list[NasaSpecies]) -> list[list[NasaSpecies]]:
    """Find all sequences of condensed species, terminating or not with a
    liquid species, that should form a single species. This usually happens
    when multiple species are present in the database to describe different
    physical state (solid/liquid) or different crystal structure. This method
    returns a list of list of species, where the inner list comprises one or
    more species that are in fact the same substance."""
    result = []
    i = 0
    n = len(specieslist)
    condensed = lambda s1, s2:         \
        s1.elements == s2.elements and \
        s1.dHf == s2.dHf           and \
        s1.polynomials != []       and \
        s2.polynomials != []

    while i < n:
        sublist = [specieslist[i]]
        j = i + 1
        while j < n and condensed(specieslist[i], specieslist[j]):
            sublist.append(specieslist[j])
            j += 1
        result.append(sublist)
        i += len(sublist)
    return result


def segment(lst:list, begin:int, length:int) -> list:
    """Return a segment slice of a list with given begin index and length"""
    return lst[begin:begin+length]


def parseFortranScientificNumber(word:str):
    """Parse a scientific number in FORTRAN format such as 1.23456D+03 and 1.23456D-11"""
    return float(word.replace("D", "e"))


def getStringBetweenColumns(line:str, begincol:int, endcol:int) -> str:
    """Return the string between given begin and end column indices using FORTRAN counting starting from 1."""
    assert begincol >= 1
    assert endcol >= begincol
    return line[begincol - 1:endcol]


def getIntBetweenColumns(line:str, begincol:int, endcol:int) -> int:
    """Return the integer between given begin and end column indices using FORTRAN counting starting from 1."""
    return int(getStringBetweenColumns(line, begincol, endcol))


def getFloatBetweenColumns(line:str, begincol:int, endcol:int) -> float:
    """Return the float between given begin and end column indices using FORTRAN counting starting from 1."""
    return parseFortranScientificNumber(getStringBetweenColumns(line, begincol, endcol))


def isCommentLine(line:str) -> bool:
    """Return True if this text line is a comment line."""
    line = line.strip()
    return len(line) > 0 and line[0] in ['!', '#']


def getNumberTextLinesForNextSpeciesBlock(lines:list[str]) -> int:
    """Return the number of string lines in the next species block.
    For example, if `lines` is:
    ~~~
    ALBr2             Gurvich,1996a pt1 p186 pt2 p149.
     2 tpis96 AL  1.00BR  2.00    0.00    0.00    0.00 0  186.7895380    -140662.125
        200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13397.875
     3.199375870D+04-7.119178970D+02 9.478258110D+00-4.875531670D-03 5.516512990D-06
    -3.340053040D-09 8.368476840D-13                -1.540591306D+04-1.742171366D+01
       1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13397.875
    -3.523782900D+05 4.671544170D+02 7.111908190D+00-5.551709200D-04 3.166301130D-07
    -5.521028330D-11 3.176725950D-15                -2.265004078D+04-2.695610360D+00
    CH4(L)            Methane. McBride,1996 pp85,93.
     0 g 6/96 C   1.00H   4.00    0.00    0.00    0.00 1   16.0424600     -89233.000
        111.643      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000
    ~~~
    then this method returns 8, i.e., the number of lines between `ALBr2` and `CH4(L)`.

    Args:
        lines (list): The text lines

    Returns:
        int: The number of lines in the next species block
    """
    recordline2 = lines[1]
    numintervals = getIntBetweenColumns(recordline2, 1, 2)
    return 3 if numintervals == 0 else 2 + 3 * numintervals


def parseFormula(formula:str) -> list[tuple]:
    assert len(formula) == 40
    pairs = []
    offset = 0
    for i in range(5):
        symbol = segment(formula, offset, 2).strip()
        coeff = float(segment(formula, offset + 2, 6))
        if coeff != 0.0:
            pairs.append((symbol, coeff))
        if coeff == 0.0 and symbol != "":
            raise RuntimeError(f"Cannot accept zero coefficient for element {symbol}.")
        offset += 8
    return pairs


def identifyAggregateState(name:str, aggregatecode:int) -> str:
    """Return the aggregate state of a species. This function returns Gas in
    case `aggregatecode` is zero. If `aggregatecode` is nonzero, then the
    aggregate state of the species is either liquid or solid. The
    identification is performed using information in `name`. If `(L)` or `(l)`
    is found in `name`, then Liquid is returned. If `name` ends with `)`, then
    Solid is returned. This allows any label to be used to denote solid or
    crystal configuration, with the most common ones being `(cr)`, `(a)`,
    `(b)`, `(c)`, `(I)`, `(II)`, `(III)`, and others. If none of these rules
    apply, then Liquid is returned, which covers species whose name does not
    contain `(L)` but the species is liquid, such as RP-1, JP-4, JP-5, IRFNA.

    Args:
        name (str): The name of the species such as Mg(cr), Mg(L), O2, NH4NO3(IV), NH4NO3(III), AL2O3(a)
        aggregatecode (int): The aggregate state code of the species (zero for gas, nonzero for condensed phases)

    Returns:
        str: The aggregate state of the species
    """
    if aggregatecode == 0:
        return "Gas"
    if name.find("(L)") != -1 or name.find("(l)") != -1:
        return "Liquid"
    if name[-1] == ')': # name ends with ) to indicate solid/crystal configuration such as (cr), (a), (b), (I), (II), (III)
        return "Solid"
    return "Liquid" # return liquid if nothing is provided (e.g., RP-1, JP-4, JP-5, IRFNA)


def createNasaPolynomial(lines: list[str]) -> NasaPolynomial:
    assert len(lines) == 3

    recordline3 = lines[0]
    recordline4 = lines[1]
    recordline5 = lines[2]

    segment = getStringBetweenColumns(recordline3,  1, 22).split()
    qvalues = getStringBetweenColumns(recordline3, 24, 63).split()

    if len(segment) != 2:
        raise RuntimeError(f"Expecting only two double values, "
        "Tmin and Tmax, within the first 22 chars of line:\n{recordline3}")

    if len(qvalues) != 8:
        raise RuntimeError(f"Expecting 8 values for q exponents "
        "(last being zero) between columns 24 and 63 of line:\n{recordline3}")

    poly = NasaPolynomial()
    poly.Tmin = float(segment[0])
    poly.Tmax = float(segment[1])
    poly.qN = getIntBetweenColumns(recordline3, 23, 23)
    poly.q1 = float(qvalues[0])
    poly.q2 = float(qvalues[1])
    poly.q3 = float(qvalues[2])
    poly.q4 = float(qvalues[3])
    poly.q5 = float(qvalues[4])
    poly.q6 = float(qvalues[5])
    poly.q7 = float(qvalues[6])
    poly.a1 = getFloatBetweenColumns(recordline4,  1, 16)
    poly.a2 = getFloatBetweenColumns(recordline4, 17, 32)
    poly.a3 = getFloatBetweenColumns(recordline4, 33, 48)
    poly.a4 = getFloatBetweenColumns(recordline4, 49, 64)
    poly.a5 = getFloatBetweenColumns(recordline4, 65, 80)
    poly.a6 = getFloatBetweenColumns(recordline5,  1, 16)
    poly.a7 = getFloatBetweenColumns(recordline5, 17, 32)
    poly.b1 = getFloatBetweenColumns(recordline5, 49, 64)
    poly.b2 = getFloatBetweenColumns(recordline5, 65, 80)

    if poly.qN != 7:
        raise RuntimeError(f"Cannot accept number of coefficients for "
        "Cp0 model different than 7. Got qN = {poly.qN}.")

    return poly


def createNasaSpecies(lines:list[str]) -> NasaSpecies:
    assert len(lines) >= 3
    assert len(lines) == getNumberTextLinesForNextSpeciesBlock(lines)

    recordline1 = lines[0]
    recordline2 = lines[1]
    recordline3 = lines[2]

    species = NasaSpecies()

    species.name           = getStringBetweenColumns(recordline1,  1, 18).strip()
    species.name           = correctNasaSpeciesName(species.name)
    species.comment        = getStringBetweenColumns(recordline1, 19, 80).strip()
    species.idcode         = getStringBetweenColumns(recordline2, 4, 9).strip()
    species.elements       = parseFormula(getStringBetweenColumns(recordline2, 11, 50))
    species.aggregatecode  = getIntBetweenColumns(recordline2, 52, 52)
    species.aggregatestate = identifyAggregateState(species.name, species.aggregatecode)
    species.speciestype    = ""
    species.molarmass      = getFloatBetweenColumns(recordline2, 53, 65)

    numintervals = getIntBetweenColumns(recordline2, 1, 2)

    species.polynomials = []
    if numintervals > 0:
        for i in range(numintervals):
            polynomial = createNasaPolynomial(segment(lines, 2 + 3*i, 3))
            species.polynomials.append(polynomial)

        species.dHf = getFloatBetweenColumns(recordline2, 66, 80)
        species.dH0 = getFloatBetweenColumns(recordline3, 66, 80)
        species.H0  = 0.0
        species.T0  = 0.0

    else:
        species.dHf = 0.0
        species.dH0 = 0.0
        species.H0  = getFloatBetweenColumns(recordline2, 66, 80)
        species.T0  = getFloatBetweenColumns(recordline3,  1, 11)

    return species


def createNasaSpeciesList(lines:list[str]) -> list[NasaSpecies]:
    specieslist:list[NasaSpecies] = []

    sublines = lines
    while sublines != []:

        # Change the type of all previous species with type ProductReactant to either Product or Reactant
        if sublines[0].strip() in ["END PRODUCTS", "END REACTANTS"]:
            newtype = "product" if sublines[0].strip() == "END PRODUCTS" else "reactant"
            for species in specieslist:
                if species.speciestype == "":
                    species.speciestype = newtype
            sublines = sublines[1:]
            continue

        # Construct a new species using the next lines of text
        numlines = getNumberTextLinesForNextSpeciesBlock(sublines)
        nextlines = segment(sublines, 0, numlines)
        new_species = createNasaSpecies(nextlines)
        specieslist.append(new_species)
        sublines = sublines[numlines:]

    return specieslist


def createFormulaString(elements:list[tuple[str, float]]) -> str:
    res:str = ""
    charge:int = 0
    for symbol, coeff in elements:
        if symbol == "E":
            charge = -int(coeff)
            continue
        coeff = "" if coeff == 1.0 else str(formatNumber(coeff))
        """"""
        res += symbol.title() + coeff
    if res == '' and charge == -1: return "e-"
    if charge ==  0: return res
    if charge ==  1: return f"{res}+"
    if charge == -1: return f"{res}-"
    if charge >=  2: return f"{res}+{abs(charge)}"
    if charge <= -2: return f"{res}-{abs(charge)}"


def createElementsString(elements:list[tuple[str, float]]) -> str:
    res = ""
    for symbol, coeff in elements:
        if symbol == "E":
            continue
        res += f"{formatNumber(coeff)}:{symbol.title()} "
        """"""
    return res.rstrip()


def getCharge(elements:list[tuple[str, float]]) -> int:
    for symbol, coeff in elements:
        if symbol == "E":
            return formatNumber(-coeff)
    return 0


def getTags(speciestype:str) -> list[str]:
    if speciestype == "":
        return "product reactant"
    return speciestype


def getStandardThermoModel(species_group:list[NasaSpecies]) -> dict:
    assert species_group != []

    species0 = species_group[0]

    data = {}
    if species0.polynomials == []:
        data["H0"] = species0.H0
        data["T0"] = species0.T0
        return {"Nasa": data}

    data["dHf"] = species0.dHf
    data["dH0"] = species0.dH0

    data["Polynomials"] = []
    for species in species_group:
        for polynomial in species.polynomials:
            child = {}
            child["State"] = species.aggregatestate
            child["Label"] = species.name
            child["Tmin"] = polynomial.Tmin
            child["Tmax"] = polynomial.Tmax
            child["a1"] = polynomial.a1
            child["a2"] = polynomial.a2
            child["a3"] = polynomial.a3
            child["a4"] = polynomial.a4
            child["a5"] = polynomial.a5
            child["a6"] = polynomial.a6
            child["a7"] = polynomial.a7
            child["b1"] = polynomial.b1
            child["b2"] = polynomial.b2
            # child["Comment"] = species.comment
            data["Polynomials"].append(child)

    return {"Nasa": data}


def serialize(species_groups:list[NasaSpecies]) -> dict:
    assert species_groups != []
    species0 = species_groups[0]
    names = [species.name for species in species_groups]
    name = identifyCommonSpeciesName(names)
    suffix = determineSpeciesNameSuffix(species_groups)
    name = name + suffix
    obj = dict()
    obj["Name"] = name
    obj["Formula"] = createFormulaString(species0.elements)
    obj["Elements"] = createElementsString(species0.elements)
    obj["Charge"] = getCharge(species0.elements)
    obj["Tags"] = getTags(species0.speciestype)
    obj["AggregateState"] = determineSpeciesAggregateState(species_groups)
    obj["StandardThermoModel"] = getStandardThermoModel(species_groups)
    if obj["Elements"] == "":  # this happens when species is e-
        obj.pop("Elements")
    if obj["Charge"] == 0.0:
        obj.pop("Charge")
    return obj


def createTextLines(file) -> list[str]:
    lines = []
    for line in file:
        line = line.rstrip()
        if line == "":
            continue
        if line == "thermo":
            file.readline() # skip the next line after thermo
            continue
        if isCommentLine(line):
            continue
        lines.append(line)
    return lines


def createDatabase(specieslist:list[NasaSpecies]) -> list:
    database = []
    species_groups = combineSpeciesBlocksWhenPossible(specieslist)
    for species_group in species_groups:
        database.append(serialize(species_group))
    return database

