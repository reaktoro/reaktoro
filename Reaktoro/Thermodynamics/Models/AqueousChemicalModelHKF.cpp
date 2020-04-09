// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "AqueousChemicalModelHKF.hpp"

// C++ includes
#include <map>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

/// The electrostatic constant \f$ \eta\f$ in the HKF model (in units of (A*cal)/mol)
const auto eta = 1.66027e+05;

/// Calculate the effective electrostatic radius of species (in units of A).
/// @param species The aqueous species instance of the ionic species
/// @return The effective electrostatic radius of the ionic species (in units of A)
auto effectiveIonicRadius(const Species& species) -> real;

/// Calculate the electrostatic Debye-Huckel parameter \f$ A_{\gamma}\f$.
/// @param T The temperature value in units of K)
/// @param P The pressure value (in units of Pa)
/// @return The electrostatic Debye-Huckel parameter \f$ A_{\gamma}\f$
auto debyeHuckelParamA(real T, real P) -> real;

/// Calculate the electrostatic Debye-Huckel parameter \f$ B_{\gamma}\f$.
/// @param T The temperature value in units of K)
/// @param P The pressure value (in units of Pa)
/// @return The electrostatic Debye-Huckel parameter \f$ B_{\gamma}\f$
auto debyeHuckelParamB(real T, real P) -> real;

/// Calculate the HKF solvation parameter \f$ b_{\mathrm{NaCl}}\f$.
/// @param T The temperature value in units of K)
/// @param P The pressure value (in units of Pa)
/// @return The HKF solvation parameter \f$ b_{\mathrm{NaCl}}\f$
auto solventParamNaCl(real T, real P) -> real;

/// Calculate the HKF solvation parameter \f$ b_{\mathrm{Na^{+}Cl^{-}}}\f$.
/// @param T The temperature value in units of K)
/// @param P The pressure value (in units of Pa)
/// @return The HKF solvation parameter \f$ b_{\mathrm{Na^{+}Cl^{-}}}\f$
auto shortRangeInteractionParamNaCl(real T, real P) -> real;

/// The interpolation data of the electrostatic Debye-Huckel parameter Agamma (in units of sqrt(kg/mol)).
const std::vector<double> Agamma_data =
{
    0.4939, 0.4871, 0.4810, 0.4755, 0.4705, 0.4617, 0.4544, 0.4427, 0.4337, 0.4266,
    0.5114, 0.5047, 0.4985, 0.4928, 0.4875, 0.4782, 0.4701, 0.4568, 0.4462, 0.4375,
    0.5354, 0.5281, 0.5213, 0.5151, 0.5094, 0.4991, 0.4901, 0.4750, 0.4628, 0.4526,
    0.5649, 0.5565, 0.5488, 0.5418, 0.5353, 0.5236, 0.5134, 0.4963, 0.4824, 0.4707,
    0.5996, 0.5897, 0.5807, 0.5725, 0.5649, 0.5515, 0.5398, 0.5202, 0.5043, 0.4910,
    0.6396, 0.6276, 0.6168, 0.6070, 0.5981, 0.5823, 0.5688, 0.5463, 0.5282, 0.5131,
    0.6855, 0.6707, 0.6573, 0.6455, 0.6348, 0.6161, 0.6003, 0.5743, 0.5536, 0.5365,
    0.7383, 0.7195, 0.7028, 0.6881, 0.6751, 0.6528, 0.6341, 0.6040, 0.5803, 0.5609,
    0.7995, 0.7753, 0.7538, 0.7355, 0.7194, 0.6925, 0.6703, 0.6352, 0.6081, 0.5861,
    0.8718, 0.8401, 0.8117, 0.7882, 0.7682, 0.7353, 0.7088, 0.6678, 0.6367, 0.6118,
    0.9596, 0.9169, 0.8783, 0.8476, 0.8221, 0.7815, 0.7498, 0.7018, 0.6662, 0.6381,
    1.0704, 1.0111, 0.9563, 0.9152, 0.8823, 0.8317, 0.7934, 0.7372, 0.6964, 0.6646,
    1.2183, 1.1325, 1.0500, 0.9934, 0.9502, 0.8865, 0.8400, 0.7739, 0.7272, 0.6915,
    1.4357, 1.3019, 1.1668, 1.0856, 1.0277, 0.9467, 0.8900, 0.8122, 0.7588, 0.7186,
    1.8233, 1.5767, 1.3188, 1.1970, 1.1175, 1.0132, 0.9438, 0.8521, 0.7910, 0.7460,
    2.2948, 2.2948, 1.5300, 1.3350, 1.2230, 1.0871, 1.0019, 0.8937, 0.8241, 0.7737,
    1.8509, 1.8509, 1.8509, 1.5114, 1.3490, 1.1700, 1.0648, 0.9373, 0.8581, 0.8019,
    2.3997, 2.3997, 2.3997, 1.7439, 1.5014, 1.2632, 1.1333, 0.9831, 0.8930, 0.8305,
    3.3844, 3.3844, 3.3844, 2.0579, 1.6877, 1.3685, 1.2078, 1.0311, 0.9291, 0.8597,
    4.7805, 4.7805, 4.7805, 2.4810, 1.9161, 1.4874, 1.2889, 1.0815, 0.9662, 0.8895,
    6.0949, 6.0949, 6.0949, 3.0235, 2.1930, 1.6213, 1.3770, 1.1344, 1.0047, 0.9202
};

/// The interpolation data of the electrostatic Debye-Huckel parameter Bgamma (in units of 10^-8 sqrt(kg/mol)/cm).
const std::vector<double> Bgamma_data =
{
    0.3253, 0.3251, 0.3250, 0.3249, 0.3248, 0.3248, 0.3248, 0.3250, 0.3253, 0.3257,
    0.3288, 0.3285, 0.3283, 0.3281, 0.3280, 0.3277, 0.3275, 0.3273, 0.3272, 0.3273,
    0.3328, 0.3325, 0.3322, 0.3319, 0.3317, 0.3312, 0.3309, 0.3304, 0.3300, 0.3298,
    0.3373, 0.3369, 0.3365, 0.3361, 0.3358, 0.3352, 0.3347, 0.3339, 0.3333, 0.3329,
    0.3421, 0.3416, 0.3411, 0.3406, 0.3402, 0.3394, 0.3388, 0.3377, 0.3369, 0.3363,
    0.3472, 0.3465, 0.3459, 0.3453, 0.3448, 0.3439, 0.3431, 0.3417, 0.3407, 0.3398,
    0.3526, 0.3517, 0.3509, 0.3502, 0.3495, 0.3484, 0.3474, 0.3458, 0.3445, 0.3435,
    0.3581, 0.3570, 0.3560, 0.3551, 0.3543, 0.3529, 0.3518, 0.3499, 0.3483, 0.3471,
    0.3639, 0.3625, 0.3612, 0.3601, 0.3591, 0.3575, 0.3561, 0.3538, 0.3521, 0.3506,
    0.3700, 0.3683, 0.3666, 0.3652, 0.3640, 0.3620, 0.3604, 0.3577, 0.3557, 0.3540,
    0.3767, 0.3744, 0.3722, 0.3705, 0.3690, 0.3666, 0.3646, 0.3615, 0.3592, 0.3573,
    0.3840, 0.3810, 0.3781, 0.3759, 0.3740, 0.3711, 0.3688, 0.3652, 0.3625, 0.3604,
    0.3925, 0.3885, 0.3845, 0.3815, 0.3792, 0.3756, 0.3729, 0.3688, 0.3657, 0.3633,
    0.4031, 0.3976, 0.3915, 0.3875, 0.3845, 0.3801, 0.3769, 0.3722, 0.3688, 0.3661,
    0.4180, 0.4096, 0.3994, 0.3940, 0.3901, 0.3848, 0.3809, 0.3755, 0.3717, 0.3687,
    0.4324, 0.4324, 0.4089, 0.4011, 0.3961, 0.3894, 0.3849, 0.3787, 0.3745, 0.3712,
    0.4209, 0.4209, 0.4209, 0.4090, 0.4024, 0.3942, 0.3889, 0.3819, 0.3771, 0.3736,
    0.4366, 0.4366, 0.4366, 0.4180, 0.4092, 0.3991, 0.3929, 0.3849, 0.3797, 0.3759,
    0.4560, 0.4560, 0.4560, 0.4282, 0.4165, 0.4042, 0.3969, 0.3879, 0.3822, 0.3780,
    0.4718, 0.4718, 0.4718, 0.4393, 0.4243, 0.4093, 0.4009, 0.3909, 0.3846, 0.3801,
    0.4780, 0.4780, 0.4780, 0.4503, 0.4322, 0.4145, 0.4049, 0.3937, 0.3869, 0.3822
};

/// The interpolation data of the HKF solvation parameter bNaCl (in units of kg/cal).
/// The values given on page 1477 (PDF page 229) of Helgeson et al. (1981), Table 29,
/// have different units than those below. The value in Table 29 are normalized by 2.303RT,
/// where R is in units of cal/mol*K).
/// Note: The values of bNaCl have been multiplied by 1.0e+07.
const std::vector<double> bNaCl_data =
{
     21.962,  22.211,  22.437,  22.643,  22.831,  23.162,  23.444,  23.901,  24.258,  24.548,
     18.081,  18.321,  18.542,  18.746,  18.934,  19.273,  19.569,  20.063,  20.461,  20.792,
     14.530,  14.783,  15.016,  15.232,  15.432,  15.794,  16.112,  16.648,  17.088,  17.458,
     11.235,  11.516,  11.775,  12.012,  12.233,  12.630,  12.979,  13.570,  14.055,  14.465,
     8.1250,  8.4490,  8.7450,  9.0150,  9.2640,  9.7100,  10.100,  10.757,  11.295,  11.749,
     5.1380,  5.5210,  5.8680,  6.1830,  6.4700,  6.9800,  7.4210,  8.1580,  8.7570,  9.2600,
     2.2140,  2.6780,  3.0950,  3.4670,  3.8040,  4.3950,  4.9000,  5.7330,  6.4020,  6.9600,
    -0.7100, -0.1370,  0.3760,  0.8260,  1.2270,  1.9190,  2.5030,  3.4490,  4.2000,  4.8210,
    -3.7030, -2.9810, -2.3360, -1.7830, -1.2980, -0.4770,  0.2010,  1.2830,  2.1270,  2.8170,
    -6.8580, -5.9300, -5.0960, -4.4020, -3.8070, -2.8240, -2.0280, -0.7870,  0.1640,  0.9320,
    -10.304, -9.0820, -7.9690, -7.0800, -6.3390, -5.1470, -4.2090, -2.7790, -1.7070, -0.8520,
    -14.247, -12.590, -11.042, -9.8730, -8.9340, -7.4760, -6.3620, -4.7090, -3.4980, -2.5470,
    -19.060, -16.716, -14.437, -12.856, -11.643, -9.8370, -8.5060, -6.5920, -5.2230, -4.1660,
    -25.556, -21.993, -18.343, -16.122, -14.524, -12.262, -10.663, -8.4390, -6.8930, -5.7180,
    -36.227, -29.865, -23.076, -19.806, -17.650, -14.783, -12.853, -10.265, -8.5190, -7.2140,
    -48.639, -48.639, -29.220, -24.100, -21.113, -17.436, -15.095, -12.080, -10.108, -8.6630,
    -38.002, -38.002, -38.002, -29.287, -25.032, -20.261, -17.410, -13.895, -11.671, -10.072,
    -52.214, -52.214, -52.214, -35.791, -29.554, -23.299, -19.818, -15.719, -13.214, -11.449,
    -76.652, -76.652, -76.652, -44.209, -34.857, -26.594, -22.338, -17.561, -14.746, -12.801,
    -110.79, -110.79, -110.79, -55.166, -41.133, -30.187, -24.986, -19.429, -16.271, -14.134,
    -143.76, -143.76, -143.76, -68.878, -48.527, -34.113, -27.778, -21.329, -17.796, -15.455
};

/// The interpolation data of the HKF short-range interaction parameter bNa+Cl- (in units of kg/mol).
/// Note: The values of bNapClm have been multiplied by 1.0e+02.
const std::vector<double> bNapClm_data =
{
    -15.448, -14.872, -14.390, -14.002, -13.708, -13.401, -13.471, -14.739, -17.512, -21.789,
    -9.7520, -9.5630, -9.4040, -9.2760, -9.1780, -9.0730, -9.0900, -9.4870, -10.370, -11.739,
    -5.6300, -5.6030, -5.5790, -5.5600, -5.5440, -5.5240, -5.5180, -5.5520, -5.6470, -5.8010,
    -2.4110, -2.4660, -2.5100, -2.5460, -2.5710, -2.5940, -2.5770, -2.4300, -2.1280, -1.6720,
     0.2440,  0.1450,  0.0630, -0.0020, -0.0510, -0.0970, -0.0750,  0.1720,  0.6910,  1.4820,
     2.5290,  2.4050,  2.3010,  2.2180,  2.1560,  2.0970,  2.1220,  2.4260,  3.0690,  4.0520,
     4.5590,  4.4210,  4.3050,  4.2120,  4.1420,  4.0740,  4.1010,  4.4380,  5.1530,  6.2460,
     6.4060,  6.2630,  6.1390,  6.0400,  5.9660,  5.8940,  5.9210,  6.2760,  7.0320,  8.1870,
     8.1190,  7.9760,  7.8480,  7.7460,  7.6700,  7.5950,  7.6230,  7.9870,  8.7630,  9.9500,
     9.7290,  9.5910,  9.4620,  9.3590,  9.2820,  9.2050,  9.2330,  9.6000,  10.384,  11.583,
     11.259,  11.130,  11.001,  10.897,  10.820,  10.743,  10.771,  11.138,  11.921,  13.119,
     12.723,  12.608,  12.479,  12.377,  12.300,  12.224,  12.251,  12.614,  13.391,  14.581,
     14.133,  14.036,  13.909,  13.807,  13.731,  13.656,  13.682,  14.041,  14.808,  15.984,
     15.496,  15.421,  15.296,  15.196,  15.121,  15.047,  15.073,  15.426,  16.182,  17.339,
     16.818,  16.771,  16.648,  16.550,  16.476,  16.403,  16.428,  16.775,  17.517,  18.654,
     18.090,  18.090,  17.969,  17.872,  17.800,  17.728,  17.753,  18.093,  18.821,  19.935,
     19.380,  19.380,  19.262,  19.167,  19.096,  19.026,  19.050,  19.383,  20.096,  21.188,
     20.645,  20.645,  20.529,  20.437,  20.367,  20.298,  20.322,  20.648,  21.346,  22.415,
     21.887,  21.887,  21.774,  21.683,  21.615,  21.547,  21.570,  21.890,  22.572,  23.618,
     23.108,  23.108,  22.997,  22.908,  22.841,  22.775,  22.798,  23.110,  23.777,  24.801,
     24.308,  24.308,  24.199,  24.113,  24.048,  23.983,  24.005,  24.310,  24.963,  25.964
};

/// The pressure coordinates of the interpolation points (in units of celsius).
const std::vector<double> temperature_range =
{
      0,  25,  50,  75,
    100, 125, 150, 175,
    200, 225, 250, 275,
    300, 325, 350, 375,
    400, 425, 450, 475, 500
};

/// The pressure coordinates of the interpolation points (in units of bar).
const std::vector<double> pressure_range =
{
       1,  250,  500,  750, 1000,
    1500, 2000, 3000, 4000, 5000
};

/// The effective electrostatic radii of ionic species (in units of angstrom).
/// This data was taken from Table 3 of Helgeson et al. (1981).
const std::map<std::string, real> effective_radii =
{
    {"H+"  , 3.08}, {"Fe+++", 3.46},
    {"Li+" , 1.64}, {"Al+++", 3.33},
    {"Na+" , 1.91}, {"Au+++", 3.72},
    {"K+"  , 2.27}, {"La+++", 3.96},
    {"Rb+" , 2.41}, {"Gd+++", 3.79},
    {"Cs+" , 2.61}, {"In+++", 3.63},
    {"NH4+", 2.31}, {"Ca+++", 3.44},
    {"Ag+" , 2.20}, {"F-"   , 1.33},
    {"Au+" , 2.31}, {"Cl-"  , 1.81},
    {"Cu+" , 1.90}, {"Br-"  , 1.96},
    {"Mg++", 2.54}, {"I-"   , 2.20},
    {"Sr++", 3.00}, {"OH-"  , 1.40},
    {"Ca++", 2.87}, {"HS-"  , 1.84},
    {"Ba++", 3.22}, {"NO3-" , 2.81},
    {"Pb++", 3.08}, {"HCO3-", 2.10},
    {"Zn++", 2.62}, {"HSO4-", 2.37},
    {"Cu++", 2.60}, {"ClO4-", 3.59},
    {"Cd++", 2.85}, {"ReO4-", 4.23},
    {"Hg++", 2.98}, {"SO4--", 3.15},
    {"Fe++", 2.62}, {"CO3--", 2.81},
    {"Mn++", 2.68}
};

/// The bilinear interpolator of the HKF parameter Agamma
BilinearInterpolator Agamma(pressure_range, temperature_range, Agamma_data);

/// The bilinear interpolator of the HKF parameter Bgamma
BilinearInterpolator Bgamma(pressure_range, temperature_range, Bgamma_data);

/// The bilinear interpolator of the HKF parameter bNaCl
BilinearInterpolator bNaCl(pressure_range, temperature_range, bNaCl_data);

/// The bilinear interpolator of the HKF parameter bNapClm
BilinearInterpolator bNapClm(pressure_range, temperature_range, bNapClm_data);

auto debyeHuckelParamA(real T, real P) -> real
{
    const auto TdegC = convertKelvinToCelsius(T);
    const auto Pbar  = convertPascalToBar(P);
    return Agamma(Pbar, TdegC);
}

auto debyeHuckelParamB(real T, real P) -> real
{
    const auto TdegC = convertKelvinToCelsius(T);
    const auto Pbar  = convertPascalToBar(P);
    return Bgamma(Pbar, TdegC);
}

auto solventParamNaCl(real T, real P) -> real
{
    const auto TdegC = convertKelvinToCelsius(T);
    const auto Pbar  = convertPascalToBar(P);
    return 1.0e-07 * bNaCl(Pbar, TdegC);
}

auto shortRangeInteractionParamNaCl(real T, real P) -> real
{
    const auto TdegC = convertKelvinToCelsius(T);
    const auto Pbar  = convertPascalToBar(P);
    return 1.0e-02 * bNapClm(Pbar, TdegC);
}

auto effectiveIonicRadius(const Species& species) -> real
{
    // Find the effective ionic radius of the species in `effective_radii`.
    // Note that `species` might have a different name convention than those
    // used in `effective_radii`. Thus, we need to check if the name of given
    // `species` is an alternative to a name in `effective_radii`.
    for(auto pair : effective_radii)
        if(isAlternativeChargedSpeciesName(species.name(), pair.first))
            return pair.second;

    // The electrical charge of the species
    const auto Zi = species.charge();

    // Estimated effective ionci radius of the species based on TOUGHREACT approach
    if(Zi == -1) return 1.81;        // based on Cl- value
    if(Zi == -2) return 3.00;        // based on rounded average of CO3-- and SO4-- values
    if(Zi == -3) return 4.20;        // based on estimation from straight line fit with charge
    if(Zi == +1) return 2.31;        // based on NH4+ value
    if(Zi == +2) return 2.80;        // based on rounded average for +2 species in the HKF table of effective ionic radii
    if(Zi == +3) return 3.60;        // based on rounded average for +3 species in the HKF table of effective ionic radii
    if(Zi == +4) return 4.50;        // based on estimaton using HKF eq. 142
    if(Zi <  -3) return -Zi*4.2/3.0; // based on linear extrapolation
    return Zi*4.5/4.0;               // based on linear extrapolation
}

} // namespace

auto aqueousChemicalModelHKF(const AqueousMixture& mixture)-> ActivityModelFn
{
    // The number of species in the mixture
    const auto num_species = mixture.numSpecies();

    // The number of charged species in the mixture
    const auto num_charged_species = mixture.numChargedSpecies();

    // The indices of the charged species
    const auto icharged_species = mixture.indicesChargedSpecies();

    // The index of the water species
    const auto iwater = mixture.indexWater();

    // The effective electrostatic radii of the charged species
    std::vector<real> effective_radii;

    // The electrical charges of the charged species only
    std::vector<double> charges;

    // The Born coefficient of the ion H+
    const auto omegaH = 0.5387e+05;

    // The natural log of 10
    const auto ln10 = std::log(10);

    // The molar mass of water
    const auto Mw = waterMolarMass;

    // The state of the aqueous mixture
    AqueousMixtureState state;

    // Collect the effective radii of the ions
    for(Index idx_ion : icharged_species)
    {
        const Species& species = mixture.species(idx_ion);
        effective_radii.push_back(effectiveIonicRadius(species));
        charges.push_back(species.charge());
    }

    // Define the chemical model function of the aqueous phase
    ActivityModelFn model = [=](ActivityProps res, real T, real P, ArrayXrConstRef x) mutable
    {
        using std::abs;
        using std::log;
        using std::log10;
        using std::pow;
        using std::sqrt;

        // Evaluate the state of the aqueous mixture
        state = mixture.state(T, P, x);

        // Auxiliary references to state variables
        const auto& I = state.Ie;
        const auto& m = state.m;

        // The square root of the ionic strength
        const auto sqrtI = sqrt(I);

        // The mole fraction of the water species and its molar derivatives
        const auto xw = x[iwater];

        // The ln and log10 of water mole fraction
        const auto ln_xw = log(xw);
        const auto log10_xw = log10(xw);

        // The alpha parameter
        const auto alpha = xw/(1.0 - xw) * log10_xw;

        // The parameters for the HKF model
        const auto A = debyeHuckelParamA(T, P);
        const auto B = debyeHuckelParamB(T, P);
        const auto bNaCl = solventParamNaCl(T, P);
        const auto bNapClm = shortRangeInteractionParamNaCl(T, P);

        // The osmotic coefficient of the aqueous phase
        real phi = {};

        // Set the activity coefficients of the neutral species to
        // water mole fraction to convert it to molality scale
        res.ln_g.fill(0.0);
//        res.ln_activity_coefficients = ln_xw;

        // Loop over all charged species in the mixture
        for(auto i = 0; i < num_charged_species; ++i)
        {
            // The index of the charged species in the mixture
            const auto ispecies = icharged_species[i];

            // The molality of the charged species and its molar derivatives
            const auto mi = m[ispecies];

            // Check if the molality of the charged species is zero
            if(mi == 0.0)
                continue;

            // The electrical charge of the charged species
            const auto z = charges[i];
            const auto z2 = z*z;

            // The effective radius of the charged species
            const auto eff_radius = effective_radii[i];

            // The Born coefficient of the current charged species
            const auto omega = eta*z2/eff_radius - z*omegaH;

            // The absolute Born coefficient of the current charged species
            const auto omega_abs = eta*z2/eff_radius;

            // The Debye-Huckel ion size parameter of the current ion as computed by Reed (1982) and also in TOUGHREACT
            const auto a = (z < 0) ?
                2.0*(eff_radius + 1.91*abs(z))/(abs(z) + 1.0) :
                2.0*(eff_radius + 1.81*abs(z))/(abs(z) + 1.0);

            // The \Lamba parameter of the HKF activity coefficient model and its molar derivatives
            const auto lambda = 1.0 + a*B*sqrtI;

            // The log10 of the activity coefficient of the charged species (in mole fraction scale) and its molar derivatives
            // This is the equation (298) in Helgeson et a. (1981) paper, page 230.
            const auto log10_gi = -(A*z2*sqrtI)/lambda + log10_xw + (omega_abs * bNaCl + bNapClm - 0.19*(abs(z) - 1.0)) * I;

            // Set the activity coefficient of the current charged species
            res.ln_g[ispecies] = log10_gi * ln10;

            // Check if the mole fraction of water is one
            if(xw != 1.0)
            {
                // The sigma parameter of the current ion and its molar derivatives
                const auto sigma = 3.0/pow(a*B*sqrtI, 3) * (lambda - 1.0/lambda - 2.0*log(lambda));

                // The psi contribution of the current ion and its molar derivatives
                const auto psi = A*z2*sqrtI*sigma/3.0 + alpha - 0.5*(omega*bNaCl + bNapClm - 0.19*(abs(z) - 1.0)) * I;

                // Update the osmotic coefficient with the contribution of the current charged species
                phi += mi * psi;
            }
        }

        // Set the activities of the solutes (molality scale)
        res.ln_a = res.ln_g + m.log();

        // Set the activity of water (in mole fraction scale)
        if(xw != 1.0) res.ln_a[iwater] = ln10 * Mw * phi;
                 else res.ln_a[iwater] = ln_xw;

        // Set the activity coefficient of water (mole fraction scale)
        res.ln_g[iwater] = res.ln_a[iwater] - ln_xw;
    };

    return model;
}

} // namespace Reaktoro
