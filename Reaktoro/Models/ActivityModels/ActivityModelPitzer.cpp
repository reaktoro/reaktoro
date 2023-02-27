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

#include "ActivityModelPitzer.hpp"

// C++ includes
#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Serialization/Models/ActivityModels.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

using std::log;
using std::sqrt;
using std::exp;
using std::abs;
using std::round;

/// The calculated values of function J0(x) on the points x[i] = 0.0 + i*0.01 in the region (0.0, 1.0)
const Vec<double> J0region1 =
{
    0.00000000e+00, 7.06207558e-05, 2.38544646e-04, 4.80518193e-04, 7.84985974e-04, 1.14430570e-03, 1.55285818e-03, 2.00625822e-03, 2.50094831e-03, 3.03396197e-03,
    3.60277406e-03, 4.20520037e-03, 4.83932706e-03, 5.50345931e-03, 6.19608287e-03, 6.91583456e-03, 7.66147916e-03, 8.43189104e-03, 9.22603924e-03, 1.00429754e-02,
    1.08818235e-02, 1.17417716e-02, 1.26220648e-02, 1.35219989e-02, 1.44409153e-02, 1.53781966e-02, 1.63332630e-02, 1.73055680e-02, 1.82945965e-02, 1.92998614e-02,
    2.03209017e-02, 2.13572800e-02, 2.24085812e-02, 2.34744104e-02, 2.45543916e-02, 2.56481662e-02, 2.67553921e-02, 2.78757422e-02, 2.90089038e-02, 3.01545772e-02,
    3.13124753e-02, 3.24823227e-02, 3.36638550e-02, 3.48568178e-02, 3.60609669e-02, 3.72760670e-02, 3.85018914e-02, 3.97382219e-02, 4.09848476e-02, 4.22415655e-02,
    4.35081791e-02, 4.47844989e-02, 4.60703414e-02, 4.73655292e-02, 4.86698908e-02, 4.99832598e-02, 5.13054751e-02, 5.26363808e-02, 5.39758253e-02, 5.53236619e-02,
    5.66797479e-02, 5.80439451e-02, 5.94161188e-02, 6.07961386e-02, 6.21838773e-02, 6.35792115e-02, 6.49820210e-02, 6.63921888e-02, 6.78096011e-02, 6.92341471e-02,
    7.06657187e-02, 7.21042108e-02, 7.35495206e-02, 7.50015482e-02, 7.64601961e-02, 7.79253690e-02, 7.93969741e-02, 8.08749206e-02, 8.23591201e-02, 8.38494861e-02,
    8.53459342e-02, 8.68483818e-02, 8.83567482e-02, 8.98709545e-02, 9.13909238e-02, 9.29165805e-02, 9.44478508e-02, 9.59846627e-02, 9.75269454e-02, 9.90746298e-02,
    1.00627648e-01, 1.02185934e-01, 1.03749423e-01, 1.05318051e-01, 1.06891756e-01, 1.08470476e-01, 1.10054152e-01, 1.11642726e-01, 1.13236139e-01, 1.14834335e-01, 1.16437258e-01
};

/// The calculated values of function J0(x) on the points x[i] = 1.0 + i*0.1 in the region (1.0, 10.0)
const Vec<double> J0region2 =
{
    1.16437258e-01, 1.32715063e-01, 1.49411766e-01, 1.66486062e-01, 1.83902904e-01, 2.01632240e-01, 2.19648056e-01, 2.37927656e-01, 2.56451090e-01, 2.75200712e-01,
    2.94160824e-01, 3.13317387e-01, 3.32657787e-01, 3.52170640e-01, 3.71845632e-01, 3.91673386e-01, 4.11645345e-01, 4.31753678e-01, 4.51991199e-01, 4.72351292e-01,
    4.92827852e-01, 5.13415235e-01, 5.34108209e-01, 5.54901912e-01, 5.75791820e-01, 5.96773717e-01, 6.17843660e-01, 6.38997965e-01, 6.60233175e-01, 6.81546049e-01,
    7.02933538e-01, 7.24392774e-01, 7.45921053e-01, 7.67515825e-01, 7.89174680e-01, 8.10895338e-01, 8.32675641e-01, 8.54513543e-01, 8.76407104e-01, 8.98354481e-01,
    9.20353920e-01, 9.42403756e-01, 9.64502399e-01, 9.86648338e-01, 1.00884013e+00, 1.03107639e+00, 1.05335582e+00, 1.07567714e+00, 1.09803916e+00, 1.12044072e+00,
    1.14288073e+00, 1.16535811e+00, 1.18787186e+00, 1.21042101e+00, 1.23300461e+00, 1.25562177e+00, 1.27827163e+00, 1.30095335e+00, 1.32366614e+00, 1.34640922e+00,
    1.36918186e+00, 1.39198333e+00, 1.41481295e+00, 1.43767006e+00, 1.46055401e+00, 1.48346418e+00, 1.50639998e+00, 1.52936082e+00, 1.55234615e+00, 1.57535543e+00,
    1.59838813e+00, 1.62144375e+00, 1.64452179e+00, 1.66762179e+00, 1.69074329e+00, 1.71388583e+00, 1.73704899e+00, 1.76023235e+00, 1.78343550e+00, 1.80665805e+00,
    1.82989961e+00, 1.85315981e+00, 1.87643830e+00, 1.89973472e+00, 1.92304873e+00, 1.94638000e+00, 1.96972820e+00, 1.99309303e+00, 2.01647419e+00, 2.03987136e+00, 2.06328427e+00
};

/// The calculated values of function J0(x) on the points x[i] = 10.0 + i*1.0 in the region (10.0, 100.0)
const Vec<double> J0region3 =
{
    2.06328427e+00, 2.29822044e+00, 2.53446209e+00, 2.77181680e+00, 3.01013128e+00, 3.24928140e+00, 3.48916529e+00, 3.72969829e+00, 3.97080934e+00, 4.21243822e+00,
    4.45453343e+00, 4.69705062e+00, 4.93995136e+00, 5.18320209e+00, 5.42677335e+00, 5.67063917e+00, 5.91477653e+00, 6.15916491e+00, 6.40378598e+00, 6.64862328e+00,
    6.89366200e+00, 7.13888874e+00, 7.38429135e+00, 7.62985879e+00, 7.87558100e+00, 8.12144876e+00, 8.36745365e+00, 8.61358792e+00, 8.85984443e+00, 9.10621661e+00,
    9.35269838e+00, 9.59928411e+00, 9.84596857e+00, 1.00927469e+01, 1.03396146e+01, 1.05865675e+01, 1.08336015e+01, 1.10807131e+01, 1.13278988e+01, 1.15751554e+01,
    1.18224798e+01, 1.20698692e+01, 1.23173210e+01, 1.25648325e+01, 1.28124014e+01, 1.30600255e+01, 1.33077027e+01, 1.35554308e+01, 1.38032082e+01, 1.40510328e+01,
    1.42989031e+01, 1.45468174e+01, 1.47947742e+01, 1.50427719e+01, 1.52908093e+01, 1.55388850e+01, 1.57869977e+01, 1.60351462e+01, 1.62833294e+01, 1.65315461e+01,
    1.67797954e+01, 1.70280762e+01, 1.72763875e+01, 1.75247285e+01, 1.77730983e+01, 1.80214960e+01, 1.82699208e+01, 1.85183719e+01, 1.87668486e+01, 1.90153502e+01,
    1.92638760e+01, 1.95124254e+01, 1.97609976e+01, 2.00095922e+01, 2.02582084e+01, 2.05068458e+01, 2.07555038e+01, 2.10041819e+01, 2.12528795e+01, 2.15015962e+01,
    2.17503315e+01, 2.19990850e+01, 2.22478561e+01, 2.24966446e+01, 2.27454499e+01, 2.29942718e+01, 2.32431097e+01, 2.34919634e+01, 2.37408324e+01, 2.39897164e+01, 2.42386152e+01
};

/// The calculated values of function J0(x) on the points x[i] = 100.0 + i*100.0 in the region (100.0, 1000.0)
const Vec<double> J0region4 =
{
    2.42386152e+01, 4.91709893e+01, 7.41388190e+01, 9.91190695e+01, 1.24105388e+02, 1.49095205e+02, 1.74087256e+02, 1.99080833e+02, 2.24075510e+02, 2.49071007e+02
};

/// The calculated values of function J0(x) on the points x[i] = 1000.0 + i*1000.0 in the region (1000.0, 10000.0)
const Vec<double> J0region5 =
{
    2.49071007e+02, 4.99046818e+02, 7.49036364e+02, 9.99030282e+02, 1.24902622e+03, 1.49902328e+03, 1.74902104e+03, 1.99901926e+03, 2.24901780e+03, 2.49901659e+03
};

/// The calculated values of function J1(x) on the points x[i] = 0.0 + i*0.01 in the region (0.0, 1.0)
const Vec<double> J1region1 =
{
    0.00000000e-00, 1.25151744e-04, 4.14749318e-04, 8.24176731e-04, 1.33185071e-03, 1.92376087e-03, 2.58976962e-03, 3.32207863e-03, 4.11444419e-03, 4.96172312e-03,
    5.85958733e-03, 6.80433336e-03, 7.79274934e-03, 8.82201875e-03, 9.88964869e-03, 1.09934151e-02, 1.21313201e-02, 1.33015585e-02, 1.45024900e-02, 1.57326178e-02,
    1.69905698e-02, 1.82750838e-02, 1.95849943e-02, 2.09192222e-02, 2.22767652e-02, 2.36566899e-02, 2.50581251e-02, 2.64802559e-02, 2.79223180e-02, 2.93835939e-02,
    3.08634081e-02, 3.23611240e-02, 3.38761406e-02, 3.54078895e-02, 3.69558327e-02, 3.85194600e-02, 4.00982869e-02, 4.16918532e-02, 4.32997208e-02, 4.49214724e-02,
    4.65567101e-02, 4.82050541e-02, 4.98661416e-02, 5.15396256e-02, 5.32251741e-02, 5.49224690e-02, 5.66312055e-02, 5.83510910e-02, 6.00818448e-02, 6.18231971e-02,
    6.35748887e-02, 6.53366701e-02, 6.71083013e-02, 6.88895510e-02, 7.06801965e-02, 7.24800229e-02, 7.42888229e-02, 7.61063965e-02, 7.79325506e-02, 7.97670984e-02,
    8.16098597e-02, 8.34606598e-02, 8.53193299e-02, 8.71857068e-02, 8.90596322e-02, 9.09409528e-02, 9.28295200e-02, 9.47251900e-02, 9.66278229e-02, 9.85372835e-02,
    1.00453440e-01, 1.02376165e-01, 1.04305334e-01, 1.06240827e-01, 1.08182528e-01, 1.10130321e-01, 1.12084096e-01, 1.14043746e-01, 1.16009167e-01, 1.17980255e-01,
    1.19956912e-01, 1.21939042e-01, 1.23926549e-01, 1.25919343e-01, 1.27917334e-01, 1.29920434e-01, 1.31928560e-01, 1.33941627e-01, 1.35959556e-01, 1.37982267e-01,
    1.40009683e-01, 1.42041730e-01, 1.44078334e-01, 1.46119422e-01, 1.48164927e-01, 1.50214778e-01, 1.52268909e-01, 1.54327256e-01, 1.56389753e-01, 1.58456339e-01, 1.60526953e-01
};

/// The calculated values of function J1(x) on the points x[i] = 1.0 + i*0.1 in the region (1.0, 10.0)
const Vec<double> J1region2 =
{
    1.60526953e-01, 1.81442105e-01, 2.02701762e-01, 2.24262453e-01, 2.46088272e-01, 2.68149211e-01, 2.90419926e-01, 3.12878818e-01, 3.35507329e-01, 3.58289398e-01,
    3.81211036e-01, 4.04259985e-01, 4.27425447e-01, 4.50697864e-01, 4.74068736e-01, 4.97530475e-01, 5.21076279e-01, 5.44700030e-01, 5.68396207e-01, 5.92159812e-01,
    6.15986308e-01, 6.39871564e-01, 6.63811811e-01, 6.87803601e-01, 7.11843772e-01, 7.35929420e-01, 7.60057870e-01, 7.84226655e-01, 8.08433496e-01, 8.32676281e-01,
    8.56953056e-01, 8.81262001e-01, 9.05601427e-01, 9.29969760e-01, 9.54365532e-01, 9.78787370e-01, 1.00323399e+00, 1.02770420e+00, 1.05219687e+00, 1.07671094e+00,
    1.10124543e+00, 1.12579939e+00, 1.15037195e+00, 1.17496229e+00, 1.19956963e+00, 1.22419322e+00, 1.24883238e+00, 1.27348644e+00, 1.29815479e+00, 1.32283683e+00,
    1.34753201e+00, 1.37223979e+00, 1.39695969e+00, 1.42169120e+00, 1.44643390e+00, 1.47118733e+00, 1.49595110e+00, 1.52072481e+00, 1.54550809e+00, 1.57030058e+00,
    1.59510196e+00, 1.61991188e+00, 1.64473005e+00, 1.66955617e+00, 1.69438996e+00, 1.71923115e+00, 1.74407948e+00, 1.76893471e+00, 1.79379660e+00, 1.81866492e+00,
    1.84353946e+00, 1.86842001e+00, 1.89330637e+00, 1.91819835e+00, 1.94309576e+00, 1.96799843e+00, 1.99290619e+00, 2.01781888e+00, 2.04273633e+00, 2.06765841e+00,
    2.09258496e+00, 2.11751584e+00, 2.14245093e+00, 2.16739009e+00, 2.19233319e+00, 2.21728012e+00, 2.24223076e+00, 2.26718500e+00, 2.29214274e+00, 2.31710386e+00, 2.34206827e+00
};

/// The calculated values of function J1(x) on the points x[i] = 10.0 + i*1.0 in the region (10.0, 100.0)
const Vec<double> J1region3 =
{
    2.34206827e+00, 2.59187348e+00, 2.84191919e+00, 3.09214678e+00, 3.34251316e+00, 3.59298607e+00, 3.84354101e+00, 4.09415916e+00, 4.34482586e+00, 4.59552961e+00,
    4.84626133e+00, 5.09701375e+00, 5.34778108e+00, 5.59855865e+00, 5.84934266e+00, 6.10013007e+00, 6.35091839e+00, 6.60170559e+00, 6.85249003e+00, 7.10327038e+00,
    7.35404555e+00, 7.60481466e+00, 7.85557702e+00, 8.10633207e+00, 8.35707936e+00, 8.60781855e+00, 8.85854938e+00, 9.10927165e+00, 9.35998525e+00, 9.61069006e+00,
    9.86138606e+00, 1.01120732e+01, 1.03627516e+01, 1.06134211e+01, 1.08640820e+01, 1.11147342e+01, 1.13653778e+01, 1.16160129e+01, 1.18666397e+01, 1.21172583e+01,
    1.23678687e+01, 1.26184711e+01, 1.28690656e+01, 1.31196523e+01, 1.33702314e+01, 1.36208030e+01, 1.38713673e+01, 1.41219243e+01, 1.43724742e+01, 1.46230171e+01,
    1.48735532e+01, 1.51240825e+01, 1.53746053e+01, 1.56251215e+01, 1.58756314e+01, 1.61261350e+01, 1.63766325e+01, 1.66271241e+01, 1.68776097e+01, 1.71280895e+01,
    1.73785636e+01, 1.76290322e+01, 1.78794953e+01, 1.81299530e+01, 1.83804055e+01, 1.86308528e+01, 1.88812950e+01, 1.91317322e+01, 1.93821645e+01, 1.96325920e+01,
    1.98830149e+01, 2.01334330e+01, 2.03838467e+01, 2.06342558e+01, 2.08846606e+01, 2.11350611e+01, 2.13854573e+01, 2.16358494e+01, 2.18862373e+01, 2.21366213e+01,
    2.23870013e+01, 2.26373774e+01, 2.28877497e+01, 2.31381183e+01, 2.33884832e+01, 2.36388444e+01, 2.38892021e+01, 2.41395563e+01, 2.43899070e+01, 2.46402544e+01, 2.48905984e+01
};

/// The calculated values of function J1(x) on the points x[i] = 100.0 + i*100.0 in the region (100.0, 1000.0)
const Vec<double> J1region4 =
{
    2.48905984e+01, 4.99141212e+01, 7.49270387e+01, 9.99355618e+01, 1.24941747e+02, 1.49946507e+02, 1.74950320e+02, 1.99953466e+02, 2.24956118e+02, 2.49958395e+02
};

/// The calculated values of function J1(x) on the points x[i] = 1000.0 + i*1000.0 in the region (1000.0, 10000.0)
const Vec<double> J1region5 =
{
    2.49958395e+02, 4.99971193e+02, 7.49977057e+02, 9.99980578e+02, 1.24998298e+03, 1.49998474e+03, 1.74998611e+03, 1.99998720e+03, 2.24998810e+03, 2.49998886e+03
};

/// The temperature points (in K) for the interpolation of Aphi from 0 to 350 °C.
const Vec<double> Aphi_temperatures =
{
    273.15, 298.15, 323.15, 348.15, 373.15, 398.15, 423.15, 448.15, 473.15, 498.15, 523.15, 548.15, 573.15, 598.15, 623.15
};

/// The pressure points (in Pa) for the interpolation of Aphi from 1 to 5000 bar.
const Vec<double> Aphi_pressures =
{
    1e5, 100e5, 200e5, 400e5, 600e5, 800e5, 1000e5, 1500e5, 2000e5, 3000e5, 4000e5, 5000e5,
};

/// The data used for the interpolation of Aphi.
/// The following equation was used to compute these data:
/// \f[
///      A^{\phi}=1.400684\cdot10^{6}\frac{1}{\rho_{w}}\left(\frac{\rho_{w}}{\epsilon_{w}T}\right)^{\frac{3}{2}},
/// \f]
/// where @f$T@f$ is temperature in K, @f$\epsilon_{w}@f$ is the dielectric constant of water, and @f$\rho_{w}@f$ its density in g/m<sup>3<\sup>
const Vec<double> Aphi_data =
{
    0.37674222, 0.39147278, 0.41031975, 0.43330198, 0.46058910, 0.49248630, 0.52954906, 0.57258689, 0.62280952, 0.68208055, 0.75343848, 0.84226817, 0.95932062, 1.13030749, 1.43881804,
    0.37511283, 0.38962104, 0.40816127, 0.43071699, 0.45741532, 0.48854946, 0.52461006, 0.56635820, 0.61497721, 0.67238823, 0.74195091, 0.83021165, 0.95211387, 1.13030749, 1.43881804,
    0.37349119, 0.38778969, 0.40603633, 0.42818435, 0.45432255, 0.48468653, 0.51967722, 0.55991097, 0.60632611, 0.66040079, 0.72461917, 0.80356970, 0.90695201, 1.06002794, 1.35961218,
    0.37031872, 0.38423960, 0.40194465, 0.42334123, 0.44845489, 0.47742645, 0.51051253, 0.54810417, 0.59077384, 0.63937254, 0.69522720, 0.76056419, 0.83951064, 0.94084255, 1.08752618,
    0.36723706, 0.38083131, 0.39805034, 0.41877218, 0.44297410, 0.47072404, 0.50217061, 0.53754295, 0.57716443, 0.62148709, 0.67116338, 0.72719329, 0.79124059, 0.86637464, 0.95905835,
    0.36424217, 0.37755559, 0.39433845, 0.41445319, 0.43784086, 0.46451346, 0.49453834, 0.52802729, 0.56513171, 0.60604478, 0.65101410, 0.70037330, 0.75461241, 0.81453466, 0.88162611,
    0.36133042, 0.37440405, 0.39079550, 0.41036291, 0.43302090, 0.45873872, 0.48752239, 0.51939798, 0.55439682, 0.59254133, 0.63383034, 0.67822265, 0.72561799, 0.77583430, 0.82858092,
    0.35439502, 0.36702034, 0.38260006, 0.40101786, 0.42215421, 0.44591174, 0.47220102, 0.50092041, 0.53193360, 0.56504392, 0.59996245, 0.63626628, 0.67334062, 0.71029469, 0.74583251,
    0.34791990, 0.36026613, 0.37522638, 0.39274277, 0.41269203, 0.43494673, 0.45937116, 0.48580434, 0.51403915, 0.54379790, 0.57470298, 0.60623985, 0.63770811, 0.66815393, 0.69627074,
    0.33621768, 0.34833639, 0.36246295, 0.37869462, 0.39694660, 0.41708490, 0.43894731, 0.46233761, 0.48701196, 0.51266165, 0.53889233, 0.56519854, 0.59093075, 0.61524994, 0.63706014,
    0.32597289, 0.33811759, 0.35176187, 0.36716097, 0.38429206, 0.40304427, 0.42326335, 0.44475696, 0.46728781, 0.49056178, 0.51421251, 0.53778186, 0.56069387, 0.58221783, 0.60141144,
    0.31695465, 0.32925197, 0.34262585, 0.35747187, 0.37383264, 0.39162945, 0.41072659, 0.43094633, 0.45206745, 0.47381721, 0.49585921, 0.51777702, 0.53905156, 0.55902802, 0.57686399
};

auto interpolate(real x, double x0, double x1, Vec<double> const& ypoints) -> real
{
    auto const N = ypoints.size();
    auto const dx = (x1 - x0)/(N - 1);

    auto const k = floor((x.val() - x0)/dx);
    auto const i = k + 2 < N ? k : N - 3;

    const auto xi0 = x0 + dx*i;
    const auto xi1 = xi0 + dx;
    const auto xi2 = xi1 + dx;
    const auto yi0 = ypoints[i];
    const auto yi1 = ypoints[i + 1];
    const auto yi2 = ypoints[i + 2];

    return interpolateQuadratic(x, xi0, xi1, xi2, yi0, yi1, yi2);
}

auto J0(real x) -> real
{
    if(0.0 <= x && x <= 1.0) return interpolate(x, 0.0, 1.0, J0region1);
    if(1.0 <= x && x <= 10.0) return interpolate(x, 1.0, 10.0, J0region2);
    if(10.0 <= x && x <= 100.0) return interpolate(x, 10.0, 100.0, J0region3);
    if(100.0 <= x && x <= 1000.0) return interpolate(x, 100.0, 1000.0, J0region4);
    if(1000.0 <= x && x <= 10000.0) return interpolate(x, 1000.0, 10000.0, J0region5);

    Exception exception;
    exception.error << "Cannot interpolate the Pitzer function J0(x) with the provided x = " << x << ".";
    exception.reason << "The value of x must be within the interval [0, 10000].";
    RaiseError(exception);

    return 0.0;
}

auto J1(real x) -> real
{
    if(0.0 <= x && x <= 1.0) return interpolate(x, 0.0, 1.0, J1region1);
    if(1.0 <= x && x <= 10.0) return interpolate(x, 1.0, 10.0, J1region2);
    if(10.0 <= x && x <= 100.0) return interpolate(x, 10.0, 100.0, J1region3);
    if(100.0 <= x && x <= 1000.0) return interpolate(x, 100.0, 1000.0, J1region4);
    if(1000.0 <= x && x <= 10000.0) return interpolate(x, 1000.0, 10000.0, J1region5);

    Exception exception;
    exception.error << "Cannot interpolate the Pitzer function J1(x) with the provided x = " << x << ".";
    exception.reason << "The value of x must be within the interval [0, 10000].";
    RaiseError(exception);

    return 0.0;
}

/// Auxiliary alias for ActivityModelParamsPitzer::CorrectionModel.
using PitzerParamCorrectionModel = ActivityModelParamsPitzer::CorrectionModel;

/// Return the Pitzer temperature-pressure correction model based on constant expression.
auto createParamCorrectionModelConstant(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    auto const& c = coefficients;

    errorif(c.size() == 0, "Cannot create the Constant temperature-dependent Pitzer parameter function with empty coefficients");

    return [=](real const& T, real const& Pbar) { return c[0]; };

    errorif(true, "Cannot create the Constant temperature-dependent Pitzer parameter function with given coefficients (which should be only one): ", str(c));
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by PHREEQC v3.
auto createParamCorrectionModelPhreeqc(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    auto const& c = coefficients;
    auto const Tr = 298.15;

    errorif(c.size() == 0, "Cannot create the Phreeqc temperature-dependent Pitzer parameter function with empty coefficients");

    if(c.size() == 1) return [=](real const& T, real const& Pbar) { return c[0]; };
    if(c.size() == 2) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr); };
    if(c.size() == 3) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr); };
    if(c.size() == 4) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr); };
    if(c.size() == 5) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr); };
    if(c.size() == 6) return [=](real const& T, real const& Pbar) { return c[0] + c[1]*(1.0/T - 1.0/Tr) + c[2]*log(T/Tr) + c[3]*(T - Tr) + c[4]*(T*T - Tr*Tr) + c[5]*(1.0/(T*T) - 1.0/(Tr*Tr)); };

    errorif(true, "Cannot create the Phreeqc temperature-dependent Pitzer parameter function with given coefficients: ", str(c));
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by He and Morse (1993) (doi: 10.1016/0016-7037(93)90137-L).
auto createParamCorrectionModelHeMorse1993(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `HeMorse1993` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Dai et al. (2013) (doi: 10.2118/164045-ms).
auto createParamCorrectionModelDai2013(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Dai2013` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Dai et al. (2014) (doi: 10.2118/169786-ms).
auto createParamCorrectionModelDai2014(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Dai2014` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Christov and Møller (2004) (see Table 5 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelChristovMoller2004(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `ChristovMoller2004` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Holmes et al. (1987) (see Table 6 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelHolmes1987(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Holmes1987` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Pitzer et al. (1984) (see Table 7 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelPitzer1984(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Pitzer1984` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Palaban and Pitzer (1987) (see Table 8 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelPalabanPitzer1987(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `PalabanPitzer1987` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Polya et al. (2001) (see Table 10 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelPolya2001(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `Polya2001` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model based on expression provided by Li and Duan (2007) (see Table 12 in 10.1016/j.chemgeo.2007.07.023).
auto createParamCorrectionModelLiDuan2007(Vec<Param> const& coefficients) -> Fn<real(real const&, real const&)>
{
    errorif(true, "Currently, the `LiDuan2007` temperature-pressure correction model for Pitzer parameters are not implemented");
}

/// Return the Pitzer temperature-pressure correction model with given coefficients and model option.
auto createParamCorrectionModel(Vec<Param> const& coefficients, PitzerParamCorrectionModel const& option) -> Fn<real(real const&, real const&)>
{
    switch(option)
    {
        case PitzerParamCorrectionModel::Constant:           return createParamCorrectionModelConstant(coefficients);
        case PitzerParamCorrectionModel::Phreeqc:            return createParamCorrectionModelPhreeqc(coefficients);
        case PitzerParamCorrectionModel::HeMorse1993:        return createParamCorrectionModelHeMorse1993(coefficients);
        case PitzerParamCorrectionModel::Dai2013:            return createParamCorrectionModelDai2013(coefficients);
        case PitzerParamCorrectionModel::Dai2014:            return createParamCorrectionModelDai2014(coefficients);
        case PitzerParamCorrectionModel::ChristovMoller2004: return createParamCorrectionModelChristovMoller2004(coefficients);
        case PitzerParamCorrectionModel::Holmes1987:         return createParamCorrectionModelHolmes1987(coefficients);
        case PitzerParamCorrectionModel::Pitzer1984:         return createParamCorrectionModelPitzer1984(coefficients);
        case PitzerParamCorrectionModel::PalabanPitzer1987:  return createParamCorrectionModelPalabanPitzer1987(coefficients);
        case PitzerParamCorrectionModel::Polya2001:          return createParamCorrectionModelPolya2001(coefficients);
        case PitzerParamCorrectionModel::LiDuan2007:         return createParamCorrectionModelLiDuan2007(coefficients);
        default:                                             return createParamCorrectionModelPhreeqc(coefficients);
    }
}

/// Used to represent an interaction parameter in the Pitzer activity model.
struct PitzerParam
{
    /// The indices of the species associated with this interaction parameter.
    const Indices ispecies;

    /// The temperature-pressure correction model for the interaction parameter.
    const Fn<real(real const&, real const&)> model;

    /// The current value of the interaction parameter since last update.
    real value;

    /// Update the current value of the interaction parameter.
    auto update(real const& T, real const& P)
    {
        value = model(T, P);
    }
};

/// Auxiliary alias for ActivityModelParamsPitzer::InteractionParamAttribs.
using PitzerInteractionParamAttribs = ActivityModelParamsPitzer::InteractionParamAttribs;

/// Return a given list of species indices sorted in ascending order of electric charge.
auto sortedSpeciesIndicesByCharge(SpeciesList const& specieslist, Indices const& ispecies) -> Indices
{
    // Lambda function that returns the charge of the ith species in `specieslist`
    auto chargefn = [&](auto i) { return specieslist[i].charge(); };

    // Lambda function that returns true if two species indices are in ascending order of electric charge.
    auto lesschargefn = [&](Index ispecies1, Index ispecies2) { return chargefn(ispecies1) < chargefn(ispecies2); };

    return sortedfn(ispecies, lesschargefn);
}

/// Convert a PitzerInteractionParamAttribs object to a PitzerParam one for a binary interaction parameter.
auto createPitzerParamBinary(SpeciesList const& specieslist, PitzerInteractionParamAttribs const& attribs) -> PitzerParam
{
    errorifnot(attribs.formulas.size() == 2, "Expecting two chemical formulas when processing a Pitzer binary interaction parameter but got formulas: ", str(attribs.formulas));

    auto const numspecies = specieslist.size();

    auto const ispecies1 = specieslist.findWithFormula(attribs.formulas[0]);
    if(ispecies1 >= numspecies)
        return {}; // species1 is not present in the list of aqueous species; return empty PitzerParam object

    auto const ispecies2 = specieslist.findWithFormula(attribs.formulas[1]);
    if(ispecies2 >= numspecies)
        return {}; // species2 is not present in the list of aqueous species; return empty PitzerParam object

    Indices ispecies = sortedSpeciesIndicesByCharge(specieslist, {ispecies1, ispecies2});

    return PitzerParam{ ispecies, createParamCorrectionModel(attribs.parameters, attribs.model) };
}

/// Convert a PitzerInteractionParamAttribs object to a PitzerParam one for a ternary interaction parameter.
auto createPitzerParamTernary(SpeciesList const& specieslist, PitzerInteractionParamAttribs const& attribs) -> PitzerParam
{
    errorifnot(attribs.formulas.size() == 3, "Expecting two chemical formulas when processing a Pitzer ternary interaction parameter but got formulas: ", str(attribs.formulas));

    auto const numspecies = specieslist.size();

    auto const ispecies1 = specieslist.findWithFormula(attribs.formulas[0]);
    if(ispecies1 >= numspecies)
        return {}; // species1 is not present in the list of aqueous species; return empty PitzerParam object

    auto const ispecies2 = specieslist.findWithFormula(attribs.formulas[1]);
    if(ispecies2 >= numspecies)
        return {}; // species2 is not present in the list of aqueous species; return empty PitzerParam object

    auto const ispecies3 = specieslist.findWithFormula(attribs.formulas[2]);
    if(ispecies3 >= numspecies)
        return {}; // species3 is not present in the list of aqueous species; return empty PitzerParam object

    Indices ispecies = sortedSpeciesIndicesByCharge(specieslist, {ispecies1, ispecies2, ispecies3});

    return PitzerParam{ ispecies, createParamCorrectionModel(attribs.parameters, attribs.model) };
}

/// Return the default value for \eq{alpha_1} parameter according to that used in PHREEQC v3 (see file pitzer.cpp under comment "Set alpha values").
auto determineDefaultAlpha1(ChemicalFormula const& formula0, ChemicalFormula const& formula1) -> Param
{
    auto const z0 = round(abs(formula0.charge()));
    auto const z1 = round(abs(formula1.charge()));
    return (z0 == 2 && z1 == 2) ? 1.4 : 2.0; // alpha1 = 2.0 for all electrolytes except the 2-2 electrolytes, for which alpha1 = 1.4 (Pitzer and Silvester, 1978; Plummer et al., 1988, page 4)
};

/// Return the default value for \eq{alpha_2} parameter according to that used in PHREEQC v3 (see file pitzer.cpp under comment "Set alpha values").
auto determineDefaultAlpha2(ChemicalFormula const& formula0, ChemicalFormula const& formula1) -> Param
{
    auto const z0 = round(abs(formula0.charge()));
    auto const z1 = round(abs(formula1.charge()));
    return (z0 == 1 && z1 == 1) || (z0 == 2 && z1 == 2) ? 12.0 : 50.0; // alpha2 = 12.0 for 1-1 and 2-2 electrolytes, alpha2 = 50.0 for 3-2 and 4-2 electrolytes (Pitzer and Silvester, 1978; Plummer et al., 1988, page 4)
};

/// Return the index of the entry in `alphas` containing a pair of formulas equivalent to the pair `formula0` and `formula1`.
auto findAlphaEntry(ChemicalFormula const& formula0, ChemicalFormula const& formula1, Vec<PitzerInteractionParamAttribs> const& alphas) -> Index
{
    return indexfn(alphas, RKT_LAMBDA(x,
        ( formula0.equivalent(x.formulas[0]) && formula1.equivalent(x.formulas[1]) ) ||
        ( formula0.equivalent(x.formulas[1]) && formula1.equivalent(x.formulas[0]) )
    ));
}

/// Return the value for \eq{alpha_1} parameter (either user-specified or determined from default value according to that used in PHREEQC v3).
auto determineAlpha1(ChemicalFormula const& formula0, ChemicalFormula const& formula1, Vec<PitzerInteractionParamAttribs> const& alphas) -> Param
{
    if(auto const idx = findAlphaEntry(formula0, formula1, alphas); idx < alphas.size())
        return alphas[idx].parameters[0];
    return determineDefaultAlpha1(formula0, formula1);
}

/// Return the value for \eq{alpha_2} parameter (either user-specified or determined from default value according to that used in PHREEQC v3).
auto determineAlpha2(ChemicalFormula const& formula0, ChemicalFormula const& formula1, Vec<PitzerInteractionParamAttribs> const& alphas) -> Param
{
    if(auto const idx = findAlphaEntry(formula0, formula1, alphas); idx < alphas.size())
        return alphas[idx].parameters[0];
    return determineDefaultAlpha2(formula0, formula1);
}

/// Determine the coefficients multiplying the terms where the lambda Pitzer
/// parameter is involved. Check equations (34-37) of Clegg (1991) (Activity
/// Coefficients in Natural Waters). Also check code in PHREEQC file pitzer.dat,
/// under comments:
///  * "Coef for Osmotic coefficient for TYPE_LAMDA"
auto determineLambdaCoeffs(double z1, double z2, Index i1, Index i2) -> Tuple<double, double, double>
{
    assert(z1 <= z2 && "Expecting species associated to binary lambda Pitzer parameter to be ordered in ascending order of charge.");

    const bool NxNx = z1 == 0.0 && z2 == 0.0 && i1 == i2;
    if(NxNx)
        return {1.0, 1.0, 0.5}; // (37.i => 2/2), (37.i => 2/2), (34.iv => 1/2)

    const bool NxNy = z1 == 0.0 && z2 == 0.0 && i1 != i2;
    if(NxNy)
        return {2.0, 2.0, 1.0}; // (37.i => 2), (37.i => 2), (34.v => 1)

    const bool NxCx = z1 == 0.0 && z2  > 0.0;
    if(NxCx)
        return {2.0, 2.0, 1.0}; // (37.iv => 2), (35.iv => 2), (34.vii => 1)

    const bool AxNx = z1  < 0.0 && z2 == 0.0;
    if(AxNx)
        return {2.0, 2.0, 1.0}; // (36.iv => 2), (37.iv => 2), (34.viii => 1)

    errorif(true, "Could not determine mu factor in osmotic coefficient equation with given triad of species with charges (", z1, ", ", z2, ").");

    return {};
}

/// Determine the coefficients multiplying the terms where the mu Pitzer
/// parameter is involved. Check equations (34-37) of Clegg (1991) (Activity
/// Coefficients in Natural Waters). Also check code in PHREEQC file pitzer.dat,
/// under comments:
///  * "Coef for Osmotic coefficient for TYPE_MU" and
///  * "Coef for gammas for TYPE_MU"
auto determineMuCoeffs(double z1, double z2, double z3, Index i1, Index i2, Index i3) -> Tuple<double, double, double, double>
{
    assert(z1 <= z2 && z2 <= z3 && "Expecting species associated to ternary mu Pitzer parameter to be ordered in ascending order of charge.");

    const bool NNN = z1 == 0.0 && z2 == 0.0 && z3 == 0.0;
    const bool NNC = z1 == 0.0 && z2 == 0.0 && z3  > 0.0;
    const bool ANN = z1  < 0.0 && z2 == 0.0 && z3 == 0.0;

    const bool NxNxNx = NNN && i1 == i2 && i2 == i3;
    if(NxNxNx)
        return {1.0, 1.0, 1.0, 1.0}; // (37.i => 3/3), (37.i => 3/3), (37.i => 3/3), (34.iv => 1)

    const bool NxNxNy = NNN && i1 == i2 && i2 != i3;
    if(NxNxNy)
        return {3.0, 3.0, 3.0, 3.0}; // (37.ii => 6/2), (37.ii => 6/2), (37.i => 3), (34.v => 3)

    const bool NyNxNx = NNN && i1 != i2 && i2 == i3;
    if(NyNxNx)
        return {3.0, 3.0, 3.0, 3.0}; // (37.i => 3), (37.ii => 6/2), (37.ii => 6/2), (34.v => 3)

    const bool NxNyNx = NNN && i1 != i2 && i1 == i3;
    if(NxNyNx)
        return {3.0, 3.0, 3.0, 3.0}; // (37.ii => 6/2), (37.i => 3), (37.ii => 6/2), (34.v => 3)

    const bool NxNyNz = NNN && i1 != i2 && i1 != i3;
    if(NxNyNz)
        return {6.0, 6.0, 6.0, 6.0}; // (37.iii => 6), (37.iii => 6), (37.iii => 6), (34.vi => 6)

    const bool NxNxCx = NNC && i1 == i2;
    if(NxNxCx)
        return {3.0, 3.0, 3.0, 3.0}; // (37.vii => 6/2), (37.vii => 6/2), (35.iv => 3), (34.vii => 3)

    const bool AxNxNx = ANN && i2 == i3;
    if(AxNxNx)
        return {3.0, 3.0, 3.0, 3.0}; // (36.iv => 3), (37.vii => 6/2), (37.vii => 6/2), (34.viii => 3)

    const bool NxNyCx = NNC && i1 != i2;
    if(NxNyCx)
        return {6.0, 6.0, 6.0, 6.0}; // (37.vii => 6), (37.vii => 6), (35.v => 6), (34.xii => 6)

    const bool AxNxNy = ANN && i2 != i3;
    if(AxNxNy)
        return {6.0, 6.0, 6.0, 6.0}; // (36.v => 6), (37.vii => 6/2), (37.vii => 6/2), (34.xiii => 6)

    errorif(true, "Could not determine mu factor in osmotic coefficient equation with given triad of species with charges (", z1, ", ", z2, ", ", z3, ").");

    return {};
}

/// The function \eq{g(x) = 2[1-(1+x)e^{-x}]/x^2} in the Pitzer model (see Eq. 11 of Plummer et al 1988).
auto G(real const& x) -> real
{
    return x == 0.0 ? x : 2.0*(1.0 - (1.0 + x)*exp(-x))/(x*x);
}

/// The function \eq{g^\prime(x) = -2\left[1-\left(1+x+\dfrac{1}{2}x^2\right)e^{-x}\right]/x^2} in the Pitzer model (see Eq. 12 of Plummer et al 1988).
auto GP(real const& x) -> real
{
    return x == 0.0 ? x : -2.0*(1.0 - (1.0 + x + 0.5*x*x)*exp(-x))/(x*x);
}

/// The function that computes electrostatic mixing effects of unsymmetrical cation-cation and anion-anion pairs \eq{^{E}\theta_{ij}(I)} and \eq{^{E}\theta_{ij}^{\prime}(I)}.
auto computeThetaValues(real const& I, real const& sqrtI, real const& Aphi, double zi, double zj) -> Pair<real, real>
{
    auto const aux = 6.0*Aphi*sqrtI;

    auto const xij = zi*zj*aux;
    auto const xii = zi*zi*aux;
    auto const xjj = zj*zj*aux;

    auto const J0ij = J0(xij);
    auto const J0ii = J0(xii);
    auto const J0jj = J0(xjj);

    auto const J1ij = J1(xij);
    auto const J1ii = J1(xii);
    auto const J1jj = J1(xjj);

    auto const thetaE = zi*zj/(4*I) * (J0ij - 0.5*J0ii - 0.5*J0jj);
    auto const thetaEP = zi*zj/(8*I*I) * (J1ij - 0.5*J1ii - 0.5*J1jj) - thetaE/I;

    return { thetaE, thetaEP };
}

/// The auxiliary type used to store computed values from the Pitzer activity model implemented by @ref Pitzer.
struct PitzerState
{
    real phiw;        ///< The osmotic coefficient of water.
    real ln_aw;       ///< The activity of water (natural log).
    ArrayXr ln_gamma; ///< The activity coefficients of the aqueous species (natural log).
};

/// The auxiliary type used to implement the Pitzer activity model.
struct PitzerModel
{
    AqueousMixture solution; ///< The aqueous solution for which this Pitzer activity model is defined.

    Vec<PitzerParam> beta0;  ///< The parameters \eq{\beta^{(0)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> beta1;  ///< The parameters \eq{\beta^{(1)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> beta2;  ///< The parameters \eq{\beta^{(2)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> Cphi;   ///< The parameters \eq{C^{\phi}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<PitzerParam> theta;  ///< The parameters \eq{\theta_{ij}(T, P)} in the Pitzer model for cation-cation and anion-anion interactions.
    Vec<PitzerParam> psi;    ///< The parameters \eq{\psi_{ijk}(T, P)} in the Pitzer model for cation-cation-anion and anion-anion-cation interactions.
    Vec<PitzerParam> lambda; ///< The parameters \eq{\lambda_{ij}(T, P)} in the Pitzer model for neutral-cation and neutral-anion interactions.
    Vec<PitzerParam> zeta;   ///< The parameters \eq{\zeta_{ijk}(T, P)} in the Pitzer model for neutral-cation-anion interactions.
    Vec<PitzerParam> mu;     ///< The parameters \eq{\mu_{ijk}(T, P)} in the Pitzer model for neutral-neutral-neutral, neutral-neutral-cation, and neutral-neutral-anion interactions.
    Vec<PitzerParam> eta;    ///< The parameters \eq{\eta_{ijk}(T, P)} in the Pitzer model for neutral-cation-cation and neutral-anion-anion interactions.

    Vec<Param> alpha1; ///< The parameters \eq{alpha_1_{ij}} associated to the parameters \eq{\beta^{(1)}_{ij}}.
    Vec<Param> alpha2; ///< The parameters \eq{alpha_2_{ij}} associated to the parameters \eq{\beta^{(1)}_{ij}}.

    using Tuples2i = Tuples<Index, Index>;                   ///< Auxiliary type for a tuple of 2 index values.
    using Tuples3d = Tuples<double, double, double>;         ///< Auxiliary type for a tuple of 3 double values.
    using Tuples4d = Tuples<double, double, double, double>; ///< Auxiliary type for a tuple of 4 double values.

    Tuples3d lambda_coeffs; ///< The coefficients multiplying the terms where the lambda Pitzer parameter is involved.
    Tuples4d mu_coeffs;     ///< The coefficients multiplying the terms where the mu Pitzer parameter is involved.

    Tuples2i thetaij; ///< The indices (i, j) of the cation-cation and anion-anion species pairs used to account for \eq{^{E}\theta_{ij}(I)} and \eq{^{E}\theta_{ij}^{\prime}(I)} contributions.
    ArrayXr thetaE;   ///< The current values of the parameters \eq{^{E}\theta_{ij}(I)} associated to the \eq{\theta_{ij}} parameters.
    ArrayXr thetaEP;  ///< The current values of the parameters \eq{^{E}\theta_{ij}^{\prime}(I)} associated to the \eq{\theta_{ij}} parameters.

    BilinearInterpolator Aphi; ///< The parameter \eq{A^\phi(T, P)} in the Pitzer model.

    /// Construct a default Pitzer object.
    PitzerModel()
    {}

    /// Construct a Pitzer object with given list of species in the aqueous solution and the parameters for the Pitzer activity model.
    // PitzerModel(SpeciesList const& specieslist, ActivityModelParamsPitzer const& params)
    PitzerModel(AqueousMixture const& solution, ActivityModelParamsPitzer const& params)
    : solution(solution)
    {
        for(auto const& entry : params.beta0)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                beta0.push_back(param);

        for(auto const& entry : params.beta1)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                beta1.push_back(param);

        for(auto const& entry : params.beta2)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                beta2.push_back(param);

        for(auto const& entry : params.Cphi)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                Cphi.push_back(param);

        for(auto const& entry : params.theta)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                theta.push_back(param);

        for(auto const& entry : params.psi)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                psi.push_back(param);

        for(auto const& entry : params.lambda)
            if(PitzerParam param = createPitzerParamBinary(solution.species(), entry); !param.ispecies.empty())
                lambda.push_back(param);

        for(auto const& entry : params.zeta)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                zeta.push_back(param);

        for(auto const& entry : params.mu)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                mu.push_back(param);

        for(auto const& entry : params.eta)
            if(PitzerParam param = createPitzerParamTernary(solution.species(), entry); !param.ispecies.empty())
                eta.push_back(param);

        for(auto const& entry : params.beta1)
            alpha1.push_back(determineAlpha1(entry.formulas[0], entry.formulas[1], params.alpha1));

        for(auto const& entry : params.beta2)
            alpha2.push_back(determineAlpha2(entry.formulas[0], entry.formulas[1], params.alpha2));

        auto const& ications = solution.indicesCations();
        auto const& ianions = solution.indicesAnions();

        for(auto i = 0; i < ications.size() - 1; ++i)
            for(auto j = i + 1; j < ications.size(); ++j)
                thetaij.emplace_back(ications[i], ications[j]);

        for(auto i = 0; i < ianions.size() - 1; ++i)
            for(auto j = i + 1; j < ianions.size(); ++j)
                thetaij.emplace_back(ianions[i], ianions[j]);

        auto const& z = solution.charges();

        for(auto const& entry : lambda)
        {
            auto const i1 = entry.ispecies[0];
            auto const i2 = entry.ispecies[1];
            lambda_coeffs.push_back(determineLambdaCoeffs(z[i1], z[i2], i1, i2));
        }

        for(auto const& entry : mu)
        {
            auto const i1 = entry.ispecies[0];
            auto const i2 = entry.ispecies[1];
            auto const i3 = entry.ispecies[2];
            mu_coeffs.push_back(determineMuCoeffs(z[i1], z[i2], z[i3], i1, i2, i3));
        }

        Aphi = BilinearInterpolator(Aphi_temperatures, Aphi_pressures, Aphi_data);
    }

    /// Update all Pitzer interaction parameters according to current temperature and pressure.
    auto updateParams(real const& T, real const& P)
    {
        // Note: Do not try to avoid the calculations before in case T and P is
        // the same as last time because T and P could be the same but the Param
        // objects in these models could be changing!

        auto const Pbar = P * 1e-5; // from Pa to bar

        for(auto& param : beta0)
            param.value = param.model(T, Pbar);

        for(auto& param : beta1)
            param.value = param.model(T, Pbar);

        for(auto& param : beta2)
            param.value = param.model(T, Pbar);

        for(auto& param : Cphi)
            param.value = param.model(T, Pbar);

        for(auto& param : theta)
            param.value = param.model(T, Pbar);

        for(auto& param : psi)
            param.value = param.model(T, Pbar);

        for(auto& param : lambda)
            param.value = param.model(T, Pbar);

        for(auto& param : zeta)
            param.value = param.model(T, Pbar);

        for(auto& param : mu)
            param.value = param.model(T, Pbar);

        for(auto& param : eta)
            param.value = param.model(T, Pbar);
    }

    /// Evaluate the Pitzer model and compute the properties of the aqueous solution.
    auto evaluate(AqueousMixtureState const& aqstate, PitzerState& pzstate)
    {
        auto const& T = aqstate.T; // in K
        auto const& P = aqstate.P; // in Pa
        auto const& M = aqstate.m; // in molal
        auto const& z = solution.charges();
        auto const& iH2O = solution.indexWater();
        auto const& icharged = solution.indicesCharged();
        auto const& Mw = waterMolarMass; // in kg/mol

        /// Update all Pitzer interaction parameters according to current temperature and pressure.
        updateParams(T, P);

        // The ionic strength of the solution and its square-root
        auto const I = aqstate.Ie;
        auto const DI = sqrt(I);

        auto const numspecies = M.size();

        auto& ln_aw  = pzstate.ln_aw = 0.0;
        auto& OSMOT  = pzstate.phiw = 0.0;
        auto& LGAMMA = pzstate.ln_gamma = zeros(numspecies);

        real CSUM = 0.0;
        real BIGZ = (M * z.abs()).sum();
        real OSUM = M.sum() - M[iH2O];

        // The Debye-Huckel coefficient Aphi evaluated at (T, P)
        auto const Aphi0 = Aphi(T, P);

        // The b parameter of the Pitzer model
	    auto const B = 1.2;

        // The F term in the Pitzer model
	    real F = -Aphi0*(DI/(1 + B*DI) + 2.0*log(1.0 + B*DI)/B);

        // The osmotic coefficient of water in the Pitzer model
        OSMOT = -Aphi0*I*DI/(1 + B*DI);

        for(auto const& param : beta0)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            LGAMMA[i0] += M[i1] * 2.0 * param.value;
            LGAMMA[i1] += M[i0] * 2.0 * param.value;
            OSMOT += M[i0] * M[i1] * param.value;
        }

        for(auto const& [i, param] : enumerate(beta1))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            F += M[i0] * M[i1] * param.value * GP(alpha1[i] * DI)/I;
            LGAMMA[i0] += M[i1] * 2.0 * param.value * G(alpha1[i] * DI);
            LGAMMA[i1] += M[i0] * 2.0 * param.value * G(alpha1[i] * DI);
            OSMOT += M[i0] * M[i1] * param.value * exp(-alpha1[i] * DI);
        }

        for(auto const& [i, param] : enumerate(beta2))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            F += M[i0] * M[i1] * param.value * GP(alpha2[i] * DI)/I;
            LGAMMA[i0] += M[i1] * 2.0 * param.value * G(alpha2[i] * DI);
            LGAMMA[i1] += M[i0] * 2.0 * param.value * G(alpha2[i] * DI);
            OSMOT += M[i0] * M[i1] * param.value * exp(-alpha2[i] * DI);
        }

        for(auto const& param : Cphi)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            auto const aux = 2.0 * sqrt(abs(z[i0] * z[i1]));

            CSUM += M[i0] * M[i1] * param.value/aux;
            LGAMMA[i0] += M[i1] * BIGZ * param.value/aux;
            LGAMMA[i1] += M[i0] * BIGZ * param.value/aux;
            OSMOT += M[i0] * M[i1] * BIGZ * param.value/aux;
        }

        for(auto const& param : theta)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            LGAMMA[i0] += 2.0 * M[i1] * param.value;
            LGAMMA[i1] += 2.0 * M[i0] * param.value;
            OSMOT += M[i0] * M[i1] * param.value;
        }

        for(auto const& [i0, i1] : thetaij)
        {
            auto const [etheta, ethetap] = computeThetaValuesInterpolation(I, DI, Aphi0, z[i0], z[i1]);

            F += M[i0] * M[i1] * ethetap;
            LGAMMA[i0] += 2.0 * M[i1] * etheta;
            LGAMMA[i1] += 2.0 * M[i0] * etheta;
            OSMOT += M[i0] * M[i1] * (etheta + I*ethetap);
        }

        for(auto const& param : psi)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        for(auto const& [i, param] : enumerate(lambda))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];

            auto const [clng0, clng1, cosm] = lambda_coeffs[i];

            LGAMMA[i0] += M[i1] * param.value * clng0;
            LGAMMA[i1] += M[i0] * param.value * clng1;
            OSMOT += M[i0] * M[i1] * param.value * cosm;
        }

        for(auto const& param : zeta)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        for(auto const& [i, param] : enumerate(mu))
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            auto const [clng0, clng1, clng2, cosm] = mu_coeffs[i];

            LGAMMA[i0] += M[i1] * M[i2] * param.value * clng0;
            LGAMMA[i1] += M[i0] * M[i2] * param.value * clng1;
            LGAMMA[i2] += M[i0] * M[i1] * param.value * clng2;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value * cosm;
        }

        for(auto const& param : eta)
        {
            auto const i0 = param.ispecies[0];
            auto const i1 = param.ispecies[1];
            auto const i2 = param.ispecies[2];

            LGAMMA[i0] += M[i1] * M[i2] * param.value;
            LGAMMA[i1] += M[i0] * M[i2] * param.value;
            LGAMMA[i2] += M[i0] * M[i1] * param.value;
            OSMOT += M[i0] * M[i1] * M[i2] * param.value;
        }

        // Finalise the calculation of the activity coefficient by adding the missing F and CSUM contributions
        for(auto i : icharged)
            LGAMMA[i] += z[i]*z[i]*F + abs(z[i])*CSUM;

        // Finalise the calculation of the activity coefficient by adding the missing F and CSUM contributions
        OSMOT = 1 + 2.0/OSUM * OSMOT;

        // Compute the activity of the water species
        ln_aw = -OSMOT * OSUM * Mw;
    }
};

} // namespace anonymous

auto createActivityModelPitzer(SpeciesList const& species, ActivityModelParamsPitzer const& params) -> ActivityModel
{
    // Create the aqueous solution
    AqueousMixture solution(species);

    // The index of water in the solution
    const Index iH2O = solution.indexWater();

    // Create the PitzerModel object responsible for evaluating the Pitzer model
    PitzerModel pzmodel(solution, params);

    // The PitzerState object that holds computed properties of the aqueous solution by the Pitzer model
    PitzerState pzstate;

    // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
    auto aqstateptr = std::make_shared<AqueousMixtureState>();
    auto aqsolutionptr = std::make_shared<AqueousMixture>(solution);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        // Evaluate the state of the aqueous solution
        auto const& aqstate = *aqstateptr = solution.state(T, P, x);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;

        // Export the aqueous solution and its state via the `extra` data member
        props.extra["AqueousMixtureState"] = aqstateptr;
        props.extra["AqueousMixture"] = aqsolutionptr;

        // Evaluate the Pitzer activity model with given aqueous state
        pzmodel.evaluate(aqstate, pzstate);

        // The mole fraction of water
        auto const xw = x[iH2O];

        // Set the activity coefficients of the solutes
        props.ln_g = pzstate.ln_gamma;

        // Set the activities of the solutes
        props.ln_a = props.ln_g + log(aqstate.m);

        // Set the activitiy of water
        props.ln_a[iH2O] = pzstate.ln_aw;

        // Set the activity coefficient of water (mole fraction scale)
        props.ln_g[iH2O] = pzstate.ln_aw - log(xw);
    };

    return fn;
}

auto ActivityModelPitzer() -> ActivityModelGenerator
{
    return ActivityModelPitzer(Params::embedded("Pitzer.yaml"));
}

auto ActivityModelPitzer(ActivityModelParamsPitzer const& params) -> ActivityModelGenerator
{
    return [=](SpeciesList const& species) { return createActivityModelPitzer(species, params); };
}

auto ActivityModelPitzer(Params const& params) -> ActivityModelGenerator
{
    auto const& data = params.data();
    errorif(!data.exists("ActivityModelParams"), "Expecting Pitzer activity model parameters in given Params object, but it lacks a `ActivityModelParams` section within which another section `Pitzer` should exist.");
    errorif(!data.at("ActivityModelParams").exists("Pitzer"), "Expecting Pitzer activity model parameters in given Params object, under the section `Pitzer`.");
    errorif(!data.at("ActivityModelParams").at("Pitzer").isDict(), "Expecting section `Pitzer` with Pitzer activity model parameters to be a dictionary.");

    ActivityModelParamsPitzer pzparams =
        data["ActivityModelParams"]["Pitzer"].as<ActivityModelParamsPitzer>();

    return ActivityModelPitzer(pzparams);
}

} // namespace Reaktoro
