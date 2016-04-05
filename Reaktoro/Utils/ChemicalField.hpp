// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// A type that contains the values of a scalar field and its derivatives.
class ChemicalField
{
public:
    /// Construct a default ChemicalField instance.
    ChemicalField();

    /// Resize the chemical field.
    auto resize(unsigned num_points, unsigned num_components) -> void;

    /// Set the value of the chemical field at the i-th point.
    auto val(Index i, double val) -> void;

    /// Return a reference to the values of the scalar chemical field.
    auto val() -> Vector&;

    /// Return a const reference to the values of the scalar chemical field.
    auto val() const -> const Vector&;

    /// Activate the calculation of derivatives w.r.t. temperature.
    auto ddt(bool active) -> void;

    /// Set the temperature derivative of the chemical field at the i-th point.
    auto ddt(Index i, double ddt) -> void;

    /// Return a reference to the derivatives of the scalar chemical field w.r.t. temperature.
    auto ddt() -> Vector&;

    /// Return a const reference to the derivatives of the scalar chemical field w.r.t. temperature.
    auto ddt() const -> const Vector&;

    /// Activate the calculation of derivatives w.r.t. pressure.
    auto ddp(bool active) -> void;

    /// Set the pressure derivative of the chemical field at the i-th point.
    auto ddp(Index i, double ddp) -> void;

    /// Return a reference to the derivatives of the scalar chemical field w.r.t. pressure.
    auto ddp() -> Vector&;

    /// Return a const reference to the derivatives of the scalar chemical field w.r.t. pressure.
    auto ddp() const -> const Vector&;

    /// Activate the calculation of derivatives w.r.t. component amounts.
    auto ddc(bool active) -> void;

    /// Set the composition derivative of the chemical field at the i-th point.
    auto ddc(Index i, const Vector& ddc) -> void;

    /// Return a const reference to the derivatives of the scalar chemical field w.r.t. component amounts.
    auto ddc() -> Matrix&;

    /// Return a const reference to the derivatives of the scalar chemical field w.r.t. component amounts.
    auto ddc() const -> const Matrix&;

private:
    /// The values of the scalar chemical field.
    Vector m_val;

    /// The sensitivity of the scalar chemical field with respect to temperature.
    Vector m_ddt;

    /// The sensitivity of the scalar chemical field with respect to pressure.
    Vector m_ddp;

    /// The sensitivity of the scalar chemical field with respect to component amounts.
    Matrix m_ddc;
};

} // namespace Reaktoro

