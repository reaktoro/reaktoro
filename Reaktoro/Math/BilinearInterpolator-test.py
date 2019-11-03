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

from reaktoro import BilinearInterpolator, VectorDouble

import pytest


def test_BilinearIntepolator():
    empty_interpolator = BilinearInterpolator()
    assert empty_interpolator.empty()

    assert not empty_interpolator.data()
    assert not empty_interpolator.xCoordinates()
    assert not empty_interpolator.yCoordinates()

    coordinates = VectorDouble([0, 1, 2, 3])
    data = VectorDouble([4, 5, 6, 7, 5, 5, 6, 7, 6, 6, 6, 7, 7, 7, 7, 7])
    empty_interpolator.setCoordinatesX(coordinates)
    empty_interpolator.setCoordinatesY(coordinates)
    empty_interpolator.setData(data)

    interpolator = BilinearInterpolator(coordinates, coordinates, data)
    assert not interpolator.empty()

    # it is needed to call data() twice to verify if the vector is not gonna be moved to interpolator_data and then destroyed
    interpolator_data = interpolator.data()
    assert all(i == j for i, j in zip(interpolator_data, data))
    assert all(i == j for i, j in zip(empty_interpolator.data(), data))
    assert interpolator.xCoordinates()
    assert interpolator.yCoordinates()
    assert interpolator(0.0, 0.0) == 4.0
    assert interpolator(0.5, 0.0) == 4.5
    assert interpolator(1.0, 1.0) == 5.0
    assert interpolator(2.0, 2.0) == 6.0
    assert interpolator(2.5, 2.5) == 6.75
    assert interpolator(2.99, 2.99) == pytest.approx(7, rel=0.1)
    assert interpolator(3.0, 3.0) == 7.0
    assert interpolator(3.0, 0.0) == interpolator(0.0, 3.0)
    assert interpolator.data()

    assert interpolator(3.0001, 1.0) == interpolator(3.0, 1.0)
    assert interpolator(1.0, 3.0001) == interpolator(1.0, 3.0)
    assert interpolator(-0.0001, 1.0) == interpolator(0.0, 1.0)
    assert interpolator(1.0, -0.0001) == interpolator(1.0, 0.0)
