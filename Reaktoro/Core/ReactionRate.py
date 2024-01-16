# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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


from reaktoro import *
from autodiff import real
from pytest import approx


def testReactionRate():

    rate = ReactionRate(1.23)
    assert rate.value() == 1.23
    assert rate.onEquationMode() == False

    rate = ReactionRate(real(2.34))
    assert rate.value() == 2.34
    assert rate.onEquationMode() == False

    rate = ReactionRate.enforce(123.4)
    assert rate.value() == 123.4
    assert rate.onEquationMode() == True

    # Testing assign method with double and real scalars
    rate.assign(1.0)
    assert rate.value() == 1.0
    assert rate.onEquationMode() == True

    rate.assign(real(2.0))
    assert rate.value() == 2.0
    assert rate.onEquationMode() == True

    # Testing arithmetic assign-operators with double and real scalars
    for seed in [1.0, real(1.0)]:
        rate = ReactionRate(2.0 * seed)

        rate += (1.0 * seed)
        assert rate.value() == approx(3.0)

        rate -= (1.0 * seed)
        assert rate.value() == approx(2.0)

        rate *= (3.0 * seed)
        assert rate.value() == approx(6.0)

        rate /= (3.0 * seed)
        assert rate.value() == approx(2.0)

    # Testing positive and negative arithmetic operators
    rate = ReactionRate(2.0)

    rate = +rate
    assert rate.value() == approx(2.0)

    rate = -rate
    assert rate.value() == approx(-2.0)

    # Testing arithmetic operators with double and real scalars on the right
    for seed in [1.0, real(1.0)]:
        rate = ReactionRate(2.0 * seed)

        rate = rate + (1.0 * seed)
        assert rate.value() == approx(3.0)

        rate = rate - (1.0 * seed)
        assert rate.value() == approx(2.0)

        rate = rate * (3.0 * seed)
        assert rate.value() == approx(6.0)

        rate = rate / (3.0 * seed)
        assert rate.value() == approx(2.0)

    # Testing arithmetic operators with double and real scalars on the left
    for seed in [1.0, real(1.0)]:
        rate = ReactionRate(2.0 * seed)

        rate = (1.0 * seed) + rate
        assert rate.value() == approx(3.0)

        rate = (1.0 * seed) - rate
        assert rate.value() == approx(2.0)

        rate = (3.0 * seed) * rate
        assert rate.value() == approx(6.0)

        rate = (3.0 * seed) / rate
        assert rate.value() == approx(2.0)

