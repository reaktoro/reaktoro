/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// Reaktor includes
#include <Reaktor/Activity/AqueousActivity.hpp>
#include <Reaktor/Activity/AqueousActivityDrummond.hpp>
#include <Reaktor/Activity/AqueousActivityDuanSun.hpp>
#include <Reaktor/Activity/AqueousActivityHKF.hpp>
#include <Reaktor/Activity/AqueousActivityIdeal.hpp>
#include <Reaktor/Activity/AqueousActivityPitzer.hpp>
#include <Reaktor/Activity/AqueousActivityRumpf.hpp>
#include <Reaktor/Activity/AqueousActivitySetschenow.hpp>
#include <Reaktor/Activity/GaseousActivity.hpp>
#include <Reaktor/Activity/GaseousActivityDuanMollerWeare.hpp>
#include <Reaktor/Activity/GaseousActivityDuanSun.hpp>
#include <Reaktor/Activity/GaseousActivityIdeal.hpp>
#include <Reaktor/Activity/GaseousActivityPengRobinson.hpp>
#include <Reaktor/Activity/GaseousActivitySpycherPruess.hpp>
#include <Reaktor/Activity/GaseousActivitySpycherReed.hpp>
#include <Reaktor/Activity/MineralActivity.hpp>
#include <Reaktor/Activity/MineralActivityIdeal.hpp>

/**
 * @defgroup Activity Activity
 *
 * Provides the models for activity calculations of aqueous, gaseous and mineral species
 */
