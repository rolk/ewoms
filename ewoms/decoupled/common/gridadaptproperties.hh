// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Markus Wolff                                 *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef EWOMS_GRIDADAPT_PROPERTIES_HH
#define EWOMS_GRIDADAPT_PROPERTIES_HH

#include <ewoms/common/basicproperties.hh>

namespace Ewoms
{
namespace Properties
{
//! Grid adaption type tag for all decoupled models.
NEW_TYPE_TAG(GridAdaptTypeTag);

//! Defines if the grid is h-adaptive
NEW_PROP_TAG( EnableGridAdapt);

//! Class defining the refinement/coarsening indicator
NEW_PROP_TAG(GridAdaptIndicator);

//! Class defining the refinement/coarsening indicator for grid initialization
NEW_PROP_TAG(GridAdaptInitializationIndicator);

//! Switch the use of initial grid adaption on/off
NEW_PROP_TAG(GridAdaptEnableInitializationIndicator);

//! Mimimum allowed level
NEW_PROP_TAG(GridAdaptMinLevel);

//! Maximum allowed level
NEW_PROP_TAG(GridAdaptMaxLevel);

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineTolerance);

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenTolerance);

//! Tolerance for refinement
NEW_PROP_TAG(GridAdaptRefineThreshold);

//! Tolerance for coarsening
NEW_PROP_TAG(GridAdaptCoarsenThreshold);

//! Time step interval for adaption
NEW_PROP_TAG(GridAdaptInterval);

//! Switch for refinement at Dirichlet BC's -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtDirichletBC);

//! Switch for refinement at Neumann BC's -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtFluxBC);

//! Switch for refinement at sources -> not used by all indicators!
NEW_PROP_TAG(GridAdaptRefineAtSource);

//no adaptive grid
SET_BOOL_PROP(GridAdaptTypeTag, EnableGridAdapt, false);

//standard setting
SET_INT_PROP(GridAdaptTypeTag, GridAdaptMinLevel, 0);
SET_INT_PROP(GridAdaptTypeTag, GridAdaptMaxLevel, 1);
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptRefineTolerance, 0.05);
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptCoarsenTolerance, 0.001);
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptRefineThreshold, 0.0);
SET_SCALAR_PROP(GridAdaptTypeTag, GridAdaptCoarsenThreshold, 0.0);
SET_INT_PROP(GridAdaptTypeTag, GridAdaptInterval, 1);
//Switch initial grid adaption off per default
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptEnableInitializationIndicator, false);

// Switch of extra refinement strategy at boundaries/sources
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptRefineAtDirichletBC, false);
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptRefineAtFluxBC, false);
SET_BOOL_PROP(GridAdaptTypeTag, GridAdaptRefineAtSource, false);
} // namespace Properties
} // namespace Ewoms


#endif
