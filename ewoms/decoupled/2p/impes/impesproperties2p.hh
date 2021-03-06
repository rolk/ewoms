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
#ifndef EWOMS_IMPES2P_PROPERTIES_HH
#define EWOMS_IMPES2P_PROPERTIES_HH

#include <ewoms/decoupled/common/impetproperties.hh>
#include <ewoms/decoupled/2p/2pproperties.hh>

/*!
 * \ingroup IMPES
 */
/*!
 * \file
 * \brief Properties related to the sequential IMPES algorithms
 */
namespace Ewoms
{

namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPES algorithms
 *
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! TypeTag for the two-phase IMPES scheme
NEW_TYPE_TAG(IMPESTwoP, INHERITS_FROM(IMPET, DecoupledTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
}
}

#include <ewoms/decoupled/2p/impes/impesproblem2p.hh>

namespace Ewoms
{
namespace Properties
{
}
}

#endif
