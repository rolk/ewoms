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
#ifndef EWOMS_TRANSPORT_PROPERTIES_2P_HH
#define EWOMS_TRANSPORT_PROPERTIES_2P_HH

#include <ewoms/decoupled/common/transportproperties.hh>
#include <ewoms/decoupled/2p/2pproperties.hh>

/*!
 * \ingroup Saturation2p
 */
/*!
 * \file
 * \brief Specifies the properties for immiscible 2p transport
 */
namespace Ewoms
{
namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for transport part of a decoupled two-phase model
NEW_TYPE_TAG(TransportTwoP, INHERITS_FROM(Transport, DecoupledTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG( CapillaryFlux );         //!< The type of the diffusive part in a transport equation
NEW_PROP_TAG( GravityFlux );        //!< The type of a convective part in a transport equation
}
}

#endif
