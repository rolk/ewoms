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
#ifndef EWOMS_TRANSPORT_PROPERTIES_HH
#define EWOMS_TRANSPORT_PROPERTIES_HH

#include "decoupledproperties.hh"

/*!
 * \ingroup Sequential
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

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(Transport, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG( TransportSolutionType);
NEW_PROP_TAG( EvalCflFluxFunction ); //!< Type of the evaluation of the CFL-condition
NEW_PROP_TAG( ImpetCflFactor );
NEW_PROP_TAG( ImpetSwitchNormals );

SET_SCALAR_PROP(Transport, ImpetCflFactor, 1.0);
SET_BOOL_PROP(Transport, ImpetSwitchNormals, false);

/*!
 * \brief Default implementation for the Vector of the transportet quantity
 *
 * This type defines the data type of the transportet quantity. In case of a
 * immiscible 2p system, this would represent a vector holding the saturation
 * of one phase.
 */
SET_PROP(Transport, TransportSolutionType)
{
    private:
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionType;

    public:
    typedef typename SolutionType::ScalarSolution type;//!<type for vector of scalar properties
};
}
}

#endif
