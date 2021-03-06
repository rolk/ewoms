// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
#ifndef EWOMS_BOUNDARYCONDITIONS_HH
#define EWOMS_BOUNDARYCONDITIONS_HH

/*!
* \file
* \copydoc Ewoms::BoundaryConditions
*/
namespace Ewoms
{
/*!
 * \ingroup Sequential
 */
/*!
* \brief Define a class containing boundary condition flags
*
*/

//! base Class that defines boundary condition flags
struct BoundaryConditions
{
    /** \brief These values are ordered according to precedence */
    enum Flags {
        couplingOutflow = -2, //!< An outflow boundary for coupled models
        couplingInflow = -1, //!< An inflow boundary for coupled models
        outflow = 0, //!< An outflow boundary
        neumann = 1, //!< Neumann boundary
        process = 2, //!< Processor boundary
        dirichlet = 3 //!< Dirichlet boundary
    };
};

/** \} */
}
#endif
