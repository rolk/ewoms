// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
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
#ifndef EWOMS_VELOCITYDEFAULT_HH
#define EWOMS_VELOCITYDEFAULT_HH

/*!
 * \file
 * \copydoc Ewoms::FVVelocityDefault
 */

#include <dune/grid/common/gridenums.hh>

#include <ewoms/decoupled/common/decoupledproperties.hh>

namespace Ewoms
{
//! \ingroup IMPET
/*! \brief Default implementation of a velocity class.
 *
 * If the velocity is reconstructed in the pressure model this default implementation is used in the transport model.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVVelocityDefault
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GridView::Intersection Intersection;


public:
    //! Constructs a FVVelocityDefault object
    /*!
     * \param problem A problem class object
     */
    FVVelocityDefault(Problem& problem)
    {}

    //! For initialization
    void initialize(){};

    //! Local velocity calculation
    void calculateVelocity(const Intersection& intersection, CellData& cellData)
    {}

    //! Local velocity calculation
    void calculateVelocityOnBoundary(const Intersection& intersection, CellData& cellData)
    {}

    //! \brief Indicates if velocity is reconstructed the transport step
    bool calculateVelocityInTransport()
    {
        return false;
    }

    /*! \brief Adds velocity output to the output file
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {}


};
}
#endif
