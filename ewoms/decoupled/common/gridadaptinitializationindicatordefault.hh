// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *   Copyright (C) 2012 by Christoph Grueninger                              *
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
#ifndef EWOMS_GRIDADAPTINITIALIZATIONINDICATORDEFAULT_HH
#define EWOMS_GRIDADAPTINITIALIZATIONINDICATORDEFAULT_HH

#include "decoupledproperties.hh"

#include <dune/common/dynvector.hh>

/*!
 * \file
 * \copydoc Ewoms::GridAdaptInitializationIndicatorDefault
 */
namespace Ewoms
{
/*!\ingroup IMPES
 * \brief  Class defining a start indicator for grid adaption
 *
 *  Uses the defined grid adaptation indicator and further accounts for sources and boundaries.
 *  Only for grid initialization!
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptInitializationIndicatorDefault
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridAdaptIndicator) GridAdaptIndicator;

public:
    static void registerParameters()
    { }

    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {}

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        return false;
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
        return false;
    }

    int maxLevel()
    {
        return maxLevel_;
    }


    /*! \brief Initializes the adaption indicator class*/
    void init()
    {}

    /*! \brief Constructs a GridGridAdaptIndicator instance
     *
     * This standard indicator is based on the saturation gradient. It checks the local gradient
     * compared to the maximum global gradient. The indicator is compared locally to a
     * refinement/coarsening threshold to decide whether a cell should be marked for refinement
     * or coarsening or should not be adapted.
     *
     * \param problem The problem object
     * \param adaptionIndicator Indicator whether a be adapted
     */
    GridAdaptInitializationIndicatorDefault(Problem& problem, GridAdaptIndicator& adaptionIndicator)
    {
        maxLevel_ = GET_PARAM(TypeTag, int, GridAdaptMaxLevel);
    }

private:
    int maxLevel_;
};
}

#endif
