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
#ifndef EWOMS_GRIDADAPTIONINDICATOR2P_HH
#define EWOMS_GRIDADAPTIONINDICATOR2P_HH

#include <ewoms/decoupled/common/impetproperties.hh>
#include <ewoms/decoupled/2p/2pproperties.hh>


/*!
 * \file
 * \copydoc Ewoms::GridAdaptionIndicator2P
 */
namespace Ewoms
{
/*!\ingroup IMPES
 * \brief  Class defining a standard, saturation dependent indicator for grid adaption
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class GridAdaptionIndicator2P
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };


public:
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, bool, GridAdaptEnableInitializationIndicator, "Enable the initialization indicator for grid adaption");

        REGISTER_PARAM(TypeTag, Scalar, GridAdaptRefineTolerance, "The maximum value of the grid adaption indicator before an element gets refined");
        REGISTER_PARAM(TypeTag, Scalar, GridAdaptCoarsenTolerance, "The minimum value of the grid adaption indicator before an element gets coarsened");
    }

    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     * This standard indicator is based on the saturation gradient.
     */
    void calculateIndicator()
    {
        // prepare an indicator for refinement
        if(indicatorVector_.size() != problem_.variables().cellDataGlobal().size())
        {
            indicatorVector_.resize(problem_.variables().cellDataGlobal().size());
        };
        indicatorVector_ = -1e100;

        Scalar globalMax = -1e100;
        Scalar globalMin = 1e100;

        ElementIterator eItEnd = problem_.gridView().template end<0>();
        // 1) calculate Indicator -> min, maxvalues
        // Schleife über alle Leaf-Elemente
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd;
                ++eIt)
        {
            // Bestimme maximale und minimale Sättigung
            // Index des aktuellen Leaf-Elements
            int globalIdxI = problem_.variables().index(*eIt);

            Scalar satI = 0.0;
            switch (saturationType_)
            {
            case Sw:
                satI = problem_.variables().cellData(globalIdxI).saturation(wPhaseIdx);
                break;
            case Sn:
                satI = problem_.variables().cellData(globalIdxI).saturation(nPhaseIdx);
                break;
            }

            globalMin = std::min(satI, globalMin);
            globalMax = std::max(satI, globalMax);

            // Berechne Verfeinerungsindikator an allen Zellen
            IntersectionIterator isItend = problem_.gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItend; ++isIt)
            {
                const typename IntersectionIterator::Intersection &intersection = *isIt;
                // Only consider internal intersections
                if (intersection.neighbor())
                {
                    // Access neighbor
                    ElementPointer outside = intersection.outside();
                    int globalIdxJ = problem_.variables().index(*outside);

                    // Visit intersection only once
                    if (eIt->level() > outside->level() || (eIt->level() == outside->level() && globalIdxI < globalIdxJ))
                    {
                        Scalar satJ = 0.;
                        switch (saturationType_)
                        {
                        case Sw:
                            satJ = problem_.variables().cellData(globalIdxJ).saturation(wPhaseIdx);
                            break;
                        case Sn:
                            satJ = problem_.variables().cellData(globalIdxJ).saturation(nPhaseIdx);
                            break;
                        }

                        Scalar localdelta = std::abs(satI - satJ);
                        indicatorVector_[globalIdxI][0] = std::max(indicatorVector_[globalIdxI][0], localdelta);
                        indicatorVector_[globalIdxJ][0] = std::max(indicatorVector_[globalIdxJ][0], localdelta);
                    }
                }
            }
        }

        Scalar globaldelta = globalMax - globalMin;

        refineBound_ = refinetol_*globaldelta;
        coarsenBound_ = coarsentol_*globaldelta;
    }

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        return (indicatorVector_[problem_.elementMapper().map(element)] > refineBound_);
    }

    /*! \brief Indicator function for marking of grid cells for coarsening
     *
     * Returns true if an element should be coarsened.
     *
     *  \param element A grid element
     */
    bool coarsen(const Element& element)
    {
        return (indicatorVector_[problem_.elementMapper().map(element)] < coarsenBound_);
    }

    /*! \brief Initializes the adaption indicator class*/
    void init()
    {
        refineBound_ = 0.;
        coarsenBound_ = 0.;
    }

    /*! \brief Constructs a GridAdaptionIndicator instance
     *
     *  This standard indicator is based on the saturation gradient. It checks the local gradient compared to the maximum global gradient.
     *  The indicator is compared locally to a refinement/coarsening threshold to decide whether a cell should be marked for refinement or coarsening or should not be adapted.
     *
     * \param problem The problem object
     */
    GridAdaptionIndicator2P (Problem& problem):
        problem_(problem)
    {
        refinetol_ = GET_PARAM(TypeTag, Scalar, GridAdaptRefineTolerance);
        coarsentol_ = GET_PARAM(TypeTag, Scalar, GridAdaptCoarsenTolerance);
    }

private:
    Problem& problem_;
    Scalar refinetol_;
    Scalar coarsentol_;
    Scalar refineBound_;
    Scalar coarsenBound_;
    ScalarSolutionType indicatorVector_;
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
};
}

#endif
