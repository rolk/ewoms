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
#ifndef EWOMS_GRIDADAPTINITIALIZATIONINDICATOR_HH
#define EWOMS_GRIDADAPTINITIALIZATIONINDICATOR_HH

#include "decoupledproperties.hh"

#include <dune/common/dynvector.hh>

/*!
 * \file
 * \copydoc Ewoms::GridAdaptInitializationIndicator
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
class GridAdaptInitializationIndicator
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, GridAdaptIndicator) GridAdaptIndicator;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    enum
    {
        refineCell = 1,
        coarsenCell = -1
    };


public:
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, bool, GridAdaptEnableInitializationIndicator, "Enable the initialization indicator for grid adaption");

        REGISTER_PARAM(TypeTag, bool, GridAdaptRefineAtDirichletBC, "Refine the grid at Dirichlet boundaries");
        REGISTER_PARAM(TypeTag, bool, GridAdaptRefineAtFluxBC, "Refine the grid at flux boundaries");
        REGISTER_PARAM(TypeTag, bool, GridAdaptRefineAtSource, "Refine the grid at regions with a source term");
    }

    /*! \brief Calculates the indicator used for refinement/coarsening for each grid cell.
     *
     */
    void calculateIndicator()
    {
        adaptionIndicator_.calculateIndicator();

        // prepare an indicator for refinement
        indicatorVector_.resize(problem_.variables().cellDataGlobal().size());

        indicatorVector_ = coarsenCell;

        if (!enableInitializationIndicator_)
            return;

        ElementIterator eItEnd = problem_.gridView().template end<0>();
        // 1) calculate Indicator -> min, maxvalues
        // Schleife über alle Leaf-Elemente
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd;
                ++eIt)
        {
            // Bestimme maximale und minimale Sättigung
            // Index des aktuellen Leaf-Elements
            int globalIdxI = problem_.variables().index(*eIt);

            int level = eIt->level();
            maxLevel_ = std::max(level, maxLevel_);

            if (level < minAllowedLevel_)
            {
                indicatorVector_[globalIdxI] = refineCell;
                continue;
            }

            if (refineAtSource_)
            {
            PrimaryVariables source(0.0);
            problem_.source(source, *eIt);
            for (int i = 0; i < numEq; i++)
            {
                if (std::abs(source[i]) > 1e-10)
                {
                    indicatorVector_[globalIdxI] = refineCell;
                    break;
                }
            }
            }

            if (indicatorVector_[globalIdxI] == refineCell)
                continue;

            // Berechne Verfeinerungsindikator an allen Zellen
            IntersectionIterator isItend = problem_.gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItend; ++isIt)
            {
                if (indicatorVector_[globalIdxI] == refineCell)
                    break;
                if (isIt->boundary())
                {
                    if (maxLevel_ != maxAllowedLevel_)
                        {
                            indicatorVector_[globalIdxI] = refineCell;
                        }
                    else
                    {
                        BoundaryTypes bcTypes;
                        problem_.boundaryTypes(bcTypes, *isIt);

                        for (int i = 0; i < numEq; i++)
                        {
                            if (bcTypes.isNeumann(i))
                            {
                                PrimaryVariables flux(0.0);
                                problem_.neumann(flux, *isIt);

                                bool fluxBound = false;
                                for (int j = 0; j < numPhases; j++)
                                {
                                if (std::abs(flux[j]) < 1e-10)
                                {
                                    indicatorVector_[globalIdxI] = coarsenCell;
                                }
                                else
                                {
                                    if (refineAtFluxBC_)
                                    {
                                        indicatorVector_[globalIdxI] = refineCell;
                                        fluxBound = true;
                                        break;
                                    }
                                    else
                                    {
                                        indicatorVector_[globalIdxI] = coarsenCell;
                                    }
                                }
                                }
                                if (fluxBound)
                                    break;
                            }
                            else if (bcTypes.isDirichlet(i))
                            {
                                if (refineAtDirichletBC_)
                                {
                                    indicatorVector_[globalIdxI] = refineCell;
                                    break;
                                }
                                else
                                {
                                    indicatorVector_[globalIdxI] = coarsenCell;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*! \brief Indicator function for marking of grid cells for refinement
     *
     * Returns true if an element should be refined.
     *
     *  \param element A grid element
     */
    bool refine(const Element& element)
    {
        int idx = problem_.elementMapper().map(element);
        if (indicatorVector_[idx] == refineCell || adaptionIndicator_.refine(element))
            return true;
        else
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
        int idx = problem_.elementMapper().map(element);
        if (maxLevel_ == maxAllowedLevel_ && indicatorVector_[idx] == coarsenCell && !adaptionIndicator_.refine(element))
            return true;
        else if (indicatorVector_[idx] == coarsenCell &&  adaptionIndicator_.coarsen(element))
            return true;
        else
            return false;
    }

    int maxLevel()
    {
        return maxLevel_;
    }

    /*! \brief Initializes the adaption indicator class */
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
    GridAdaptInitializationIndicator(Problem& problem, GridAdaptIndicator& adaptionIndicator):
        problem_(problem), adaptionIndicator_(adaptionIndicator), maxLevel_(0)
    {
        minAllowedLevel_ = GET_PARAM(TypeTag, int, GridAdaptMinLevel);
        maxAllowedLevel_ = GET_PARAM(TypeTag, int, GridAdaptMaxLevel);
        enableInitializationIndicator_ = GET_PARAM(TypeTag, bool, GridAdaptEnableInitializationIndicator);
        refineAtDirichletBC_ = GET_PARAM(TypeTag, bool, GridAdaptRefineAtDirichletBC);
        refineAtFluxBC_ = GET_PARAM(TypeTag, bool, GridAdaptRefineAtFluxBC);
        refineAtSource_ = GET_PARAM(TypeTag, bool, GridAdaptRefineAtSource);
    }

private:
    Problem& problem_;
    GridAdaptIndicator& adaptionIndicator_;
    Dune::DynamicVector<int> indicatorVector_;
    int maxLevel_;
    int minAllowedLevel_;
    int maxAllowedLevel_;
    bool enableInitializationIndicator_;
    bool refineAtDirichletBC_;
    bool refineAtFluxBC_;
    bool refineAtSource_;
};
}
#endif
