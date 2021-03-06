// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
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
#ifndef EWOMS_GRAVITYPART_HH
#define EWOMS_GRAVITYPART_HH

#include <ewoms/decoupled/2p/transport/fv/convectivepart.hh>
#include "fvtransportproperties2p.hh"

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

/*!
 * \file
 * \copydoc Ewoms::GravityPart
 */

namespace Ewoms
{
/*!\ingroup FVSaturation2p
 * \brief  Class for defining the gravity term of  a two-phase flow saturation equation
 *
 * Defines the gravity term of the form
 *
 * \f[
 * \bar \lambda \boldsymbol K \, (\rho_n - \rho_w) \, g \, \text{grad} \, z,
 * \f]
 *
 * where \f$ \bar \lambda = \lambda_w f_n = \lambda_n f_w \f$ and \f$ \lambda \f$ is a phase mobility and \f$ f \f$ a phase fractional flow function,
 * \f$ \boldsymbol K \f$ is the intrinsic permeability, \f$ \rho \f$ is a phase density and  \f$ g \f$ is the gravity constant.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class GravityPart: public ConvectivePart<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView)GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

public:

    /*! \brief Returns convective term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           saturation of current element
     *  \param satJ           saturation of neighbor element
     */
    void getFlux(DimVector& flux, const Intersection& intersection, const Scalar satI, const Scalar satJ) const
    {
        ElementPointer element = intersection.inside();

        int globalIdxI = problem_.variables().index(*element);
        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        int indexInInside = intersection.indexInInside();

        //get lambda_bar = lambda_n*f_w
        Scalar lambdaWI = 0;
        Scalar lambdaNWI = 0;
        Scalar lambdaWJ = 0;
        Scalar lambdaNWJ = 0;

        if (preComput_)
        {
            lambdaWI=cellDataI.mobility(wPhaseIdx);
            lambdaNWI=cellDataI.mobility(nPhaseIdx);
        }
        else
        {
            lambdaWI = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satI);
            lambdaWI /= viscosity_[wPhaseIdx];
            lambdaNWI = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satI);
            lambdaNWI /= viscosity_[nPhaseIdx];
        }

        Scalar potentialW = cellDataI.fluxData().potential(wPhaseIdx, indexInInside);
        Scalar potentialNW = cellDataI.fluxData().potential(nPhaseIdx, indexInInside);

        DimMatrix meanPermeability(0);
        DimVector distVec(0);
        Scalar lambdaW =  0;
        Scalar lambdaNW = 0;

        if (intersection.neighbor())
        {
            // access neighbor
            ElementPointer neighborPointer = intersection.outside();

            int globalIdxJ = problem_.variables().index(*neighborPointer);
            CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

            distVec = neighborPointer->geometry().center() - element->geometry().center();

            // get permeability
            problem_.spatialParams().meanK(meanPermeability,
                    problem_.spatialParams().intrinsicPermeability(*element),
                    problem_.spatialParams().intrinsicPermeability(*neighborPointer));

            //get lambda_bar = lambda_n*f_w
            if (preComput_)
            {
                lambdaWJ=cellDataJ.mobility(wPhaseIdx);
                lambdaNWJ=cellDataJ.mobility(nPhaseIdx);
            }
            else
            {
                lambdaWJ = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*neighborPointer), satJ);
                lambdaWJ /= viscosity_[wPhaseIdx];
                lambdaNWJ = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*neighborPointer), satJ);
                lambdaNWJ /= viscosity_[nPhaseIdx];
            }

            lambdaW = (potentialW >= 0) ? lambdaWI : lambdaWJ;
            lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
            lambdaNW = (potentialNW >= 0) ? lambdaNWI : lambdaNWJ;
            lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;
        }
        else
        {
            // get permeability
            problem_.spatialParams().meanK(meanPermeability,
                    problem_.spatialParams().intrinsicPermeability(*element));

            distVec = intersection.geometry().center() - element->geometry().center();

            //calculate lambda_n*f_w at the boundary
            lambdaWJ = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satJ);
            lambdaWJ /= viscosity_[wPhaseIdx];
            lambdaNWJ = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satJ);
            lambdaNWJ /= viscosity_[nPhaseIdx];

            //If potential is zero always take value from the boundary!
            lambdaW = (potentialW > 0) ? lambdaWI : lambdaWJ;
            lambdaNW = (potentialNW > 0) ? lambdaNWI : lambdaNWJ;
        }

        // set result to K*grad(pc)
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();
        Scalar dist = distVec.two_norm();
        //calculate unit distVec
        distVec /= dist;
        Scalar areaScaling = (unitOuterNormal * distVec);

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        Scalar scalarPerm = permeability.two_norm();

        Scalar scalarGravity = problem_.gravity() * distVec;

        flux = unitOuterNormal;

        // set result to f_w*lambda_n*K*grad(pc)
        flux *= lambdaW*lambdaNW/(lambdaW+lambdaNW) * scalarPerm * (density_[wPhaseIdx] - density_[nPhaseIdx]) * scalarGravity * areaScaling;
    }
    /*! \brief Constructs a GravityPart object
     *
     *  \param problem A problem class object
     */
    GravityPart (Problem& problem)
    : ConvectivePart<TypeTag>(problem), problem_(problem), preComput_(GET_PROP_VALUE(TypeTag, PrecomputedConstRels))
    {}

    //! For initialization
    void initialize()
    {
        ElementIterator element = problem_.gridView().template begin<0> ();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(*element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(*element));
        fluidState.setTemperature(problem_.temperature(*element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);
    }

private:
    Problem& problem_; //problem data
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object, if preCompute = false new mobilities will be taken (for implicit Scheme)
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
};}

#endif
