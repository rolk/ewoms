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
#ifndef EWOMS_CAPILLARYDIFFUSION_HH
#define EWOMS_CAPILLARYDIFFUSION_HH

#include <ewoms/decoupled/2p/transport/fv/diffusivepart.hh>
#include "fvtransportproperties2p.hh"

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

/*!
 * \file
 * \copydoc Ewoms::CapillaryDiffusion
 */
namespace Ewoms
{
/*!\ingroup FVSaturation2p
 * \brief  Class for defining the diffusive capillary pressure term of a saturation equation
 *
 * Defines the diffusive capillary pressure term of the form
 *
 * \f[
 * \bar \lambda \boldsymbol K \text{grad} \, p_c,
 * \f]
 *
 * where \f$ \bar \lambda = \lambda_w f_n = \lambda_n f_w \f$ and \f$ \lambda \f$ is a phase mobility and \f$ f \f$ a phase fractional flow function,
 * \f$ \boldsymbol K \f$ is the intrinsic permeability and \f$ p_c = p_c(S_w) \f$ the capillary pressure.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class CapillaryDiffusion: public DiffusivePart<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
      typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

      typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
      typedef typename SpatialParams::MaterialLaw MaterialLaw;

      typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
      typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

      typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

      typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
      typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
      typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        pressEqIdx = Indices::pressEqIdx
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

public:

    //! Returns capillary diffusion term
    /*! Returns capillary diffusion term for current element face
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           saturation of current element
     *  \param satJ           saturation of neighbor element
     *  \param pcGradient     gradient of capillary pressure between element I and J
     */
    void getFlux (DimVector& flux, const Intersection& intersection, Scalar satI, Scalar satJ, const DimVector& pcGradient) const
    {
        ElementPointer element = intersection.inside();
        // get global coordinate of cell center
        const GlobalPosition& globalPos = element->geometry().center();

        int globalIdxI = problem_.variables().index(*element);
        CellData& CellDataI = problem_.variables().cellData(globalIdxI);

        // get geometry type of face
        //Dune::GeometryType faceGT = isIt->geometryInInside().type();

        Scalar temperature = problem_.temperature(*element);
        Scalar referencePressure = problem_.referencePressure(*element);

        //get lambda_bar = lambda_n*f_w
        Scalar mobBar = 0;
        Scalar mobilityWI = 0;
        Scalar mobilityNWI = 0;

        if (preComput_)
        {
            mobilityWI = CellDataI.mobility(wPhaseIdx);
            mobilityNWI = CellDataI.mobility(nPhaseIdx);
        }
        else
        {
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, referencePressure);
            fluidState.setPressure(nPhaseIdx, referencePressure);
            fluidState.setTemperature(temperature);
            mobilityWI = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satI);
            mobilityWI /= FluidSystem::viscosity(fluidState, wPhaseIdx);
            mobilityNWI = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satI);
            mobilityNWI /= FluidSystem::viscosity(fluidState, nPhaseIdx);
        }

        DimMatrix meanPermeability(0);

        if (intersection.neighbor())
        {
            // access neighbor
            ElementPointer neighborPointer = intersection.outside();

            int globalIdxJ = problem_.variables().index(*neighborPointer);
            CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

            // distance vector between barycenters
            DimVector distVec = globalPosNeighbor - globalPos;

            // compute distance between cell centers
            Scalar dist = distVec.two_norm();

            DimVector unitDistVec(distVec);
            unitDistVec /= dist;

            // get permeability
            problem_.spatialParams().meanK(meanPermeability,
                    problem_.spatialParams().intrinsicPermeability(*element),
                    problem_.spatialParams().intrinsicPermeability(*neighborPointer));


            Scalar mobilityWJ = 0;
            Scalar mobilityNWJ = 0;
            //get lambda_bar = lambda_n*f_w
            if(preComput_)
            {
                mobilityWJ = cellDataJ.mobility(wPhaseIdx);
                mobilityNWJ = cellDataJ.mobility(nPhaseIdx);
            }
            else
            {
                FluidState fluidState;
                fluidState.setPressure(wPhaseIdx, referencePressure);
                fluidState.setPressure(nPhaseIdx, referencePressure);
                fluidState.setTemperature(temperature);

                mobilityWJ = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*neighborPointer), satJ);
                mobilityWJ /= FluidSystem::viscosity(fluidState, wPhaseIdx);
                mobilityNWJ = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*neighborPointer), satJ);
                mobilityNWJ /= FluidSystem::viscosity(fluidState, nPhaseIdx);
            }
            Scalar mobilityWMean = 0.5*(mobilityWI + mobilityWJ);
            Scalar mobilityNWMean = 0.5*(mobilityNWI + mobilityNWJ);
            mobBar = mobilityWMean*mobilityNWMean/(mobilityWMean+mobilityNWMean);
         }//end intersection with neighbor
        else
        {
            BoundaryTypes bcTypes;
            problem_.boundaryTypes(bcTypes, intersection);
            if (bcTypes.isNeumann(pressEqIdx))
            {
                PrimaryVariables priVars;
                problem_.neumann(priVars, intersection);
                if (priVars[wPhaseIdx] == 0)
                {
                    flux = 0;
                    return;
                }
            }
            // get permeability
            problem_.spatialParams().meanK(meanPermeability,
                    problem_.spatialParams().intrinsicPermeability(*element));

            Scalar mobilityWJ = 0;
            Scalar mobilityNWJ = 0;

            //calculate lambda_n*f_w at the boundary
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, referencePressure);
            fluidState.setPressure(nPhaseIdx, referencePressure);
            fluidState.setTemperature(temperature);
            mobilityWJ = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satJ);
            mobilityWJ /= FluidSystem::viscosity(fluidState, wPhaseIdx);
            mobilityNWJ = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satJ);
            mobilityNWJ /= FluidSystem::viscosity(fluidState, nPhaseIdx);

            Scalar mobWMean = 0.5 * (mobilityWI + mobilityWJ);
            Scalar mobNWMean = 0.5 * (mobilityNWI + mobilityNWJ);

            mobBar = mobWMean * mobNWMean / (mobWMean + mobNWMean);
        }

        // set result to K*grad(pc)
        meanPermeability.mv(pcGradient, flux);

        // set result to f_w*lambda_n*K*grad(pc)
        flux *= mobBar;
    }

    /*! \brief Constructs a CapillaryDiffusion object
     *
     *  \param problem A problem class object
     */
    CapillaryDiffusion (Problem& problem)
    : DiffusivePart<TypeTag>(problem), problem_(problem), preComput_(GET_PROP_VALUE(TypeTag, PrecomputedConstRels))
    {}

private:
    Problem& problem_;//problem data
    const bool preComput_;//if preCompute = true the mobilities are taken from the variable object, if preCompute = false new mobilities will be taken (for implicit Scheme)
};
}

#endif
