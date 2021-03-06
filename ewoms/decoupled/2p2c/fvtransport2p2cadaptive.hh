// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Benjamin Faigle                                   *
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
#ifndef EWOMS_FVTRANSPORT2P2C_ADAPTIVE_HH
#define EWOMS_FVTRANSPORT2P2C_ADAPTIVE_HH

#include <dune/grid/common/gridenums.hh>
#include <ewoms/decoupled/2p2c/2p2cadaptiveproperties.hh>
#include <ewoms/decoupled/2p2c/fvtransport2p2c.hh>
#include <ewoms/common/math.hh>
#include <ewoms/parallel/gridcommhandles.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

/*!
 * \file
 * \copydoc Ewoms::FVTransport2P2CAdaptive
 */

namespace Ewoms
{
//! Compositional Transport step in a Finite Volume discretization
/*! \ingroup Adaptive2p2c
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2CAdaptive : public FVTransport2P2C<TypeTag>
{
    typedef FVTransport2P2C<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;


    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;


    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld,
        NumPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx=Indices::contiWEqIdx, contiNEqIdx=Indices::contiNEqIdx
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;
    typedef Dune::FieldVector<Scalar,dim+1> TransmissivityMatrix;
    typedef Dune::FieldVector<Scalar, NumPhases> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    //! Acess function for the current problem
    Problem& problem()
    {return problem_;};

public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes);

    void getMpfaFlux(Dune::FieldVector<Scalar, 2>&, Dune::FieldVector<Scalar, 2>&,
            const IntersectionIterator&, CellData&);

    //! Constructs a FVTransport2P2C object
    /*!
     * Currently, the miscible transport scheme can not be applied with a global pressure / total velocity
     * formulation.
     *
     * \param problem a problem class object
     */
    FVTransport2P2CAdaptive(Problem& problem) : FVTransport2P2C<TypeTag>(problem),
        problem_(problem)
    {
        enableMPFA = GET_PARAM(TypeTag,bool,EnableMultiPointFluxApproximation);
    }

    virtual ~FVTransport2P2CAdaptive()
    {
    }

    static void registerParameters()
    {
        ParentType::registerParameters();

        REGISTER_PARAM(TypeTag, bool, EnableMultiPointFluxApproximation, "Specify whether the multi-point flux appoximation scheme ought to be used instead of two-point flux approximation");
        REGISTER_PARAM(TypeTag, bool, MpfaEnableSecondHalfEdge, "Specify whether the flux for the second half-edge of the MPFA should be calculated");
    }

protected:
    Problem& problem_;

    bool enableMPFA;
    static const int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
};

//! \brief Calculate the update vector and determine timestep size
/*!
 *  This method calculates the update vector \f$ u \f$ of the discretized equation
 *  \f[
       C^{\kappa , new} = C^{\kappa , old} + u,
 *  \f]
 *  where \f$ u = \sum_{\gamma} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{\gamma} \f$,
 *  \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{\gamma} \f$ is the face area of face \f$ \gamma \f$.
 *
 *  In addition to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition. This method uses a standard \a Tpfa method for regular fluxes,
 *  and a \a mpfa can be used near hanging nodes.
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param[out] dt Time step size \f$\mathrm{[s]}\f$
 *  \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2CAdaptive<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet)
{
    this->impet_ = impet;
    // initialize dt very large
    dt = 1E100;
    this->averagedFaces_ = 0;

    // resize update vector and set to zero
    int size_ = problem_.gridView().size(0);
    updateVec.resize(GET_PROP_VALUE(TypeTag, NumComponents));
    updateVec[wCompIdx].resize(size_);
    updateVec[nCompIdx].resize(size_);
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;
    //also resize PV vector if necessary
    if(int(this->totalConcentration_.size()) != size_)
    {
        this->totalConcentration_[wCompIdx].resize(size_);
        this->totalConcentration_[nCompIdx].resize(size_);
        // copy data //TODO: remove this, remove PM Transport Vector!!
        // loop through all elements
        for (int i = 0; i< problem().gridView().size(0); i++)
        {
            CellData& cellDataI = problem().variables().cellData(i);
            for(int compIdx = 0; compIdx < GET_PROP_VALUE(TypeTag, NumComponents); compIdx++)
            {
                this->totalConcentration_[compIdx][i]
                        = cellDataI.totalConcentration(compIdx);
            }
        }
    }

    // Cell which restricts time step size
    int restrictingCell = -1;

    Dune::FieldVector<Scalar, 2> entries(0.), timestepFlux(0.);
    // compute update vector
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get cell infos
        int globalIdxI = problem().variables().index(*eIt);
//        const GlobalPosition globalPos = eIt->geometry().center();
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        // some variables for time step calculation
        double sumfactorin = 0;
        double sumfactorout = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // handle interior face
            if (isIt->neighbor())
            {
                if (enableMPFA && isIt->outside()->level()!=eIt.level())
                    getMpfaFlux(entries, timestepFlux, isIt, cellDataI);
                else
                    this->getFlux(entries, timestepFlux, *isIt, cellDataI);
            }

            //     Boundary Face
            if (isIt->boundary())
            {
                this->getFluxOnBoundary(entries, timestepFlux, *isIt, cellDataI);
            }

            // add to update vector
            updateVec[wCompIdx][globalIdxI] += entries[wCompIdx];
            updateVec[nCompIdx][globalIdxI] += entries[nCompIdx];

            // for time step calculation
            sumfactorin += timestepFlux[0];
            sumfactorout += timestepFlux[1];

        }// end all intersections

        /***********     Handle source term     ***************/
        PrimaryVariables q(NAN);
        problem().source(q, *eIt);
        updateVec[wCompIdx][globalIdxI] += q[contiWEqIdx];
        updateVec[nCompIdx][globalIdxI] += q[contiNEqIdx];

        // account for porosity
        sumfactorin = std::max(sumfactorin,sumfactorout)
                        / problem().spatialParams().porosity(*eIt);

        if ( 1./sumfactorin < dt)
        {
            dt = 1./sumfactorin;
            restrictingCell= globalIdxI;
        }
    } // end grid traversal

#if HAVE_MPI
    // communicate updated values
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ElementMapper ElementMapper;
    typedef GridCommHandleSum<Dune::FieldVector<Scalar, 1>, Dune::BlockVector<Dune::FieldVector<Scalar, 1> >, ElementMapper, /*commCodim=*/0> DataHandle;
    for (int i = 0; i < updateVec.size(); i++)
    {
        DataHandle dataHandle(updateVec[i], problem_.variables().elementMapper());
        problem_.gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
    }
    dt = problem_.gridView().comm().min(dt);
#endif

    if(impet)
    {
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "
                <<dt * GET_PARAM(TypeTag, Scalar, ImpetCflFactor)<< std::endl;
        if(this->averagedFaces_ != 0)
            Dune::dinfo  << " Averageing done for " << this->averagedFaces_ << " faces. "<< std::endl;
    }
    return;
}

//! Compute flux over an irregular interface using a \a mpfa method
/** A mpfa l-method is applied to calculate fluxes near hanging nodes, using:
 * \f[
      - \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
        \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha}
    \f]
 *
 * We provide two options: Calculating the flux expressed by twice the flux
 * through the one unique interaction region on the hanging node if one
 * halfedge is stored (eg on boundaries). Or using the second interaction
 * region covering neighboring cells.
 *
 * \param fluxEntries The flux entries, mass influx from cell \f$j\f$ to \f$i\f$.
 * \param timestepFlux flow velocities for timestep estimation
 * \param intersectionIterator Iterator of the intersection between cell I and J
 * \param cellDataI The cell data for cell \f$i\f$
 */
template<class TypeTag>
void FVTransport2P2CAdaptive<TypeTag>::getMpfaFlux(Dune::FieldVector<Scalar, 2>& fluxEntries, Dune::FieldVector<Scalar, 2>& timestepFlux,
                                        const IntersectionIterator& intersectionIterator, CellData& cellDataI)
{
    fluxEntries = 0.;
    timestepFlux = 0.;
    // cell information
    ElementPointer elementI= intersectionIterator->inside();
    int globalIdxI = problem().variables().index(*elementI);

    // get position
    const GlobalPosition globalPos = elementI->geometry().center();
    const GlobalPosition& gravity_ = problem().gravity();
    // cell volume, assume linear map here
    Scalar volume = elementI->geometry().volume();

    // get values of cell I
    Scalar pressI = problem().pressureModel().pressure(globalIdxI);
    Scalar pcI = cellDataI.capillaryPressure();
    DimMatrix K_I(problem().spatialParams().intrinsicPermeability(*elementI));

    PhaseVector SmobI(0.);
    SmobI[wPhaseIdx] = std::max((cellDataI.saturation(wPhaseIdx)
                            - problem().spatialParams().materialLawParams(*elementI).Swr())
                            , 1e-2);
    SmobI[nPhaseIdx] = std::max((cellDataI.saturation(nPhaseIdx)
                                - problem().spatialParams().materialLawParams(*elementI).Snr())
                            , 1e-2);

    Scalar densityWI (0.), densityNWI(0.);
    densityWI= cellDataI.density(wPhaseIdx);
    densityNWI = cellDataI.density(nPhaseIdx);

    PhaseVector potential(0.);

    // access neighbor
    ElementPointer neighborPointer = intersectionIterator->outside();
    int globalIdxJ = problem().variables().index(*neighborPointer);
    CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

    // neighbor cell center in global coordinates
    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosNeighbor - globalPos;
    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    // phase densities in neighbor
    Scalar densityWJ (0.), densityNWJ(0.);
    densityWJ = cellDataJ.density(wPhaseIdx);
    densityNWJ = cellDataJ.density(nPhaseIdx);

    // average phase densities with central weighting
    double densityW_mean = (densityWI + densityWJ) * 0.5;
    double densityNW_mean = (densityNWI + densityNWJ) * 0.5;

    double pressJ = problem().pressureModel().pressure(globalIdxJ);
    Scalar pcJ = cellDataJ.capillaryPressure();

    // determine potentials for upwind
        /** get geometrical Info, transmissibility matrix */
        GlobalPosition globalPos3(0.);
        int globalIdx3=-1;
        TransmissivityMatrix T(0.);
        IntersectionIterator additionalIsIt = intersectionIterator;
        TransmissivityMatrix additionalT(0.);

        int halfedgesStored
            = problem().variables().getMpfaData(*intersectionIterator, additionalIsIt, T, additionalT, globalPos3, globalIdx3);
        if (halfedgesStored == 0)
            halfedgesStored = problem().pressureModel().computeTransmissibilities(intersectionIterator,additionalIsIt, T,additionalT, globalPos3, globalIdx3 );

        // acess cell 3 and prepare mpfa
        Scalar press3 = problem().pressureModel().pressure(globalIdx3);
        CellData& cellData3 = problem().variables().cellData(globalIdx3);
        Scalar pc3 = cellData3.capillaryPressure();
        Scalar temp1 = globalPos * gravity_;
        Scalar temp2 = globalPosNeighbor * gravity_;
        Scalar temp3 = globalPos3 * gravity_;
        if(pressureType==pw)
        {
            potential[wPhaseIdx] += (pressI-temp1*densityW_mean) * T[2]
                            +(pressJ-temp2*densityW_mean) * T[0]
                            +(press3- temp3*densityW_mean) * T[1];
            potential[nPhaseIdx] += (pressI+pcI-temp1*densityNW_mean) * T[2]
                            +(pressJ+pcJ-temp2*densityNW_mean) * T[0]
                            +(press3+pc3- temp3*densityNW_mean) * T[1];
            // second half edge, if there is one
            if(halfedgesStored == 2)
            {
                int AdditionalIdx = problem().variables().index(*(additionalIsIt->outside()));
                CellData& cellDataAdditional = problem().variables().cellData(AdditionalIdx);
                potential[wPhaseIdx] += (pressI-temp1*densityW_mean) * additionalT[2]
                                +(pressJ-temp2*densityW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                    -(additionalIsIt->outside()->geometry().center()*gravity_*densityW_mean)
                                  ) * additionalT[1];
                potential[nPhaseIdx] += (pressI+pcI-temp1*densityNW_mean) * additionalT[2]
                                +(pressJ+pcJ-temp2*densityNW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                    + cellDataAdditional.capillaryPressure()
                                    -(additionalIsIt->outside()->geometry().center()*gravity_*densityNW_mean)
                                  ) * additionalT[1];
            }
        }
        else if(pressureType==pn)
        {
            potential[wPhaseIdx] += (pressI-pcI-temp1*densityW_mean) * T[2]
                          + (pressJ-pcJ-temp2*densityW_mean) * T[0]
                          + (press3-pc3- temp3*densityW_mean) * T[1];
            potential[nPhaseIdx] += (pressI-temp1*densityNW_mean) * T[2]
                          + (pressJ-temp2*densityNW_mean) * T[0]
                          + (press3-temp3*densityNW_mean) * T[1];
            // second half edge, if there is one
            if(halfedgesStored == 2)
            {
                int AdditionalIdx = problem().variables().index(*(additionalIsIt->outside()));
                CellData& cellDataAdditional = problem().variables().cellData(AdditionalIdx);

                potential[wPhaseIdx] += (pressI-pcI-temp1*densityW_mean) * additionalT[2]
                                +(pressJ-pcJ-temp2*densityW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                   - cellDataAdditional.capillaryPressure()
                                   -(additionalIsIt->outside()->geometry().center()*gravity_*densityW_mean)
                                  ) * additionalT[1];
                potential[nPhaseIdx] += (pressI-temp1*densityNW_mean) * additionalT[2]
                                +(pressJ-temp2*densityNW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                    -(additionalIsIt->outside()->geometry().center()*gravity_*densityNW_mean))
                                    * additionalT[1];
            }
        }
    //end of mpfa specific stuff

    // determine upwinding direction, perform upwinding if possible
    Dune::FieldVector<bool, NumPhases> doUpwinding(true);
    PhaseVector lambda(0.);
    for(int phaseIdx = 0; phaseIdx < NumPhases; phaseIdx++)
    {
        int contiEqIdx = 0;
        if(phaseIdx == wPhaseIdx)
            contiEqIdx = contiWEqIdx;
        else
            contiEqIdx = contiNEqIdx;

        if(!this->impet_ or !this->restrictFluxInTransport_) // perform a strict uwpind scheme
        {
            if(potential[phaseIdx] > 0.)
            {
                lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                cellDataI.setUpwindCell(intersectionIterator->indexInInside(), contiEqIdx, true);
            }
            else if(potential[phaseIdx] < 0.)
            {
                lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
                cellDataI.setUpwindCell(intersectionIterator->indexInInside(), contiEqIdx, false);
            }
            else
            {
                doUpwinding[phaseIdx] = false;
                cellDataI.setUpwindCell(intersectionIterator->indexInInside(), contiEqIdx, false);
                cellDataJ.setUpwindCell(intersectionIterator->indexInOutside(), contiEqIdx, false);
            }
        }
        else // Transport after PE with check on flow direction
        {
            bool cellIwasUpwindCell;
            //get the information from smaller (higher level) cell, as its IS is unique
            if(elementI->level()>neighborPointer->level())
                cellIwasUpwindCell = cellDataI.isUpwindCell(intersectionIterator->indexInInside(), contiEqIdx);
            else // reverse neighbors information gathered
                cellIwasUpwindCell = !cellDataJ.isUpwindCell(intersectionIterator->indexInOutside(), contiEqIdx);

            if (potential[phaseIdx] > 0. && cellIwasUpwindCell)
                lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
            else if (potential[phaseIdx] < 0. && !cellIwasUpwindCell)
                lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
            // potential direction does not coincide with that of P.E.
            else if(this->restrictFluxInTransport_ == 2)   // perform central averageing for all direction changes
                doUpwinding[phaseIdx] = false;
            else    // i.e. restrictFluxInTransport == 1
            {
               //check if harmonic weithing is necessary
                if (!cellIwasUpwindCell && cellDataJ.mobility(phaseIdx) != 0.) // check if outflow induce neglected (i.e. mob=0) phase flux
                    lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                else if (cellIwasUpwindCell && cellDataI.mobility(phaseIdx) != 0.) // check if inflow induce neglected phase flux
                    lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
                else
                    doUpwinding[phaseIdx] = false;
            }

            // do not perform upwinding if so desired
            if(!doUpwinding[phaseIdx])
            {
                //a) no flux if there wont be any flux regardless how to average/upwind
                if(cellDataI.mobility(phaseIdx)+cellDataJ.mobility(phaseIdx)==0.)
                {
                    potential[phaseIdx] = 0;
                    continue;
                }

                //b) perform harmonic averageing
                fluxEntries[wCompIdx] -= potential[phaseIdx] / volume
                        * harmonicMean(cellDataI.massFraction(phaseIdx, wCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                                cellDataJ.massFraction(phaseIdx, wCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
                fluxEntries[nCompIdx] -= potential[phaseIdx] / volume
                        * harmonicMean(cellDataI.massFraction(phaseIdx, nCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                                cellDataJ.massFraction(phaseIdx, nCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
                // c) timestep control
                // for timestep control : influx
                timestepFlux[0] += std::max(0.,
                        - potential[phaseIdx] / volume
                          * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx)));
                // outflux
                timestepFlux[1] += std::max(0.,
                        potential[phaseIdx]  / volume
                        * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx))/SmobI[phaseIdx]);

                //d) output (only for one side)
                this->averagedFaces_++;
                #if DUNE_MINIMAL_DEBUG_LEVEL < 3
                // verbose (only for one side)
                if(globalIdxI > globalIdxJ)
                    Dune::dinfo << "harmonicMean flux of phase" << phaseIdx <<" used from cell" << globalIdxI<< " into " << globalIdxJ
                    << " ; TE upwind I = "<< cellDataI.isUpwindCell(intersectionIterator->indexInInside(), contiEqIdx) << " but pot = "<< potential[phaseIdx] <<  " \n";
                #endif

                //e) stop further standard calculations
                potential[phaseIdx] = 0;
            }
        }
    }

    // calculate and standardized velocity
    double velocityJIw = std::max((-lambda[wPhaseIdx] * potential[wPhaseIdx]) / volume, 0.0);
    double velocityIJw = std::max(( lambda[wPhaseIdx] * potential[wPhaseIdx]) / volume, 0.0);
    double velocityJIn = std::max((-lambda[nPhaseIdx] * potential[nPhaseIdx]) / volume, 0.0);
    double velocityIJn = std::max(( lambda[nPhaseIdx] * potential[nPhaseIdx]) / volume, 0.0);

    // for timestep control
    timestepFlux[0] += velocityJIw + velocityJIn;

    double foutw = velocityIJw/SmobI[wPhaseIdx];
    double foutn = velocityIJn/SmobI[nPhaseIdx];
    if (std::isnan(foutw) || std::isinf(foutw) || foutw < 0) foutw = 0;
    if (std::isnan(foutn) || std::isinf(foutn) || foutn < 0) foutn = 0;
    timestepFlux[1] += foutw + foutn;

    fluxEntries[wCompIdx] +=
          velocityJIw * cellDataJ.massFraction(wPhaseIdx, wCompIdx) * densityWJ
        - velocityIJw * cellDataI.massFraction(wPhaseIdx, wCompIdx) * densityWI
        + velocityJIn * cellDataJ.massFraction(nPhaseIdx, wCompIdx) * densityNWJ
        - velocityIJn * cellDataI.massFraction(nPhaseIdx, wCompIdx) * densityNWI;
    fluxEntries[nCompIdx] +=
          velocityJIw * cellDataJ.massFraction(wPhaseIdx, nCompIdx) * densityWJ
        - velocityIJw * cellDataI.massFraction(wPhaseIdx, nCompIdx) * densityWI
        + velocityJIn * cellDataJ.massFraction(nPhaseIdx, nCompIdx) * densityNWJ
        - velocityIJn * cellDataI.massFraction(nPhaseIdx, nCompIdx) * densityNWI;
    return;
}

}
#endif
