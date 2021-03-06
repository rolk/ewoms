// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Benjamin Faigle                              *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
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
#ifndef EWOMS_FVTRANSPORT2P2C_MULTIPHYSICS_HH
#define EWOMS_FVTRANSPORT2P2C_MULTIPHYSICS_HH

#include <ewoms/decoupled/2p2c/fvtransport2p2c.hh>
#include <ewoms/parallel/gridcommhandles.hh>

#include <dune/common/fvector.hh>

/*!
 * \file
 * \copydoc Ewoms::FVTransport2P2CMultiPhysics
 */

namespace Ewoms
{
//! Compositional Transport Step in a Finite Volume discretization
/*!
 * \ingroup multiphysics
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *
 * The model domain is automatically divided into a single-phase and a two-phase domain. As the flux computation is relatively cheap,
 * the same method is used for the real transport step independently of the subdomain.
 * The pressure equation does not need any derivatives in simple
 * subdomains, therefore in the transport estimate step inter-cell fluxes in the simple subdomain are omitted.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2CMultiPhysics : public FVTransport2P2C<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;


    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, 2> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    //! Acess function for the current problem
    Problem& problem()
    {return this->problem_;}
public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false);

    //! Constructs a FVTransport2P2CMultiPhysics object
    /*!
     * \param problem a problem class object
     */
    FVTransport2P2CMultiPhysics(Problem& problem) : FVTransport2P2C<TypeTag>(problem)
    {}

    virtual ~FVTransport2P2CMultiPhysics()
    {     }
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
 *  employing a CFL condition.
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param[out] dt Time step size \f$\mathrm{[s]}\f$
 *  \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2CMultiPhysics<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet)
{
    // initialize dt very large
    dt = 1E100;
    // store if we do update Estimate for flux functions
    this->impet_ = impet;
    this->averagedFaces_ = 0.;

    // resize update vector and set to zero
    updateVec.resize(GET_PROP_VALUE(TypeTag, NumComponents));
    updateVec[wCompIdx].resize(problem().gridView().size(0));
    updateVec[nCompIdx].resize(problem().gridView().size(0));
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    PhaseVector entries(0.), timestepFlux(0.);
    // compute update vector
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get cell infos
        int globalIdxI = problem().variables().index(*eIt);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        if(impet or cellDataI.subdomain()==2)   // estimate only necessary in subdomain
        {
            // some variables for time step calculation
            double sumfactorin = 0;
            double sumfactorout = 0;

            // run through all intersections with neighbors and boundary
            IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            {

                /****** interior face   *****************/
                if (isIt->neighbor())
                    this->getFlux(entries, timestepFlux, *isIt, cellDataI);

                /******  Boundary Face   *****************/
                if (isIt->boundary())
                    this->getFluxOnBoundary(entries, timestepFlux, *isIt, cellDataI);

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
            updateVec[wCompIdx][globalIdxI] += q[Indices::contiWEqIdx];
            updateVec[nCompIdx][globalIdxI] += q[Indices::contiNEqIdx];

            // account for porosity in fluxes for time-step
            sumfactorin = std::max(sumfactorin,sumfactorout)
                            / problem().spatialParams().porosity(*eIt);

            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= globalIdxI;
            }
        }
    } // end grid traversal

#if HAVE_MPI
    // communicate updated values
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ElementMapper ElementMapper;
    typedef GridCommHandleSum<Dune::FieldVector<Scalar, 1>, Dune::BlockVector<Dune::FieldVector<Scalar, 1> >, ElementMapper, /*commCodim=*/0> DataHandle;
    for (int i = 0; i < updateVec.size(); i++)
    {
        DataHandle dataHandle(updateVec[i], problem().variables().elementMapper());
        problem().gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
    }
    dt = problem().gridView().comm().min(dt);
#endif

    if(impet)
    {
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "<<dt * GET_PARAM(TypeTag, Scalar, ImpetCflFactor)<< std::endl;
        if(this->averagedFaces_ != 0)
            Dune::dinfo  << " Averageing done for " << this->averagedFaces_ << " faces. "<< std::endl;
    }
    return;
}
}
#endif
