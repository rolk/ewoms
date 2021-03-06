// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
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
/*!
 * \file
 * \copydoc Ewoms::MimeticGroundwaterEquationLocalStiffness
 */
#ifndef EWOMS_MIMETICGROUNDWATER_HH
#define EWOMS_MIMETICGROUNDWATER_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

#include <ewoms/decoupled/common/boundaryconditions.hh>
#include "localstiffness.hh"
#include <ewoms/decoupled/2p/diffusion/diffusionproperties2p.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include<iomanip>

namespace Ewoms
{
/*!
 *  \ingroup MimeticPressure2p
 */
/*!
 * \brief compute local stiffness matrix for conforming finite elements for diffusion equation
 *
 */

//! A class for computing local stiffness matrices
/*! A class for computing local stiffness matrix for the diffusion equation
 * \f{eqnarray*}{
 * \text{div} \boldsymbol v_{total} &=& q \\
 * \boldsymbol v_{total} &=& - \lambda \boldsymbol K \text{grad} p
 * \f}
 * where \f$ p = p_D \f$ on \f$ \Gamma_{Dirichlet} \f$, and \f$ \boldsymbol v_{total} \cdot  \boldsymbol n  = q_N \f$ on \f$ \Gamma_{Neumann} \f$.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class MimeticGroundwaterEquationLocalStiffness: public LocalStiffness<TypeTag, 1>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView)GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // grid types
    enum
    {
        dim = GridView::dimension
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        pressEqIdx = Indices::pressEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum
    {   m=1};
    enum
    {   size=2*dim};

    //! Constructor
    MimeticGroundwaterEquationLocalStiffness (Problem& problem,
            bool levelBoundaryAsDirichlet, const GridView& gridView,
            bool procBoundaryAsDirichlet=true)
    : problem_(problem), gridView_(gridView)
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
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, paramCache, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, paramCache, nPhaseIdx);
    }

    //! assemble local stiffness matrix for given element and order
    /*! On exit the following things have been done:
     - The stiffness matrix for the given entity and polynomial degree has been assembled and is
     accessible with the mat() method.
     - The boundary conditions have been evaluated and are accessible with the bc() method
     - The right hand side has been assembled. It contains either the value of the essential boundary
     condition or the assembled source term and neumann boundary condition. It is accessible via the rhs() method.
     \param[in]  element a codim 0 entity reference
     \param[in]  k order of CR basis (only k = 1 is implemented)
     */
    void assemble (const Element& element, int k=1)
    {
        unsigned int numFaces = element.template count<1>();
        this->setcurrentsize(numFaces);

        // clear assemble data
        for (unsigned i=0; i<numFaces; i++)
        {
            this->b[i] = 0;
            this->bctype[i].reset();
            for (unsigned j=0; j<numFaces; j++)
            this->A[i][j] = 0;
        }

        assembleV(element,k);
        assembleBC(element,k);
    }

    // TODO/FIXME: this is only valid for linear problems where
    // the local stiffness matrix is independend of the current
    // solution. We need to implement this properly, but this
    // should at least make the thing compile...
    typedef Dune::FieldVector<Scalar, m> VBlockType;
    void assemble(const Element &cell, const Dune::BlockVector<VBlockType>& localSolution, int orderOfShapeFns = 1)
    {
        assemble(cell, orderOfShapeFns);
    }

    //! assemble only boundary conditions for given element
    /*! On exit the following things have been done:
     - The boundary conditions have been evaluated and are accessible with the bc() method
     - The right hand side contains either the value of the essential boundary
     condition or the assembled neumann boundary condition. It is accessible via the rhs() method.
     \param[in]  element a codim 0 entity reference
     \param[in]  k order of CR basis
     */
    void assembleBoundaryCondition (const Element& element, int k=1)
    {
        unsigned int numFaces = element.template count<1>();
        this->setcurrentsize(numFaces);

        // clear assemble data
        for (unsigned i=0; i<numFaces; i++)
        {
            this->b[i] = 0;
            this->bctype[i].reset();
        }

        assembleBC(element,k);
    }

    void assembleElementMatrices(const Element& element, Dune::FieldVector<Scalar,2*dim>& faceVol,
            Dune::FieldMatrix<Scalar,2*dim,2*dim>& W, Dune::FieldVector<Scalar,2*dim>& c,
            Dune::FieldMatrix<Scalar,2*dim,2*dim>& Pi, Scalar& dinv, Dune::FieldVector<Scalar,2*dim>& F, Scalar& qmean)
    {
        unsigned int numFaces = element.template count<1>();
        this->setcurrentsize(numFaces);

        // get global coordinate of cell center
        const Dune::FieldVector<Scalar,dim>& centerGlobal = element.geometry().center();

        int globalIdx = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        // eval diffusion tensor, ASSUMING to be constant over each cell
        Dune::FieldMatrix<Scalar,dim,dim> K(0);
        K = problem_.spatialParams().intrinsicPermeability(element);

        K *= (cellData.mobility(wPhaseIdx) + cellData.mobility(nPhaseIdx));

        // cell volume
        Scalar volume = element.geometry().volume();

        // build the matrices R and ~N
        Dune::FieldMatrix<Scalar,2*dim,dim> R(0), N(0);

        //       std::cout << "element " << elemId << ": center " << centerGlobal << std::endl;;

        typedef typename GridView::IntersectionIterator IntersectionIterator;

        int i = -1;
        IntersectionIterator endit = gridView_.iend(element);
        for (IntersectionIterator it = gridView_.ibegin(element); it!=endit; ++it)
        {
            // local number of facet
            i = it->indexInInside();

            const Dune::FieldVector<Scalar,dim>& faceGlobal = it->geometry().center();
            faceVol[i] = it->geometry().volume();

            // get normal vector
            const Dune::FieldVector<Scalar,dim>& unitOuterNormal = it->centerUnitOuterNormal();

            N[i] = unitOuterNormal;

            for (int k = 0; k < dim; k++)
            // move origin to the center of gravity
            R[i][k] = faceVol[i]*(faceGlobal[k] - centerGlobal[k]);
        }

        //      std::cout << "N =\dim" << N;
        //      std::cout << "R =\dim" << R;

        // proceed along the lines of Algorithm 1 from
        // Brezzi/Lipnikov/Simonicini M3AS 2005
        // (1) orthonormalize columns of the matrix R
        Scalar norm = R[0][0]*R[0][0];
        for (unsigned i = 1; i < numFaces; i++)
        norm += R[i][0]*R[i][0];
        norm = std::sqrt(norm);
        for (unsigned i = 0; i < numFaces; i++)
        R[i][0] /= norm;
        Scalar weight = R[0][1]*R[0][0];
        for (unsigned i = 1; i < numFaces; i++)
        weight += R[i][1]*R[i][0];
        for (unsigned i = 0; i < numFaces; i++)
        R[i][1] -= weight*R[i][0];
        norm = R[0][1]*R[0][1];
        for (unsigned i = 1; i < numFaces; i++)
        norm += R[i][1]*R[i][1];
        norm = std::sqrt(norm);
        for (unsigned i = 0; i < numFaces; i++)
        R[i][1] /= norm;
        if (dim == 3)
        {
            Scalar weight1 = R[0][2]*R[0][0];
            Scalar weight2 = R[0][2]*R[0][1];
            for (unsigned i = 1; i < numFaces; i++)
            {
                weight1 += R[i][2]*R[i][0];
                weight2 += R[i][2]*R[i][1];
            }
            for (unsigned i = 0; i < numFaces; i++)
            R[i][1] -= weight1*R[i][0] + weight2*R[i][1];
            norm = R[0][2]*R[0][2];
            for (unsigned i = 1; i < numFaces; i++)
            norm += R[i][2]*R[i][2];
            norm = std::sqrt(norm);
            for (unsigned i = 0; i < numFaces; i++)
            R[i][2] /= norm;
        }
        //      std::cout << "~R =\dim" << R;

        // (2) Build the matrix ~D
        Dune::FieldMatrix<Scalar,2*dim,2*dim> D(0);
        for (unsigned s = 0; s < numFaces; s++)
        {
            Dune::FieldVector<Scalar,2*dim> es(0);
            es[s] = 1;
            for (unsigned k = 0; k < numFaces; k++)
            {
                D[k][s] = es[k];
                for (unsigned i = 0; i < dim; i++)
                {
                    D[k][s] -= R[s][i]*R[k][i];
                }
            }
        }

        Scalar traceK = K[0][0];
        for (unsigned i = 1; i < dim; i++)
        traceK += K[i][i];
        D *= 2.0*traceK/volume;
        //      std::cout << "u~D =\dim" << D;

        // (3) Build the matrix W = Minv
        Dune::FieldMatrix<Scalar,2*dim,dim> NK(N);
        NK.rightmultiply (K);
        for (unsigned i = 0; i < numFaces; i++)
        {
            for (unsigned j = 0; j < numFaces; j++)
            {
                W[i][j] = NK[i][0]*N[j][0];
                for (unsigned k = 1; k < dim; k++)
                W[i][j] += NK[i][k]*N[j][k];
            }
        }

        W /= volume;
        W += D;
        //       std::cout << "W = \dim" << W;
        //       std::cout << D[2][2] << ", " << D[2][0] << ", " << D[2][3] << ", " << D[2][1] << std::endl;
        //       std::cout << D[0][2] << ", " << D[0][0] << ", " << D[0][3] << ", " << D[0][1] << std::endl;
        //       std::cout << D[3][2] << ", " << D[3][0] << ", " << D[3][3] << ", " << D[3][1] << std::endl;
        //       std::cout << D[1][2] << ", " << D[1][0] << ", " << D[1][3] << ", " << D[1][1] << std::endl;

        // Now the notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
        // The matrix W developed so far corresponds to one element-associated
        // block of the matrix B^{-1} there.

        // Corresponding to the element under consideration,
        // calculate the part of the matrix C coupling velocities and element pressures.
        // This is just a row vector of size numFaces.
        // scale with volume
        for (unsigned i = 0; i < numFaces; i++)
        c[i] = faceVol[i];

        // Set up the element part of the matrix \Pi coupling velocities
        // and pressure-traces. This is a diagonal matrix with entries given by faceVol.
        for (unsigned i = 0; i < numFaces; i++)
        Pi[i][i] = faceVol[i];

        // Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
        Dune::FieldVector<Scalar,2*dim> Wc(0);
        W.umv(c, Wc);
        dinv = 1.0/(c*Wc);

        // Calculate the element part of the matrix F = Pi W c^T which is a column vector.
        F = 0;
        Pi.umv(Wc, F);
        //      std::cout << "Pi = \dim" << Pi << "c = " << c << ", F = " << F << std::endl;

        PrimaryVariables sourceVec(0.0);
        problem_.source(sourceVec, element);
        qmean = volume*(sourceVec[wPhaseIdx]/density_[wPhaseIdx] + sourceVec[nPhaseIdx]/density_[nPhaseIdx]);
    }

private:
    void assembleV (const Element& element, int k=1)
    {
        unsigned int numFaces = element.template count<1>();
        this->setcurrentsize(numFaces);

        // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
        // The matrix W developed here corresponds to one element-associated
        // block of the matrix B^{-1} there.
        Dune::FieldVector<Scalar,2*dim> faceVol(0);
        Dune::FieldMatrix<Scalar,2*dim,2*dim> W(0);
        Dune::FieldVector<Scalar,2*dim> c(0);
        Dune::FieldMatrix<Scalar,2*dim,2*dim> Pi(0);
        Dune::FieldVector<Scalar,2*dim> F(0);
        Scalar dinv;
        Scalar qmean;
        this->assembleElementMatrices(element, faceVol, W, c, Pi, dinv, F, qmean);

        // Calculate the element part of the matrix Pi W Pi^T.
        Dune::FieldMatrix<Scalar,2*dim,2*dim> PiWPiT(W);
        PiWPiT.rightmultiply(Pi);
        PiWPiT.leftmultiply(Pi);

        // Calculate the element part of the matrix F D^{-1} F^T.
        Dune::FieldMatrix<Scalar,2*dim,2*dim> FDinvFT(0);
        for (unsigned i = 0; i < numFaces; i++)
            for (unsigned j = 0; j < numFaces; j++)
                FDinvFT[i][j] = dinv*F[i]*F[j];

        // Calculate the element part of the matrix S = Pi W Pi^T - F D^{-1} F^T.
        for (unsigned i = 0; i < numFaces; i++)
            for (unsigned j = 0; j < numFaces; j++)
                this->A[i][j] = PiWPiT[i][j] - FDinvFT[i][j];

        // Calculate the source term F D^{-1} f
        // NOT WORKING AT THE MOMENT
        Scalar factor = dinv*qmean;
        for (unsigned i = 0; i < numFaces; i++)
            this->b[i] = F[i]*factor;

        //        std::cout << "faceVol = " << faceVol << std::endl << "W = " << std::endl << W << std::endl
        //              << "c = " << c << std::endl << "Pi = " << std::endl << Pi << std::endl
        //              << "dinv = " << dinv << std::endl << "F = " << F << std::endl;
        //        std::cout << "dinvF = " << dinvF << ", q = " << qmean
        //             << ", b = " << this->b[0] << ", " << this->b[1] << ", " << this->b[2] << ", " << this->b[3] << std::endl;
    }

    void assembleBC (const Element& element, int k=1)
    {
        // evaluate boundary conditions via intersection iterator
        typedef typename GridView::IntersectionIterator IntersectionIterator;

        IntersectionIterator endIsIt = gridView_.iend(element);
        for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != endIsIt; ++isIt)
        {
            int faceIndex = isIt->indexInInside();
            if (isIt->boundary())
            {
                problem_.boundaryTypes(this->bctype[faceIndex], *isIt);
                PrimaryVariables boundValues(0.0);

                if (this->bctype[faceIndex].isNeumann(pressEqIdx))
                {
                    problem_.neumann(boundValues, *isIt);
                    Scalar J = (boundValues[wPhaseIdx]/density_[wPhaseIdx] + boundValues[nPhaseIdx]/density_[nPhaseIdx]);
                    this->b[faceIndex] -= J * isIt->geometry().volume();
                }
                else if (this->bctype[faceIndex].isDirichlet(pressEqIdx))
                {
                    problem_.dirichlet(boundValues, *isIt);
                    this->b[faceIndex] = boundValues[pressureIdx];
                }
            }
        }
    }

private:
    Problem& problem_;
    const GridView& gridView_;

    Scalar density_[numPhases];
};

/** \} */
}
#endif
