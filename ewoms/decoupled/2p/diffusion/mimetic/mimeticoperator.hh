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
/*!
 * \file
 * \copydoc Ewoms::MimeticOperatorAssembler
 */
#ifndef EWOMS_MIMETICOPERATOR_HH
#define EWOMS_MIMETICOPERATOR_HH

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/geometry/type.hh>

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<ewoms/decoupled/common/boundaryconditions.hh>

#include"croperator.hh"

namespace Ewoms
{
/*!
 * \ingroup MimeticPressure2p
 * \brief Levelwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  \tparam Scalar The field type used in the elements of the stiffness matrix
  \tparam GridView The grid view of the simulation grid
*/
template<class Scalar, class GridView>
class MimeticOperatorAssembler : public CROperatorAssembler<Scalar, GridView>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum {dim=GridView::dimension};
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,2*dim> > VType;
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > PType;
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;

public:

    MimeticOperatorAssembler (const GridView& gridView)
    : CROperatorAssembler<Scalar, GridView>(gridView), elementMapper(gridView)
    {}

    template<class LocalStiffness, class Vector>
    void calculatePressure (LocalStiffness& loc, Vector& u,
                            VType& velocity, PType& pressure)
    {
        // run over all level elements
        Iterator eendit = this->gridView_.template end<0>();
        for (Iterator it = this->gridView_.template begin<0>(); it!=eendit; ++it)
        {
            unsigned int numFaces = it->template count<1>();

            int elemId = elementMapper.map(*it);

            // get local to global id map and pressure traces
            Dune::FieldVector<Scalar,2*dim> pressTrace(0);
            for (unsigned k = 0; k < numFaces; k++)
            {
                pressTrace[k] = u[this->faceMapper_.map(*it, k, 1)];
            }

            // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
            // The matrix W developed here corresponds to one element-associated
            // block of the matrix B^{-1} there.
            Dune::FieldVector<Scalar,2*dim> faceVol(0);
            Dune::FieldMatrix<Scalar,2*dim,2*dim> W(0);
            Dune::FieldVector<Scalar,2*dim> c(0);
            Dune::FieldMatrix<Scalar,2*dim,2*dim> Pi(0);
            Dune::FieldVector<Scalar,2*dim> F(0);
            Scalar dinv = 0;
            Scalar qmean = 0;
            loc.assembleElementMatrices(*it, faceVol, W, c, Pi, dinv, F, qmean);

            pressure[elemId] = dinv*(qmean + (F*pressTrace));

            Dune::FieldVector<Scalar,2*dim> v(0);
            for (int i = 0; i < 2*dim; i++)
                for (int j = 0; j < 2*dim; j++)
                    v[i] += W[i][j]*faceVol[j]*(pressure[elemId] - pressTrace[j]);

            velocity[elemId] = v;
        }
    }

private:
    ElementMapper elementMapper;
};
}
#endif
