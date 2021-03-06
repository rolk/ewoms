// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Benjamin Faigle                              *
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
 * \copydoc Ewoms::Test2P2CSpatialParams
 */
#ifndef TEST_2P2C_SPATIAL_PARAMS_HH
#define TEST_2P2C_SPATIAL_PARAMS_HH

#include <ewoms/decoupled/2p2c/2p2cproperties.hh>
#include <ewoms/decoupled/spatialparams/fvspatialparams.hh>
#include <ewoms/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dune/common/fmatrix.hh>

namespace Ewoms
{
//forward declaration
template<class TypeTag>
class Test2P2CSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test2P2CSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Test2P2CSpatialParams, SpatialParams, Ewoms::Test2P2CSpatialParams<TypeTag>);

// Set the material law
SET_PROP(Test2P2CSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    //    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
    typedef LinearMaterial<Scalar>         RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the sequential 2p2c test
 */
template<class TypeTag>
class Test2P2CSpatialParams : public FVSpatialParams<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid)     Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
        return constPermeability_;
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }


    Test2P2CSpatialParams(const GridView& gridView) : FVSpatialParams<TypeTag>(gridView),
            constPermeability_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0);
        materialLawParams_.setSnr(0);

//        // parameters for the Brooks-Corey Law
//        // entry pressures
//        materialLawParams_.setPe(10000);
//
//        // Brooks-Corey shape parameters
//        materialLawParams_.setLambda(2);

        // parameters for the linear
        // entry pressures function

        materialLawParams_.setEntryPC(0);
        materialLawParams_.setMaxPC(10000);

        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-12;
        }
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix constPermeability_;

};

} // end namespace
#endif
