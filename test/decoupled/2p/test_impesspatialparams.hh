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
/*!
 * \file
 * \copydoc Ewoms::TestIMPESSpatialParams
 */
#ifndef TEST_IMPES_SPATIAL_PARAMS_HH
#define TEST_IMPES_SPATIAL_PARAMS_HH

#include <ewoms/decoupled/spatialparams/fvspatialparams.hh>
#include <ewoms/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <ewoms/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Ewoms
{

//forward declaration
template<class TypeTag>
class TestIMPESSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TestIMPESSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TestIMPESSpatialParams, SpatialParams, Ewoms::TestIMPESSpatialParams<TypeTag>);

// Set the material law
SET_PROP(TestIMPESSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 *
 * \ingroup IMPETtests
 * \brief spatial parameters for the sequential 2p test
 */
template<class TypeTag>
class TestIMPESSpatialParams: public FVSpatialParams<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef FVSpatialParams<TypeTag> ParentType;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;


public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;


    Scalar intrinsicPermeability (const Element& element) const
    {
        return 1.0e-7;
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
//    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    const MaterialLawParams& materialLawParams(const Element& element) const
    {
            return materialLawParams_;
    }


    TestIMPESSpatialParams(const GridView& gridView)
    : ParentType(gridView)
    {
        // residual saturations
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.2);

//        // parameters for the Brooks-Corey Law
//        // entry pressures
        materialLawParams_.setPe(0);
//        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(2);

        // parameters for the linear
        // entry pressures function
//        materialLawParams_.setEntryPC(0);
//        materialLawParams_.setMaxPC(0);
    }

private:
    MaterialLawParams materialLawParams_;
};

} // end namespace
#endif
