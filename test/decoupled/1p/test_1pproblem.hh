// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
 *                                                                            *
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
 * \copydoc Ewoms::TestProblemOneP
 */
#ifndef EWOMS_TEST_1P_PROBLEM_HH
#define EWOMS_TEST_1P_PROBLEM_HH

#include "test_1pspatialparams.hh"

#include <ewoms/io/cubegridcreator.hh>
#include <ewoms/material/fluidsystems/liquidphase.hh>
#include <ewoms/material/components/unit.hh>
#include <ewoms/decoupled/1p/diffusion/fv/fvpressureproperties1p.hh>
#include <ewoms/decoupled/1p/diffusion/diffusionproblem1p.hh>
#include <ewoms/decoupled/common/fv/fvvelocity.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/common/fvector.hh>

namespace Ewoms {

template<class TypeTag>
class TestProblemOneP;

//////////
// Specify the properties
//////////
namespace Properties {
NEW_TYPE_TAG(TestProblemOneP, INHERITS_FROM(FVPressureOneP));

NEW_PROP_TAG(Delta);

// set the GridCreator property
SET_TYPE_PROP(TestProblemOneP, GridCreator, CubeGridCreator<TypeTag>);

// Set the grid type
SET_PROP(TestProblemOneP, Grid)
{
        typedef Dune::YaspGrid<2> type;
//    typedef Dune::SGrid<2, 2> type;
};

// Set the wetting phase
SET_PROP(TestProblemOneP, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Ewoms::LiquidPhase<Scalar, Ewoms::Unit<Scalar> > type;
};

// Set the spatial parameters
SET_TYPE_PROP(TestProblemOneP, SpatialParams, Ewoms::TestOnePSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TestProblemOneP, EnableGravity, false);

//Set the problem
SET_TYPE_PROP(TestProblemOneP, Problem, Ewoms::TestProblemOneP<TypeTag>);

SET_INT_PROP(TestProblemOneP, LinearSolverVerbosity, 1);

// define the properties required by the cube grid creator
SET_SCALAR_PROP(TestProblemOneP, DomainSizeX, 1.0);
SET_SCALAR_PROP(TestProblemOneP, DomainSizeY, 1.0);
SET_SCALAR_PROP(TestProblemOneP, DomainSizeZ, 0.0);

SET_INT_PROP(TestProblemOneP, CellsX, 30);
SET_INT_PROP(TestProblemOneP, CellsY, 30);
SET_INT_PROP(TestProblemOneP, CellsZ, 0);

SET_SCALAR_PROP(TestProblemOneP, Delta, 1e-6);
SET_SCALAR_PROP(TestProblemOneP, EndTime, 0);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the decoupled one-phase model.
 */
template<class TypeTag>
class TestProblemOneP: public DiffusionProblem1P<TypeTag >
{
    typedef DiffusionProblem1P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;


public:
    TestProblemOneP(TimeManager &timeManager) :
        ParentType(GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView()), velocity_(*this)
    {
        delta_ = GET_PARAM(TypeTag, Scalar, Delta);

        this->spatialParams().setDelta(delta_);
    }

    static void registerParameters()
    {
        ParentType::registerParameters();

        REGISTER_PARAM(TypeTag, Scalar, Delta, "The epsilon value used to calculate the derivative of the permeability.");
    }

    /*!
    * \name Problem parameters
    */
    // \{

    /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
    const char *name() const
    {
        return "test_1p";
    }

    bool shouldWriteRestartFile() const
    { return false; }

    void addOutputVtkFields()
    {
        velocity_.calculateVelocity();
        velocity_.addOutputVtkFields(this->resultWriter());
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}

    //! Returns the reference pressure for evaluation of constitutive relations
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
    }

    //!source term [kg/(m^3 s)]
    void source(PrimaryVariables &values, const Element& element) const
        {
        values = 0;

        values = integratedSource_(element, 4);
        }

    /*!
    * \brief Returns the type of boundary condition.
    *
    * BC can be dirichlet (pressure) or neumann (flux).
    */
    void boundaryTypes(BoundaryTypes &bcType,
            const Intersection& intersection) const
    {
        bcType.setAllDirichlet();
    }

    //! return dirichlet condition  (pressure, [Pa])
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values = exact(globalPos);
    }


    //! return neumann condition  (flux, [kg/(m^2 s)])
    void neumann(PrimaryVariables &values, const Intersection& intersection) const
        {
        values = 0;
        }

private:
    Scalar exact (const GlobalPosition& globalPos) const
    {
        double pi = 4.0*std::atan(1.0);

        return (std::sin(pi*globalPos[0])*std::sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
        {
        Dune::FieldVector<Scalar,dim> grad(0);
        double pi = 4.0*std::atan(1.0);
        grad[0] = pi*std::cos(pi*globalPos[0])*std::sin(pi*globalPos[1]);
        grad[1] = pi*std::cos(pi*globalPos[1])*std::sin(pi*globalPos[0]);

        return grad;
        }

    Scalar integratedSource_(const Element& element, int integrationPoints) const
    {
        Scalar source = 0.;
        LocalPosition localPos(0.0);
        GlobalPosition globalPos(0.0);
        Scalar halfInterval = 1.0/double(integrationPoints)/2.;
        for (int i = 1; i <= integrationPoints; i++)
        {
            for (int j = 1; j <= integrationPoints; j++)
            {
                localPos[0] = double(i)/double(integrationPoints) - halfInterval;
                localPos[1] = double(j)/double(integrationPoints) - halfInterval;
                globalPos = element.geometry().global(localPos);
                source += 1./(integrationPoints*integrationPoints) * evaluateSource_(globalPos);
            }
        }

        return source;
    }

    Scalar evaluateSource_(const GlobalPosition& globalPos) const
    {
        Scalar temp = temperatureAtPos(globalPos);
        Scalar referencePress = referencePressureAtPos(globalPos);

        Scalar pi = 4.0 * std::atan(1.0);
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        Scalar dpdx = pi * std::cos(pi * x) * std::sin(pi * y);
        Scalar dpdy = pi * std::sin(pi * x) * std::cos(pi * y);
        Scalar dppdxx = -pi * pi * std::sin(pi * x) * std::sin(pi * y);
        Scalar dppdxy = pi * pi * std::cos(pi * x) * std::cos(pi * y);
        Scalar dppdyx = dppdxy;
        Scalar dppdyy = dppdxx;
        Scalar kxx = (delta_* x*x + y*y)/(x*x + y*y);
        Scalar kxy = -(1.0 - delta_) * x * y / (x*x + y*y);
        Scalar kyy = (x*x + delta_*y*y)/(x*x + y*y);
        Scalar dkxxdx = 2 * x * y*y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkyydy = 2 * x*x * y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkxydx = (1.0 - delta_) * y * (x*x - y*y) /((x*x + y*y) * (x*x + y*y));
        Scalar dkxydy = (1.0 - delta_) * x * (y*y - x*x) /((x*x + y*y) * (x*x + y*y));

        Scalar fx = dkxxdx * dpdx + kxx * dppdxx + dkxydx * dpdy + kxy * dppdyx;
        Scalar fy = dkxydy * dpdx + kxy * dppdxy + dkyydy * dpdy + kyy * dppdyy;

        return -(fx + fy) / Fluid::viscosity(temp, referencePress) * Fluid::density(temp, referencePress);
    }

    double delta_;
    Ewoms::FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity) > velocity_;
};
} //end namespace

#endif
