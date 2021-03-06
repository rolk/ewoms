// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
 *   Copyright (C) 2010-2012 by Benjamin Faigle                              *
 *   Copyright (C) 2011 by Michael Sinsbeck                                  *
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
#ifndef EWOMS_DECOUPLED_PROPERTIES_HH
#define EWOMS_DECOUPLED_PROPERTIES_HH

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/basicproperties.hh>
#include <ewoms/decoupled/common/gridadaptproperties.hh>

#include <dune/common/fvector.hh>

/*!
 * \ingroup Sequential
 */
/*!
 * \file
 * \brief Base file for properties related to sequential (decoupled) models
 */
namespace Ewoms
{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! Create a type tag for all decoupled models
NEW_TYPE_TAG(DecoupledModel, INHERITS_FROM(NumericModel, GridAdaptTypeTag));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

//! Property tag for types associated with the solution of the PDE.
//! This means vectors of primary variables, solution functions on the
//! grid, and elements, and shape functions.
NEW_PROP_TAG( SolutionTypes);
NEW_PROP_TAG( PrimaryVariables);
NEW_PROP_TAG( Indices);

NEW_PROP_TAG( Grid); //!< The type of the DUNE grid
NEW_PROP_TAG( GridView); //!< The type of the grid view

NEW_PROP_TAG( Problem); //!< The type of the problem
NEW_PROP_TAG( Model); //!< The type of the discretizations
NEW_PROP_TAG( PressureModel ); //!< The type of the discretization of a pressure model
NEW_PROP_TAG( TransportModel ); //!< The type of the discretization of a transport model
NEW_PROP_TAG( Velocity ); //!< The type velocity reconstruction
NEW_PROP_TAG( NumEq ); //!< Number of equations in the system of PDEs
NEW_PROP_TAG( NumPhases); //!< Number of phases in the system
NEW_PROP_TAG( NumComponents); //!< Number of components in the system
NEW_PROP_TAG( Variables); //!< The type of the container of global variables
NEW_PROP_TAG( CellData );//!< Defines data object to be stored
NEW_PROP_TAG( TimeManager );  //!< Manages the simulation time
NEW_PROP_TAG( BoundaryTypes ); //!< Stores the boundary types of a single degree of freedom
NEW_PROP_TAG( MaxIntersections ); //!< Gives maximum number of intersections of an element and neighboring elements
NEW_PROP_TAG( VtkOutputLevel); //! VtkOutputLevel is equal to zero only primary variables are written, all available quantities are written if it's larger than 0.
}
}

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bvector.hh>

#include <ewoms/common/timemanager.hh>
#include <ewoms/decoupled/common/boundarytypes.hh>
#include <ewoms/decoupled/common/boundaryconditions.hh>

namespace Ewoms
{

template<class TypeTag>
class GridAdaptInitializationIndicatorDefault;

template<class TypeTag>
class VariableClass;

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_PROP(DecoupledModel, GridView)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

public:
    typedef typename Grid::LeafGridView type;
};

//! Default number of intersections for quadrilaterals
SET_PROP(DecoupledModel, MaxIntersections)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum
    {
        dim = GridView::dimension
    };
public:
    static const int value = 2*dim;
};

/*!
 * \brief Specifies the types which are assoicated with a solution.
 *
 * This means shape functions, solution vectors, etc.
 */
SET_PROP(DecoupledModel, SolutionTypes)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, Variables) Variables;

    enum
    {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        maxIntersections = GET_PROP_VALUE(TypeTag, MaxIntersections)
    };

    template<int dim>
    struct VertexLayout
    {
        bool contains (Dune::GeometryType gt) const
        {   return gt.dim() == 0;}
    };

    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt) const
        {   return gt.dim() == dim;}
    };

public:
    /*!
     * \brief Mapper for the grid view's vertices.
     */
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, VertexLayout> VertexMapper;

    /*!
     * \brief Mapper for the grid view's elements.
     */
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, ElementLayout> ElementMapper;

    /*!
     * \brief The type of a solution at a fixed time.
     *
     * This defines the primary and secondary variable vectors at each degree of freedom.
     */
    typedef Dune::FieldVector<Scalar, numEq> PrimaryVariables;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarSolution;//!<type for vector of scalars
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numComponents> ComponentProperty;//!<type for vector of phase properties
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numPhases> PhaseProperty;//!<type for vector of phase properties
    typedef Dune::FieldVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> >, numPhases> FluidProperty;//!<type for vector of fluid properties: Vector[element][phase]
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, numPhases>, maxIntersections > > PhasePropertyElemFace;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef Dune::BlockVector<Dune::FieldVector<Dune::FieldVector<Scalar, dim>, maxIntersections > > DimVecElemFace;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars
};

SET_TYPE_PROP(DecoupledModel,  Variables, VariableClass<TypeTag>);

SET_TYPE_PROP(DecoupledModel,  PrimaryVariables, typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables);

//! Set the default type for the time manager
SET_TYPE_PROP(DecoupledModel, TimeManager, Ewoms::TimeManager<TypeTag>);

//! By default, write out everything
SET_INT_PROP(DecoupledModel, VtkOutputLevel, 2);

SET_SCALAR_PROP(DecoupledModel, InitialTimeStepSize, 0.0);

/*!
 * \brief Boundary types at a single degree of freedom.
 */
SET_PROP(DecoupledModel, BoundaryTypes)
{ private:
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
public:
    typedef Ewoms::BoundaryTypes<numEq>  type;
};

//Set default class for adaptation initialization indicator
SET_TYPE_PROP(GridAdaptTypeTag,  GridAdaptInitializationIndicator, GridAdaptInitializationIndicatorDefault<TypeTag>);
}
}

#include "gridadaptinitializationindicatordefault.hh"

#endif
