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
/*!
 * \file
 * \copydoc Ewoms::Adaptive2p2c
 */
#ifndef EWOMS_TEST_ADAPTIVE_2P2C_PROBLEM_HH
#define EWOMS_TEST_ADAPTIVE_2P2C_PROBLEM_HH

#include <dune/grid/alugrid.hh>
// #include <dune/grid/uggrid.hh>

#include <ewoms/io/cubegridcreator.hh>

#include <ewoms/common/math.hh>
#include <ewoms/decoupled/2p2c/2p2cadaptiveproperties.hh>
#include <ewoms/decoupled/2p2c/2p2cproblem.hh>
#include <ewoms/decoupled/2p/impes/gridadaptionindicator2p.hh>
#include <ewoms/material/fluidmatrixinteractions/2p/linearmaterial.hh>

#include <ewoms/material/fluidsystems/h2oairfluidsystem.hh>
#include <ewoms/material/fluidsystems/h2on2fluidsystem.hh>

#include "test_dec2p2c_spatialparams.hh"

#include <dune/common/fvector.hh>
#include <string>

namespace Ewoms
{

template<class TypeTag>
class Adaptive2p2c;

namespace Properties
{
NEW_TYPE_TAG(Adaptive2p2c, INHERITS_FROM(DecoupledTwoPTwoCAdaptive,Test2P2CSpatialParams));

NEW_PROP_TAG(OutputInterval);
NEW_PROP_TAG(SimulationName);

// Set the grid type
SET_PROP(Adaptive2p2c, Grid)
{
#if HAVE_ALUGRID
    typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> type;
// #elif HAVE_UG
//     typedef Dune::UGGrid<2> type;
#endif

};
// set the GridCreator property
SET_TYPE_PROP(Adaptive2p2c, GridCreator, CubeGridCreator<TypeTag>);

// Set the problem property
SET_PROP(Adaptive2p2c, Problem)
{
    typedef Ewoms::Adaptive2p2c<TTAG(Adaptive2p2c)> type;
};

SET_SCALAR_PROP(Adaptive2p2c, EndTime, 3e3);
SET_SCALAR_PROP(Adaptive2p2c, InitialTimeStepSize, 200);
SET_INT_PROP(Adaptive2p2c, OutputInterval, 1);
SET_STRING_PROP(Adaptive2p2c, SimulationName, "test_adaptive2p2c");

// Select fluid system
SET_PROP(Adaptive2p2c, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Ewoms::FluidSystems::H2OAir<Scalar, Ewoms::SimpleH2O<Scalar> > type;
};

// Specify indicator
SET_TYPE_PROP(Adaptive2p2c, GridAdaptIndicator, GridAdaptionIndicator2P<TypeTag>);

// Enable gravity
SET_BOOL_PROP(Adaptive2p2c, EnableGravity, true);
SET_INT_PROP(Adaptive2p2c,
        BoundaryMobility,
        GET_PROP_TYPE(TypeTag, Indices)::permDependent);
SET_BOOL_PROP(Adaptive2p2c, EnableCapillarity, true);
SET_INT_PROP(Adaptive2p2c, PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureNW);

SET_INT_PROP(Adaptive2p2c, GridAdaptMinLevel, 0);
SET_INT_PROP(Adaptive2p2c, GridAdaptMaxLevel, 1);
SET_BOOL_PROP(Adaptive2p2c, EnableMultiPointFluxApproximation, true);
SET_BOOL_PROP(Adaptive2p2c, MpfaEnableSecondHalfEdge, true);

SET_SCALAR_PROP(Adaptive2p2c, ImpetCflFactor, 0.8);
SET_SCALAR_PROP(Adaptive2p2c, ImpetErrorTermLowerBound, 0.2);
SET_SCALAR_PROP(Adaptive2p2c, ImpetErrorTermUpperBound, 0.9);
SET_SCALAR_PROP(Adaptive2p2c, ImpetErrorTermFactor, 0.5);
SET_BOOL_PROP(Adaptive2p2c, ImpetEnableVolumeIntegral, true);
SET_BOOL_PROP(Adaptive2p2c, ImpetRestrictFluxInTransport, true);

SET_SCALAR_PROP(Adaptive2p2c, DomainSizeX, 10.0);
SET_SCALAR_PROP(Adaptive2p2c, DomainSizeY, 10.0);
SET_SCALAR_PROP(Adaptive2p2c, DomainSizeZ, 0.0);

SET_SCALAR_PROP(Adaptive2p2c, CellsX, 10);
SET_SCALAR_PROP(Adaptive2p2c, CellsY, 10);
SET_SCALAR_PROP(Adaptive2p2c, CellsZ, 0);
}

/*!
 * \ingroup Adaptive2p2cs
 * \ingroup IMPETtests
 *
 * \brief The problem file for the semi-implicit two-phase two-component test with grid-adaption.
 */
template<class TypeTag = TTAG(Adaptive2p2c)>
class Adaptive2p2c: public IMPETProblem2P2C<TypeTag>
{
typedef IMPETProblem2P2C<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;


// boundary typedefs
typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

enum
{
    dim = GridView::dimension,
    dimWorld = GridView::dimensionworld
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename Grid::Traits::template Codim<0>::EntityPointer ElementPointer;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
Adaptive2p2c(TimeManager &timeManager) :
    ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView()),
    debugWriter_(this->gridView(), "gridAfterAdapt")
{
    this->setGrid(GridCreator::grid());
    std::string s = GET_PARAM(TypeTag, std::string, SimulationName);
    this->setName(s.c_str());
    this->setOutputInterval(GET_PARAM(TypeTag, int, OutputInterval));
    // initialize the tables of the fluid system
    FluidSystem::init();
}

    static void registerParameters()
    {
        ParentType::registerParameters();

        REGISTER_PARAM(TypeTag, std::string, SimulationName, "The name of the simulation used for writing result files");
        REGISTER_PARAM(TypeTag, int, OutputInterval, "The number of time steps between two VTK files get written to disk");
    }

//void preTimeStep()
//{
//    ParentType::preTimeStep();
//            // use second writer
//            debugWriter_.gridChanged();
//            // write
//            debugWriter_.beginWrite(this->timeManager().time());
//            //write stuff out
//            typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::ScalarSolution ScalarSolutionType;
//            typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;
//            int size = this->gridView().size(0);
//            ScalarSolutionType *pressureW = debugWriter_.allocateManagedBuffer (size);
//            ScalarSolutionType *pressureN = debugWriter_.allocateManagedBuffer (size);
//            ScalarSolutionType *totalConcentration1 = debugWriter_.allocateManagedBuffer (size);
//            ScalarSolutionType *totalConcentration2 = debugWriter_.allocateManagedBuffer (size);
//            for (int i = 0; i < size; i++)
//            {
//                CellData& cellData = this->variables().cellData(i);
//                (*pressureW)[i] = cellData.pressure(wPhaseIdx);
//                (*pressureN)[i] = cellData.pressure(nPhaseIdx);
//                (*totalConcentration1)[i] = cellData.massConcentration(wPhaseIdx);
//                (*totalConcentration2)[i] = cellData.massConcentration(nPhaseIdx);
//            }
//            debugWriter_.attachCellData(*pressureW, "wetting pressure");
//            debugWriter_.attachCellData(*pressureN, "nonwetting pressure");
//            debugWriter_.attachCellData(*totalConcentration1, "C^w from cellData");
//            debugWriter_.attachCellData(*totalConcentration2, "C^n from cellData");
//            debugWriter_.endWrite();
//            return;
//}

/*!
 * \name Problem parameters
 */
// \{

bool shouldWriteRestartFile() const
{
    return false;
}

//! Returns the temperature within the domain.
/*! This problem assumes a temperature of 10 degrees Celsius.
 * \param globalPos The global Position
 */
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::referencePressureAtPos()
 */
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e6;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::boundaryTypesAtPos()
 */
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > this->bboxMax()[0]-1E-6 || globalPos[0] < 1e-6)
        bcTypes.setAllDirichlet();
    else
        // all other boundaries
        bcTypes.setAllNeumann();
}

/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::boundaryFormulation()
 */
const void boundaryFormulation(typename Indices::BoundaryFormulation &bcFormulation, const Intersection& intersection) const
{
    bcFormulation = Indices::concentration;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::dirichletAtPos()
 */
void dirichletAtPos(PrimaryVariables &bcValues, const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    // Dirichlet for pressure equation
    bcValues[Indices::pressureEqIdx] = (globalPos[0] < 1e-6) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);

    // Dirichlet values for transport equations
    bcValues[Indices::contiWEqIdx] = 1.;
    bcValues[Indices::contiNEqIdx] = 1.- bcValues[Indices::contiWEqIdx];

}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::neumannAtPos()
 */
void neumannAtPos(PrimaryVariables &neumannValues, const GlobalPosition& globalPos) const
{
    this->setZero(neumannValues, Indices::contiWEqIdx);
}
/*!
 * \copydoc Ewoms::IMPETProblem::source()
 */
void source(PrimaryVariables &values, const Element &element)
{
    this->setZero(values, Indices::contiWEqIdx);
    ElementPointer father(element);
    // access level 1 entity
    while (father->level() != this->gridAdapt().getMinLevel())
    {
        father = father->father();
    }
    GlobalPosition globalPos = father->geometry().center();
    if (fabs(globalPos[0] - 4.8) < 0.5 && fabs(globalPos[1] - 4.8) < 0.5)
        values[Indices::contiNEqIdx] = 0.0001;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::initialFormulation()
 */
const void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::concentration;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::initConcentrationAtPos()
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1.0;
}

private:
Grid grid_;
Ewoms::VtkMultiWriter<GridView> debugWriter_;
};
} //end namespace

#endif
