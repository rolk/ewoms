// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::Co2FingerProblem
 */
#ifndef EWOMS_CO2_FINGER_PROBLEM_HH
#define EWOMS_CO2_FINGER_PROBLEM_HH

#include <ewoms/models/immiscible/immisciblemodel.hh>
#include <ewoms/linear/parallelamgbackend.hh>
#include "co2octanefluidsystem.hh"

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class Co2FingerProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Co2FingerBaseProblem);

// Set the grid type
#if HAVE_DUNE_ALUGRID
// use dune-alugrid if available
#warning "using dune-alugrid. adaptive grid refinement will be available, but parallelism won't"
SET_TYPE_PROP(Co2FingerBaseProblem,
              Grid,
              Dune::ALUGrid</*dim=*/2,
                            /*dimWorld=*/2,
                            Dune::cube,
                            Dune::nonconforming>);
#else
#warning "using YaspGrid. parallism will be available, but adaptive grid refinement won't"
SET_TYPE_PROP(Co2FingerBaseProblem, Grid, Dune::YaspGrid<2>);
#endif

// declare the CO2 finger problem specific property tags
NEW_PROP_TAG(FluidSystemPressureLow);
NEW_PROP_TAG(FluidSystemPressureHigh);
NEW_PROP_TAG(FluidSystemNumPressure);
NEW_PROP_TAG(FluidSystemTemperatureLow);
NEW_PROP_TAG(FluidSystemTemperatureHigh);
NEW_PROP_TAG(FluidSystemNumTemperature);

NEW_PROP_TAG(MaxDepth);
NEW_PROP_TAG(Temperature);
NEW_PROP_TAG(SimulationName);

// Set the problem property
SET_TYPE_PROP(Co2FingerBaseProblem, Problem,
              Ewoms::Co2FingerProblem<TypeTag>);

// Set fluid configuration
SET_PROP(Co2FingerBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Ewoms::Co2OctaneFluidSystem<Scalar> type;
};

// Set the material Law
SET_PROP(Co2FingerBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffMaterialLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Set the thermal conduction law
SET_PROP(Co2FingerBaseProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::SomertonThermalConductionLaw<FluidSystem, Scalar> type;
};

// set the energy storage law for the solid phase
SET_TYPE_PROP(Co2FingerBaseProblem, SolidEnergyLaw,
              Opm::ConstantSolidHeatCapLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(Co2FingerBaseProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(Co2FingerBaseProblem, EnableGravity, true);

// set the defaults for the problem specific properties
SET_SCALAR_PROP(Co2FingerBaseProblem, FluidSystemPressureLow, 3e7);
SET_SCALAR_PROP(Co2FingerBaseProblem, FluidSystemPressureHigh, 4e7);
SET_INT_PROP(Co2FingerBaseProblem, FluidSystemNumPressure, 100);
SET_SCALAR_PROP(Co2FingerBaseProblem, FluidSystemTemperatureLow, 290);
SET_SCALAR_PROP(Co2FingerBaseProblem, FluidSystemTemperatureHigh, 500);
SET_INT_PROP(Co2FingerBaseProblem, FluidSystemNumTemperature, 100);

SET_SCALAR_PROP(Co2FingerBaseProblem, MaxDepth, 2500);
SET_SCALAR_PROP(Co2FingerBaseProblem, Temperature, 293.15);
SET_STRING_PROP(Co2FingerBaseProblem, SimulationName, "co2finger");

// The default for the end time of the simulation
SET_SCALAR_PROP(Co2FingerBaseProblem, EndTime, 1e4);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(Co2FingerBaseProblem, InitialTimeStepSize, 250);

// The default DGF file to load
SET_STRING_PROP(Co2FingerBaseProblem, GridFile, "co2finger.dgf");

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 */
template <class TypeTag>
class Co2FingerProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { OctaneIdx = FluidSystem::OctaneIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLaw) ThermalConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;
    typedef typename ThermalConductionLaw::Params ThermalConductionLawParams;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2FingerProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;

        pressureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureLow);
        pressureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureHigh);
        nPressure_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumPressure);

        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);

        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        // intrinsic permeabilities
        K_ = this->toDimMatrix_(1e-13);

        // porosities
        porosity_ = 0.3;

        // residual saturations
        materialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        materialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        materialParams_.setEntryPressure(1e4);
        materialParams_.setLambda(2.0);

        materialParams_.finalize();

        // parameters for the somerton law thermal conduction
        computeThermalCondParams_(thermalCondParams_, porosity_);

        // assume constant heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow,
                             "The lower temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh,
                             "The upper temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumTemperature,
                             "The number of intervals between the lower and "
                             "upper temperature");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureLow,
                             "The lower pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureHigh,
                             "The upper pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumPressure,
                             "The number of intervals between the lower and "
                             "upper pressure");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName)
            << "_" << Model::name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";
        oss << "_" << Model::discretizationName();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context OPM_UNUSED,
                                           unsigned spaceIdx OPM_UNUSED,
                                           unsigned timeIdx OPM_UNUSED) const
    { return K_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context  OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    {
        return porosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context OPM_UNUSED,
                                               unsigned spaceIdx OPM_UNUSED,
                                               unsigned timeIdx OPM_UNUSED) const
    {
        return materialParams_;
    }

    /*!
     * \brief Return the parameters for the heat storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyLawParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context OPM_UNUSED,
                            unsigned spaceIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED) const
    {
        return thermalCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (onLeftBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate[contiCO2EqIdx] = -1e-3; // [kg/(m^3 s)]

            typedef Opm::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;
            FluidState fs;
            fs.setSaturation(gasPhaseIdx, 1.0);
            const auto& pg =
                context.intensiveQuantities(spaceIdx, timeIdx).fluidState().pressure(gasPhaseIdx);
            fs.setPressure(gasPhaseIdx, Toolbox::value(pg));
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));

            // impose an forced inflow boundary condition for pure CO2
            values.setMassRate(massRate);

            if (enableEnergy) {
                typename FluidSystem::template ParameterCache<Scalar> paramCache;
                paramCache.updatePhase(fs, gasPhaseIdx);
                Scalar h = FluidSystem::template enthalpy<FluidState, Scalar>(fs, paramCache, gasPhaseIdx);
                values.setEnthalpyRate(massRate[contiCO2EqIdx] * h);
            }
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // const auto& matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        // values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    {
        rate = Scalar(0.0);
        //const auto& pos = context.pos(spaceIdx, time);
    }

    //! \}

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar densityL = FluidSystem::Octane::liquidDensity(temperature_, Scalar(1e5));
        Scalar depth = maxDepth_ - pos[dim - 1];
        Scalar pl = 1e5 - densityL * this->gravity()[dim - 1] * depth;

        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
        fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        //////
        // set composition of the liquid phase
        //////
        fs.setMoleFraction(liquidPhaseIdx, CO2Idx, 0.005);
        fs.setMoleFraction(liquidPhaseIdx, OctaneIdx,
                           1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/liquidPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);
    }

    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onRightBoundary_(pos); }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(liquidPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    DimMatrix K_;

    Scalar porosity_;

    MaterialLawParams materialParams_;

    ThermalConductionLawParams thermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    unsigned nTemperature_;
    unsigned nPressure_;

    Scalar temperatureLow_, temperatureHigh_;
    Scalar pressureLow_, pressureHigh_;
};
} // namespace Ewoms

#endif
