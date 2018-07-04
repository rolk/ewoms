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
 * \copydoc Ewoms::Co2OctaneFluidSystem
 */
#ifndef OPM_CO2_OCTANE_FLUID_SYSTEM_HPP
#define OPM_CO2_OCTANE_FLUID_SYSTEM_HPP

#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>

#include <opm/material/IdealGas.hpp>
#include "octane.hh"
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/binarycoefficients/H2O_CO2.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/material/common/Exceptions.hpp>

#include <iostream>
#include <cassert>

namespace Ewoms {
/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with water and nitrogen as components.
 */
template <class Scalar>
class Co2OctaneFluidSystem
    : public Opm::BaseFluidSystem<Scalar, Co2OctaneFluidSystem<Scalar> >
{
    typedef Co2OctaneFluidSystem<Scalar> ThisType;
    typedef Opm::BaseFluidSystem<Scalar, ThisType> Base;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = Opm::NullParameterCache<Evaluation>;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the liquid phase
    static const int liquidPhaseIdx = 0;
    //! Index of the gas phase
    static const int gasPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {
            "liquid",
            "gas"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        return
            (phaseIdx == gasPhaseIdx)
            ? true
            :Octane::liquidIsCompressible();// the water component decides for the liquid phase...
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        return
            (phaseIdx == gasPhaseIdx)
            ? Octane::gasIsIdeal() && CO2::gasIsIdeal() // let the components decide
            : false; // not a gas
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! The component index of water
    static const int OctaneIdx = 0;
    //! The component index of molecular nitrogen
    static const int CO2Idx = 1;

    //! The component for pure water
    typedef Ewoms::Octane<Scalar> Octane;

    //! The component for pure nitrogen
    typedef Opm::SimpleCO2<Scalar> CO2;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        static const char* name[] = {
            Octane::name(),
            CO2::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == OctaneIdx)
            ? Octane::molarMass()
            : (compIdx == CO2Idx)
            ? CO2::molarMass()
            : 1e30;
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::criticalTemperature()
            : (compIdx == CO2Idx)
            ? CO2::criticalTemperature()
            : 1e30;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::criticalPressure()
            : (compIdx == CO2Idx)
            ? CO2::criticalPressure()
            : 1e30;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::acentricFactor()
            : (compIdx == CO2Idx)
            ? CO2::acentricFactor()
            : 1e30;
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * If a tabulated Octane component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/50,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/50);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (Octane::isTabulated) {
            Octane::init(tempMin, tempMax, nTemp,
                         pressMin, pressMax, nPress);
        }
    }

    /*!
     * \copydoc BaseFluidSystem::density
     */
#warning "TODO: hardcore thermodynamics"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        LhsEval sumMoleFrac = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoleFrac += Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, compIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            // assume ideal mixture where each molecule occupies the same volume regardless
            // of whether it is water or nitrogen.
            const LhsEval& clOctane = Octane::liquidDensity(T, p)/Octane::molarMass();

            const auto& xlOctane = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, OctaneIdx));
            const auto& xlCO2 = Opm::decay<LhsEval>(fluidState.moleFraction(liquidPhaseIdx, CO2Idx));

            return clOctane*(Octane::molarMass()*xlOctane + CO2::molarMass()*xlCO2)/sumMoleFrac;
        }

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        // assume ideal mixture: steam and nitrogen don't "distinguish" each other
        const auto& xgOctane = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, OctaneIdx));
        const auto& xgCO2 = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, CO2Idx));
        const auto& rho_gOctane = Octane::gasDensity(T, p*xgOctane);
        const auto& rho_gCO2 = CO2::gasDensity(T, p*xgCO2);
        return (rho_gOctane + rho_gCO2)/Opm::max(1e-5, sumMoleFrac);
    }

    //! \copydoc BaseFluidSystem::viscosity
#warning "TODO: hardcore thermodynamics"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            // assume pure water for the liquid phase
            return Octane::liquidViscosity(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 20001, p. 9.21/22
         */
        LhsEval muResult = 0;
        const LhsEval mu[numComponents] = {
            Octane::gasViscosity(T, Octane::vaporPressure(T)),
            CO2::gasViscosity(T, p)
        };

        LhsEval sumx = 0.0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumx += Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, compIdx));
        sumx = Opm::max(1e-10, sumx);

        for (unsigned i = 0; i < numComponents; ++i) {
            LhsEval divisor = 0;
            for (unsigned j = 0; j < numComponents; ++j) {
                LhsEval phiIJ = 1 + Opm::sqrt(mu[i]/mu[j]) * std::pow(molarMass(j)/molarMass(i), 1/4.0);
                phiIJ *= phiIJ;
                phiIJ /= std::sqrt(8*(1 + molarMass(i)/molarMass(j)));
                divisor +=
                    Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, j))
                    /sumx*phiIJ;
            }
            muResult +=
                Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, i))
                /sumx*mu[i]/divisor;
        }
        return muResult;
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
#warning "TODO: hardcore thermodynamics"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            if (compIdx == OctaneIdx)
                return Octane::vaporPressure(T)/p;
            return Opm::BinaryCoeff::H2O_CO2::henry<Scalar>(T)/p;
        }

        assert(phaseIdx == gasPhaseIdx);

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
#warning "TODO (currently unused)"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned /*compIdx*/)

    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // liquid phase
        if (phaseIdx == liquidPhaseIdx)
            return Opm::BinaryCoeff::H2O_CO2::liquidDiffCoeff(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);
        return Opm::BinaryCoeff::H2O_CO2::gasDiffCoeff(T, p);
    }

    //! \copydoc BaseFluidSystem::enthalpy
#warning "TODO (currently unused)"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        Opm::Valgrind::CheckDefined(T);
        Opm::Valgrind::CheckDefined(p);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            // TODO: correct way to deal with the solutes???
            return Octane::liquidEnthalpy(T, p);
        }

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        // assume ideal mixture: Molecules of one component don't discriminate between
        // their own kind and molecules of the other component.
        const auto& XgOctane = Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, OctaneIdx));
        const auto& XgCO2 = Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, CO2Idx));

        LhsEval hOctane = XgOctane*Octane::gasEnthalpy(T, p);
        LhsEval hCO2 = XgCO2*CO2::gasEnthalpy(T, p);
        return hOctane + hCO2;
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
#warning "TODO (currently unused)"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        if (phaseIdx == liquidPhaseIdx) // liquid phase
            return Octane::liquidThermalConductivity(T, p);

        // gas phase
        assert(phaseIdx == gasPhaseIdx);

        // return the sum of the partial conductivity of Nitrogen and Steam
        const auto& xOctane = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, OctaneIdx));
        const auto& xCO2 = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));

        // Assuming Raoult's, Daltons law and ideal gas in order to obtain the
        // partial pressures in the gas phase
        const auto& lambdaCO2 = CO2::gasThermalConductivity(T, p*xCO2);
        const auto& lambdaOctane = Octane::gasThermalConductivity(T, p*xOctane);

        return lambdaCO2 + lambdaOctane;
    }

    //! \copydoc BaseFluidSystem::heatCapacity
#warning "TODO (currently unused)"
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval heatCapacity(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& xAlphaOctane = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, OctaneIdx));
        const auto& xAlphaCO2 = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));
        const auto& XAlphaOctane = Opm::decay<LhsEval>(fluidState.massFraction(phaseIdx, OctaneIdx));
        const auto& XAlphaCO2 = Opm::decay<LhsEval>(fluidState.massFraction(phaseIdx, CO2Idx));

        if (phaseIdx == liquidPhaseIdx)
            return Octane::liquidHeatCapacity(T, p);

        assert(phaseIdx == gasPhaseIdx);

        // for the gas phase, assume ideal mixture
        LhsEval c_pCO2;
        LhsEval c_pOctane;

        c_pCO2 = CO2::gasHeatCapacity(T, p*xAlphaCO2);
        c_pOctane = Octane::gasHeatCapacity(T, p*xAlphaOctane);

        // mingle both components together. this assumes that there is no "cross
        // interaction" between both flavors of molecules.
        return XAlphaOctane*c_pOctane + XAlphaCO2*c_pCO2;
    }
};

} // namespace Opm

#endif
