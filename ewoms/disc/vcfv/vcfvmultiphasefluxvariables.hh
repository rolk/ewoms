// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 *
 * \copydoc Ewoms::VcfvMultiPhaseFluxVariables
 */
#ifndef EWOMS_VCFV_MULTIPHASE_FLUX_VARIABLES_HH
#define EWOMS_VCFV_MULTIPHASE_FLUX_VARIABLES_HH

#include <ewoms/common/spline.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/fvector.hh>

#include "vcfvproperties.hh"

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(EnableSmoothUpwinding);
NEW_PROP_TAG(MaterialLaw);
NEW_PROP_TAG(EnableGravity);
NEW_PROP_TAG(VelocityModule);
}

/*!
 * \ingroup VcfvModel
 *
 * \brief This class calculates the pressure potential gradients and
 *        the filter velocities for multi-phase flow in porous media
 */
template <class TypeTag>
class VcfvMultiPhaseFluxVariables
    : public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { useTwoPointGradients = GET_PROP_VALUE(TypeTag, UseTwoPointGradients) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;
    typedef typename VelocityModule::VelocityFluxVariables VelocityFluxVariables;

public:
    /*!
     * \brief Register all run-time parameters for the flux variables.
     */
    static void registerParameters()
    {
        VelocityModule::registerParameters();
    }

    /*!
     * \brief Update the flux variables for a given sub-control-volume-face.
     *
     * \param elemCtx Reference to the current element context.
     * \param scvfIdx The local index of the sub-control-volume face for which the flux variables should be calculated.
     * \param timeIdx The index used by the time discretization.
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];
        insideScvIdx_ = scvf.i;
        outsideScvIdx_ = scvf.j;

        extrusionFactor_ =
            (elemCtx.volVars(insideScvIdx_, timeIdx).extrusionFactor()
             + elemCtx.volVars(outsideScvIdx_, timeIdx).extrusionFactor()) / 2;
        Valgrind::CheckDefined(extrusionFactor_);
        assert(extrusionFactor_ > 0);

        // compute the pressure potential gradients
        calculateGradients_(elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(potentialGrad_);

        // determine the upstream indices. since this is a semi-smooth
        // non-linear solver, make upstream only look at the
        // evaluation point for the upstream decision
        const auto &evalFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);
        if (&evalFluxVars == this) {
            // we _are_ the evaluation point. Check whether the
            // pressure potential is in the same direction as the face
            // normal or in the opposite one
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                    upstreamScvIdx_[phaseIdx] = insideScvIdx_;
                    downstreamScvIdx_[phaseIdx] = outsideScvIdx_;
                    continue;
                }

                if (potentialGrad(phaseIdx) * scvf.normal > 0) {
                    upstreamScvIdx_[phaseIdx] = outsideScvIdx_;
                    downstreamScvIdx_[phaseIdx] = insideScvIdx_;
                }
                else {
                    upstreamScvIdx_[phaseIdx] = insideScvIdx_;
                    downstreamScvIdx_[phaseIdx] = outsideScvIdx_;
                }
            }
        }
        else {
            // we are *not* the evaluation point. in this case, we
            // just take the up-/downstream indices from the
            // evaluation point.
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                upstreamScvIdx_[phaseIdx] = evalFluxVars.upstreamIndex(phaseIdx);
                downstreamScvIdx_[phaseIdx] = evalFluxVars.downstreamIndex(phaseIdx);
            }
        }

        VelocityFluxVariables::calculateVelocities_(elemCtx, scvfIdx, timeIdx);
    }


    /*!
     * \brief Update the flux variables for a given boundary face.
     *
     * \param context Reference to the current execution context.
     * \param bfIdx The local index of the boundary face for which the flux variables should be calculated.
     * \param timeIdx The index used by the time discretization.
     * \param fluidState The FluidState on the domain boundary.
     * \param paramCache The FluidSystem's parameter cache.
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context,
                        int bfIdx,
                        int timeIdx,
                        const FluidState &fluidState,
                        typename FluidSystem::ParameterCache &paramCache)
    {
        int scvIdx = context.insideScvIndex(bfIdx, timeIdx);
        insideScvIdx_ = scvIdx;
        outsideScvIdx_ = scvIdx;

        extrusionFactor_ = context.volVars(bfIdx, timeIdx).extrusionFactor();
        Valgrind::CheckDefined(extrusionFactor_);
        assert(extrusionFactor_ > 0);

        calculateBoundaryGradients_(context, bfIdx, timeIdx, fluidState, paramCache);
        VelocityFluxVariables::calculateBoundaryVelocities_(context, bfIdx, timeIdx, fluidState, paramCache);
    }

    /*!
     * \brief Returns th extrusion factor for the sub-control volume face
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase
     *        at the face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const DimVector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return the local index of the control volume which is on
     *        the "inside" of the sub-control volume face.
     */
    short insideIndex() const
    { return insideScvIdx_; }

    /*!
     * \brief Return the local index of the control volume which is on
     *        the "outside" of the sub-control volume face.
     */
    short outsideIndex() const
    { return outsideScvIdx_; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    short upstreamIndex(int phaseIdx) const
    { return upstreamScvIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    short downstreamIndex(int phaseIdx) const
    { return downstreamScvIdx_[phaseIdx]; }

    /*!
     * \brief Return the weight of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar upstreamWeight(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief Return the weight of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar downstreamWeight(int phaseIdx) const
    { return 1.0 - upstreamWeight(phaseIdx); }

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx,
                             int timeIdx)
    {
        const auto &fvElemGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.subContVolFace[scvfIdx];

        if (useTwoPointGradients) {
            const auto &fluidStateI = elemCtx.volVars(insideScvIdx_, timeIdx).fluidState();
            const auto &fluidStateJ = elemCtx.volVars(outsideScvIdx_, timeIdx).fluidState();
            const auto &scvI = fvElemGeom.subContVol[insideScvIdx_];
            const auto &scvJ = fvElemGeom.subContVol[outsideScvIdx_];

            // the "normalized normal" of the scvf divided by the
            // distance of the centers of the two adjacent SCVs
            DimVector n = scvf.normal;
            n /= scvf.normal.two_norm();

            // distance between the centers of the two SCVs
            Scalar dist = 0;
            for (int i = 0; i < dimWorld; ++ i) {
                Scalar tmp = scvI.global[i] - scvJ.global[i];
                dist += tmp*tmp;
            }
            dist = std::sqrt(dist);

            // calculate the pressure gradient
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                    potentialGrad_[phaseIdx] = 0;
                    continue;
                }

                potentialGrad_[phaseIdx] = n;
                potentialGrad_[phaseIdx] *= (fluidStateJ.pressure(phaseIdx) - fluidStateI.pressure(phaseIdx))/dist;
            }
        }
        else {
            // reset all gradients to 0
            for (int phase = 0; phase < numPhases; ++phase) {
                potentialGrad_[phase] = Scalar(0);
            }
            Valgrind::CheckDefined(potentialGrad_);

            // calculate gradients
            for (int scvIdx = 0;
                 scvIdx < elemCtx.numScv();
                 scvIdx ++)
            {
                // FE gradient at vertex
                const DimVector &feGrad = scvf.grad[scvIdx];
                const auto &fluidState = elemCtx.volVars(scvIdx, timeIdx).fluidState();
                Valgrind::CheckDefined(feGrad);

                // compute sum of pressure gradients for each phase
                for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
                {
                    if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                        continue;

                    // the pressure gradient
                    DimVector tmp(feGrad);
                    tmp *= fluidState.pressure(phaseIdx);
                    Valgrind::CheckDefined(tmp);
                    potentialGrad_[phaseIdx] += tmp;
                    Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
                }
            }
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            DimVector g(elemCtx.problem().gravity(elemCtx, insideScvIdx_, timeIdx));
            g += elemCtx.problem().gravity(elemCtx, outsideScvIdx_, timeIdx);
            g /= 2;
            Valgrind::CheckDefined(g);

            const auto &fluidStateIn = elemCtx.volVars(insideScvIdx_, timeIdx).fluidState();
            const auto &fluidStateOut = elemCtx.volVars(outsideScvIdx_, timeIdx).fluidState();
            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = fluidStateIn.saturation(phaseIdx);
                Scalar SJ = fluidStateOut.saturation(phaseIdx);
                Scalar rhoI = fluidStateIn.density(phaseIdx);
                Scalar rhoJ = fluidStateOut.density(phaseIdx);
                Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                if (fI + fJ == 0)
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fI = fJ = 0.5;
                Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);
                Valgrind::CheckDefined(density);

                // make gravity acceleration a force
                DimVector f(g);
                f *= density;

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
                if (!std::isfinite(potentialGrad_[phaseIdx].two_norm())) {
                    DUNE_THROW(NumericalProblem,
                               "Non finite potential gradient for phase '"
                               << FluidSystem::phaseName(phaseIdx) << "'");
                }
            }
        }
    }

    template <class Context, class FluidState>
    void calculateBoundaryGradients_(const Context &context,
                                     int bfIdx,
                                     int timeIdx,
                                     const FluidState &fluidState,
                                     const typename FluidSystem::ParameterCache &paramCache)
    {
        const auto &fvElemGeom = context.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.boundaryFace[bfIdx];

        const auto &elemCtx = context.elemContext();
        const auto &insideScv = elemCtx.fvElemGeom(timeIdx).subContVol[insideScvIdx_];

        const auto &fluidStateI = elemCtx.volVars(insideScvIdx_, timeIdx).fluidState();
        const auto &fluidStateJ = fluidState;

        // the "normalized normal" of the scvf divided by the
        // distance of the centers of the two adjacent SCVs
        DimVector n = scvf.normal;
        n /= scvf.normal.two_norm();

        // distance between the center of the SCV and center of the boundary face
        DimVector distVec = scvf.ipGlobal;
        distVec -= context.element().geometry().global(insideScv.localGeometry->center());
        Scalar dist = distVec * n;

        // if the following assertation triggers, the center of the
        // center of the interior SCV was not inside the element!
        assert(dist > 0);

        // calculate the pressure gradient using two-point gradient
        // appoximation
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                potentialGrad_[phaseIdx] = 0;
                continue;
            }

            potentialGrad_[phaseIdx] = n;
            potentialGrad_[phaseIdx] *= (fluidStateJ.pressure(phaseIdx) - fluidStateI.pressure(phaseIdx))/dist;
        }

        fluidStateI.checkDefined();

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            DimVector g(context.problem().gravity(context.elemContext(), insideScvIdx_, timeIdx));

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate volumetric gravity acceleration force
                DimVector f(g);
                f *= fluidStateI.density(phaseIdx);

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
            }
        }
    }

#if 0
    void handleSmoothUpwinding_(const ElementContext &elemCtx,
                                int scvfIdx,
                                int timeIdx,
                                int phaseIdx,
                                const DimVector &normal)
    {
        const VolumeVariables &up = elemCtx.volVars(upstreamIndex(phaseIdx), timeIdx);
        const VolumeVariables &dn = elemCtx.volVars(downstreamIndex(phaseIdx), timeIdx);

        // first, calculate the component of the "prelimary velocity"
        // which is parallel to the normal of the sub-control volume
        // face
        DimVector parV(normal);
        parV *= normal * filterVelocity_[phaseIdx];
        Scalar x = parV.two_norm();
        assert(x >= 0);

        if (x == 0.0) {
            filterVelocity_[phaseIdx] = 0.0;
            return;
        }

        // slopes of the velocity in upwind, downwind direction and at
        // x = 0
        Scalar mUp = up.mobility(phaseIdx);
        Scalar mDn = dn.mobility(phaseIdx);
        Scalar m0 = Ewoms::harmonicMean(mUp, mDn);

        Scalar maxMob = std::max(mUp, mDn);
        if (maxMob < 1e-8) {
            filterVelocity_[phaseIdx] = 0.0;
            return;
        }

        // put the mean viscosity and permeanbility in
        // relation to the viscosity of water at
        // approximatly 20 degrees Celsius.
        const Scalar pGradRef = 100; // [Pa/m]
        const Scalar mobRef = 1.0/1e-3; // [1 / (Pa s)]
        const Scalar KRef = 1e-12; // [m^2] = approx 1 Darcy
        const Scalar vRef = mobRef * KRef * pGradRef; // [m/s]

        // calculate the velocity below which smooth upwinding should
        // kick in. For this, we assume that the reference medium is
        // fully saturated with liquid water.
        Scalar eps = vRef / maxMob;
        if (x > eps) {
            // we only do tricks if x is below the epsilon
            // value. Here, this is not the case, so we use
            // the full upwinding scheme
            filterVelocity_[phaseIdx] *= mUp;
        }
        else {
            // interpolate between zero and epsilon using a cubic
            // spline
            Ewoms::Spline<Scalar> sp(/*x0=*/0.0,
                                     /*x1=*/eps,
                                     /*y0=*/0.0,
                                     /*y1=*/mUp*eps,
                                     /*m0=*/m0,
                                     /*m1=*/mUp);

            // set the length of the velocity component which is
            // parallel to the face normal to the one from the spline,
            // and do not modify the perpendicular part
            DimVector perpV = filterVelocity_[phaseIdx];
            perpV -= parV;
            perpV *= mUp;

            parV *= sp.eval(x)/parV.two_norm();

            filterVelocity_[phaseIdx] = perpV;
            filterVelocity_[phaseIdx] += parV;
        }
    }
#endif

    // extrusion factor for the sub-control volume face
    Scalar extrusionFactor_;

    // local indices of the inside and the outside sub-control volumes
    short insideScvIdx_;
    short outsideScvIdx_;

    short upstreamScvIdx_[numPhases];
    short downstreamScvIdx_[numPhases];

    // pressure potential gradients of all phases
    DimVector potentialGrad_[numPhases];
};

} // end namepace

#endif
