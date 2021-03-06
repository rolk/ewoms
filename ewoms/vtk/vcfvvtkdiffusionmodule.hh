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
 * \copydoc Ewoms::VcfvVtkDiffusionModule
 */
#ifndef EWOMS_VCFV_VTK_DIFFUSION_MODULE_HH
#define EWOMS_VCFV_VTK_DIFFUSION_MODULE_HH

#include "vcfvvtkoutputmodule.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

namespace Ewoms
{
namespace Properties
{
// create new type tag for the VTK output of the quantities for molecular diffusion
NEW_TYPE_TAG(VtkDiffusion);

// create the property tags needed for the diffusion module
NEW_PROP_TAG(VtkWriteTortuosities);
NEW_PROP_TAG(VtkWriteDiffusionCoefficients);
NEW_PROP_TAG(VtkWriteEffectiveDiffusionCoefficients);

// set default values for what quantities to output
SET_BOOL_PROP(VtkDiffusion, VtkWriteTortuosities, false);
SET_BOOL_PROP(VtkDiffusion, VtkWriteDiffusionCoefficients, false);
SET_BOOL_PROP(VtkDiffusion, VtkWriteEffectiveDiffusionCoefficients, false);
}

/*!
 * \ingroup VcfvVtk
 *
 * \brief VTK output module for quantities which make sense for models which
 *        incorperate molecular diffusion.
 *
 * This module deals with the following quantities:
 * - Molecular diffusion coefficients of all components in all fluid phases
 * - Effective molecular diffusion coefficients of the porous medium of all components in all fluid phases
 */
template<class TypeTag>
class VcfvVtkDiffusionModule : public VcfvVtkOutputModule<TypeTag>
{
    typedef VcfvVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename ParentType::PhaseComponentBuffer PhaseComponentBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;
    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

public:
    VcfvVtkDiffusionModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, bool, VtkWriteTortuosities, "Include the tortuosity for each phase in the VTK output files");
        REGISTER_PARAM(TypeTag, bool, VtkWriteDiffusionCoefficients, "Include the molecular diffusion coefficients in the VTK output files");
        REGISTER_PARAM(TypeTag, bool, VtkWriteEffectiveDiffusionCoefficients, "Include the effective molecular diffusion coefficients the medium in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (tortuosityOutput_()) this->resizePhaseBuffer_(tortuosity_);
        if (diffusionCoefficientOutput_()) this->resizePhaseComponentBuffer_(diffusionCoefficient_);
        if (effectiveDiffusionCoefficientOutput_()) this->resizePhaseComponentBuffer_(effectiveDiffusionCoefficient_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        const auto &vertexMapper = elemCtx.problem().vertexMapper();
        const auto &elem = elemCtx.element();

        for (int i = 0; i < elemCtx.numScv(); ++i) {
            int I = vertexMapper.map(elem, i, dim);
            const auto &volVars = elemCtx.volVars(i, /*timeIdx=*/0);

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (tortuosityOutput_()) tortuosity_[phaseIdx][I] = volVars.tortuosity(phaseIdx);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (diffusionCoefficientOutput_()) diffusionCoefficient_[phaseIdx][compIdx][I] = volVars.diffusionCoefficient(phaseIdx, compIdx);
                    if (effectiveDiffusionCoefficientOutput_()) effectiveDiffusionCoefficient_[phaseIdx][compIdx][I] = volVars.effectiveDiffusionCoefficient(phaseIdx, compIdx);
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (tortuosityOutput_()) this->commitPhaseBuffer_(writer, "tortuosity", tortuosity_);
        if (diffusionCoefficientOutput_()) this->commitPhaseComponentBuffer_(writer, "diffusionCoefficient", diffusionCoefficient_);
        if (effectiveDiffusionCoefficientOutput_()) this->commitPhaseComponentBuffer_(writer, "effectiveDiffusionCoefficient", effectiveDiffusionCoefficient_);
    }

private:
    static bool tortuosityOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteTortuosities); }

    static bool diffusionCoefficientOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteDiffusionCoefficients); }

    static bool effectiveDiffusionCoefficientOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteEffectiveDiffusionCoefficients); }

    PhaseBuffer tortuosity_;
    PhaseComponentBuffer diffusionCoefficient_;
    PhaseComponentBuffer effectiveDiffusionCoefficient_;
};

}

#endif
