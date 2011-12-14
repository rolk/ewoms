// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 */
#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include "2p2cproperties.hh"
#include "2p2clocalresidual.hh"
#include "2p2cproblem.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{{\bf D}_{\alpha, pm}^\kappa \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mass fraction of, e.g., air in the wetting phase \f$X^a_w\f$ is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^a_w<X^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^w_n<X^w_{n,max})\f$</li>
 * </ul>
 */
template<class TypeTag>
class TwoPTwoCModel: public BoxModel<TypeTag>
{
    typedef TwoPTwoCModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases,

        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        setSwitched_(false);
    }

    /*!
     * \brief Returns a string with the model's human-readable name
     */
    std::string name() const
    { return "2p2c"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        if (pvIdx == Indices::pressureIdx) {
            if (formulation == plSg)
                oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
            else
                oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/1);
        }
        else if (pvIdx == Indices::switchIdx)
            oss << "switch";
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    std::string eqName(int eqIdx) const
    { 
        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "continuity^" << FluidSystem::componentName(/*compIdx=*/1);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \param dest Contains the storage of each component for one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(PrimaryVariables &dest, int phaseIdx)
    {
        dest = 0;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);

            this->localResidual().addPhaseStorage(dest, elemCtx, /*timeIdx=*/0, phaseIdx);
        };

        this->gridView_().comm().sum(dest);
    }

    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();
        setSwitched_(false);
    };

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressureIdx == pvIdx)
            return std::min(1.0/this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx], 1.0);
        return 1;
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();
        setSwitched_(false);
    }

    /*!
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param vert The vertex
     */
    void serializeEntity(std::ostream &outStream, const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);

        int vertIdx = this->dofMapper().map(vert);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);

        outStream << this->solution(/*timeIdx=*/0)[vertIdx].phasePresence() << " ";
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param vert The vertex
     */
    void deserializeEntity(std::istream &inStream, const Vertex &vert)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vert);

        // read phase presence
        int vertIdx = this->dofMapper().map(vert);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << vertIdx);

        int tmp;
        inStream >> tmp;
        this->solution(/*timeIdx=*/0)[vertIdx].setPhasePresence(tmp);
        this->solution(/*timeIdx=*/1)[vertIdx].setPhasePresence(tmp);
    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     *
     * \param curGlobalSol The current global solution
     */
    void updatePhasePresence_(SolutionVector &curGlobalSol)
    {
        bool wasSwitched = false;

        std::vector<bool> visited(curGlobalSol.size(), false);

        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            bool fvElemGeomUpdated = false;
            int numScv = elemIt->template count<dim>();
            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, scvIdx, dim);

                if (visited[globalIdx])
                    continue;
                visited[globalIdx] = true;

                if (!fvElemGeomUpdated) {
                    fvElemGeomUpdated = true;
                    elemCtx.updateFVElemGeom(*elemIt);
                }

                // compute the volume variables of the current
                // sub-control volume
                elemCtx.updateScvVars(curGlobalSol[globalIdx],
                                       scvIdx,
                                       /*timeIdx=*/0);
                const VolumeVariables &volVars = elemCtx.volVars(scvIdx, /*timeIdx=*/0);

                const GlobalPosition &globalPos = elemCtx.pos(scvIdx, /*timeIdx=*/0);
                if (primaryVarSwitch_(curGlobalSol[globalIdx],
                                      volVars,
                                      globalIdx,
                                      globalPos))
                {
                    this->jacobianAssembler().markVertexRed(globalIdx);
                    wasSwitched = true;
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        if (this->gridView_().comm().size() > 1)
            wasSwitched = this->gridView_().comm().max(wasSwitched);

        setSwitched_(wasSwitched);
    }

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    };

    /*!
     * \brief Set whether there was a primary variable switch after in
     *        the last timestep.
     */
    void setSwitched_(bool yesno)
    { switchFlag_ = yesno; }

    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(PrimaryVariables &priVars,
                           const VolumeVariables &volVars,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        int phasePresence = priVars.phasePresence();
        int newPhasePresence = phasePresence;

        Scalar switchTol = 0.0;
        if (priVars.wasSwitched())
            switchTol = 0.02;

        // check if a primary var switch is necessary
        if (phasePresence == gPhaseOnly)
        {
            // calculate mole fraction in the hypothetic liquid phase
            Scalar xll = volVars.fluidState().moleFraction(/*phaseIdx=*/0, /*compIdx=*/0);
            Scalar xlg = volVars.fluidState().moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);

            // if the sum of the mole fractions would be larger than
            // 100%, liquid phase appears
            if (xll + xlg > 1.0 + switchTol)
            {
                // liquid phase appears
                std::cout << "liquid phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xll + xlg: "
                          << xll + xlg << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pgSl)
                    priVars[switchIdx] = 0.0;
                else if (formulation == plSg)
                    priVars[switchIdx] = 1.0;
            };
        }
        else if (phasePresence == lPhaseOnly)
        {
            // calculate fractions of the partial pressures in the
            // hypothetic gas phase
            Scalar xgl = volVars.fluidState().moleFraction(/*phaseIdx=*/1, /*compIdx=*/0);
            Scalar xgg = volVars.fluidState().moleFraction(/*phaseIdx=*/1, /*compIdx=*/1);

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xgl + xgg > 1.0 + switchTol)
            {
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", x_gl + x_gg: "
                          << xgl + xgg << std::endl;
                std::cout.flush();
                newPhasePresence = bothPhases;
                if (formulation == pgSl)
                    priVars[switchIdx] = 0.999;
                else if (formulation == plSg)
                    priVars[switchIdx] = 0.001;
            }
        }
        else if (phasePresence == bothPhases)
        {
            if (volVars.fluidState().saturation(/*phaseIdx=*/1) <= 0.0 - switchTol)
            {
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", Sg: "
                          << volVars.fluidState().saturation(/*phaseIdx=*/1)
                          << " X_0^1: " << volVars.fluidState().massFraction(/*phaseIdx=*/0, /*compIdx=*/1)
                          << std::endl;
                newPhasePresence = lPhaseOnly;

                priVars[switchIdx] =
                    volVars.fluidState().massFraction(/*phaseIdx=*/0, /*compIdx=*/1);
            }
            else if (volVars.fluidState().saturation(/*phaseIdx=*/0) <= 0.0 - switchTol)
            {
                // liquid phase disappears
                std::cout << "Liquid phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", Sl: "
                          << volVars.fluidState().saturation(/*phaseIdx=*/0)
                          << " X_1^0: " << volVars.fluidState().massFraction(/*phaseIdx=*/1, /*compIdx=*/0)
                          << std::endl;
                newPhasePresence = gPhaseOnly;

                priVars[switchIdx] =
                    volVars.fluidState().massFraction(/*phaseIdx=*/1, /*compIdx=*/0);
            }
        }

        priVars.setPhasePresence(newPhasePresence);
        priVars.setSwitched(phasePresence != newPhasePresence);

        return phasePresence != newPhasePresence;
    }

protected:
    // parameters given in constructor
    bool switchFlag_;
};

}

#include "2p2cpropertydefaults.hh"

#endif
