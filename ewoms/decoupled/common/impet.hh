// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Markus Wolff                                 *
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
#ifndef EWOMS_IMPET_HH
#define EWOMS_IMPET_HH

/*!
 * \file
 * \copydoc Ewoms::IMPET
 */

#include "impetproperties.hh"

#include <iostream>

namespace Ewoms
{
/*!
 * \ingroup IMPET
 * \brief IMplicit Pressure Explicit Transport (IMPET) scheme for the solution of weakly coupled diffusion-transport formulations.
 *
 * The model implements the decoupled equations of two-phase flow.
 * These equations can be derived from the two-phase flow equations shown for the two-phase VCVF discretization (TwoPVcfvModel).
 * The first equation to solve is a pressure equation of elliptic character. The second one is a transport equation (e.g. for saturation, concentration,...),
 * which can be hyperbolic or parabolic.
 *
 * As the equations are only weakly coupled they do not have to be solved simultaneously
 * but can be solved sequentially. First the pressure equation is solved implicitly,
 * second the transport equation can be solved explicitly. This solution procedure is called IMPES algorithm
 * (IMplicit Pressure Explicit Saturation) for immiscible flow or IMPEC algorithm
 * (IMplicit Pressure Explicit Concentration) for miscible flow.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class IMPET
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

public:
    typedef typename SolutionTypes::ScalarSolution SolutionType;

    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, Scalar, ImpetErrorTermFactor, "The weight of the calculated error term");
        REGISTER_PARAM(TypeTag, Scalar, ImpetErrorTermLowerBound, "The lower limit of the error term for IMPET");
        REGISTER_PARAM(TypeTag, Scalar, ImpetErrorTermUpperBound, "The upper limit of the error term for IMPET");
        REGISTER_PARAM(TypeTag, int, ImpetIterationFlag, "Enable multiple iterations per time step");
        REGISTER_PARAM(TypeTag, int, ImpetIterationNumber, "The number of iterations done per time step");
        REGISTER_PARAM(TypeTag, Scalar, ImpetMaximumDefect, "The maximum residual tolerated");
        REGISTER_PARAM(TypeTag, Scalar, ImpetRelaxationFactor, "The relaxation factor for the update");
    }

    //! Set initial solution and initialize parameters
    void initialize()
    {
        //initial values of transported quantity
        problem_.transportModel().initialize();
        //call function with true to get a first initialisation of the pressure field
        problem_.pressureModel().initialize();

        //update constitutive functions
        problem_.pressureModel().updateMaterialLaws();

        Dune::dinfo << "IMPET: initialization done." << std::endl;

        return;
    }

    //! Calculate the update.
    /*!
     *  \param  t         time
     *  \param dt         time step size
     *  \param updateVec  vector for the update values
     *
     *  Calculates the new pressure and velocity and determines the time step size
     *  and the update of the transported quantity for the explicit time step.
     */
    void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec)
    {
        if (iterFlag_) // only needed if iteration has to be done
        {
            bool converg = false;
            int iter = 0;
            int iterTot = 0;

            // the method is valid for any transported quantity.
            TransportSolutionType transValueOldIter;
            problem_.transportModel().getTransportedQuantity(transValueOldIter);
            TransportSolutionType updateOldIter(transValueOldIter);
            updateOldIter = 0;
            TransportSolutionType transportedQuantity(transValueOldIter);
            TransportSolutionType updateHelp(transValueOldIter);
            TransportSolutionType updateDiff(transValueOldIter);

            while (!converg)
            {
                iter++;
                iterTot++;

                problem_.pressureModel().update();

                //calculate defect of transported quantity
                problem_.transportModel().update(t, dt, updateVec, true);

                updateHelp = updateVec;
                problem_.transportModel().getTransportedQuantity(transportedQuantity);
                transportedQuantity += (updateHelp *= (dt * cFLFactor_));
                transportedQuantity *= omega_;
                transValueOldIter *= (1 - omega_);
                transportedQuantity += transValueOldIter;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                transValueOldIter = transportedQuantity;
                updateOldIter = updateVec;
                Dune::dinfo << " defect = " << dt * updateDiff.two_norm() / transportedQuantity.two_norm();
                // break criteria for iteration loop
                if (iterFlag_ == 2 && dt * updateDiff.two_norm() / transportedQuantity.two_norm() <= maxDefect_)
                {
                    converg = true;
                }
                else if (iterFlag_ == 1 && iter > nIter_)
                {
                    converg = true;
                }

                if (iterFlag_ == 2 && transportedQuantity.infinity_norm() > (1 + maxDefect_))
                {
                    converg = false;
                }
                if (!converg && iter > nIter_)
                {
                    converg = true;
                    std::cout << "Nonlinear loop in IMPET.update exceeded nIter = " << nIter_ << " iterations."
                            << std::endl;
                    std::cout << transportedQuantity.infinity_norm() << std::endl;
                }
            }
            // outputs
            if (iterFlag_ == 2)
                std::cout << "Iteration steps: " << iterTot << std::endl;
            std::cout.setf(std::ios::scientific, std::ios::floatfield);
        }
        else
        {
            //update pressure
            problem_.pressureModel().update();

            //calculate defect of transported quantity
            problem_.transportModel().update(t, dt, updateVec, true);
        }

        dt *= cFLFactor_;
    }

    //! \brief Write data files
    /*!
     *  calls the output methods of both models, pressure and transport.
     *  \param writer the current VTKwriter
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem_.pressureModel().addOutputVtkFields(writer);
        problem_.transportModel().addOutputVtkFields(writer);

        return;
    }

    //! Constructs an IMPET object
    /*!
     * \param prob Problem
     */
    IMPET(Problem& prob) :
            problem_(prob)
    {
        cFLFactor_ = GET_PARAM(TypeTag, Scalar, ImpetCflFactor);
        iterFlag_ = GET_PARAM(TypeTag, int, ImpetIterationFlag);
        nIter_ = GET_PARAM(TypeTag, int, ImpetIterationNumber);
        maxDefect_ = GET_PARAM(TypeTag, Scalar, ImpetMaximumDefect);
        omega_ = GET_PARAM(TypeTag, Scalar, ImpetRelaxationFactor);
    }

private:
    Problem& problem_; //!< object of type Problem including problem
    Scalar cFLFactor_;
    int iterFlag_;
    int nIter_;
    Scalar maxDefect_;
    Scalar omega_;
};
}
#endif
