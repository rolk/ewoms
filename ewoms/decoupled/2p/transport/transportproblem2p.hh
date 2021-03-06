// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
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
 * \copydoc Ewoms::TransportProblem2P
 */
#ifndef EWOMS_TRANSPORTPROBLEM_2P_HH
#define EWOMS_TRANSPORTPROBLEM_2P_HH

#include "transportproperties2p.hh"
#include <ewoms/decoupled/common/onemodelproblem.hh>
#include <ewoms/decoupled/common/fv/fvvelocitydefault.hh>

#include <dune/common/fvector.hh>

namespace Ewoms
{
namespace Properties
{
// Set the model properties
SET_PROP(TransportTwoP, Model)
{
    typedef typename GET_PROP_TYPE(TypeTag, TransportModel) type;
};
SET_TYPE_PROP(TransportTwoP, Velocity, FVVelocityDefault<TypeTag>); //this Property should be set by the pressure model, only for a pure transport it is set here for the transportproblem!
}

/*!
 * \ingroup Saturation2p
 * \ingroup IMPETproblems
 * \brief  Base class for a decoupled two-phase transport problem
 *
 * \tparam TypeTag The problem Type Tag
 */
template<class TypeTag>
class TransportProblem2P : public OneModelProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef OneModelProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution Solution;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;


    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld>      GlobalPosition;

    // private!! copy constructor
    TransportProblem2P(const TransportProblem2P&)
    {}

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TransportProblem2P(TimeManager& timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
        gravity_(0)
    {
        cFLFactor_ = GET_PARAM(TypeTag, Scalar, ImpetCflFactor);

        newSpatialParams_ = true;
        spatialParams_ = new SpatialParams(gridView);

        gravity_ = 0;
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim - 1] = - 9.81;
    }

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param spatialParams SpatialParams instantiation
     */
    TransportProblem2P(TimeManager &timeManager, const GridView &gridView, SpatialParams &spatialParams)
        : ParentType(timeManager, gridView),
        gravity_(0),spatialParams_(spatialParams)
    {
        cFLFactor_ = GET_PARAM(TypeTag, Scalar, ImpetCflFactor);

        newSpatialParams_ = false;

        gravity_ = 0;
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim - 1] = - 9.81;
    }

    ~TransportProblem2P()
    {
        if (newSpatialParams_)
        {
        delete spatialParams_;
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     *
     */
    Scalar temperature(const Element& element) const
    {
        return this->asImp_().temperatureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position of the center of an element
     *
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a temperatureAtPos() method.");
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param element The element
     *
     */
    Scalar referencePressure(const Element& element) const
    {
        return this->asImp_().referencePressureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param globalPos The position of the center of an element
     *
     */
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a referencePressureAtPos() method.");
    }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return *spatialParams_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

    void timeIntegration()
    {
        // allocate temporary vectors for the updates
        Solution updateVector;

        Scalar t = this->timeManager().time();
        Scalar dt = 1e100;

        // obtain the first update and the time step size
        this->model().update(t, dt, updateVector);

        //make sure t_old + dt is not larger than tend
        dt = std::min(dt*cFLFactor_, this->timeManager().episodeMaxTimeStepSize());
        this->timeManager().setTimeStepSize(dt);

        // explicit Euler: Sat <- Sat + dt*N(Sat)
        this->model().updateTransportedQuantity(updateVector);
    }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Ewoms::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;

    // material properties
    SpatialParams*  spatialParams_;
    bool newSpatialParams_;

    Scalar cFLFactor_;
};

}

#endif
