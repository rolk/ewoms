// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
 * \copydoc Ewoms::FVSpatialParamsOneP
 */
#ifndef EWOMS_FV_SPATIAL_PARAMS_ONE_P_HH
#define EWOMS_FV_SPATIAL_PARAMS_ONE_P_HH

#include <ewoms/decoupled/common/decoupledproperties.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/math.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms
{
// forward declation of property tags
namespace Properties
{
NEW_PROP_TAG( SpatialParams);
}


/*!
 * \ingroup SpatialParams
 */

/*!
 * \brief The base class for spatial parameters of problems using the
 *        fv method.
 */
template<class TypeTag>
class FVSpatialParamsOneP
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) Implementation;

    enum
    {
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

public:
    FVSpatialParamsOneP(const GridView &gv)
    {
    }

    ~FVSpatialParamsOneP()
    {
    }

    /*!
     * \brief Averages the intrinsic permeability (Scalar).
     * \param K1 intrinsic permeability of the first element
     * \param K2 intrinsic permeability of the second element
     */
    Scalar meanK(Scalar K1, Scalar K2) const
    {
        const Scalar K = Ewoms::harmonicMean(K1, K2);
        return K;
    }

    /*!
     * \brief Averages the intrinsic permeability (Scalar).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first element
     * \param K2 intrinsic permeability of the second element
     */
    void meanK(Tensor &result, Scalar K1, Scalar K2) const
    {
        const Scalar K = Ewoms::harmonicMean(K1, K2);
        for (int i = 0; i < dimWorld; ++i)
        {
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = 0;
            result[i][i] = K;
        }
    }

    /*!
     * \brief Averages the intrinsic permeability (Tensor).
     * \param result averaged intrinsic permeability
     * \param K1 intrinsic permeability of the first element
     * \param K2 intrinsic permeability of the second element
     */
    void meanK(Tensor &result, const Tensor &K1, const Tensor &K2) const
    {
        // entry-wise harmonic mean at the main diagonal and arithmetic mean at the off-diagonal
        for (int i = 0; i < dimWorld; ++i)
        {
            result[i][i] = harmonicMean(K1[i][i], K2[i][i]);
            for (int j = 0; j < dimWorld; ++j)
            {
                if (i != j)
                {
                    result[i][j] = 0.5 * (K1[i][j] + K2[i][j]);
                }
            }
        }
    }

    /*!
     * \brief Dummy function that can be used if only one value exist (boundaries).
     * \param result intrinsic permeability
     * \param K intrinsic permeability of the element
     */
    void meanK(Tensor &result, Scalar K) const
    {
        for (int i = 0; i < dimWorld; ++i)
        {
            for (int j = 0; j < dimWorld; ++j)
                result[i][j] = 0;
            result[i][i] = K;
        }
    }

    /*!
     * \brief Dummy function that can be used if only one value exist (boundaries).
     * \param result intrinsic permeability
     * \param K intrinsic permeability of the element
     */
    void meanK(Tensor &result, const Tensor &K) const
    {
        result = K;
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability.
     *
     * \return intrinsic (absolute) permeability
     * \param element The element
     */
    const Tensor& intrinsicPermeability (const Element& element) const
    {
        return asImp_().intrinsicPermeabilityAtPos(element.geometry().center());
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability.
     *
     * \return intrinsic (absolute) permeability
     * \param globalPos The position of the center of the element
     */
    const Tensor& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a intrinsicPermeabilityAtPos() method.");
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \return porosity
     * \param element The element
     */
    Scalar porosity(const Element& element) const
    {
        return asImp_().porosityAtPos(element.geometry().center());
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \return porosity
     * \param globalPos The position of the center of the element
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a porosityAtPos() method.");
    }

private:
    Implementation &asImp_()
    {
        return *static_cast<Implementation*> (this);
    }

    const Implementation &asImp_() const
    {
        return *static_cast<const Implementation*> (this);
    }
};

} // namespace Ewoms

#endif
