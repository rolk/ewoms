// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
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
#ifndef EWOMS_CELLDATA1P_HH
#define EWOMS_CELLDATA1P_HH

#include "1pproperties.hh"
#include "fluxData1p.hh"

/*!
 * \file
 * \copydoc Ewoms::CellData1P
 */

namespace Ewoms
{
template<class TypeTag>
class FluxData1P;

/*!
 * \ingroup OnePhase
 */
//! Class including data of one grid cell.
/*! The variables of one-phase flow, which are the pressure as well as additional data assigned to cell-cell interfaces, so-called flux-data, are stored.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class CellData1P
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef FluxData1P<TypeTag> FluxData;

private:
    Scalar pressure_;
    FluxData fluxData_;

public:

    //! Constructs a CellData1P object
    CellData1P() :
        pressure_(0.0)
    {
    }
    //! Returns the flux data of the cell
    FluxData& fluxData()
    {
        return fluxData_;
    }
    //! Returns the flux data of the cell
    const FluxData& fluxData() const
    {
        return fluxData_;
    }

    ////////////////////////////////////////////////////////////
    // functions returning primary variables
    ////////////////////////////////////////////////////////////

    //!\brief Returns the cell pressure
    Scalar pressure()
    {
        return pressure_;
    }
    //!\brief Returns the cell pressure
    Scalar pressure() const
    {
        return pressure_;
    }
    //!Sets the cell pressure
    void setPressure(Scalar press)
    {
        pressure_ = press;
    }
};

}
#endif
