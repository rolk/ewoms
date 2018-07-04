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
 * \copydoc Opm::Octane
 */
#ifndef OPM_OCTANE_HPP
#define OPM_OCTANE_HPP

#include <opm/material/components/Component.hpp>
#include <opm/material/IdealGas.hpp>

#include <opm/material/common/Unused.hpp>

namespace Ewoms {
/*!
 * \ingroup Components
 *
 * \brief A simple representation of linear octane
 *
 * \tparam Scalar The type used for scalar values
 */
#warning "TODO: check that the values are consistent with reference CIPR simulator."
template <class Scalar>
class Octane : public Opm::Component<Scalar, Octane<Scalar> >
{
    typedef Opm::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for the iso-octane.
     */
    static const char* name()
    { return "Octane"; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of iso-octane.
     */
    static Scalar molarMass()
    { return 0.11423; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of steam at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        // Assume an ideal gas
        return molarMass()*IdealGas::molarDensity(temperature, pressure);
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& /*temperature*/,
                                   const Evaluation& /*pressure*/)
    {
        return 1e-05; // TODO: this is roughly the viscosity of air at standard conditions
    }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, const Evaluation& density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief Rough estimate of the density of oil \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { return 703.0; }

    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& temperature OPM_UNUSED)
    {
        return 1.47e3; // at 20 degC
    }

    /*!
     * \brief Rough estimate of the viscosity of oil in \f$\mathrm{[Pa*s]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { return 0.005; }

    /*!
     * \brief The enthalpy of iso-octane at a given pressure and temperature \f$\mathrm{[J/kg]}\f$.
     *
     * We simply use the value of iso-octane here.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& pressure)
    {
        return liquidHeatCapacity(temperature, pressure)*temperature; // [J/kg]
    }

    /*!
     * \brief Specific isobaric heat capacity \f$[J/(kg K)]\f$ of liquid iso-octane.
     *
     * We simply use the value of iso-octane here.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature OPM_UNUSED,
                                         const Evaluation& pressure OPM_UNUSED)
    {
        return 240.0/molarMass();
    }

    /*!
     * \brief The enthalpy of iso-octane at a given pressure and temperature \f$\mathrm{[J/kg]}\f$.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                     const Evaluation& pressure)
    {
        return gasHeatCapacity(temperature, pressure)*temperature; // [J/kg]
    }

    /*!
     * \brief Specific isobaric heat capacity \f$[J/(kg K)]\f$ of liquid iso-octane.
     *
     * We simply use the value of iso-octane here.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature OPM_UNUSED,
                                      const Evaluation& pressure OPM_UNUSED)
    {
        return 240/molarMass();
    }


    /*!
     * \brief Specific heat conductivity of liquid TCE \f$\mathrm{[W/(m K)]}\f$.
     *
     * \todo The value returned here is a guess which does not necessarily correspond to reality in any way!
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 0.14; // value for 0 degC
    }
};

} // namespace Opm

#endif
