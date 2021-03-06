// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
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
 * \brief Source file for a test of the diffusion model
 */
#include "config.h"
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/gridinfo.hh>

#include "test_diffusionproblem.hh"
#include "resultevaluation.hh"

////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << "usage: " << progname << " #refine [delta]\n";
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(FVVelocity2PTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        static const int dim = Grid::dimension;
        typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

        TimeManager::registerParameters();
        END_PARAM_REGISTRATION;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc != 2 && argc != 3)
            usage(argv[0]);

        int numRefine;
        std::istringstream(argv[1]) >> numRefine;

        double delta = 1e-3;
        if (argc == 3)
            std::istringstream(argv[2]) >> delta;

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        Dune::FieldVector<int,dim> N(1);
        GlobalPosition L(0.0);
        GlobalPosition H(1.0);
        Grid grid(N,L,H);
        grid.globalRefine(numRefine);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        Dune::Timer timer;
        bool consecutiveNumbering = true;

        typedef GET_PROP_TYPE(TTAG(FVVelocity2PTestProblem), Problem) FVProblem;
        FVProblem fvProblem(grid.leafView(), delta);
        fvProblem.setName("fvdiffusion");
        timer.reset();
        fvProblem.init();
        fvProblem.calculateFVVelocity();
        double fvTime = timer.elapsed();
        fvProblem.writeOutput();
        Ewoms::ResultEvaluation fvResult;
        fvResult.evaluate(grid.leafView(), fvProblem, consecutiveNumbering);

        typedef GET_PROP_TYPE(TTAG(FVMPFAOVelocity2PTestProblem), Problem) MPFAOProblem;
        MPFAOProblem mpfaProblem(grid.leafView(), delta);
        mpfaProblem.setName("fvmpfaodiffusion");
        timer.reset();
        mpfaProblem.init();
        double mpfaTime = timer.elapsed();
        mpfaProblem.writeOutput();
        Ewoms::ResultEvaluation mpfaResult;
        mpfaResult.evaluate(grid.leafView(), mpfaProblem, consecutiveNumbering);

        typedef GET_PROP_TYPE(TTAG(MimeticPressure2PTestProblem), Problem) MimeticProblem;
        MimeticProblem mimeticProblem(grid.leafView(), delta);
        mimeticProblem.setName("mimeticdiffusion");
        timer.reset();
        mimeticProblem.init();
        double mimeticTime = timer.elapsed();
        mimeticProblem.writeOutput();
        Ewoms::ResultEvaluation mimeticResult;
        mimeticResult.evaluate(grid.leafView(), mimeticProblem, consecutiveNumbering);

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        std::cout.precision(2);
        std::cout << "\t error press \t error grad\t sumflux\t erflm\t\t uMin\t\t uMax\t\t time" << std::endl;
        std::cout << "2pfa\t " << fvResult.relativeL2Error << "\t " << fvResult.ergrad << "\t " << fvResult.sumflux
                        << "\t " << fvResult.erflm << "\t " << fvResult.uMin << "\t " << fvResult.uMax << "\t " << fvTime << std::endl;
        std::cout << "mpfa-o\t " << mpfaResult.relativeL2Error << "\t " << mpfaResult.ergrad << "\t " << mpfaResult.sumflux
                        << "\t " << mpfaResult.erflm << "\t " << mpfaResult.uMin << "\t " << mpfaResult.uMax << "\t " << mpfaTime << std::endl;
        std::cout << "mimetic\t " << mimeticResult.relativeL2Error << "\t " << mimeticResult.ergrad << "\t " << mimeticResult.sumflux
                        << "\t " << mimeticResult.erflm << "\t " << mimeticResult.uMin << "\t " << mimeticResult.uMax << "\t " << mimeticTime << std::endl;



        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
