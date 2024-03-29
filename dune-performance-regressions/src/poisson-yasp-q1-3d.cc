//
// Created by abeltluk on 5/12/23.
//
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include "dune/performance-regressions/poisson-problem.hh"

//===============================================================
//===============================================================
// Solve the Poisson equation
//           - \Delta u = f in \Omega,
//                    u = g on \partial\Omega_D
//  -\nabla u \cdot \nu = j on \partial\Omega_N
//===============================================================
//===============================================================

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
    try{
        //Maybe initialize Mpi
        Dune::MPIHelper::instance(argc, argv);

        // YaspGrid Q1 3D test
        {
            // make grid
            Dune::FieldVector<double,3> L(1.0);
            std::array<int,3> N(Dune::filledArray<3,int>(1));
            Dune::YaspGrid<3> grid(L,N);
            grid.globalRefine(3);

            // get view
            typedef Dune::YaspGrid<3>::LeafGridView GV;
            const GV& gv=grid.leafGridView();

            // make finite element map
            typedef GV::Grid::ctype DF;
            typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
            FEM fem(gv);

            // solve problem
            poisson<GV,
                    FEM,
                    Dune::PDELab::ConformingDirichletConstraints,
                    Dune::CGSolver,
                    PreconditionMarker::SeqILU>(gv,fem,"poisson_yasp_Q1_3d",2);
        }

        // test passed
        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}