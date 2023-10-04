//
// Created by abeltluk on 5/12/23.
//
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>

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

#if HAVE_ALBERTA
        {
      // make grid
      AlbertaUnitSquare grid;
      grid.globalRefine(8);

      // get view
      typedef AlbertaUnitSquare::LeafGridView GV;
      const GV& gv=grid.leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_Alberta_Pk_2d",q);
    }
#endif

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