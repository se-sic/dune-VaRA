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

        // YaspGrid Q1 2D test
        {
            // make grid
            Dune::FieldVector<double,2> L(1.0);
            std::array<int,2> N(Dune::filledArray<2,int>(1));
            Dune::YaspGrid<2> grid(L,N);
            grid.globalRefine(3);

            // get view
            typedef Dune::YaspGrid<2>::LeafGridView GV;
            const GV& gv=grid.leafGridView();

            // make finite element map
            typedef GV::Grid::ctype DF;
            typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
            FEM fem(gv);

            // solve problem
            poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q1_2d",2);
        }

        // YaspGrid Q2 2D test
        {
            // make grid
            Dune::FieldVector<double,2> L(1.0);
            std::array<int,2> N(Dune::filledArray<2,int>(1));
            Dune::YaspGrid<2> grid(L,N);
            grid.globalRefine(3);

            // get view
            typedef Dune::YaspGrid<2>::LeafGridView GV;
            const GV& gv=grid.leafGridView();

            // make finite element map
            typedef GV::Grid::ctype DF;
            typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
            FEM fem(gv);

            // solve problem
            poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q2_2d",2);
        }

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
            poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q1_3d",2);
        }

        // YaspGrid Q2 3D test
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
            typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
            FEM fem(gv);

            // solve problem
            poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_yasp_Q2_3d",2);
        }

        // UG Pk 2D test
#if HAVE_DUNE_UGGRID
        {
      // make grid
      std::shared_ptr<Dune::UGGrid<2> > grid(TriangulatedUnitSquareMaker<Dune::UGGrid<2> >::create());
      grid->globalRefine(4);

      // get view
      typedef Dune::UGGrid<2>::LeafGridView GV;
      const GV& gv=grid->leafGridView();

      // make finite element map
      typedef GV::Grid::ctype DF;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_UG_Pk_2d",q);
    }
#endif

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

#if HAVE_DUNE_ALUGRID
        {
      using ALUType = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
      auto alugrid = Dune::StructuredGridFactory<ALUType>::createSimplexGrid(Dune::FieldVector<ALUType::ctype, 2>(0.0), Dune::FieldVector<ALUType::ctype, 2>(1.0), std::array<uint,2>{1u, 1u});
      alugrid->globalRefine(4);

      // get view
      using GV = ALUType::LeafGridView;
      auto gv = alugrid->leafGridView();

      // make finite element map
      typedef ALUType::ctype DF;
      const int k=3;
      const int q=2*k;
      typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,k> FEM;
      FEM fem(gv);

      // solve problem
      poisson<GV,FEM,Dune::PDELab::ConformingDirichletConstraints>(gv,fem,"poisson_ALU_Pk_2d",q);
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
