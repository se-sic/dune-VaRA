// -*- tab-width: 4; indent-tabs-mode: nil -*-
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include<dune/grid/yaspgrid.hh>
#include<dune/grid/common/gridfactory.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#include <dune/grid/albertagrid/gridfactory.hh>
#endif
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid/uggridfactory.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

class YaspUnitSquare : public Dune::YaspGrid<2>
{
public:

    YaspUnitSquare () : Dune::YaspGrid<2>(Dune::FieldVector<double,2>(1.0),
                                          {{1,1}},
                                          std::bitset<2>(false),
                                          0)
    {}
};

#if HAVE_DUNE_ALUGRID
class ALUUnitSquare : public Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>
{
public:
  ALUUnitSquare () : Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>(GRIDSDIR "/2dsimplex.alu") {}
};

#endif //HAVE_DUNE_ALUGRID


#if HAVE_ALBERTA
#  if ALBERTA_DIM == 2
class AlbertaLDomain : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaLDomain () : Dune::AlbertaGrid<2,2>(GRIDSDIR "/ldomain.al") {}
};

class AlbertaUnitSquare : public Dune::AlbertaGrid<2,2>
{
public:
  AlbertaUnitSquare () : Dune::AlbertaGrid<2,2>(GRIDSDIR "/2dgrid.al") {}
};

class AlbertaReentrantCorner : public Dune::GridPtr<Dune::AlbertaGrid<2,2> >
{
public:
  AlbertaReentrantCorner()
    : Dune::GridPtr<Dune::AlbertaGrid<2,2> >(GRIDSDIR "/2dreentrantcorner.dgf")
  { }
};
#  endif //ALBERTA_DIM == 2
#endif //HAVE_ALBERTA

template<typename Grid>
class TriangulatedLDomainMaker {
    static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
    static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
    static std::unique_ptr<Grid> create() {
        Dune::GridFactory<Grid> gf;

        Dune::FieldVector<typename Grid::ctype, 2> pos;
        pos[0] =-1.0;  pos[1] =-1.0; gf.insertVertex(pos);
        pos[0] = 0.0;  pos[1] =-1.0; gf.insertVertex(pos);
        pos[0] =-1.0;  pos[1] = 0.0; gf.insertVertex(pos);
        pos[0] = 0.0;  pos[1] = 0.0; gf.insertVertex(pos);
        pos[0] = 1.0;  pos[1] = 0.0; gf.insertVertex(pos);
        pos[0] =-1.0;  pos[1] = 1.0; gf.insertVertex(pos);
        pos[0] = 0.0;  pos[1] = 1.0; gf.insertVertex(pos);
        pos[0] = 1.0;  pos[1] = 1.0; gf.insertVertex(pos);

        auto type = Dune::GeometryTypes::triangle;
        std::vector<unsigned int> vid(3);
        vid[0] = 0;  vid[1] = 1;  vid[2] = 2; gf.insertElement(type, vid);
        vid[0] = 2;  vid[1] = 1;  vid[2] = 3; gf.insertElement(type, vid);
        vid[0] = 2;  vid[1] = 3;  vid[2] = 5; gf.insertElement(type, vid);
        vid[0] = 5;  vid[1] = 3;  vid[2] = 6; gf.insertElement(type, vid);
        vid[0] = 3;  vid[1] = 4;  vid[2] = 6; gf.insertElement(type, vid);
        vid[0] = 6;  vid[1] = 4;  vid[2] = 7; gf.insertElement(type, vid);

        return std::unique_ptr<Grid>(gf.createGrid());
    }
};

//////////////////////////////////////////////////////////////////////
//
// UnitTriangle
//

template<typename Grid>
class UnitTriangleMaker {
    static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
    static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
    static std::unique_ptr<Grid> create() {
        Dune::GridFactory<Grid> gf;
        Dune::FieldVector<typename Grid::ctype, 2> pos;

        pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);

        auto type = Dune::GeometryTypes::triangle;
        std::vector<unsigned int> vid(3);

        vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);

        return std::unique_ptr<Grid>(gf.createGrid());
    }
};

#if HAVE_DUNE_ALUGRID
template<>
class UnitTriangleMaker<Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> > {
  typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> Grid;
public:
  static std::unique_ptr<Grid> create() {
    return std::unique_ptr<Grid>(new Grid(GRIDSDIR "/2dtriangle.alu"));
  }
};
#endif // HAVE_DUNE_ALUGRID

//////////////////////////////////////////////////////////////////////
//
// TriangulatedUnitSquare
//

template<typename Grid>
class TriangulatedUnitSquareMaker {
    static_assert(Grid::dimension == 2, "Dimension of grid must be 2");
    static_assert(Grid::dimensionworld == 2, "Dimension of world must be 2");
public:
    static std::unique_ptr<Grid> create() {
        Dune::GridFactory<Grid> gf;
        Dune::FieldVector<typename Grid::ctype, 2> pos;

        pos[0] = 0; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 1; gf.insertVertex(pos);

        auto type = Dune::GeometryTypes::triangle;
        std::vector<unsigned int> vid(3);

        vid[0] = 0; vid[1] = 1; vid[2] = 2; gf.insertElement(type, vid);
        vid[0] = 1; vid[1] = 2; vid[2] = 3; gf.insertElement(type, vid);

        return std::unique_ptr<Grid>(gf.createGrid());
    }
};

#if HAVE_DUNE_ALUGRID
template<>
class TriangulatedUnitSquareMaker<Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> > {
  typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> Grid;
public:
  static std::unique_ptr<Grid> create() {
    return std::unique_ptr<Grid>(new Grid(GRIDSDIR "/2dsimplex.alu"));
  }
};
#endif //HAVE_DUNE_ALUGRID

//////////////////////////////////////////////////////////////////////
//
// UnitTetrahedron
//

template<typename Grid>
class UnitTetrahedronMaker {
    static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
    static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
    static std::unique_ptr<Grid> create() {
        Dune::GridFactory<Grid> gf;
        Dune::FieldVector<typename Grid::ctype, 3> pos;

        pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);

        auto type = Dune::GeometryTypes::tetrahedron;
        std::vector<unsigned int> vid(4);

        vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 3; gf.insertElement(type, vid);

        return std::unique_ptr<Grid>(gf.createGrid());
    }
};

//////////////////////////////////////////////////////////////////////
//
// TriangulatedUnitCube
//

// Minimal triangulation with 5 tets, does contain unit tet
// AlbertaSimplexGrid<3,3> cannot refine this, see Flyspry#569
template<typename Grid>
class TriangulatedUnitCubeMaker {
    static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
    static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
    static std::unique_ptr<Grid> create() {
        Dune::GridFactory<Grid> gf;
        Dune::FieldVector<typename Grid::ctype, 3> pos;

        pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
        pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
        pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

        auto type = Dune::GeometryTypes::tetrahedron;
        std::vector<unsigned int> vid(4);

        // tet at vertex 0
        vid[0] = 0; vid[1] = 1; vid[2] = 2; vid[3] = 4; gf.insertElement(type, vid);
        // tet at vertex 3
        vid[0] = 1; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
        // central tet
        vid[0] = 1; vid[1] = 2; vid[2] = 4; vid[3] = 7; gf.insertElement(type, vid);
        // tet at vertex 5
        vid[0] = 1; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
        // tet at vertex 6
        vid[0] = 2; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);

        return std::unique_ptr<Grid>(gf.createGrid());
    }
};

#if HAVE_ALBERTA
# if ALBERTA_DIM == 3
#  ifndef ALLOW_ALBERTA_MINIMAL_TRIANGULATED_CUBE
// AlbertaSimplexGrid<3,3> cannot refine the minimal triangulated cube, see
// Flyspry#569.  If you want to use it nevertheless, define
// ALLOW_ALBERTA_MINIMAL_TRIANGULATED_CUBE before including gridexamples.hh.
// specialize the template to make any attempt to use the create() method fail.
template<>
class TriangulatedUnitCubeMaker<Dune::AlbertaGrid<3,3> >
{};
#  endif //ALLOW_ALBERTA_MINIMAL_TRIANGULATED_CUBE
# endif //ALBERTA_DIM == 3
#endif //HAVE_ALBERTA

//////////////////////////////////////////////////////////////////////
//
// KuhnTriangulatedUnitCubeMaker
//

// Kuhn triangulation with 6 tets, does not contain unit tet, all tets have
// (0,7) as a common edge
template<typename Grid>
class KuhnTriangulatedUnitCubeMaker {
    static_assert(Grid::dimension == 3, "Dimension of grid must be 3");
    static_assert(Grid::dimensionworld == 3, "Dimension of world must be 3");
public:
    static std::unique_ptr<Grid> create() {
        Dune::GridFactory<Grid> gf;

        int fake_argc = 0;
        char **fake_argv = NULL;

        if(Dune::MPIHelper::instance(fake_argc, fake_argv).rank() == 0) {
            Dune::FieldVector<typename Grid::ctype, 3> pos;

            pos[0] = 0; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
            pos[0] = 1; pos[1] = 0; pos[2] = 0; gf.insertVertex(pos);
            pos[0] = 0; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
            pos[0] = 1; pos[1] = 1; pos[2] = 0; gf.insertVertex(pos);
            pos[0] = 0; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
            pos[0] = 1; pos[1] = 0; pos[2] = 1; gf.insertVertex(pos);
            pos[0] = 0; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);
            pos[0] = 1; pos[1] = 1; pos[2] = 1; gf.insertVertex(pos);

            auto type = Dune::GeometryTypes::tetrahedron;
            std::vector<unsigned int> vid(4);

            vid[0] = 0; vid[1] = 1; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
            vid[0] = 0; vid[1] = 1; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
            vid[0] = 0; vid[1] = 4; vid[2] = 5; vid[3] = 7; gf.insertElement(type, vid);
            vid[0] = 0; vid[1] = 4; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
            vid[0] = 0; vid[1] = 2; vid[2] = 6; vid[3] = 7; gf.insertElement(type, vid);
            vid[0] = 0; vid[1] = 2; vid[2] = 3; vid[3] = 7; gf.insertElement(type, vid);
        }

        std::unique_ptr<Grid> gp(gf.createGrid());
        gp->loadBalance();
        return gp;
    }
};

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================

template<typename GV, typename RF>
class PoissonModelProblem
{
    typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
    typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

    //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
    static constexpr bool permeabilityIsConstantPerCell()
    {
        return true;
    }


    //! tensor diffusion coefficient
    typename Traits::PermTensorType
    A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
        typename Traits::PermTensorType I;
        for (std::size_t i=0; i<Traits::dimDomain; i++)
            for (std::size_t j=0; j<Traits::dimDomain; j++)
                I[i][j] = (i==j) ? 1 : 0;
        return I;
    }

    //! velocity field
    typename Traits::RangeType
    b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
        typename Traits::RangeType v(0.0);
        return v;
    }

    //! sink term
    typename Traits::RangeFieldType
    c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
        return 0.0;
    }

    //! source term
    typename Traits::RangeFieldType
    f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
        const auto& xglobal = e.geometry().global(x);
        if (xglobal[0]>0.25 && xglobal[0]<0.375 && xglobal[1]>0.25 && xglobal[1]<0.375)
            return 50.0;
        else
            return 0.0;
    }

    //! boundary condition type function
    BCType
    bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
        typename Traits::DomainType xglobal = is.geometry().global(x);

        if (xglobal[1]<1E-6 || xglobal[1]>1.0-1E-6)
        {
            return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        }
        if (xglobal[0]>1.0-1E-6 && xglobal[1]>0.5+1E-6)
        {
            return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
        }
        return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    }

    //! Dirichlet boundary condition value
    typename Traits::RangeFieldType
    g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
    {
        typename Traits::DomainType xglobal = e.geometry().global(x);
        xglobal -= 0.5;
        return exp(-xglobal.two_norm2());
    }

    //! Neumann boundary condition
    typename Traits::RangeFieldType
    j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
        typename Traits::DomainType xglobal = is.geometry().global(x);

        if (xglobal[0] > 1.0 - 1E-6 && xglobal[1] > 0.5 + 1E-6) {
            return -5.0;
        } else {
            return 0.0;
        }
    }

    //! outflow boundary condition
    typename Traits::RangeFieldType
    o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
    {
        return 0.0;
    }
};

//===============================================================
// Problem setup and solution
//===============================================================

// generate a P1 function and output it
template<typename GV, typename FEM, typename CON>
void poisson (const GV& gv, const FEM& fem, std::string filename, int q)
{
    using Dune::PDELab::Backend::native;

    // constants and types
    typedef typename FEM::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

    // make function space
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,
    Dune::PDELab::ISTL::VectorBackend<> > GFS;
    GFS gfs(gv,fem);
    gfs.name("solution");

    // make model problem
    typedef PoissonModelProblem<GV,R> Problem;
    Problem problem;

    // make constraints map and initialize it from a function
    typedef typename GFS::template ConstraintsContainer<R>::Type C;
    C cg;
    cg.clear();
    Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gv,problem);
    Dune::PDELab::constraints(bctype,gfs,cg);

    // make local operator
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FEM> LOP;
    LOP lop(problem);

#ifdef OLD_BACKEND
    typedef Dune::PDELab::ISTLMatrixBackend MBE;
  MBE mbe;
#else
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    MBE mbe(27); // 27 is too large / correct for all test cases, so should work fine
#endif

    // make grid operator
    typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,
            MBE,
            double,double,double,
            C,C> GridOperator;
    GridOperator gridoperator(gfs,cg,gfs,cg,lop,mbe);

    // make coefficent Vector and initialize it from a function
    // There is some weird shuffling around here - please leave it in,
    // it's there to test the copy constructor and assignment operator of the
    // matrix wrapper
    typedef typename GridOperator::Traits::Domain DV;
    DV x0(gfs,Dune::PDELab::Backend::unattached_container());
    {
        DV x1(gfs);
        DV x2(x1);
        x2 = 0.0;
        x0 = x1;
        x0 = x2;
    }

    // initialize DOFs from Dirichlet extension
    typedef Dune::PDELab::ConvectionDiffusionDirichletExtensionAdapter<Problem> G;
    G g(gv,problem);
    Dune::PDELab::interpolate(g,gfs,x0);
    Dune::PDELab::set_nonconstrained_dofs(cg,0.0,x0);


    // represent operator as a matrix
    // There is some weird shuffling around here - please leave it in,
    // it's there to test the copy constructor and assignment operator of the
    // matrix wrapper
    typedef typename GridOperator::Traits::Jacobian M;
    M m;
    {
        Dune::Timer patternTimer;
        M m1(gridoperator);
        std::cout << "pattern creation:" << patternTimer.elapsed() << std::endl;
#ifndef OLD_BACKEND
        std::cout << m1.patternStatistics() << std::endl;
#endif
        M m2(m1);
        m2 = 0.0;
        m = m1;
        m = m2;
    }
    gridoperator.jacobian(x0,m);
    //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

    // evaluate residual w.r.t initial guess
    typedef typename GridOperator::Traits::Range RV;
    RV r(gfs);
    r = 0.0;
    gridoperator.residual(x0,r);

    // make ISTL solver
    Dune::MatrixAdapter<typename M::Container,typename DV::Container,typename RV::Container> opa(native(m));
    //ISTLOnTheFlyOperator opb(gridoperator);
    Dune::SeqSSOR<typename M::Container,typename DV::Container,typename RV::Container> ssor(native(m),1,1.0);
    Dune::SeqILU<typename M::Container,typename DV::Container,typename RV::Container> ilu0(native(m),1.0);
    Dune::Richardson<typename DV::Container,typename RV::Container> richardson(1.0);

//   typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<M,
//     Dune::Amg::FirstDiagonal> > Criterion;
//   typedef Dune::SeqSSOR<M,V,V> Smoother;
//   typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
//   SmootherArgs smootherArgs;
//   smootherArgs.iterations = 2;
//   int maxlevel = 20, coarsenTarget = 100;
//   Criterion criterion(maxlevel, coarsenTarget);
//   criterion.setMaxDistance(2);
//   typedef Dune::Amg::AMG<Dune::MatrixAdapter<M,V,V>,V,Smoother> AMG;
//   AMG amg(opa,criterion,smootherArgs,1,1);

    Dune::CGSolver<typename DV::Container> solvera(opa,ilu0,1E-10,5000,2);
    // FIXME: Use ISTLOnTheFlyOperator in the second solver again
    Dune::CGSolver<typename DV::Container> solverb(opa,richardson,1E-10,5000,2);
    Dune::InverseOperatorResult stat;

    // solve the jacobian system
    r *= -1.0; // need -residual
    DV x(gfs,0.0);
    solvera.apply(native(x),native(r),stat);
    x += x0;

    // output grid function with VTKWriter
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
    vtkwriter.write(filename,Dune::VTK::ascii);
}
