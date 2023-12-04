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

#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include "gridmakers.hh"

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

namespace PreconditionMarker {
    class SeqSSOR;
    class SeqILU;
    class Richardson;
}




template<typename PrecMarker>
struct PreconditionerMaker {
    template <typename MatrixTy,
            typename DomainContainerTy,
            typename RangeContainerTy>
    static auto create(MatrixTy m);
};

template<>
struct __attribute__((feature_variable("SeqSSOR"))) PreconditionerMaker<PreconditionMarker::SeqSSOR> {
    template <typename MatrixTy,
              typename DomainContainerTy,
              typename RangeContainerTy>
    static auto create(MatrixTy m) {
        using Dune::PDELab::Backend::native;
        return Dune::SeqSSOR<typename MatrixTy::Container,
                             DomainContainerTy,
                             RangeContainerTy>(native(m),1,1.0);
    };
};

template<>
struct __attribute__((feature_variable("SeqILU"))) PreconditionerMaker<PreconditionMarker::SeqILU> {
    template <typename MatrixTy,
              typename DomainContainerTy,
              typename RangeContainerTy>
    static auto create(MatrixTy m) {
        using Dune::PDELab::Backend::native;
        return Dune::SeqILU<typename MatrixTy::Container,
                DomainContainerTy,
                RangeContainerTy>(native(m),1.0);
    };
};

template<>
struct __attribute__((feature_variable("Richardson"))) PreconditionerMaker<PreconditionMarker::Richardson> {
    template <typename MatrixTy,
              typename DomainContainerTy,
              typename RangeContainerTy>
    static auto create(MatrixTy m) {
        return Dune::Richardson<DomainContainerTy,
        RangeContainerTy>(1.0);
    };
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
template<typename GridViewTy,
         typename FiniteElementMapTy,
         typename ConstraintsTy,
         template<typename> typename SolverTy = Dune::CGSolver,
         typename PreconditionerTy = PreconditionMarker::SeqILU>
void poisson (const GridViewTy& gv,
              const FiniteElementMapTy& fem,
              std::string filename,
              int q)
{
    using Dune::PDELab::Backend::native;

    // constants and types
    typedef typename FiniteElementMapTy::Traits::FiniteElementType::Traits::
    LocalBasisType::Traits::RangeFieldType R;

    // make function space
    typedef Dune::PDELab::GridFunctionSpace<GridViewTy,FiniteElementMapTy,ConstraintsTy,
    Dune::PDELab::ISTL::VectorBackend<> > GFS;
    GFS gfs(gv,fem);
    gfs.name("solution");

    // make model problem
    typedef PoissonModelProblem<GridViewTy,R> Problem;
    Problem problem;

    // make constraints map and initialize it from a function
    typedef typename GFS::template ConstraintsContainer<R>::Type C;
    C cg;
    cg.clear();
    Dune::PDELab::ConvectionDiffusionBoundaryConditionAdapter<Problem> bctype(gv,problem);
    Dune::PDELab::constraints(bctype,gfs,cg);

    // make local operator
    typedef Dune::PDELab::ConvectionDiffusionFEM<Problem,FiniteElementMapTy> LOP;
    LOP lop(problem);

    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    MBE mbe(27); // 27 is too large / correct for all test cases, so should work fine

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
        std::cout << m1.patternStatistics() << std::endl;
        M m2(m1);
        m2 = 0.0;
        m = m1;
        m = m2;
    }
    gridoperator.jacobian(x0,m);

    // evaluate residual w.r.t initial guess
    typedef typename GridOperator::Traits::Range RV;
    RV r(gfs);
    r = 0.0;
    gridoperator.residual(x0,r);

    // make ISTL solver
    Dune::MatrixAdapter<typename M::Container,typename DV::Container,typename RV::Container> opa(native(m));

    auto prec = PreconditionerMaker<PreconditionerTy>::template create<M,typename DV::Container,typename RV::Container>(m);

    SolverTy<typename DV::Container> solver(opa, prec, 1E-10, 5000, 2);
    Dune::InverseOperatorResult stat;

    // solve the jacobian system
    r *= -1.0; // need -residual
    DV x(gfs,0.0);
    solver.apply(native(x), native(r), stat);
    x += x0;

    // output grid function with VTKWriter
    Dune::VTKWriter<GridViewTy> vtkwriter(gv, Dune::VTK::conforming);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,x);
    vtkwriter.write(filename,Dune::VTK::ascii);
}
