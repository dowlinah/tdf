#include <dolfin.h>
#include <iostream>
using namespace dolfin;
#include "forms.h"

#include "TLV_file.hpp"

class u_D_expr : public Expression {
    public:
        void eval(
            Array<double>& values,
            const Array<double>& x
        ) const {
            values[0] = 1+(x[0]*x[0]) - 2*(x[1]*x[1]);
        }
};

class Boundary : public SubDomain {
public:
    bool inside(
        const Array<double>& x, 
        bool on_boundary
    ) const override {
        return on_boundary;
    }
};  

int main(int argc, char** argv) {
    dolfin::init(argc,argv);

    auto mesh = std::make_shared<UnitSquareMesh>(128,128);
    auto V = std::make_shared<forms::FunctionSpace>(mesh);

    auto u_D = std::make_shared<u_D_expr>();
    auto boundary = std::make_shared<Boundary>();
    DirichletBC bc(V, u_D, boundary);

    auto f = std::make_shared<Constant>(-6.0);
    forms::LinearForm L(V);
    L.f = f;
    forms::BilinearForm a(V,V);

    Function u(V);
    solve( a==L, u, bc);

    std::string ofname = "solution.h5";
    auto ofh = HDF5File(mesh->mpi_comm(), ofname, "w");
    ofh.write(u,"solution");
    ofh.close();

    TLVFile out_tlv = TLVFile(mesh->mpi_comm(), "solution.tlv", "meta");
    out_tlv.save_metadata = true;
    out_tlv.write(u);

}
