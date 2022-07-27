#include "TLV_file.hpp"

#include <dolfin/io/HDF5Utility.h>

#include <fstream>

void TLVFile::encodeHeader() {
    header.reset();
    header |= (0b01110000 & (stored_type << 4));
    header |= (0b00001111 & label);
    header_set = true;
}  

unsigned TLVFile::get_size(
        const std::string& fname
) {
    std::ifstream is;
    is.open(fname, std::ios::in | std::ios::binary );
    if( !is.good() ) {
        return false;
    }
    char raw_head = 0x0;
    is.read( reinterpret_cast<char*>(&raw_head), sizeof(char) );
    std::int64_t L = -1;
    is.read( reinterpret_cast<char*>(&L), sizeof L );
    return unsigned(L);
}

bool TLVFile::write(const dolfin::Function& u) {
    using namespace dolfin;

    int mpirank;
    MPI_Comm_rank(_mpi_comm.comm(),&mpirank);

    if( save_metadata ) {
        const Mesh& mesh = *u.function_space()->mesh();
        const GenericDofMap& dofmap = *u.function_space()->dofmap();

        const std::size_t tdim = mesh.topology().dim();
        std::vector<dolfin::la_index> cell_dofs;
        std::vector<std::size_t> x_cell_dofs;
        const std::size_t n_cells = mesh.topology().ghost_offset(tdim);

        std::vector<std::size_t> local_to_global_map;
        dofmap.tabulate_local_to_global_dofs(local_to_global_map);

        for(std::size_t i=0; i != n_cells; ++i) {
            x_cell_dofs.push_back(cell_dofs.size());
            auto cell_dofs_i = dofmap.cell_dofs(i);
            for( Eigen::Index j = 0; j < cell_dofs_i.size(); ++j ) {
                auto p = cell_dofs_i[j];
                cell_dofs.push_back(local_to_global_map[p]);
            }
        }

        std::size_t offset = MPI::global_offset(_mpi_comm.comm(), cell_dofs.size(), true);
        std::transform(x_cell_dofs.begin(), x_cell_dofs.end(), x_cell_dofs.begin(),
                std::bind2nd(std::plus<std::size_t>(), offset));

        const bool mpi_io = _mpi_comm.size() > 1 ? true : false;
        std::vector<std::int64_t> global_size(1, MPI::sum(_mpi_comm.comm(), cell_dofs.size()));
        write_mpi_vec(cell_dofs,global_size,metaname+"_cell_dofs");
        
        if(_mpi_comm.rank() == _mpi_comm.size() -1) {
            x_cell_dofs.push_back(global_size[0]);
        }
        global_size[0] = mesh.num_entities_global(tdim)+1;
        write_mpi_vec(x_cell_dofs,global_size,metaname+"_x_cell_dofs");

        std::vector<std::size_t> cells(mesh.topology().global_indices(tdim).begin(),
                mesh.topology().global_indices(tdim).begin() + n_cells);
        global_size[0] = mesh.num_entities_global(tdim);
        write_mpi_vec(cells,global_size,metaname+"_cells");
    }

    // write vector HDF5:158
    std::vector<double> local_data;
    u.vector()->get_local(local_data);
    std::pair<std::size_t, std::size_t> local_range = u.vector()->local_range();
    const std::vector<std::int64_t> global_size_vec(1, u.vector()->size());
    if(mpirank == 0) { // setup blank TLVField in file.
        std::vector<double> initv(2,0);
        set_type(initv);
        std::ofstream ofinitv;
        ofinitv.open(ofname, std::ios::out | std::ios::binary);
        encodeHeader();
        char hea = header.to_ullong();
        ofinitv.write(reinterpret_cast<char*>(&hea),sizeof hea);
        std::int64_t s = global_size_vec[0];
        ofinitv.write(reinterpret_cast<char*>(&s),sizeof s);
        // ofinitv.write(reinterpret_cast<char*>(&initv[0]),s*sizeof(double));
        ofinitv.close();
    }
    MPI_Barrier(_mpi_comm.comm());
    std::pair<std::size_t,std::size_t> range = u.vector()->local_range();
    set_type(local_data);
    write_range(ofname,range.first,(range.second-range.first),local_data);
    
    MPI_Barrier(_mpi_comm.comm());
    return true;
}

bool TLVFile::read(dolfin::Function& u) {
    using namespace dolfin;
    const Mesh& mesh = *u.function_space()->mesh();
    const GenericDofMap& dofmap = *u.function_space()->dofmap();

    const std::size_t num_global_cells = get_size(metaname+"_cells");
    const std::pair<std::size_t, std::size_t> cell_range =
        MPI::local_range(_mpi_comm.comm(),
                num_global_cells);

    std::vector<std::size_t> input_cells;
    std::vector<std::size_t> x_cell_dofs;
    std::vector<dolfin::la_index> input_cell_dofs;
    if( !metadata_cached ) {
        read_range(
            metaname+"_cells",
            cell_range.first,
            cell_range.second,
            input_cells
        );

        read_range(
            metaname+"_x_cell_dofs",
            cell_range.first,
            cell_range.second+1,
            x_cell_dofs
        );

        read_range(
            metaname+"_cell_dofs",
            x_cell_dofs.front(),
            x_cell_dofs.back(),
            input_cell_dofs
        );
    }
    if( cache_metadata && !metadata_cached ) {
        cached_input_cells.resize(input_cells.size());
        for( unsigned i = 0; i < input_cells.size(); ++i ) {
            cached_input_cells[i] = input_cells[i];
        }
        cached_x_cell_dofs.resize(x_cell_dofs.size());
        for( unsigned i = 0; i < x_cell_dofs.size(); ++i ) {
            cached_x_cell_dofs[i] = x_cell_dofs[i];
        }
        cached_input_cell_dofs.resize(input_cell_dofs.size());
        for( unsigned i = 0; i < input_cell_dofs.size(); ++i ) {
            cached_input_cell_dofs[i] = input_cell_dofs[i];
        }
        metadata_cached = true;
    }

    GenericVector& x = *u.vector();
    
    const std::size_t num_global_dofs = get_size(ofname);
    const std::pair<dolfin::la_index, dolfin::la_index> 
        input_vector_range = MPI::local_range(_mpi_comm.comm(), num_global_dofs);

    std::vector<double> input_values;
    read_range(
        ofname,
        input_vector_range.first,
        input_vector_range.second,
        input_values
    );

    if( metadata_cached ) {
        HDF5Utility::set_local_vector_values(_mpi_comm.comm(), x, mesh, cached_input_cells,
                cached_input_cell_dofs, cached_x_cell_dofs, input_values, input_vector_range, dofmap);
    } else { 
        HDF5Utility::set_local_vector_values(_mpi_comm.comm(), x, mesh, input_cells,
                input_cell_dofs, x_cell_dofs, input_values, input_vector_range, dofmap);
    }

}
