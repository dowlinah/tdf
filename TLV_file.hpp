#ifndef _TLV_FILE_HPP_
#define _TLV_FILE_HPP_

#include <dolfin.h>
#include <iostream>
#include <fstream>

class TLVFile {
    public:
        TLVFile(MPI_Comm comm, const std::string& name, const std::string& mname="") :
            _mpi_comm(comm)
        {
            ofname = name;
            metaname = mname;
            save_metadata = false;
            cache_metadata = false;
            metadata_cached = false;
        }

        bool write(const dolfin::Function& u);
        bool read(dolfin::Function& u);
        bool save_metadata;
        bool metadata_cached, cache_metadata;

    private:
        std::vector<std::size_t> cached_input_cells;
        std::vector<std::size_t> cached_x_cell_dofs;
        std::vector<dolfin::la_index> cached_input_cell_dofs;

        std::string ofname;
        std::string metaname;
        dolfin::MPI::Comm _mpi_comm;
        std::bitset<8> header;
        unsigned stored_type;
        unsigned label;
        bool header_set;

        std::vector<dolfin::la_index> cell_dofs;
        std::vector<std::size_t> x_cell_dofs;

        void encodeHeader();
        void set_type(std::vector<unsigned>& x)    {stored_type = 0;}
        void set_type(std::vector<int>& x)         {stored_type = 1;}
        void set_type(std::vector<float>& x)       {stored_type = 2;}
        void set_type(std::vector<double>& x)      {stored_type = 3;}
        void set_type(std::vector<char>& x)        {stored_type = 4;}
        void set_type(std::vector<std::size_t>& x) {stored_type = 5;}
        unsigned get_size(const std::string& fname);

        template<typename T>
        void write_mpi_vec(
                std::vector<T>& dat, 
                std::vector<std::int64_t>& global_size,
                const std::string& fname
        );

        template<class T>
        bool read_range(
                const std::string& fname,
                const unsigned& start,
                const unsigned& end,
                std::vector<T>& out
        );

        template<class T>
        bool write_range(
                const std::string& fname,
                const unsigned& start,
                const unsigned& end,
                std::vector<T>& out
        );

};
        
template<typename T>
void TLVFile::write_mpi_vec(
        std::vector<T>& dat, 
        std::vector<std::int64_t>& global_size,
        const std::string& fname
) {
    int mpirank;
    MPI_Comm_rank(_mpi_comm.comm(),&mpirank);
    if(mpirank == 0) { // setup blank TLVField in file.
        std::vector<T> initv(1,0);
        set_type(initv);
        std::ofstream ofinit;
        ofinit.open(fname, std::ios::out | std::ios::binary);
        encodeHeader();
        char hea = header.to_ullong();
        ofinit.write(reinterpret_cast<char*>(&hea),sizeof hea);
        std::int64_t s = global_size[0];
        ofinit.write(reinterpret_cast<char*>(&s),sizeof s);
        ofinit.close();
    }
    MPI_Barrier(_mpi_comm.comm());
    std::size_t num_local_items = 1;
    for(std::size_t i= 1; i < global_size.size(); i++ ) {
        num_local_items *= global_size[i];
    }
    num_local_items = dat.size()/num_local_items;
    const std::size_t offset = dolfin::MPI::global_offset(_mpi_comm.comm(), num_local_items, true);

    std::pair<std::size_t,std::size_t> range(offset, offset+num_local_items);
    set_type(dat);
    write_range(fname,range.first,range.second-range.first,dat);
    MPI_Barrier(_mpi_comm.comm());
}

template<class T>
bool TLVFile::read_range(
        const std::string& fname,
        const unsigned& start,
        const unsigned& end,
        std::vector<T>& out
) {
    std::ifstream is;
    is.open(fname, std::ios::in | std::ios::binary );
    if( !is.good() ) {
        return false;
    }
    char raw_head = 0x0;
    is.read( reinterpret_cast<char*>(&raw_head), sizeof(char) );
    header.reset();
    header |= raw_head;
    stored_type = (header>>4).to_ullong();
    label = ((header<<4)>>4).to_ullong();
    header_set = true;
    std::int64_t L = -1;
    is.read( reinterpret_cast<char*>(&L), sizeof L );
    unsigned in_u;
    int in_int;
    float in_float;
    double in_double;
    char in_char;
    std::size_t in_st;
    if(end > L) {
        return false;
    }
    unsigned l_to_read = end-start;
    out.clear();
    out.resize(l_to_read);
    std::size_t sizes[] = {sizeof in_u, sizeof in_int, sizeof in_float, sizeof in_double, sizeof in_char, sizeof in_st};
    unsigned offset=sizeof(char) + sizeof(std::int64_t);
    offset += start*sizes[stored_type];
    is.seekg(offset);
    is.read(reinterpret_cast<char*>(&out[0]),l_to_read*sizes[stored_type]);
    return true;
}

template<class T>
bool TLVFile::write_range(
        const std::string& fname,
        const unsigned& start,
        const unsigned& len,
        std::vector<T>& out
) {
    std::fstream os;
    os.open(fname, std::ios::out | std::ios::in | std::ios::binary);
    std::size_t sizes[] = {sizeof(unsigned), sizeof(int), sizeof(float), sizeof(double), sizeof(char), sizeof(std::size_t)};
    unsigned offset=sizeof(char) + sizeof(std::int64_t);
    offset += start*sizes[stored_type];
    os.seekp(offset);
    for( unsigned i = 0; i < out.size(); ++i ) {
        switch(stored_type) {
            case 0: // unsigned
                os.write(reinterpret_cast<char*>(&out[i]),sizeof(unsigned));
                break;
            case 1: // int
                os.write(reinterpret_cast<char*>(&out[i]),sizeof(int));
                break;
            case 2: // float
                os.write(reinterpret_cast<char*>(&out[i]),sizeof(float));
                break;
            case 3: // double
                os.write(reinterpret_cast<char*>(&out[i]),sizeof(double));
                break;
            case 4: // char
                os.write(reinterpret_cast<char*>(&out[i]),sizeof(char));
                break;
            case 5: // size_t
                os.write(reinterpret_cast<char*>(&out[i]),sizeof(std::size_t));
                break;
        }
    }
    os.close();
    return true;
}

#endif

