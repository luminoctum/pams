#ifndef MODULEBASE
#define MODULEBASE
#include <map>
#include <Eigen/Dense>
#include "utility.hh"
#include "netcdf.hh"

class ModuleBase{
public:
    std::string filename;
    std::map<std::string, Eigen::ArrayXXf> ncVar;
    long current;
    int nrows, ncols;
    float xlen, ylen;
    float dx, dy;
    float start, end, step;
    int frame;
    Eigen::ArrayXXf buffer;
public:
    ModuleBase(std::string _fname = "dynamics.nc"){
        filename = _fname;
        NcFile dataFile(_fname.c_str(), NcFile::ReadOnly);
        current = dataFile.get_dim("time")->size() - 1;
        if (!dataFile.is_valid()){ ASSERT_FILE_NOT_FOUND(_fname); }
        for (int i = 0; i < dataFile.num_vars(); i++){
            NcVar *data = dataFile.get_var(i);
            long *edges = data->edges();
            switch (data->num_dims()){
                case 1:
                    buffer.resize(1, edges[0]);
                    data->get(&buffer(0, 0), edges[0]);
                    break;
                case 2:
                    buffer.resize(edges[1], edges[0]);
                    data->get(&buffer(0, 0), edges[0], edges[1]);
                    break;
                case 3:
                    buffer.resize(edges[2], edges[1]);
                    data->set_cur(current, 0, 0);
                    data->get(&buffer(0, 0), 1, edges[1], edges[2]);
                    break;
            }
            ncVar[data->name()] = buffer;
        }
        nrows   = dataFile.get_att("nx")->as_int(0);
        ncols   = dataFile.get_att("ny")->as_int(0);
        xlen    = dataFile.get_att("xlen")->as_float(0);
        ylen    = dataFile.get_att("ylen")->as_float(0);
        start   = ncVar["time"](current);
        end     = dataFile.get_att("end")->as_float(0);
        step    = dataFile.get_att("step")->as_float(0);
        frame   = dataFile.get_att("frame")->as_int(0);
        dx      = xlen / (nrows - 1);
        dy      = ylen / (ncols - 1);
    }
    friend std::ostream& operator<< (std::ostream &os, const ModuleBase &other){
        os << "========== Experiment file name : " << other.filename << " ==========" << std::endl
            << "Number of Grids in X: " << other.nrows << std::endl
            << "Number of Grids in Y: " << other.ncols << std::endl
            << "Total length in X: " << other.xlen / 1000. << " km" << std::endl
            << "Total length in Y: " << other.ylen / 1000. << " km" << std::endl
            << "Grid size in X: " << other.dx / 1000. << " km" << std::endl
            << "Grid size in Y: " << other.dy / 1000. << " km" <<std::endl
            << "Time start: " << other.start << " s" << std::endl
            << "Time end: " << other.end << " s" << std::endl
            << "Time step: " << other.step << " s" << std::endl
            << "Times per frame: " << other.frame << std::endl;
        return os;
    }
    template <typename StateType>
    void ncwrite(const StateType &state, float time){
        current++;
        NcFile dataFile(filename.c_str(),NcFile::Write);
        for (size_t i = 0; i < state.size(); i++){
            if (state[i].scale_by_mass){
                buffer = state[i].main() / state[0].main_stag(state[i].stag);
            } else{
                buffer = state[i].main();
            }
            dataFile.get_var(state[i].name.c_str())->put_rec(&buffer(0, 0), current);
            //dataFile.get_var(state[i].name.c_str())->put_rec(&state[i].value(0, 0), current);
        }
        dataFile.get_var("time")->put_rec(&time, current);
    }
};

#endif
