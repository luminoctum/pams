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
public:
    ModuleBase(std::string _fname = "dynamics.nc", long _current = 0){
        current = _current;
        filename = _fname;
        NcFile dataFile(_fname.c_str(), NcFile::ReadOnly);
        if (!dataFile.is_valid()){ ASSERT_FILE_NOT_FOUND(_fname); }
        Eigen::ArrayXXf buffer;
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
                    data->set_cur(_current, 0, 0);
                    data->get(&buffer(0, 0), 1, edges[1], edges[2]);
                    break;
            }
            ncVar[data->name()] = buffer;
        }
        nrows   = dataFile.get_att("nx")->as_int(0);
        ncols   = dataFile.get_att("ny")->as_int(0);
        xlen    = dataFile.get_att("xlen")->as_float(0);
        ylen    = dataFile.get_att("ylen")->as_float(0);
        start   = dataFile.get_att("start")->as_float(0);
        end     = dataFile.get_att("end")->as_float(0);
        step    = dataFile.get_att("step")->as_float(0);
        frame   = dataFile.get_att("frame")->as_int(0);
        dx      = xlen / (nrows - 1);
        dy      = ylen / (ncols - 1);
    }
    template <typename StateType>
    void ncwrite(const StateType &state, float time){
        current++;
        NcFile dataFile(filename.c_str(),NcFile::Write);
        for (size_t i = 0; i < state.size(); i++)
            dataFile.get_var(state[i].name.c_str())->put_rec(&state[i].value(0, 0), current);
        dataFile.get_var("time")->put_rec(&time, current);
    }
};

#endif
