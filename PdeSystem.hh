#ifndef PDESYSTEM
#define PDESYSTEM
#include <iostream>
#include <iomanip>
#include <ads/timer.h>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include "netcdf.hh"
#include "Include.hh"

struct ncconfig{
    std::string fname;
    long current;
};

class PdeSystem{
public:
	PdeSystem(std::string control_file = "control.in"){
		timer.tic();
		load_domain(control_file);
		load_parameter(control_file);
		ncfile.current = 0;
		ncfile.fname = "dynamics.nc";
        load_nc_file();
	}

	friend std::ostream& operator<< (std::ostream &os, const PdeSystem &other){
		os << "Name : " << other.name << std::endl
			<< "Number of Grids in X: " << other.nrows << std::endl
			<< "Number of Grids in Y: " << other.ncols << std::endl
			<< "Total length in X: " << other.xlen / 1000. << " km" << std::endl
			<< "Total length in Y: " << other.ylen / 1000. << " km" << std::endl
			<< "Grid size in X: " << other.dx / 1000. << " km" << std::endl
			<< "Grid size in Y: " << other.dy / 1000. << " km" <<std::endl
			<< "Time start: " << other.tbegin << " s" << std::endl
			<< "Time end: " << other.tend << " s" << std::endl
			<< "Time step: " << other.dt << " s" << std::endl
			<< "Times per frame: " << other.frame << std::endl;
		os << "Parameters: " << std::endl;
		for (std::map<std::string, float>::const_iterator it = other.sp.begin(); 
				it != other.sp.end(); it++) 
			os << it->first << " = " << it->second << std::endl;
		return os;
	}

	/* define differential equations */
	virtual void operator() (const State &, State &, float t) = 0;
    /* define boundary conditions for variables, needed to update halo */
    virtual void set_boundary_conditions(){}
    /* update diagnostic variables */
	virtual void update(float){}
    /* print out running information on screen */
	virtual void observe(float t){
		long ostep = std::floor(t / dt + 0.5);
		if (ostep % frame != 0) return;
		ncwrite(ostep * dt);
		std::cout 
			<< std::setw(8) << std::left
			<< ostep
			<< std::setw(15) << std::left << std::setprecision(5)
			<< ostep * dt
			<< std::setw(15) << std::left << std::setprecision(5)
			<< timer.toc() << std::endl;
	}
    /* read initial values from nc file */
	void init_variables(){
		for (size_t i = 0; i < attr.size(); i++) var.push_back(gp[attr[i].name]);
	}
	int rows(){ return nrows; }
	int cols(){ return ncols; }
	float start(){ return tbegin; }
	float end(){ return tend; }
	float step(){ return dt; }

protected:
	/* Read time and domain */
	void load_domain(std::string file){
		std::ifstream infile;
		std::string str;

		infile.open(file.c_str(), std::ios::in);
		if (!infile) {ASSERT_FILE_NOT_FOUND(file);}
		else{
			infile >> str; getline(infile, name);
			sclocate(infile, "Time and domain");
			infile
				>> str >> nrows
				>> str >> ncols
				>> str >> xlen
				>> str >> ylen
				>> str >> tbegin
				>> str >> tend
				>> str >> dt
				>> str >> frame;
			infile.close();
			dx = xlen / (nrows - 1);
			dy = ylen / (ncols - 1);
		}
	}

	/* Read Parameter */
	void load_parameter(std::string file){
		std::ifstream infile;
		qstring str;
		float value;

		infile.open(file.c_str(), std::ios::in);
		if (!infile) {ASSERT_FILE_NOT_FOUND(file);}
		else {
			sclocate(infile, "Parameters");
			infile >> str >> value;
			while (!str.empty()){
				sp[str] = value;
				infile >> str >> value; 
			} 
			infile.close();
		}
	}

    /* Read initial condition stored in nc file */
	void load_nc_file(){
		NcFile dataFile(ncfile.fname.c_str(), NcFile::ReadOnly);
		if (!dataFile.is_valid()){
			std::cerr << "Cannot open file: " << ncfile.fname << std::endl;
			exit(-1);
		}
		for (size_t i = 0; i < dataFile.num_vars(); i++){
			NcVar *data = dataFile.get_var(i);
			Grid temp;
			long *edges = data->edges();
			switch (data->num_dims()){
				case 1:
					temp.resize(1, edges[0]);
					data->get(&temp(0, 0), edges[0]);
					break;
				case 2:
					temp.resize(edges[1], edges[0]);
					data->get(&temp(0, 0), edges[0], edges[1]);
					break;
				case 3:
					temp.resize(edges[2], edges[1]);
					data->set_cur(0, 0, 0);
					data->get(&temp(0, 0), 1, edges[1], edges[2]);
					break;
			}
			gp[data->name()] = temp;
		}
	} 

    /* write to nc file */
	virtual void ncwrite(float t){
		NcFile dataFile(ncfile.fname.c_str(),NcFile::Write);
		for (size_t i = 0; i < attr.size(); i++)
			dataFile.get_var(attr[i].name.c_str())->put_rec(&var[i](0, 0), ncfile.current);
		dataFile.get_var("time")->put_rec(&t, ncfile.current);
		ncfile.current++;
	}
    
protected:
	std::string name;
	int nrows, ncols;
	float xlen, ylen;
	float dx, dy;
	float tbegin, tend, dt;
	int frame;
	ncconfig ncfile;
    std::vector<Attribute> attr;
	std::map<std::string, float> sp; // scalar parameter
	std::map<std::string, Grid> gp; // grid parameter
	ads::Timer timer;
};

#endif
