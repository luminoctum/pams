#! /usr/bin/env python2.7
from pylab import *
from pycli.ncio.ncfile import *

class InitBase:
    var_file    = 'variables.lst'
    nc_file     = 'dynamics.nc'
    id          = '0'
    var         = {}
    def __init__(self, setups):
        self.setups = setups
        self.name   = setups['name']
        self.nx     = setups['nx']
        self.ny     = setups['ny']
        self.xlen   = setups['xlen']
        self.ylen   = setups['ylen']
        self.start  = setups['start']
        self.end    = setups['end']
        self.step   = setups['step']
        self.frame  = setups['frame']
        self.xaxis  = linspace(0, self.xlen, self.nx)
        self.yaxis  = linspace(0, self.ylen, self.ny)
    def initialize(self):
        self.set_variables()
        self.write_ncfile()
    def set_variables(self): pass
    def write_ncfile(self):
        varlist = genfromtxt(self.var_file,
                dtype = 'string',
                delimiter = '"',
                skip_header = 1,
                usecols = [1, 3, 5, 7, 9])
        file = ncfile(self.nc_file, 'w')
        dims = {}
        for i in range(varlist.shape[0]):
            if self.id in varlist[i, 2]:
                if varlist[i, 1][0] == 'd':
                    if varlist[i, 0] == 'time':
                        file.add_dim(varlist[i, 0], [], varlist[i, 3], varlist[i, 4])
                    else:
                        file.add_dim(varlist[i, 0], self.var[varlist[i, 0]], varlist[i, 3], varlist[i, 4])
                    dims[varlist[i, 1][1]] = varlist[i, 0]
                else:
                    dim_name, dim_length = (), ()
                    for key in varlist[i, 1]: 
                        dim_name = dim_name + (dims[key],)
                        if dims[key] != 'time':
                            dim_length = dim_length + (len(self.var[dims[key]]),)
                    if varlist[i, 0] in self.var.keys():
                        file.add_var(varlist[i, 0], dim_name, self.var[varlist[i, 0]].T, varlist[i, 3], varlist[i, 4])
                    else:
                        file.add_var(varlist[i, 0], dim_name, zeros(dim_length), varlist[i, 3], varlist[i, 4])
        file.time[:] = 0
        file.add_atts(self.setups)
