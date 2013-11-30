#! /usr/bin/env python2.7
from init_base import *
from pycli.ode.gridop import *

class InitShallowWater(InitBase):
    id = '1'
    def balanced_mean_flow(self):
        sigma   = 0.5E6
        x0, y0 = self.xlen * 4/5., self.ylen/2.
        X, Y = meshgrid(self.xaxis, self.yaxis)
        X, Y = X.T, Y.T
        phi0 = 20000;
        dy  = self.ylen / (self.ny - 1)
        uwind = ones((self.nx + 1, self.ny)) * 10
        uwindx = tohalf(uwind, axis = 0)
        phi = zeros((self.nx, self.ny))
        for i in range(self.ny):
            phi[:, i] = - trapz(uwindx[:, :i], dx = dy, axis = 1) * self.setups['f0'];
        #phi = phi0 - array([0.05/3. * self.setups['beta'] * self.yaxis * self.yaxis * self.yaxis
        #        for x in self.xaxis])
        phi = phi + phi0 + phi0 * 1 * exp(-((X - x0)**2 + (Y - y0)**2)/(2. * sigma * sigma))
        return phi, uwind
    def set_variables(self):
        sigma   = 0.5E6
        x0, y0 = self.xlen * 4/5., self.ylen/2.
        X, Y = meshgrid(self.xaxis, self.yaxis)
        X, Y = X.T, Y.T
        phi = 10000 * (1 + 2 * exp(-((X - x0)**2 + (Y - y0)**2)/(2. * sigma * sigma)))
        #phi, uwind = self.balanced_mean_flow();
        model.var   =   {
                'phi'           :   phi,
                'west_east'     :   self.xaxis,
                'south_north'   :   self.yaxis,
                'west_eastb'    :   self.xaxisb,
                'south_northb'  :   self.yaxisb,
                'tracer'        :   X
                #'uwind'         :   uwind
        }

if __name__ == '__main__':
    setups  =   {
            'name'    :   'negative beta',
            'nx'      :   400,
            'ny'      :   100,
            'xlen'    :   40.E6,
            'ylen'    :   10.E6,
            'start'   :   0.,
            'end'     :   360000.,
            'step'    :   600.,
            'frame'   :   2,
            'f0'      :   1.E-4,
            'beta'    :   1.E-11,
    }
    model = InitShallowWater(setups)
    model.initialize()
