#! /usr/bin/env python2.7
from init_base import *
from pycli.ode.gridop import *

class InitShallowWater(InitBase):
    id = '1'
    def set_variables(self):
        sigma   = 0.5E6
        x0, y0 = self.xlen * 4/5., self.ylen/2.
        X, Y = meshgrid(self.xaxis, self.yaxis)
        phi = 10000 * (1 + 2 * exp(-((X - x0)**2 + (Y - y0)**2)/(2. * sigma * sigma)))
        model.var   =   {
                'phi'           :   phi.T,
                'west_east'     :   self.xaxis,
                'south_north'   :   self.yaxis,
                'west_eastb'    :   tohalf(self.xaxis, ext = 'both'),
                'south_northb'  :   tohalf(self.yaxis, ext = 'both'),
                'tracer'        :   X.T
        }

if __name__ == '__main__':
    setups  =   {
            'name'    :   'negative beta',
            'nx'      :   400,
            'ny'      :   200,
            'xlen'    :   40.E6,
            'ylen'    :   20.E6,
            'start'   :   0.,
            'end'     :   2880000.,
            'step'    :   600.,
            'frame'   :   2,
            'f0'      :   1.E-4,
            'beta'    :   -1.E-11,
    }
    model = InitShallowWater(setups)
    model.initialize()
