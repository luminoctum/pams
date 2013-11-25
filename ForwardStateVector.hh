#ifndef Forward
#define Forward
#include "ShallowWater.hh"
#include <vector>

typedef PatchVariable<Grid, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> MainPatch;
typedef PatchVariable<StagGridx, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> StagPatchx;
typedef PatchVariable<StagGridy, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> StagPatchy;

template<template<int> class System>
class ForwardStateVector : public System<O_INTERP>{
public:
    ForwardStateVector() : System<O_INTERP>(){}
    void operator()(StagPatchx &uflux, StagPatchy &vflux, MainPatch &phi){
        #pragma omp parallel for
        for (int i = 0; i < NTILES_IN_X * NTILES_IN_Y; i++){
            this->advection(uflux(i), vflux(i), phi(i));
            //this->self_advection(phi(i), uflux(i), vflux(i));
            this->gradient(phi(i), uflux(i), vflux(i));
            this->coriolis(this->f, uflux(i), vflux(i));
        }
        uflux.update(1.);
        vflux.update(1.);
        phi.update(1.);
        uflux.setLeftRightZero();
        vflux.setBottomTopZero();
    }
};

#endif
