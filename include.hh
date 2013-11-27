#ifndef INCLUDE
#define INCLUDE

/* define domain decomposition */
#define     NTILES_IN_X     1
#define     NTILES_IN_Y     1
/* define interpolation order */
#define     O_INTERP        2
/* typedef of major variables */
typedef PatchVariable<TileVariable, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> Patch;
typedef std::vector<Patch> StateVector;
typedef ShallowWater<O_INTERP, NTILES_IN_X, NTILES_IN_Y> ForwardModel;


#endif
