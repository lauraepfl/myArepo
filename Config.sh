
#!/bin/bash            # this line only there to enable syntax highlighting in this file

#VORONOI
#CHUNKING
PERIODIC 
COOLING
#MESHRELAX_DENSITY_IN_INPUT		 #input density instead of masses
#--------------------------------------- Basic operation mode of code
REFLECTIVE_Z=2                           # in-/outflow boundary conditions in x direction
#REFLECTIVE_Y=2                           # in-/outflow boundary conditions in y direction
#REFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction
#LONG_Z=0.0
READ_MASS_AS_DENSITY_IN_INPUT
BARROS
#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization
#REFINEMENT_VOLUME_LIMIT       # Limit the volume of cells and the maximum volume difference between neighboring cels
#JEANS_REFINEMENT              # Refinement criterion to ensure Jeans stability of cells
#REFINEMENT_HIGH_RES_GAS
#VORONOI_DYNAMIC_UPDATE	
#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)
#----------------------------------------

#--------------------------------------- Magnetohydrodynamics
MHD                           # Master switch for magnetohydrodynamics
MHD_POWELL                    
EXTERNALGRAVITY
MHD_SEEDFIELD
#---------------------------------------
#RIEMANN_HLLC
RIEMANN_HLLD
#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_CENTER_OF_MASS                    # output centers of cells
OUTPUT_PRESSURE
OUTPUT_VOLUME
#OUTPUTCOOLRATE
#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
#DEBUG                                    # enables core-dumps, should this be standard?

#####
#### a seconda della massa dissolve o merge le celle 
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS
MESHRELAX