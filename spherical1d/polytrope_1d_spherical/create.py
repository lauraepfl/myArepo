""" @package examples/polytrope_1d_spherical/create.py
Code that creates 1d polytrope test problem;
supposed to be as siple as possible and self-contained

created by Rainer Weinberger, last modified 04.03.2019
"""


""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.argv[1])
print("examples/polytrope_1d_spherical/create.py: creating ICs in directory " + simulation_directory)


""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(100.0)
NumberOfCells = IntType(256)

## initial state
G = FloatType(1.0)
density_0 = FloatType(1.0)
#PolytropicIndex = FloatType(1.0)
velocity_0 = FloatType(0.0)
gamma = 1.666  ## note: this has to be consistent with the parameter settings for Arepo
Temperature = 1e+4
kb = 1.38e-16
mH = 1.67e-24
pressure = kb*Temperature*density_0 #density in cgs changer
scale_v =  1e+5
scale_l   =  3.086e+18
scale_m    = 1.98e+33
scale_t = scale_l/scale_v
scale_p = scale_m/(scale_l*scale_t**2)
scale_d = 6.807e-23
scale_E = scale_m*scale_v**2
pressure_0 = pressure/scale_p


""" set up grid: uniform radial 1d grid """
## spacing
dx = Boxsize / FloatType(NumberOfCells)
## position of first and last cell
pos_first, pos_last = FloatType(0.1) + FloatType(0.5) * dx, Boxsize - FloatType(0.5) * dx

## set up grid
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = np.linspace(pos_first, pos_last, NumberOfCells, dtype=FloatType)
Volume = FloatType(4.) * np.pi / FloatType(3.) * ( (Pos[:,0]+FloatType(0.5)*dx)**3 - (Pos[:,0]-FloatType(0.5)*dx)**3 )


""" set up hydrodynamical quantitites """
#Rscale = FloatType(0.8) * Boxsize / np.pi
Radius = Pos[:,0]
#theta = np.sinc(Radius/Rscale/np.pi) ## np.sinc(x) = np.sin(np.pi x)/(np.pi x) or 1 if x=0
#theta[theta<1.e-12] = FloatType(1.e-12)
#K = Rscale**2 * FloatType(4.0) * np.pi * G / ( (PolytropicIndex+FloatType(1.0) ) * density_0**(FloatType(1.0)/PolytropicIndex - FloatType(1.0) ) )
## explicit n=1 case
#K = Rscale**2 * FloatType(2.0) * np.pi * G

Density = np.zeros([NumberOfCells], dtype=FloatType)
Pressure = np.zeros([NumberOfCells], dtype=FloatType)


## hydrodynamic quantities
Density.fill(density_0) # * theta**PolytropicIndex
#print(Density)
Velocity = np.zeros([NumberOfCells, 3], dtype=FloatType)
Velocity[:,0] = velocity_0
Pressure.fill(pressure_0) #K * Density**gamma
Uthermal = Pressure / Density / (gamma - FloatType(1.0) )

## density in mass filed in input
Mass = Density
ncells = 0
TotalEnergy = 1e+51
TotalEnergy_0 = TotalEnergy/scale_E
TotVolume = 4./3.*np.pi*((Boxsize/10)**3)

def Volume(i) :
    if i == 0 :
            v = 4./3.*np.pi*(Radius[i]**3)
    else :
        v = 4./3.*np.pi*(Radius[i]**3 - Radius[i-1]**3)

    return v


for i,r in enumerate(Radius) :
    if r < Boxsize/10 :
        ncells += 1

for i,r in enumerate(Radius) :                      #initialize high energy in a small part of space (1/10 of boxsize)
    if r < Boxsize/10 :
       
       Uthermal[i] += TotalEnergy_0/ncells#/Volume(i)

print(Uthermal.sum()*scale_E)#*4./3.*np.pi*(Boxsize)**3)
#quit()

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(FilePath, 'w')    # open/cr#eate file

## create hdf5 groups
header = IC.create_group("Header")    # create header group
part0 = IC.create_group("PartType0")    # create particle group for gas cells

## write header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## write cell data
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Mass)
part0.create_dataset("Velocities", data=Velocity)
part0.create_dataset("InternalEnergy", data=Uthermal)

## close file
IC.close()


""" normal exit """
sys.exit(0) 
