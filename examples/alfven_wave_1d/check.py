""" @package ./examples/alfven_wave_1d/check.py
Code that checks results of 1d Alfven wave propagation problem

created by Alessandro Stenghel and Federico Marinacci, 
last modified 13.7.2020 -- comments welcome
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os     # file specific calls

simulation_directory = str(sys.argv[1])
print("alfven_wave_1d: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32 # integer type

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

""" open initial conditiions to get parameters """
try:
    data = h5py.File(simulation_directory + "/ics.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = np.int32(data["Header"].attrs["NumPart_Total"][0])

""" maximum L1 error after two propagations """
DeltaMaxAllowed = 1e-4 * FloatType(NumberOfCells)**-2

""" initial state -- copied from create.py """
density_0 = FloatType(1.0)
velocity_0 = FloatType(0.0)
pressure_0 = FloatType(1.0)
gamma = FloatType(5.0) / FloatType(3.0)
gamma_minus_one = gamma - FloatType(1.0)
delta = FloatType(1e-6)    # relative velocity perturbation
uthermal_0 = pressure_0 / density_0 / gamma_minus_one
bfield_0= FloatType(1.0)
k_z = FloatType(2.0*np.pi)
omega = bfield_0*k_z/np.sqrt(density_0) 

""" 
    loop over all output files; need to be at times when analytic
    solution equals the initial conditions
"""
i_file = 0
status = 0
error_data = []
while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    
    """ get simulation data """
    ## simulation data
    time = FloatType(data['Header'].attrs['Time'])
    Pos = np.array(data["PartType0"]["CenterOfMass"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    Bfield = np.array(data["PartType0"]["MagneticField"], dtype = FloatType)/FloatType(np.sqrt(4.0*np.pi))
    Volume = Mass / Density
    Pressure = gamma_minus_one*Density*Uthermal

    """ calculate analytic solution at cell positions """
    Density_ref = np.full(Pos.shape[0], density_0, dtype=FloatType)
    Velocity_ref = np.zeros(Pos.shape, dtype=FloatType)
    Pressure_ref = np.full(Pos.shape[0], pressure_0, dtype=FloatType)
    Bfield_ref = np.zeros(Pos.shape, dtype=FloatType)
    
    ## perturbations
    Velocity_ref[:,0] = velocity_0
    Velocity_ref[:,1] = delta*np.sin(k_z*Pos[:,0]-omega*time)
    Velocity_ref[:,2] = delta*np.cos(k_z*Pos[:,0]-omega*time)
    Bfield_ref[:,0] = bfield_0
    Bfield_ref[:,1] = -k_z*bfield_0/omega*Velocity_ref[:,1]
    Bfield_ref[:,2] = -k_z*bfield_0/omega*Velocity_ref[:,2]

    """ compare data """
    ## density
    abs_delta_dens = np.abs(Density - Density_ref)
    L1_dens = np.average(abs_delta_dens, weights=Volume)
    
    ## velocity 
    abs_delta_vel_y = np.abs(Velocity - Velocity_ref)[:,1]
    L1_vel_y = np.average(abs_delta_vel_y, weights=Volume)

    abs_delta_vel_z = np.abs(Velocity - Velocity_ref)[:,2]
    L1_vel_z = np.average(abs_delta_vel_z, weights=Volume)

    ## magnetic field 
    abs_delta_bfield_y = np.abs(Bfield -Bfield_ref)[:,1]
    L1_bfield_y = np.average(abs_delta_bfield_y, weights=Volume)

    abs_delta_bfield_z = np.abs(Bfield -Bfield_ref)[:,2]
    L1_bfield_z = np.average(abs_delta_bfield_z, weights=Volume)
    
    ## pressure
    abs_delta_pressure = np.abs(Pressure-Pressure_ref)
    L1_pressure = np.average(abs_delta_pressure, weights=Volume)

    """ printing results """
    print("alfven_wave_1d: L1 error of " + filename +":")
    print("\t density: %g" % L1_dens)
    print("\t velocity y: %g" % L1_vel_y)
    print("\t velocity z: %g" % L1_vel_z)
    print("\t magnetic field y: %g" % L1_bfield_y)
    print("\t magnetic field z: %g" % L1_bfield_z)
    print("\t pressure: %g" % L1_pressure)
    print("\t tolerance: %g for %d cells" % (DeltaMaxAllowed, NumberOfCells) )
    
    error_data.append(np.array([time, L1_dens, L1_vel_y, L1_vel_z, \
                      L1_bfield_y, L1_bfield_z, L1_pressure], dtype=FloatType))
    
    
    """ criteria for failing the test """
    if L1_dens > DeltaMaxAllowed or \
       L1_vel_y > DeltaMaxAllowed or L1_vel_z > DeltaMaxAllowed or \
       L1_bfield_y > DeltaMaxAllowed or L1_bfield_z > DeltaMaxAllowed or \
       L1_pressure > DeltaMaxAllowed:
         status = 1
    
    if makeplots and i_file >= 0:
      if not os.path.exists( simulation_directory+"/plots" ):
        os.mkdir( simulation_directory+"/plots" )
      
      # only import matplotlib if needed
      import matplotlib.pyplot as plt

      # plot density
      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      dx = Boxsize / FloatType(Pos.shape[0])
      ax.plot( Pos[:,0], Density_ref , 'k', lw=0.7, label="Analytical solution" )
      ax.plot( Pos[:,0], Density, 'o-r', mec='r', mfc="None", label="Arepo" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Density" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_dens), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/density_%02d.pdf"%(i_file) )
      plt.close(f)
      
      # plot pressure
      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      ax.plot( Pos[:,0], Pressure_ref , 'k', lw=0.7, label="Analytical solution" )
      ax.plot( Pos[:,0], Pressure, 'o-r', mec='r', mfc="None", label="Arepo" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Pressure" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_pressure), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/pressure_%02d.pdf"%(i_file) )
      plt.close(f)

      # plot velocities
      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      ax.plot( Pos[:,0], Velocity_ref[:,1] , ':k', lw=0.7, label="Analytical solution y" )
      ax.plot( Pos[:,0], Velocity[:,1] , 'o-m', mec='m', mfc="None", label="Arepo v y" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Velocity y" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_vel_y), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/velocityy_%02d.pdf"%(i_file) )
      plt.close(f)

      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      ax.plot( Pos[:,0], Velocity_ref[:,2] , ':k', lw=0.7, label="Analytical solution z" )
      ax.plot( Pos[:,0], Velocity[:,2], 'o-c', mec='c', mfc="None", label="Arepo v z" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Velocity z" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_vel_z), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/velocityz_%02d.pdf"%(i_file) )
      plt.close(f)

      #plot Bfields
      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      ax.plot( Pos[:,0], Bfield_ref[:,1] , ':k', lw=0.7, label="Analytical solution y" )
      ax.plot( Pos[:,0], Bfield[:,1], 'o-m', mec='m', mfc="None", label="Arepo B y" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Magnetic Field y" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_bfield_y), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/bfieldy_%02d.pdf"%(i_file) )
      plt.close(f)

      plt.rcParams['text.usetex'] = True
      f = plt.figure( figsize=(3.5,3.5) )
      ax = plt.axes( [0.19, 0.12, 0.75, 0.75] )
      ax.plot( Pos[:,0], Bfield_ref[:,2] , ':k', lw=0.7, label="Analytical solution z" )
      ax.plot( Pos[:,0], Bfield[:,2], 'o-c', mec='c', mfc="None", label="Arepo B z" )
      ax.set_xlim( 0, Boxsize )
      ax.set_xlabel( "x" )
      ax.set_ylabel( "Magnetic Field z" )
      ax.legend( loc='upper right', frameon=False, fontsize=8 )
      ax.set_title( "$\mathrm{alfven\_wave\_1d:}\ \mathrm{N}=%d,\ \mathrm{L1}=%4.1e$" % (NumberOfCells,L1_bfield_z), loc='right', size=8 )
      plt.ticklabel_format( axis='y', style='sci', scilimits=(0,0) )
      f.savefig( simulation_directory+"plots/bfieldz_%02d.pdf"%(i_file) )
      plt.close(f)
    
    i_file += 1

#save L1 errors
np.savetxt(simulation_directory+"/error_%d.txt"%NumberOfCells, np.array(error_data))

""" normal exit """
sys.exit(status) 
