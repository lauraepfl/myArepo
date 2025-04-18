""" @package ./examples/polytrope_1d_spherical/check.py
Code that checks results of 1d polytrope test problem

created by Rainer Weinberger, last modified 04.03.2019
"""

""" load libraries """
import sys    # needed for exit codes
import numpy as np    # scientific computing package
import h5py    # hdf5 format
import os      # file specific calls
import matplotlib.pyplot as plt    ## needs to be active for plotting!
#plt.rcParams['text.usetex'] = True

makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.argv[1])
print("examples/polytrope_1d_spherical/check.py: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(-1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = IntType(data["Header"].attrs["NumPart_Total"][0])

## maximum L1 error after one propagation; based on experience
DeltaMaxAllowed = 0.001 * FloatType(NumberOfCells/256.0)**-2


""" Initial conditions as reference """
i_snap = 0
directory = simulation_directory+"/output/"
filename = "snap_%03d.hdf5" % (i_snap)
#data = h5py.File(directory+filename, "r")
#Pos_ref = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)[:,0]
#Density_ref = np.array(data["PartType0"]["Density"], dtype = FloatType)
#Velocity_ref = np.array(data["PartType0"]["Velocities"], dtype = FloatType)[:,0]
Uthermal_0 = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)

#Accel_ref = np.array(data["PartType0"]["Acceleration"], dtype = FloatType)[:,0]
#GradPress_ref = np.array(data["PartType0"]["PressureGradient"], dtype = FloatType)[:,0]
UnitTime = 0.977 #in MY
gamma = 1.666  ## note: this has to be consistent with the parameter settings for Arepo
Temperature = 1e+4
kb = 1.38e-16
mH = 1.67e-24
scale_v =  1e+5
scale_l   =  3.086e+18
scale_m    = 1.98e+33
scale_t = scale_l/scale_v
scale_p = scale_m/(scale_l*scale_t**2)
scale_d = 6.807e-23
scale_E = scale_m*scale_v**2
CGSTime = UnitTime*31.6e+12
CGSMass = 1.988e+33
CGSlenght = 3.086e+18
CGSVelocity = CGSlenght/CGSTime
ergs = 1e-7
meters = 1e-3
kg = 1e-3

#print(Uthermal_0.sum()*scale_E)
#quit()



time = []
peak_positions = []
peak_Rad = []
time_step = 0.01  # Remplace par ton vrai pas de temps
T=1e+4
kb=1.38e-23
RadialMomentum = np.empty((0,i_snap))
totRadMom = np.zeros(20)
totUtherm = np.zeros(20)
totUkin = np.zeros(20)
totUtot = np.zeros(20)
Uthermal_all = np.empty((0,i_snap))
Ukinetic_all = np.empty((0,i_snap))
Utot_all = (Uthermal_all + Ukinetic_all)
Ufond =  kb*T #Uthermal_ref[Uthermal_ref.size-1]
eta = 2.026
E = 1e+51
density_0 = 1.0*scale_d #repris de create.py QUOI METTRE???
rAna = np.zeros(20)
MassAll = np.zeros(20)
vAna = np.zeros(20)
RadAna = np.zeros(20)




while True:
    i_snap +=1
    """ compare to snapshots """
    filename = "snap_%03d.hdf5" % (i_snap)
    print(f"Trying to open: {directory+filename}")


    try:
        data = h5py.File(directory+filename, "r")
    except:
        break

    Pos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)[:,0]
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)[:,0]
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    Volume = []
    mu = 14/23*mH
    Temperature = 2/3*mu*mH*Uthermal/kb/Density

    for i,r in enumerate(Pos) :
      if i == 0 :
        Volume.append(4./3.*np.pi*(Pos[i]**3))
      else :
        Volume.append(4./3.*np.pi*(Pos[i]**3 - Pos[i-1]**3))

    n = Density*CGSMass/Volume/CGSlenght**3/mu/mH
    Mass = Density*Volume
    RadMom = Velocity*Mass
    Ukinetic = 0.5*Velocity**2*Mass
    Pressure = n*Temperature*kb #Temperature*kb*Density/mu/mH MAUVAISE VALEUR

    totRadMom[i_snap-1] = RadMom.sum()
    totUtherm[i_snap-1] = (Uthermal.sum() - Ufond)*CGSMass*CGSlenght**2/CGSTime**2
    totUkin[i_snap-1] = Ukinetic.sum()*CGSMass*CGSlenght**2/CGSTime**2
    totUtot[i_snap-1] = (Uthermal.sum() - Ufond + Ukinetic.sum())*CGSMass*CGSlenght**2/CGSTime**2

    t= time_step*i_snap

    ################################# Solution analytique #####################################
    rAna[i_snap-1] = (eta*E/scale_E/density_0*scale_d)**0.2*t #pas cgs
    MassAll[i_snap-1] = Mass.sum()/scale_m
   # if MassAll.size == 0:
    #  MassAll = Mass[np.newaxis, :]  # Crée une matrice 2D avec une seule ligne
    #else:
    #  MassAll = np.vstack((MassAll, Mass))  # Ajoute une nouvelle ligne


    
    
    Accel = np.array(data["PartType0"]["Acceleration"], dtype = FloatType)[:,0]
    GradPress = np.array(data["PartType0"]["PressureGradient"], dtype = FloatType)[:,0]

     # Find peak of density
    peak_index = np.argmax(Density)
    peak_position = Pos[peak_index]

    # Find peak of Radial Momentum
    lastpeak = RadMom[::-1] #inverse la liste
    peakRad_index = np.argmax(lastpeak)
    RadPos = Pos[256 - peakRad_index]

    # time of the snap
    time.append(i_snap * time_step * UnitTime)
    peak_positions.append(peak_position)
    peak_Rad.append(RadPos)


    if RadialMomentum.size == 0:
      RadialMomentum = RadMom[np.newaxis, :]  # Crée une matrice 2D avec une seule ligne
    else:
      RadialMomentum = np.vstack((RadialMomentum, RadMom))  # Ajoute une nouvelle ligne

    if Uthermal_all.size == 0:
      Uthermal_all = Uthermal[np.newaxis, :]  # Crée une matrice 2D avec une seule ligne
    else:
      Uthermal_all = np.vstack((Uthermal_all, Uthermal))  # Ajoute une nouvelle ligne

    if Ukinetic_all.size == 0:
      Ukinetic_all = Ukinetic[np.newaxis, :]  # Crée une matrice 2D avec une seule ligne
    else:
      Ukinetic_all = np.vstack((Ukinetic_all, Ukinetic))  # Ajoute une nouvelle ligne



    
    """ plots """
    if makeplots:
        fig, ax = plt.subplots(ncols=1, nrows=4, sharex=True, figsize=np.array([6.9,6.0]) )
        fig.subplots_adjust(left = 0.13, bottom = 0.09,right = 0.98, top = 0.98)
        
        ax[0].plot(Pos, Density, 'b', label='evolved profile')
        ax[0].plot(Pos, Density, 'r+', label='Arepo cells')
       # ax[0].plot(Pos_ref, Density_ref, 'k--', label='initial profile', lw=0.7)
        ax[0].set_ylabel(r"density")
        #ax[0].fill_between([0.1,0.7],[0.0,0.0],[1.01,1.01], color='k',alpha=0.2)
        #ax[0].set_ylim( 0., 1. )
        ax[0].set_yscale("log")
        ax[0].set_xscale("log")

        
        ax[1].plot(Pos, Temperature, 'b')
       # ax[1].plot(Pos_ref, [1.0]*Pos_ref.shape[0], 'k--', lw=0.7)
        ax[1].set_ylabel(r"Temperature [K]")
        #ax[1].set_ylim([0.99,1.01])
        #ax[1].fill_between([0.1,0.7],[0.99,0.99],[1.01,1.01], color='k',alpha=0.2)
        
        ax[2].plot(Pos, Velocity, 'b')
       # ax[2].plot(Pos_ref, [0.0]*Pos_ref.shape[0], 'k--', lw=0.7)
        ax[2].set_ylabel(r"velocity")
        #ax[2].set_ylim([-0.01,0.01])
        #ax[2].fill_between([0.1,0.7],[-0.01,-0.01],[0.01,0.01], color='k',alpha=0.2)
        
        ax[3].plot(Pos, Pressure, 'b')
        #ax[3].plot(Pos_ref[:-1], Accel_ref[:-1] - GradPress_ref[:-1] / Density_ref[:-1], 'k--', lw=0.7)
       # ax[3].set_ylim([-0.5,0.1])#
        
        ax[3].set_ylabel(r"Pressure")
        ax[3].set_xlabel(r"radius")
       # ax[3].set_xlim([0.0,0.8])
       # ax[3].fill_between([0.1,0.7],[-0.5,-0.5],[0.1,0.1], color='k',alpha=0.2)
        
        ax[0].legend(loc=1, frameon=False)
        
        fig.align_ylabels(ax[:])
        
        if not os.path.exists( simulation_directory+"/plots" ):
          os.mkdir( simulation_directory+"/plots" )
        fig.savefig(simulation_directory+"/plots/profiles_%03d.png"%i_snap)
        plt.close(fig)
    
    
    """ compare to ICs by interpolating to IC positions """
    
  #  delta_dens = np.interp(Pos_ref, Pos, Density) - Density_ref
  #  delta_vel = np.interp(Pos_ref, Pos, Velocity) - Velocity_ref
  #  delta_uthermal = np.interp(Pos_ref, Pos, Uthermal) - Uthermal_ref
    
    ## L1 norm
   # i_check, = np.where(Pos_ref < 0.7)
    
   # L1_dens = np.average( np.abs(delta_dens[i_check]) )
   # L1_vel = np.average( np.abs(delta_vel[i_check]) )
   # L1_uthermal = np.average( np.abs(delta_uthermal[i_check]) )
    
    
    """ printing results """
   # print("examples/polytrope_1d_spherical/check.py: L1 error of " + filename +":")
   # print("\t density: %g" % L1_dens)
   # print("\t velocity: %g" % L1_vel)
   # print("\t specific internal energy: %g" % L1_uthermal)
   # print("\t tolerance: %g for %d cells" % (DeltaMaxAllowed, NumberOfCells) )
    
    """ criteria for failing the test """
  #  if L1_dens > DeltaMaxAllowed or L1_vel > DeltaMaxAllowed or L1_uthermal > DeltaMaxAllowed:
  #      sys.exit(1)
########################Solution Analytique######################
vAna = np.gradient(rAna, t)
RadAna = vAna*MassAll*scale_m
    

####################################### plots ######################################

# front wave evolution
fig, ax = plt.subplots()
ax.plot(time, peak_positions, marker='x', linestyle='-', color='b')
ax.plot(time, rAna, marker='x', linestyle='-', color='r')
ax.set_xlabel("Time [MY]")
ax.set_ylabel("front position [pc]")
ax.grid(True)
fig.savefig(simulation_directory+"/plots/front")

# Radial Momentum evolution
fig, ax = plt.subplots()
for i in range(RadialMomentum.shape[0]):
    ax.plot(Pos, RadialMomentum[i, :], label=f'Time step {i}')

ax.set_xlabel("Position [pc]")
ax.set_ylabel("Radial Momentum []")
#ax.set_title("Case with cooling")
ax.legend()
ax.grid(True)
fig.savefig(simulation_directory+"/plots/RadialMomentum")

# Peak radial momentum evolution
fig, ax = plt.subplots()
ax.plot(time, peak_Rad, marker='x', linestyle='-', color='b')
ax.set_xlabel("Time [MY]")
ax.set_ylabel("Radial Momentum peak position [pc]")
#ax.set_title("Case with cooling")
ax.grid(True)
fig.savefig(simulation_directory+"/plots/RadialPeak")

# Thermal energy evolution
fig, ax = plt.subplots()
for i in range(Uthermal_all.shape[0]):
    ax.plot(Pos, Uthermal_all[i, :], label=f'Time step {i}')

ax.set_xlabel("Position [pc]")
ax.set_ylabel("Thermal Energy [ergs]")
#ax.set_title("Case with cooling")
ax.legend()
ax.grid(True)
fig.savefig(simulation_directory+"/plots/Uthermal")

# Kinetic energy evolution
fig, ax = plt.subplots()
for i in range(Ukinetic_all.shape[0]):
    ax.plot(Pos, Ukinetic_all[i, :], label=f'Time step {i}')

ax.set_xlabel("Position [pc]")
ax.set_ylabel("Kinetic Energy [ergs]")
#ax.set_title("Case with cooling")
ax.legend()
ax.grid(True)
fig.savefig(simulation_directory+"/plots/UKinetic")

# Total energy evolution
fig, ax = plt.subplots()
for i in range(Ukinetic_all.shape[0]):
    ax.plot(Pos, Ukinetic_all[i, :], label=f'Time step {i}')

ax.set_xlabel("Position [pc]")
ax.set_ylabel("Total Energy [ergs]")
#ax.set_title("Case with cooling")
ax.legend()
ax.grid(True)
fig.savefig(simulation_directory+"/plots/Tot_Energy")

#  Total radial momentum
fig, ax = plt.subplots()
ax.plot(time, totRadMom, marker='x', linestyle='-', color='b')
ax.plot(time, RadAna, marker='x', linestyle='-', color='r')
ax.set_xlabel("Time [MY]")
ax.set_ylabel("Total Radial Momentum")
ax.set_xscale("log")
ax.set_yscale("log")
ax.grid(True)
fig.savefig(simulation_directory+"/plots/TotalRad")

#  Total Thermal energy
fig, ax = plt.subplots()
ax.plot(time, totUtherm, marker='x', linestyle='-', color='b')
ax.set_xlabel("Time [MY]")
ax.set_ylabel("Total Thermal Energy")
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.grid(True)
fig.savefig(simulation_directory+"/plots/TotalUth")

#  Total Kinetic Energy
fig, ax = plt.subplots()
ax.plot(time, totUkin, marker='x', linestyle='-', color='b')
ax.set_xlabel("Time [MY]")
ax.set_ylabel("Total Kinetic Energy")
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.grid(True)
fig.savefig(simulation_directory+"/plots/TotalUkin")

#  Total Energy
fig, ax = plt.subplots()
ax.plot(time, totUtot, marker='x', linestyle='-', color='b')
ax.set_xlabel("Time [MYs]")
ax.set_ylabel("Total Energy [ergs]")
#ax.set_xscale("log")
#ax.set_yscale("log")
ax.grid(True)
fig.savefig(simulation_directory+"/plots/TotalU")



""" normal exit """
sys.exit(0) 
