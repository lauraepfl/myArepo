""" @package ./examples/mhd_shocktube_1d/check.py
Code that checks results of 1d mhd shocktube problem

created by Rainer Weinberger, last modified 03.07.2020
"""

""" load libraries """
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False


simulation_directory = str(sys.argv[1])
print("mhd_shocktube_1d: checking simulation output in directory " + simulation_directory) 

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
CellsPerDimension = np.sqrt(NumberOfCells) ## 2d sim


""" loop over all output files """
status = 0
i_file = 1
while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/output/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    
    """ get simulation data """
    
    x = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    vel = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    rho = np.array(data["PartType0"]["Density"], dtype = FloatType)
    u = np.array(data['PartType0']['InternalEnergy'], dtype = FloatType)
    B = np.array(data["PartType0"]["MagneticField"], dtype = FloatType)
    
    absB = np.sqrt(B[:,0]*B[:,0] + B[:,1]*B[:,1] + B[:,2]*B[:,2]) / np.sqrt(4. * np.pi)
    
    
    verificationData = np.array([x[:,0], rho, vel[:,0], u, absB]).T
    
    checkData = np.loadtxt(simulation_directory+'/alpha_3.dat')
    checkData[:,0] += 1.0
    rho_ref = np.interp(x[:,0], checkData[:,0], checkData[:,1])
    vel_ref = np.interp(x[:,0], checkData[:,0], checkData[:,2])
    u_ref = np.interp(x[:,0], checkData[:,0], checkData[:,8]/checkData[:,1]*3./2.)
    absB_ref = np.interp(x[:,0], checkData[:,0], np.sqrt(checkData[:,5]**2+checkData[:,6]**2+checkData[:,7]**2) )
    
    
    delta_rho = rho - rho_ref
    delta_vel = vel[:,0] - vel_ref
    delta_u = u - u_ref
    delta_B = absB - absB_ref
    
    if makeplots:
        fig, ax = plt.subplots(4, sharex=True, figsize=np.array([6.9,6.0]))
        fig.subplots_adjust(left = 0.09, bottom = 0.09,right = 0.98, top = 0.98)
        
        ax[0].plot(x[:,0], rho, 'b', zorder=3, label="Arepo")
        ax[0].plot(checkData[:,0], checkData[:,1], c='k', zorder=1, lw=0.7, label="Exact solution")
        ax[0].set_ylabel(r'density')
        
        ax[1].plot(x[:,0], vel[:,0], 'b', zorder=3)
        ax[1].plot(checkData[:,0], checkData[:,2], c='k', zorder=1, lw=0.7)
        ax[1].set_ylabel(r'vel')
        
        ax[2].plot(x[:,0], u, 'b', zorder=3)
        ax[2].plot(checkData[:,0], checkData[:,8]/checkData[:,1]*3./2., c='k', zorder=1, lw=0.7)
        ax[2].set_ylabel(r'u$_{th}$')
        
        ax[3].plot(x[:,0], absB, 'b', zorder=3)
        ax[3].plot(checkData[:,0], np.sqrt(checkData[:,5]**2+checkData[:,6]**2+checkData[:,7]**2), c='k', zorder=1, lw=0.7)
        ax[3].set_ylabel(r'$\left|\mathbf{B}\right|$')
        ax[3].set_xlabel(r'position')
        ax[3].set_xlim([-0.1,2.6])
        
        
        fig.align_ylabels(ax[:])
        ax[0].legend( loc='upper right', frameon=False, fontsize=8 )
        
        if not os.path.exists( simulation_directory+"/plots" ):
          os.mkdir( simulation_directory+"/plots" )
        
        fig.savefig(simulation_directory+'/plots/figure_%03d.pdf'%i_file)
    

    
    res_scaling = 300. / np.float64(len(rho))
    tolerance_rho = 0.04 * res_scaling
    tolerance_vel = 0.15 * res_scaling
    tolerance_u = 0.2 * res_scaling
    tolerance_B = 0.03 * res_scaling
    
    i_compare = np.where((x<0.6) | (x>0.9))[0]
    
    if np.std(delta_rho[i_compare]) > tolerance_rho:
        status += 1
    if np.std(delta_vel[i_compare]) > tolerance_vel:
        status += 1
    if np.std(delta_u[i_compare]) > tolerance_u:
        status += 1
    if np.std(delta_B[i_compare]) > tolerance_B:
        status +=1
    
    print('standard deviations of absolute error and tolerance (density, velocity, int. energy, magnetic field):')
    print(np.std(delta_rho[i_compare]), tolerance_rho)
    print(np.std(delta_vel[i_compare]), tolerance_vel)
    print(np.std(delta_u[i_compare]), tolerance_u)
    print(np.std(delta_B[i_compare]), tolerance_B)
    
    i_file += 1
    
exit(status)
