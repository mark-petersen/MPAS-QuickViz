'''
plot global stats
August 2016, Mark Petersen, LANL
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

wd = '/lustre/scratch5/turquoise/mpeterse/runs/s08'
del2or4=2
if del2or4==2:
    abc = 'rstu'
    titleText = 'MPAS-Ocean u_t = nu2 del2 u'
elif del2or4==4:
    abc = 'vwxy'
    titleText = 'MPAS-Ocean u_t = -nu4 del4 u'
Nyc = '8421'
fileName = '/analysis_members/globalStats.0001-01-01_00.00.00.nc'
var = 'normalVelocityMax'
lt = '-:'

# this plot only:
Ly = 16.0*30.0e3
nu2 = 1000.0
nu4 = 1.2e11
kyv = [1,2,3,4]
pi2 = 2.0* np.pi
col = 'kbgr'
   
for j in range(len(abc)):
    ncDS = Dataset(wd+abc[j]+fileName,'r') 
    x = ncDS.variables['daysSinceStartOfSim']
    y = ncDS.variables[var]
    if del2or4==2:
        sol = np.exp(- nu2 * (kyv[j]*pi2/Ly)**2 *x[:]*86400 )
    elif del2or4==4:
        sol = np.exp(- nu4 * (kyv[j]*pi2/Ly)**4 *x[:]*86400 )

    plt.plot(x[:], y[:]/y[0], col[j], label='MPAS Kmax/'+Nyc[j])
    plt.plot(x[:], sol[:]/sol[0], '--'+col[j], label='Sol Kmax/'+Nyc[j])
    plt.xlabel('time, days')
    plt.ylabel(var)
    #plt.grid(True,which="both",ls="dotted")
    plt.legend()
    ncDS.close()

plt.title(titleText)
plt.savefig('heat_eqn_del{}.png'.format(del2or4))

