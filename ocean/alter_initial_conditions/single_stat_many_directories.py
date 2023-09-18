'''
plot global stats in time
2023, Mark Petersen, LANL
'''

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from datetime import date

wd = '/lustre/scratch5/mpeterse/runs/'

del2or4=2
if del2or4==2:
    casename = 's10'; abc = 'abcd'
    titleText = 'MPAS-Ocean u_t = nu2 del2 u, orig. form.'
    casename = 's09'; abc = 'mnop'
    titleText = 'MPAS-Ocean u_t = nu2 del2 u, stencil form.'
    casename = 's10'; abc = 'efgh'
    titleText = 'MPAS-Ocean u_t = nu2 del2 u, half size domain'
elif del2or4==4:
    casename = 's09'; abc = 'qrst'
    titleText = 'MPAS-Ocean u_t = -nu4 del4 u'
Nyc = '8421'
fileName = '/analysis_members/globalStats.0001-01-01_00.00.00.nc'
var = 'normalVelocityMax'
lt = '-:'

# this plot only:
ny = 16.0
dc = 30.0e3
Ly = ny*dc*np.sqrt(3)/2
nu2 = 1000.0
nu4 = 1.2e11
kyv = [1,2,4,8]
pi2 = 2.0* np.pi
col = 'kbgr'
   
for j in range(len(abc)):
    ncDS = Dataset(wd+casename+abc[j]+fileName,'r') 
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
plt.figtext(.1, .01, casename+abc)
plt.figtext(.9, .01, date.today().strftime("%m/%d/%y") )

plt.savefig('heat_eqn_del{}_'.format(del2or4)+casename+abc+'.png')

