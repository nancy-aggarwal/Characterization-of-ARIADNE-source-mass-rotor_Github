# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
# from IPython import get_ipython

# %%
from scipy.io import loadmat
import numpy as np
from scipy.optimize import minimize
from datetime import datetime
now = datetime.now

import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import os
import pickle

import json

from scipy.stats import chi2


# %%
SaveFitFigs = True
SaveFitData = True
dpiN = 700
dark_plots = False
if dark_plots:
    dark='darkbg/'
    q = mpl.rc_params_from_file('matplotlibrc_dark')
else:
    dark = 'whitebg/'
    mpl.rcParams.update(mpl.rcParamsDefault)
SavePlotDir_Exp2  = '../Results/Test/'+dark+'FittingFigs/'
SaveDataDir_Exp2  = '../Results/Test/'+'Pickles/'
LoadDataDir_Exp2 = SaveDataDir_Exp2 # The other notebook stored the pickle in the same folder
if SaveFitFigs:
    if not os.path.exists(SavePlotDir_Exp2):
        os.makedirs(SavePlotDir_Exp2)
if SaveFitData:
    if not os.path.exists(SaveDataDir_Exp2):
        os.makedirs(SaveDataDir_Exp2)


# %%
if dark_plots:
    mpl.rcParams.update(q)
    # get_ipython().run_line_magic('matplotlib', 'inline')
    mpl.rcParams.update({
                    #'legend.borderpad': 0.3,
                    #'legend.borderaxespad': 0.25,
#                     'legend.columnspacing': 0.6,
#                     'legend.handlelength': 0.7,
                    #'legend.handleheight': 0.4,
                    #'legend.handletextpad': 0.2,
#                     'legend.labelspacing': 0.45,
#                     'text.usetex': True,
                    'font.size':13,
                    })
else:
    # get_ipython().run_line_magic('matplotlib', 'inline')
    # mpl.rcParams.update(mpl.rcParamsDefault)
    font = {
    #    'weight' : 'normal',
       'size'   : 15,
       'family': 'Times New Roman'}
    plt.rc('font', **font)
#     mpl.rcParams.update({'font.family':'serif'})


# %%
# get_ipython().run_line_magic('load_ext', 'autoreload')


# %%
from B_calc_script import TopFunctionOneExpAnyFreq


# %%
# get_ipython().run_line_magic('autoreload', '2')


# %%


# %% [markdown]
# # Load data

# %%
Exp2_data_filename = LoadDataDir_Exp2+'Exp2_cut_averaged_data.pk'


# %%
with open(Exp2_data_filename,'rb') as file_obj:
    Exp2_data_cut = pickle.load(file_obj)

# %% [markdown]
# ## Load parameters ##

# %%
with open('../Params/Exp2_dimensions_and_locations.json', 'r') as fp:
    params_dims_locs = json.load(fp)


# %%
params_dims_locs


# %%
rtr_dims = params_dims_locs['rotor_dims']


# %%
Exp2_AW_sensor_loc =params_dims_locs['AW_location']
string_to_parse = params_dims_locs['AW_location']['location']
Exp2_AW_sensor_loc['location'] = eval(string_to_parse.replace('rotor_dims','rtr_dims').replace('D_wheel_sensor','params_dims_locs[\'D_wheel_sensor\']'))


# %%
Exp2_AV_sensor_loc =params_dims_locs['AV_location']
string_to_parse = params_dims_locs['AV_location']['location']
Exp2_AV_sensor_loc['location'] = eval(string_to_parse.replace('rotor_dims','rtr_dims').replace('D_wheel_sensor','params_dims_locs[\'D_wheel_sensor\']'))


# %%
with open('../Params/Params_3sources.pk','rb') as filehandle:
    Exp2_Opt_Params_3_sources = pickle.load(filehandle)

# %%


# %% [markdown]
# # Calculate fitted field, chi, and plot #

# %%
Exp2_settings = {
    'rotor dimensions':rtr_dims,
    'sensor locations':{
        'AW':Exp2_AW_sensor_loc,
        'AV':Exp2_AV_sensor_loc},
    'bar location':0,
#     'DC shifts':[DC_shift_AVx,DC_shift_AVy,DC_shift_AWy,DC_shift_AWz]
#    'deltaB':1 #picoTesla
}


# %%

nu = 5

Exp2_data = {
    'theta':np.concatenate([Exp2_data_cut['theta avg'][nu]
    # ,360+Exp2_data_cut['theta avg'][nu]
    ]), 
    #theta positive for ac, negative for clockwise
    'B':{
        'AW':{
            'Z':np.concatenate([
                # Exp2_data_cut['AW']['Z avg'][nu],
                                Exp2_data_cut['AW']['Z avg'][nu]['B']]),
            'Y':np.concatenate([
                # Exp2_data_cut['AW']['Y avg'][nu],
                                Exp2_data_cut['AW']['Y avg'][nu]['B']])
        },
        'AV':{
#             'X':np.concatenate([Exp2_data_cut['AV']['Z avg'][nu]+20,Exp2_data_cut['AV']['Z avg'][nu]+20]),
            'X':np.concatenate([
                # Exp2_data_cut['AV']['X avg'][nu],
                                Exp2_data_cut['AV']['X avg'][nu]['B']]),
            'Y':np.concatenate([
                # Exp2_data_cut['AV']['Y avg'][nu],
                                Exp2_data_cut['AV']['Y avg'][nu]['B']])
#             'Y':np.concatenate([-Exp2_data_cut['AV']['Y avg'][nu]-70,-Exp2_data_cut['AV']['Y avg'][nu]-70])
        }        
    },
    'error in B':{
        'AW':{
            'Z':Exp2_data_cut['AW']['Z avg'][nu]['sigma'],
            'Y':Exp2_data_cut['AW']['Y avg'][nu]['sigma']
        },
        'AV':{
            'X':Exp2_data_cut['AV']['X avg'][nu]['sigma'],
            'Y':Exp2_data_cut['AV']['Y avg'][nu]['sigma']
        }
    }
            }


# %%
# nowtext = now().strftime("%Y%m%d%H%M")

nowtext = '_15font'
fitplotfilename = SavePlotDir_Exp2+'FittedData_{}Hz'.format(nu)+nowtext+'.png'
fitdatafilename = SaveDataDir_Exp2+'FittedData_{}Hz'.format(nu)+nowtext+'.pk'


Exp2_optimization_settings = {
    'print':True,
    'number of sources':3,
    'location dimensions':3,
    'moment dimensions':3,
    'location coordinate system':'polar',
    'moment coordinate system':'polar',
    'chi tolerance':10,
    'optimize DC shifts':True,
    'optimize bar location':True
}

Exp2_plot_settings = {
    'plot':True,
#     'memo':'{} Hz (AV X&Y inverted)'.format(nu),
    'memo':'{} Hz'.format(nu),
    'doubleplot':False,
    'saveplot':SaveFitFigs,
    'dpi':dpiN,
    'figname':fitplotfilename
}
Exp2_save_settings ={
    'save fit data':SaveFitData,
    'fit data filename':fitdatafilename
}

Exp2_all_settings = {
    'experiment settings':Exp2_settings,
    'data':Exp2_data,
    'optimization settings':Exp2_optimization_settings,
    'plot settings':Exp2_plot_settings,
    'save settings':Exp2_save_settings
    
}

Exp2_Opt_Params = Exp2_Opt_Params_3_sources
E_opt = TopFunctionOneExpAnyFreq(Exp2_Opt_Params,Exp2_all_settings)


# %%
fitplotfilename

# %% [markdown]
# # Get $\chi^2$ from Error Function ##

# %%
N_points = 4*len(Exp2_data_cut['theta avg'][nu])
N_params = len(Exp2_Opt_Params)


# %%
chi2_opt = E_opt*N_points*N_points
dof = N_points-N_params


# %%
sf_opt = chi2.sf(chi2_opt,dof)


# %%
print('Error function value is {}'.format(E_opt))
print('Number of points is ',N_points)
print("$\chi^2$ is {}".format(chi2_opt))
print('Number of parameters is ',N_params)
print('degrees of freedom is ',dof)
print("$\chi^2$/dof is {}".format(chi2_opt/dof))
print("Survival fraction is {}".format(sf_opt))


# %%



# %%



# %%



# %%



# %%



