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
from copy import deepcopy

import json

from scipy.stats import chi2

# %%
SaveFitFigs = False
SaveFitData = False
dpiN = 1000
dark_plots = True
n_sig = 8
n_print_sigfigs = 2
if dark_plots:
    dark='darkbg/'
    q = mpl.rc_params_from_file('matplotlibrc_dark')
else:
    dark = 'whitebg/'
    mpl.rcParams.update(mpl.rcParamsDefault)
SavePlotDir_Exp2  = '../Results/2021-12-21_twosigfigs/Exp2/'+dark+'FittingFigs/'
SaveDataDir_Exp2  = '../Results/2021-12-20/Exp2/'+'Pickles/'
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
    # %matplotlib inline
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
    # %matplotlib inline
    # mpl.rcParams.update(mpl.rcParamsDefault)
    font = {
    #    'weight' : 'normal',
       'size'   : 15,
       'family': 'Times New Roman'}
    plt.rc('font', **font)
#     mpl.rcParams.update({'font.family':'serif'})

# %%
# %load_ext autoreload

# %%
from B_calc_script import TopFunctionOneExpAnyFreq
from B_calc_script import signif

# %%
# %autoreload 2

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
rtr_dims = deepcopy(params_dims_locs['rotor_dims'])

for key in rtr_dims:
    rtr_dims[key] = signif(rtr_dims[key],n_sig)

# %%
Exp2_AW_sensor_loc =params_dims_locs['AW_location']
string_to_parse = params_dims_locs['AW_location']['location']
Exp2_AW_sensor_loc['location'] = eval(string_to_parse.replace('rotor_dims','rtr_dims').replace('D_wheel_sensor','params_dims_locs[\'D_wheel_sensor\']'))

# %%
Exp2_AV_sensor_loc =params_dims_locs['AV_location']
string_to_parse = params_dims_locs['AV_location']['location']
Exp2_AV_sensor_loc['location'] = eval(string_to_parse.replace('rotor_dims','rtr_dims').replace('D_wheel_sensor','params_dims_locs[\'D_wheel_sensor\']'))

# %% [markdown]
# # Calculate fitted field, chi, and plot #

# %%
Exp2_settings = {
    'rotor dimensions':rtr_dims,
    'sensor locations':{
        'AW':Exp2_AW_sensor_loc,
        'AV':Exp2_AV_sensor_loc},
    # 'bar location':0,
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
nowtext = '3sources'
fitplotfilename = SavePlotDir_Exp2+'FittedData_{}Hz_'.format(nu)+nowtext+'.png'
fitdatafilename = SaveDataDir_Exp2+'FittedData_{}Hz_'.format(nu)+nowtext+'.pk'

Exp2_plot_settings = {
    'plot':True,
#     'memo':'{} Hz (AV X&Y inverted)'.format(nu),
    'memo':'{} Hz anticlockwise'.format(nu),
    'saveplot':SaveFitFigs,
    'figname':fitplotfilename,
    'dpi' : dpiN,
    'doubleplot':False,
    'print sigfigs':n_print_sigfigs
}

Exp2_optimization_settings = {
    'print':True,
    'number of sources':3,
    'location dimensions':3,
    'moment dimensions':3,
    'location coordinate system':'polar',
    'moment coordinate system':'polar',
    'chi tolerance':n_sig+1,
    'optimize DC shifts':True,
    'optimize bar location':True,
    'significant figures':n_sig
}
Exp2_save_settings ={
    'save fit data':False,
}
Exp2_all_settings = {
    'experiment settings':Exp2_settings,
    'data':Exp2_data,
    'optimization settings':Exp2_optimization_settings,
    'plot settings':Exp2_plot_settings,
    'save settings':Exp2_save_settings
    
}

Exp2_Opt_Params_3_sources = [ 13.34880956,   8.93204632,   4.54077966,  13.86327989,
       -18.5934633 , -21.08123588,  15.02431798, 274.11433097,
        -4.89142693,  12.91915001, -39.7559461 ,  21.92261841,
        15.14663565,  93.49895702,  -3.9010842 ,  37.87557816,
       -30.16727982, -33.31713104, -34.93498927,  24.45478623,
         7.14475312, -34.82079508,  93.56467221]

chi = TopFunctionOneExpAnyFreq(Exp2_Opt_Params_3_sources,Exp2_all_settings)

# %%
fitplotfilename

# %%
initial = []
initial.append(Exp2_Opt_Params_3_sources)

for ip in range(len(Exp2_Opt_Params_3_sources)):
#     print(ip)
    arr_dum = Exp2_Opt_Params_3_sources.copy()
    arr_dum[ip] = Exp2_Opt_Params_3_sources[ip]*1.1
#     print(arr_dum)
    initial.append(arr_dum)

# %%



# %%
Exp2_all_settings['plot settings']['plot'] = 0
Exp2_all_settings['optimization settings']['print'] = 0
# opt_res_Exp2 = minimize(TopFunctionOneExpAnyFreq,Exp2_Opt_Params,Exp2_all_settings,
#                         method = 'Nelder-Mead',options = {'maxfev': 10000,'maxiter': 10000,'adaptive':True,'fatol':.01})
# opt_res_Exp2 = minimize(TopFunctionOneExpAnyFreq,Exp2_Opt_Params,Exp2_all_settings,
#                         method = 'Powell',options = {'maxfev': 1000,'maxiter': 1000,'ftol':1})
opt_res_Exp2 = minimize(TopFunctionOneExpAnyFreq,Exp2_Opt_Params_3_sources,Exp2_all_settings,
                        method = 'Nelder-Mead',options = {'maxfev': 30000,'maxiter': 20000,'fatol':1e-2,'xatol' : 1e-3
                                                         ,'initial_simplex':initial
                                                         })

# %%
Exp2_all_settings['plot settings']['plot'] = True
Exp2_plot_settings['saveplot']=SaveFitFigs
Exp2_all_settings['optimization settings']['print'] = 1
Exp2_all_settings['save settings']['save fit data'] = SaveFitData
Exp2_all_settings['save settings']['fit data filename'] = fitdatafilename
print(opt_res_Exp2['x'])

# nowtext = now().strftime("%Y%m%d%H%M")
Exp2_plot_settings['figname']=SavePlotDir_Exp2+'FittedData_{}Hz'.format(nu)+nowtext+'optimized'+'.png'
chi_opt = TopFunctionOneExpAnyFreq(opt_res_Exp2['x'],Exp2_all_settings)

# %%
data_to_write = opt_res_Exp2['x']

# %%
data_to_write

# %%
if SaveFitData:
    with open(SaveDataDir_Exp2+'FittedDipoles_{}Hz_'.format(nu)+nowtext+'.pk','wb') as file_obj:
        pickle.dump(data_to_write,file_obj)
    with open(SaveDataDir_Exp2+'FittedDipoles_{}Hz_'.format(nu)+nowtext+'.json','w',encoding = 'utf-8') as file_obj:
        json.dump(data_to_write.tolist(),file_obj)
    

# %%



