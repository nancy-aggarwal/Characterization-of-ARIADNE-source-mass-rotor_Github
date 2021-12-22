"""Code for calculating magnetic field from given dipole moments,
and  then calculating chi^2 from a given data.

Give some description here...
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from copy import deepcopy
import pickle

npfuncs = {}
npfuncs['pi'] = np.pi
npfuncs['abs'] = np.abs
npfuncs['sum'] = np.sum
npfuncs['array'] = np.array
npfuncs['size'] = np.size
npfuncs['sqrt'] = np.sqrt
npfuncs['sindeg'] = lambda thetadeg: np.sin(np.pi*thetadeg/180)
npfuncs['cosdeg'] = lambda thetadeg: np.cos(np.pi*thetadeg/180)
#**************************************************************************
#                                Function List
#**************************************************************************
#                                 signif
#                             FieldAtAnyLocation
#                          TopFunctionOneExpAnyFreq
#                                   Parser
#                                  calc_chi
#                               PlotCalcField
#                               PlotEverything
#                                GetAxisNumber
#                                   CalcB
#                          PlotSetupAndMakeLocation
#                             PlotCalcOnlyFields
#                                 PlotFields
#                                CheckValid
#                                GetCartesian
#                                 GetPolar
#                            ReadAndParseLocation
#**************************************************************************

#``````````````````````````````````````````````````````````````````````````
#``````````````````````````````````````````````````````````````````````````
#                                signif
#``````````````````````````````````````````````````````````````````````````
#``````````````````````````````````````````````````````````````````````````


def signif(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


#``````````````````````````````````````````````````````````````````````````
#``````````````````````````````````````````````````````````````````````````
#                         FieldAtAnyLocation
#``````````````````````````````````````````````````````````````````````````
#``````````````````````````````````````````````````````````````````````````

def FieldAtAnyLocation(calculation_parameters_in,settings):
    """
    This is the top level function

    settings:
        * 'experiment settings':
            - 'rotor dimensions'
            - 'sensor locations'
                + '3He'
                    ~ 'location'
            - 'bar location' (optional)
        * 'data'
            - 'theta'
        * 'optimization settings'
            - 'number of sources'
            - 'optimize bar location'
            - 'location dimensions'
            - 'moment dimensions'
            - 'location coordinate system'
            - 'moment coordinate system'
            - 'significant figures'
        * 'plot settings'
            - 'plot'
            - 'saveplot'
            - 'figname'
            - 'memo'
            - 'doubleplot'
            - 'print sigfigs'
    Enter all lengths in mm, all angles in degrees, all magnetic moments in
     1e-11 Am^2, all fields and DC offsets in pT

    """
    # print('In top level function')
    # calculation_parameters = list(np.around(np.array(calculation_parameters_in),2))
    calculation_parameters = list(signif(calculation_parameters_in,settings['optimization settings']['significant figures']))
    # optimization_parameters = optimization_parameters_in
    calculation_dict = deepcopy(settings['experiment settings'])
    data_dict         = {}
    data_dict['theta'] = deepcopy(settings['data']['theta'])
    data_dict['calculated field'] = {}

    #========================================================================
    # Call parser
    #========================================================================
    ParseSuccess = Parser(calculation_parameters,settings,calculation_dict)
    """ Now calculation_dict should have following keywords:
        * 'number of sources'
        * 'rotor dimensions'
        * 'sensor locations'
            - 'AW'
                + 'location'
        * 'source locations'
        * 'source moments'
        * 'DC shift'
                + 'AW'
                    ~ 'Z'
        * 'bar location'
    """
    if ParseSuccess!=1:
        if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
            print('Parsing was not successful, errorcode {}'.format(ParseSuccess))
        return 1000

    #=========================================================================
    # Calculate field
    #=========================================================================
    for sensorname in calculation_dict['sensor locations']:
        sensor_loc = calculation_dict['sensor locations'][sensorname]
        data_dict['calculated field'][sensorname] = {}
        B_AC_3_axis = 0
        # print(sensor_loc)
        for i_source in range(calculation_dict['number of sources']):
            source_m = calculation_dict['source moments'][i_source]
            source_loc = calculation_dict['source locations'][i_source]
            # print('Calculating field for sensor {}, source {}'\
            #  .format(sensorname,i_source))
            # print(source_m)
            theta_rotation_deg = data_dict['theta']+\
                calculation_dict['bar location']
            B_i_3axis = CalcB(source_loc,source_m,sensor_loc,theta_rotation_deg)
            # CalcB takes locations in millimeters, degrees;
            # moments in 1e-11 Am^2
            # returns fields in picoTesla
            B_AC_3_axis += B_i_3axis
        for axis in ['X','Y','Z']:
            axnum = GetAxisNumber(axis)
            data_dict['calculated field'][sensorname][axis] = B_AC_3_axis[axnum,:]
    
    
    #=========================================================================
    # Call PlotCalcField
    #=========================================================================
    if 'plot settings' in settings:
        if settings['plot settings']['plot']:
            [fig,side_ax,front_ax,data_ax]=\
             PlotCalcField(calculation_dict,data_dict,settings['plot settings'])
            if 'saveplot' in settings['plot settings']:
                    if settings['plot settings']['saveplot']:
                        fig.savefig(settings['plot settings']['figname'],bbox_inches='tight',dpi = settings['plot settings']['dpi'])

    if 'save settings' in settings:
        if 'save fit data' in settings['save settings']:
            if settings['save settings']['save fit data']:
                with open(settings['save settings']['fit data filename'],'wb') as file_obj:
                    pickle.dump(data_dict,file_obj)
    return data_dict['calculated field']


#``````````````````````````````````````````````````````````````````````````
#``````````````````````````````````````````````````````````````````````````
#                         TopFunctionOneExpAnyFreq
#``````````````````````````````````````````````````````````````````````````
#``````````````````````````````````````````````````````````````````````````


def TopFunctionOneExpAnyFreq(optimization_parameters_in,settings):
    """
    This is the top level function

    settings:
        * 'experiment settings':
            - 'rotor dimensions'
            - 'sensor locations'
                + 'AW'
                    ~ 'location'
            - 'bar location' (optional)
        * 'data'
            - 'theta'
            - 'B'
                + 'AW'
                    ~ 'Z'
        * 'optimization settings'
            - 'number of sources'
            - 'optimize bar location'
            - 'location dimensions'
            - 'moment dimensions'
            - 'location coordinate system'
            - 'moment coordinate system'
            - 'significant figures'
        * 'plot settings'
            - 'plot'
            - 'saveplot'
            - 'figname'
            - 'memo'
            - 'doubleplot'
            - 'print sigfigs'
    Enter all lengths in mm, all angles in degrees, all magnetic moments in
     1e-11 Am^2, all fields and DC offsets in pT

    """
    # print('In top level function')
    # optimization_parameters = list(np.around(np.array(optimization_parameters_in),2))

    optimization_parameters = list(signif(optimization_parameters_in,settings['optimization settings']['significant figures']))
    # print(optimization_parameters)

    # optimization_parameters = optimization_parameters_in
    optimization_dict = deepcopy(settings['experiment settings'])
    data_dict         = {}
    data_dict['theta'] = deepcopy(settings['data']['theta'])
    data_dict['measured field']  = deepcopy(settings['data'])
    data_dict['calculated field'] = {}

    #========================================================================
    # Call parser
    #========================================================================
    ParseSuccess = Parser(optimization_parameters,settings,optimization_dict)
    """ Now optimiztion_dict should have following keywords:
        * 'number of sources'
        * 'rotor dimensions'
        * 'sensor locations'
            - 'AW'
                + 'location'
        * 'source locations'
        * 'source moments'
        * 'DC shift'
                + 'AW'
                    ~ 'Z'
        * 'bar location'
    """
    if ParseSuccess!=1:
        if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
            print('Parsing was not successful, errorcode {}'.format(ParseSuccess))
        return 1000

    #=========================================================================
    # Calculate chi
    #=========================================================================
    chi_sq_red,n_points_total = calc_chi(optimization_dict,data_dict)
    if 'print' in settings['optimization settings']:
        if settings['optimization settings']['print']:
            print("Total number of points are {}".format(n_points_total))
    # chi_sq_round = round(chi_sq_red,6)
    chi_sq_round = chi_sq_red
    """
    After this function, data_dict['calculated field']
    and data_dict['reduced chi squared']
    should be populated and can be used by the plotter
    """
    if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
        print('Reduced chi sq = {}'.format(chi_sq_round))
    #=========================================================================
    # Call PlotEverything
    #=========================================================================
    if 'plot settings' in settings:
        if settings['plot settings']['plot']:
            [fig,side_ax,front_ax,data_ax]=\
             PlotEverything(optimization_dict,data_dict,settings['plot settings'])
            if 'saveplot' in settings['plot settings']:
                    if settings['plot settings']['saveplot']:
                        fig.savefig(settings['plot settings']['figname'],bbox_inches='tight',dpi = settings['plot settings']['dpi'])

    if 'save settings' in settings:
        if 'save fit data' in settings['save settings']:
            if settings['save settings']['save fit data']:
                with open(settings['save settings']['fit data filename'],'wb') as file_obj:
                    pickle.dump(data_dict,file_obj)
    return chi_sq_round

#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
#                                   Parser
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
def Parser(array_to_parse,settings,optimization_dict):
    """
    Parser takes the same format dictionary settings as the Toplevel function.
    It then populates optimization_dict with values from array_to_parse
    It checks whether the user wants to optimize the location of the bar at
    the starting, then verifies that the number of parameters are exactly
    right.
    It then populates the dictionary that contain source locations, moments,
    DC shifts, and all other settings for the calculations. This dictionary is
    used in the rest of the code.
    Finally, it checks whether the proposed sources lie inside the rotor.
    """

    #========================================================================
    # Decide whether to optimize for bar location or extract it
    #========================================================================
    if ('optimize bar location' in settings['optimization settings']):
        if settings['optimization settings']['optimize bar location']:
            # print("Extracting bar location from parameter list.")
            n_bar = 1
            bar_loc = array_to_parse[-1]
            if 'print' in settings['optimization settings']:
                if settings['optimization settings']['print']:
                    print('Bar location: {}°'.format(bar_loc))
        else:
            if 'print' in settings['optimization settings']:
                if settings['optimization settings']['print']:
                    print("Expecting bar location in 'experiment setting'.")
            n_bar = 0
            bar_loc = settings['experiment settings']['bar location']
    else:
        if 'print' in settings['optimization settings']:
            if settings['optimization settings']['print']:
                print("Expecting bar location in 'experiment setting'.")
        n_bar = 0
        bar_loc = settings['experiment settings']['bar location']

    if bar_loc<0:
        if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
            print('Bar must be greater than 0 degrees')
        return -1
    if bar_loc>180:
        if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
            print('Bar must be less than 180 degrees')
        return -2
    optimization_dict['bar location'] = bar_loc

    #=========================================================================
    # DC shifts
    #=========================================================================
    
    if ('optimize DC shifts' in settings['optimization settings']):
        if settings['optimization settings']['optimize DC shifts']:
            n_DC = 0
            optimization_dict['DC shift']={}
            sensornames = list(settings['data']['B'].keys())
            # sensornames = list(optimization_dict['sensor locations'].keys())
            sensornames.sort(reverse=True)
            n_sensor = 0
            DC_arr = []
            for sensor in sensornames:
                optimization_dict['DC shift'][sensor]={}
                axes = list(settings['data']['B'][sensor].keys())
                axes.sort(reverse=True)
                # print(axes)
                for axis in axes:
                    DC_shift = array_to_parse[-1-n_bar-n_DC]
                    n_DC+=1
                    optimization_dict['DC shift'][sensor][axis] = DC_shift
                    DC_arr.append(DC_shift)
                n_sensor+=1
            optimization_dict['number of sensors'] = n_sensor
        else:
            if 'print' in settings['optimization settings']:
                if settings['optimization settings']['print']:
                    print("Expecting DC shifts in 'experiment setting'.")
            n_DC = 0
            i_DC = 0
            DC_arr_read = settings['experiment settings']['DC shifts']
            n_DC_arr_read = len(DC_arr_read)
            optimization_dict['DC shift']={}
            sensornames = list(settings['data']['B'].keys())
            sensornames.sort(reverse=True)
            n_sensor = 0
            DC_arr = []
            for sensor in sensornames:
                optimization_dict['DC shift'][sensor]={}
                axes = list(settings['data']['B'][sensor].keys())
                axes.sort(reverse=True)
                # print(axes)
                for axis in axes:
                    DC_shift = DC_arr_read[n_DC_arr_read-1-i_DC]
                    i_DC+=1
                    optimization_dict['DC shift'][sensor][axis] = DC_shift
                    DC_arr.append(DC_shift)
                n_sensor+=1
            optimization_dict['number of sensors'] = n_sensor
    else:
        if 'print' in settings['optimization settings']:
            if settings['optimization settings']['print']:
                print("Expecting DC shifts in 'experiment setting'.")
        n_DC = 0
        i_DC = 0
        DC_arr_read = settings['experiment settings']['DC shifts']
        n_DC_arr_read = len(DC_arr_read)
        # print(DC_arr_read)
        # print(n_DC_arr_read)
        optimization_dict['DC shift']={}
        sensornames = list(settings['data']['B'].keys())
        sensornames.sort(reverse=True)
        n_sensor = 0
        DC_arr = []
        for sensor in sensornames:
            optimization_dict['DC shift'][sensor]={}
            axes = list(settings['data']['B'][sensor].keys())
            axes.sort(reverse=True)
            # print(axes)
            for axis in axes:
                DC_shift = DC_arr_read[n_DC_arr_read-1-i_DC]
                i_DC+=1
                optimization_dict['DC shift'][sensor][axis] = DC_shift
                DC_arr.append(DC_shift)
            n_sensor+=1
        optimization_dict['number of sensors'] = n_sensor
    if 'print' in settings['optimization settings']:
        if settings['optimization settings']['print']:
            print('DC shifts: {} pT'.format(DC_arr))

    #=========================================================================
    # Check number of parameters
    #=========================================================================
    n_sources = settings['optimization settings']['number of sources']
    loc_dim   = settings['optimization settings']['location dimensions']
    m_dim     = settings['optimization settings']['moment dimensions']

    anticipated_size = n_sources*(loc_dim+m_dim) + n_DC + n_bar

    if len(array_to_parse)<anticipated_size:
        if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
            print('insufficient parameters. Expected number for {} sensor axes,\
             {} bar location, {}  sources, {} location dimensions,\
             {} moment dimensions is {}'.format(n_DC,n_bar,n_sources,loc_dim,\
             m_dim,anticipated_size))
        return -3
    elif len(array_to_parse)>anticipated_size:
        if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
            print('too many parameters. Expected number for {} sensor axes,\
             {} bar location, {}  sources, {} location dimensions,\
             {} moment dimensions is {}'.format(n_DC,n_bar,n_sources,loc_dim,\
             m_dim,anticipated_size))
        return -4
    # else:
        # print('Correct number of parameters, proceeding to extract rest')

    optimization_dict['number of sources'] = n_sources

    #========================================================================
    # Parse source stuff populate dict if the sources lie inside rotor
    #========================================================================
    source_locations = []
    source_moments = []
    i_item = 0
    i_source = 0
    # now first loc_dim + m_dim elements are for first source, and so on
    while i_item+1 < n_sources*(loc_dim+m_dim):
        if loc_dim>2:
            src_loc = {
             'coordinate':
             settings['optimization settings']['location coordinate system'],
             'location':array_to_parse[i_item:i_item+loc_dim]}
        else:
            src_loc = {
             'coordinate':
             settings['optimization settings']['location coordinate system'],
             'location':list(array_to_parse[i_item:i_item+loc_dim])+[0]}
#             print(src_loc)
        val = CheckValid(src_loc,settings['experiment settings']\
         ['rotor dimensions'])
        if not val:
            if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
                print('Source %0d not in rotor'%(i_source+1))
            return -5
        if src_loc['coordinate']!='cartesian':
            if src_loc['location'][1]>360:
                if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
                    print('Source %0.0d has theta bigger than 2pi'%(i_source+1))
                return -6
            if src_loc['location'][1]<0:
                if 'print chi' not in settings['optimization settings'] or settings['optimization settings']['print chi']:
                    print('Source %0.0d has theta negative'%(i_source+1))
                return -7

        source_locations.append(src_loc)
        if m_dim>2:
            src_m = {
             'coordinate':settings['optimization settings']\
             ['moment coordinate system'],
             'moment':array_to_parse[i_item+loc_dim:i_item+loc_dim+m_dim]}
        else:
            src_m = {
             'coordinate':settings['optimization settings']\
             ['moment coordinate system'],
             'moment':list(array_to_parse\
             [i_item+loc_dim:i_item+loc_dim+m_dim])+[0]}
        source_moments.append(src_m)
        i_item += loc_dim+m_dim
        i_source +=1
        if 'print' in settings['optimization settings']:
            if settings['optimization settings']['print']:
                print('Source {} is located at {}'.format(i_source,src_loc))
                print('Source {} has moment {}'.format(i_source,src_m))
    optimization_dict['source locations']=source_locations
    optimization_dict['source moments']=source_moments
    # print(optimization_dict['source locations'])

    if 'chi tolerance' in settings['optimization settings']:
        optimization_dict['chi tolerance']=settings['optimization settings']\
         ['chi tolerance']
    # print('Parsing done')
    return 1

#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
#                                 calc_chi
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
def calc_chi(optimization_dict,data_dict):
    """
    This function takes the optimization_dict that has information about the
    experiment, and calculates the magnetic field.
    Optimization_dict should have following:
        * 'number of sources'
        * 'rotor dimensions'
        * 'sensor locations'
            - 'AW'
                + 'location'
        * 'source locations'
        * 'source moments'
        * 'DC shift'
                + 'AW'
                    ~ 'Z'
        * 'bar location'

    data_dict has:
        * 'theta'
        * 'measured field'
            - 'AW'
                + 'Y'
        * 'calculated field' - empty, will be filled in this function
    """
    chi_sq_total = 0
    n_points_total = 0
    for sensorname in optimization_dict['sensor locations']:
        sensor_loc = optimization_dict['sensor locations'][sensorname]
        data_dict['calculated field'][sensorname] = {}
        B_AC_3_axis = 0
        # print(sensor_loc)
        for i_source in range(optimization_dict['number of sources']):
            source_m = optimization_dict['source moments'][i_source]
            source_loc = optimization_dict['source locations'][i_source]
            # print('Calculating field for sensor {}, source {}'\
            #  .format(sensorname,i_source))
            # print(source_m)
            theta_rotation_deg = data_dict['theta']+\
             optimization_dict['bar location']
            B_i_3axis = CalcB(source_loc,source_m,sensor_loc,theta_rotation_deg)
            # CalcB takes locations in millimeters, degrees;
            # moments in 1e-11 Am^2
            # returns fields in picoTesla
            B_AC_3_axis += B_i_3axis

        for sensoraxis in data_dict['measured field']['B'][sensorname]:
            axnum = GetAxisNumber(sensoraxis)
            B_AC_this_axis = B_AC_3_axis[axnum,:]
            B_DC_this_axis = optimization_dict['DC shift'][sensorname]\
             [sensoraxis]
            B_total_this_axis = B_AC_this_axis+B_DC_this_axis
            data_dict['calculated field'][sensorname][sensoraxis] = \
             B_total_this_axis
            B_measured_this_axis = data_dict['measured field']['B'][sensorname]\
             [sensoraxis]
            chi_sq_array = np.power((B_total_this_axis  -\
             B_measured_this_axis)/data_dict['measured field']['error in B'][sensorname][sensoraxis],2)
            chi_sq_this_axis = np.sum(chi_sq_array)
            chi_sq_total +=chi_sq_this_axis
            n_points_total += np.size(chi_sq_array)
    chi_sq_reduced = chi_sq_total/(n_points_total**2)
    if 'chi tolerance' in optimization_dict:
        chi_sq_reduced = np.around(chi_sq_reduced,decimals = optimization_dict['chi tolerance'])
    data_dict['reduced chi squared'] = chi_sq_reduced
    return chi_sq_reduced,n_points_total

#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
#                                PlotCalcField
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
def PlotCalcField(calculation_dict,data_dict,plot_settings):
    # Make textstring to print about dipole locations

    n_sources = calculation_dict['number of sources']

    
    # textstring = """                                    

    # Dipole #                         Position               \
    #                      Moment"""
    # for i_source in range(n_sources):
    #     source_loc = deepcopy(
    #      calculation_dict['source locations'][i_source])
    #     source_loc_polar = GetPolar(source_loc)
    #     source_m = calculation_dict['source moments'][i_source]
    # #   source_m_polar = GetPolar(source_m)
    #     strloop = """
    #     %0.0f             (%0.2f mm, %0.1f°, %0.2f mm)\
    #        (%0.1f,%0.1f,%0.1f) $\\times$ 1e-11 A$\mathrm{m}^2$  """\
    #     %(i_source+1,*source_loc_polar['location'],*source_m['moment'])
    #     textstring = textstring+strloop
    #     textstringE = """ $E = $ %0.3e """\
    # %(data_dict['reduced chi squared'])

    textstringD = """ Dipole #                         Position               \
                         Moment"""
    for i_source in range(n_sources):
        source_loc = deepcopy(
         calculation_dict['source locations'][i_source])
        source_loc_polar = GetPolar(source_loc)
        source_m = calculation_dict['source moments'][i_source]
    #   source_m_polar = GetPolar(source_m)
        fitparamslistsigfigs = signif([i_source+1]+source_loc_polar['location']+source_m['moment'],plot_settings['print sigfigs'])
        fitparamslist = []
        for number in fitparamslistsigfigs:
            if int(number) == float(number):
               number = int(number)
            else:
               number = float(number)
            fitparamslist.append(number)
        # source_m_polar = GetPolar(source_m)
        strloop = """
        {}             ({} mm, {}°, {} mm)\
           ({},{},{}) $\\times$ $10^{{-11}}$ Am$^2$  """\
        .format(*fitparamslist)
        textstringD = textstringD+strloop
    
    [fig,side_ax,front_ax,data_ax_spec] = PlotSetupAndMakeLocation(calculation_dict)
    # fig.text(0.1,-0.32,textstring)
    # fig.text(0.4,-0.3,textstring, 
    #     bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    if plot_settings['print sigfigs']<3:
        fig.text(0.4,-0.3,textstringD, 
            bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    else:
            fig.text(0.4,-0.3,textstringD, 
                bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    # fig.text(0.4,-0.32,textstringD, 
    #     bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    if 'memo' in plot_settings:
            titlestring = plot_settings['memo']    
            fig.suptitle(titlestring,fontsize=16,weight='bold')
    PlotCalcOnlyFields(data_dict,fig,data_ax_spec,plot_settings)
    # fig.tight_layout()


    return fig,side_ax,front_ax,data_ax_spec

#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
#                                PlotEverything
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
def PlotEverything(optimization_dict,data_dict,plot_settings):
    # Make textstring to print about dipole locations

    n_sources = optimization_dict['number of sources']

    # textstring = """                                      $E = $ %0.3e

    # Dipole #                         Position               \
    #                      Moment"""\
    # %(data_dict['reduced chi squared'])
    # for i_source in range(n_sources):
    #     source_loc = deepcopy(
    #      optimization_dict['source locations'][i_source])
    #     source_loc_polar = GetPolar(source_loc)
    #     source_m = optimization_dict['source moments'][i_source]
    # #   source_m_polar = GetPolar(source_m)
    #     strloop = """
    #     %0.0f             (%0.2f mm, %0.1f°, %0.2f mm)\
    #        (%0.1f,%0.1f,%0.1f) $\\times$ 1e-11 A$\mathrm{m}^2$  """\
    #     %(i_source+1,*source_loc_polar['location'],*source_m['moment'])
    #     textstring = textstring+strloop


    textstringE = """E$^\dagger$ = {} """.format(signif(data_dict['reduced chi squared'],plot_settings['print sigfigs']))
    textstringD = """ Dipole #                        Position$^\dagger$               \
                  Moment$^\dagger$ """
    for i_source in range(n_sources):
        source_loc = deepcopy(
         optimization_dict['source locations'][i_source])
        source_loc_polar = GetPolar(source_loc)
        source_m = optimization_dict['source moments'][i_source]
        fitparamslistsigfigs = signif([i_source+1]+source_loc_polar['location']+source_m['moment'],plot_settings['print sigfigs'])
        fitparamslist = []
        for number in fitparamslistsigfigs:
            if int(number) == float(number):
               number = int(number)
            else:
               number = float(number)
            fitparamslist.append(number)
        # source_m_polar = GetPolar(source_m)
        strloop = """
        {}             ({} mm, {}°, {} mm)\
           ({},{},{}) $\\times$ $10^{{-11}}$ Am$^2$  """\
        .format(*fitparamslist)
        textstringD = textstringD+strloop


    [fig,side_ax,front_ax,data_ax_spec] = PlotSetupAndMakeLocation(optimization_dict)
    # fig.text(0.1,-0.32,textstring)
    # fig.text(0.4,-0.3,textstring, 
    #     bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    if plot_settings['print sigfigs']<3:
        fig.text(0.37,-0.17,textstringE, 
            bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
        fig.text(0.5,-0.3,textstringD, 
            bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    else:
            fig.text(0.35,-0.17,textstringE, 
                bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
            fig.text(0.46,-0.3,textstringD, 
                bbox=dict(facecolor='none', edgecolor='gray', boxstyle='round'))
    if 'memo' in plot_settings:
        titlestring = plot_settings['memo']
        fig.suptitle(titlestring,fontsize=16,weight='bold')
    PlotFields(data_dict,fig,data_ax_spec,plot_settings)
    # fig.tight_layout()


    return fig,side_ax,front_ax,data_ax_spec

#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
#                                 GetAxisNumber
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
def GetAxisNumber(axisname):
    if axisname=='X' or axisname=='x':
        axnumber = 0
    elif axisname=='Y' or axisname=='y':
        axnumber = 1
    elif axisname=='Z' or axisname=='z':
        axnumber = 2
    else:
        print('Axis can be x, y, or z. Using z as default')
        return 2
    return axnumber
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
#                                    CalcB
#`````````````````````````````````````````````````````````````````````````````
#`````````````````````````````````````````````````````````````````````````````
def CalcB(source_loc,source_m,sensor_loc,theta_rotation_deg):
    # This function will take a single source, a single sensor (everything 3D),
    # and give a 3 dimensional magnetic field at the sensor due to the source
    # as a function of rotation
    # Enter all lengths in mm, all angles in degrees, all magnetic moments in
    # 1e-11 SI, all DC offsets in pT
    mu0_by_4pi = 1e-7
    # First get source parameters
    sensor_cartesian = GetCartesian(sensor_loc)['location']
    source_polar = GetPolar(source_loc)
    source_r = source_polar['location'][0]
    source_theta = source_polar['location'][1]
    source_z = source_polar['location'][2]
    # Now we make a rotation operator and multiply it to m and source location
    theta_rotation_plus_source = theta_rotation_deg + source_theta
    R_theta = np.array([[[npfuncs['cosdeg'](th),-npfuncs['sindeg'](th),0],
                         [npfuncs['sindeg'](th),npfuncs['cosdeg'](th),0],
                         [0,0,1]]
                                  for th in theta_rotation_plus_source])
    source_m_cartesian_theta = \
     np.matmul(R_theta,source_m['moment']).transpose()
    source_r_cartesian_theta = \
     np.matmul(R_theta,[source_r,0,source_z]).transpose()
#     print(sensor_cartesian)
#     print(source_m_cartesian_theta)
#     print(source_r_cartesian_theta)
    d_vec = (source_r_cartesian_theta.transpose()-sensor_cartesian).transpose()

    # Convert lengths to meters, and moments to SI
    d_vec = 1e-3*d_vec
    source_m_cartesian_theta = 1e-11*source_m_cartesian_theta

    #print(d_vec.shape)
    #plt.plot(theta,m_theta[0,:])
    #plt.plot(theta,m_theta[1,:])
    #plt.plot(theta,m_theta[2,:])
    mdotd = np.diagonal(np.matmul(d_vec.transpose(),source_m_cartesian_theta))
   # print(mdotd.shape)
    comp1 = 3*np.divide(np.multiply(d_vec,mdotd),
     np.linalg.norm(d_vec,axis=0)**2)
#    print(comp1.shape)
    B_vec = mu0_by_4pi/np.linalg.norm(d_vec,axis=0)**3 *\
     (comp1-source_m_cartesian_theta)
#    print(B_vec.shape)

    return B_vec*1e12

#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                          PlotSetupAndMakeLocation
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def PlotSetupAndMakeLocation(input_dict):

    ## What all is needed in input_dict:
    # Source locations
    # Source moments
    # DC offsets
    # Sensor locations
    # Sensor names
    # Sensor axis
    # Bar location for single freq

    #======================================================================
    # Make arrays that will be used to draw the rotor
    #======================================================================
    th = np.arange(0,360,1)
    rotor_out_x = input_dict['rotor dimensions']['outer_radius']*\
     npfuncs['cosdeg'](th)
    rotor_out_y = input_dict['rotor dimensions']['outer_radius']*\
     npfuncs['sindeg'](th)
    rotor_in_x = input_dict['rotor dimensions']['inner_radius']*\
     npfuncs['cosdeg'](th)
    rotor_in_y = input_dict['rotor dimensions']['inner_radius']*\
     npfuncs['sindeg'](th)
    rotor_hole_x = input_dict['rotor dimensions']['hole_radius']*\
     npfuncs['cosdeg'](th)
    rotor_hole_y = input_dict['rotor dimensions']['hole_radius']*\
     npfuncs['sindeg'](th)

    r = input_dict['rotor dimensions']['inner_radius']
    bar_theta = np.deg2rad(input_dict['bar location'])
    a = input_dict['rotor dimensions']['bar_width']/2
    phi_pos = bar_theta + np.arcsin(a/r)
    phi_neg = bar_theta - np.arcsin(a/r)

    # Make bar at the angle bar location
    bar_theta_x_up = [r*np.cos(phi_pos),-r*np.cos(phi_neg)]
    bar_theta_y_up = [r*np.sin(phi_pos),-r*np.sin(phi_neg)]
    bar_theta_x_down = [r*np.cos(phi_neg),-r*np.cos(phi_pos)]
    bar_theta_y_down = [r*np.sin(phi_neg),-r*np.sin(phi_pos)]

    rotor_out_z = \
     input_dict['rotor dimensions']['height']/2*np.ones(len(rotor_out_x))
    vert_z = np.arange(-input_dict['rotor dimensions']['height']/2,
     input_dict['rotor dimensions']['height']/2,
     input_dict['rotor dimensions']['height']/20)
    vert_x =\
     input_dict['rotor dimensions']['outer_radius']*np.ones(len(vert_z))
    vert_y =\
     input_dict['rotor dimensions']['outer_radius']*np.ones(len(vert_z))

    # rtr_color = 'deepskyblue'
    rtr_color = 'slategray'

    #====================================================================
    # Decide figure size and subplot sizes
    #====================================================================
    max_x = max(rotor_out_x)
    min_x = min(rotor_out_x)
    max_y = max(rotor_out_y)
    min_y = min(rotor_out_y)
    max_z = max(vert_z)+0.5
    min_z = min(vert_z)-0.5
    # Sources are inside the rotor by default, so only need to check the
    # sensors to decide plot boundaries
    for sensor in input_dict['sensor locations']:
        sensor_cartesian = GetCartesian(input_dict['sensor locations'][sensor])
        max_x = max(max_x,sensor_cartesian['location'][0])
        min_x = min(min_x,sensor_cartesian['location'][0])
        max_y = max(max_y,sensor_cartesian['location'][1])
        min_y = min(min_y,sensor_cartesian['location'][1])
        max_z = max(max_z,sensor_cartesian['location'][2])
        min_z = min(min_z,sensor_cartesian['location'][2])

    x_size = max_x - min_x
    y_size = max_y - min_y
    z_size = max_z - min_z

    total_height = 0.1*(y_size)+0.2 # space for annotation on top
    total_width = 0.2*(x_size + z_size )+2.8 # space for axes ticks

    #===================================================================
    # Initialize figure
    #===================================================================
    # fig, [side_view,front_view,fields] = plt.subplots(nrows=1, ncols=3
    #  ,gridspec_kw={'width_ratios': [z_size/x_size,1,1.3]}
    #  ,figsize = (total_width,total_height)
    #  )
    fig = plt.figure(figsize = (total_width,total_height))
    [rotor_spec,fields_spec] = gridspec.GridSpec(nrows=1,ncols=2,wspace=0.2)
    # fields = plt.Subplot(fig,fields_spec)
    [side_view_spec,front_view_spec] = gridspec.GridSpecFromSubplotSpec(nrows=1,ncols=2
     ,width_ratios= [z_size/x_size,1]
     ,subplot_spec=rotor_spec
     ,wspace=0.1
     )
    side_view = plt.Subplot(fig,side_view_spec)
    front_view = plt.Subplot(fig,front_view_spec
    ,sharey = side_view
    )
    plt.setp(front_view.get_yticklabels(), visible=False)
    # fig.add_subplot(fields)
    fig.add_subplot(side_view)
    fig.add_subplot(front_view)

    # side_view.get_shared_y_axes().join(side_view,front_view)
    # front_view.set_yticklabels([])

    #==================================================================
    # Plot rotor and sensors and Sources
    #==================================================================

    side_view.set_xlabel('z (mm)')
    front_view.set_xlabel('x (mm)')
    side_view.set_ylabel('y (mm)')

    side_view.plot(-rotor_out_z,rotor_out_y,color = rtr_color)
    side_view.plot(rotor_out_z,rotor_out_y,color = rtr_color)
    side_view.plot(vert_z,vert_y,vert_z,-vert_y,color = rtr_color)

    front_view.plot(rotor_out_x,rotor_out_y,
     label = 'rotor boundary',color = rtr_color)
    front_view.plot(rotor_in_x,rotor_in_y,color = rtr_color)
    front_view.plot(rotor_hole_x,rotor_hole_y,color = rtr_color)

    front_view.plot(bar_theta_x_up,bar_theta_y_up,color = rtr_color)
    front_view.plot(bar_theta_x_down,bar_theta_y_down,color = rtr_color)

    # numberofdots = input_dict['number of sources'] + input_dict['number of sensors']

    # cmap = plt.get_cmap('Accent')
    # colors = [cmap(i_col) for i_col in np.linspace(0, 1, numberofdots)]

    # colors = ['magenta','purple','red','violet','tomato', 'darkviolet']
    # n_sensors = len(input_dict['sensor locations'].keys())
    # colors = plt.get_cmap('tab10',input_dict['number of sensors'])
    colors = plt.get_cmap('Set1')
    i_col = 0
    for sensor in input_dict['sensor locations']:
        sensor_cartesian = GetCartesian(input_dict['sensor locations'][sensor])
        front_view.scatter(sensor_cartesian['location'][0],
         sensor_cartesian['location'][1],label = 'sensor {}'.format(sensor),
         # color = colors[i_col]
         color = colors(i_col)
         ,alpha = 0.6
         ,marker = 'o',s = 60)
        side_view.scatter(sensor_cartesian['location'][2],
         sensor_cartesian['location'][1]
         # ,color = colors[i_col]
         ,color = colors(i_col)
         ,alpha = 0.6
         ,marker='o',
         s=60)
        i_col +=1

    # colors = ['yellowgreen','olivedrab','darkgoldenrod','lime','olive']
    i_col = 0
    # cmap = plt.get_cmap('jet')
    # colors = [cmap(i_col) for i_col in np.linspace(0, 1, input_dict['number of sources'])]
    # colors = plt.get_cmap('Dark2',input_dict['number of sources'])
    colors = plt.get_cmap('Accent')
    for i_source in range(input_dict['number of sources']):
        source_polar     = GetPolar(
         input_dict['source locations'][i_source])
        source_polar_shifted = deepcopy(source_polar)
        source_polar_shifted['location'][1] = source_polar_shifted['location'][1]\
         + input_dict['bar location']
        source_cartesian = GetCartesian(source_polar_shifted)
        front_view.scatter(source_cartesian['location'][0],
         source_cartesian['location'][1],label = 'source ' + str(i_source+1),
         # color = colors[i_col]
         color = colors(i_col)
         ,marker='*',s=140)
        side_view.scatter(source_cartesian['location'][2],
         source_cartesian['location'][1]
         # ,color = colors[i_col]
         ,color = colors(i_col)
         ,marker='*',s=140)
        i_col +=1


    front_view.legend(loc='upper left',bbox_to_anchor=(-0.65,-0.2),ncol=2)
    front_view.axis('equal')
    side_view.axis('equal')

    return fig,side_view,front_view,fields_spec

#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                                PlotCalcOnlyFields
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def PlotCalcOnlyFields(input_dict,fig,data_spec,plot_settings):

    # colors = ['red','green','blue','magenta','orange','brown','pink','olive']


    [blah,data_plots_space_spec] = gridspec.GridSpecFromSubplotSpec(nrows=1,ncols=2
     ,subplot_spec=data_spec
     ,width_ratios = [0.001,0.999]
     ,wspace = 0.15
     )
    plt.Subplot(fig,blah)
    fig.add_subplot(blah,frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.ylabel('Magnetic Field (pT)')
    n_plot = 3
    # for sensor in input_dict['calculated field']:
    # for axis in input_dict['calculated field'][sensor]:
    #         n_plot +=1

    # cmap = plt.get_cmap('tab10')
    # colors = [cmap(i_col) for i_col in np.linspace(0, 1, n_plot)]

    data_plots_specs = gridspec.GridSpecFromSubplotSpec(nrows=n_plot,ncols=1
     ,subplot_spec=data_plots_space_spec
     ,hspace = 0
     )
    i_plot = n_plot-1
    theta_plot = input_dict['theta']/180

    # sensor_colors = plt.get_cmap('tab10',len(input_dict['calculated field'].keys()))
    sensor_colors = plt.get_cmap('Set1')
    
    for axis in ['X','Y','Z']:
        if i_plot == n_plot-1:
            ax_bottom = plt.Subplot(fig,data_plots_specs[i_plot])
            ax_i = ax_bottom
        else:
            ax_i = plt.Subplot(fig,data_plots_specs[i_plot]
            ,sharex = ax_bottom
            )
            plt.setp(ax_i.get_xticklabels(), visible=False)
        fig.add_subplot(ax_i)
        # h = ax_i.plot(theta_plot,input_dict['measured field']['B'][sensor][axis]
        # # ,color = colors[i_plot]
        # ,color = "C{}".format(i_plot)
        # # ,alpha = 0.8
        # ,linewidth = 2
        # )
        icol = 0
        for sensor in input_dict['calculated field']:
            ax_i.plot(theta_plot,input_dict['calculated field'][sensor][axis]
                # ,color = h[-1].get_color()
                ,color = sensor_colors(icol)
                ,linewidth = 3
                ,alpha = 0.6
                # ,linestyle = '--'
                )
            if 'doubleplot' in plot_settings:
                if plot_settings['doubleplot']:
                    # ax_i.plot(theta_plot+2,input_dict['measured field']['B'][sensor][axis]
                    # # ,color = colors[i_plot]
                    # ,color = "C{}".format(i_plot)
                    # # ,alpha = 0.8
                    # ,linewidth = 2
                    # )

                    ax_i.plot(theta_plot+2,input_dict['calculated field'][sensor][axis]
                    # ,color = h[-1].get_color()
                    ,color = sensor_colors(icol)
                    ,linewidth = 3
                    ,alpha = 0.6
                    # ,linestyle = '--'
                    )
            icol +=1

            

        ax_i.set_ylabel('{}'.format(axis))
        plt.grid()
        # if i_plot!=3:
        #      plt.setp(ax_data_plots[i_plot].get_xticklabels(), visible=False)
        i_plot-=1

    # lgd = pltax.legend(bbox_to_anchor=(1.35, 1.3))
    # pltax.legend(loc='best')
    # pltax.grid()
    # xmin,xmax = pltax.get_xlim()
    # ymin,ymax = pltax.get_ylim()
    # xloc = xmin - 1*(xmax-xmin)
    # yloc = ymax + 0.05*(ymax-ymin)
    # pltax.text(xloc,yloc,textstring)
    ax_bottom.set_xlabel('$\\theta/\pi$')
    # pltax.set_ylabel('Magnetic field (pT)')


#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                                PlotFields
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def PlotFields(input_dict,fig,data_spec,plot_settings):

    # colors = ['red','green','blue','magenta','orange','brown','pink','olive']


    [blah,data_plots_space_spec] = gridspec.GridSpecFromSubplotSpec(nrows=1,ncols=2
     ,subplot_spec=data_spec
     ,width_ratios = [0.001,0.999]
     ,wspace = 0.15
     )
    plt.Subplot(fig,blah)
    fig.add_subplot(blah,frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.ylabel('Magnetic Field (pT)')
    n_plot = 0
    for sensor in input_dict['calculated field']:
        for axis in input_dict['calculated field'][sensor]:
            n_plot +=1

    # cmap = plt.get_cmap('tab10')
    # colors = [cmap(i_col) for i_col in np.linspace(0, 1, n_plot)]

    data_plots_specs = gridspec.GridSpecFromSubplotSpec(nrows=n_plot,ncols=1
     ,subplot_spec=data_plots_space_spec
     ,hspace = 0
     )
    i_plot = n_plot-1
    theta_plot = input_dict['theta']/180
    for sensor in input_dict['calculated field']:
        for axis in input_dict['calculated field'][sensor]:
            if i_plot == n_plot-1:
                ax_bottom = plt.Subplot(fig,data_plots_specs[i_plot])
                ax_i = ax_bottom
            else:
                ax_i = plt.Subplot(fig,data_plots_specs[i_plot]
                ,sharex = ax_bottom
                )
                plt.setp(ax_i.get_xticklabels(), visible=False)
            fig.add_subplot(ax_i)
            h = ax_i.plot(theta_plot,input_dict['measured field']['B'][sensor][axis]
            # ,color = colors[i_plot]
            ,color = "C{}".format(i_plot)
            # ,alpha = 0.8
            ,linewidth = 2
            )
            ax_i.plot(theta_plot,input_dict['calculated field'][sensor][axis]
             # ,color = h[-1].get_color()
             ,color = [0.4,0.4,0.4]
             ,linewidth = 3
             ,alpha = 0.6
             # ,linestyle = '--'
             )
            if 'doubleplot' in plot_settings:
                if plot_settings['doubleplot']:
                    ax_i.plot(theta_plot+2,input_dict['measured field']['B'][sensor][axis]
                    # ,color = colors[i_plot]
                    ,color = "C{}".format(i_plot)
                    # ,alpha = 0.8
                    ,linewidth = 2
                    )

                    ax_i.plot(theta_plot+2,input_dict['calculated field'][sensor][axis]
                    # ,color = h[-1].get_color()
                    ,color = [0.4,0.4,0.4]
                    ,linewidth = 3
                    ,alpha = 0.6
                    # ,linestyle = '--'
                    )

            

            ax_i.set_ylabel('{} {}'.format(sensor, axis))
            plt.grid()
            # if i_plot!=3:
            #      plt.setp(ax_data_plots[i_plot].get_xticklabels(), visible=False)
            i_plot-=1

    # lgd = pltax.legend(bbox_to_anchor=(1.35, 1.3))
    # pltax.legend(loc='best')
    # pltax.grid()
    # xmin,xmax = pltax.get_xlim()
    # ymin,ymax = pltax.get_ylim()
    # xloc = xmin - 1*(xmax-xmin)
    # yloc = ymax + 0.05*(ymax-ymin)
    # pltax.text(xloc,yloc,textstring)
    ax_bottom.set_xlabel('$\\theta/\pi$')
    # pltax.set_ylabel('Magnetic field (pT)')

#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                                CheckValid
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def CheckValid(source_location_in,rotor_dimensions):
    # First check if source lies inside the volume:
    valid_basic = 0
    valid = 0
    source_location = GetPolar(source_location_in)

    if len(source_location['location'])>2: # In cases when z is being used
        if ((2*npfuncs['abs'](source_location['location'][2])<=\
         rotor_dimensions['height']) and
         (source_location['location'][0]<=rotor_dimensions['outer_radius'])and
         (source_location['location'][0]>=rotor_dimensions['hole_radius'])):
            # This now means the thing is inside the rotor boundary
            valid_basic = 1
            # print('Source matches z and r')
#            print(valid)
    else: # In cases when z is not being used
        if source_location['location'][0]<=rotor_dimensions['outer_radius']:
            valid_basic = 1
            # print('no z given, source matches r')
#     print(source_location)
    # remove negative r
    if source_location['location'][0]<0:
        valid_basic  = 0
        # print('r negative')
    # remove theta outside 2pi
    if source_location['location'][1]>360:
        valid_basic = 0
        # print('theta greater than $2\pi$')
    if source_location['location'][1]<0:
        # print('theta negative')
        valid_basic = 0

    # print('Basic validity  = %0.0f'%(valid_basic))
    # Now let's check if the source is in the rotor's anulus or bar
    if valid_basic:
        # Check if its in the anulus
        # print(source_location['location'][0])
        # print(rotor_dimensions['inner_radius'])
        if source_location['location'][0]>=rotor_dimensions['inner_radius']:
            valid = 1
            # print("inside anulus")
        # If not, check if its in the bar
        elif npfuncs['abs'](2*source_location['location'][0]*\
         npfuncs['sindeg'](source_location['location'][1]))<=\
         rotor_dimensions['bar_width']:
            valid = 1
    #         print("in the bar")
    # print('Final validity  = %0.0f'%(valid))
    return valid

#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                                GetCartesian
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def GetCartesian(loc_in):
    loc = ReadAndParseLocation(loc_in)
    loc_cartesian = {'coordinate':'cartesian'}
    if loc['coordinate'] == 'spherical':
        r = loc['location']['r']
        theta = loc['location']['theta']
        if 'phi' in loc['location']:
            phi = loc['location']['phi']
            loc_cartesian['location'] =\
             r([npfuncs['cosdeg'](theta)*npfuncs['cosdeg'](phi),
             npfuncs['sindeg'](theta)*npfuncs['cosdeg'](phi),
             npfuncs['sindeg'](phi)])
        else:
            loc_cartesian['location'] = r([npfuncs['cosdeg'](theta),
                                    npfuncs['sindeg'](theta)])
    elif loc['coordinate'] == 'polar':
        r = loc['location']['r']
        theta = loc['location']['theta']
        if 'z' in  loc['location']:
            z = loc['location']['z']
            loc_cartesian['location'] = ([r*npfuncs['cosdeg'](theta),
                                    r*npfuncs['sindeg'](theta),
                                    z])
        else:

            loc_cartesian['location'] = ([r*npfuncs['cosdeg'](theta),
                                    r*npfuncs['sindeg'](theta)])
    elif loc['coordinate'] == 'cartesian':
        if 'z' in loc['location']:
            loc_cartesian['location'] = [loc['location']['x'],
                                         loc['location']['y'],
                                         loc['location']['z']]
        else:
            loc_cartesian['location'] =\
             [loc['location']['x'],loc['location']['y']]

    else:
        print(['Location coordinate must be either spherical (r,theta,phi),',\
        ' polar (r,theta,z), or cartesian (x,y,z)'])
    return loc_cartesian

#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                                 GetPolar
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def GetPolar(loc_in):
    loc = ReadAndParseLocation(loc_in)
    loc_polar={'coordinate':'polar'}
    if loc['coordinate'] == 'spherical':
        if 'phi' in loc['location']:
            loc_polar['location'] =\
             [loc['location']['r']*npfuncs['cosdeg'](loc['location']['phi']),\
             loc['location']['theta'],\
             loc['location']['r']*npfuncs['sindeg'](loc['location']['phi'])]
        else:
            loc_polar['location'] = \
             [loc['location']['r'],loc['location']['theta']]
    elif loc['coordinate'] =='polar':
        if 'z' in loc['location']:
            loc_polar['location'] = [loc['location']['r'],
                         loc['location']['theta'],
                         loc['location']['z']]
        else:
            loc_polar['location'] = \
             [loc['location']['r'],loc['location']['theta']]
    elif loc['coordinate'] == 'cartesian':
        if 'z' in loc['location']:
            loc_polar['location'] =\
             [npfuncs['sqrt'](loc['location']['x']**2+loc['location']['y']**2),\
             npfuncs['atan2'](loc['location']['y'],loc['location']['x']),\
             loc['location']['z']]
        else:
            loc_polar['location'] =\
             [npfuncs['sqrt'](loc['location']['x']**2+loc['location']['y']**2),\
             npfuncs['atan2'](loc['location']['y'],loc['location']['x'])]
    else:
        print(['Location coordinate must be either spherical (r,theta,phi)',\
        ', polar (r,theta,z), or cartesian (x,y,z)'])
    return loc_polar

#..........................................................................
#``````````````````````````````````````````````````````````````````````````
#                            ReadAndParseLocation
#..........................................................................
#``````````````````````````````````````````````````````````````````````````
def ReadAndParseLocation(loc):
    # Give a dict containing coordinate system and a tuple of values for the
    # vector in that coordinate system
    # Works in 3 D if its a 3-tuple
    # This function is used by other functions and most times not directly
    # by the user.
    # The user should be able to input all locations in a dict with coordinate
    #  type and location as a tuple.
    # print(loc)
    if loc['coordinate'] == 'spherical':
        loc_return = {'coordinate':'spherical'}
        r = loc['location'][0]
        theta = loc['location'][1]
        loc_return['location'] = {'r':r, 'theta':theta}
        if len(loc['location'])>2:
            phi = loc['location'][2]
            loc_return['location']['phi'] = phi
    elif loc['coordinate'] =='polar':
        loc_return = {'coordinate':'polar'}
        r = loc['location'][0]
        theta = loc['location'][1]
        loc_return['location'] = {'r':r, 'theta':theta}
        if len(loc['location'])>2:
            z = loc['location'][2]
            loc_return['location']['z'] = z
    elif loc['coordinate'] == 'cartesian':
        loc_return = {'coordinate':'cartesian'}
        x = loc['location'][0]
        y = loc['location'][1]
        loc_return['location'] = {'x':x, 'y':y}
        if len(loc['location'])>2:
            z = loc['location'][2]
            loc_return['location']['z'] = z
    else:
        print(['Location coordinate must be either spherical (r,theta,phi),',
        ' polar (r,theta,z), or cartesian (x,y,z)'])
    return loc_return
