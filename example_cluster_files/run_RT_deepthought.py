#!/bin/env python

"""
Runs the RT code within the context of a deepthought job.

author: @arjunsavel
"""

from RT import *
import march
from tqdm import tqdm
import os
import pdb
import sys
from sklearn.model_selection import ParameterGrid

# the only thing that gets incremented is already_done — the last entry
# in the slurm script.


def construct_phases():
    """
    Helper function. Returns 50 phase points normalized between -15.936 and 15.936 as strings. Negative
    phases are cast as 360 + phase.
    """

    def custom_piecewise(x, condition1, condition2, condition3, func1, func2, func3):
        """
        Another helper function. Makes a custom piecewise function based on three functions
        and their conditions.
        """
        y = np.empty(np.shape(x)[0])
        y[condition1] = func1(x[condition1])
        y[condition2] = func2(x[condition2] - np.max(x[condition1])) + np.max(y[condition1])
        y[condition3] = func3(x[condition3] - np.max(x[condition2])) + np.max(y[condition2])
        return y
    
    def normalize_to_pm_val(val):
        """
        Normalizes to +/- whatever value is input.
        """
        return 2 * val * (phase_evals - np.amin(phase_evals)) / (np.amax(phase_evals) - np.amin(phase_evals)) - val
    
    phase_x = np.arange(1, 53, 1, dtype=np.float)

    phase_evals = custom_piecewise(phase_x, (phase_x < 20), (phase_x >=20) & (phase_x <= 30), 
                  (phase_x > 30),
                lambda x: .3*x, 
                  lambda x: x, 
                  lambda x: .3 * x)

    max_x_coord = 15.936
    phase_evals = normalize_to_pm_val(phase_evals, max_x_coord)

    phase_evals[phase_evals < 0.0] = 360 + phase_evals[phase_evals < 0.0]
    
    return phase_evals.astype(str)[1:-1]



def run_RT_deepthought():
    """
    Main function for running the RT code within deepthought.

    When called with main: the first val is the number of the chunk,
    the second is whether doppler is on or off, 
    and the third is how many spectra have already been computed.
    """

    os.chdir('RT_3D_Transmission_Code')


    input_val = sys.argv[1]
    doppler = bool(eval(sys.argv[2])) # 1 is Doppler on, 0 is doppler off

    # the "already done" parameter must be updated on each job array submission
    already_done = eval(sys.argv[3]) 
    spectrum_number = (eval(input_val) + already_done) // 113 # wrap at 113 because I have 113 chunks.
    wave_chunk = (eval(input_val) + already_done) % 113


    ######################################### Set up parameter grid ###################################
    drag_timescales = ['1e7']
    clouds_list = [False]
    phases = construct_phases()
    models = ['_deep_rcb']
    condensations = [False]
    condense_species = ['Fe']

    parameter_list = list(ParameterGrid({'drag_timescale': drag_timescales,
       'model': models,
        'condensation':condensations,
        'condense_species':condense_species,
        'clouds':clouds_list,
                                  'phase':phases}))

    ##################################### End parameter grid setup ####################################
    
    
    # now access parameters for this spectrum!
    parameters = parameter_list[spectrum_number]
    drag_timescale = parameters['drag_timescale']
    phase = parameters['phase']
    model = parameters['model']
    condensation = parameters['condensation']
    clouds = parameters['clouds']
    condense_species = parameters['condense_species']


    ################# Below is mostly specific to Arjun's grid / C code wrapper #######################


    # Select TP file and decide whether we need to create a new one
    if eval(phase):
        t_p_filename = f'init_t_p_3D_WASP76{model}_{drag_timescale}.dat_phase_{phase}_inc_0.0.txt'
    else:
        t_p_filename = f't_p_3D_WASP76{model}_{drag_timescale}.dat'
    t_p_path = 'Data/t_p_profiles/' + t_p_filename

    # create TP profile if it doesn't exist
    if not os.path.exists(t_p_path): 
      os.chdir('Data/t_p_profiles')
      planet_name = f't_p_3D_WASP76{model}_{drag_timescale}.dat'
      altitudes, latitudes, longitudes = march.read_t_p_file(planet_name)
      outdir = os.getcwd()
      phase_list = [eval(phase)]
      march.rotate_planet(phase_list, 
                          planet_name, 
                          latitudes, 
                          longitudes, 
                          outdir, 
                          parallel=False)
      os.chdir('../../')
    if eval(phase):
        altitudes, latitudes, longitudes = march.read_t_p_file(t_p_path)
        a = 0.0330 * 1.496e+11  # meters
        xcen, ycen = march.get_planet_coords(eval(phase), a)
        print(phase)
        march.planet_intensity(xcen, ycen,
                       altitudes,
                       latitudes,
                       "Data/",
                       eval(phase))

    # input is all now changed in run_single_chunk 
    # run simulation
    further_outdir = f'spectrum_number_{spectrum_number}'
    obj = RT(doppler=doppler,
           condensation=condensation,
           change_elems=change_elems,
           further_outdir=further_outdir,
           phase=eval(phase),
             TP_file=t_p_path,
           clouds=clouds)


    """
    create a "further outdir" so that multiple simulations can
    run the same chunk at the same time without overwriting
    one another.
    """
    further_outdir = f'spectrum_number_{spectrum_number}'
    obj.run_single_chunk(wave_chunk)

    """
    cleaning will have to happen later.
    Settings for each file can be replicated
    with the same parameter grid.
    """                
if ( __name__ == '__main__' ):
  run_RT_deepthought()
