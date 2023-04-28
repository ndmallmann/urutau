#Programmer: Nicolas Dullius Mallmann
#Name: starlight_call
#Status: Development

#Modules section
import glob as gl
import numpy as np
import os.path as pth
import threading as th
import subprocess as sp

#Function to create file grids
def run_starlight(executable_path, synthesized_spectra_path, executable, bases_dir, \
                  observed_spectra_path, mask_path, seed, llow_sn, lupp_sn, \
                  olsyn_ini, olsyn_fin, odlsyn, fscale_chi2, fit_fxk, iserrspecavailable, \
                  isflagspecavailable, starlight_config, master_base_file, mask_file,
                  reddening_law, velocity_recession, velocity_dispersion, \
                  spectra_extension, output_extension, processors):

    #Constants section
    _SPECTRA_FINDER_ = pth.join(observed_spectra_path, "*{0}".format(spectra_extension))
    _FORMAT_ = "{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n{7}\n{8}\n{9}\n{10}\n{11}\n{12}\n{13}\n{14}\n"
    _LINE_FORMAT_ = "{0} {1} {2} {3} {4} {5} {6} {7}\n"
    _NAME_FORMAT_ = "grid_{0}.{1}"
    _GRID_EXTENSION_ = "inp"
    _FIRST_ = 0
    _BAR_ = "/"
    #End of constants section

    #Change the directories to absolute values
    executable_path = pth.abspath(pth.expanduser(executable_path)) + _BAR_
    synthesized_spectra_path = pth.abspath(pth.expanduser(synthesized_spectra_path)) + _BAR_
    bases_dir = pth.abspath(pth.expanduser(bases_dir)) + _BAR_
    observed_spectra_path = pth.abspath(pth.expanduser(observed_spectra_path)) + _BAR_
    mask_path = pth.abspath(pth.expanduser(mask_path)) + _BAR_
    grids_path = observed_spectra_path

    #Get list of spectra
    list_of_spectra = gl.glob(_SPECTRA_FINDER_)



    #Separate the list into a lot of lists
    if processors > 1:
        size_list = len(list_of_spectra)
        division_list = np.arange(1,processors)*(size_list/processors)
        list_of_spectra = np.split(list_of_spectra, division_list.astype(int))
    else:
        list_of_spectra = [list_of_spectra]

    list_of_runners = []

    #Create all the grids files
    for i in range(0, len(list_of_spectra)):
            # Opening grid file
            file_name = _NAME_FORMAT_.format(i+1, _GRID_EXTENSION_)
            file_path = pth.join(grids_path, file_name)
            if not pth.exists(file_path):
                save=open(file_path, 'w')

                # Writing the common parameters
                line = _FORMAT_.format(len(list_of_spectra[i]), bases_dir, observed_spectra_path, \
                                           mask_path, synthesized_spectra_path, seed, llow_sn, lupp_sn, \
                                           olsyn_ini, olsyn_fin, odlsyn, fscale_chi2, fit_fxk, \
                                           iserrspecavailable, isflagspecavailable)
                save.write(line)


                # Write the rest of the file (line with the input and output spectra names)
                for j in list_of_spectra[i]:
                    j = pth.basename(j)
                    save.write(_LINE_FORMAT_.format(j, starlight_config, master_base_file, mask_file, \
                                                    reddening_law, velocity_recession, velocity_dispersion, \
                                                    j.replace(spectra_extension, output_extension)))

                # Close grid file
                save.close()

                #Create a consumer to run the starlight
                list_of_runners.append(runner(executable_file_path = executable_path, grid_path = file_path, executable_file_name = executable))
                list_of_runners[-1].start()


    #Wait for all the runners to end
    for runner_element in list_of_runners:
        runner_element.join()

    
    
#Create a class of consumers to run the starlight
class runner(th.Thread):
    def __init__(self, executable_file_path, grid_path, executable_file_name):
        th.Thread.__init__(self)
        self.executable_path = executable_file_path
        self.grid_path = grid_path
        self.executable = executable_file_name

    #Run starlight
    def run(self):

        #Open a grid file to be used as a standard input for the starlight
        grid = open(self.grid_path)

        #Mount the arguments
        arguments = (pth.join(self.executable_path, self.executable))

        #Run starlight
        prog = sp.Popen(args = arguments, stdin = grid, cwd = self.executable_path)

        #Wait for the program to end
        prog.wait()

