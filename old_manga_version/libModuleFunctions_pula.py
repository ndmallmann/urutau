#Programmer: Nicolas Dullius Mallmann
#Functions to use with the control

import os
import numpy as np
import inspect as ins
import libManga as me
import libConfig as lc
import libCreateMegacube as cm
import astropy.io.fits as fits
import libStarlightCaller as sc

def extraction(config, job_order, job_list):
    """
    Wrapper for extraction function in libMaNGA.
    """

    #Get the parameters
    work_path = config["master"]["root_path"]

    input_data = config["extraction"]["input"]
    drp = config["extraction"]["drp_file_path"]
    stn = config["extraction"]["signal_to_noise"]
    wdw = config["extraction"]["window_range"]
    ext = config["extraction"]["spectra_extension"]
    cor = config["extraction"]["correct_extinction"]


    ########### PARTE MODIFICADA (04/28/2020) ##########################
    #Ignores this run if the MEGACUBE already exists
    plateifu = os.path.basename(input_data).split("-LINCUBE")[0]
    plateifu = plateifu.split("manga-")[-1]
    mega_output_file = "manga-" + plateifu + "-MEGA.fits"
    mega_output_path = os.path.join(work_path, plateifu, "Megacube", mega_output_file)

    print("****** VERIFICAR!!!!!")
    print("****** -ROOT PATH: {0}".format(work_path))
    print("****** -PLATEIFU: {0}".format(plateifu))
    print("****** -CAMINHO DO MEGACUBE: {0}".format(mega_output_path))
    print("****** -ARQUIVO ORIGINAL: {0}".format(input_data))

    if os.path.exists(mega_output_path):
        print("--- THE FILE EXISTS, TRYING TO END RUN")
        #Add input variable to (or modify variables of) the next job
        if (job_order + 1) < len(job_list):
            next_job = job_list[job_order + 1]
            config.set(next_job, "input", [extracted_path, input_data])
        return

    ########### FIM DA PARTE MODIFICADA (04/28/2020) ###################

    #Run function
    extracted_path = me.extraction(input_data, drp, work_path, stn, wdw, ext, cor)

    #Add input variable to (or modify variables of) the next job
    if (job_order + 1) < len(job_list):
        next_job = job_list[job_order + 1]
        config.set(next_job, "input", [extracted_path, input_data])

    #OPTIONAL: Save used configuration
    config.write(os.path.join(os.path.dirname(extracted_path), "configEXT.ini"))



def synthesis(config, job_order, job_list):
    """
    Wrapper for the synthesis function in libStarlightCaller.
    """

    #Get the parameters
    param = {}
    param["executable_path"] = config["synthesis"]["executable_path"]
    directory = os.path.join(config["master"]["root_path"],
                os.path.basename(os.path.basename(os.path.split(config["synthesis"]["input"][0].rstrip("/"))[0])))
    param["synthesized_spectra_path"] = os.path.join(directory, "Sinthesized")
    param["executable"] = config["synthesis"]["executable"]
    param["bases_dir"] = config["synthesis"]["bases_dir"]
    param["observed_spectra_path"] = config["synthesis"]["input"][0]
    param["mask_path"] = config["synthesis"]["mask_path"]    
    param["seed"] = config["synthesis"]["seed"]
    param["llow_sn"] = config["synthesis"]["llow_sn"]
    param["lupp_sn"] = config["synthesis"]["lupp_sn"]
    param["olsyn_ini"] = config["synthesis"]["olsyn_ini"]
    param["olsyn_fin"] = config["synthesis"]["olsyn_fin"]
    param["odlsyn"] = config["synthesis"]["odlsyn"]
    param["fscale_chi2"] = config["synthesis"]["fscale_chi2"]
    param["fit_fxk"] = config["synthesis"]["fit_fxk"]
    param["iserrspecavailable"] = config["synthesis"]["iserrspecavailable"]
    param["isflagspecavailable"] = config["synthesis"]["isflagspecavailable"]
    param["starlight_config"] = config["synthesis"]["starlight_config"]
    param["master_base_file"] = config["synthesis"]["master_base_file"]
    param["mask_file"] = config["synthesis"]["mask_file"]
    param["reddening_law"] = config["synthesis"]["reddening_law"]
    param["velocity_recession"] = config["synthesis"]["velocity_recession"]
    param["velocity_dispersion"] = config["synthesis"]["velocity_dispersion"]
    param["spectra_extension"] = config["synthesis"]["spectra_extension"]
    param["output_extension"] = config["synthesis"]["output_extension"]
    param["processors"] = config["synthesis"]["processors"]


    ########### PARTE MODIFICADA (04/28/2020) ##########################

    #Try to create the directory
    if os.path.exists(param["synthesized_spectra_path"]):
        # Directory already exists (skipping starlight run)
        #Add input variable to (or modify variables of) the next job
        if (job_order + 1) < len(job_list):
            next_job = job_list[job_order + 1]
            config.set(next_job, "input", [param["synthesized_spectra_path"], config["synthesis"]["input"][1]])
        return
    else:
        os.makedirs(param["synthesized_spectra_path"])

    ########### PARTE MODIFICADA (04/28/2020) ##########################

    #Run the function
    sc.run_starlight(**param)

    #Add input variable to (or modify variables of) the next job
    if (job_order + 1) < len(job_list):
        next_job = job_list[job_order + 1]
        config.set(next_job, "input", [param["synthesized_spectra_path"], config["synthesis"]["input"][1]])

    #OPTIONAL: Save used configuration
    config.write(os.path.join(os.path.dirname(param["synthesized_spectra_path"]), "configSYNTH.ini"))

def analysis(config, job_order, job_list):

    #Get the parameters
    c = 299792.458
    H0 = 73.0

    synth_dir = config["analysis"]["input"][0] + "/"
    cube_file = config["analysis"]["input"][1]
    plateifu = fits.getheader(cube_file, 0)["plateifu"]

    directory = os.path.join(os.path.join(config["master"]["root_path"], os.path.basename(os.path.split(synth_dir.rstrip("/"))[0])))
    output_file_path = os.path.join(directory, "Megacube")

    #Try to create the directory
    try:
        os.makedirs(output_file_path)
    except:
        pass

    output_file_name = os.path.join(output_file_path, "manga-" + plateifu + "-MEGA.fits")
    config.set("analysis", "SaveFilePath", output_file_name)

    if os.path.exists(output_file_name):
        return

    arq = fits.getdata(config["analysis"]["drp_file_path"], 1)
    zdist = arq[np.where(arq["plateifu"] == plateifu)]["nsa_zdist"][0]
    config.set("analysis", "gal_dist", (c*zdist)/H0)

    cube_name = os.path.basename(cube_file)
    cubes_dir = os.path.dirname(cube_file) + "/"

    config.set("analysis", "cubes_dir", cubes_dir)
    config.set("analysis", "original_cube", cube_file)
    config.set("analysis", "synth_dir", synth_dir)
    config.set("analysis", "plateifu", plateifu)

    config.write(os.path.join(directory, "configANL.ini"))

    #Run the function
    cm.CuboStarlightFits(config, "analysis")
