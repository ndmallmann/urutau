#Programmer: Nicolas Dullius Mallmann

import math, sys, os
import copy as cp
import glob as gl
import noisy as ns
import numpy as np
import scipy as sp
import astropy.io.fits as fits
import ReddeningCorrections as rc

def extract_point(cube_file_handler, extract_path = "./", point = [0,0], redshift = 0.0, extension = "spec", correct_extinction = False):
    """
        Extract a single point of a MaNGA datacube and saves the data in an ASCII file.

        Parameters:
            cube_file_handler       -> astropy.io.fits file handler of a MaNGA cube (LINCUBE)
            extract_path            -> relative or absolute path to extract the ascii file
            point                   -> list with the X and Y coordinates of a point of the cube
            redshift                -> redshift used to correct the wavelenght of the data
            extension               -> extension of the extracted spectrum
            correct_extinction      -> if set to True, apply ccm correction

        Warnings:
            -This function do not test if the cube exists;
            -This function do not test if the extraction path directories exist;
            -This function do not test if the point is inside the cube limits.
    """

    #Constants section
    _PRIMARY_INDEX_ = "PRIMARY"
    _FLUX_INDEX_ = "FLUX"
    _IVAR_INDEX_ = "IVAR"
    _MASK_INDEX_ = "MASK"
    _WAVE_INDEX_ = "WAVE"
    _PLATEIFU_ = "PLATEIFU"
    _FORMAT_ = ("%.2f","%.3e", "%.3e", "%i")
    _REDSHIFT_CORRECTION_ = 1.0/(1.0+redshift)
    _EXTRACT_PREFIX_ = cube_file_handler[_FLUX_INDEX_].header[_PLATEIFU_]
    _NAME_ = "{0}_{1}_{2}.{3}"
    _IGNORE_ = "ignore"
    _EBV_ = "ebvgal"
    _WARN_ = "warn"
    _SPACES_ = "   "
    _X_ = 0
    _Y_ = 1

    #Reads all the relevant data (flux, inverse variance, mask and wave) for the point (X,Y)
    flux = cube_file_handler[_FLUX_INDEX_].data[:, point[_Y_], point[_X_]]
    ivar = cube_file_handler[_IVAR_INDEX_].data[:, point[_Y_], point[_X_]]
    mask = cube_file_handler[_MASK_INDEX_].data[:, point[_Y_], point[_X_]]
    wave = cube_file_handler[_WAVE_INDEX_].data

    #Get all the good flux data (non-negative values)
    good_flux = (flux >= 0.0)

    #Get all the good inverse variance data (positives)
    good_ivar = (ivar > 0.0)

    #Ignore warning about division by zero (don't worry the extraction function handles the division by zero)
    np.seterr(divide=_IGNORE_)

    #Calculates the error (if good ivar, return the square root of the inverse of ivar, else, return zero)
    error = np.where(good_ivar, np.sqrt(1.0/ivar), 0.0)

    #Restore the warnings about division by zero
    np.seterr(divide=_WARN_)

    #Correct for extinction if requested
    if correct_extinction:
        ebv = cube_file_handler[_PRIMARY_INDEX_].header[_EBV_]
        w, f, err = rc.deredden(np.column_stack([wave, flux, error]), "ccm", ebv)
        wave = w
        flux = f
        error = err

    #Updates the mask data with the good flux data and good ivar data positions
    mask = np.where(good_flux * good_ivar, mask, 1.0)

    #Redshift correction
    wave = wave*_REDSHIFT_CORRECTION_

    #New wave, flux, error and mask data
    new_wave = np.arange(start = math.ceil(wave[0]), stop = math.ceil(wave[-1]), step = 1.0)
    new_flux = sp.interp(new_wave, wave, flux)
    new_error = sp.interp(new_wave, wave, error)
    new_mask = np.ceil(sp.interp(new_wave, wave, mask))
    

    #Count the number of useful points
    n = 0
    for num in new_mask:
        if num == 0:
            n+=1

    #If at least 1000 points aren't crap, save the spec ## RR changed to 500
    #if n < 500: print('too less useful points: ',n)
    if n > 500:
        #Organize all the data to save afterwards
        data = np.column_stack((new_wave, new_flux, new_error, new_mask))

        #Create the output file name and its path to save the data
        file_name = _NAME_.format(_EXTRACT_PREFIX_, point[_Y_], point[_X_], extension)
        save_path = os.path.join(extract_path, file_name)

        #Save all the relevant data
        np.savetxt(save_path, data, fmt=_FORMAT_, delimiter=_SPACES_)


def extraction(cube_file_path, drp_file_path, extracted_spectra_path = "./", signal_to_noise = 10.0, window_range = [5650, 5750], spectra_extension = "spec", correct_extinction = False):
    """
        Extract MaNGA datacube to be used with the Starlight Code (by Cid Fernandes)

        This function extracts all points of one cube file with a minimal signal to noise and saves them in separate ascii files.
        This function uses a MaNGA drpall file to get the object's redshift.

        Parameters:
            cube_file_path          -> MaNGA datacube's path (LINCUBE)
            drp_file_path           -> drp datacube's path
            extract_spectra_path    -> relative or absolute path to extract the ascii file
            signal_to_noise         -> minimal signal to noise to extract a point
            window_range            -> window of wavelenght to be used to calculate the signal to noise mask
            spectra_extension       -> extension of the extracted spectra
            correct_extinction      -> if set to True, apply ccm correction

        Warnings:
            -This function do not test if the cube exists;
            -This function do not test if the drpall exists;
            -This function do not test if the extraction path directories exist;
            -This function do not test if the window_range is good or not.
    """

    print(cube_file_path,signal_to_noise)
    #Constants section
    _SN_MASK_EXTENSION_NAME_ = "S2N_MASK"
    _REDSHIFT_ = "nsa_z"
    _PRIMARY_HDU_ = "PRIMARY"
    _PLATEIFU_ = "PLATEIFU"
    _BAD_MASK_VALUE_ = 0.0
    _BUTTERWORTH_ORDER_ = 3
    _BUTTERWORTH_RANGE_ = 0.3
    _MODE_NOISE_ = "r"
    _DRP_TABLE_ = 1
    _FIRST_ = 0
    #End of the contants section

    #Create a file handler for the cube
    cube_file = fits.open(cube_file_path)
    cube_aux = fits.open(cube_file_path)

    #Pass the butterworth filter
    cube_aux["FLUX"].data = _butterworth_cube(cube_aux, _BUTTERWORTH_ORDER_, _BUTTERWORTH_RANGE_)

    #Reads the plate and ifu number
    cube_plateifu = cube_file[_PRIMARY_HDU_].header[_PLATEIFU_]

    extracted_spectra_path = os.path.join(extracted_spectra_path, cube_plateifu, "Extracted")
    if os.path.exists(extracted_spectra_path):
        #return extracted_spectra_path
        pass
    else:
        os.makedirs(extracted_spectra_path)

    #Reads the drp file to get extra info about the cube
    drp_data = fits.getdata(drp_file_path, _DRP_TABLE_)

    #Get the index for the info about the cube, the cube info and its redshift
    cube_info_index = np.where(drp_data[_PLATEIFU_] == cube_plateifu)[_FIRST_]
    cube_info = drp_data[cube_info_index]
    cube_redshift = cube_info[_REDSHIFT_] 

    #Open a noise object to calculate the mask of signal to noise
    noise = ns.noise(cube_aux)

    #Calculate the mask of the signal to noise ratio
    noise_mask = noise.S_TO_N_CUBE(do_mask = True, window = window_range, SN_tresh = signal_to_noise, mode = _MODE_NOISE_)
    noise_mask = noise_mask[_SN_MASK_EXTENSION_NAME_].data[_FIRST_]

    #Get all the points with signal to noise above the threshold
    list_of_points = np.column_stack(np.where(noise_mask > _BAD_MASK_VALUE_))
#    print(noise_mask)
    #Use the function to extract all the points
    for point_to_extract in list_of_points:
        extract_point(cube_file_handler = cube_aux, extract_path = extracted_spectra_path, point = [point_to_extract[1], point_to_extract[0]], redshift = cube_redshift, extension = spectra_extension, correct_extinction = correct_extinction)

    #Append butterworth extension to original cube and save it inside the extracted_spectra_path
    hdu = fits.hdu.ImageHDU()
    hdu.data = cube_aux["FLUX"].data
    hdu.header.append(("AUTHOR", "Nicolas D. Mallmann"))
    hdu.header.append(("FILTER", "Butterworth - Highpass"))
    hdu.header.append(("FORDER", 3))
    hdu.header.append(("FRANGE", 0.3))
    hdu.header.append(("PLATEIFU", cube_plateifu))
    hdu.header.append(("EXTNAME", "FLUX_B"))

    cube_file.append(hdu)

    path_cube = os.path.join(extracted_spectra_path, os.path.basename(cube_file_path).rstrip(".fits") + "_BW3.fits")
    cube_file.writeto(path_cube,overwrite=True)

    return extracted_spectra_path


def _butterworth_cube(manga_cube, N, W):

    z,y,x = manga_cube["FLUX"].data.shape
    circular = np.zeros((x,y))

    for i in range(0, x):
        for j in range(0, y):
            if x%2 == 0:
                if y%2 == 0:
                    circular[i,j] = np.sqrt(((i-x/2+0.5)**2. + (j-y/2+0.5)**2.))
                else:
                    circular[i,j] = np.sqrt(((i-x/2+0.5)**2. + (j-y/2)**2.))
            else:
                if y%2 == 0:
                    circular[i,j] = np.sqrt(((i-x/2)**2. + (j-y/2+0.5)**2.))
                else:
                    circular[i,j] = np.sqrt(((i-x/2)**2. + (j-y/2)**2.))

    circular = circular/circular.max()
    H = np.zeros_like(circular)
    H = 1. - 1./(1. + (W/circular)**(2.*N))
    aux = np.roll(H, - int(y/2), axis=0)
    H = np.roll(aux, - int(x/2), axis=1)
    data = cp.deepcopy(manga_cube["FLUX"].data)

    for k in range(0, z):
        temp = np.fft.fft2(data[k,:,:])
        temp = temp*H
        data[k,:,:] = np.fft.ifft2(temp).real

    return data
