#PROGRAMMER: Rogerio Riffel
#MODIFIED BY: Nicolas D. Mallmann

import re, glob

import noisy as ns
import numpy as np
import os.path as pth
import libConfig as conf
import astropy.io.fits as pf

from pylab import *

from ReadStarlightParameters import starlightPars
from ReadStarlightParameters import popVectors
from ReadStarlightParameters import StSyntesis


#This is just a list o constants used by the function to read the data
#You can change the string VALUES if you need (they MUST correspond to the OPTION NAMES inside the configuration object/file)

_original_cube_ = "original_cube"
_Olsyn_ini_ = "Olsyn_ini"
_Olsyn_fin_ = "Olsyn_fin"
_DeltaLamb_ = "DeltaLamb"
_StOutFromList_ = "StOutFromList"
_DataList_ = "DataList"
_BaseFile_ = "BaseFile"
_gal_D_ = "gal_dist"
_outdir_ = "synth_dir"
_GalaxyShortName_ = "plateifu"
_extension_ = "extension"
_IsAGNComp_ = "IsAGNComp"
_OnlnyFC_ = "OnlnyFC"
_BinFCVec_ = "BinFCVec"
#_BinFCVecLab_ = "BinFCVecLab"
_BinHDVec_ = "BinHDVec"
#_BinHDVecLab_ = "BinHDVecLab"
_BinPopVec_ = "BinPopVec"
#_BinPopVecLab_ = "BinPopVecLab"
#_BinPopVecMassLab_ = "BinPopVecMassLab"
_BinSFR_ = "BinSFR"
#_BinSFRLabs_ = "BinSFRLabs"
_NormFac_ = "NormFac"
_filein_ = "FileInput"
_save_file_path_ = "SaveFilePath"
_signal_noise_list_ = "SN_List"
_window_range_ = "WindowRange"
_young_Mage_ = "YoungMage"
_doYoungMage_ = "doYoungMage"

###########################################################################################
def GetCubesPars(config, category_name = "data"):

    #Read the configuration file
    data = config[category_name]

    #Making Header
    headerKeys = []
    headerComments = []
    if data[_doYoungMage_]:
        extraheaderkey=[['Av','Optical extinction'],\
        ['Mage_L','Mean age light weigthed'],\
        ['Mage_M','Mean age mass weigthed'],\
        ['MZ_L','Mean metalicity light weigthed'],\
        ['MZ_M','Mean metalicity mass weigthed'],\
        ['Mage_L'+str(data[_young_Mage_]/1E9)+'G','Mage_L for ages less than: '+str(data[_young_Mage_]/1E9)+'G'],\
        ['Mage_M'+str(data[_young_Mage_]/1E9)+'G','Mage_M for ages less than: '+str(data[_young_Mage_]/1E9)+'G'],\
        ['MZ_L'+str(data[_young_Mage_]/1E9)+'G','MZ_L for ages less than: '+str(data[_young_Mage_]/1E9)+'G'],\
        ['MZ_M'+str(data[_young_Mage_]/1E9)+'G','MZ_M for ages less than: '+str(data[_young_Mage_]/1E9)+'G'],\
        ['Mstar','Present mass in stars (M_sun, M* from starligh)'],\
        ['Mpross','Mass that has been processed in stars (~2xMstar)' ],\
        ['F_Norm','Normalization flux in input units'],\
        ['Sigma_star','Stellar dispersion (see starlight manual)'],\
        ['vrot_star','Stellar rotation (see starlight manual)'],\
        ['Adev','Precentage mean deviation (see manual)'],\
        ['ChiSqrt','ChiSqrt/Nl_eff (see manual)'],\
        ['SNR','SNR on normalization window']]
    else:
        extraheaderkey=[['Av','Optical extinction'],\
        ['Mage_L','Mean age light weigthed'],\
        ['Mage_M','Mean age mass weigthed'],\
        ['MZ_L','Mean metalicity light weigthed'],\
        ['MZ_M','Mean metalicity mass weigthed'],\
        ['Mstar','Present mass in stars (M_sun, M* from starligh)'],\
        ['Mpross','Mass that has been processed in stars (~2xMstar)' ],\
        ['F_Norm','Normalization flux in input units'],\
        ['Sigma_star','Stellar dispersion (see starlight manual)'],\
        ['vrot_star','Stellar rotation (see starlight manual)'],\
        ['Adev','Precentage mean deviation (see manual)'],\
        ['ChiSqrt','ChiSqrt/Nl_eff (see manual)'],\
        ['SNR','SNR on normalization window']]


    #Control variables
    failFile=0
    failFileAGN=0
    failFileHD=0
    failFileSFR=0
    failFileHDMass=0


    #Add extra data if the object contains an AGN
    if data[_IsAGNComp_]:
            aux=data[_BinFCVec_]
            for w in aux:
                headerKeys.append('FC'+str(w))
                headerComments.append('Featureless cont. F_l '+str(aux[w][0])+' & '+str(aux[w][1]))
                failFileAGN += 1
            if not data[_OnlnyFC_]:
                aux=data[_BinHDVec_]
                for w in aux:
                    headerKeys.append('BB'+str(w))
                    headerComments.append('Planck function T: '+str(aux[w][0])+' & '+str(aux[w][1]))
                    failFileHD += 1
    aux=data[_BinPopVec_]
    for w in aux:
           headerKeys.append(w+'_light')
           headerComments.append('Light binned pop. vec '+ '{:.1E}'.format(aux[w][0])+' < age <= '+'{:.1E}'.format(aux[w][1]))
           failFile += 1
    for w in aux:
           headerKeys.append(w+'_mass')
           headerComments.append('Mass binned pop. vec '+ '{:.1E}'.format(aux[w][0])+' < age <= '+'{:.1E}'.format(aux[w][1]))
           failFile += 1
    
    aux=data[_BinSFR_]
    for w in aux:
           headerKeys.append(w)
           headerComments.append('SFR for '+ '{:.1E}'.format(aux[w][0])+' < age <= '+'{:.1E}'.format(aux[w][1]))
           failFileSFR += 1
    for w in np.arange(0,len(extraheaderkey)):
           headerKeys.append(extraheaderkey[w][0])
           headerComments.append(extraheaderkey[w][1])

    try:
        #Get the starlight parameters from the file
        pars = starlightPars(data[_filein_])

       # print('doing file: {0} computed with Starlight version {1}'.format(data[_filein_], pars[-1]))

        #Getting some paramenters of interest
        gd = (3.08567758E24) * float(data[_gal_D_])

        if pars[-1]=='V4':
            chi2 = pars[0][0]
            fobs_norm = pars[0][1]
            adev = pars[0][2]
            Mini_tot = pars[0][4]
            Mcor_tot = pars[0][5]
            av = pars[0][6]
            SNR = pars[0][13]
            SumPopVecs = pars[0][15]
            Sigma_star=pars[0][16]
            v_rot=pars[0][12]

            #Making some calculations for the SFR and total mass, see STARLIGHT manual, page: 22
            #Using 1 Mpc = 3.08567758E24 cm
            Mini_t = data[_NormFac_] * (Mini_tot * 4.0*np.pi) * ((gd**2.) / 3.826E33)
            Mcor_t = data[_NormFac_] * (Mcor_tot * 4.0*np.pi) * ((gd**2.) / 3.826E33)

        elif pars[-1]=='V5':
            chi2= pars[0][0]
            fobs_norm= pars[0][3]
            adev= pars[0][6]
            av= pars[0][10]
            SNR=pars[0][16]
            SumPopVecs=pars[0][18]
            Sigma_star=pars[0][19]
            v_rot=pars[0][15]
            #Making some calculations for the SFR and total mass, see STARLIGHT manual, page: 22
            #Using 1 Mpc = 3.08567758E24 cm
            #IMPORTANT: For V5 the values of Mini_t and Mcor_t are taken from the output file
            if float(pars[0][5]) != float(data[_gal_D_]):
                # Testing if distance is OK. If not, using starlight input distance.
                print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                print ("The distance informed in the input file and in the config files are not equal:")
                print ("Distance in Starlight: {0} Distance in ConfigFile: {1}".format(pars[0][5], gal_D))
                print ("For consistency I'm setting value to Starlight output = {0}".format(pars[0][5]))
                print ("Check this value or do not consider all values depending on distance.")
                print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                config.set(category_name, _gal_D_, pars[0][5])
            Mini_t = pars[0][8]
            Mcor_t =pars[0][9]

        #Getting the population vectors paramenters
        popVecs = popVectors(data[_filein_])

        #Separating the population vectors parameters in arrays (normalizing to 100%)
        popx = popVecs[0][:,0] * (100./ popVecs[0][:,0].sum()) #ADD THIS SUM TO A HEADER
        popMi = popVecs[0][:,1] * (100./ popVecs[0][:,1].sum()) #ADD THIS SUM TO A HEADER
        popm = popVecs[0][:,2] * (100./ popVecs[0][:,2].sum())
        popage = popVecs[0][:,3]
        popZ = popVecs[0][:,4]
        components = popVecs[1]

        #Getting the observed and synthetic spectra
        sintRes = StSyntesis(data[_filein_])
        l_obs = sintRes[:,0]
        f_obs = sintRes[:,1] * fobs_norm * data[_NormFac_] 
        f_syn = sintRes[:,2] * fobs_norm * data[_NormFac_]
        wei = sintRes[:,3]


        #########################################################
        #                                                       #
        #    Extracting the Bined population Vectors            #
        #                                                       #
        #########################################################
        PopBin=[]
        aux=data[_BinPopVec_]
        for bins in aux:
            PopBini = (np.sum(popx[(popage > aux[bins][0]) & (popage <= aux[bins][1])]))
            PopBin.append(PopBini)
       # to add the population vectos of the mass fractions below the light fraction
        for mass in  aux:
            PopBiniM = (np.sum(popm[(popage > aux[mass][0]) & (popage <= aux[mass][1])]))
            PopBin.append(PopBiniM)

        #########################################################
        #                                                       #
        #    Extracting mean AGE and Mean Metallicities         #
        #                                                       #
        #########################################################
        l=0
        sumvec=[]
        for component in components:
            if ('BB_') not in component:
                if ('Power') not in component:
                    sumvec.append(l)
            l += 1

        #Here it was renormalized in order to the stellar fraction sum 100%, since only it gives mean SP age.
        meanAgex = np.sum(((100*(popx[sumvec])/np.sum(popx[sumvec]))/100)*np.log10(popage[sumvec]))

        #Here it was renormalized in order to the stellar fraction sum 100%
        meanZx= np.sum(((100*(popx[sumvec])/np.sum(popx[sumvec]))/100)*popZ[sumvec])

        #Here it was renormalized in order to the stellar fraction sum 100%
        meanAgem= np.sum(((100*(popm[sumvec])/np.sum(popm[sumvec]))/100)*np.log10(popage[sumvec]))

        #Here it was renormalized in order to the stellar fraction sum 100%
        meanZm=sum(((100*(popm[sumvec])/np.sum(popm[sumvec]))/100)*popZ[sumvec])
        if data[_doYoungMage_]:
            q=0
            sumvec1g=[]
            for component in components:
                if ('BB_') not in component:
                    if ('Power') not in component:
                        if (popage[q] <= data[_young_Mage_]):
                             sumvec1g.append(q)
                q=q+1
            meanAgex1g=sum(((100*(popx[sumvec1g])/sum(popx[sumvec1g]))/100)*np.log10(popage[sumvec1g])) # here it was renormalized in order to the stellar fraction sum 100%, since only it gives mean SP age.
            meanZx1g= sum(((100*(popx[sumvec1g])/sum(popx[sumvec1g]))/100)*popZ[sumvec1g]) # here it was renormalized in order to the stellar fraction sum 100%

            meanAgem1g= sum(((100*(popm[sumvec1g])/sum(popm[sumvec1g]))/100)*np.log10(popage[sumvec1g])) # here it was renormalized in order to the stellar fraction sum 100%
            meanZm1g=sum(((100*(popm[sumvec1g])/sum(popm[sumvec1g]))/100)*popZ[sumvec1g]) # here it was renormalized in order to the stellar fraction sum 100%


        #########################################################
        #                                                       #
        #            Extracting SFR paramenters                 #
        #                                                       #
        #########################################################

        SFR=[]
        aux=data[_BinSFR_]
        for bins in aux:
            SFR_t = (np.sum(popMi[(popage > aux[bins][0]) & \
                (popage <= aux[bins][1] )]))
            SF = aux[bins][1] - aux[bins][0]
            SFR.append(( (SFR_t/100) * Mini_t)/SF)

        if data[_IsAGNComp_]:
                #########################################################
                #                                                       #
                #              Extracting AGN Components                #
                #                                                       #
                #########################################################

                #Searching for the FC components.
                #Here I'm using the string vector with the component_j names searching for those with Power in the name, wich means the FC
                #Note here that the power is multiplied by 100.

                ind = 0
                nu = np.zeros(0)
                fcind = []
                for comp in components:
                    if 'Power' in comp:
                        ptmp = re.split('\_', comp)
                        #print(comp)
                        #print(re.sub('\'','',ptmp[1]))
                        nui = float(re.sub('\'','',ptmp[1]))/100.
                        nu = np.append(nu,nui)
                        fcind.append(ind)
                    ind += 1
#                print('passou')
                FCBin=[]
                aux=data[_BinFCVec_]
                for fc in aux:
                    FCBini = (np.sum(popx[fcind][(nu > aux[fc][0]) & (nu <= aux[fc][1])]))
                    FCBin.append(FCBini)

                #Searching for the hot dust components, using the same idea as for the FC
                #i = hdind
                
                if not data[_OnlnyFC_]:
                    ind = 0
                    BB = np.zeros(0)
                    hdind = []
                    for comp in components:
                        if 'BB_' in comp:
                            ptmp = re.split('\_', comp)
                            BBi = float(ptmp[1])
                            BB = np.append(BB,BBi)
                            hdind.append(ind)
                        ind += 1
                    HDBin=[]
                    aux= data[_BinHDVec_]
                    for Temp in aux:
                        HDBini = (np.sum(popx[hdind][(BB > aux[Temp][0]) & (BB <= aux[Temp][1])]))
                    HDBin.append(HDBini)
    except: 
        print ("===>>> Problem with input file: {0} <<<=== :(  ".format(data[_filein_]))
        if data[_IsAGNComp_]:
            FCBin = np.zeros(failFileAGN)
            if not data[_OnlnyFC_]:
               HDBin = np.zeros(failFileHD)
        PopBin = np.zeros(failFile)
        SFR = np.zeros(failFileSFR)
        av = np.nan
        meanAgex = np.nan
        meanAgem = np.nan
        meanZx = np.nan
        meanZm = np.nan
        Mcor_t = np.nan
        Mini_t = np.nan
        fobs_norm = np.nan
        adev = np.nan
        chi2 = np.nan
        SNR  = np.nan
        meanAgex1g = np.nan
        meanAgem1g = np.nan
        meanZx1g = np.nan
        meanZm1g = np.nan
        Sigma_star = np.nan
        v_rot = np.nan
        l_obs = np.zeros(int((1+(data[_Olsyn_fin_] - data[_Olsyn_ini_])/data[_DeltaLamb_])))
        f_obs = np.zeros(int((1+(data[_Olsyn_fin_] - data[_Olsyn_ini_])/data[_DeltaLamb_])))
        f_syn = np.zeros(int((1+(data[_Olsyn_fin_] - data[_Olsyn_ini_])/data[_DeltaLamb_])))
        wei = np.zeros(int((1+(data[_Olsyn_fin_] - data[_Olsyn_ini_])/data[_DeltaLamb_])))
        Nbase = re.split(' ', open(data[_BaseFile_]).readline())[0]
        popx = np.zeros(int(Nbase))
        popm = np.zeros(int(Nbase))

    #Append all the relevant data inside the Pops vector
    Pops=[]

    #Special maps calculated for AGNs
    if data[_IsAGNComp_]:
        for v in FCBin:
         Pops.append(v)
         if not data[_OnlnyFC_]:
            for v in HDBin:
              Pops.append(v)

    for v in PopBin:
     Pops.append(v)
    for s in SFR:
     Pops.append(s)
    Pops.append(av)
    Pops.append(meanAgex)
    Pops.append(meanAgem)
    Pops.append(meanZx)
    Pops.append(meanZm)
    if data[_doYoungMage_]:
        Pops.append(meanAgex1g)
        Pops.append(meanAgem1g)
        Pops.append(meanZx1g)
        Pops.append(meanZm1g)
    Pops.append(Mcor_t)
    Pops.append(Mini_t)
    Pops.append(fobs_norm)
    Pops.append(Sigma_star)
    Pops.append(v_rot)
    Pops.append(adev)
    Pops.append(chi2)
    Pops.append(SNR)

    Nbase=int(re.split(' ',open(data[_BaseFile_]).readline())[0])
    try:
        Nbaseextra=int(re.split(' ',open(data[_extrabase_]).readline())[0])
    except:
        Nbaseextra=0

    if  Nbaseextra > Nbase:
        basetest=Nbaseextra
    else:
        basetest=Nbase

    if (len(popx) < basetest):
        x=np.zeros(basetest - len(popx))
        popx=np.append(popx,x)
        popm=np.append(popm,x)
    #print(np.array(Pops))
    return headerKeys,headerComments,np.array(Pops),l_obs,f_obs,f_syn,wei,popx,popm

###########################################################################################
def makeCubeFits(config, category_name = "data"):

        #Read configuration
        data = config[category_name]

        #Get the list of files
        #Read the list from a file
        if data[_StOutFromList_]:
            Table=open(data[_DataList_])
            inputfiles = Table.readlines()
        #Create the list from a directory and the extension
        elif not data[_StOutFromList_]:
            print(pth.join(data[_outdir_], "{0}_*_*.{1}".format(data[_GalaxyShortName_], data[_extension_].lstrip("."))))

            inputfiles=glob.glob(pth.join(data[_outdir_], "{0}_*_*.{1}".format(data[_GalaxyShortName_], data[_extension_].lstrip("."))))
            print(pth.join(data[_outdir_], "{0}_*_*.{1}".format(data[_GalaxyShortName_], data[_extension_].lstrip("."))))
            for i in range(0, len(inputfiles)):
                inputfiles[i] = "{0} {1}".format(inputfiles[i], data[_gal_D_])
        #Get the first file as a reference
#### MODIFICADO #####
        try:
            print(inputfiles[0])
    #             headerExtra=False
    #        except:
    #             tmp_z,tmp_y,tmp_x = np.shape(pf.getdata(glob.glob(pth.join(data[_original_cube_]))[0]))
    #             inputfiles=asarray([])
    #             y,x=indices((tmp_y,tmp_x))
    #             temp=np.column_stack([ravel(y),ravel(x)])
    #             for pix in temp:
    #               outname=data[_outdir_]+data[_GalaxyShortName_]+'_'+str(pix[0]+1)+'_'+str(pix[1]+1)+'.'+data[_extension_]+' 0.0'
    #               inputfiles=np.append(inputfiles,outname)
    ##             print(len(inputfiles))
    ##             for i in range(0, len(inputfiles)):
    ##                inputfiles[i] = inputfiles[i]' 0.0'
    ##                print(inputfiles[i])
    #             print(inputfiles[0])
    #             headerExtra=True
    #### MODIFICADO #####
            topars = re.split("\s+", re.sub("(^\s+)|(\s+(\n)?$)", "", inputfiles[0]))

            #MAYBE IT WONT BE USED (CHANGE THE GETCUBESPARS TO USE TOPARS AS FILEIN AND GAL_D
            config.set(category_name, _filein_, topars[0])
            config.set(category_name, _gal_D_, topars[1])

            #Read some of the common starlight parameters using any of the input files
            headerKeysPars,headerCommentsPars,pars,l,fob,fsyn,wei,popx,popm=GetCubesPars(config, category_name)

            #Adding the used SSPs to the headerKeysPars
            SSPsUsedHeader=np.loadtxt(data[_BaseFile_], usecols=(0,1,2), dtype='str', skiprows=1)
    #### MODIFICADO #####
    #        if headerExtra:
    #            headerKeysPars.append('----- Synthesis was not possible  -----')        
    #### MODIFICADO #####
            headerKeysPars.append('----- SSPs used in the Synthesis, ages, metallicities  -----')
            headerCommentsPars.append(' . ')
            baseAges = []
            baseZs = []
            for k in SSPsUsedHeader:
                headerKeysPars.append(k[0]+',  '+k[1]+',  '+k[2])
                headerCommentsPars.append(' SSPs, age, Z')
                baseAges.append(float(k[1]))
                baseZs.append(float(k[2]))
            #Get the size of the original data
            (z,y,x) = np.shape(pf.getdata(data[_original_cube_]))

            #Creating empty data
            cuboPopBins=np.zeros((len(pars),y,x))
            cuboPopX=np.zeros((len(popx),y,x))
            cuboPopM=np.zeros((len(popm),y,x))
            cuboFobs=np.zeros((len(fob),y,x))
            cuboFsyn=np.zeros((len(fsyn),y,x))
            cuboWei=np.zeros((len(wei),y,x))
            
            #creating a MEF for AGES and Metallicities
            AgeMetal=np.zeros((len(baseZs),2))
            AgeMetal[:,0] = baseAges
            AgeMetal[:,1] = baseZs

            #List of elements "name distance"
            for filein_i in inputfiles:
               # print(filein_i)
                #Get the name and the galaxy distance by removing any spaces and end line characters
                (filein,gal_D) = re.split("\s+", re.sub("(^\s+)|(\s+(\n)?$)", "", filein_i))

                config.set(category_name, _filein_, filein)
                config.set(category_name, _gal_D_, gal_D)

                #Get the name of the file without the extension
                outname = re.sub("[\.][^\.]+$", "", filein)

                #Get the values of y and x from the file name
                y,x=re.split('_',outname)[-2:]

                #Makes a correction to the x and y values (first pixel of the extraction is considered x=1 and y=1)
                #y = int(y) - 1
                #x = int(x) - 1
                y = int(y)
                x = int(x)

                #Run the function to get the values of the synthesis
                CubePars=GetCubesPars(config, category_name)

                #Put the values inside the data to be saved # updated
                cuboPopBins[:, y, x]=CubePars[2]
                cuboFobs[:, y, x]=CubePars[4]
                cuboFsyn[:, y, x]=CubePars[5]
                cuboWei[:, y, x]=CubePars[6]
                cuboPopX[:, y, x]=CubePars[7]
                cuboPopM[:, y, x]=CubePars[8]
            headerKeysPars_r= np.column_stack((headerKeysPars,headerCommentsPars))
            return headerKeysPars_r,cuboPopBins, cuboFobs,cuboFsyn,cuboWei,cuboPopX,cuboPopM,AgeMetal
        except:
           Fail=True
           return Fail

###########################################################################################
#   output_name = "cuboMega.fits"
#   clobber = Flase

def CuboStarlightFits(config, category_name = "data"):

    #Read configuration file
    data = config[category_name]

    #See if the file already exists
    if not pth.exists(data[_save_file_path_]):

        #Run function to create the data
        cuboFinal=makeCubeFits(config, category_name)
        # see if the synthesis has failed or not.
        if cuboFinal == True:
           logFile = open('LogFile.log', 'a')
           logFile.write('Cube has Failed '+glob.glob(pth.join(data[_original_cube_]))[0]+'\n')
           logFile.close()
           print('Cube has Failed '+glob.glob(pth.join(data[_original_cube_]))[0]+'\n'+' LogFile.log updated \n')
        else:
            #Load original fits file
            saveCube=pf.open(data[_original_cube_])

            n_obj = ns.noise(saveCube)
            data_noise_masks = []
            for i in data[_signal_noise_list_]: 
                noise_mask = n_obj.S_TO_N_CUBE(do_mask = True, window = data[_window_range_], SN_tresh = i, mode = "r")
                noise_mask = (noise_mask["S2N_MASK"].data[0]/100. == 0.0)*1.0
                data_noise_masks.append(noise_mask)


            #Adding all the HDUs to the new cube
            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Ages & Mettalicities in the Base'))
            hdr.append(('EXTNAME', 'BaseAgeMetal'))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[7], header=hdr)
            saveCube.append(hdu)

            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Synthesis Parameters & Bined Population Vectors'))
            hdr.append(('EXTNAME', 'PoPBins'))
            for h in np.arange(0,(len(cuboFinal[0]))):
                   hdr.append(('DATA'+str(h), cuboFinal[0][h][0],cuboFinal[0][h][1]))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[1], header=hdr)
            saveCube.append(hdu)
    #        print(hdr)

            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Population Vectors Not Bined in Light Fractions'))
            hdr.append(('EXTNAME', 'PoPVecsL'))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[5], header=hdr)
            saveCube.append(hdu)

            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Population Vectors Not Bined in Mass Fractions'))
            hdr.append(('EXTNAME', 'PoPVecsM'))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[6], header=hdr)
            saveCube.append(hdu)

            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Observed Flux (STARLIGHT output) in input units (erg/cm-2/s/A)'))
            hdr.append(('EXTNAME', 'OBSFLUX'))
            hdr.append(('CRPIX3',1))
            hdr.append(('CRVAL3',data[_Olsyn_ini_]))
            hdr.append(('CD3_3',data[_DeltaLamb_]))
            hdr.append(('WAT3_001', 'wtype=linear axtype=wave' ))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[2], header=hdr)
            saveCube.append(hdu)

            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Synthetic Flux (STARLIGHT output) in input units  (erg/cm-2/s/A)'))
            hdr.append(('EXTNAME', 'FLXSYN'))
            hdr.append(('CRPIX3',1))
            hdr.append(('CRVAL3',data[_Olsyn_ini_]))
            hdr.append(('CD3_3',data[_DeltaLamb_]))
            hdr.append(('WAT3_001', 'wtype=linear axtype=wave' ))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[3], header=hdr)
            saveCube.append(hdu)

            hdr = pf.Header()
            hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
            hdr.append(('OBJECT','Weight used in the fits (STARLIGHT output)'))
            hdr.append(('EXTNAME', 'WEIGHT'))
            hdr.append(('CRPIX3',1))
            hdr.append(('CRVAL3',data[_Olsyn_ini_]))
            hdr.append(('CD3_3',data[_DeltaLamb_]))
            hdr.append(('WAT3_001', 'wtype=linear axtype=wave' ))
            hdu = pf.hdu.ImageHDU(data=cuboFinal[4], header=hdr)
            saveCube.append(hdu)

            s = 0
            for i in data[_signal_noise_list_]:
                hdr = pf.Header()
                hdr.append(('AUTHORS', 'Rogerio Riffel (riffel@ufrgs.br) & Nicolas Dullius Mallmann (nicolas.mallmann@ufrgs.br)'))
                hdr.append(('OBJECT', 'Signal/Noise mask maps'))
                hdr.append(('EXTNAME', 'SN_MASKS_{0}'.format(i)))
                hdu = pf.hdu.ImageHDU(data=data_noise_masks[s], header=hdr)
                saveCube.append(hdu)
                s = s + 1

            saveCube.writeto(data[_save_file_path_])
