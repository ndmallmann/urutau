#Programmer: Jaderson Schimoia
#Modified by: Nicolas Dullius Mallmann
#Name: noisy
#Status: Developing

import numpy
import pylab
import astropy.io.fits as pyfits

"""
DESCRIPTION: THIS SCRIPT CALCULATES THE S/N RATIO OF THE MaNGA DATA CUBES
AND ADD THE RESULT AS A NEW EXTENSION TO THE FILE
"""


class noise(object):

    def __init__(self, manga_cube):

        #MODIFIED
        if not isinstance(manga_cube, pyfits.HDUList):
            hdulist = pyfits.open(manga_cube)
        else:
            hdulist = manga_cube
        #hdulist is a pyfits object
        self.header=hdulist[0].header



        self.flux=hdulist['FLUX'].data
        #gets the data flux of the 'FLUX' extension in the datacube
        self.flux_header=hdulist['FLUX'].header
        #gets the header of the data flux extension from the data cube


        #gets the NAXIS sizes from the flux header
        self.naxis1=self.flux_header['NAXIS1']
        self.naxis2=self.flux_header['NAXIS2']
        self.naxis3=self.flux_header['NAXIS3']

        self.ivar=hdulist['IVAR'].data
        #gets the ivar data from the 'IVAR' extension in the datacube
        self.ivar_header=hdulist['IVAR'].header
        #gets the header of the 'IVAR' extension in the data cube
        self.mask=hdulist['MASK'].data
        #gets the mask data from the datacube
        self.mask_header=hdulist['MASK'].header
        #gets the mask header data from the datacube


        #save the name (or instance) of the input file

        self.name=manga_cube


        #This part of the function writes the wavelength
        #The dispersion is linear in this case

        pix_ref=self.flux_header['CRPIX3']
        #reference pixel of the wavelength
        wave_ref=self.flux_header['CRVAL3']
        #starting wavelength at the reference pixel
        wave_delt=self.flux_header['CD3_3']
        #delta wavelength in linear dispersion
        wave_size=self.flux_header['NAXIS3']
        #size of the spectral dimension

        #print wave_size debugging
        index=numpy.arange(wave_size)+1
        self.wave=numpy.zeros(wave_size) + wave_ref + (index-pix_ref)*wave_delt
        #self.wave is the wavelength array of the cube

        #MODIFIED
        #hdulist.close()




    def splot(self,x,y):
        """
        plot the spectrum for a given (x,y) position
        """
        #plot the spectrum for a given (x,y) position

        #check for mask values
        mask_index=numpy.where(self.mask[:,x,y]>0)[0]


        pylab.plot(self.wave,self.flux[:,x,y])
        pylab.fill_between(self.wave,numpy.zeros(numpy.size(self.wave)),self.mask[:,x,y],color='red',alpha=0.3)
        pylab.ylim(0,1.15*numpy.max(self.flux[:,x,y]))


    def s_to_n(self,x,y,window=None, plot=False):
        """
        Calculates the signal-to-noise spectrum at a given (x,y) positon, with the ollowing steps
        1) Check for good and bad pixels using the MASK extension
        2) Badd pixels are flaged with S/N = 0
        3) Foo good pixels valvulate the flux uncertainty using the inverse of the variance extension:
            flux_er=1./sqrt(ivar)
        4) Divide the flux by the flux_er pixel-to-pixel to generate the signal-to-noise spectrum
        5) Return the signal-to-noise spectrum at the (x,y) position

        """
        #plot the signal-to-noise ration spectrum for a given (x,y) position




        #first check using mask if the x,y pixel contains only bad pixels:

        mask_good=numpy.where(self.mask[:,x,y]==0)[0]
        #checking good pixels


        #print mask_good debugging

        mask_bad=numpy.where(self.mask[:,x,y]>0)[0]
        #checking bad pixels


        #print mask_bad debuggin

        #create an empty array to containg the S/N spectrum
        S_to_N=numpy.zeros(numpy.size(self.wave))


        #print self.naxis3  #debugging
        if numpy.size(mask_good)>0:
            #If we have good pixels, then calculate statistics


            #S_to_N at the good pixels:

            S_to_N[mask_good]=(self.flux[mask_good,x,y])/(1./numpy.sqrt(self.ivar[mask_good,x,y]))



            #S_to_N at the bad pixels:

            if numpy.size(S_to_N[mask_bad])!=0:
                S_to_N[mask_bad]=0

            #print("(X,Y)=("+str(x)+","+str(y)+")   MEAN S/N RATIO:"+str(int(numpy.mean(S_to_N))))

            if plot==True:
                pylab.plot(self.wave,S_to_N)
                pylab.fill_between(self.wave,numpy.zeros(numpy.size(self.wave)),self.mask[:,x,y],color='red',alpha=0.3)
                pylab.ylim(0,1.15*numpy.max(S_to_N))

        else:
            #if we don't have good pixels, thens S/N == 0
            #print("WARNING: ONLY BAD SPECTRAL PIXELS AT THIS X,Y POSITION")
            pass


        return S_to_N


    def S_TO_N_CUBE(self,do_mask=False,window=None,SN_tresh=10, mode='r'):
        """
        Calcultes the signal-to-noise cube using the s_to_n definition.
        """

        l_index=0
        u_index=numpy.size(self.wave)

        print(do_mask,SN_tresh,'SNR')
        if window:
            #CORRECTION MADE: 
            #BEFORE: l_index=numpy.where(self.wave<=window[0])[0][0]
            l_index=numpy.where(self.wave<=window[0])[0][-1]
            u_index=numpy.where(self.wave>=window[1])[0][0]


        #=========================================
        #This part will be executed if a mask for good signal-to-noise pixels is desired

        if do_mask==True:
            tmp_sn_mask=numpy.ndarray([1,self.naxis1,self.naxis2])


        #first create an empty cube from the naxis sizes:
        tmp_cube=numpy.ndarray([self.naxis3,self.naxis1,self.naxis2])

        #double loop over x and y to calculate the S/N spectra
        for x in range(self.naxis1):
            for y in range(self.naxis2):
                tmp_cube[:,x,y]=self.s_to_n(x,y)


        # calculate the statistics for the specific wavelength range:
                if do_mask==True:
                    #print l_index, u_index
                    #print tmp_cube[:,x,y]
                    tmp_sn_mask[0,x,y]=numpy.mean(tmp_cube[:,x,y][l_index:u_index])

        #for the signal-to-noise cube I'll keep the header of the FLUX extension
        S_TO_N_header=self.flux_header.copy()


        #first update a few keywords
        S_TO_N_header['EXTNAME']=('S2N')


        #MODIFIED
        if isinstance(self.name, pyfits.HDUList):
            hdulist=self.name
        else:
            hdulist=pyfits.open(self.name)

        #extract all extension's names keywords
        extname=numpy.array([])
        for i in hdulist[1::]:
            extname=numpy.append(extname,i.header['EXTNAME'])


        if do_mask==True:
            I=numpy.where(tmp_sn_mask>=SN_tresh)
            MASK=tmp_sn_mask.copy()
            MASK=MASK*0
            MASK[I]=100


            #create a new header extension for the mean SN_IMAGE
            #I will use the IVAR header as a base for the SN_IMAGE extension
            #sn_image_header=self.ivar_header.copy()
            sn_image_header=self.flux_header.copy()
            #sn_image_header['EXTNAME']=('S2N_IMAGE','Signal-to-noise image')
            sn_image_header['EXTNAME']=('S2N_IMAGE')

            #if window:
            #    key_window=str(window)
            #else:
            #    key_window='all_spectrum'
            #sn_image_header.update('WINDOW',key_window,'wavelength window range to calculate the mean signal-to-noise')
            #sn_image_header.update('HDUCLAS3','S2N_RATIO')




            #sn_mask_header=self.ivar_header.copy()
            sn_mask_header=self.flux_header.copy()
            #sn_mask_header['EXTNAME']=('S2N_MASK','Signal-to-noise mask')
            sn_mask_header['EXTNAME']=('S2N_MASK')



            #There are three options for output in this definition:
            #r: return
            #a: append
            #u: update


            #check whixh option is desired by the user:
            #print('*** Output of the SIGNAL-TO-NOISE extension ***')
            while ((mode!='r')&(mode!='a')&(mode!='u')):
                #print("WARNING: Output option is not a valid option!")
                mode=raw_input('Please enter a valid option:\nr: return\n\: append\nu: update\n: ')

            if mode=='a':
                #check for existence of the extension:
                if 'S2N' in extname:
                    #print('WARNING: signal-to-noise already exists!')
                    mode=raw_input('a: append anyway!\nu: update existing signal-to-noise extension\n: ')


            if mode=='u':
                #check for existence of the extension:
                if not 'S2N' in extname:
                    #print('WARNING: there is no signal-to-noise extension to update!')
                    #print('Setting mode=a for appending')
                    mode='a'


            #Now it's time for output the results of the calculations:
            if mode=='a':
                #print('Appending SIGNAL-TO-NOISE extensions:')
                #appending
                #S2N = tmp_cube
                #S2N_IMAGE = tmp_sn_mask
                #S2N_MASK = MASK
                pyfits.append(self.name,tmp_cube,header=S_TO_N_header)
                #print('S2N append-> Done!')
                pyfits.append(self.name,tmp_sn_mask,header=sn_image_header)
                #print('S2N_IMAGE append-> Done!')
                pyfits.append(self.name,MASK,sn_mask_header)
                #print('S2N_MASK append-> Done!')


            if mode=='u':
                #print('Updating signal-to-noise extensions')
                pyfits.update(self.name,data=tmp_cube,ext=hdulist.index_of('S2N'),header=S_TO_N_header)
                #print('S2N update-> Done!')
                pyfits.update(self.name,data=tmp_sn_mask,ext=hdulist.index_of('S2N_IMAGE'),header=sn_image_header)
                #print('S2N_IMAGE update-> Done!')
                pyfits.update(self.name,data=MASK,ext=hdulist.index_of('S2N_MASK'),header=sn_mask_header)
                #print('S2N_MASK update-> Done!')


            if mode=='r':
                #***
                #***
                #VERIFY WHY THERE IS NO NAME FOR EXTENSIONS HERE!!!!
                #print('Returning signal-to-noise CUBE!')
                #create the cube that will be returned
                new_cube=pyfits.HDUList()
                new_cube.append(pyfits.PrimaryHDU()) #Primary HDU
                new_cube.append(pyfits.ImageHDU())  # S2N extension
                new_cube.append(pyfits.ImageHDU())  #S2N_IMAGE extension
                new_cube.append(pyfits.ImageHDU())  #S2N_MASK extension

                new_cube[0].header=hdulist[0].header
                new_cube[1].data,new_cube[1].header=tmp_cube,S_TO_N_header
                new_cube[2].data,new_cube[2].header=tmp_sn_mask,sn_image_header
                new_cube[3].data,new_cube[3].header=MASK,sn_mask_header

                new_cube[1].header["EXTNAME"] = 'S2N'
                new_cube[2].header["EXTNAME"] = 'S2N_IMAGE'
                new_cube[3].header["EXTNAME"] = 'S2N_MASK'
                return new_cube
                #will return the entyre cube in which the prymary extension is the same of the



            #pyfits.update(self.name,tmp_sn_mask,ext='S2N_IMAGE')




            #pyfits.open('signaltonoise_cube.fits',tmp_sn_mask,header=sn_image_header)

            #falta anexar o treshhold for the signal-to-noise
            #pyfits.append('signaltonoise_cube.fits',MASK,header=sn_mask_header)

            #pyfits.update(self.name,MASK,ext='S2N_MASK')

            #MODIFIED
            #hdulist.close()

        #S_TO_N_CUBE.writeto('signaltonoise_cube.fits')

        #return S_TO_N_CUBE


