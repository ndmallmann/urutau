import numpy as np
import os, glob

def starlightPars(filename):
        '''
        This routine retuns some important parameters from starlight V4 output.
        usage: parms=starlightPars(filename)
           parms[0] = array with parameters
           parms[1] = array with the list of parameters extracted.
           parms[2] = Beginig of the population vectos fractions
           parms[3] = End of the population vectos fractions.
           parms[4] = Beginig of the syntheis results
           The last 4 parameters are used by popVectors and retur the Starlight version
        '''

# To add more parameters only increase the ParmList array
       
        cnt=1
        pars=np.zeros(0)
        Table=open(filename)
        filecont = Table.readlines()
        if 'PANcMExStarlight' in filecont[1]:
                version='V5'
                ParList=['[chi2/Nl_eff', 'chi2]','flux_unit]','fobs_norm','Lobs_norm','LumDistInMpc','[adev (%)]','[Lum_tot (Lsun/A if distance & flux_unit are Ok...)]',
                        '[Mini_tot (Msun  if distance & flux_unit are Ok...)]','[Mcor_tot (Msun  if distance & flux_unit are Ok...)]','[AV_min  (mag)]',
                        '[sum-of-x (%)]','[l_norm (A) - for base]','[llow_norm (A) - window for f_obs]','[lupp_norm (A) - window for f_obs]',
                        '[v0_min  (km/s)]','[S/N in S/N window]','[N_base]','[sum-of-x (%)]','[vd_min  (km/s)]']
                for par in ParList:
                    for line in filecont:
                        if par in line:
                           tline = line.split()
                           if par=='chi2]':
                               parm=float(tline[0:][1])
                               pars=np.append(pars,parm)
                           elif par=='flux_unit]':
                               parm=float(tline[0:][1])
                               pars=np.append(pars,parm)
                           elif par=='Lobs_norm':
                               parm=float(tline[0:][1])
                               pars=np.append(pars,parm)

                           elif par=='LumDistInMpc':
                               parm=float(tline[0:][2])
                               pars=np.append(pars,parm)
                           else:
                               parm=float(tline[0:][0])
                               pars=np.append(pars,parm)

                for line in filecont:
                      if('x_j(%)' in line):
                           popin=cnt
                      if('## Synthesis Results' in line):
                           popend=cnt-2
                      if('## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei' in line):
                           syntin=cnt+1
                      cnt=cnt+1

       
        if 'StarlightChains_v04' in filecont[1]:
            version='V4'
            ParList=['[chi2/Nl_eff]','[fobs_norm (in input units)]','[adev (%)]','[Flux_tot (units of input spectrum!)]',
                    '[Mini_tot (???)]','[Mcor_tot (???)]','[AV_min  (mag)]','[YAV_min (mag)]','[sum-of-x (%)]',
                    '[l_norm (A) - for base]','[llow_norm (A) - window for f_obs]','[lupp_norm (A) - window for f_obs]',
                    '[v0_min  (km/s)]','[S/N in S/N window]','[N_base]','[sum-of-x (%)]','[vd_min  (km/s)]']       
            for par in ParList:
                for line in filecont:
                    if par in line:
                       tline = line.split()
                       parm=float(tline[0:][0])
                       pars=np.append(pars,parm)
            for line in filecont:
                  if('x_j(%)' in line):
                       popin=cnt
                  if('## Synthesis Results' in line):
                       popend=cnt-2
                  if('## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei' in line):
                       syntin=cnt+1
                  cnt=cnt+1
              
              
        return pars,ParList,popin,popend,syntin,version

def popVectors(filein):
        '''
        This routine is used to read the population vectors table from Starlight Output. It loads 
        x_j(%)      Mini_j(%)     Mcor_j(%)     age_j(yr)     Z_j      (L/M)_j     and j, Components_j to an array. 
        Usage:
        x=popVectors(filein)
        Returns:
        x[0] - all float components
        x[1] - str component_j
        '''
        x=starlightPars(filein)
        linebeg=int(x[2])
        linefin=int(x[3])
        tmp=open(filein)
        tmpcont=tmp.readlines()
        tempcut=tmpcont[linebeg:linefin]
        save=open(filein+'.temp','w+')
        for line in tempcut:
             save.write(line)
        save.close()
        pop=np.loadtxt(filein+'.temp',usecols=(1,2,3,4,5,6,0))
        popComps=np.loadtxt(filein+'.temp',usecols=(9,),dtype='str')
        os.system("rm -rf  %s"%(filein+'.temp'))
        return pop,popComps

#def StSyntesis(filein):
#        '''
#        This routine loads the Synthetic spectrum (Best Model): l_obs f_obs f_syn wei
#        Usage
#        x=StSyntesis(filein)
#        '''
#        x=starlightPars(filein)
#        skip=int(x[4])
#        # Fror compatibility with new numpy version
#        try:
#                spec=np.genfromtxt(filein,skiprows=skip,missing_values='*********')
#        except:
#                spec=np.genfromtxt(filein,skip_header=skip,missing_values='*********')
#        return spec

def StSyntesis(filein):
        '''
        This routine loads the Synthetic spectrum (Best Model): l_obs f_obs f_syn wei
        Usage
        x=StSyntesis(filein)
        '''
        x=starlightPars(filein)
        skip=int(x[4])
        # Fror compatibility with new numpy version
        try:
                spec=np.genfromtxt(filein,skiprows=skip,missing_values='*',filling_values='-1.00')
        except:
                spec=np.genfromtxt(filein,skip_header=skip,missing_values='*',filling_values='-1.00')
        return spec






