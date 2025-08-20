
import numpy as np
import sys
import h5py as h5

#Quick fudge to make import from ../Scripts work
sys.path.append('../../common_code/')


from PostProcessingScripts import * 
# from formation_channels import * 

import pandas as pd
from astropy import units as u
from astropy import constants as const

channelList = ['classic', 'stable B no CEE', 'vii',  'immediate CE',  r'double-core CE', 'other', 'vi']




def add_formation_channel_H5file(pathToData, DCOtype, BPSmodelName, redshifts):

    DCOname = DCOname_dict[DCOtype]

    # path for files 
    path  = pathToData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'

    # read in data 
    fdata = h5.File(path, 'r+')

    # delete the group if it already exists (so that we overwrite it)
    # if fdata.get("formationchannel_z_rates"):
    #     del f['dataset_name']
    if "formationchannel_z_rates" in fdata.keys():
        print("formation channel rates group already exists... deleting it now")
        del fdata["formationchannel_z_rates"]

    Rfc = fdata.create_group("formationchannel_z_rates")

    fDCO = fdata['doubleCompactObjects']
    channels = fdata['doubleCompactObjects']["formaton channel"][...].squeeze()
    # channels = fdata['doubleCompactObjects']['formaton channel'][...].squeeze()

    # add a line with the redshifts 
    channel_no = fdata["formationchannel_z_rates"].create_dataset(u"redshift", data=redshifts)
    channel_no.attrs[u"units"] = u"#"    



    # add all 
    # for ind_mssfr, mssfr in enumerate(['312']): # temp fix
    for ind_mssfr, mssfr in enumerate(MSSFRnameslist[:]):
        print('adding all for mssfr ', mssfr)
        # print(fdata.keys())
        fparam_key = 'weights_intrinsicPerRedshift'
        fractions_z = np.zeros_like(redshifts)
        for ind_z, redshift in enumerate(redshifts):
            redshift = np.round(redshift,4)
            weightheader = 'w_' + mssfr + '_z_' +  str(redshift)
            weights_ = fdata[fparam_key][weightheader][...].squeeze()  
            if np.sum(weights_)!=0:
                fractions_z[ind_z] = np.sum(weights_)
            del weights_
            # else:
            #     fractions_z = np.zeros_like(redshifts)

        # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
        # Rdetxyz = Rdet.create_dataset(u"w_%s"%mssfrName, data=weightsDet)
        # Rdetxyz.attrs[u"units"] = u"yr^{-1}"
        # append total intrinsic rate at redshift 0 
        fc_header_name = "rate_" + "total" + "_" + mssfr 
        channel_no = fdata["formationchannel_z_rates"].create_dataset(u"%s"%fc_header_name, data=fractions_z)
        channel_no.attrs[u"units"] = u"Gpc^-3 yr^-1"  
        del fractions_z
    
    for nrC, Channel in enumerate(channelList): 
        fc_ind_wanted = dictFormationChannelIndex[Channel]
        mask_fc = (channels==fc_ind_wanted)
        # del channels
        for ind_mssfr, mssfr in enumerate(MSSFRnameslist[:]):
            print('mssfr ', mssfr)
            if np.sum(mask_fc)>1:
                fparam_key = 'weights_intrinsicPerRedshift'
                fractions_z = np.zeros_like(redshifts)
                for ind_z, redshift in enumerate(redshifts):
                    redshift = np.round(redshift,4)
                    weightheader = 'w_' + mssfr + '_z_' +  str(redshift)
                    weights_ = fdata[fparam_key][weightheader][...].squeeze()  
                    if np.sum(weights_)!=0:
                        fractions_z[ind_z] = np.sum(weights_[mask_fc])/np.sum(weights_)
                    del weights_
                        
            else:
                fractions_z = np.zeros_like(redshifts)

            # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
            # Rdetxyz = Rdet.create_dataset(u"w_%s"%mssfrName, data=weightsDet)
            # Rdetxyz.attrs[u"units"] = u"yr^{-1}"
            # append total intrinsic rate at redshift 0 
            fc_header_name = "fraction_" + Channel + "_" + mssfr 
            channel_no = fdata["formationchannel_z_rates"].create_dataset(u"%s"%fc_header_name, data=fractions_z)
            channel_no.attrs[u"units"] = u"#"    




    fdata.close() 

    print()
    print('-----------------------------------------------')
    print('completed')
    print('added formation channels rates to their own file :) YES!')


#, 'BBH' 'G', #'A', 'B', 'C', 'D', 'E', 'F',  'G',  'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'



def obtain_redshiftsruns(pathData = '/Volumes/SimonsFoundation/DataDCO/'):
    BPSmodelName='A'
    DCOtype='BNS'
    path_ = '/Volumes/SimonsFoundation/DataDCO/' + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
    fdata = h5.File(path, 'r')
    redshifts = fdata['redshifts']['redshift'][...].squeeze()
    fdata.close()
    return redshifts 

##### PLOT LVKM1 

pathData='/Volumes/SimonsFoundation/DataDCO/'
redshifts_runs = obtain_redshiftsruns(pathData = pathData)

#  
for DCOtype in ['BBH']:
    for BPSmodelName in ['K']:

        pathToData = '/Volumes/SimonsFoundation/DataDCO/'


        print()
        print('-------------------------------------------------------')
        print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
        print('-------------------------------------------------------')
        # print('which is in this directory', pathToDataWithoutCOMPASname)
        # print('optimistic = ', optimistic)

        add_formation_channel_H5file(pathToData=pathToData, DCOtype=DCOtype, BPSmodelName=BPSmodelName, redshifts=redshifts_runs)




