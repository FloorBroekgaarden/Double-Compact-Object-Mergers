
import numpy as np
import sys
import h5py as h5

#Quick fudge to make import from ../Scripts work
sys.path.append('../../common_code/')


from PostProcessingScripts import * 
from formation_channels import * 

import pandas as pd
from astropy import units as u
from astropy import constants as const






def add_formation_channel_H5file(pathToData, DCOtype, BPSmodelName):

    DCOname = DCOname_dict[DCOtype]

    # path for files 
    path_dir = pathToData
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'

    # read in data 
    fdata = h5.File(path, 'r+')


    # set optimistic true if that is the variation (H) 
    # OPTIMISTIC=False
    # if (bps_model=='F') or (bps_model=='K'):
    #     OPTIMISTIC=True
    #     print('doing optimistic version of %s'%alphabetDirDict[bps_model])

    # immediateRLOFAfterCEE = fdata["commonEnvelopes"]["immediateRLOFAfterCEE"][...].squeeze()
    # immediateRLOFAfterCEEmask = (immediateRLOFAfterCEE==1)
    # CErandomSeed = fdata["commonEnvelopes"]["randomSeed"][...].squeeze()
    
    # seeds_withImmediateRLOFAfterCEE = np.in1d(fdata['doubleCompactObjects']['seed'][...].squeeze(), CErandomSeed[immediateRLOFAfterCEEmask])
    # print('systems with withImmediateRLOFAfterCEE',np.sum(seeds_withImmediateRLOFAfterCEE), ' which is a fraction', np.sum(seeds_withImmediateRLOFAfterCEE)/len(seeds_withImmediateRLOFAfterCEE))



    # delete the group if it already exists (so that we overwrite it)
    if "formaton channel" in fdata['doubleCompactObjects'].keys():
        print("formation channel rates group already exists... deleting it now")
        del fdata['doubleCompactObjects']["formaton channel"]

    seeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
    channels = identify_formation_channels(seeds=seeds, file=fdata)
    print(fdata.keys())
    # fdata.close()

    #######
    # fdata = h5.File(path 'r+')

    # Rdet = fdata.create_group("weights_detected")
    fDCO = fdata['doubleCompactObjects']


    # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
    # Rdetxyz = Rdet.create_dataset(u"w_%s"%mssfrName, data=weightsDet)
    # Rdetxyz.attrs[u"units"] = u"yr^{-1}"
    # append total intrinsic rate at redshift 0 
    channel_no = fdata['doubleCompactObjects'].create_dataset(u"formaton channel", data=channels)
    channel_no.attrs[u"units"] = u"#"    




    fdata.close() 

    print()
    print('-----------------------------------------------')
    print('completed')
    print('added formation channels numbers to DCO file :) YES!')


#, 'BBH' 'G', #'A', 'B', 'C', 'D', 'E', 'F',  'G',  'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'

for DCOtype in ['BBH']:
    for BPSmodelName in ['K']:

        pathToData = '/Volumes/SimonsFoundation/DataDCO/'


        print()
        print('-------------------------------------------------------')
        print(' now at ', DCOtype, 'adding FC CHANNELS for modelname', BPSmodelName)
        print('-------------------------------------------------------')
        # print('which is in this directory', pathToDataWithoutCOMPASname)
        # print('optimistic = ', optimistic)

        add_formation_channel_H5file(pathToData=pathToData, DCOtype=DCOtype, BPSmodelName=BPSmodelName)




