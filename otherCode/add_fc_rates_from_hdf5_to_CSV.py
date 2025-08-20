
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





def initalize_formationChannels(DCOname):

    headerDict_intrinsic = { 5:'All intrinsic (z=0) [Gpc^-3 yr^-1]',  0:'channel V intrinsic (z=0)',  1:'channel I intrinsic (z=0)', 2:'channel II intrinsic (z=0)',3:'channel III intrinsic (z=0)', 4:'channel IV intrinsic (z=0)'}
    stringgg =  'Formation_Channels_Local_Rates'
    headerDict_intrinsicOther = {'other_withCE':'channel V with CE intrinsic (z=0)', 'other_withoutCE':'channel V without CE intrinsic (z=0)'}


    for MSSFR in MSSFRnameslist:
        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '_xyz_' + MSSFR + '.csv' 

        NAMES = []
        for Channel in range(len(headerDict_intrinsic)):
            namez0 = headerDict_intrinsic[Channel]      
            NAMES.append(namez0)

        # add special 'other channels'
        for Channel in ['other_withCE', 'other_withoutCE']:
            namez0 = headerDict_intrinsicOther[Channel]      
            NAMES.append(namez0)            


        datas=[]
        for ii in range(len(headerDict_intrinsic)):
            datas.append(np.zeros_like(BPSnameslist))    
        # add special 'other channels'  
        for ii in range(len(['other_withCE', 'other_withoutCE'])):
            datas.append(np.zeros_like(BPSnameslist)) 
        
        df = pd.DataFrame(data=datas, index=NAMES, columns=BPSnameslist).T
        df.columns =   df.columns.map(str)
        df.index.names = ['formation channel']
        df.columns.names = ['model']
        df.to_csv(writePath)


    return 



def writeToRatesFile_FormationChannels(pathData='/Volumes/SimonsFoundation/DataDCO/', DCOtype='BHNS'):
    """ DCOType = ['BBH', 'BHNS', BNS]"""

    headerDict_intrinsic = { 5:'All intrinsic (z=0) [Gpc^-3 yr^-1]',  'other':'channel V intrinsic (z=0)',  'classic':'channel I intrinsic (z=0)', 'stable B no CEE':'channel II intrinsic (z=0)','immediate CE':'channel III intrinsic (z=0)', \
    r'double-core CE':'channel IV intrinsic (z=0)', 'other_withCE':'channel V with CE intrinsic (z=0)', 'other_withoutCE':'channel V without CE intrinsic (z=0)'}
    # headerDict_intrinsicOther = {'other_withCE':'channel V with CE intrinsic (z=0)', 'other_withoutCE':'channel V without CE intrinsic (z=0)'}

    DCOname = DCOname_dict[DCOtype]
    which_z_ind=0 # gives lowest redshift. 
    stringgg =  'Formation_Channels_Local_Rates'
    fparam_key = "formationchannel_z_rates"


    for ind_L, MSSFRname in enumerate(MSSFRnameslist[:]):
        print('at MSSFR: ', MSSFRname)
        #iterate over the formation channels 
        for ind_c, whichChannel in enumerate(['classic', 'stable B no CEE',  'immediate CE',  r'double-core CE', 'other', 'other_withCE', 'other_withoutCE']): #'vii', , 'vi'
            merger_ratio_z = np.zeros(len(BPSnameslist))
            #iterate over stellar evolution models 
            # for ind_m, BPSmodelName in enumerate(["I"]):
            for ind_m, BPSmodelName in enumerate(BPSnameslist):
                full_data_path  = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
                fparam_key = "formationchannel_z_rates"
                header = "fraction_" + whichChannel + "_" + MSSFRname
                fdata = h5.File(full_data_path,'r')
                fc_fraction = fdata[fparam_key][header][...].squeeze()[which_z_ind]  # [which_z_ind = 0] gives at lowest redshift
                # add the case A MT systems that initially were split off 
                if whichChannel=='classic':
                    header = "fraction_" + 'vi' + "_" + MSSFRname
                    fc_fraction += fdata[fparam_key][header][...].squeeze()[which_z_ind]  # [which_z_ind = 0] gives at lowest redshift  
                elif whichChannel=='stable B no CEE':
                    header = "fraction_" + 'vii' + "_" + MSSFRname
                    fc_fraction += fdata[fparam_key][header][...].squeeze()[which_z_ind]  # [which_z_ind = 0] gives at lowest redshift  
                merger_ratio_z[ind_m] = fc_fraction
                fdata.close()
                
            writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '_xyz_' + MSSFRname + '.csv' 
            df = pd.read_csv(writePath, index_col=0)                
            column_name = headerDict_intrinsic[whichChannel] 
            df[column_name] = merger_ratio_z
            df.to_csv(writePath)

        # add total rate: 
        totals_z = np.zeros(len(BPSnameslist))
        for ind_m, BPSmodelName in enumerate(BPSnameslist):
            full_data_path  = pathData + alphabetDirDict[BPSmodelName] +'/COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            header = "rate_" + "total" + "_" + MSSFRname
            fdata = h5.File(full_data_path,'r')
            totals_z[ind_m] =  fdata[fparam_key][header][...].squeeze()[which_z_ind] # [which_z_ind = 0] gives at lowest redshift
            fdata.close()

        writePath = '/Users/floorbroekgaarden/Projects/GitHub/DCO_FormationChannels/dataFiles/data_Fig_1/'  + stringgg + '_'  + DCOname + '_xyz_' + MSSFRname + '.csv' 
        df = pd.read_csv(writePath, index_col=0)                
        column_name = 'All intrinsic (z=0) [Gpc^-3 yr^-1]' 
        df[column_name] = totals_z
        df.to_csv(writePath)


    return
####







INITIALIZE_FormationChannels = True #False #True #True

if INITIALIZE_FormationChannels==True:
    for DCOtype in ['BNS', 'BHNS','BBH']:
    # for DCOtype in ['BHNS','BBH', 'BNS']:
        DCOname = DCOname_dict[DCOtype]
        initalize_formationChannels(DCOname)
    print('done')

runFormationChannels=True#False #True# #False

if runFormationChannels ==True:
    # for DCOtype in ['BHNS', 'BNS', 'BBH']:
    for DCOtype in [ 'BNS', 'BHNS','BBH']:
        print('at DCOtype =', DCOtype)
        writeToRatesFile_FormationChannels(DCOtype=DCOtype)
        print('done with ', DCOtype)





