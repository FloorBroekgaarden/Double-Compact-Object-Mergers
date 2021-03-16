import cv2
import os
import string
import os
import moviepy.video.io.ImageSequenceClip
from PIL import Image
import glob





def makeMovie_rates(whichRate='intrinsic', fps=.4, duration=300):
	'''
	whichRate = 'intrinsic' or 'observed'
	fps=0.4, frames per second
	duration = duration of the movie 
	'''




	GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
	MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
	SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']


	MSSFRnameslist = []
	MSSFRnameslist.append('000') # add phenomenological 



	for ind_SFR, SFR in enumerate(SFRs):
		ind_x = ind_SFR + 1
		for ind_GSMF, GSMF in enumerate(GSMFs):
			ind_y = ind_GSMF + 1
			for ind_MZ, MZ in enumerate(MZs):
				ind_z = ind_MZ + 1
				
				MSSFRnameslist.append('%s%s%s'%(ind_x, ind_y, ind_z))



                
  
	image_folder = '../plottingCode/Fig_2/supplementary_material/'

	images = []
	if whichRate=='intrinsic':
		for ind_m, SFRD_model in enumerate(MSSFRnameslist):
			images.append(image_folder +   'Rates_intrinsic_'  + SFRD_model + '.png')
	elif whichRate=='observed':
		for ind_m, SFRD_model in enumerate(MSSFRnameslist):
			images.append(image_folder +   'Rates_observed_'  + SFRD_model + '.png')


	image_files = images
	clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
	clip.write_videofile(image_folder+'movie_'+ 'rates_' + whichRate + '.mp4')


	# make also gif:
 
	# Create the frames
	frames = []
	# imgs = glob.glob("*.png")
	for i in images:
	    new_frame = Image.open(i)
	    frames.append(new_frame)
	 
	# Save into a GIF file that loops forever
	frames[0].save(image_folder+'gif_'+ 'rates_' + whichRate +  '.gif', format='GIF',
	               append_images=frames[1:],
	               save_all=True,
	               duration=duration, loop=0)



	return 



makeMovie_intrinsicRates=True



# Run rhis using python 3!! 

if makeMovie_intrinsicRates==True:
	makeMovie_rates(whichRate='intrinsic')

