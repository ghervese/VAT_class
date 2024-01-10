# VAT_class
A python class with GUI to automatize the download of .fits corresponding to a set of tiles for specific deep space object and different surveys available with SkyView
                                            VAT - Virtual Astrophotographer Tool                                                                               V0.1 by Guillaume Hervé-Secourgeon // herve-guillaume[at]orange.fr                               
----------------------------------------------------------------------------------------------------------------

This Python class based on astropy and astroquery is dedicated to generate a set of fits images based on NASA's SkyView open observatory.
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
0 - The prerequisites for this Python class that are to be installed are the following classes:\n
      -astroquery.skyview
      -astropy units
      -astropy.coordinates import SkyCoord
      -astropy.coordinates import Angle
      -astropy.coordinates import position_angle
      -astropy.wcs import WCS
      -astropy.visualization import astropy_mpl_style
      -matplotlib.pyplot
      -numpy

&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
This class provides a GUI to prepare a set of .fits files that can be post-processed with the prefered softwares of the user.
1 - In a first frame :


 The user can define the object of interest (ROI) within the Messier, New General Catalog or International Catalog.
 The user can define the field of view (FOV) ot the region of interest in degree.
 On overview is exposed after having selected the object and the radius of the FOV.
 The region of interest is a square area in this first version.
 The metadata are downloaded for the considered ROI.
2 - In a second frame :

 The user can specify the expected resolution in arcseconds / pixel
 Based on that option a set of tiles is proposed considering a definition of 2400 px X 2400 px for each tile.
 The tiles are superimposed over the overview of the ROI
 The user can also select among the available set of data that are accessible within the diferrent serveys hosted on SkyView server
 The user specify the targeted directory where the data will be uploaded
 The format of the .fits is the following : 
 - For the tile centered on the selected object : 'center_'+survey_name+'_' + object_name + '.fits'
 - For the other tiles : 'tile#_'+survey_name+'_' + object_name + '.fits'
 A short report is exposed, in the prompt, on the screen that summarizes the .fits files that have been downloaded and their location on the computer with the the full path.
