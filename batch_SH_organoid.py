from skimage import measure
import pandas as pd
from skimage.filters import threshold_otsu, rank
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np

import pyshtools as pysh

from skimage.io import imread
from skimage.io import imsave
from skimage import filters
from skimage import morphology
from pyclesperanto_prototype import imshow
import pyclesperanto_prototype as cle
import matplotlib.pyplot as plt

### import function https://github.com/gioda/FeARLesS/blob/main/Fearless/utils.py
import pyshtools
import gc
import numpy as np
import os
from vedo import printc, spher2cart, probePoints

#from vedo import Points
import shutil
from sys import exit
#from vedo import spher2cart
#from vedo import printc, spher2cart, probePoints

from sys import argv, exit
import numpy as np
import pyshtools

from vedo import printc, load, spher2cart, mag, ProgressBar, Points, write
from vedo import *

from skimage.io import imread
from skimage.io import imsave
#from skimage import filters
#from skimage import morphology
#from pyclesperanto_prototype import imshow
from skimage import measure
import pandas as pd
from skimage.filters import threshold_otsu, rank
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from vedo import printc, load, ProgressBar
#from utils import pathExists, voxelIntensity
import numpy as np
import gc
#import matplotlib
#matplotlib.use('Agg')  # Use Agg backend for non-interactive plotting
import matplotlib.pyplot as plt

#receiving input argument
input_image = sys.argv[1]
path_results = sys.argv[2]



path_results = path_results + '/'
if not os.path.exists(path_results):
    os.makedirs(path_results)

#pathExists(path_results)
lmax = 20
N = 500          # number of grid intervals on the unit sphere
rmax = 1400
x0 = [0, 0, 0]  # set object at this position
xLimb = [-200, 0, 200]
cutOrigin = [150, 0, 0]
deg_fit = 6
radiusDiscretisation = 50
N = 200
expo = 1.0

printc('lmax =', lmax, 'N =', N, 'deg_fit =', deg_fit, c='y')

from itertools import chain
col_names =  [['image', 'cyst_index', 'cyst_size', 'cyst_r'], ['power_per_l_m'+i for i in map(str, range(lmax+1))], ['power_per_l'+i for i in map(str, range(lmax+1))]]
col_names = list(chain(*col_names))
#col_names

#CPoutDir = "../images_data/CPouts2"
for nm in os.listdir(CPoutDir):
    if "_C4_CystMask.tiff" in nm:
        #print(nm)
        df = pd.DataFrame(columns = col_names)
        #df
        fileName = nm.replace('_C4_CystMask.tiff','')

        outDir = os.path.join(path_results, str(fileName + '/'))
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        #files += [nm]
        print(fileName)
        #make_composite_C5(newName, ImageDir, outDir)
        mask = imread(os.path.join(CPoutDir, str(fileName + "_C4_CystMask.tiff"))) # mask segmeted cysts
        C3 = imread(os.path.join(WTDir, str(fileName + "_C3.tif"))) # FoxA2, nuclei
        C2 = imread(os.path.join(WTDir, str(fileName + "_C2.tif"))) # Pax6, nuclei
        C4 = imread(os.path.join(WTDir, str(fileName + "_C4.tif"))) # Dapi, nuclei

        labels_mask, nb_cyst = measure.label(mask, return_num = True)

        # for each orgnoid
        for i in range(nb_cyst+1):
            if i > 0:

                cyst_select = np.where(labels_mask == i)
                cyst_index = i

                if cyst_select[0].size > 100: # orgnoid must > 100
                    # find the boundary of the orgnoid
                    [z, rows, columns] = np.where(labels_mask == cyst_index)
                    z1 = min(z)
                    z2 = max(z)
                    row1 = min(rows)
                    row2 = max(rows)
                    col1 = min(columns)
                    col2 = max(columns)

                    # crop the orgnoid and save foxa2 and dapi
                    xx = np.where(labels_mask == cyst_index, C3, 0)
                    newImage_foxa2 = xx[z1:z2, row1:row2, col1:col2]

                    xx = np.where(labels_mask == cyst_index, C4, 0)
                    newImage_dapi = xx[z1:z2, row1:row2, col1:col2]

                    xx = np.where(labels_mask == cyst_index, labels_mask, 0)
                    newImage_cyst = xx[z1:z2, row1:row2, col1:col2]

                    imsave(os.path.join(outDir, str(fileName + "_cyst_2_channel_FoxA2_" + str(cyst_index)+ ".tif")), newImage_foxa2)
                    imsave(os.path.join(outDir, str(fileName + "_cyst_2_channel_dapi_"+ str(cyst_index) + ".tif")), newImage_dapi)

                    del xx, newImage_dapi, newImage_cyst, newImage_foxa2, z1, z2, row1, row2, col1, col2

                    # reload the dapi and foxa2 using vedo
                    vol_foxa2 = load(os.path.join(outDir, str(fileName + "_cyst_2_channel_FoxA2_" + str(cyst_index)+ ".tif")))
                    vol_dapi = load(os.path.join(outDir, str(fileName + "_cyst_2_channel_dapi_"+ str(cyst_index) + ".tif")))

                    # the center and radius of organoid
                    pos = vol_dapi.center()
                    rmax = vol_dapi.diagonalSize()/2

                    # Compute voxel intensities
                    scalars = []
                    for th in np.linspace(0, np.pi, N, endpoint=False):
                        for ph in np.linspace(0, 2*np.pi, N, endpoint=False):

                            # compute sample points
                            p = spher2cart(rmax, th, ph)

                            samplePointsTmp = []
                            # making discretization more dense away from the center
                            p_tmp = p / (radiusDiscretisation-1)**expo

                            for j in range(radiusDiscretisation):
                                SP = pos + p_tmp * (j**expo)
                                samplePointsTmp.append(SP)

                            # compute intensities
                            pb = probePoints(vol_foxa2, samplePointsTmp)

                            scalarsTmp = pb.getPointArray()

                            scalars.append(scalarsTmp.max())

                            del pb, scalarsTmp, samplePointsTmp, p, p_tmp
                            #gc.collect()


                    matrixOfIntensities = np.array(scalars)
                    #coeff = matrixOfIntensities
                    #matrixOfIntensities
                    #len(matrixOfIntensities)
                    formattedcoeff = np.reshape(matrixOfIntensities, (N, N))

                    # remove the mean to reduce the l=0 contribution
                    formattedcoeff_norm = formattedcoeff - formattedcoeff.mean()

                    SH = pyshtools.SHGrid.from_array(formattedcoeff_norm, grid = "DH")
                    clm = SH.expand(lmax_calc= lmax)
                    #grid_reco = clm.expand(lmax=lmax).to_array()  # cut "high frequency" components
                    grid = clm.expand(lmax = 20)

                    Plot_intermediate = True
                    if Plot_intermediate:
                        fig, ax = plt.subplots(1, 1)
                        ax.imshow(formattedcoeff, extent=(0, N, 0, N))
                        ax.set(xlabel='theta', ylabel='phi')
                        #plt.savefig(os.path.join(outDir, str(fileName + "_FoxA2_intensity.pdf")))
                        plt.savefig(os.path.join(outDir, str(fileName + "_FoxA2_intensity_" + str(cyst_index)+".pdf")))

                        plt.clf()
                        plt.close()
                        gc.collect()

                        fig, ax = clm.plot_spectrum(show=False, lmax = 10, convention = 'power', unit = "per_lm")
                        plt.savefig(os.path.join(outDir, str(fileName + "_SH_power_per_lm_" + str(cyst_index)+".pdf")))
                        plt.clf()
                        plt.close()
                        gc.collect()

                        fig, ax = clm.plot_spectrum(show=False, lmax = 10, convention = 'power', unit = "per_l")
                        plt.savefig(os.path.join(outDir, str(fileName + "_SH_power_per_l_" + str(cyst_index)+ ".pdf")))
                        plt.clf()
                        plt.close()
                        gc.collect()

                        fig, ax = clm.plot_spectrum2d(show=False, lmax = 10)
                        plt.savefig(os.path.join(outDir, str(fileName + "_SH_power_2d_per_l_m_"+ str(cyst_index)+".pdf")))
                        plt.clf()
                        plt.close()
                        gc.collect()


                        fig, ax = grid.plot(show=False)
                        plt.savefig(os.path.join(outDir, str(fileName + "_FoxA2Intensity_restored_with20Components_"+ str(cyst_index)+".pdf")))
                        plt.clf()
                        plt.close()
                        gc.collect()

                        plt.close('all')
                        gc.collect();

                    del scalars, matrixOfIntensities, formattedcoeff
                    #gc.collect()

                    # save the SH power
                    df.loc[len(df)] =  list(chain([fileName, cyst_index, cyst_select[0].size, rmax],
                                              clm.spectrum(lmax = lmax, convention = 'power', unit = "per_lm"),
                                              clm.spectrum(lmax = lmax, convention = 'power', unit = "per_l")))

                    del SH, clm, grid, cyst_select, vol_foxa2, vol_dapi
                    gc.collect()
                    print('gc collection done for cyst ' + str(i))

        df.to_csv(os.path.join(outDir, "Condion_eachOrgnoid_SHspectrum_power.csv"), index=True, header=True)
        del mask, C2, C3, C4, df, labels_mask, nb_cyst
        gc.collect()
        print('gc collection done for image' + fileName)
