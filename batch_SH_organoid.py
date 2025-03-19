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

def samplePoints(vol, expo, N, radiusDiscretisation):
    """Compute sample points."""
    pos = vol.center()
    rmax = vol.diagonalSize()/2

    samplePoints = []
    for th in np.linspace(0, np.pi, N, endpoint=False):
        for ph in np.linspace(0, 2*np.pi, N, endpoint=False):

            # compute sample points
            p = spher2cart(rmax, th, ph)
            # making discretization more dense away from the center
            p_tmp = p / (radiusDiscretisation-1)**expo
            for j in range(radiusDiscretisation):
                SP = pos + p_tmp * (j**expo)
                samplePoints.append(SP)

    del vol
    return np.array(samplePoints)

def confirm(message):
    """
    Ask user to enter Y or N (case-insensitive).

    :return: True if the answer is Y.
    :rtype: bool
    """
    answer = ""
    while answer not in ["y", "n"]:
        answer = input(message).lower()
    return answer == "y"


def pathExists(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        printc("Directory ", path, " Created ", c='green')
    else:
        printc("Directory ", path, " already exists", c='red')
        if confirm("Should I delete the folder and create a new one [Y/N]? "):
            shutil.rmtree(path)
            os.makedirs(path, exist_ok=True)
            printc("Directory ", path, " Created ", c='green')
        else:
            exit()


def voxelIntensity(vol, expo, N, radiusDiscretisation):
    """Compute voxel intensities."""
    pos = vol.center()
    rmax = vol.diagonalSize()/2

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
            pb = probePoints(vol, samplePointsTmp)

            del samplePointsTmp

            # making the intensities growing outside the volume according to the gradient
            scalarsTmp = pb.getPointArray()
            nonz = np.nonzero(scalarsTmp)[0]
            if len(nonz) > 2:
                lastNoZeroId = nonz[-1]  # find the last value != 0
                secondlastNoZeroId = nonz[-2]
                # find the last value != 0
                lastNoZero = scalarsTmp[lastNoZeroId]
                secondlastNoZero = scalarsTmp[secondlastNoZeroId]
                dx = lastNoZero - secondlastNoZero

            for i in range(lastNoZeroId+1, len(scalarsTmp)):
                scalarsTmp[i] = scalarsTmp[i-1] + dx
            scalars.append(scalarsTmp.tolist())

            del pb, scalarsTmp

    del vol
    gc.collect()

    # return allIntensitiesMatrix
    return np.array(scalars).reshape((N * N, radiusDiscretisation))


def forwardTransformation(matrixOfIntensities, N, lmax):

    ##############################################

    coeff = matrixOfIntensities

    ##############################################
    # SPHARNM
    allClm = np.zeros((matrixOfIntensities.shape[1], 2, lmax, lmax))
    for j in range(allClm.shape[0]):
        formattedcoeff = np.reshape(coeff[:, j], (N, N))
        SH = pyshtools.SHGrid.from_array(formattedcoeff)
        clm = SH.expand()

        allClm[j, :, :, :] = clm.to_array(lmax=lmax - 1)

    del formattedcoeff, clm, matrixOfIntensities

    return allClm


def inverseTransformations(allClm, allIntensitiesShape, N, lmax):
    """Make inverse SPHARM."""
    from scipy.interpolate import griddata

    aSH_recoMatrix = np.zeros((allIntensitiesShape[0], allIntensitiesShape[1]))

    for j in range(allClm.shape[0]):
        # inverse SPHARM coefficients
        clmCoeffs = pyshtools.SHCoeffs.from_array(allClm[j, :, :, :])
        SH_reco = clmCoeffs.expand(lmax=lmax - 1)
        # grid_reco.plot()
        aSH_reco = SH_reco.to_array()

        ##############################
        pts1 = []
        ll = []
        for ii, long in enumerate(np.linspace(0, 360, num=aSH_reco.shape[1], endpoint=True)):
            for jj, lat in enumerate(np.linspace(90, -90, num=aSH_reco.shape[0], endpoint=True)):
                th = np.deg2rad(90 - lat)
                ph = np.deg2rad(long)
                p = spher2cart(aSH_reco[jj][ii], th, ph)
                pts1.append(p)
                ll.append((lat, long))

        radii = aSH_reco.T.ravel()

        # make a finer grid
        n = N * 1j
        l_min, l_max = np.array(ll).min(axis=0), np.array(ll).max(axis=0)
        grid = np.mgrid[l_max[0]:l_min[0]:n, l_min[1]:l_max[1]:n]
        grid_x, grid_y = grid
        agrid_reco_finer = griddata(ll, radii, (grid_x, grid_y), method='cubic')
        ##############################

        formatted_aSH_reco = np.reshape(agrid_reco_finer, (N * N))

        aSH_recoMatrix[:, j] = formatted_aSH_reco

    del formatted_aSH_reco, agrid_reco_finer, grid_x, grid_y, grid

    return aSH_recoMatrix

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
#import pyclesperanto_prototype as cle
#import matplotlib.pyplot as plt

#from scipy.interpolate import griddata


lmax = 20
N = 500          # number of grid intervals on the unit sphere
rmax = 1400
x0 = [0, 0, 0]  # set object at this position
xLimb = [-200, 0, 200]
cutOrigin = [150, 0, 0]
deg_fit = 6

CPoutDir = "/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/test_WTd6"
WTDir = "/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/test_WTd6"

#CPoutDir = "/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/CPouts2"
#ImageDir = "/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/d4_10x_Pax6KO_WTchim"

#DataPath = '/Users/jingkui.wang/workspace/imp/image_analysis/S-BIAD441/limbs/limbs-noFlank/'
path_results = 'res/' + 'testScript_organoid_WTd6' + '/'

#pathExists(path_results)

printc('lmax =', lmax, 'N =', N, 'deg_fit =', deg_fit, c='y')

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

radiusDiscretisation = 50
N = 200
FFTexpansion = radiusDiscretisation
expo = 1.0
#print(outDir)

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
        
