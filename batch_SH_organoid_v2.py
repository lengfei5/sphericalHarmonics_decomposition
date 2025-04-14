

import numpy as np
import pandas as pd

import pyshtools as pysh

### import function https://github.com/gioda/FeARLesS/blob/main/Fearless/utils.py
import os

#from vedo import Points
import shutil
from sys import argv, exit

from vedo import printc, load, spher2cart, mag, ProgressBar, Points, write, probePoints
from vedo import *

from scipy import ndimage as ndi

from skimage.io import imread
from skimage.io import imsave
from skimage import measure
from skimage.filters import rank
from skimage.filters import threshold_otsu as filter_otsu
#from skimage import filters
from skimage import morphology
from skimage.morphology import disk
from skimage.feature import peak_local_max
from skimage.filters import gaussian

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

from pyclesperanto_prototype import imshow
import pyclesperanto_prototype as cle
import napari_segment_blobs_and_things_with_membranes as nsbatwm
import napari_simpleitk_image_processing as nsitk




import gc
from itertools import chain

#receiving input argument
import sys
input_image = sys.argv[1]
path_results = sys.argv[2]


print("input image -- " + input_image)
inputDir = os.path.dirname(input_image)
name_image = os.path.basename(input_image)

print("input directory -- " + inputDir)


path_results = path_results + '/'

if not os.path.exists(path_results):
    os.makedirs(path_results)

print("result directory -- " + path_results)



#pathExists(path_results)
N = 300
lmax = 100
radiusDiscretisation = 50
expo = 1.0

printc('lmax =', lmax, 'N =', N, c='y')

#col_names =  [['image', 'cyst_index', 'cyst_size', 'cyst_r'], ['power_per_lm_'+i for i in map(str, range(lmax+1))], ['power_per_l_'+i for i in map(str, range(lmax+1))]]
#col_names = list(chain(*col_names))
col_names =  [['image', 'cyst_index', 'cyst_size', 'cyst_r', 'nb_local_max', 'nb_segmentedObj'],
              ['power_per_lm_'+i for i in map(str, range(lmax+1))],
              ['power_per_l_'+i for i in map(str, range(lmax+1))],
              ['power_per_dlogl_'+i for i in map(str, range(lmax+1))]
             ]
col_names = list(chain(*col_names))


if "_C4_CystMask.tiff" in name_image:
    #print(nm)
    df = pd.DataFrame(columns = col_names)
    #df
    fileName = name_image.replace('_C4_CystMask.tiff','')

    outDir = os.path.join(path_results, str(fileName + '/'))
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    print("-------------")
    print("file name -- " + fileName)
    print("output directory -- " + outDir)
    print("-------------")
    #files += [nm]
    #print(fileName)
    #make_composite_C5(newName, ImageDir, outDir)
    mask = imread(os.path.join(inputDir, str(fileName + "_C4_CystMask.tiff"))) # mask segmeted cysts
    C3 = imread(os.path.join(inputDir, str(fileName + "_C3.tif"))) # FoxA2, nuclei
    #C2 = imread(os.path.join(inputDir, str(fileName + "_C2.tif"))) # Pax6, nuclei
    C4 = imread(os.path.join(inputDir, str(fileName + "_C4.tif"))) # Dapi, nuclei

    labels_mask, nb_cyst = measure.label(mask, return_num = True)

    # for each orgnoid
    for i in range(nb_cyst+1):
        if i > 0:

            cyst_select = np.where(labels_mask == i)
            cyst_index = i

            if cyst_select[0].size > 500: # orgnoid must > 100
                # find the boundary of the orgnoid
                [z, rows, columns] = np.where(labels_mask == cyst_index)
                z1 = min(z)
                z2 = max(z)
                row1 = min(rows)
                row2 = max(rows)
                col1 = min(columns)
                col2 = max(columns)

                df = pd.DataFrame(columns = col_names)

                # crop the orgnoid and save foxa2 and dapi
                xx = np.where(labels_mask == cyst_index, C3, 0)
                newImage_foxa2 = xx[z1:z2, row1:row2, col1:col2]

                xx = np.where(labels_mask == cyst_index, C4, 0)
                newImage_dapi = xx[z1:z2, row1:row2, col1:col2]

                #xx = np.where(labels_mask == cyst_index, labels_mask, 0)
                #newImage_cyst = xx[z1:z2, row1:row2, col1:col2]

                imsave(os.path.join(outDir, str(fileName + "_cyst_2_channel_FoxA2_" + str(cyst_index)+ ".tif")), newImage_foxa2)
                imsave(os.path.join(outDir, str(fileName + "_cyst_2_channel_dapi_"+ str(cyst_index) + ".tif")), newImage_dapi)

                del xx, newImage_dapi, newImage_foxa2, z1, z2, row1, row2, col1, col2

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
                        nonz = np.nonzero(scalarsTmp)[0]
                        scalarsTmp_nonz = scalarsTmp[nonz]

                        scalars.append(scalarsTmp_nonz.max())

                        del pb, scalarsTmp, samplePointsTmp, p, p_tmp
                        #gc.collect()


                matrixOfIntensities = np.array(scalars)
                formattedcoeff = np.reshape(matrixOfIntensities, (N, N))

                ## smoothing the spherical grid for local maximum detection
                #formattedcoeff = gaussian(formattedcoeff, sigma=2, preserve_range=True)
                preprocessed = gaussian(formattedcoeff, sigma=4, preserve_range=True)
                ## find the cutoff of FoxA2
                binary = nsbatwm.threshold_otsu(preprocessed)

                #imshow(binary)
                split_objects = nsbatwm.split_touching_objects(binary, sigma=3.5)
                #imshow(split_objects)

                #touching_labels = nsitk.touching_objects_labeling(binary, sigma=6)

                mask_grid, nb_labels = measure.label(split_objects, return_num = True)

                cutoff_foxa2 = filter_otsu(preprocessed)
                print(cutoff_foxa2)

                coordinates = peak_local_max(preprocessed, threshold_abs=cutoff_foxa2, min_distance=10)

                fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
                ax = axes.ravel()

                ax[0].imshow(formattedcoeff)
                ax[0].axis('off')
                ax[0].set_title('Original')

                ax[1].imshow(preprocessed, cmap=plt.cm.gray)
                ax[1].axis('off')
                ax[1].set_title('processed image')

                ax[2].imshow(formattedcoeff, cmap=plt.cm.gray)
                ax[2].autoscale(False)
                ax[2].plot(coordinates[:, 1], coordinates[:, 0], 'r.')
                ax[2].axis('off')
                ax[2].set_title('Peak local max')

                ax[3].imshow(split_objects)
                ax[3].autoscale(False)
                #ax[3].plot(coordinates[:, 1], coordinates[:, 0], 'r.')
                ax[3].axis('off')
                ax[3].set_title('watershed segmentation')

                fig.tight_layout()

                plt.savefig(os.path.join(outDir, str(fileName + "_" + str(cyst_index) +  "locaMaximum_watershedSegmentation.pdf")))
                #plt.show()


                # remove the mean to reduce the l=0 contribution
                #formattedcoeff_norm = formattedcoeff - formattedcoeff.mean()
                formattedcoeff_norm = (formattedcoeff - formattedcoeff.mean())/np.var(formattedcoeff)

                SH = pysh.SHGrid.from_array(formattedcoeff_norm, grid = "DH")
                clm = SH.expand(lmax_calc= lmax)
                grid = clm.expand(lmax = lmax)

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

                    fig, ax = clm.plot_spectrum(show=False, lmax = 20, convention = 'power', unit = "per_l")
                    plt.savefig(os.path.join(outDir, str(fileName + "_SH_power_per_l_" + str(cyst_index)+ ".pdf")))
                    plt.clf()
                    plt.close()
                    gc.collect()

                    fig, ax = clm.plot_spectrum(show=False, lmax = 20, convention = 'power', unit = "per_lm")
                    plt.savefig(os.path.join(outDir, str(fileName + "_SH_power_per_lm_" + str(cyst_index)+".pdf")))
                    plt.clf()
                    plt.close()
                    gc.collect()

                    fig, ax = clm.plot_spectrum(show=False, lmax = 20, convention = 'power', unit = "per_dlogl")
                    plt.savefig(os.path.join(outDir, str(fileName + "_SH_power_per_dlogl_" + str(cyst_index)+ ".pdf")))
                    plt.clf()
                    plt.close()
                    gc.collect()

                    plt.close('all')
                    gc.collect();

                #coefs = clm.spectrum(lmax = 50, convention = 'power', unit = "per_l")

                del scalars, matrixOfIntensities, formattedcoeff
                #gc.collect()

                # save the SH power
                df.loc[len(df)] =  list(chain([fileName, cyst_index, cyst_select[0].size, rmax, len(coordinates), nb_labels],
                                          clm.spectrum(lmax = lmax, convention = 'power', unit = "per_lm"),
                                          clm.spectrum(lmax = lmax, convention = 'power', unit = "per_l"),
                                 clm.spectrum(lmax = lmax, convention = 'power', unit = "per_dlogl")))


                df.to_csv(os.path.join(outDir, str("SHspectrum_power_" + fileName + "_" + str(cyst_index) + ".csv")), index=True, header=True)

                del SH, clm, grid, cyst_select, vol_foxa2, vol_dapi, df
                gc.collect()
                print('gc collection done for cyst ' + str(i))

    #df.to_csv(os.path.join(outDir, "Condion_eachOrgnoid_SHspectrum_power.csv"), index=True, header=True)
    del mask, C3, C4, labels_mask, nb_cyst
    gc.collect()
    print('gc collection done for image' + fileName)

else:
    print("input image format is correct...")
    sys.exit(0)
