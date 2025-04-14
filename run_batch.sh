##########
## run SH python script image by image
## the script should be run in the image_env conda environment
##########

inputDir="/Volumes/groups/tanaka/People/current/Hannah/RawImagesVienna/M44_OlympusSD/211208-11_HM1-extended-compTC/Isotropic/Isotropic_stainA_RAd2_forSH"

outDir="res/test_organoid_wt_timepoints_v1"

for file in ${inputDir}/*_C4_CystMask.tiff;
do
    echo $file

    python batch_SH_organoid_v2.py $file $outDir
    
    #break
    
done

	
	   
	    


