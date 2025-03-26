##########
## run SH python script image by image
## the script should be run in the image_env conda environment
##########

inputDir="/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/d4_WT"
outDir="res/test_organoid_WTd4_v2"

for file in ${inputDir}/*_C4_CystMask.tiff;
do
    echo $file
    python batch_SH_organoid.py $file $outDir
    
    #break
    
done

	
	   
	    


