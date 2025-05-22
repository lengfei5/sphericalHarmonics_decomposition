##########
## run SH python script image by image
## the script should be run in the image_env conda environment
##########

inputDir="/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/d3-4-6_10x_PKO-FKO_TetOnF-TetOnP_forSH"

outDir="res/test_organoid_KOKO_TetOnTetOn_timepoints_v1"

for file in ${inputDir}/*_C5_CystMask.tiff;
do
    echo $file

    python batch_SH_organoid_koko_tetonteton.py $file $outDir
    
    #break
    
done

	
	   
	    


