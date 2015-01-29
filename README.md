# Mammut
So, I want a record of the scripts used at various points in my project analyzing gene expression information, and particularly it's divergence during evolution in maize. However, the project has been going on a long time, overwhich my programming style has developed a fair amount, and I don't have time now to make up for a the mess before. Therefore, this is a sadly unorganized compilation of scripts used for the maize paralog project. I'll try a drop a note about usage, but no promises. 

###general deconvolution scripts
-LinRegRNAseqData_11.R
-numerize_11.R

###scripts used to analyze enzyme and metabolite data

*function for 'deconvolution' error tolerant version
-linregsuperflex2.R
*running the normalization of enzymes from milli U activity/ mg FW to deconvolution/significanct testing, etc.
-pubnorm.R
*function to get mean and standard deviation for summed up tissues
-calc_absSE.R
