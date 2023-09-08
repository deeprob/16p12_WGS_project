


Instruction for running smoove (a wrapper for lumpy) are available at https://github.com/brentp/smoove


Smoove is installed at /data3/software/bin
It requires that other software is available on the PATH
so in every bash script, I set PATH to /data3/software/bin and 
/data3/software/lumpy-sv/bin 

Lumpy is installed at /data3/software/lumpy-sv/bin
For the regenotyping step, smoove also uses a software called duphold
I installed it at /data5/software/duphold from https://github.com/brentp/duphold/releases/tag/v0.2.3





# There are 4 steps
	1_smoove.sh is the initial smoove run. 
	2_union.sh combines all of the vcfs from step 1
	3_regonotype.sh regenotypes at every step
	4_paste.sh is an optional step that creates a single cohort vcf 



# there is filtering information at https://github.com/brentp/smoove/issues/130



