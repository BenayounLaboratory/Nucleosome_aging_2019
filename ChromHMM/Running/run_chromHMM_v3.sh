# 2016-12-14

java -Xmx8G -jar /Users/benayoun/Softwares/ChromHMM_1.12/ChromHMM.jar BinarizeBed \
	/Users/benayoun/Softwares/ChromHMM_1.12/CHROMSIZES/mm9.txt \
	/Volumes/MyBook_3/BD_aging_project/Public_datasets/LICR_Datasets/ChromHMM/CHROM_HMM_RUN_v3/Input_Dir \
	LiverHeartCerebellumNPCs_fileable3.txt \
	/Volumes/MyBook_3/BD_aging_project/Public_datasets/LICR_Datasets/ChromHMM/CHROM_HMM_RUN_v3/Output_bin_dir

java -Xmx8G -jar /Users/benayoun/Softwares/ChromHMM_1.12/ChromHMM.jar LearnModel -i OBLiverHeartCerebellumNPCs_mm9_HMM -p 2 \
	/Volumes/MyBook_3/BD_aging_project/Public_datasets/LICR_Datasets/ChromHMM/CHROM_HMM_RUN_v3/Output_bin_dir \
	/Volumes/MyBook_3/BD_aging_project/Public_datasets/LICR_Datasets/ChromHMM/CHROM_HMM_RUN_v3/learned_model_OBLiverHeartCerebellumNPCs \
	10 \
	mm9

