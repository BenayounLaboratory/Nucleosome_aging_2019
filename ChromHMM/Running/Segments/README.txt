# 2016-12-22

# separate segments into different files for intersection analysis

parseSegments.pl Cerebellum Cerebellum_10_OBLiverHeartCerebellumNPCs_mm9_HMM_segments.bed
parseSegments.pl Heart Heart_10_OBLiverHeartCerebellumNPCs_mm9_HMM_segments.bed
parseSegments.pl Liver Liver_10_OBLiverHeartCerebellumNPCs_mm9_HMM_segments.bed
parseSegments.pl NPCs NPCs_10_OBLiverHeartCerebellumNPCs_mm9_HMM_segments.bed
parseSegments.pl OlfactoryBulb OlfactoryBulb_10_OBLiverHeartCerebellumNPCs_mm9_HMM_segments.bed


E1    =>  Low Signal
E2    =>  Active Enhancer
E3    =>  Weak/poised enhancer
E4    =>  Inactive poised/promoter
E5    =>  Polycomb repressed
E6    =>  Low Signal
E7    =>  Insulator
E8    =>  Weak/poised enhancer
E9    =>  Active Promoter
E10   =>  Flanking Promoter
