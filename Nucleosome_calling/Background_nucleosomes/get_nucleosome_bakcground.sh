DPOSDIR='/Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Nucleosome_analysis/DANPOS'

cat $DPOSDIR/Heart_H3_positionning_aging_1e-30_2016-06-02/Heart_Heart_29m-Heart_Heart_3m.positions.integrative.xls                           | cut -f 1,2,3 > Heart_H3_nuc_background.bed
cat $DPOSDIR/OB_H3_positionning_aging_1e-30_2016-06-02/OlfactoryBulb_OB_29m-OlfactoryBulb_OB_3m.positions.integrative.xls                    | cut -f 1,2,3 > OlfactoryBulb_H3_nuc_background.bed
cat $DPOSDIR/Cerebellum_H3_positionning_aging_1e-30_2016-06-02/Cerebellum_Cerebellum_29m-Cerebellum_Cerebellum_3m.positions.integrative.xls  | cut -f 1,2,3 > Cerebellum_H3_nuc_background.bed
cat $DPOSDIR/NPCs_H3_positionning_aging_1e-30_2016-06-02/NPCs_NPCs_29m-NPCs_NPCs_3m.positions.integrative.xls                                | cut -f 1,2,3 > NPCs_H3_nuc_background.bed
cat $DPOSDIR/Liver_H3_positionning_aging_1e-30_2016-06-02/Liver_Liver_29m-Liver_Liver_3m.positions.integrative.xls                           | cut -f 1,2,3 > Liver_H3_nuc_background.bed



# 12644712 Cerebellum_H3_nuc_background.bed
# 12553854 Heart_H3_nuc_background.bed
# 12877928 Liver_H3_nuc_background.bed
# 11920670 NPCs_H3_nuc_background.bed
# 12017141 OlfactoryBulb_H3_nuc_background.bed
