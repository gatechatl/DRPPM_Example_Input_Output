drppm -RemoveColumnsFromMatrix GSE98638_HCC.TCell.S5063.count.txt 0 GSE98638_HCC.TCell.S5063.count_filt1col.txt
drppm -MergeGeneName GSE98638_HCC.TCell.S5063.count_filt1col.txt MAX GSE98638_HCC.TCell.S5063.count_filt1col_max.txt
drppm -RemoveZeroCountGenes GSE98638_HCC.TCell.S5063.count_filt1col_max.txt 1 GSE98638_HCC.TCell.S5063.count_filt1col_max_zero.txt


drppm -RemoveColumnsFromMatrix GSE98638_HCC.TCell.S5063.TPM.txt 0 GSE98638_HCC.TCell.S5063.TPM_filt1col.txt
