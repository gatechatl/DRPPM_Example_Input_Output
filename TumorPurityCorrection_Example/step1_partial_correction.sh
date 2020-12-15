drppm -GenerateScriptForPartiallingOut matrix.txt 1 tumor_purity_correction_script.r matrix_purity_corrected.txt

R --vanilla < tumor_purity_correction_script.r
