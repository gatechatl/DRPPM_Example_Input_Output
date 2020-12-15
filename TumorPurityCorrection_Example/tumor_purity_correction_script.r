data = read.table("matrix.txt", sep="\t", header = T, row.names=1);
r2 = resid(m2 <- lm(HSC.Prog_Normal.derived.combined_PMID_30827681 ~ BM_Blast, data = data))
HSC.Prog_Normal.derived.combined_PMID_30827681 = coef(m2)["(Intercept)"] + r2
r3 = resid(m3 <- lm(GMP_Normal.derived.combined_PMID_30827681 ~ BM_Blast, data = data))
GMP_Normal.derived.combined_PMID_30827681 = coef(m3)["(Intercept)"] + r3
r4 = resid(m4 <- lm(Myeloid_Normal.derived.combined_PMID_30827681 ~ BM_Blast, data = data))
Myeloid_Normal.derived.combined_PMID_30827681 = coef(m4)["(Intercept)"] + r4
r5 = resid(m5 <- lm(HSC.Prog.like_Tumor.derived.combined_PMID_30827681 ~ BM_Blast, data = data))
HSC.Prog.like_Tumor.derived.combined_PMID_30827681 = coef(m5)["(Intercept)"] + r5
r6 = resid(m6 <- lm(GMP.like_Tumor.derived.combined_PMID_30827681 ~ BM_Blast, data = data))
GMP.like_Tumor.derived.combined_PMID_30827681 = coef(m6)["(Intercept)"] + r6
r7 = resid(m7 <- lm(Myeloid.like_Tumor.derived.combined_PMID_30827681 ~ BM_Blast, data = data))
Myeloid.like_Tumor.derived.combined_PMID_30827681 = coef(m7)["(Intercept)"] + r7
r8 = resid(m8 <- lm(HSC.like_Tumor.derived.per.cell.type_PMID_30827681 ~ BM_Blast, data = data))
HSC.like_Tumor.derived.per.cell.type_PMID_30827681 = coef(m8)["(Intercept)"] + r8
out = cbind(HSC.Prog_Normal.derived.combined_PMID_30827681,GMP_Normal.derived.combined_PMID_30827681,Myeloid_Normal.derived.combined_PMID_30827681,HSC.Prog.like_Tumor.derived.combined_PMID_30827681,GMP.like_Tumor.derived.combined_PMID_30827681,Myeloid.like_Tumor.derived.combined_PMID_30827681,HSC.like_Tumor.derived.per.cell.type_PMID_30827681)
write.table(out, file = "matrix_purity_corrected.txt", sep="\t", col.names=NA)
