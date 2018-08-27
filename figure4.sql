
/* 65 genes AFF + CRISPRi */
select * from (select o.gene, crispri_aln_pval, crispri_aln_rho, crispri_zol_pval, crispri_zol_rho
 from aff_multiple_case_only_variants_v5 o
join (select a.gene, a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho from crispri.zol_digested_copy z
join crispri.aln_copy a
on z.gene=a.gene where ((a.`r avg3 Mann-Whitney p-value` <= 0.05 and a.`r avg3 Mann-Whitney p-value` !='') or (z.`r avg3 Mann-Whitney p-value` <= 0.05 and z.`r avg3 Mann-Whitney p-value` !='')) and (abs(a.`r avg3 average phenotype of strongest 3`) > 0.15 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.15) group by a.gene) c
on
o.gene=c.gene
group by o.gene) c
join aff_all_variants_v3 vv
on c.gene=vv.gene
group by vv.gene


/* 16 genes AFF + CRISPRi + gene expression */
select o.symbol as gene, o.logFC as o_fold, o.`AveExpr` as o_avg, o.`P.Value` as o_pval, c.logFC as c_fold, c.`AveExpr` as c_avg, c.`P.Value` as c_pval,
a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho, v.*
 from cheng_onj o
join cheng_cancer_her2_1000vs2500days_copy c
on o.symbol = c.symbol
join crispri.zol_digested_copy z
on o.symbol=z.gene
join crispri.aln_copy a
on o.symbol=a.gene
join aff_multiple_case_only_variants_v5 v
on o.symbol = v.gene
where c.symbol !='.' and o.symbol !='.' and ((a.`r avg3 Mann-Whitney p-value` <= 0.05 and a.`r avg3 Mann-Whitney p-value` !='') or (z.`r avg3 Mann-Whitney p-value` <= 0.05 and z.`r avg3 Mann-Whitney p-value` !='')) and (c.`P.Value` <= 0.05 or o.`adj.P.Val` <= 0.05) and (abs(a.`r avg3 average phenotype of strongest 3`) > 0.15 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.15)
 order by abs(crispri_aln_rho) desc, abs(crispri_zol_rho) desc

/* the query to select the top 8 genes */
select o.symbol as gene, o.logFC as o_fold, o.`AveExpr` as o_avg, o.`P.Value` as o_pval, c.logFC as c_fold, c.`AveExpr` as c_avg, c.`P.Value` as c_pval,
a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho, v.*
 from cheng_onj o
join cheng_cancer_her2_1000vs2500days_copy c
on o.symbol = c.symbol
join crispri.zol_digested_copy z
on o.symbol=z.gene
join crispri.aln_copy a
on o.symbol=a.gene
join aff_all_variants_v3 v
on o.symbol = v.gene
where ((o.logFC < 0 and c.logFC < 0) or (o.logFC > 0 and c.logFC > 0) ) and c.symbol !='.' and o.symbol !='.' and ((a.`r avg3 Mann-Whitney p-value` <= 0.05 and a.`r avg3 Mann-Whitney p-value` !='') or (z.`r avg3 Mann-Whitney p-value` <= 0.05 and z.`r avg3 Mann-Whitney p-value` !='')) and ((a.`r avg3 average phenotype of strongest 3` < 0 and z.`r avg3 average phenotype of strongest 3` < 0) or (a.`r avg3 average phenotype of strongest 3` > 0 and z.`r avg3 average phenotype of strongest 3` > 0)) and v.`only_in_cases` = 1 and c.`P.Value` <= 0.05 and o.`adj.P.Val` <= 0.05 and (abs(a.`r avg3 average phenotype of strongest 3`) > 0.1 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.1)
group by o.symbol


/* incorporating haploid screen, which doesn't change any results */

select o.symbol as gene, o.logFC as o_fold, o.`AveExpr` as o_avg, o.`P.Value` as o_pval, c.logFC as c_fold, c.`AveExpr` as c_avg, c.`P.Value` as c_pval,
a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho, k.name as kbm7_gene, v.*
 from cheng_onj o
join cheng_cancer_her2_1000vs2500days_copy c
on o.symbol = c.symbol
join crispri.zol_digested_copy z
on o.symbol=z.gene
join crispri.aln_copy a
on o.symbol=a.gene
left join (
		select m.name from 
		(select gene_id, a.* 
		from tbone.aln_kbm7_hits a 
		join morpheome.aliases m
		on a.name=m.name) a
		left join `morpheome`.aliases m
		on m.gene_id=a.gene_id where m.type = "NCBI_official_symbol"
	
	) k
on 
o.symbol=k.name
join aff_all_variants_v3 v
on o.symbol = v.gene
where ((o.logFC < 0 and c.logFC < 0) or (o.logFC > 0 and c.logFC > 0) ) and c.symbol !='.' and o.symbol !='.' and ((a.`r avg3 Mann-Whitney p-value` <= 0.05 and a.`r avg3 Mann-Whitney p-value` !='') or (z.`r avg3 Mann-Whitney p-value` <= 0.05 and z.`r avg3 Mann-Whitney p-value` !='') or k.name is not null) and ((a.`r avg3 average phenotype of strongest 3` < 0 and z.`r avg3 average phenotype of strongest 3` < 0) or (a.`r avg3 average phenotype of strongest 3` > 0 and z.`r avg3 average phenotype of strongest 3` > 0) or k.name is not null) and v.`only_in_cases` = 1 and c.`P.Value` <= 0.05 and o.`adj.P.Val` <= 0.05 and (abs(a.`r avg3 average phenotype of strongest 3`) > 0.1 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.1) 
group by o.symbol