
/* select all 2357 genes from ALN or ZOL CRISPRi screen with absolute value rho phenotype > 0.1 and p value <= 0.05  */
select a.gene, a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho  
from crispri.`aln` a
join 
crispri.`zol_digested_copy` z
on a.gene=z.gene
where 
(abs(a.`r avg3 average phenotype of strongest 3`) > 0.1 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.1) and (a.`r avg3 Mann-Whitney p-value` <= 0.05 or z.`r avg3 Mann-Whitney p-value` <= 0.05) and a.gene not like "%pseudo_%" and z.gene not like "%pseudo_%"

/* select all 774 statistically significant differentially expressed genes in both ONJ and DTC patients */
select o.symbol as gene, o.logFC as o_fold, o.`AveExpr` as o_avg, o.`P.Value` as o_pval from tbone.cheng_onj o where o.symbol !='.' and o.`P.Value` <= 0.05 group by symbol) o
join
(select  c.symbol as gene, c.logFC as c_fold, c.`AveExpr` as c_avg, c.`P.Value` as c_pval from tbone.cheng_cancer_her2_1000vs2500days c where c.symbol !='.' and c.`P.Value` <= 0.05 group by symbol) c
on o.gene=c.gene
group by o.gene
 

/* 165 genes CRISPRi-AFF intersection */
select * from (select * from aff_multiple_case_only_variants_v5) AFF
join
(select a.gene, a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho  
from crispri.`aln_copy` a
join 
crispri.`zol_digested_copy` z
on a.gene=z.gene
where 
(abs(a.`r avg3 average phenotype of strongest 3`) > 0.1 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.1) and (a.`r avg3 Mann-Whitney p-value` <= 0.05 or z.`r avg3 Mann-Whitney p-value` <= 0.05) and a.gene not like "%pseudo_%" and z.gene not like "%pseudo_%") CRISPRI
on AFF.gene=CRISPRI.gene
group by AFF.gene


/* 68 genes ONJ/DTC gene expression-AFF intersection */
select * from (select * from aff_multiple_case_only_variants_v5) AFF
join 
(select o.* from 
(select o.symbol as gene, o.logFC as o_fold, o.`AveExpr` as o_avg, o.`P.Value` as o_pval from tbone.cheng_onj o where o.symbol !='.' and o.`P.Value` <= 0.05 group by symbol) o
join
(select  c.symbol as gene, c.logFC as c_fold, c.`AveExpr` as c_avg, c.`P.Value` as c_pval from tbone.cheng_cancer_her2_1000vs2500days c where c.symbol !='.' and c.`P.Value` <= 0.05 group by symbol) c
on o.gene=c.gene
group by o.gene) RNA
on AFF.gene=RNA.gene 
group by AFF.gene

/* 5 genes ONJ/DTC gene expression-AFF-CRISPRi intersection */
select * from (select * from (select o.*, c_fold, c_avg, c_pval from 
(select o.symbol as gene, o.logFC as o_fold, o.`AveExpr` as o_avg, o.`P.Value` as o_pval from tbone.cheng_onj o where o.symbol !='.' and o.`P.Value` <= 0.05 group by symbol) o
join
(select  c.symbol as gene, c.logFC as c_fold, c.`AveExpr` as c_avg, c.`P.Value` as c_pval from tbone.cheng_cancer_her2_1000vs2500days c where c.symbol !='.' and c.`P.Value` <= 0.05 group by symbol) c
on o.gene=c.gene
group by o.gene) RNA
join
(select a.gene as crispri_gene, a.`r avg3 Mann-Whitney p-value` as crispri_aln_pval, a.`r avg3 average phenotype of strongest 3` as crispri_aln_rho,
z.`r avg3 Mann-Whitney p-value` as crispri_zol_pval, z.`r avg3 average phenotype of strongest 3` as crispri_zol_rho  
from crispri.`aln` a
join 
crispri.`zol_digested_copy` z
on a.gene=z.gene
where 
(abs(a.`r avg3 average phenotype of strongest 3`) > 0.1 or abs(z.`r avg3 average phenotype of strongest 3`) > 0.1) and (a.`r avg3 Mann-Whitney p-value` <= 0.05 or z.`r avg3 Mann-Whitney p-value` <= 0.05) and a.gene not like "%pseudo_%" and z.gene not like "%pseudo_%") CRISPRI
on RNA.gene=CRISPRI.crispri_gene
join
(select gene as aff_gene, patient1,patient2,patient3,patient4 from aff_multiple_case_only_variants_v5) AFF
on AFF.aff_gene=RNA.gene
group by RNA.gene) a
join aff_all_variants_v3_copy b
on b.gene like concat( '%', a.gene, '%')  