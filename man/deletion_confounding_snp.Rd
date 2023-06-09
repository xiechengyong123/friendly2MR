% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.R
\name{deletion_confounding_snp}
\alias{deletion_confounding_snp}
\title{调用phenoscanner包，查找工具变量相关的混杂因素}
\usage{
deletion_confounding_snp(
  confound = NULL,
  exposure_dat = NULL,
  query_gene = NULL,
  query_region = NULL,
  catalogue = "GWAS",
  pvalue = 5e-08,
  proxies = "None",
  r2 = 0.8,
  build = 37,
  write = TRUE,
  save_path = "MR_ivs_confounding.csv"
)
}
\arguments{
\item{confound}{设定的混杂因素，可以是一个或者多个，例如：c("Body mass index","Coronary heart disease")}

\item{exposure_dat}{TwosampleMR包中函数extract_instruments或者harmonise_data返回的数据框}

\item{query_gene}{设置phenoscanner函数搜索的基因，绝大部分情况下用不到该参数，默认：NULL}

\item{query_region}{设置phenoscanner函数搜索的基因组范围，绝大部分情况下用不到该参数，默认：NULL}

\item{catalogue}{设置phenoscanner函数搜索的目录（None, GWAS, eQTL, pQTl, mQTL, methQTL），默认：catalogue = GWAS}

\item{pvalue}{设置工具变量与混杂因素的筛选阈值，，默认：pvalue = 5e-08}

\item{proxies}{设置phenoscanner函数搜索的代理人种（None, AFR, AMR, EAS, EUR, SAS），默认：proxies = None}

\item{r2}{设置工具变量与其他可能与混杂因素相关连SNP的连锁不平衡阈值，大于该阈值的代理SNP得以保留，默认：r2 = 0.8}

\item{build}{设置基因组的版本（37，38），默认：build = 37}

\item{write}{逻辑值,将工具变量相关的混杂因素导出，默认：write=T 导出结果文件}

\item{save_path}{设置结果的保存路径和名称}
}
\value{
返回TwosampleMR包中exposure格式的数据框，若PhenoScanner网站中有筛选到工具变量SNP的混杂因素，则将该SNP的mr_keep.exposure/mr_keep状态设置为FALSE(⑅˃◡˂⑅)
}
\description{
调用phenoscanner包，查找工具变量相关的混杂因素
}
