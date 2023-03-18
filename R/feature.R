#' 查询多个暴露对一个结局的阳性结果
#'
#' @export
#' @return data frame
#查询多个暴露对一个结局的阳性结果
find_anyexposur_outcome<-function(exposure=exposure,outcome=outcome,write=F,p1 = 5e-08,clump = TRUE,p2 = 5e-08,r2= 0.001,kb = 10000,LD=0.8){
  #设置MR_result_combine为空向量  
  MR_result_combine<-c()
  start=Sys.time()
  #对每一个暴露进行循环分析
  for (id in 1:length(exposure)) {
    exposure_id=exposure[id]
    
    #计算还剩下几个暴露
    remain=length(exposure)-id
    print(paste0("当前分析的是第",id,"个暴露,","还剩下",remain,"个暴露"))
    
    #提取IV
    exposure_dat <- extract_instruments(
      exposure_id,
      p1 = p1,
      clump = clump,
      p2 = p2,
      r2 = r2,
      kb = kb,
      access_token = NULL
    )
    
    #若暴露P值太小，则无法进行MR分析
    if (is.null(exposure_dat)) {
      res<-cbind(exposure_id,outcome,"筛选暴露P值太小，无法进行MR分析")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #从结局中获取有效的工具变量
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = outcome,rsq =LD ,access_token = NULL)
    
    #若在结局中无有效的IV，则需重新设置LD，然后在查找合适的工具变量
    if (is.null(outcome_dat)) {
      res<-cbind(exposure_id,outcome,"需要重新设定LD阈值，否则无法进行MR分析")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #协调和矫正数据
    dat <- harmonise_data(exposure_dat, outcome_dat)
    
    
    #执行MR分析
    res <- mr(dat)
    
    if (sum(dat$mr_keep)==0) {
      res=cbind(exposure_id,outcome,"初筛结果无阳性")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #若初筛过程中的5种方法中有显著的结果则进行记录
    if (sum(res$pval<0.05)>=1) {
      res<-cbind(exposure_id,outcome,"初筛结果阳性")
      MR_result_combine<-rbind(MR_result_combine,res)
    }else{
      res<-cbind(exposure_id,outcome,"初筛结果无阳性")
      MR_result_combine<-rbind(MR_result_combine,res)
    }
    
  }
  
  #将MR_result_combine设置为数据狂格式然后重新设置列名
  MR_result_combine=data.frame(MR_result_combine)
  colnames(MR_result_combine)<-c("id.exposure","id.outcome","conclusion")
  if (write) {
    write.csv(MR_result_combine,"MR多对一初筛.csv",row.names = F)
  }
  
  end=Sys.time()
  print(end-start)
  return(MR_result_combine)
  
}

#' 查询一个暴露对多个结局的阳性结果
#' @export
#' @return data frame
#查询一个暴露对多个结局的阳性结果
find_exposur_anyoutcome<-function(exposure=exposure,outcome=outcome,write=F,p1 = 5e-08,clump = TRUE,p2 = 5e-08,r2= 0.001,kb = 10000,LD=0.8){
  #设置MR_result_combine为空向量
  MR_result_combine<-c()
  starttime=Sys.time()
  #对每一个暴露进行循环分析
  for (id in 1:length(outcome)) {
    outcome_id=outcome[id]
    
    #计算还剩下几个暴露
    remain=length(outcome)-id
    print(paste0("当前分析的是第",id,"个结局,","还剩下",remain,"个结局"))
    
    #提取IV
    exposure_dat <- extract_instruments(
      exposure,
      p1 = p1,
      clump = clump,
      p2 = p2,
      r2 = r2,
      kb = kb,
      access_token = NULL
    )
    
    #若暴露P值太小，则无法进行MR分析
    if (is.null(exposure_dat)) {
      res<-cbind(exposure,outcome_id,"筛选暴露P值太小，无法进行MR分析")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #从结局中获取有效的工具变量
    outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = outcome_id,rsq =LD ,access_token = NULL)
    
    #若在结局中无有效的IV，则需重新设置LD，然后在查找合适的工具变量
    if (is.null(outcome_dat)) {
      res<-cbind(exposure,outcome_id,"需要重新设定LD阈值，否则无法进行MR分析")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #协调和矫正数据
    dat <- harmonise_data(exposure_dat, outcome_dat)
    
    #执行MR分析
    res <- mr(dat)
    
    if (sum(dat$mr_keep)==0) {
      res=cbind(exposure,outcome_id,"初筛结果无阳性")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #若初筛过程中的5种方法中有显著的结果则进行记录
    if (sum(res$pval<0.05)>=1) {
      res<-cbind(exposure,outcome_id,"初筛结果阳性")
      MR_result_combine<-rbind(MR_result_combine,res)
    }else{
      res<-cbind(exposure,outcome_id,"初筛结果无阳性")
      MR_result_combine<-rbind(MR_result_combine,res)
    }
    
  }
  
  #将MR_result_combine设置为数据狂格式然后重新设置列名
  MR_result_combine=data.frame(MR_result_combine)
  colnames(MR_result_combine)<-c("id.exposure","id.outcome","conclusion")
  if (write) {
    write.csv(MR_result_combine,"MR一对多初筛.csv",row.names = F)
  }
  
  endtime=Sys.time()
  print(endtime-starttime)
  
  return(MR_result_combine)
  
}


#' 从工具变量中插补缺失的SNP(1000基因组数据)
#' @export
#' @return data frame
#查询一个暴露对多个结局的阳性结果
#从工具变量中插补缺失的SNP(1000基因组数据)
snp_add_eaf=function(exposure=NULL,outcome=NULL,pop="EUR"){
  
  if(!is.null(exposure)&!is.null(outcome)){
    stop("只能输入暴露或者结局，二选其一")
  }
  
  if (!is.null(exposure)) {
    cat("###############################\n#########暴露数据eaf插补#########\n###############################")
    for(i in 1:nrow(exposure)){
      if (is.na(exposure$eaf.exposure[i])) {
        #the alternative allele frequencies and the LD scores for each variant, calculated for each super population separately
        snp_info=as.data.frame(ieugwasr::afl2_rsid(exposure$SNP[i]))
        if (nrow(snp_info)==0) {
          next()
        }
        
        if (pop=="EUR") {
          eaf=snp_info$AF.EUR  
        }else if (pop=="AFR") {
          eaf=snp_info$AF.AFR  
        }else if (pop=="AMR") {
          eaf=snp_info$AF.AMR
        }else if (pop=="EAS") {
          eaf=snp_info$AF.EAS
        }else if (pop=="SAS"){
          eaf=snp_info$AF.SAS
        }else{
          cat("pop 输入错误,只能输入：c(EUR,AFR,AMR,EAS,SAS)")
        }
        if (is.null(eaf)) {
          next()
        }
        exposure$eaf.exposure[i]=eaf
      }
    }
    return(exposure)
  } 
  
  
  if (!is.null(outcome)) {
    cat("###############################\n#########结局数据eaf插补#########\n###############################")
    for(i in 1:nrow(outcome)){
      if (is.na(outcome$eaf.outcome[i])) {
        #the alternative allele frequencies and the LD scores for each variant, calculated for each super population separately
        snp_info=as.data.frame(ieugwasr::afl2_rsid(outcome$SNP[i]))
        if (nrow(snp_info)==0) {
          next()
        }
        
        if (pop=="EUR") {
          eaf=snp_info$AF.EUR  
        }else if (pop=="AFR") {
          eaf=snp_info$AF.AFR  
        }else if (pop=="AMR") {
          eaf=snp_info$AF.AMR
        }else if (pop=="EAS") {
          eaf=snp_info$AF.EAS
        }else if (pop=="SAS"){
          eaf=snp_info$AF.SAS
        }else{
          cat("pop 输入错误,只能输入：(EUR,AFR,AMR,EAS,SAS)中的一种")
        }
        if (is.null(eaf)) {
          next()
        }
        outcome$eaf.outcome[i]=eaf
        
      }
    }
    return(outcome)
  }     
  
}

#' 定义函数：查找工具变量相关的混杂因素 
#'
#' @export
#' @return data frame
#定义函数：查找工具变量相关的混杂因素
deletion_confounding_snp = function(confound = NULL,
                                    exposure_dat = NULL,
                                    query_gene = NULL,
                                    query_region = NULL,
                                    catalogue = "GWAS",
                                    pvalue = 1e-05,
                                    proxies = "None",
                                    r2 = 0.8,
                                    build = 37,
                                    write=TRUE,
                                    filepath="MR_ivs_confounding.csv") {
  #将多个混杂因素进行合并，以|为分隔符
  confound = paste(confound, collapse = "|")
  
  #设置MR_result_combine为空向量  
  MR_result_combine<-c()
  #循环进行phenoscanner函数查找与query_snp相关的GWAS结果
  for (i in 1:length(exposure_dat$SNP)) {
    #获取SNP的rsID号
    query_snp = exposure_dat$SNP[i]
    #显示进程程度
    cat(
      "#################################################################################\n"
    )
    num = length(exposure_dat$SNP) - i
    print = paste0("查询第", i, "个工具变量的混杂因素", "，还剩下", num, "个工具变量")
    cat(paste0("###############", print, "########################"))
    cat(
      "\n#################################################################################\n"
    )
    
    #phenoscanner函数查找与query_snp相关的混杂因素
    phenoscanner_data <-
      phenoscanner::phenoscanner(
        snpquery = query_snp,
        genequery = query_gene,
        regionquery = query_region,
        catalogue = catalogue,
        pvalue = pvalue,
        proxies = proxies,
        r2 = r2,
        build = build
      )
    #获取工具变量SNP的信息
    phe_snps = phenoscanner_data$snps
    
    #获取phenoscanner函数查找结果
    phe_results = phenoscanner_data$results
    
    #若PhenoScanner网站中有筛选到工具变量SNP的混杂因素，则将该SNP的mr_keep.exposure/mr_keep状态设置为FALSE
    if (sum(("mr_keep" %in% colnames(exposure_dat)) == 1)) {
      if (sum(grepl(confound, phe_results$trait)) >= 1) {
        exposure_dat$mr_keep[i] =as.logical("FALSE")
        res<-data.frame(cbind(query_snp,phe_results$trait,phe_results$p))
        res<-res[grepl(confound,res$V2),]
        MR_result_combine<-rbind(MR_result_combine,res)
      }
    } else {
      if (sum(grepl(confound, phe_results$trait)) >= 1) {
        exposure_dat$mr_keep.exposure[i] = as.logical("FALSE")
        res<-data.frame(cbind(query_snp,phe_results$trait,phe_results$p))
        res<-res[grepl(confound,res$V2),]
        MR_result_combine<-rbind(MR_result_combine,res)
      }
    }
  }
  #将MR_result_combine设置为数据框格式然后重新设置列名
  MR_result_combine=data.frame(MR_result_combine)
  colnames(MR_result_combine)<-c("工具变量","混杂因素","P值")
  if (write) {
    write.csv(MR_result_combine,file = filepath,row.names = F)
  }  
  #返回暴露/矫正后的数据
  return(exposure_dat)
  
}

