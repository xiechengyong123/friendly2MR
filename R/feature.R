#' 调用TwoSampleMR包来查询多个暴露对一个结局的阳性结果
#' @param exposure 暴露的id号，可以是多个，例如：c("ieu-a-2","ieu-a-7")
#' @param outcome  结局的id号，有且只能是一个
#' @param write 逻辑值，将初筛结果导出到MR多对一初筛.csv文件中，默认：write=T 导出该文件
#' @param p1 筛选工具变量的显著性阈值，默认：p1 = 5e-08
#' @param clump 逻辑值，决定暴露中的工具变量是否进行去连锁不平衡，默认：clump = TRUE
#' @param p2 TwoSampleMR官网给出的解释是：Secondary clumping threshold，但经我查阅源代码发现这个参数压根没用到(⑅˃◡˂⑅)
#' @param r2 暴露中工具变量去连锁不平衡时的r2阈值，默认：r2= 0.001
#' @param kb 暴露中工具变量去连锁不平衡时的步长也叫窗口（windows）的阈值，默认：kb = 10000
#' @param LD 结局中工具变量的连锁不平衡阈值，即：若暴露中不存在某个工具变量，则寻找该工具变量的代理SNP，LD就是筛选阈值，只有大于该阈值才可作为代理SNP(⑅˃◡˂⑅)
#' @export
#' @return 多个暴露对一个结局的数据框
find_multiexposure_outcome=function(exposure=exposure,outcome=outcome,write=T,p1 = 5e-08,clump = TRUE,p2 = 5e-08,r2= 0.001,kb = 10000,LD=0.8){
  .find_multiexposure_outcome_fun(exposure,outcome,write,p1,clump,p2,r2,kb,LD)
}

if(T){
.find_multiexposure_outcome_fun=function (exposure,outcome,write,p1,clump,p2,r2,kb,LD) 
  UseMethod(".find_multiexposure_outcome_fun")
}

#查询多个暴露对一个结局的阳性结果
.find_multiexposure_outcome_fun.default=function(exposure,outcome,write,p1,clump,p2,r2,kb,LD){
  #设置MR_result_combine为空数据框  
  columns= c("id.exposure","id.outcome","conclusion")
  MR_result_combine = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(MR_result_combine) = columns

  #判断结局的个数
  if(length(outcome)>1){
    stop("结局id只能有一个")
  }

  start=Sys.time()
  #对每一个暴露进行循环分析
  for (id in 1:length(exposure)) {
    exposure_id=exposure[id]
    
    #计算还剩下几个暴露
    remain=length(exposure)-id
    print(paste0("当前分析的是第",id,"个暴露,","还剩下",remain,"个暴露"))

    #提取IV
    exposure_dat <- TwoSampleMR::extract_instruments(
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
    outcome_dat <- TwoSampleMR::extract_outcome_data(snps=exposure_dat$SNP, outcomes = outcome,rsq =LD ,access_token = NULL)
    
    #若在结局中无有效的IV，则需重新设置LD，然后在查找合适的工具变量
    if (is.null(outcome_dat)) {
      res<-cbind(exposure_id,outcome,"需要重新设定LD阈值，否则无法进行MR分析")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #协调和矫正数据
    dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)
    
    #执行MR分析
    res <- TwoSampleMR::mr(dat)
    
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
  
  #将MR_result_combine设置为数据框然后重新设置列名
  MR_result_combine=data.frame(MR_result_combine)
  colnames(MR_result_combine)<-c("暴露","结局","结果")
  if (write) {
    write.csv(MR_result_combine,"MR多对一初筛.csv",row.names = F)
  }
  
  end=Sys.time()
  print(end-start)
  return(MR_result_combine)
  
}



#' 调用TwoSampleMR包来查询一个暴露对多个结局的阳性结果
#' @param exposure 暴露的id号，可以是多个，例如：c("ieu-a-2","ieu-a-7")
#' @param outcome  结局的id号，有且只能是一个
#' @param write 逻辑值，将初筛结果导出到MR多对一初筛.csv文件中,默认：write=T 导出该文件
#' @param p1 筛选工具变量的显著性阈值，默认：p1 = 5e-08
#' @param clump 逻辑值，决定暴露中的工具变量是否进行去连锁不平衡，默认：clump = TRUE
#' @param p2 TwoSampleMR官网给出的解释是：Secondary clumping threshold，但经我查阅源代码发现这个参数压根没用到(⑅˃◡˂⑅)
#' @param r2 暴露中工具变量去连锁不平衡时的r2阈值，默认：r2= 0.001
#' @param kb 暴露中工具变量去连锁不平衡时的步长也叫窗口（windows）的阈值，默认：kb = 10000
#' @param LD 结局中工具变量的连锁不平衡阈值，即：若暴露中不存在某个工具变量，则寻找该工具变量的代理SNP，LD就是筛选阈值，只有大于该阈值才可作为代理SNP(⑅˃◡˂⑅)
#' @export
#' @return 一个暴露对多个结局的数据框
find_exposure_multioutcome=function(exposure=exposure,outcome=outcome,write=T,p1 = 5e-08,clump = TRUE,p2 = 5e-08,r2= 0.001,kb = 10000,LD=0.8){
  .find_exposure_multioutcome_fun(exposure,outcome,write,p1,clump,p2,r2,kb,LD)
}

if(T){
.find_exposure_multioutcome_fun=function (exposure,outcome,write,p1,clump,p2,r2,kb,LD) 
  UseMethod(".find_exposure_multioutcome_fun")
}

#查询一个暴露对多个结局的阳性结果
.find_exposure_multioutcome_fun.default=function (exposure,outcome,write,p1,clump,p2,r2,kb,LD){
  #设置MR_result_combine为空数据框  
  columns= c("暴露","结局","结果")
  MR_result_combine = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(MR_result_combine) = columns

  #判断暴露的个数
  if(length(exposure)>1){
    stop("暴露id只能有一个")
  }
  
  starttime=Sys.time()
  #对每一个暴露进行循环分析
  for (id in 1:length(outcome)) {
    outcome_id=outcome[id]
    
    #计算还剩下几个暴露
    remain=length(outcome)-id
    print(paste0("当前分析的是第",id,"个结局,","还剩下",remain,"个结局"))
    
    #提取IV
    exposure_dat <- TwoSampleMR::extract_instruments(
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
    outcome_dat <- TwoSampleMR::extract_outcome_data(snps=exposure_dat$SNP, outcomes = outcome_id,rsq =LD ,access_token = NULL)
    
    #若在结局中无有效的IV，则需重新设置LD，然后在查找合适的工具变量
    if (is.null(outcome_dat)) {
      res<-cbind(exposure,outcome_id,"需要重新设定LD阈值，否则无法进行MR分析")
      MR_result_combine<-rbind(MR_result_combine,res)
      next()
    }
    
    #协调和矫正数据
    dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)
    
    #执行MR分析
    res <- TwoSampleMR::mr(dat)
    
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
  colnames(MR_result_combine)<-c("暴露","结局","结果")
  if (write) {
    write.csv(MR_result_combine,"MR一对多初筛.csv",row.names = F)
  }
  
  endtime=Sys.time()
  print(endtime-starttime)
  
  return(MR_result_combine)
  
}


#' 调用ieugwasr包对工具变量中缺失的效应等位基因频率（eaf）进行填充，数据来源于1000基因组数据
#' @param exposure TwosampleMR包读入暴露的数据框
#' @param outcome TwosampleMR包读入结局的数据框
#' @param pop 参考1000基因组数据中的人种
#' @export
#' @return 返回TwosampleMR包中exposure或outcome格式的数据框；注意：有的SNP可能会出现eaf填充失败，这是因为1000基因组数据中缺乏该SNP的eaf信息(⑅˃◡˂⑅)
find_snp_add_eaf=function(exposure=NULL,outcome=NULL,pop="EUR"){
  .snp_add_eaf_fun(exposure,outcome,pop)
}

if(T){
.snp_add_eaf_fun=function(exposure,outcome,pop)
  UseMethod(".snp_add_eaf_fun")
}

#填充eaf
.snp_add_eaf_fun.default=function(exposure,outcome,pop){
  #如果exposure和outcome都有或都无数据则停止运行
  if(!is.null(exposure)&!is.null(outcome)){
    stop("只能输入暴露或者结局，二选其一")
  }else if (is.null(exposure)&is.null(outcome)) {
    stop("请输入暴露或者结局，二选其一")
  }
  #如果exposure不为空，则对exposure中的eaf数据进行插补
  if (!is.null(exposure)) {
    #判断exposure中SNP的eaf是否存在NA，若不存在则停止程序
    if(sum(is.na(exposure$eaf.exposure))==0){
      stop("exposure中SNP的eaf不存在NA，无需填充")
    }
    cat("###################################################################################\n")
    cat("##################################开始暴露数据eaf插补#################################\n")
    cat("###################################################################################\n")

    for(i in 1:nrow(exposure)){
      if (is.na(exposure$eaf.exposure[i])) {
        cat(paste0("#########################","当前填充的是第",i,"个SNP：",exposure$SNP[i],"#########################\n"))
        #如果该SNP的eaf缺失则填充SNP的eaf
        snp_info=as.data.frame(ieugwasr::afl2_rsid(exposure$SNP[i]))
        if (nrow(snp_info)==0) {
          message(paste0("第",i,"个SNP：",exposure$SNP[i],"填充失败"))
          next()
        }
        #根据pop参数选择人种
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
          stop("pop 输入错误,只能输入：c(EUR,AFR,AMR,EAS,SAS)")
        }
        if (is.null(eaf)) {
          message(paste0("第",i,"个SNP：",exposure$SNP[i],"填充失败"))
          next()
        }
        exposure$eaf.exposure[i]=eaf
      }
    }
    return(exposure)
  }

  #如果outcome不为空，则对outcome中的eaf数据进行插补
  if (!is.null(outcome)) {
    #判断outcome中SNP的eaf是否存在NA，若不存在则停止程序
    if(sum(is.na(outcome$eaf.outcome))==0){
      stop("outcome中SNP的eaf不存在NA，无需填充")
    }
    cat("###################################################################################\n")
    cat("##################################开始结局数据eaf插补#################################\n")
    cat("###################################################################################\n")
    for(i in 1:nrow(outcome)){
      if (is.na(outcome$eaf.outcome[i])) {
        cat(paste0("#########################","当前填充的是第",i,"个SNP：",outcome$SNP[i],"#########################\n"))
        #如果该SNP的eaf缺失则填充SNP的eaf
        snp_info=as.data.frame(ieugwasr::afl2_rsid(outcome$SNP[i]))
        if (nrow(snp_info)==0) {
          message(paste0("第",i,"个SNP：",outcome$SNP[i],"填充失败"))
          next()
        }
        #根据pop参数选择人种
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
          stop("pop 输入错误,只能输入：(EUR,AFR,AMR,EAS,SAS)中的一种")
        }
        if (is.null(eaf)) {
          message(paste0("第",i,"个SNP：",outcome$SNP[i],"填充失败"))
          next()
        }
        outcome$eaf.outcome[i]=eaf
      }
    }
    return(outcome)
  }
}


#' 调用phenoscanner包，查找工具变量相关的混杂因素 
#' @param confound 设定的混杂因素，可以是一个或者多个，例如：c("Body mass index","Coronary heart disease")
#' @param exposure_dat TwosampleMR包中函数extract_instruments或者harmonise_data返回的数据框
#' @param query_gene  设置phenoscanner函数搜索的基因，绝大部分情况下用不到该参数，默认：NULL
#' @param query_region  设置phenoscanner函数搜索的基因组范围，绝大部分情况下用不到该参数，默认：NULL
#' @param catalogue 设置phenoscanner函数搜索的目录（None, GWAS, eQTL, pQTl, mQTL, methQTL），默认：catalogue = GWAS
#' @param pvalue  设置工具变量与混杂因素的筛选阈值，，默认：pvalue = 5e-08
#' @param proxies 设置phenoscanner函数搜索的代理人种（None, AFR, AMR, EAS, EUR, SAS），默认：proxies = None
#' @param r2  设置工具变量与其他可能与混杂因素相关连SNP的连锁不平衡阈值，大于该阈值的代理SNP得以保留，默认：r2 = 0.8
#' @param build 设置基因组的版本（37，38），默认：build = 37
#' @param write 逻辑值,将工具变量相关的混杂因素导出，默认：write=T 导出结果文件
#' @param save_path 设置结果的保存路径和名称
#' 
#' @export
#' @return 返回TwosampleMR包中exposure格式的数据框，若PhenoScanner网站中有筛选到工具变量SNP的混杂因素，则将该SNP的mr_keep.exposure/mr_keep状态设置为FALSE(⑅˃◡˂⑅)
deletion_confounding_snp = function(confound = NULL,
                                    exposure_dat = NULL,
                                    query_gene = NULL,
                                    query_region = NULL,
                                    catalogue = "GWAS",
                                    pvalue = 5e-08,
                                    proxies = "None",
                                    r2 = 0.8,
                                    build = 37,
                                    write=TRUE,
                                    save_path="MR_ivs_confounding.csv"){
  .deletion_confounding_snp_fun(confound,exposure_dat,query_gene,query_region,catalogue,pvalue,proxies,r2,build,write,save_path)
}

if(T){
.deletion_confounding_snp_fun=function(confound,exposure_dat,query_gene,query_region,catalogue,pvalue,proxies,r2,build,write,save_path)
  UseMethod(".deletion_confounding_snp_fun")
}

#查找工具变量相关的混杂因素
.deletion_confounding_snp_fun.default = function(confound,exposure_dat,query_gene,query_region,catalogue,pvalue,proxies,r2,build,write,save_path) {
  #将多个混杂因素进行合并，以|为分隔符
  confound = paste(confound, collapse = "|")

  #设置MR_result_combine为空数据框  
  columns= c("工具变量","混杂因素","P值")
  MR_result_combine = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(MR_result_combine) = columns

  #循环进行phenoscanner函数查找与query_snp相关的GWAS结果
  for (i in 1:length(exposure_dat$SNP)) {
    #获取SNP的rsID号
    query_snp = exposure_dat$SNP[i]
    #显示运行程度
    cat("#################################################################################\n")
    num = length(exposure_dat$SNP) - i
    print = paste0("查询第", i, "个工具变量的混杂因素", "，还剩下", num, "个工具变量")
    cat(paste0("#####################", print, "#####################"))
    cat("\n#################################################################################\n")
    
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
  #将MR_result_combine设置为数据框格式并重新命名
  MR_result_combine=data.frame(MR_result_combine)
  colnames(MR_result_combine) = columns

  if (write) {
    write.csv(MR_result_combine,file = save_path,row.names = F)
  }  
  #返回暴露/矫正后的数据
  return(exposure_dat)
}


#' 调用ieugwasr包，在本地电脑进行去连锁不平衡，此功能可用于工具变量较多，在线进行clump时总是报错的情况
#' @param dat TwosampleMR包读入暴露数据的后的数据框
#' @param clump_kb  暴露中工具变量去连锁不平衡时的步长也叫窗口（windows）的阈值，默认：kb = 10000
#' @param clump_r2  暴露中工具变量去连锁不平衡时的r2阈值，默认：r2 = 0.001
#' @param clump_p1  暴露中工具变量去连锁不平衡时的显著性阈值, 默认：clump_p1 = 1
#' @param bfile 1000基因组数据中作为参考的plink格式人种数据，下载地址：http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz，默认：欧洲人群
#' @export
#' @return 返回TwosampleMR包中exposure格式的数据框
clump_data_local=function(dat=exposure_dat,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 1,bfile = "D:/生物信息学习/生信自研课程/孟德尔随机化/data/1kg.v3/EUR"){
  .clump_data_local_fun(dat,clump_kb,clump_r2,clump_p1,bfile)
}

if(T){
.clump_data_local_fun=function(dat,clump_kb,clump_r2,clump_p1,bfile)
  UseMethod(".clump_data_local_fun")
}

#定义函数在本地进行去连锁不平衡
.clump_data_local_fun=function(dat,clump_kb,clump_r2,clump_p1,bfile,plink=plinkbinr::get_plink_exe()){
  pval_column <- "pval.exposure"
  if (!is.data.frame(dat)) {
    stop("输入数据需得是从TwosampleMR包返回的数据框")
  }
  if ("pval.exposure" %in% names(dat) & "pval.outcome" %in% 
      names(dat)) {
    message("存在pval.exposure和pval.outcome两列。使用pval.exposure进行去连锁不平衡")
  }
  else if (!"pval.exposure" %in% names(dat) & "pval.outcome" %in% 
           names(dat)) {
    message("pval.exposure列不存在，使用pval.outcome列进行去连锁不平衡")
    pval_column <- "pval.outcome"
  }
  else if (!"pval.exposure" %in% names(dat)) {
    message("pval.exposure列不存在, 所有工具变量的clumping p-value设置为0.99")
    dat$pval.exposure <- 0.99
  }
  else {
    pval_column <- "pval.exposure"
  }
  if (!"id.exposure" %in% names(dat)) {
    dat$id.exposure <- random_string(1)
  }
  d <- data.frame(rsid = dat$SNP, pval = dat[[pval_column]], 
                  id = dat$id.exposure)
  
  out <- ieugwasr::ld_clump_local(dat = d, clump_kb = clump_kb, clump_r2 = clump_r2, 
                            clump_p = clump_p1, bfile = bfile, plink=plink)
  keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, 
                                                     out$id)
  return(dat[keep, ])
}



#' 调用ukbb数据库的参考数据variants.tsv.gz对ukbb中的GWAS数据批量化进行SNP位点注释
#' @param ukbb_gwas_dir ukbb格式的GWAS数据所在文件夹
#' @param save_dir  设置输出数据的文件夹
#' @export
#' @return 返回一个或多个SNP位点注释后的ukbb的GWAS数据
batch_ukbb_chrpos2rsID=function(ukbb_gwas_dir=NULL,save_dir=NULL){
  .batch_ukbb_chrpos2rsID_fun(ukbb_gwas_dir,save_dir)
}

if(T){
.batch_ukbb_chrpos2rsID_fun=function(ukbb_gwas_dir,save_dir)
  UseMethod(".batch_ukbb_chrpos2rsID_fun")
}

.batch_ukbb_chrpos2rsID_fun=function(ukbb_gwas_dir,save_dir){
  #记录开始时间
  start_time = Sys.time()

  #导入参考数据集
  if(!"ukbb_ref" %in% ls(envir = .GlobalEnv, all.names = T)){
    #rm(ukbb_ref, envir = .GlobalEnv)
    file_ukbb_ref_Rdata=system.file("data", # 子目录
                                    "ukbb_ref.Rdata", # 文件名
                                    package="friendly2MR", 
                                    mustWork=TRUE)
    #加载ukbb的SNP参考数据集
    load(file = file_ukbb_ref_Rdata)
  }

  #列出所有的ukbb数据集
  files<-list.files(ukbb_gwas_dir)
  #提取.bgz结尾的文件
  files<-files[grepl("\\.bgz$",files)]
  #显示所有.bgz结尾的文件
  print(files)
  
  #修改文件扩展名
  for (f in files){
  newname<-sub(".bgz",'.gz',f)
  file.rename(paste0(ukbb_gwas_dir,f),paste0(ukbb_gwas_dir,newname))
  }

  #列出所有的以.gz为结尾的ukbb数据集
  ukbb_gwas_dir_file<-list.files(ukbb_gwas_dir)
  all_ukbb_gwas_file<-ukbb_gwas_dir_file[grepl("\\.gz$",ukbb_gwas_dir_file)]

  #循环处理数据
  for(i in 1:length(all_ukbb_gwas_file)){

    #显示ukbb数据的处理数量
    cat("################################正在进行UKBB的chrpos2rsID处理#####################\n")
    cat(paste0("#############当前处理的是第",i,"个数据集：",all_ukbb_gwas_file[i],"###########\n"))
    cat(paste0("##############################","还剩下：",length(all_ukbb_gwas_file)-i,"个数据集##################################\n"))

    #设置读入文件
    gwas_ukb=paste0(ukbb_gwas_dir,"/",all_ukbb_gwas_file[i])
    #读入数据
    ukbb = as.data.frame(data.table::fread(gwas_ukb))

    #合并数据
    merge=dplyr::left_join(ukbb,ukbb_ref,by=c("variant"="variant"))
    #将variant列拆分为"CHR","BP","other_allele","effect_allele"
    merge=tidyr::separate(merge, variant, c("CHR","BP","other_allele","effect_allele"),sep=":")

    #计算效应等位基因频率
    merge$eaf<-merge[,"AC"]/(merge[,"n_complete_samples"]*2)

    #提取合适的列,c("CHR","BP","rsid","effect_allele","other_allele","eaf","beta","se","pval","n_complete_samples")
    merge<-merge[,c("CHR","BP","rsid","effect_allele","other_allele","eaf","beta","se","pval","n_complete_samples")]
    #列名重命名
    colnames(merge)<-c("CHR","BP","SNP","effect_allele","other_allele","eaf","beta","se","pval","N")

    #去除SNP列中非rsID号的位点
    merge=merge[grep("rs",merge$SNP),]%>%na.omit()
    head(merge)

    single_gwas_name=stringr::str_sub(all_ukbb_gwas_file[i], start = 1L, end = -7L)
    #保存数据
    if (!is.null(save_dir)) {
      data.table::fwrite(
        as.data.frame(merge),
        file = paste0(save_dir,"/",single_gwas_name,"txt"),
        sep = "\t",
        row.names = F
      )
    }else if (is.null(save_dir)){
      data.table::fwrite(
      as.data.frame(merge),
      file = paste0(ukbb_gwas_dir,"/",single_gwas_name,"txt"),
      sep = "\t",
      row.names = F
    )}else{stop("save_dir参数设置的文件保存路径不正确，请重新输入！！")}

  }

  end_time = Sys.time()
  print(end_time-start_time)
  #清空merge,ukbb变量
  rm(merge,ukbb)
}


#' 调用epigraphdb数据库的数据,对多个暴露和/或多个结局进行快速筛选
#' @param exposure  暴露的id号或名称，可以是多个，例如：c("ieu-a-2","ieu-a-7")，c("Body mass index","Coronary heart disease")，如果exposure = NUL，则纳入全部的数据
#' @param outcome 结局的id号或名称，可以是多个，例如：c("ieu-a-2","ieu-a-7")，c("Body mass index","Coronary heart disease")，如果outcome = NUL，则纳入全部的数据
#' @param pval_threshold 暴露和结局的孟德尔随机化阈值，默认：pval_threshold = 1e-05
#' @param write 逻辑值，是否将输出结果保存到文件中，默认：write = T 导出结果文件
#' @param save_path 筛查结果的保存路径和文件名称，默认：save_path = "find_multiexposure_multioutcome_epigraphdb.csv"
#' @export
#' @return  epigraphdb数据库中快速筛查多个暴露和/或多个结局的筛选结果
find_multiexposure_multioutcome_epigraphdb=function(exposure = NULL,outcome = NULL,pval_threshold = 1e-05,write=T,save_path="find_multiexposure_multioutcome_epigraphdb.csv"){
  .find_multiexposure_multioutcome_epigraphdb_fun(exposure,outcome,pal_threshold,write,save_path)
}

if(T){
.find_multiexposure_multioutcome_epigraphdb_fun=function(exposure,outcome,pal_threshold,write,save_path)
  UseMethod(".find_multiexposure_multioutcome_epigraphdb_fun")
}

# 筛选多个暴露和多个结局中的显著结果 -------------------------------------------------------
.find_multiexposure_multioutcome_epigraphdb_fun=function(exposure,outcome,pal_threshold,write,save_path){
  start=Sys.time()
  cat("##########################开始筛查######################################\n")
  #设置mr_fast_result_all
  mr_fast_result_all<-data.frame()

  #设置ieu的id号保存路径
  ao_file_RData=system.file("data", # 子目录
                            "IEU_ao.RData", # 文件名
                            package="friendly2MR", 
                            mustWork=TRUE)
  #首先生成IEU的性状信息：
  if(!file.exists(ao_file_RData)){
    ao <- TwoSampleMR::available_outcomes()
    save(ao,file = ao_file_RData)
  }
  load(file = ao_file_RData)
  
  exposure_length=length(exposure)
  outcome_length=length(outcome)
  
  for (i in 1:exposure_length) {
    for (j in 1:outcome_length) {
      exposure_id = exposure[i]
      outcome_id = outcome[j]
      # #判断暴露id和性状是否在数据库中
      # if (sum((exposure_id %in% ao$trait),(exposure_id %in% ao$id))>0) {
      # }else{stop("\n##########################该暴露不在数据库中######################################\n")}
      # # #判断结局id和性状是否在数据库中
      # if (sum((outcome_id %in% ao$trait),(outcome_id %in% ao$id))>0) {
      # }else{stop("\n##########################该结局不在数据库中######################################\n")}
      #首先判断是ieu的id还是性状名称：
      if (!is.null(exposure_id)) {
        if (exposure_id %in% ao$id) {
          exposure_trait <- ao[ao$id == exposure_id, ]$trait
        } else{
          exposure_trait = exposure_id
        }
      } else{
        exposure_trait = exposure_id
      }
      
      if (!is.null(outcome)) {
        if (outcome_id %in% ao$id) {
          outcome_trait <- ao[ao$id == outcome_id, ]$trait
        } else{
          outcome_trait = outcome_id
        }
      } else{
        outcome_trait = outcome_id
      }
      
      mr_fast_result = epigraphdb::mr(
        exposure_trait = exposure_trait,
        pval_threshold = pval_threshold,
        outcome_trait = outcome_trait,
        mode = "table"
      )
      
      #根据ieu号进行筛选：
      # mr_fast_result1<-data.frame()
      if (!is.null(exposure_id)) {
        if (exposure_id %in% ao$id) {
          mr_fast_result <-
            mr_fast_result[mr_fast_result$exposure.id == exposure_id, ]
        }else{mr_fast_result=mr_fast_result}
      }

      if (!is.null(outcome_id)) {
        if (outcome_id %in% ao$id) {
          mr_fast_result <- mr_fast_result[mr_fast_result$outcome.id == outcome_id, ]
        }else{mr_fast_result=mr_fast_result}
      }
      # 多个结果合并
      mr_fast_result_all<-unique(rbind(mr_fast_result,mr_fast_result_all))
      
    } 
  }
  #判断是否保存  
  if (write) {
    write.csv(mr_fast_result_all,file=save_path,row.names = F)
  }
  cat("##########################筛查结束######################################\n")
  end=Sys.time()
  print(end-start)
  return(mr_fast_result_all)
}




#' 将从IEU官网下载到本地的gwasvcf文件转换为TwosampleMR的暴露或结局格式
#' @param dat 从IEU官网下载到本地的gwasvcf文件，例如：IEU-a-2.vcf.gz
#' @param type  指定该数据作为暴露还是结局，例如：exposure，outcomes，默认：type = "exposure"
#' @param p1  筛选工具变量的显著性阈值，默认：p1 = 5e-08
#' @param clump 逻辑值，决定暴露中的工具变量是否进行去连锁不平衡，默认：clump = TRUE
#' @param clump_kb 暴露中工具变量去连锁不平衡时的步长也叫窗口（windows）的阈值，默认：kb = 10000
#' @param clump_r2 暴露中工具变量去连锁不平衡时的r2阈值，默认：clump_r2 = 0.001
#' @param pop 去连锁不平衡时用到的1000人基因组人种，默认：pop = "EUR"
#' @param write 逻辑值，是否将输出结果保存到文件中，默认：write = T 导出结果文件
#' @param save_path 筛查结果的保存路径和文件名称，默认：save_path = "TwosampleMR_formate.txt"
#' @export
#' @return  epigraphdb数据库中快速筛查多个暴露和/或多个结局的筛选结果
gwasvcf2TwosampleMR_local=function(dat=dat,type = "exposure",p1=5e-08,clump=T,clump_kb = 10000,clump_r2 = 0.001,pop = "EUR",write=T,save_path ="TwosampleMR_formate.txt"){
  .gwasvcf2TwosampleMR_local_fun(dat,type,p1,clump,clump_kb,clump_r2,pop,write,save_path)
}

if(T){
.gwasvcf2TwosampleMR_local_fun=function(dat,type,p1,clump,clump_kb,clump_r2,pop,write,save_path)
  UseMethod(".gwasvcf2TwosampleMR_local_fun")
}

.gwasvcf2TwosampleMR_local_fun<-function(dat,type,p1,clump,clump_kb,clump_r2,pop,write,save_path){
  #1.判断是否存在dat文件
  #2.若dat是文件路径，则直接将本地的gwasvcf文件转换为TwosampleMR格式
  start=Sys.time()
  if(file.exists(dat)){
    gwasvcf_file=VariantAnnotation::readVcf(dat)
    my_gwas_data=gwasglue::gwasvcf_to_TwoSampleMR(gwasvcf_file,type = type)
  }

  if (is.null(p1)) {
    stop("###################请输入合适的p1以筛选显著的SNP作为工具变量###################\n")
  }
  if(type == "exposure"){
    my_gwas_data<-subset(my_gwas_data,pval.exposure<p1)
  }
  else if (type == "outcomes") {
    my_gwas_data<-subset(my_gwas_data,pval.outcomes<p1)
  }
  else{stop("###############################输入的type不正确###############################\n")}
  #判断是否进行去来连锁不平衡
  if (clump) {
     my_gwas_data <- TwoSampleMR::clump_data(my_gwas_data,clump_kb=clump_kb,clump_r2=clump_r2,pop=pop)
  }
  
  #判断是否需要将文件导出
  if (write) {
    write.table(my_gwas_data,file =save_path, quote = FALSE,sep = "\t",row.names = F)
  }

  end=Sys.time()
  print(end-start)
  cat("###################################转换完成###################################\n")
  return(my_gwas_data)
}



