# friendly2MR

## **写给观众老爷们的话：**

    首先，这是一个非常简约（俗称：渣渣）的R包，主要包含了孟德尔随机化中的一些简单辅助函数，目前实现了以下功能：

    1.查询多个暴露对一个结局的阳性结果：find_anyexposur_outcome()：一个一个跑，有点慢，凑合着跑吧！

    2.查询一个暴露对多个结局的阳性结果:find_exposur_anyoutcome()：同样，一个一个跑，有点慢，凑合着跑吧！

    3.填充工具变量中插补SNP的eaf：主要是调用1000基因组数据做的填充，与原暴露或结局中GWAS中SNP的eaf或有所不同，谨慎使用！

    4.查找工具变量相关的混杂因素：这个主要是调用了[phenoscanner](http://www.phenoscanner.medschl.cam.ac.uk/)函数，但存在SNP在该数据库中找不到的情况，此种情况下，将默认无混杂因素

## 安装和加载

```
if (!requireNamespace("remotes", quietly = TRUE))install.packages("remotes")
remotes::install_github("xiechengyong123/friendly2MR")
library(friendly2MR)
```


## 使用链接

[b站视屏链接](https://www.bilibili.com/video/BV1Lk4y1h7j8/?spm_id_from=333.999.0.0&vd_source=559aa6843f51710f9b5e95a85661a0f3)

## 后续计划

* **模块：特色功能**
* **模块：双样本单变量的MR分析（一键生成）**
* **模块：双样本多变量的MR分析（一键生成）**
* **模块：双向因果效应（一键生成）**
* **模块：中介MR（一键生成）**
* **模块：共定位**
* **模块：计算PRS**
* **模块：非线性孟德尔随机化**
* **模块：GSMR分析功能**
* **模块：药物靶点的孟德尔分析（SMR）**
* **模块：SMR结合QTL分析**
* **模块：添加TWAS（FUSION）分析功能**

## 与我联系

**微信：**

![1679155331649](https://github.com/xiechengyong123/friendly2MR/blob/main/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20230319000148.jpg)

## 更新日志:

1.添加了四个函数：分别是batch_ukbb_chrpos2rsID，clump_data_local，find_multiexposure_multioutcome_epigraphdb，gwasvcf2TwosampleMR_local；

    batch_ukbb_chrpos2rsID：调用ukbb数据库的参考数据variants.tsv.gz对ukbb中的GWAS数据批量化进行SNP位点注释；

    clump_data_local：调用ieugwasr包，在本地电脑进行去连锁不平衡，此功能可用于工具变量较多，在线进行clump时总是报错的情况；

    find_multiexposure_multioutcome_epigraphdb：调用epigraphdb数据库的数据,对多个暴露和/或多个结局进行快速筛选；

    gwasvcf2TwosampleMR_local：将从IEU官网下载到本地的gwasvcf文件转换为TwosampleMR的暴露或结局格式；

2.修改函数名称：原包中的find_anyexposur_outcome和find_exposur_anyoutcome函数名称进行更改：

    find_anyexposur_outcome-->find_multiexposure_outcome;

    find_exposur_anyoutcome-->find_exposure_multioutcome;

3.添加了函数和参数的注释，便于观众老爷们的理解。

4.文件中包含了ukbb的参考文件，太大了无法上传，只能存在我的百度云盘，链接：https://pan.baidu.com/s/1fdQu2rluE2fd6sgXcHIX0g
提取码：4sc4。

5.**该R包只分享10个人且30天内有效，先到先得，过时无效**。

----------------------------------------------------------2023年-4月-9日-----------------------------------------------------
