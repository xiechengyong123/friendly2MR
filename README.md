# friendly2MR

## **写给观众老爷们的话：**

    首先，这是一个非常简约（俗称：渣渣）的R包，主要包含了孟德尔随机化中的一些简单辅助函数，目前实现了以下功能：

    1.查询多个暴露对一个结局的阳性结果：find_anyexposur_outcome()：一个一个跑，有点慢，凑合着跑吧！

    2.查询一个暴露对多个结局的阳性结果:find_exposur_anyoutcome()：同样，一个一个跑，有点慢，凑合着跑吧！

    3.填充工具变量中插补SNP的eaf：主要是调用1000基因组数据做的填充，**与原暴露或结局中GWAS中SNP的eaf或有所不同，谨慎使用！**

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
* **模块：模块：双样本多变量的MR分析（一键生成）**
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

微信：

![1679155331649](https://github.com/xiechengyong123/friendly2MR/blob/main/%E5%BE%AE%E4%BF%A1%E5%9B%BE%E7%89%87_20230319000148.jpg)
