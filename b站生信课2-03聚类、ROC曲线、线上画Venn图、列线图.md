<a id="mulu">目录</a>
<a href="#mulu" class="back">回到目录</a>
<style>
    .back{width:40px;height:40px;display:inline-block;line-height:20px;font-size:20px;background-color:lightyellow;position: fixed;bottom:50px;right:50px;z-index:999;border:2px solid pink;opacity:0.3;transition:all 0.3s;color:green;}
    .back:hover{color:red;opacity:1}
    img{vertical-align:bottom;}
</style>

<!-- @import "[TOC]" {cmd="toc" depthFrom=3 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [一致性聚类与无监督聚类](#一致性聚类与无监督聚类)
- [ROC曲线](#roc曲线)
- [Venn图绘制](#venn图绘制)
- [列线图](#列线图)

<!-- /code_chunk_output -->

<!-- 打开侧边预览：f1->Markdown Preview Enhanced: open...
只有打开侧边预览时保存才自动更新目录 -->

写在前面：本篇教程来自b站课程[TCGA及GEO数据挖掘入门必看](https://www.bilibili.com/video/BV1b34y1g7RM) P18-P21

### 一致性聚类与无监督聚类
需要数据：多因素cox回归结果（其实也可以是lasso/单因素cox回归结果）、tpm表达矩阵
需要包：`ConsensusClusterPlus`
``` r
if(!require("ConsensusClusterPlus", quietly = T))
{
  library("BiocManager");
  BiocManager::install("ConsensusClusterPlus");
  library("ConsensusClusterPlus");
}
```
**读取数据**：选出肿瘤组的样本，取出筛选后基因的在各样本中的表达量
``` r
# tpm表达矩阵
data <- read.table("C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\TCGA_LUSC_TPM.txt", check.names = F, row.names = 1, sep = '\t', header = T);
dimnames <- list(rownames(data), colnames(data));
data <- matrix(as.numeric(as.matrix(data)), nrow = nrow(data), dimnames = dimnames);
# 选出肿瘤组的样本
group <- sapply(strsplit(colnames(data), '\\-'), "[", 4);
group <- sapply(strsplit(group, ''), "[", 1);
data <- data[, group==0];
# 筛选基因
multi_cox_gene <- read.table("C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\multiCox.txt", check.names = F, row.names = 1, sep = '\t', header = T);
data <- data[rownames(multi_cox_gene), ];
```
![一致性聚类与无监督聚类1](./md-image/一致性聚类与无监督聚类1.png){:width=250 height=250}
**对样品进行聚类分型**，使用`ConsensusClusterPlus`函数：
- `maxK`最大的K值，形成一系列梯度
- `pItem`选择百分之多少的样本重复抽样
- `pfeature`选择百分之多少的基因重复抽样
- `reps`重复抽样的数目，可以先设置为100，结果不错再设置为1000（这样结果更严谨）
- `clusterAlg`聚类算法，取值："hc"/"pam"/"km"
- `distanc`距离矩阵算法，取值："pearson"/"spearman"/"euclidean"
- `title`输出结果的文件夹名字，包含输出的图片等
- `seed`随机种子，用于固定结果
- `plot`输出图片的格式
``` r
res <- ConsensusClusterPlus(
  data,
  maxK = 9,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  title = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\ConsensusClusterPlus",
  clusterAlg = "pam",
  distance = "euclidean",
  seed = 123,
  plot = "png"
);
```
**结果分析**：
- 第一张图标识相关度与颜色的关系（图例）：1是非常相关（蓝色），0是不相关（白色）
  ![一致性聚类与无监督聚类2](./md-image/一致性聚类与无监督聚类2.png){:width=300 height=300}
- 002-009：每个k值（分成了多少组）对应的聚类结果
  ![一致性聚类与无监督聚类3](./md-image/一致性聚类与无监督聚类3.png){:width=400 height=400}
  红框部分的不同颜色代表不同的组，它下面的2*2个方块代表每组的差异
  评判标准（以k=2为例）：
  - 组内的差异小（右上和左下两个块足够蓝）
  - 组间的差异大（右下和左上两个块足够白）
  - 每组的样本数不能过小（不能小于总样本的10%），可以通过红框中颜色占比看出
  
  可以看到k=2的图是符合标准的
- 010CDF值：
  ![一致性聚类与无监督聚类4](./md-image/一致性聚类与无监督聚类4.png){:width=400 height=400}
  曲线在x=0.1~0.9的变化越小的越好。可以看出k=2（红色线）符合标准
- 011CDF值变化：
  ![一致性聚类与无监督聚类5](./md-image/一致性聚类与无监督聚类5.png){:width=400 height=400}
  一般选取折线拐点（曲线变化趋势改变最大的点）附近的点，该图中拐点为x=3

综合上面的分析，我们选取k=2的结果
**根据上面的分组，对样本进行分型**：
``` r
clu_num <- 2;  # 分成几组（k值）
clu <- res[[clu_num]][["consensusClass"]];  # 聚类结果（分组信息）
clu <- as.data.frame(clu);
colnames(clu) <- c("cluster");
letter <- LETTERS[1:10];  # 每组的名称，这里是ABCD大写字母
uniq_clu <- levels(factor(clu$cluster));  # 原来每组的名称
clu$cluster <- letter[match(clu$cluster, uniq_clu)];  # 将每组名称改成我们刚才定义的大写字母
clu_save <- rbind(ID = colnames(clu), clu);
write.table(clu_save, file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\cluster.txt", sep = '\t', row.names = F, quote = F);
```
![一致性聚类与无监督聚类6](./md-image/一致性聚类与无监督聚类6.png){:width=300 height=300}
可以看到样品被分为了AB两组
### ROC曲线
![ROC曲线1](./md-image/ROC曲线1.png){:width=300 height=300}
横坐标1-Specificity是特异性（假阳性概率），纵坐标Sentivity是敏感性（真阳性概率）
AUC指曲线下的面积
- AUC=1：是完美的分类器
- 0.5< AUC <1：优于随机猜测，数值越大越好
- AUC=0.5：等同于随机猜测，没有预测价值
- AUC<0.5：比随机猜测差。但如果是反向猜测，该模型也可能优于随机猜测

``` r
if(!require("pROC", quietly = T))
{
  library("BiocManager");
  BiocManager::install("pROC");
  library(pROC);
}
library(survival);
library(survminer);
if(!require("timeROC", quietly = T))
{
  library("BiocManager");
  install.packages('listenv');
  install.packages('parallelly');
  BiocManager::install("timeROC");
  library(timeROC);
}
```
**读取并处理生存信息，与风险得分合并**：
``` r
# 读取生存信息
library("readxl");
library(tidyverse);
cli <- read_excel("C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\clinical.xlsx");
cli <- cli[, c("survival_time", "vital_status", "days_to_birth", "gender", "T", "N", "M", "stage_event", "anatomic_neoplasm_subdivision", "bcr_patient_barcode")];
# 处理生存信息
cli <- column_to_rownames(cli, "bcr_patient_barcode");  # 更改行名为样本名
cli$time <- as.numeric(cli$survival_time)/365;  # 存活天数用年表示
cli$state <- ifelse(cli$vital_status=='Alive', 0, 1);  # 0表示存活，1表示死亡
cli$Age <- round(as.numeric(cli$days_to_birth)/(-365));  # 年龄用年表示
cli$Gender <- ifelse(cli$gender=="MALE", 0, 1);  # 性别用01表示
cli$`T` <- substr(cli$`T`, 1, 1);
cli$`N` <- substr(cli$`N`, 1, 1);
cli$`M` <- substr(cli$`M`, 1, 1);  # TNM列只取第一个字符
cli$`T` <- gsub('X', NA, cli$`T`);
cli$`N` <- gsub('X', NA, cli$`N`);
cli$`M` <- gsub('X', NA, cli$`M`);  # TNM列将X替换为NA
cli$stage_event <- ifelse(grepl('X', cli$stage_event), NA, cli$stage_event);  # X变NA
cli$stage_event <- ifelse(grepl('IV', cli$stage_event), "4", cli$stage_event);  # IV变4
cli$stage_event <- ifelse(grepl('III', cli$stage_event), "3", cli$stage_event);  # III变3
cli$stage_event <- ifelse(grepl('II', cli$stage_event), "2", cli$stage_event);  # II变2
cli$stage_event <- ifelse(grepl('I', cli$stage_event), "1", cli$stage_event);  # I变1
cli$Stage <- as.numeric(cli$stage_event);
cli$`T` <- as.numeric(cli$`T`);
cli$`N` <- as.numeric(cli$`N`);
cli$`M` <- as.numeric(cli$`M`);  # TNM、Stage列转为数值型
cli$subdivision <- ifelse(  # 开头L->1  R->2
  startsWith(cli$anatomic_neoplasm_subdivision, "L"),
  1,
  ifelse(
    startsWith(cli$anatomic_neoplasm_subdivision, "R"),
    2,
    NA
  )
);
# 读取风险得分
risk <- read.table("C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\risk.txt", check.names = F, row.names = 1, sep = '\t', header = T);
# 合并
same_sample <- intersect(row.names(risk), row.names(cli));
risk <- risk[same_sample, ];
cli <- cli[same_sample, ];
rt <- cbind(
  cli[, c("time", "state")],
  riskScore = risk[, c("riskScore")],
  cli[, c("Age", "Gender", "T", "N", "M", "Stage", "subdivision")]
);
# 保存生存信息
library(writexl);
write_xlsx(cbind(ID = row.names(cli), cli[, c("time", "state", "Age", "Gender", "T", "N", "M", "Stage", "subdivision")]), "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\new_clinical.xlsx");
```
![ROC曲线2](./md-image/ROC曲线2.png){:width=180 height=180}
**ROC分析并绘图**：
``` r
# 可以是state/是否为肿瘤组~基因表达量/风险得分
roc1 <- roc(rt$state ~ rt$riskScore);  
pdf(file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\ROC.riskscore.pdf", width = 5, height = 5);
bioCol = c("DarkOrchid", "Orange2", "MediumSeaGreen", "NavyBlue", "#8B668B", "#FF4500", "#135612", "#561214");
plot(
  roc1,
  print.auc = T,
  col = bioCol,
  legacy.axes = T
);
dev.off();
```
![ROC曲线3](./md-image/ROC曲线3.png){:width=300 height=300}
另一种ROC分析--**`timeROC`时间依赖型生存曲线**：
- `T`事件时间
- `delta`事件状态（删失数据值为0）
- `marker`一个标记值，值越大，事件越可能发生，此处使用风险得分。如果使用的数据值越小越可能发生，则可以在前面加负号
- `other_markers`协变量（矩阵形式输入）
- `cause`所关注的时间结局，一般为1（死亡）
- `weighting`计算方法，默认使用KM模型，还可以是"cox"cox模型、"aalen"additive Aalen模型
- `times`想计算ROC曲线的时间节点
- `ROC`是否保存sensitivities的specificties值（默认为T）
- `iid`是否保存置信区间（默认为F，样本量大时很耗时间）

``` r
roc_rt <- timeROC(
  T = risk$time,
  delta = risk$state,
  marker = risk$riskScore,
  cause = 1,
  weighting = 'aalen',
  times = c(1, 3, 5),
  ROC = T
);
pdf(file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\ROC.all.pdf", width = 5, height = 5);
plot(roc_rt, time = 1, col = bioCol[1], title = F, lwd = 4);
plot(roc_rt, time = 3, col = bioCol[2], title = F, lwd = 4, add = T);  # 在前一条线上继续添加
plot(roc_rt, time = 5, col = bioCol[3], title = F, lwd = 4, add = T);  # 在前一条线上继续添加
legend(
  'bottomright',
  c(paste0('AUC at 1 year: ', sprintf("%.03f", roc_rt$AUC[1])),
    paste0('AUC at 3 years: ', sprintf("%.03f", roc_rt$AUC[2])),
    paste0('AUC at 5 years: ', sprintf("%.03f", roc_rt$AUC[3]))),
  col = bioCol[1:3],
  lwd = 4,
  bty = 'n',
  title = "All set"
);
dev.off();
```
![ROC曲线4](./md-image/ROC曲线4.png){:width=300 height=300}
**其它临床特征的ROC曲线**：
``` r
pre_time <- 5;  # 预测年限
# 先使用风险得分作roc预测
roc_rt <- timeROC(
  T = risk$time,
  delta = risk$state,
  marker = risk$riskScore,
  cause = 1,
  weighting = 'aalen',
  times = c(pre_time),
  ROC = T
);
pdf(file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\cliROC.all.pdf", width = 5.5, height = 5.5);
# 风险得分的roc曲线
plot(roc_rt, time = pre_time, col = bioCol[1], title = F, lwd = 4);
# 
abline(0, 1);
auc_text <- c(paste0("Risk", ", AUC=", sprintf("%.3f", roc_rt$AUC[2])));
# 再使用临床数据的其它列作roc预测
for (i in 4:ncol(rt)) {
  roc_rt <- timeROC(
    T = rt$time,
    delta = rt$state,
    marker = rt[, i],
    cause = 1,
    weighting = 'aalen',
    times = c(pre_time),
    ROC = T
  );
  plot(roc_rt, time = pre_time, col = bioCol[i-2], title = F, lwd = 4, add = T);
  auc_text <- c(auc_text, paste0(colnames(rt)[i], ", AUC=", sprintf("%.3f", roc_rt$AUC[2])));
}
legend("bottomright", auc_text, lwd = 4, bty = 'n', title = "All set", col = bioCol[1:(ncol(rt)-1)]);
dev.off();
```
![ROC曲线5](./md-image/ROC曲线5.png){:width=300 height=300}
可以看到使用风险得分进行roc预测的准确度明显大于其它临床特征
### Venn图绘制
除了使用代码绘制，还可以使用[在线网站](https://bioinfogp.cnb.csic.es/tools/venny/index.html)
将4种差异表达分析结果基因输入，并调整相关设置
![Venn图绘制](./md-image/Venn图绘制.png){:width=400 height=400}
之后右键图片->`将图像另存为`
除此之外，点击图片中的各区域，左下角results就可列出这部分交叉的基因
![Venn图绘制2](./md-image/Venn图绘制2.png){:width=350 height=350}
### 列线图
使用`regplot`包的版本为0.2，[下载](https://cran.r-project.org/src/contrib/Archive/regplot/regplot_0.2.tar.gz)
``` r
library(survival);
if(!require("regplot", quietly = T))
{
  install.packages("regplot");
  library("regplot");
}
library(rms);
```
**读取风险得分和临床数据，合并**：
``` r
# 风险得分
risk <- read.table("C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\risk.txt", check.names = F, row.names = 1, sep = '\t', header = T);
# 临床数据
library("readxl");
library(tibble);
cli <- read_excel("C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\new_clinical.xlsx");
cli <- column_to_rownames(cli, "ID");
# 删除包含NA的行
cli <- cli[apply(cli, 1, function(x)any(is.na(match(NA, x)))),  ];
# 合并
same_sample <- intersect(rownames(risk), rownames(cli));  # 共同样本名
risk <- risk[same_sample, ];
cli <- cli[same_sample, ];  # 过滤
rt <- cbind(risk[, c("time", "state", "risk")], cli[, c("Age", "Gender", "T", "N", "M", "Stage", "subdivision")]);  # 合并
```
![列线图2](./md-image/列线图2.png){:width=150 height=150}
**画列线图**：
``` r
# 因为我们这里有生存时间和生存状态，所以用cox
# 如果只是用来预测二分类（如正常/肿瘤），就用logistic回归
res.cox <- coxph(Surv(time, state) ~ ., data = rt);
pdf(file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\nomogram.pdf", width = 10, height = 11);
# 画图
nom1 <- regplot(
  res.cox,
  title = "",
  points = T,
  droplines = T,
  observation = rt[50, ],
  failtime = c(1, 3, 5),
  prfail = F
);
dev.off();
```
![列线图3](./md-image/列线图3.png){:width=700 height=700}
更高版本的`regplot`包画出的列线图：
![列线图1](./md-image/列线图1.png){:width=600 height=600}
图左侧标明了画图过程中使用的变量，图右侧的points区域是每个变量相应能获得多少分，根据所有变量分数之和可以获得总体的打分`Total Points`，这个总体打分可以预测每个样本存活1/3/5年的概率
**列线图打分**：
``` r
nomo_risk <- predict(res.cox, data = rt, type = "risk");
rt <- cbind(Nomogram = nomo_risk, risk);
outTab <- rbind(ID = colnames(rt), rt);
write.table(outTab, file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\nomoRisk.txt", sep = '\t', row.names = F, quote = F);
```
![列线图4](./md-image/列线图4.png){:width=150 height=150}
**校准曲线**：
``` r
pdf(file = "C:\\Users\\WangTianHao\\Documents\\GitHub\\R-for-bioinformatics\\b站生信课03\\save_data\\calibration.pdf", width = 6, height = 6);
# 1年
f <- cph(Surv(time, state) ~ Nomogram, x = T, y = T, surv = T, data = rt, time.inc = 1);
cal <- calibrate(f, cmethod = 'KM', method = "boot", u = 1, m = (nrow(rt)/3), B = 1000);
plot(
  cal, 
  xlim = c(0, 1), 
  ylim = c(0, 1), 
  xlab = "Nomogram-predicted OS (%)", 
  ylab = "Observed OS (%)",
  lwd = 3,
  col = "Firebrick2",
  sub = F
);
# 3年
f <- cph(Surv(time, state) ~ Nomogram, x = T, y = T, surv = T, data = rt, time.inc = 3);
cal <- calibrate(f, cmethod = 'KM', method = "boot", u = 3, m = (nrow(rt)/3), B = 1000);
plot(
  cal, 
  xlim = c(0, 1), 
  ylim = c(0, 1), 
  xlab = "Nomogram-predicted OS (%)", 
  ylab = "Observed OS (%)",
  lwd = 3,
  col = "MediumSeaGreen",
  sub = F,
  add = T
);
# 5年
f <- cph(Surv(time, state) ~ Nomogram, x = T, y = T, surv = T, data = rt, time.inc = 5);
cal <- calibrate(f, cmethod = 'KM', method = "boot", u = 5, m = (nrow(rt)/3), B = 1000);
plot(
  cal, 
  xlim = c(0, 1), 
  ylim = c(0, 1), 
  xlab = "Nomogram-predicted OS (%)", 
  ylab = "Observed OS (%)",
  lwd = 3,
  col = "NavyBlue",
  sub = F,
  add = T
);
# 图例
legend(
  "bottomright",
  c("1-year", "3-year", "5-year"),
  col = c("Firebrick2", "MediumSeaGreen", "NavyBlue"),
  lwd = 3,
  bty = 'n'
);
dev.off();
```
![列线图5](./md-image/列线图5.png){:width=400 height=400}
横坐标是这个列线图预测的总体生存率，纵坐标是实际生存率。当预测==实际（即曲线趋近对角线）时，预测效果较好
