<a id="mulu">目录</a>
<a href="#mulu" class="back">回到目录</a>
<style>
    .back{width:40px;height:40px;display:inline-block;line-height:20px;font-size:20px;background-color:lightyellow;position: fixed;bottom:50px;right:50px;z-index:999;border:2px solid pink;opacity:0.3;transition:all 0.3s;color:green;}
    .back:hover{color:red;opacity:1}
    img{vertical-align:bottom;}
</style>

<!-- @import "[TOC]" {cmd="toc" depthFrom=3 depthTo=6 orderedList=false} -->

<!-- code_chunk_output -->

- [ggplot2进阶](#ggplot2进阶)
    - [更多画子图的方式](#更多画子图的方式)
      - [cowplot包](#cowplot包)
      - [gridExtra包](#gridextra包)
      - [ggExtra包](#ggextra包)
    - [在图中写公式或统计信息](#在图中写公式或统计信息)
      - [如何写公式](#如何写公式)
      - [调整公式位置--hjust/vjust/angle](#调整公式位置-hjustvjustangle)
      - [希腊字符](#希腊字符)
      - [更多公式写法](#更多公式写法)
    - [先计算再做图--画图函数的计算方法](#先计算再做图-画图函数的计算方法)
    - [position参数--以柱形图为例](#position参数-以柱形图为例)
      - [堆叠柱状图(stacked bars)](#堆叠柱状图stacked-bars)
      - [并列柱状图](#并列柱状图)
    - [主题](#主题)
- [ggplot2实战](#ggplot2实战)
    - [swiss](#swiss)
    - [iris](#iris)
    - [箱型图](#箱型图)
    - [二维密度图](#二维密度图)
    - [facet子图](#facet子图)
    - [极坐标图(polar图)](#极坐标图polar图)
    - [多图组合](#多图组合)
      - [cowplot包](#cowplot包-1)
      - [gridExtra包](#gridextra包-1)
      - [拓展：patchwork包](#拓展patchwork包)
    - [拓展：lattice包](#拓展lattice包)
    - [与factor配合使用](#与factor配合使用)

<!-- /code_chunk_output -->

<!-- 打开侧边预览：f1->Markdown Preview Enhanced: open...
只有打开侧边预览时保存才自动更新目录 -->

### ggplot2进阶
##### 更多画子图的方式
###### cowplot包
``` r
if (!require("cowplot")){ 
  install.packages("cowplot");
  library("cowplot");
} 
```
``` r
cowplot::plot_grid(
    子图对象1, 子图对象2, ... , 
    labels=c(图1名称, 图2名称, ...), 
    ncol = 总列数, 
    nrow = 总行数,
    byrow = TRUE,  # 第二个子图默认向右画
    ...
)
```
注意：`byrow`决定的是子图对象的位置，不决定`labels`的位置，`labels`都是向右延伸（第二个labels都在第一个labels的右面）
``` r
sp <- ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl)))+
  geom_point(size=2.5);
bp <- ggplot(diamonds, aes(clarity, fill = cut)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=70, vjust=0.5));
cowplot::plot_grid(sp, bp, sp, labels=c("sp", "bp", "sp"), ncol = 2, nrow = 2);
cowplot::plot_grid(sp, bp, sp, labels = c("sp", "bp", "sp"), ncol = 2, nrow = 2, byrow = F);
```
![cowplot包1](./md-image/cowplot包1.png){:width=400 height=400}
可以看到图的顺序发生了改变，而左上角标签的位置不变

---

**指定每个图的大小和位置**：`draw_plot(plot, x = 0, y = 0, width = 1, height = 1)`其中xy以及宽高取值均为0-1间，表示占整个图的百分比，坐标以左下角为原点
``` r
# 使用上面画的sp bp
plot.iris <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) + 
  geom_point() + 
  facet_grid(. ~ Species) + 
  stat_smooth(method = "lm") +
  background_grid(major = 'y', minor = "none") + # add thin horizontal lines 
  panel_border();
(plot <- 
  ggdraw() +  # 相当于创建一个大的画板
  draw_plot(plot.iris, x=0, y=.5, width=1, height=.5) +
  draw_plot(sp, 0, 0, .5, .5) +
  draw_plot(bp, .5, 0, .5, .5));  # 分别画在大画板的哪个位置，以及高度宽度分别是多少
```
![cowplot包2](./md-image/cowplot包2.png){:width=400 height=400}
**为每个子图设置标签**：`draw_plot_label(c(标签名), c(标签x轴位置), c(标签y轴位置))`
还可以设置字体颜色color、字体尺寸size等
``` r
plot +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15);
```
表示"A"在x=0 y=1的位置，"B"在x=0 y=0.5的位置，"C"在x=0.5 y=0.5的位置
![cowplot包3](./md-image/cowplot包3.png){:width=400 height=400}
###### gridExtra包
``` r
if (!require("gridExtra")){ 
  install.packages("gridExtra");
  library("gridExtra");
} 
```
**先创建子图对象**：
``` r
df <- ToothGrowth
df$dose <- as.factor(df$dose)
bp <- ggplot(df, aes(x=dose, y=len, color=dose)) +
  geom_boxplot() + 
  theme(legend.position = "none") + 
  labs( tag = "A");
dp <- ggplot(df, aes(x=dose, y=len, fill=dose)) +
  geom_dotplot(binaxis='y', stackdir='center')+
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")+
  theme(legend.position = "none") + 
  labs( tag = "B");
vp <- ggplot(df, aes(x=factor(dose), y=len)) +
  geom_violin()+
  geom_boxplot(width=0.1) + 
  labs( tag = "C");
sc <- ggplot(df, aes(x=dose, y=len, color=dose, shape=dose)) +
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "none") +
  theme_gray() + 
  labs( tag = "D");
```
`grid.arrange(子图对象1, 子图对象2, ... , ncol, nrow)`指定行列数：
``` r
grid.arrange(bp, dp, vp, sc, ncol=2, nrow =2);
```
![gridExtra包1](./md-image/gridExtra包1.png){:width=400 height=400}
使用`layout_matrix`参数：由一个矩阵确定绘图区域的每个格子都画哪些图：
``` r
grid.arrange(bp, dp, vp, sc, ncol = 2, 
             layout_matrix = cbind(c(1,1,1), c(2,3,4)));
```
![gridExtra包2](./md-image/gridExtra包2.png){:width=400 height=400}
- `c(1,1,1)`是第一列画哪些图，全是1就全画第一张图
- `c(2,3,4)`是第二列画第2 3 4张图

`cbind(c(1,1,1), c(2,3,4))`矩阵的样子：
```
     [,1] [,2]
[1,]    1    2
[2,]    1    3
[3,]    1    4
```
可以这样理解：将整个画板按n行m列的矩阵分成n行m列的格子，画板的每个格子与矩阵每个位置的值相对应——矩阵中为1的位置都画第一张图，为2的位置画第二张图，...
**另一个例子**：
``` r
grid.arrange(bp, dp, vp, sc, ncol = 3, 
             layout_matrix = cbind(c(1,1), c(2,2), c(3,4)));
```
![gridExtra包3](./md-image/gridExtra包3.png){:width=400 height=400}

---

**gridExtra包与图例**：
实现效果：让两张图共享一个图例，图例显示在两张图的正上方
思路：把图例画成一张图，使用`grid.arrange`让第一行是图例，第二行是两张图
``` r
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
```
该函数接收一个图，返回这个图中的图例
``` r
# 第一张子图（有图例）
bp <- ggplot(df, aes(x=dose, y=len, color=dose)) +
  geom_boxplot() + 
  labs(tag = "A") +
  theme(legend.position = "top");  # 让图例在图的上方
# 第二张子图（无图例）
vp <- ggplot(df, aes(x=dose, y=len, color=dose)) +
  geom_violin()+ 
  geom_boxplot(width=0.1) + 
  labs( tag = "B") +
  theme(legend.position="none");  # 去掉图例
# 将第一张子图的图例提取出来
legend <- get_legend(bp);
# 再删去第一张子图的图例
bp2 <- bp + theme(legend.position="none");
```
绘图：
``` r
grid.arrange(legend, bp2, vp,  ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,1), c(2,3)),
             widths = c(2.7, 2.7), heights = c(0.2, 2.5));
```
![gridExtra包4](./md-image/gridExtra包4.png){:width=300 height=300}
其中widths和heights指定了第一行和第二行的宽高，第一行因为是图例，所以高度较小
`rbind(c(1,1), c(2,3))`矩阵的样子：
```
     [,1] [,2]
[1,]    1    1
[2,]    2    3
```
即第一行的两个格子全是第一张图的，第二行分别是第2/3张图
###### ggExtra包
用于向已有图中添加边缘直方图(marginal histograms)，展示数据的分布状况
``` r
if (!require("ggExtra")){ 
  install.packages("ggExtra");
  library("ggExtra");
} 
```
一个基本图：
``` r
(piris <- ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
  geom_point());
```
![ggExtra包1](./md-image/ggExtra包1.png){:width=300 height=300}
添加边缘直方图：
``` r
ggMarginal(piris, groupColour = TRUE, groupFill = TRUE);
```
![ggExtra包2](./md-image/ggExtra包2.png){:width=300 height=300}
##### 在图中写公式或统计信息
###### 如何写公式
先看一个例子：
``` r
m = lm(Fertility ~ Education, swiss);  # 回归拟合分析
c = cor.test( swiss$Fertility, swiss$Education );  # 相关性分析
eq <- substitute(  # 写公式
  atop( 
    paste( italic(y), " = ",  a + b %.% italic(x), sep = ""),
    paste( italic(r)^2, " = ", r2, ", ", italic(p)==pvalue, sep = "" ) 
  ),
  list(
    a = as.vector( format(coef(m)[1], digits = 2) ),
    b =  as.vector( format(coef(m)[2], digits = 2) ),
    r2 =  as.vector( format(summary(m)$r.squared, digits = 2) ),
    pvalue =  as.vector( format( c$p.value , digits = 2) ) 
  )
); 
eq <- as.character(as.expression(eq));  # 先把eq转成公式，再变成字符串，以便写入图中
ggplot(swiss, aes( x = Education,  y = Fertility ) ) +  # 画图
  geom_point( shape = 20 ) +  # 散点图
  geom_smooth( se = T ) +  # 拟合曲线
  geom_text( data = NULL,  # 写文字
             aes( x = 30, y = 80, label= eq, hjust = 0, vjust = 1),  # 文字的内容和位置
             size = 4, parse = TRUE, inherit.aes=FALSE);
```
![在图中写公式或统计信息1](./md-image/在图中写公式或统计信息1.png){:width=300 height=300}
- `atop(<equation_1> , <equation_2>)`将两个公式上下放置，返回一个`<equation>`。其中`italic`函数将普通字符转为斜体加粗的公式变量形式
- `substitute(<equation>, list(变量名=值,...))`将公式中的变量名替换为数值，如上例中就是替换了`a` `b` `r2` `pvalue`
- `geom_text`
  - `data`是否需要继承前面的画图数据`swiss`
  - `aes(x, y, label)`标签的`x`/`y`位置及内容`label`
  - `size`文字大小
  - `family`字体
  - `nudge_x`/`nudge_y`设置文字距原坐标点的距离
  - `check_overlap=T`设置不画与同一层中的上一个文本重叠的文本
  - `parse=T`将字符串表示成公式形式
  - `inherit.aes=F`不继承之前的aes设置

**公式的写法1**：
``` r
paste( italic(y), " = ",  a + b %.% italic(x), sep = "")
paste( italic(r)^2, " = ", r2, ", ", italic(p)==pvalue, sep = "" )
```
![在图中写公式或统计信息2](./md-image/在图中写公式或统计信息2.png){:width=200 height=200}
**公式的写法2**（不使用paste函数拼接）：
``` r
italic(y) == a + b %.% italic(x)
italic(r)^2~"="~r2*","~italic(p)==pvalue
```
![在图中写公式或统计信息3](./md-image/在图中写公式或统计信息3.png){:width=300 height=300}
在这种方法中，引号两边必须有`*`或`~`字符，`~`表示空格，`*`表示什么都没有，`~~`表示两个空格
**公式中的代数负号**：
![在图中写公式或统计信息8](./md-image/在图中写公式或统计信息8.png){:width=500 height=500}
###### 调整公式位置--hjust/vjust/angle
[参考文章](https://juejin.cn/post/7130987741747085349)
- `hjust/vjust`调整水平/垂直方向的位置
  - `hjust`取值通常为0-1，0是左对齐（默认）、0.5是居中、1是右对齐
  - `vjust`取正值就是垂直向下移动，负值是向上移动
- `angle`使文本绕中心点旋转一定角度，取值为整数（正--逆时针、负--顺时针），单位为°

这三个参数可被用于各种需要写文字的地方，如坐标轴标签、文本注释等等
**例1**：
``` r
df = data.frame(team=c('The Amazing Amazon Anteaters',
                       'The Rowdy Racing Raccoons',
                       'The Crazy Camping Cobras'),
                points=c(14, 22, 11));
ggplot(data=df, aes(x=team, y=points)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=90)) ;
```
![在图中写公式或统计信息4](./md-image/在图中写公式或统计信息4.png){:width=400 height=400}
这里使用了`angle=90`让x轴标签变为竖向，避免堆叠
可以使用hjust和vjust参数来调整x轴标签，使其与x轴上的刻度线更紧密地排列：
``` r
ggplot(data=df, aes(x=team, y=points)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1));
```
![在图中写公式或统计信息5](./md-image/在图中写公式或统计信息5.png){:width=400 height=400}
可以看到x轴标签向左/上移动了一些
**例2**：
``` r
df <- data.frame(player=c('Brad', 'Ty', 'Spencer', 'Luke', 'Max'),
                 points=c(17, 5, 12, 20, 22),
                 assists=c(4, 3, 7, 7, 5));
ggplot(df) +
  geom_point(aes(x=points, y=assists)) + 
  geom_text(aes(x=points, y=assists, label=player));
```
![在图中写公式或统计信息6](./md-image/在图中写公式或统计信息6.png){:width=300 height=300}
将文字向下移动使更容易阅读：
``` r
ggplot(df) +
  geom_point(aes(x=points, y=assists)) + 
  geom_text(aes(x=points, y=assists, label=player), vjust=1.2);
```
![在图中写公式或统计信息7](./md-image/在图中写公式或统计信息7.png){:width=300 height=300}
###### 希腊字符
使用`geom_text(aes(label="alpha"), parse=T)`
希腊字符的英文写法：
![在图中写公式或统计信息9](./md-image/在图中写公式或统计信息9.png){:width=300 height=300}
如何画出这个图？
思路：先画出4x6个点，再对点添加文本
准备数据：
``` r
greeks <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta",
            "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu",
            "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma",
            "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega");
dat <- data.frame( x = rep( 1:6, 4 ), y = rep( 4:1, each = 6), greek = greeks );  # rep是将数组复制多少次  
```
![在图中写公式或统计信息10](./md-image/在图中写公式或统计信息10.png){:width=300 height=300}
即每个文本对应的xy坐标
绘图：
``` r
plot2 <- 
  ggplot( dat, aes(x=x,y=y) ) + 
  geom_point(size = 0) +
  geom_text( aes( x, y + 0.1, label = tolower( greek ) ), size = 10, parse = T ) +
  geom_text( aes( x, y - 0.1, label = tolower( greek ) ), size = 5 );
```
注意两个`geom_text`的`label`都相同，但一个是写希腊符号，另一个是写普通英文，区别是`parse = T`参数，它控制是否要对字符串进行公式化转换
###### 更多公式写法
**例1**：分数、根号、指数
``` r
eq <- expression(
  paste(
    frac(1, sigma*sqrt(2*pi)), 
    " ",
    plain(e)^{frac(-(x-mu)^2, 2*sigma^2)}
  )
);
```
![更多公式写法1](./md-image/更多公式写法1.png){:width=150 height=150}
**例2**：`bquote`函数
``` r
x <- 1.24;
y <- 0.6;
ex <- bquote(
  .(
    parse(
        text = paste( 
          "observed (", 
          "italic(R)^2==",
          x,  
          "^bold(", x, "), 
          n == ", y, 
          ")", 
          sep = "  " 
        )
    )
  ) 
);
```
![更多公式写法2](./md-image/更多公式写法2.png){:width=100 height=100}
**例3**：`ggtitle()`指定标题不用写`parse = T`
``` r
x_mean <- 1.5;
x_sd <- 1.2;
ex <- substitute(
    paste(
      X[i], 
      " ~ N(", mu, "=", m, ", ", 
      sigma^2, "=", s2, ")"
    ),
    list(m = x_mean, s2 = x_sd^2)
);
ggplot( data.frame( x = rnorm(100, x_mean, x_sd) ), aes( x ) ) +
    geom_histogram( binwidth=0.5 ) +  # 直方图
    ggtitle(ex);  # 添加标题
```
![更多公式写法3](./md-image/更多公式写法3.png){:width=300 height=300}
##### 先计算再做图--画图函数的计算方法
先看一个例子：
``` r
# 准备数据
grades2 <- read_delim( file = "data/grades2.txt", delim = "\t",
                       quote = "", col_names = T);
knitr::kable( grades2 );
```
![默认计算方法1](./md-image/默认计算方法1.png){:width=200 height=200}
现在要画出每位学生及格的课程数，使用filter筛选行即可
``` r
ggplot( grades2 %>% dplyr::filter( grade >= 60 ), 
        aes( name ) ) +
  geom_bar();
```
![默认计算方法2](./md-image/默认计算方法2.png){:width=300 height=300}
仔细观察上面的函数，我们只是筛选了行并给出了x轴取值，并没有分组和计算行数，但`geom_bar()`仍给出了正确的结果
这是因为`geom_bar()`中有一个默认的参数`stat = "count"`，它表示画图函数对数据的统计方法，柱状图默认是按x轴的数据分组并计算行数
以上函数实际相当于：
``` r
# 先做统计
cnt <- grades2 %>% 
  group_by( name ) %>% 
  summarise( cnt = sum( grade >= 60 )  );
ggplot( cnt, aes( x = name, y = cnt ) ) +
  geom_bar( stat = "identity" );  # 取消默认统计方法，按指定的df画图
```
其它画图函数的默认计算方法：
- `geom_bar` : count
- `geom_boxplot` : boxplot
- `geom_count` : sum
- `geom_density` : density
- `geom_histogram` : bin
- `geom_quantile` : quantile
##### position参数--以柱形图为例
###### 堆叠柱状图(stacked bars)
应用场景：宏基因组多样本物种丰度图
![堆叠柱状图1](./md-image/堆叠柱状图1.png){:width=300 height=300}
方法：为`geom_bar()`指定参数`position = "stack"`
例：
``` r
# 准备数据
speabu <-read_tsv( file = "data/mock_species_abundance.txt"  );
head( speabu );
```
![堆叠柱状图2](./md-image/堆叠柱状图2.png){:width=150 height=150}
相同id的画一个柱子，柱子堆叠的每部分颜色不同（根据genus列来取），柱子高度由abundance列决定
``` r
ggplot( speabu, aes( x = id, y = abundance, fill = genus ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "black", width = 0.2 );
```
![堆叠柱状图3](./md-image/堆叠柱状图3.png){:width=300 height=300}
注意fill是填充色，而color是边框色

---

**需求1：指定Genus展示顺序**
思路：使用factor对Genus列进行重排
``` r
speabu$genus <- factor( speabu$genus, 
                        levels = rev( c( "Enterobacteriaceae", "Lachnospiraceae", "Bacteroidaceae", "Lactobacillaceae", "Clostridiaceae", 
                        "Ruminococcaceae", "Prevotellaceae", "Erysipelotrichaceae", "Streptococcaceae", "Enterococcaceae", "Other" ) ) );
ggplot( speabu, aes( x = id, y = abundance, fill = genus ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "black", width = 0.8 );
```
![堆叠柱状图4](./md-image/堆叠柱状图4.png){:width=300 height=300}

---

需求2：按丰度中值大小排序
使用reorder函数：`reorder(返回结果列, 顺序决定列, 排序方法)`
``` r
speabu$genus <- reorder( speabu$genus, speabu$abundance, median );
ggplot( speabu, aes( x = id, y = abundance, fill = genus ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "white", width = 0.8 );
```
![堆叠柱状图5](./md-image/堆叠柱状图5.png){:width=300 height=300}

---

**需求3：显示数值**
即每个堆叠部分的占比（abundance值）
``` r
# 先计算显示位置
speabu <- speabu %>% 
  arrange( id, desc( factor( genus ) ) ) %>%  # 按id和genus排序
  group_by( id ) %>%  # 分组，让相同id的为一个柱子，每组单独计算位置
  mutate( ypos = cumsum( abundance ) - abundance / 2 );  # 累加，再减去一半的高度使文字居中
ggplot( speabu, aes( x = id, y = abundance, fill = genus ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "black", width = 0.8 ) +
  geom_text( aes( y = ypos, label = paste( abundance, "%", sep = "" ) ), color = "white" );
```
![堆叠柱状图6](./md-image/堆叠柱状图6.png){:width=300 height=300}
###### 并列柱状图
即让每个柱子不堆叠，而是相邻排列
![不堆叠的柱状图1](./md-image/不堆叠的柱状图1.png){:width=300 height=300}
方法：为`geom_bar()`指定参数`position = "dodge"`或`position=position_dodge()`
``` r
ggplot( speabu, aes( x = id, y = abundance, fill = genus ) ) + 
  geom_bar( stat = "identity", position = "dodge", color = "white", width = 0.8 );
```

---

如何为并列的柱状图显示数值：使用`position_dodge`位置调整文本位置
``` r
# 准备数据
df2 <- data.frame(supp=rep(c("VC", "OJ"), each=3),
               dose=rep(c("D0.5", "D1", "D2"),2),
              len=c(6.8, 15, 33, 4.2, 10, 29.5));
```
![不堆叠的柱状图2](./md-image/不堆叠的柱状图2.png){:width=150 height=150}
``` r
ggplot( df2, aes(x=factor(dose), y=len, fill=supp)) +
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_text(aes(label=len), vjust=1.6, color="black", size=3.5,
            position = position_dodge(0.9))
```
![不堆叠的柱状图3](./md-image/不堆叠的柱状图3.png){:width=150 height=150}
通过修改`position_dodge`里面的参数，可以调整文字和柱子的位置

---

position的其它取值：
- `position = position_identity()`：在指定位置，不改变
- `position = position_jitter()`：随机往别的地方移动使点不重叠
- `position = position_nudge()`：平移

[更多关于position系列函数](https://blog.csdn.net/weixin_54000907/article/details/120108707)
不同的画图函数有不同的默认position取值，比如柱形图是堆叠、箱型图是相邻
``` r
ggplot(ToothGrowth, aes(x=factor( dose ), y=len, fill=supp)) +
  geom_boxplot();
```
![不堆叠的柱状图4](./md-image/不堆叠的柱状图4.png){:width=300 height=300}
##### 主题
一般分为两种调整主题的方式：
- `theme(...)`自定义各个元素的样式
- `theme_xxx()`直接使用已经定制好的内容，包括`theme_bw`、`theme_linedraw`、`theme_light`、`theme_dark`、`theme_minimal`、`theme_classic`、`theme_void`和默认主题`theme_gray`等

``` r
ggplot(ToothGrowth, aes(x=factor( dose ), y=len, fill=supp)) +
  geom_boxplot() + scale_fill_brewer( palette = "Paired" ) + theme_classic();
```
![主题1](./md-image/主题1.png){:width=300 height=300}

---

`theme()`函数：可用通式`theme(主题.部件=element_类型()/具体值)`表示
- 主题：`plot`整幅图、`axis`坐标轴、`legend`图例、`panel`面板、`facet`子图
- 部件：`title`名字，坐标轴名字、`line`线，坐标轴的xy轴, `text`标签，坐标轴刻度的数字、`ticks`坐标轴刻度的小线条、`background`背景、`position`位置、...
- 类型：`rect`区域、`line`线条、`text`文本

其中部件要和类型一致。如部件为title、text等文字相关的元素，那么类型处就应为text
[更多关于theme函数的应用](https://blog.csdn.net/zty0104/article/details/119646934)

---

除了theme函数，还有一些函数用于设置图样式，如`labs`：
``` r
labs(
  x = "<x label>",  # x轴标签
  y = "<y label>",  # y轴标签
  tag = "图标签",  # 常用于为每个子图添加对应标签
  colour = "<legend title>",  # 与aes里的colour配合使用
  fill = "<legend title>",  # 与aes里的fill配合使用
  shape = "<legend title>",  # 与aes里的shape配合使用
  ...,
)
```
`colour/fill/shape`参数是为了给图例加标签，如果图在aes中按fill填充色进行区分并画图例，就使用`fill`参数给图例加标签
``` r
ggplot(ToothGrowth, aes(x=factor( dose ), y=len, fill=supp)) +
  geom_boxplot() + 
  scale_fill_brewer( palette = "Paired" ) + 
  labs( fill = "图例", x = "Dose (mg)" ) + 
  theme( legend.position = "top" );  # 调整图例位置为顶部
```
![主题2](./md-image/主题2.png){:width=300 height=300}
`labs`还可以同时为多个图例指定名称：
``` r
# 创建基础图
df <- ToothGrowth;
df$dose <- as.factor(df$dose);
sc <- ggplot(df, aes(x=dose, y=len, color=dose, shape=dose)) +
  geom_jitter(position=position_jitter(0.2))+
  theme(legend.position = "none") +
  theme_gray();
library("gridExtra");
grid.arrange(
  sc + labs( tag = "A" ), 
  sc + labs( colour = "Dose (mg)" , tag = "B" ),
  sc + labs( shape = "Dose (mg)" , tag = "C" ),
  sc + labs( colour = "Dose (mg)", shape = "Dose (mg)", tag = "D" ),
  ncol=4, nrow =1
);
```
![主题3](./md-image/主题3.png){:width=300 height=300}
### ggplot2实战
##### swiss
用`swiss`数据做图
每道题都有R基础作图函数和ggplot2两个版本
**1.用直方图显示Catholic列的分布情况**
基础作图函数：
``` r
hist(
  x = swiss$Catholic,  # 数据
  main = "Catholic",  # 标题
  xlab = "Catholic"  # x轴标签
);
```
![swiss1](./md-image/swiss1.png){:width=300 height=300}
ggplot2：
``` r
ggplot(
  swiss,  # 使用的数据集
  aes( x = Catholic )  # 标明x轴的取自哪列
  ) +
  geom_histogram();  # 画直方图
```
![swiss2](./md-image/swiss2.png){:width=300 height=300}
**2.用散点图显示Eduction与Fertility的关系**，同时将表示两者关系的线性公式、相关系数和p值画在图的空白处
因为要标注公式、相关系数和p值，所以首先要进行计算：
``` r
m <- lm(Fertility ~ Education,swiss);  # 线性回归
cor_val <- cor(swiss$Fertility, swiss$Education);  # 相关系数
p_val <- cor.test(swiss$Fertility,swiss$Education)$p.value;  # p值
a <- coef(m)[1];  # 回归线y=ax+b中的a
b <- coef(m)[2];  # 回归线y=ax+b中的b
```
基础作图函数：
``` r
with(swiss,plot(Education, Fertility, col = "blue"))  # 绘制基础散点图
abline(  # 添加拟合曲线
  lm(swiss$Fertility ~ swiss$Education),  # 回归线
  col="red"  # 线颜色
);
legend(  # 将文字作为图例进行添加
  "topright",  # 添加在右上角
  legend = paste(  # 图例内容
    "cor =", round(cor_val,2),  # 相关系数--保留2位小数
    "\n",
    "p =", format(p_val, scientific = T, digits = 2),  # p值--科学计数法，保留2位小数
    "\n",
    "y =",format(a ,digits = 2),"+",format(b ,digits = 2),"x",  # 公式，系数保留2位小数
    "\n"
  )
);
```
![swiss3](./md-image/swiss3.png){:width=300 height=300}
ggplot2：
``` r
ggplot(swiss, aes(x = Education, y = Fertility)) +
  geom_point(color = "blue", size = 3) +  # 绘制基础散点图
  geom_smooth(method = "lm", color = "red") +  # 添加拟合曲线
  labs(title = "Education~Fertility", x = "Education", y = "Fertility") +  # 标题和xy轴
  # 可以用之前讲过的as.expression(eq)方法添加，这里列举另一种添加方法
  annotate(  # 将文字作为注释进行添加
    "text",  # 声明注释的类型为文本
    x = 20, y = 35,  # 添加注释的位置
    label = paste(  # 注释内容--与基础作图函数中的相同
      "cor =", round(cor_val,2),
      "\n",
      "p =", format(p_val, scientific = T, digits = 2),
      "\n",
      "y =",format(a ,digits = 2),"+",format(b ,digits = 2),"x"
    )
  );
```
![swiss4](./md-image/swiss4.png){:width=300 height=300}
##### iris
用`iris`数据做图
每道题都有R基础作图函数和ggplot2两个版本
**1.用箱型图显示不同species的Sepal.Length的分布情况**
即x轴为species，y轴为Sepal.Length
基础作图函数：
``` r
boxplot(Sepal.Length ~ Species, data = iris);
```
![iris1](./md-image/iris1.png){:width=300 height=300}
ggplot2：
``` r
ggplot(iris, aes(x = Species, y = Sepal.Length))+  # 指定数据集和xy轴
  geom_boxplot();  # 画箱型图
```
![iris2](./md-image/iris2.png){:width=400 height=400}
**2.用散点图显示Sepal.Length和Petal.Length之间的关系**，并根据species为点添加颜色，同时显示图例标明哪种颜色对应哪种species
基础作图函数：
因为基础作图函数不能自动根据某列分组指定颜色，需要我们手动建立物种与颜色的对应，方法：给vector中的元素命名，使其成为一个类似于字典的数据结构
``` r
colors <- c(
  "setosa" = "black",
  "versicolor" = "red",
  "virginica" = "green"
);  # 物种与颜色的对应
with(
  iris,
  plot(  # 散点图
    Sepal.Length, Petal.Length,  # 数据
    col = colors[Species],  # 按物种分配颜色
    pch = 19,  # 设置颜色为填充色，而不是边缘颜色
    xlab = "Sepal.Length",
    ylab = "Petal.Length"  # 设置xy轴标题
  )
);
legend(
  "topleft",  # 设置图例位置
  legend = levels(iris$Species),  # 设置图例中的标签来自iris$Species
  col = colors,  # 按之前的对应设置颜色
  pch = 19,  # 颜色为填充色
  cex = 0.8,  # 调整图例大小
  x.intersp = 0.3, 
  y.intersp = 0.3,  # 图例每个标签的xy间距
  text.width = 0.3  # 文本宽度
);
```
![iris3](./md-image/iris3.png){:width=400 height=400}
注：`colors[Species]`是在画图时按`colors`这个对应关系给物种标注颜色。如果只写`col=colors`，效果等效于`colors <- c("red","green","blue")`（每个点循环使用这3个颜色）
ggplot2：
``` r
ggplot(
  iris,  # 数据
  aes(
    x = Sepal.Length,
    y = Petal.Length,  # 指定xy轴数据
    color = Species  # 按物种分颜色
  ),
  labs(x = Sepal.Length,y = Petal.Length)  # 设置xy轴标签
) +
  geom_point() +  # 画散点图
  theme(legend.position = "top");  # 设置图例位置
```
![iris4](./md-image/iris4.png){:width=400 height=400}
##### 箱型图
使用ggplot2：显示`starwars`中身高`height`与性别`gender`的关系
要求：
- 去掉height为NA的数据
- 用ggsignif包计算两种性别的身高是否有显著区别，并在图上显示
- 将此图的结果保存为变量`p1`

思路：使用`dplyr::filter`筛选行，条件为没有na的height和gender；之后使用`geom_boxplot`创建箱型图，并使用`geom_signif`函数进行分析
``` r
starwars_noNA <- starwars %>%  # 获取没有NA的数据
  dplyr::filter(
    is.na(height)==FALSE & is.na(gender)==FALSE
  );
base_boxplot <- 
  ggplot(  # 传入数据
    starwars_noNA,
    aes(x = gender, y = height)
  ) +
  geom_boxplot() +  # 箱型图
  labs(x = "gender", y = "height");  # xy轴标签
library(ggsignif);
p1 <- base_boxplot +
  geom_signif(
    comparisons = list(c("feminine", "masculine")),  # 设置需要比较的组
    map_signif_level = T  # 将图像上方的p值以*代替
  );
p1;
```
![箱型图](./md-image/箱型图.png){:width=400 height=400}
##### 二维密度图
使用ggplot2
用二维密度图和散点图显示`iris`中`Sepal.Length`列和`Sepal.Width`之间的关系，同时按`Species`分组
- 二维密度图：使用`geom_density2d`函数
- 将此图的结果保存为变量`p2`

``` r
p2 <- ggplot(  # 传入数据
  iris,
  aes(
    x = Sepal.Length,
    y = Sepal.Width,  # xy轴数据
    color = Species,  # 颜色按物种分组
    shape = Species  # 点形状也按物种分组
  )
) +
  geom_density2d() +  # 画二维密度图
  geom_point() +  # 画散点图
  labs(title = "IRIS");  # 添加标题
p2;
```
![二维密度图](./md-image/二维密度图.png){:width=400 height=400}
##### facet子图
使用ggplot2
**1.画散点图，显示`mtcars`中`wt`和`mpg`之间的关系，同时按`cyl`将数据分组画多个子图**
![facet子图1](./md-image/facet子图1.png){:width=400 height=400}
- 画出拟合曲线、图例、各图的名称
- 图中的`468`组为所有数据合在一起的结果
- 将此图的结果保存为变量`p3`

思路：新建一个数据集`mtcars_new`，由`mtcars`和`mtcars_468`合并而来。`mtcars_new`就是把`mtcars`中的`cyl`列值全改成`468`，因为`facet_wrap`函数要按`cyl`分组，这样才能分为4、6、8、468四组，其中468组包含全部mtcars全部数据
``` r
rm(mtcars);  # 重置mtcars数据集
mtcars_old <- mtcars;  # 保存原mtcars
mtcars$cyl <- 468;  # 更改列值
mtcars_468 <- mtcars;  # 保存新mtcars
mtcars_new <- dplyr::bind_rows(mtcars_old, mtcars_468)#合并
p3 <- ggplot(  # 传入数据
  mtcars_new,
  aes(
    x = wt, y = mpg,
    col = factor(cyl)  # 颜色按cyl分组
  )
) +
  geom_point() +  # 画散点图
  geom_smooth() +  # 添加拟合曲线
  scale_color_manual(  # 手动设置颜色
    breaks = c("4", "6", "8", "468"),
    values = c("red", "green", "blue", "purple")
    # breaks与values一一相对（cyl=4的颜色是red）
  ) +
  labs(x = "Weight",y = "MPG") +  # 添加xy轴标签
  facet_wrap(  # 画子图
    .~cyl,  # 按cyl分成多列
    ncol = 2,  # 列数为2
    scales = "free",  # 设置xy轴不统一
    dir = "h"  # 按水平方向排布，即第一行是46第二行是8，若是v第一行就是48
  );
p3;
```
![facet子图2](./md-image/facet子图2.png){:width=500 height=500}
**2.画点线图，显示`airquality`中`Wind`和`Temp`之间的关系，同时按`Month`将数据分组画多个子图**
- 画出拟合曲线、图例、各图的名称
- 子图按2行3列组织

``` r
ggplot(  # 传入数据
  airquality,
  aes(
    x = Wind,y = Temp,
    col = factor(Month)  # 按月份分组
  )
) +
  geom_line() +  # 折线图
  geom_point() +  # 散点图
  geom_smooth() +  # 添加拟合曲线
  labs(x = "Wind", y = "Temp") +  # 添加xy轴标签
  facet_wrap(  # 画子图
    .~Month,  # 按Month分成多列
    ncol = 3,  # 列数为3
    scales = "free",  # 设置xy轴不统一
    dir = "h"  # 按水平方向排布
  );
```
![facet子图3](./md-image/facet子图3.png){:width=500 height=500}
##### 极坐标图(polar图)
使用`mtcars`的`mpg`列
- 先按`cyl`分组，由`cyl`确定每个柱子的颜色
- 每个组内按`mpg`排序
- 将此图的结果保存为变量`p4`

![polar图1](./md-image/polar图1.png){:width=300 height=300}
先看看我们要使用的数据：
``` r
rm(mtcars)  # 重置mtcars
mtcars <- mtcars %>% select(cyl, mpg);
```
![polar图2](./md-image/polar图2.png){:width=300 height=300}
可以看到行名是汽车名称，即结果图中的文字内容
**数据初处理**：
- 为了方便后续增添文本，将行名单独成一列，使用`rownames_to_column`函数
- 分组并排序

``` r
data <- mtcars %>%
  rownames_to_column() %>%
  arrange(cyl,mpg);
```
![polar图3](./md-image/polar图3.png){:width=600 height=600}
`arrange(cyl,mpg)`：以cyl为主序，mpg为次序进行排序。即先排cyl，再在其基础上排mpg。可以看成先按cyl分组，再在组内排序mpg
注：不能写`group_by(cyl)`再`arrange(mpg)`，因为group_by是针对summary这类汇总函数的，这样排序后结果是只按mpg排的，不能体现分组
``` r
data <- mtcars %>%
    group_by(cyl) %>%
    arrange(mpg);
```
![polar图4](./md-image/polar图4.png){:width=600 height=600}
**向量化索引列**：将rowname列（汽车名称）改成向量形式
``` r
data$rowname <- factor(
  data$rowname,
  levels = data$rowname
);
```
这样做是为了固定每行的顺序，防止画图函数使用默认排序方法将顺序变乱。如果不加这一段代码：
![polar图5](./md-image/polar图5.png){:width=300 height=300}
**计算文本角度**：使用`angle = 90 - 360 * (id - 0.5) / n()`，其中id是行索引（是第几行）
``` r
data <- data %>%
  mutate(
    id = row_number(),  # 行索引
    angle = 90-360*(id-0.5)/n()  # 文本角度
  );
```
![polar图6](./md-image/polar图6.png){:width=300 height=300}
**画图**：
- 传入数据：xy轴分别是汽车名称和`mpg`，填充颜色fill为`cyl`
- 使用`geom_bar`画基础柱状图，注意因为我们这里直接是按`cyl`的值确定柱子高度，所以要加上`stat = "identity"`
- 使用`coord_polar`将x轴变为圆形，这是画任何polar图必须的一步
- 添加文本，注意添加`aes(x=1:length(data$rowname)`为文本设置基础坐标
- 其它微调
  - `scale_y_continuous`修改y轴范围（即起始值和结束值），它们分别对应中间空白圆大小、柱子长短
  - 不显示主网格线（平行xy轴的那些刻度线）
  - 不显示刻度和xy轴刻度文本
  - 不显示xy轴名称
  - 去掉图背景的阴影（详见下面的函数）

``` r
p4 <- ggplot(  # 传入数据
  data,
  aes(
    x = rowname,
    y = mpg,
    fill = factor(cyl)  # 指定填充颜色
  )
) +
  geom_bar(stat = "identity") +  # 基础柱状图
  coord_polar(theta = "x",start = 0) +  # 调整坐标轴
  geom_text(  # 添加文本
    aes(
      x = 1:length(data$rowname),  # 基础坐标
      label = data$rowname,  # 文本内容
      angle = data$angle,  # 角度
      hjust = -0.1  # 标签距中心距离 
    )
  ) +
  scale_y_continuous(  # 修改y轴范围
    expand = c(0.01, max(data$mpg)/3)
  ) +
  theme(  # 不显示主网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(  # 不显示刻度和xy轴刻度文本
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x = "", y = "") +  # 不显示xy轴名称
  theme(  # 去掉图背景的阴影
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  );
p4;
```
![polar图8](./md-image/polar图8.png){:width=400 height=400}
一个细节：传入数据时为什么要用`fill = factor(cyl)`？
如果fill设置颜色时不加factor，使用`fill = cyl`：
![polar图7](./md-image/polar图7.png){:width=300 height=300}
这是因为画图函数将cyl看成了连续变量，而我们想让它作为离散变量，因此用factor处理

---

**代码汇总**：
``` r
rm(mtcars);
data <- mtcars %>%
  rownames_to_column() %>%
  arrange(cyl,mpg);
data$rowname <- factor(data$rowname, levels = data$rowname);
data <- data %>%
  mutate(id = row_number(),
         angle = 90-360*(id-0.5)/n());
p4 <- ggplot(data, aes(x = rowname, y = mpg,
                       fill = factor(cyl))) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "x",start = 0) +
  geom_text(aes(x = 1:length(data$rowname),
                label = data$rowname,
                angle = data$angle,
                hjust = -0.1)) +
  scale_y_continuous(expand = c(0.01, max(data$mpg)/3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank()) +
  labs(x = "", y = "") +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA));
p4;
```
##### 多图组合
###### cowplot包
用`cowplot::ggdraw`将p1、p2和p3按下面的方式组合在一起
![多图组合1](./md-image/多图组合1.png){:width=300 height=300}
[cowplot包的具体使用](#cowplot包)
关键：以**左下角**为原点，右、上为正方向，为每张图指定的坐标也是以左下角为准
- A(p3)：(0, 0)，宽为0.5，高为1
- B(p1)：(0.5, 0.5)，宽为0.5，高为0.5
- C(p2)：(0.5, 0)，宽为0.5，高为0.5
- 标签：ABC坐标分别为(0, 1)、(0.5, 1)、(0.5, 0.5)

``` r
library(cowplot);
ggdraw() +  # 创建画板
  draw_plot(p3, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(p1, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(p2, 0.5, 0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0.5), c(1, 1, 0.5), size = 10);  # 标签
```
![多图组合2](./md-image/多图组合2.png){:width=400 height=400}
###### gridExtra包
用`gridExtra::grid.arrange`函数将p1、p2和p4按下面的方式组合在一起
![多图组合3](./md-image/多图组合3.png){:width=300 height=300}
要求为子图加上ABC标签
[gridExtra包的具体使用](#gridextra包)
可以看到分为两列，第一行为第一张图，第二行是剩下的两张图。构建的画图矩阵应为：
```
     [,1] [,2]
[1,]    1    1
[2,]    2    3
```
``` r
library(gridExtra)
grid.arrange(
  p1+labs(tag = "A"),
  p2+labs(tag = "B"),
  p4+labs(tag = "C"),  # 先为每个图添标签，再组合在一起
  nrow = 2,
  layout_matrix = rbind(c(1,1),c(2,3))  # 画图矩阵
);
```
![多图组合4](./md-image/多图组合4.png){:width=400 height=400}
###### 拓展：patchwork包
用patchwork包中的相关函数将p1、p2、p3和p4 按下面的方式组合在一起
![多图组合5](./md-image/多图组合5.png){:width=400 height=400}
要求为子图加上ABCD标签
[patchwork包的使用](https://www.jianshu.com/p/73057774b4cb)
对于此题重要使用两个操作符：
- `p1 / p2`：将p1和p2竖直堆叠
- `p1 | p2`：将p1和p2并列放置
- 可以使用小括号来区分组合优先级

``` r
if (!require("patchwork")){ 
  install.packages("patchwork");
  library("patchwork");
} 
p1 <- p1+labs(tag = "A");
p3 <- p3+labs(tag = "B");
p4 <- p4+labs(tag = "C");
p2 <- p2+labs(tag = "D");  # 先为每个图添标签，再组合在一起
p1 / (p3 | p4) / p2;
```
![多图组合6](./md-image/多图组合6.png){:width=500 height=500}
##### 拓展：lattice包
[更多关于lattice包](https://www.jianshu.com/p/5885aaeda6c1)
lattice包内置于R，无需额外安装
它主要提供了绘制**网格图形**的方法。网格图形能够展示变量的分布或变量之间的关系，每幅图代表一个或多个变量的各个水平
以**散点图矩阵**(Scatter Plot Matrix)为例：
``` r
library(lattice);
lattice::splom( mtcars[c(1,3,4,5,6)] );
```
![lattice包](./md-image/lattice包.png){:width=400 height=400}
它可以在一张图中展示所有列间的关系。例如左上角第一个格子，它的纵坐标对应着`wt`，横坐标对应着`mpg`
##### 与factor配合使用
任务：
- 使用`readr`包中的函数读取`mouse_genes_biomart_sep2018.txt`
- 选取常染色体（1~19）和性染色体（X、Y）的基因
- 画各染色体上基因长度中值的箱型图，分别按
  - **染色体序号**排列(1, 2, 3,..., X, Y)
  - **基因长度的中值**排列（从短到长）

要求：
- 调整y轴范围，让每个箱子都完整显示
- 调整y轴范围时，不能使用`ylim()`，因为该函数会去除超出范围的值，使中值计算有偏差

解决方法：使用`coord_cartesian(ylim = c(min, max))`放大指定区域，这种方法只影响图形展示，不影响内部数据的值，而`ylim`函数会移除不在指定范围内的数据
画图思路：使用`geom_boxplot`函数，x轴为染色体名称（先排序），y轴为基因长度。在上一篇笔记中说过，只要将x轴数据设为factor形式，画箱型图函数便会自动将数据分组统计
**读取数据**：
``` r
mouse.genes <- read.delim( 
  file = "data\\mouse_genes_biomart_sep2018.txt",
  sep = "\t", header = T, stringsAsFactors = T 
);
mouse.genes %>% sample_n(5);
```
![与factor配合使用1](./md-image/与factor配合使用1.png){:width=150 height=150}
其中`Transcript.length..including.UTRs.and.CDS.`列是基因长度，`Chromosome.scaffold.name`是基因位置（在哪条染色体上）
**数据初处理**：筛选出在指定染色体的基因
``` r
chromosome <- c(as.character(1:19), "X", "Y");  # 需要哪些染色体的基因
mouse.chr <- mouse.genes %>% 
  filter(Chromosome.scaffold.name %in% chromosome);  # 筛选行
head(mouse.chr);
```
![与factor配合使用2](./md-image/与factor配合使用2.png){:width=150 height=150}
**第一个图**：按染色体序号排列
方法：将`Chromosome.scaffold.name`列转为factor，levels为染色体序号顺序
``` r
# 转为factor
x1 <- factor(mouse.chr$Chromosome.scaffold.name, levels = chromosome);
# 画图
ggplot(
  mouse.chr,  # 设定画图所用的数据集
  aes(  # 设定x,y轴对应的数据
    x = x1,
    y = Transcript.length..including.UTRs.and.CDS.
  )
) +
  geom_boxplot() +  # 画箱型图
  coord_cartesian(ylim = c(0,20000)) +  # 缩放y轴
  labs(  # 设置xy轴名称和标题
    y = "Transcript length (including UTRs and CDS)",
    x = "Chromosome/scaffold name",
    title = "按染色体序号排列"
  );
```
![与factor配合使用3](./md-image/与factor配合使用3.png){:width=400 height=400}
**第二个图**：按基因长度中值排序
方法：使用`reorder(返回结果列, 顺序决定列, 排序方法)`函数指定factor顺序，该函数在[堆叠柱状图(stacked bars)](#堆叠柱状图stacked-bars)中说过，这里不再说明
``` r
# 转为factor
x2 <- reorder(
  mouse.chr$Chromosome.scaffold.name,  # 返回结果列
  mouse.chr$Transcript.length..including.UTRs.and.CDS.,  # 顺序决定列
  median  # 排序方法
);
# 画图
ggplot(
  mouse.chr,
  aes(
    x = x2,
    y = Transcript.length..including.UTRs.and.CDS.
  )
) +  
  geom_boxplot() +
  coord_cartesian(ylim = c(0,20000)) +
  labs(
    y = "Transcript length (including UTRs and CDS)",
    x = "Chromosome/scaffold name",
    title = "按基因长度中值排列"
  );
```
![与factor配合使用4](./md-image/与factor配合使用4.png){:width=400 height=400}

---

**代码汇总**：
``` r
mouse.genes <- read.delim( 
  file = "..\\data\\mouse_genes_biomart_sep2018.txt",
  sep = "\t", header = T, stringsAsFactors = T );
chromosome <- c(as.character(1:19), "X", "Y"); 
mouse.chr <- mouse.genes %>% 
  filter(Chromosome.scaffold.name %in% chromosome);
# 按染色体序号排列↓
x1 <- factor(mouse.chr$Chromosome.scaffold.name, levels = chromosome);
ggplot(mouse.chr,
       aes(x = x1,
           y = Transcript.length..including.UTRs.and.CDS.)) +
  labs(y = "Transcript length (including UTRs and CDS)",
       x = "Chromosome/scaffold name",
       title = "按染色体序号排列") +
  geom_boxplot() +
  coord_cartesian(ylim = c(0,20000));
# 按中值排列↓
x2 <- reorder(mouse.chr$Chromosome.scaffold.name, mouse.chr$Transcript.length..including.UTRs.and.CDS., median);
ggplot(mouse.chr,
       aes(x = x2,
           y = Transcript.length..including.UTRs.and.CDS.)) +
  labs(y = "Transcript length (including UTRs and CDS)",
       x = "Chromosome/scaffold name",
       title = "按基因长度中值排列") +
  geom_boxplot() +
  coord_cartesian(ylim = c(0,20000));
```
