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

<!-- /code_chunk_output -->

<!-- 打开侧边预览：f1->Markdown Preview Enhanced: open...
只有打开侧边预览时保存才自动更新目录 -->

### ggplot2进阶
##### 更多画子图的方式
###### cowplot包
```
if (!require("cowplot")){ 
  install.packages("cowplot");
  library("cowplot");
} 
```
```
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
```
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
```
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
```
plot +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15);
```
表示"A"在x=0 y=1的位置，"B"在x=0 y=0.5的位置，"C"在x=0.5 y=0.5的位置
![cowplot包3](./md-image/cowplot包3.png){:width=400 height=400}
