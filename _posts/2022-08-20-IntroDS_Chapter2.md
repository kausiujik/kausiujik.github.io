---
layout: post
title: "第二章——数据预处理"
subtitle: "Data Preprocessing"
author: "Watthu"
header-style: text
mathjax: true
tags:
  - 数据科学导引
---

## 2.1 特征编码

在处理数据分析问题时, 针对数值型的数据, 我们利用相关数据分析算法可以很轻松地实现对其的分析; 但在某些问题中, 获取的原始数据会包含非数值型的特征, 而且它们往往是离散型特征, 这就需要我们对其进行**特征编码(feature embedding)**.

### 2.1.1 数字编码

**数字编码**是指从0开始赋予特征的每一个取值一个整数.

对于**等级型特征**, 按照特征取值从小到大进行整数编码可以保证编码后的数据保留原有次序关系. 

对于**名义型特征**, 按上述方法进行数字编码会人为引入对特征的主观偏好, 对数据分析具有误导性, 往往采用下面所讲的**独热(One-Hot)编码**.

### 2.2.2 独热(One-Hot)编码

**独热(One-Hot)编码**是指将包含 $K$ 个取值的离散型特征转换成 $K$ 个二元特征(取值为0或1的特征).

课本中给了一个例子, "汽车品牌"={路虎, 吉利, 奥迪, 大众, 奔驰}, 包含5个不同的值, 将其编码为5个特征 $f_1$, $f_2$, $f_3$, $f_4$, $f_5$, 与相应的品牌一一对应, 转换后的特征的取值如下表.

| 原始特征取值 |  $f_1$  |  $f_2$  |  $f_3$  |  $f_4$  |  $f_5$  |
| :-: | :-: | :-: | :-: | :-: | :-: |
| 路虎 | 1 | 0 | 0 | 0 | 0 |
| 吉利 | 0 | 1 | 0 | 0 | 0 |
| 奥迪 | 0 | 0 | 1 | 0 | 0 |
| 大众 | 0 | 0 | 0 | 1 | 0 |
| 奔驰 | 0 | 0 | 0 | 0 | 1 |

在线性回归模型中, 对名义型特征进行One-Hot编码的效果通常比数字编码的效果要好; One-Hot编码对包含离散型特征的分类模型的效果有很好的提升. 这也是这种编码被广泛使用的原因之一.

但是, One-Hot编码也有缺点, 一个最明显的劣势是使特征维度显著增多, 另外是增加了特征之间的相关性, 比如按上表方式编码后的5个特征存在如下线性关系: 

<center>
$f_1+f_2+f_3+f_4+f_5=1.$
</center>

这会影响线性回归等模型的效果.

### 2.2.3 哑变量编码

**哑变量编码**是指将包含 $K$ 个取值的离散型特征转换成 $K-1$ 个二元特征.

这种编码方式解决了独热编码的相关性问题, 转换后的特征的取值如下表.

| 原始特征取值 |  $f_1$  |  $f_2$  |  $f_3$  |  $f_4$  | 
| :-: | :-: | :-: | :-: | :-: |
| 路虎 | 1 | 0 | 0 | 0 |
| 吉利 | 0 | 1 | 0 | 0 |
| 奥迪 | 0 | 0 | 1 | 0 |
| 大众 | 0 | 0 | 0 | 1 |
| 奔驰 | 0 | 0 | 0 | 0 |

## 2.2 缺失值处理

**数据缺失**是指在数据采集、传输和处理等过程中, 由于某些原因导致数据不完整的情况.

### 2.2.1 删除法

**删除法**通过删除包含缺失值的数据, 来得到一个完整的数据子集. 数据的删除可以从以下两个角度进行:

 - 删除样本: 删除存在数据缺失的样本. 适用于某些样本有多个特征存在缺失值, 且存在缺失值的样本占整个数据集样本数量的比例不高的情形.
 - 删除特征: 当某个特征缺失值较多时, 且该特征对数据分析的目标影响不大时, 可以将该特征剔除.
 
删除法以减少数据来换取信息的完整, 丢失了大量隐藏在这些被删除数据中的信息. 在一些实际场景下, 数据的采集成本高且缺失值无法避免, 删除法可能会造成大量的资源浪费.

### 2.2.2 均值填补

对于存在缺失值的某一个特征, **均值填补法**首先计算该特征中非缺失值的平均值或众数, 然后使用平均值或众数来代替缺失值. 对于连续型特征, 通常使用平均值进行填补; 对于离散型特征, 通常使用众数进行填补.

均值填补法会使得数据过分集中在平均值或众数上, 导致特征的方差被低估. 此外, 由于完全忽略特征之间的相关性, 均值填补法会大大弱化特征之间的相关性. 

在实际应用中, 可以根据一定的辅助特征, 将数据集分成多组, 然后在每一组数据上分别使用均值填补.

### 2.2.3 随机填补

**随机填补**是在均值填补的基础上加上随机项, 通过增加缺失值的随机性来改善缺失值分布过于集中的缺陷. 常见的方法包括 <u>贝叶斯Bootstrap方法</u> 和 <u>近似贝叶斯Bootstrap方法</u>.

 - 贝叶斯Bootstrap方法:
   - Step1: 从均匀分布 $U(0,1)$ 中随机抽取 $k-1$ 个随机数, 并进行升序排序记为 $a_1,a_2,\dots,a_{k-1}$ ;
   - Step2: 对 $(n-k)$ 个缺失值, 分别从非缺失值 $\\{f_1,f_2,\dots,f_k\\}$ 中以概率 $\\{a_1,a_2,\dots,a_{k-1},1-a_{k-1}\\}$ 采样一个值进行填补.

 - 近似贝叶斯Bootstrap方法:
   - Step1: 从 $k$ 个非缺失值 $\\{f_1,f_2,\dots,f_k\\}$ 中有放回地抽取 $k$ 个值建立一个新的大小为 $k$ 的集合 $F$ ;
   - Step2: 对 $(n-k)$ 个缺失值, 分别从 $F$ 中随机抽取一个值进行填补.

### 2.2.4 基于模型的填补

**基于模型的方法**将缺失特征 $f$ 作为预测目标. 将数据集中其他特征或其子集作为输入特征, 通过特征 $f$ 的非缺失值构造训练集, 训练分类或回归模型, 然后使用构建的模型来预测特征 $f$ 的缺失值.

在实际使用时, 需要采用模型评估方法对模型的预测性能进行评估. 此外, 这种方法会增大特征之间的相关性.

### 2.2.5 其他缺失值处理方法

除了上述处理方法外, 缺失值处理的方法还包括**哑变量方法**和**EM算法**等.

 - 哑变量方法: 对于离散型特征, 如果存在缺失值, 可以将缺失值作为一个单独的取值处理.
 - EM算法: 利用不完整的信息实现概率模型的参数化估计.

## 2.3 数据标准化

在不少学习算法中, 目标函数往往假设其特征均值在0附近且方差齐次, 这就要求数据标准化. 而在另外一些数据分析场景下, 我们需要计算样本之间的相似度, 如果特征之间的量纲差异太大, 相似度评估结果也随之受到影响.

### 2.3.1 Z-score标准化

**Z-score标准化**通过对数据集中的每一个样本进行处理, 使得处理后的数据具有固定均值和标准差.

假设特征 $f$ 的取值集合为 $\\{f_1,f_2,\dots,f_n\\}$ , 则特征取值 $f_i$ 经过Z-score标准化后的取值 $f_i'$ 为

<center>
$f_i'=\dfrac{f_i-\mu}{\sigma},$
</center>

其中 $\mu=\dfrac{1}{n}\sum\limits_{i=1}^nf_i$ 为特征 $f$ 的平均值, $\sigma=\sqrt{\dfrac{1}{n}\sum\limits_{i=1}^n(f_i-\mu)^2}$ 为特征 $f$ 的标准差.

经过Z-score标准化后的特征能够直观反映每一个取值距离平均值的标准差距离, 从而理解特征的整体分布情况.

当数据中存在离群值(outlier)时, 为了降低离群值的影响, 可以将Z-score标准化方法中的标准差替换成平均绝对偏差 $s=\dfrac{1}{n}\sum\limits_{i=1}^n\|f_i-\mu\|$.

Z-score标准化适用于特征的最大值或最小值未知和样本分布非常分散的情况.

### 2.3.2 Min-Max标准化

**Min-Max标准化**通过对特征作线性变换, 使得转换后特征的取值分布在 $[0,1]$ 区间内.

假设特征 $f$ 的取值集合为 $\\{f_1,f_2,\dots,f_n\\}$, 则特征取值 $f_i$ 经过Min-Max标准化后的取值 $f_i'$ 为

<center>
$f_i'=\dfrac{f_i-f_{\min}}{f_{\max}-f_{\min}},$
</center>

其中 $f_{\min}$ 为特征 $f$ 取值的最小值, $f_{\max}$ 为特征 $f$ 取值的最大值.

进一步地, 如果希望将特征 $f$ 线性映射到任意区间 $[a,b]$, 则Min-Max标准化的方法为

<center>
$f_i'=\dfrac{b-a}{f_{\max}-f_{\min}}(f_i-f_{\min})+a.$
</center>

Min-Max标准化适用于需要将特征取值简单地线性映射到某一区间中的情形, 其不足之处在于当数据集中有新数据加入时, 特征的最大值或最小值会发生变化, 从而需要重新计算全部数据特征的标准化取值. 此外, 数据集中存在离群值时, Min-Max标准化后的效果也较差.

### 2.3.3 小数定标标准化

**小数定标(demical scaling)标准化**通过移动特征取值的小数点位置来进行标准化, 使得标准化后特征取值的绝对值总是小于1. 小数点移动的位数取决于特征取值的最大绝对值大小.

假设特征 $f$ 的取值集合为 $\\{f_1,f_2,\dots,f_n\\}$, 则特征取值 $f_i$ 经过小数定标标准化后的取值 $f_i'$ 为

<center>
$f_i'=\dfrac{f_i}{10^j},$
</center>

其中 $j$ 是满足 $\max\\{f_1',f_2',\dots,f_n'\\}<1$ 的最小整数.

小数定标标准化适用于特征取值比较分散, 尤其是特征取值分布在多个数量级的情况. 但这种方法也存在许多缺点, 如果特征取值分布集中在某几个量级上, 则小数定标标准化的特征取值也会集中在某几个值附近, 不利于后续数据分析时的样本区分. 此外, 当数据集中有新数据加入时, 需要重新确定小数点移动位数; 小数定标标准化的效果也会受到离群值的影响.

### 2.3.4 逻辑斯蒂标准化

**逻辑斯蒂(Logistic)标准化**利用逻辑斯蒂函数的特性, 将特征取值映射到 $[0,1]$ 区间内, 其中逻辑斯蒂函数为

<center>
$\sigma(x)=\dfrac{1}{1+\text{e}^{-x}}.$
</center>

它将数据从实数域光滑映射到 $[0,1]$ 区间.

假设特征 $f$ 的取值集合为 $\\{f_1,f_2,\dots,f_n\\}$, 则特征取值 $f_i$ 经过逻辑斯蒂标准化后的取值 $f_i'$ 为

<center>
$f_i'=\dfrac{1}{1+\text{e}^{-f_i}}.$
</center>

逻辑斯蒂标准化适用于特征取值分布相对比较集中地分布于0两侧的情况. 如果特征取值分散且均远离0, 那么标准化后的特征取值会集中于0或1附近, 造成原始特征的分布及取值间关系被改变. 因此在使用逻辑斯蒂标准化之前, 需要首先分布原始特征取值的分布情况.

### 2.3.5 不同标准化方法的对比

标准化处理会改变原始特征的取值, 因此在实际操作中需要保存所使用的标准化方法的参数, 以便对后续的数据进行统一的标准化处理.

 - Z-score标准化: 无需计算特征的最大值和最小值, 需要记录特征的均值和标准差.
 - Min-Max标准化: 需要保留原始数据间的关系.
 - 小数定标标准化: 数据分布离散且存在多个数量级.
 - 逻辑斯蒂标准化: 对分布离散并且远离零点的数据处理效果不佳.

## 2.4 特征离散化

**特征离散化(data discretization)**是指将连续型特征转换成离散型特征的过程. 特征的离散化过程是将连续型特征的取值范围划分成若干区间段(bin), 然后使用区间段代替落在该区间段的特征取值. 区间段之间的分割点称为切分点(cut point), 由切分点分割出来的子区间段的个数, 称为元数(arity).

如何在数据信息损失尽量少的前提下, 尽可能减少元数, 是特征离散化要追求的目标. 一般来说, 特征离散化方法分为有监督和无监督两种.

 - **无监督离散化**: 不参考预测特征 $y$, 直接根据特征本身的分布特性进行离散化处理. 常见的包括等距离散化、等频离散化和基于聚类的离散化等.
 - **有监督离散化**: 利用参考数据集中的预测特征 $y$, 将连续型特征进行离散化处理. 常见的包括信息增益离散化、卡方离散化和CAIM离散化等.

特征离散化方法一般分为下面四步进行:

 1. 特征排序: 对连续型特征的取值进行升序或降序排列, 以减少离散化的运算开销.
 2. 切分点选择: 根据给定的评价准则, 合理选择切分点. 常用的有基于信息增益或基于统计量的准则.
 3. 区间段分割或者合并: 基于选择好的切分点, 对现有的区间段进行分割或者合并, 得到新的区间段.
 4. 在生成的新区间段上重复第1～3步, 直到满足终止条件. 我们可以预先设定元数 $k$, 作为简单的终止判断标准, 也可以设定复杂的判断函数.

### 2.4.1 等距离散化

**等距离散化(equal width discretization)**根据连续型特征的取值, 将其均匀地划分成$k$个区间, 每个区间的宽度均相等, 然后将特征的取值划入对应的区间从而完成特征离散化.

等距离散化是一种简单的划分连续型特征的方法, 它对输入的数据质量要求高, 无法解决特征存在离群值的问题.

### 2.4.2 等频离散化

与等距离散化不同, **等频离散化**尽量使得离散化后的每一个区间内的样本量均衡, 以避免出现区间段中的样本量出现严重的不均衡.

等距离散化和等频离散化是最简单易用的两种无监督特征离散化方法.

 - 等距离散化: 倾向于将数据不均匀地划分入离散化区间中, 对离群值较为敏感.
 - 等频离散化: 将数据变换为均匀分布, 有时会将同样或接近的样本划分入不同的区间, 容易使得相邻区间段内的数据具有相似的特性.

### 2.4.3 聚类离散化

**聚类离散化**将样本划分到不同的类或簇(cluster)中, 聚类的结果是同一个簇中的样本有很大的相似性, 不同簇间的样本有很大的差异性.

基于聚类分析的离散化方法主要包括以下三个步骤:

 1. 对于需要离散化的连续型特征, 采用聚类算法(如K-means、EM算法等), 把样本依据该特征的分布划分成相应的簇或类.
 2. 在聚类结果的基础上, 基于特定的策略, 决定是否对簇进行进一步分裂或合并.
 3. 在最终确定划分的簇之后, 确定切分点以及区间个数.

在整个聚类的过程中, 需要事先确定簇的个数以及描述样本之间的距离计算方式. 如何选定簇的个数也会影响聚类算法的效果, 从而影响特征的离散化.

### 2.4.4 信息增益离散化

**信息增益离散化**灵感源自决策树模型建立时基于信息增益的评价标准, 它具体分为以下几个步骤:

 1. 对连续型特征进行排序.
 2. 把特征的每一个取值认为是候选分裂节点(切分点), 计算出对应的熵. 随后, 选择熵最小的取值作为正式的切分点, 将原来的区间一分为二.
 3. 递归处理第2步中得到的两个新区间段, 直到每个区间段内特征的类别一样为止.
 4. 合并相邻的、类的熵值为0且特征类别相同的区间段, 并重新计算新区间段类的熵值.
 5. 重复执行第4步, 直到满足终止条件. 终止条件可以是限制决策树的深度, 或者叶结点的个数等.
 
在众多的决策树算法中, ID3和C4.5是最常用的基于信息增益来进行特征选择分类的算法. 将这两种决策树算法应用于连续型特征离散化的核心是针对单个特征构建决策树建模, 然后根据决策树模型中结点分裂的阈值对特征进行划分.

### 2.4.5 卡方离散化

信息增益离散化是一种自顶向下的分裂策略, 而**卡方离散化**是一种自底向上的合并策略. 最常用的基于卡方检验的离散化方法叫做ChiMerge方法, 它通过卡方检验来判断相邻区间是否需要进行合并, 这要求区间内特征取值的类别要独立于区间. ChiMerge离散化的过程如下:

 1. 将连续型特征的每个取值看作是一个单独的区间段, 并对值进行排序.
 2. 针对每对相邻的区间段, 计算卡方统计量的值, 将卡方值最小或者低于设定阈值的相邻区间段合并在一起. 这里卡方统计量为
 <center>
 $\displaystyle\chi^2=\sum\limits_{i=1}^k\sum\limits_{j=1}^C\dfrac{(A_{ij}-E_{ij})^2}{E_{ij}},$
 </center>
 其中 $A_{ij}$ 为第 $i$ 区间段内类别为 $j$ 的样本个数, $k$ 为比较的区间个数, $C$ 为类别个数, 且有 $E_{ij}=\dfrac{1}{n}\sum\limits_{q=1}^CA_{iq}\sum\limits_{p=1}^kA_{pj}$ 为 $A_{ij}$ 的期望频数. 
 3. 对于新的区间段, 递归进行1和2两步, 直到满足终止条件.
 
ChiMerge离散化在每次迭代的过程中只能合并两个区间($k=2$), 最终是通过阈值的设定来控制元数的大小. 如果数据集中样本量较大, 算法的开销会较大.

### 2.4.6 类别属性相互依赖最大化

**类别属性相互依赖最大化(Class-attribute Interdependence Maximization, CAIM)**也是一种基于熵的特征离散化方法, 它采取自顶向下的策略, 通过选择切分点 $p$, 把特征的取值空间划分为 $f\leq p$ 和 $f>p$ 两个子区间段, 其中衡量切分点优劣的度量方法是类别属性的相互依赖程度.

假设某个连续型特征有 $n$ 个取值, $C$ 个类别. 现将特征划分为 $k$ 个子区间段, 其集合记为 $D=\\{[d_0,d_1],(d_1,d_2],\dots,(d_{k-1},d_k]\\}$, 其中 $d_0$ 和 $d_k$ 分别为特征的最小值和最大值. 设 $n_{i\cdot}$ 表示属于类别 $i$ 的样本个数, $n_{\cdot j}$ 表示落在区间段 $(d_{j-1},d_j]$ 的样本个数, $n_{ij}$ 表示在区间段 $(d_{j-1},d_j]$ 内的且属于类别 $i$ 的样本个数. 可以计算一个评价准则来评价当前离散化的好坏:

<center>
$\displaystyle\text{CAIM}=\dfrac{1}{N}\sum\limits_{j=1}^k\dfrac{M_j^2}{n_{\cdot j}},$
</center>

其中 $M_j=\max\\{n_{1j},n_{2j},\dots,n_{Cj}\\}$. CAIM取值区间为 $(0,1]$, 值越大说明类和离散区间的相互依赖程度越大, 也就说明当前离散化效果较好.

CAIM离散化分为以下几个步骤:

 1. 对进行离散化的特征进行升序排序, 确定取值区间的最小值 $d_0$ 和最大值 $d_k$, 初始化划分策略 $D=\\{[d_0,d_k]\\}$.
 2. 把区间内的每个值当作候选切分点, 计算把区间二分后的CAIM值, 并选取CAIM值最高的点作为切分点, 并更新 $D$.
 3. 对于 $D$ 中的每个区间段, 重复第2步的计算, 直到满足终止条件.
 
CAIM只关注区间内拥有最多样本的类别与特征之间的关系, 忽略同一区间内其他类别的信息. CAIM最终生成的离散化区间个数往往与样本的类别个数接近.

### 2.4.7 小结

 - 离散化本质是将连续型数据分段, 因此数据中的离群值会直接划入相应的区间段中, 进而增强模型对于数据离群值的鲁棒性.
 - 离散化后的特征, 其取值均转化为有明确含义的区间号, 相对于原始的连续型数据来说, 含义更加明确, 从而使得数据的可解释性更强, 模型更易于理解和应用.
 - 将连续化特征离散化后, 特征的取值大大减少, 既减少了数据集对于存储空间的需求, 又能减少模型训练的实际计算量, 从而可以提升模型训练的计算效率.
 
更多离散化方法描述及其在有监督学习任务中的详细综述和对比: [Garcia S, Luengo J, Sáez J A, et al. A survey of discretization techniques: Taxonomy and empirical analysis in supervised learning[J]. IEEE transactions on Knowledge and Data Engineering, 2012, 25(4): 734-750.](https://ieeexplore.ieee.org/abstract/document/6152258)

## 2.5 离群值检测

**离群值(outlier)**也称为异常值, 是指一个数据集中那些明显偏离数据集中其他样本的值. 离群值的出现主要是因为自然变异、数据测量和收集的误差以及人工操作失误等. 离群值检测可以作为数据预处理的一个步骤, 通过它来发现数据中的噪音, 为数据分析提供高质量的数据.

### 2.5.1 基于统计的方法

在统计学中, 比较常用的检测方法是**拉依达准则(Pauta criterion)**. 假定数据集服从正态分布, 并对数据集计算标准差. 在给定置信度的情况下, 如果有超过均值三倍标准差以上, 我们认为含有该误差的样本即为离群值.

拉依达准则只能在测量次数较多时使用. Z-score标准化可以用来辅助识别离群值, 通常将满足Z-score绝对值大于3的样本认为是离群值.

### 2.5.2 基于近邻的方法

如果一个样本是离群值, 那么它离其他样本的距离会比较远. 由此可以对样本的离群程度进行量化, 分数由它与 $k$ 个最近邻的距离决定, 分数的取值区间为 $[0,+\infty)$. 基于K-means算法的离群值检测方法分为以下几个步骤:

 1. 计算每一个样本与其最近的 $k$ 个近邻样本的距离, 放入到集合 $C$ 中.
 2. 对 $C$ 中的所有元素进行降序排序.
 3. 根据给定的距离阈值, 选取 $C$ 中大于给定阈值的距离所对应的样本作为离群值.
 
基于K-means算法的离群值检测简单有效, 比较常用, 但它只能发现全局的离群值, 无法发现局部离群值.

下面介绍一种基于局部密度的离群值检测方法, 称作**局部离群因子算法(Local Outlier Factor, LOF)**. 该算法给数据集中每一个样本计算一个局部离群分数, 它通过样本邻居的密度和当前样本的密度相对值来计算, 最终根据分数大小来识别离群样本.

首先, 记 $d(x_1,x_2)$ 为样本 $x_1$ 和样本 $x_2$ 之间的距离, $d_k(x)$ 为样本 $x$ 到其第 $k$ 个近邻样本的距离, $N_k(x)$ 表示样本 $x$ 的 $k$ 个近邻样本的集合. 定义样本 $x_1$ 到样本 $x_2$ 的**可达距离(reachability distance)**为

<center>
$\displaystyle\text{rd}_k(x_1,x_2)=\max\{d_k(x_2),d(x_1,x_2)\}.$
</center>

基于可达距离, 进一步定义**局部可达密度(local reachability density)**, 它是样本与其近邻样本的平均可达距离的倒数:

<center>
$\displaystyle\text{lrd}_k(x)=\left(\dfrac{1}{k}\sum\limits_{y\in N_k(x)}\text{rd}_k(x,y)\right)^{-1}.$
</center>

利用样本及其近邻样本的局部可达密度, 可以计算样本的局部离群因子:

<center>
$\displaystyle\text{lof}_k(x)=\dfrac{1}{k}\sum\limits_{y\in N_k(x)}\dfrac{\text{lrd}_k(y)}{\text{lrd}_k(x)}.$
</center>

进一步我们可以根据以下规则对样本是否离群做判断: 如果 $\text{lof}_k(x)$ 接近于1, 代表样本 $x$ 与其近邻样本的局部可达密度相似, 可以认为是正常样本; 如果 $\text{lof}_k(x)<1$, 代表样本 $x$ 的局部可达密度大于其近邻样本, 可以认为是正常样本; 如果 $\text{lof}_k(x)>1$, 代表样本 $x$ 的局部可达密度小于其近邻样本, 可能为离群值.

无论是K-means方法还是LOF算法, 都需要计算数据集中样本之间的距离. 当样本数量很大且维度很大时, 这两种方法的计算代价都很大.

### 2.5.3 小结

当我们从数据集中发现离群值后, 不能简单地认为是噪音而删除, 而是要结合具体的业务需求去考察离群值存在的合理性以及离群值产生的原因, 从而能够调整现有的模型, 提高模型的解释能力.

## 2.6 其他预处理方法

数据预处理是数据分析的重要前置步骤, 其结果直接决定了数据分析的效果. 下面再介绍一些其他的预处理方法.

 - 中心化(centering): 让数据集中特征的均值为0. 例如在部分降维算法(如PCA)运行前就要对数据进行中心化转换.
 - 分布转换: 很多模型都假设特征服从正态分布, 如果实际的特征不服从正态分布可能需要进行转换.
 
在实际的大数据项目中, 可能还需要对不同类型和来源的数据进行整合, 可能出现例如数据不一致、数据冗余和数据冲突等情况.