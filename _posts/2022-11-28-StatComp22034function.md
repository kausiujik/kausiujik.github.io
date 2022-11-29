---
layout: post
title: "StatCompWithR: 统计计算函数细节"
subtitle: "Statistical Computing with R (details about functions)"
author: "Watthu"
header-style: text
mathjax: true
tags:
  - 统计计算
---

本文仅用于阐明《统计计算》课程`StatComp22034` R包中自制函数的背景知识与相关理论推导。

## Individualized Regression For Clustering

### 模型概述

记 $y_i\in\mathbb R$ 为第 $i$ 个个体的响应变量（response）观测值，$x_i\in\mathbb R$ 为第 $i$ 个个体异质性的协变量（covariate）观测值，$z_i\in\mathbb R$ 为第 $i$ 个个体的群体同质性的协变量观测值。考虑一般的回归模型
<center>
$y_i=x_i\beta_i+z_i\alpha+\varepsilon_i,\quad i=1,\dots,N.$
</center>

其中 $\beta_i\in\mathbb R$ ($i=1,\dots,N$) 就衡量了不同个体的异质性，而 $\alpha\in\mathbb R$ 则是群体共享参数，$\varepsilon_i$ 是满足零均值，同方差的随机误差。

记 $\pmb\beta=(\beta_1,\dots,\beta_N)'$，$\pmb\gamma=(\gamma_1,\dots,\gamma_K)'$。定义
<center>
$D(\pmb\beta,\pmb\gamma)=\left\|\big(\min\{d(\beta_1,\gamma_1),\dots,d(\beta_1,\gamma_K)\},\dots,\min\{d(\beta_N,\gamma_1),\dots,d(\beta_N,\gamma_K)\}\big)'\right\|_0.$
</center>

为了实现对个体的聚类，考虑目标函数为
<center>
$Q(\pmb\beta,\alpha,\pmb\gamma)=\dfrac{1}{2}\left\|\pmb y-\pmb X\pmb\beta-\alpha\pmb z\right\|_2^2+\lambda D(\pmb\beta,\pmb\gamma),$
</center>

其中 $\pmb y=(y_1,\dots,y_N)'$，$\pmb X=\text{diag}(x_1,\dots,x_N)$ 是 $N\times N$ 矩阵，$\pmb z=(z_1,\dots,z_N)'$，$\lambda>0$ 是惩罚系数。

### 参数优化

进一步, 可以使用ADMM算法解决这一问题. 改写优化问题为
<center>
$\begin{array}{rl}
\min\limits_{\pmb\beta,\alpha,\pmb\delta,\pmb\gamma} & \dfrac{1}{2}\|\pmb y-\pmb X\pmb\beta-\alpha\pmb z\|_2^2+\lambda D(\pmb\delta,\pmb\gamma),\\
\text{s.t.}~ &\pmb\beta-\pmb\delta=\pmb 0.
\end{array}$
</center>

则该问题的增广拉格朗日函数为
<center>
$L_\rho(\pmb\beta,\alpha,\pmb\delta,\pmb\gamma;\pmb\nu)=\dfrac{1}{2}\|\pmb y-\pmb X\pmb\beta-\alpha\pmb z\|_2^2+\lambda D(\pmb\delta,\pmb\gamma)+\dfrac{\rho}{2}\left\|\pmb\beta-\pmb\delta+\dfrac{1}{\rho}\pmb\nu\right\|_2^2.$
</center>

其中 $\rho>0$ 是事先设定的固定参数，$\pmb\nu\in\mathbb R^N$是拉格朗日乘子向量.

ADMM算法的好处在于把这个函数分成好几个容易求解的部分进行计算。具体来说，

1. 设定初值 $\pmb\beta^{(0)}$，$\alpha^{(0)}$，再设定初值 $\pmb\delta^{(0)}$，$\pmb\gamma^{(0)}$ 和 $\pmb\nu^{(0)}$。
2. 对第 $r+1$（$r=0,1,\dots$）次迭代，按下述式子更新诸参数的值：

    <center>
    $\begin{array}{rl}
    \left(\begin{array}{c}
    \pmb\beta^{(r+1)}\\ \alpha^{(r+1)}
    \end{array}\right)=&\!\!\!\!\underset{\pmb\beta,\alpha}{\arg\min} ~L_\rho(\pmb\beta,\pmb\alpha,\pmb\delta^{(r)},\pmb\gamma^{(r)};\pmb\nu^{(r)}),\\
    \left(\begin{array}{c}
    \pmb\delta^{(r+1)}\\ \pmb\gamma^{(r+1)}
    \end{array}\right)=&\!\!\!\!\underset{\pmb\delta,\pmb\gamma}{\arg\min} ~L_\rho(\pmb\beta^{(r+1)},\alpha^{(r+1)},\pmb\delta,\pmb\gamma;\pmb\nu^{(r)}),\\
    \pmb\nu^{(r+1)}=&\!\!\!\!\pmb\nu^{(r)}+\rho(\pmb\beta^{(r+1)}-\pmb\delta^{(r+1)}).
    \end{array}$
    </center>
    
3. 重复迭代直至算法收敛或迭代次数过多。

    <center>
    $\dfrac{1}{N}\|\pmb\beta^{(r+1)}-\pmb\beta^{(r)}\|_2+\|\alpha^{(r+1)}-\alpha^{(r)}\|_2+\dfrac{1}{K}\|\pmb\gamma^{(r+1)}-\pmb\gamma^{(r)}\|_2<\varepsilon_1$ 或 $\|\pmb r^{(r+1)}-\pmb r^{(r)}\|_2<\varepsilon_2$，
    </center>
    
    其中 $\pmb r^{(r)}=\pmb\beta^{(r)}-\pmb\delta^{(r)}$。

#### $\pmb\beta$ 和 $\alpha$ 的更新

舍去与最小化函数无关的项，即
<center>
$\left(\begin{array}{c}
\pmb\beta^{(r+1)}\\ \alpha^{(r+1)}
\end{array}\right)=\underset{\pmb\beta,\alpha}{\arg\min}~\dfrac{1}{2}\|\pmb y-\pmb X\pmb\beta-\alpha\pmb z\|_2^2+\dfrac{\rho}{2}\left\|\pmb\beta-\pmb\delta+\dfrac{1}{\rho}\pmb\nu\right\|_2^2.$
</center>

即求解方程
<center>
$\left\{\begin{array}{l}
-\pmb X'(\pmb y-\pmb X\pmb\beta^{(r+1)}-\alpha^{(r+1)}\pmb z)+\rho(\pmb\beta^{(r+1)}-\pmb\delta^{(r)})+\pmb\nu^{(r)}=\pmb 0,\\
-\pmb Z'(\pmb y-\pmb X\pmb\beta^{(r+1)}-\alpha^{(r+1)}\pmb z)=\pmb 0.
\end{array}\right.$
</center>

该方程有显式解
<center>
$\left(\begin{array}{c}
\pmb\beta^{(r+1)}\\ \alpha^{(r+1)}
\end{array}\right)=\left[\widetilde{\pmb X}'\widetilde{\pmb X}+\left(\begin{array}{cc}
\rho\pmb I_N & \pmb 0_N\\
\pmb 0_N' & 0
\end{array}\right)\right]^{-1}\left[\widetilde{\pmb X}'\pmb y+\left(\begin{array}{c}
\rho\pmb\delta^{(r)}-\pmb\nu^{(r)}\\
0
\end{array}\right)\right],$
</center>

其中 $\widetilde{\pmb X}=(\pmb X,\pmb z)$ 是 $N\times(N+1)$ 矩阵，$\pmb I_N$ 是 $N\times N$ 单位矩阵。

#### $\pmb\delta$ 和 $\pmb\gamma$ 的更新

舍去与最小化函数无关的项，即
<center>
$\left(\begin{array}{c}
\pmb\delta^{(r+1)}\\ \pmb\gamma^{(r+1)}
\end{array}\right)=\underset{\pmb\delta,\pmb\gamma}{\arg\min}~\lambda D(\pmb\delta,\pmb\gamma)+\dfrac{\rho}{2}\left\|\pmb\beta-\pmb\delta+\dfrac{1}{\rho}\pmb\nu\right\|_2^2.$
</center>

在给定 $\pmb\gamma$ 的条件下，注意到两部分对 $\pmb\delta$ 的计算是分离的，因此计算可知
<center>
$\delta_i^{(r+1)}=\left\{\begin{array}{ll}
\gamma_{k_0}, & \exists k_0\in\{1,\dots,K\},\left\vert\gamma_{k_0}-\beta_i^{(r+1)}-\rho^{-1}\nu_i^{(r)}\right\vert\leq\dfrac{2\lambda}{\rho},\\
\beta_i^{(r+1)}+\rho^{-1}\nu_i^{(r)}, & \forall k\in\{1,\dots,K\},\left\vert\gamma_{k}-\beta_i^{(r+1)}-\rho^{-1}\nu_i^{(r)}\right\vert>\dfrac{2\lambda}{\rho}.
\end{array}\right.$
</center>

注意这里若有多个 $k_0$ 满足条件，则取使得不等号左边最小的那个 $k_0$。

而当 $\pmb\delta$ 完成更新后，基于最小化 $D(\pmb\delta,\pmb\gamma)$ 可知，$\pmb\gamma^{(r+1)}$ 为 $\pmb\delta^{(r+1)}$ 的 $K$-mediods 结果。

## EM Algorithm for Repeated Experiments

### 模型概述

记 $y_{ij}$ 为第 $i$ $(i=1,\dots,n)$ 个个体在第 $j$ $(j=1,\dots,t)$ 次重复实验中的反应时间，构建贝叶斯层次模型如下：
<center>
$\begin{array}{rl}
	y_{ij}\mid\alpha_i,z_{ij},\beta,\sigma_1^2,\sigma_2^2\sim&\!\!\!\! N(\alpha_i+\beta z_{ij},(1-z_{ij})\sigma_1^2+z_{ij}\sigma_2^2),\\
	\alpha_i\mid s_i,\mu,\gamma,\tau_1^2,\tau_2^2\sim&\!\!\!\! N(\mu+\gamma s_i,(1-s_i)\tau_1^2+s_i\tau_2^2),\\
	z_{ij}\mid\lambda,s_i\sim&\!\!\!\!\text{Bernoulli}(\lambda s_i),
\end{array}$
</center>

其中 $s_i=1$ 当且仅当第 $i$ 个个体患有精神分裂症，这些个体的平均反应时间比正常个体有 $\gamma>0$ 的偏移；$z_{ij}=1$ 当且仅当第 $i$ 个个体在第 $j$ 次重复实验中病症发作，发病率为 $\lambda$，从而平均以 $\beta>0$ 的弛豫时间才作出反应；$\alpha_i$ 是每个个体的平均反应时间，$\sigma_1^2$，$\sigma_2^2$ 分别是正常个体和病症发作个体反应时间的方差，$\tau_1^2$，$\tau_2^2$ 分别是正常个体和患病个体平均反应时间的方差.

记 $\pmb\theta=(\pmb\alpha',\beta,\gamma,\mu,\lambda,\sigma_1^2,\sigma_2^2,\tau_1^2,\tau_2^2)$ 为模型的未知参数，其中 $\pmb\alpha=(\alpha_1,\dots,\alpha_n)'$，而 $z_{ij}$ 为隐变量。于是模型在完全数据下的似然函数为
<center>
$\begin{array}{rl}
L_c(\pmb\theta\mid\pmb y,\pmb s,\pmb z)&\propto\displaystyle\prod_{i=1}^n\prod_{j=1}^t\left[f(y_{ij};\alpha_i,\sigma_1^2)\right]^{1-z_{ij}}\left[f(y_{ij};\alpha_i+\beta,\sigma_2^2)\right]^{z_{ij}}\\
&\displaystyle\quad\cdot\prod_{i=1}^n\left[f(\alpha_i;\mu,\tau_1^2)\right]^{1-s_i}\left[f(\alpha_i;\mu+\gamma,\tau_2^2)\right]^{s_i}\cdot\prod_{i=1}^n\prod_{j=1}^t\left[\lambda^{z_{ij}}(1-\lambda)^{1-z_{ij}}\right]^{s_i},
\end{array}$
</center>

其中 $f(\cdot\mid\mu,\sigma^2)$ 是正态分布 $N(\mu,\sigma^2)$ 的密度函数.

### EM算法

考虑到模型中存在缺失数据 $\pmb z$，因此使用E-M算法如下。首先计算模型在完全数据下的对数似然函数
<center>
$\begin{array}{rl}
l_c(\pmb\theta\mid\pmb y,\pmb s,\pmb z)=&\!\!\!\!\displaystyle\sum_{i=1}^n\sum_{j=1}^t\left\{(z_{ij}-1)\left[\frac{\ln\sigma_1^2}{2}+\dfrac{(y_{ij}-\alpha_i)^2}{2\sigma_1^2}\right]-z_{ij}\left[\frac{\ln\sigma_2^2}{2}+\frac{(y_{ij}-\alpha_i-\beta)^2}{2\sigma_2^2}\right]\right\}\\
&\!\!\!\!\displaystyle\quad+\sum_{i=1}^n\left\{(s_i-1)\left[\frac{\ln\tau_1^2}{2}+\dfrac{(\alpha_i-\mu)^2}{2\tau_1^2}\right]-s_i\left[\frac{\ln\tau_2^2}{2}+\dfrac{(\alpha_i-\mu-\gamma)^2}{2\tau_2^2}\right]\right\}\\
&\!\!\!\!\displaystyle\quad+\sum_{i=1}^n\sum_{j=1}^ts_i\left[z_{ij}\ln\lambda+(1-z_{ij})\ln(1-\lambda)\right]+Const,
\end{array}$
</center>

#### E步

注意到 $z_{ij}\mid\lambda,s_i\sim\text{Bernoulli}(\lambda s_i)$，因此
<center>
$\begin{array}{rl}
Q(\pmb\theta\mid\pmb\theta^{(r)})=&\!\!\!\!\displaystyle\mathbb E_{\pmb\theta^{(r)}}[l_c(\pmb\theta\mid\pmb y,\pmb s,\pmb z)\mid\pmb y,\pmb s]\\
=&\!\!\!\!\displaystyle\sum_{i=1}^n\sum_{j=1}^t\left\{(\xi_{ij}^{(r)}-1)\left[\frac{\ln\sigma_1^2}{2}+\dfrac{(y_{ij}-\alpha_i)^2}{2\sigma_1^2}\right]-\xi_{ij}^{(r)}\left[\frac{\ln\sigma_2^2}{2}+\frac{(y_{ij}-\alpha_i-\beta)^2}{2\sigma_2^2}\right]\right\}\\
&\!\!\!\!\displaystyle\quad+\sum_{i=1}^n\left\{(s_i-1)\left[\frac{\ln\tau_1^2}{2}+\dfrac{(\alpha_i-\mu)^2}{2\tau_1^2}\right]-s_i\left[\frac{\ln\tau_2^2}{2}+\dfrac{(\alpha_i-\mu-\gamma)^2}{2\tau_2^2}\right]\right\}\\
&\!\!\!\!\displaystyle\quad+\sum_{i=1}^n\sum_{j=1}^ts_i\left[\xi_{ij}^{(r)}\ln\lambda+(1-\xi_{ij}^{(r)})\ln(1-\lambda)\right],
\end{array}$
</center>

其中
<center>
$\xi_{ij}^{(r)}=\mathbb E_{\pmb\theta^{(r)}}(z_{ij}\mid y_{ij},s_i)=\dfrac{\lambda^{(r)}f(y_{ij};\alpha_i^{(r)}+\beta^{(r)},(\sigma_2^{(r)})^{2})\cdot s_i}{\lambda^{(r)}f(y_{ij};\alpha_i^{(r)}+\beta^{(r)},(\sigma_2^{(r)})^{2})+(1-\lambda^{(r)})f(y_{ij};\alpha_i^{(r)},(\sigma_1^{(r)})^{2})}.$
</center>

#### M步

这样E步就完成了，接下来通过坐标下降法来依次更新 $\pmb\theta$ 的分量，通过最大化 $Q(\pmb\theta\mid\pmb\theta^{(r)})$ 可以得到相应参数的更新公式如下：

- 更新 $\lambda$：
    <center>
    $\displaystyle\lambda^{(r+1)}=\sum\limits_{i=1}^n\sum\limits_{j=1}^t\xi_{ij}^{(r)}\bigg/\sum\limits_{i=1}^n\sum\limits_{j=1}^ts_i.$
    </center>
    
- 更新 $\pmb\alpha$：
    <center>
    $\alpha_i^{(r+1)}=\dfrac{\sum\limits_{j=1}^t\Big(\frac{(1-\xi_{ij}^{(r)})y_{ij}}{(\sigma_1^{(r)})^2}+\frac{\xi_{ij}^{(r)}(y_{ij}-\beta^{(r)})}{(\sigma_2^{(r)})^2}\Big)+\frac{(1-s_i)\mu^{(r)}}{(\tau_1^{(r)})^2}+\frac{s_i(\mu^{(r)}+\gamma^{(r)})}{(\tau_2^{(r)})^2}}{\sum\limits_{j=1}^t\Big(\frac{1-\xi_{ij}^{(r)}}{(\sigma_1^{(r)})^2}+\frac{\xi_{ij}^{(r)}}{(\sigma_2^{(r)})^2}\Big)+\frac{1-s_i}{(\tau_1^{(r)})^2}+\frac{s_i}{(\tau_2^{(r)})^2}},\quad i=1,\dots,n.$
    </center>
    
- 更新 $\beta,\sigma_1^2,\sigma_2^2$：
    <center>
    $\displaystyle\beta^{(r+1)}=\sum\limits_{i=1}^n\sum\limits_{j=1}^t\xi_{ij}^{(r)}(y_{ij}-\alpha_i^{(r+1)})\bigg/\sum\limits_{i=1}^n\sum\limits_{j=1}^t\xi_{ij}^{(r)},$
    
    $\displaystyle(\sigma_1^{(r+1)})^2=\sum\limits_{i=1}^n\sum\limits_{j=1}^t(1-\xi_{ij}^{(r)})(y_{ij}-\alpha_i^{(r+1)})^2\bigg/\sum\limits_{i=1}^n\sum\limits_{j=1}^t(1-\xi_{ij}^{(r)}),$
    
    $\displaystyle(\sigma_2^{(r+1)})^2=\sum\limits_{i=1}^n\sum\limits_{j=1}^t\xi_{ij}^{(r)}(y_{ij}-\alpha_i^{(r+1)}-\beta^{(r+1)})^2\bigg/\sum\limits_{i=1}^n\sum\limits_{j=1}^t\xi_{ij}^{(r)}.$
    </center>
    
- 更新 $\mu,\gamma,\tau_1^2,\tau_2^2$：
    <center>
    $\mu^{(r+1)}=\dfrac{\sum\limits_{i=1}^n(1-s_i)\alpha_i^{(r+1)}}{\sum\limits_{i=1}^n(1-s_i)},\quad\gamma^{(r+1)}=\dfrac{\sum\limits_{i=1}^ns_i\alpha_i^{(r+1)}}{\sum\limits_{i=1}^ns_i}-\mu^{(r+1)},$
    
    $(\tau_1^{(r+1)})^2=\dfrac{\sum\limits_{i=1}^n(1-s_i)(\alpha_i^{(r+1)}-\mu^{(r+1)})^2}{\sum\limits_{i=1}^n(1-s_i)},\quad(\tau_2^{(r+1)})^2=\dfrac{\sum\limits_{i=1}^ns_i(\alpha_i^{(r+1)}-\mu^{(r+1)}-\gamma^{(r+1)})^2}{\sum\limits_{i=1}^ns_i}.$
    </center>

## 声明

本文内容未经许可，严禁转载与挪用！