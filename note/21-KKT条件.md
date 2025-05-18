# 约束问题的最优性条件

## 一阶必要条件

### Lagrange 函数

$\mathcal{L}:R^{n+m}\to R$$$ \begin{equation}
    \nonumber
        \begin{split}
            \mathcal{L}(x,\lambda,\mu)&=f(x)-\lambda^Tg_I(x)-\mu^Th_E(x)\\ 
            &=f(x)-\sum_{i\in I}\lambda_ig_i(x)-\sum_{j\in E}\mu_jh_j(x)
        \end{split}
\end{equation} $$

这个函数$\mathcal{L}$就是**Lagrange函数**，$\lambda,\mu$就是**Lagrange乘子**

下面还有$$ \nabla_x\mathcal{L}(x,\lambda,\mu)=\nabla f(x)-\sum_{i\in I}\lambda_i\nabla g_i(x)-\sum_{j\in E}\mu_j\nabla h_j(x) $$$$ \nabla _x^2\mathcal{L}(x,\lambda,\mu)=\nabla^2f(x)-\sum_{i\in I}\lambda_i\nabla ^2g_i(x)-\sum_{j\in E}\mu_j\nabla^2h_j(x) $$

> Lagrange函数的形式中，对于优化以及等式约束，其可以看作求解 $\nabla_x\mathcal{L}=0,\nabla_\mu\mathcal{L}=0$
> 
> 对于不等式约束，其$\lambda\geq 0$主要是给违反约束添加了惩罚。

### KKT条件

设$x^*$是一个局部最优解，如果ACQ成立，即$$ \text{SFD}(x^*,D)=\text{LFD}(x^*,D) $$，则存在Lagrange乘子向量$\lambda^*,\mu^*$使得$$ \begin{equation}
    \nonumber
    \left\{
        \begin{split}
            \nabla_x\mathcal{L}(x^*,\lambda^*,\mu^*)=\nabla f(x^*)-\sum_{i\in I}\lambda_i^*\nabla g_i(x^*)-\sum_{j\in E}\mu_j^*\nabla h_j(x^*)=0\\
            \nabla_\mu\mathcal{L}(x^*,\lambda^*,\mu^*)=0\iff h_j(x)=0,\quad j\in E\\
            g_i(x^*)\geq 0,\quad \lambda^*\geq 0,\quad \lambda_i^*g_i(x^*)=0,\quad i\in I
        \end{split}
    \right.
\end{equation} $$

这就是满足问题的**一阶必要条件**，称为**KKT条件**

满足KKT条件的点$x^*$称为**KKT点**

> 其中，$\lambda^*g_i(x^*)=0,i\in I$ 称为**互补松弛条件**
> 这就是使用Lagrange函数的思想改进

> 需要注意，**并不是所有条件下，局部最优点都是KKT点**
>
> **必须还要满足ACQ条件才可**

### KKT条件的向量形式

$$ H(x,\lambda,\mu)=\begin{pmatrix}
    \nabla_x\mathcal{L}(x,\lambda,\mu)\\
    h_E(x)\\
    \min\{ \lambda_I,g_i(x) \}
\end{pmatrix}=0 $$

