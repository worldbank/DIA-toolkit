
# Distributional Impact Analysis: Toolkit and Illustrations of Impacts Beyond the Average Treatment Effect

Program evaluations often focus on average treatment effects. However, average treatment effects miss important aspects of policy evaluation, such as the impact on inequality and whether treatment harms some individuals. A growing literature develops methods to evaluate such issues by examining the distributional impacts of programs and policies. This toolkit reviews methods to do so, focusing on their application to randomized control trials. This repository contains the programs to implement these methods and illustrates select methods using data from two randomized evaluations.

The toolkit emphasizes two strands of the literature: estimation of impacts on outcome distributions and estimation of the distribution of treatment impacts. We discuss extensions to conditional treatment effect heterogeneity, that is, to analyses of how treatment impacts vary with observed characteristics, offer advice on inference, testing, and power calculations, which are important when implementing distributional analyses in practice.

## Questions of Interest

In terms of impact on the outcome distributions, consider questions such as:

* **Does microfinance boost average incomes?** 
To answer this question, we would estimate the average treatment effect, E(Y_1 )-E(Y_0 ), or the average treatment effect on the treated, E(Y_1│T=1)-E(Y_0│T=1).
* **Does hospital regulation raise minimum levels of patient safety?** 
To answer this question, we would estimate the treatment effect on the minimum, min⁡〖Y_1 〗- min⁡〖Y_0 〗, or on a quantile in the left tail, such as q_1 (0.1)-q_0 (0.1).
* **Does education reform decrease dispersion of student’s test scores?** 
To answer this question, we could estimate the treatment effect on a measure of inequality, such as the variance, var⁡〖(Y_1)〗-var⁡〖(Y_0)〗, or the Gini index.


Analyses of the distribution of treatment effects can answer questions like:

* **What proportion of students benefit from an educational reform?**
For this question, we would compute P⁡〖(Y_1>Y_0 )=P(Δ>0)〗.
* **Are the improvements in average patient outcomes from health facility inspections driven by a few people who benefit considerably? Formally: is there significant skewness in effects?**
For these questions, we would compute E{[((Y_1-Y_0-μ))⁄σ]^3 }, where μ=E(Y_1-Y_0 ) and σ^2=E[(Y_1-Y_0-μ )^2 ].
* **What is the median impact of a microfinance program? More generally, what are the quantiles of the impact distribution, like the minimum or maximum program impact?**
For this question, we would compute the relevant quantile of treatment effects.


