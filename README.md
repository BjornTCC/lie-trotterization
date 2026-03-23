# A codebase for lie-algebraic accelerated trotter steps

The optimal trotter formula which is correct to order $t^4$, is:

$$
    e^{t(A+B)} \approx e^{t \lambda_1 A / 2} e^{t \lambda_1 B} e^{t (\lambda_1 + \lambda_2) A / 2} e^{t \lambda_2 B}e^{t (\lambda_2 + \lambda_1) A / 2} e^{t \lambda_1 B} e^{t \lambda_1 A / 2}
$$

with $\lambda_2 = -\sqrt[3]{2} \lambda_1 = -\sqrt[3]{2} / (2 -\sqrt[3]{2}) $.

However, assuming access to the commutators $[A,[A,B]]$ and $[B,[A,B]]$, we can get a shorter formulae:

$$
    e^{t(A+B)} \approx e^{t A / 2} e^{t B/2} e^{\frac{t^3}{24} [A,[A,B]]} e^{\frac{t^3}{12} [B,[A,B]]} e^{t  B / 2} e^{t / 2}
$$