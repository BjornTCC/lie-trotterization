import numpy as np

def one_d_optimizer(f: callable, lb: float, ub: float, points: int = 10000, strict_positive: bool = False) -> tuple[float, float]:
    xs = np.linspace(lb,ub, points)
    ys = []
    for x in xs:
        y = f(x)
        if y < 0 and strict_positive:
            ys.append(np.inf)
        else:
            ys.append(y)
    return xs[np.argmin(ys)], np.min(ys)