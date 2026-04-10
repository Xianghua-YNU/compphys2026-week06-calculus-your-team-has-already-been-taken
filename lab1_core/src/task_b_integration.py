import math
import numpy as np

def debye_integrand(x: float) -> float:
    if abs(x) < 1e-12:
        return 0.0
    ex = math.exp(x)
    return (x**4) * ex / ((ex - 1.0) ** 2)


def trapezoid_composite(f, a: float, b: float, n: int) -> float:
    h = (b-a)/n
    x_series = np.arange(a ,b+h ,h)
    w = np.ones(n+1)*h
    w[0] = h/2
    w[-1] = h/2
    integ = 0
    for i in range(n+1):
        integ += w[i]*f(x_series[i])
    return integ
    # TODO B1: 实现复合梯形积分
    raise NotImplementedError("TODO B1")


def simpson_composite(f, a: float, b: float, n: int) -> float:
    if n % 2 == 0:
        print('n是偶数，可以进行')
    else:
        raise ValueError('n是奇数或者其他类型的数字，不可以进行计算')
    
    integ = 0
    h = (b-a)/n
    x_series = np.arange(a ,b+h ,h)
    k = np.arange(0 ,n/2 ,dtype=int) 
    w = h*np.ones(n+1 ,dtype=int)
    w[2*k] = 2*h/3
    w[2*k+1] = 4*h/3
    w[0] = h/3
    w[-1] = h/3
    for i in range(n+1):
        integ += w[i]*f(x_series[i])
    return integ
    # TODO B2: 实现复合 Simpson 积分，并检查 n 为偶数
    raise NotImplementedError("TODO B2")


def debye_integral(T: float, theta_d: float = 428.0, method: str = "simpson", n: int = 200) -> float:
    # TODO B3: 计算 Debye 积分 I(theta_d/T)
    y = theta_d/T
    if method  == 'simpson':
        return simpson_composite(debye_integrand, 0 ,y , n)
    else:
        return trapezoid_composite(debye_integrand,0 , y, n)
    raise NotImplementedError("TODO B3")


