import math
import numpy as np
import matplotlib.pyplot as plt

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

# # 误差计算与检查
# E_h_simp = debye_integral(80 ,428.0 ,'simpson' ,200)
# E_h_2_simp = debye_integral(80 ,428.0 ,'simpson' ,100)
# E_simp = (E_h_simp - E_h_2_simp)/(2**4-1)
# print(f'辛普森方法的误差在n=200的时候为{E_simp}')

# E_h_tr = debye_integral(80 ,428.0 ,'tr'  ,200)
# E_h_2_tr = debye_integral(80 ,428.0 , 'tr' ,100)
# E_tr = (E_h_tr - E_h_2_tr)/(2**2-1)
# print(f'tr方法的误差在n=200的时候为{E_tr}')


# plot the result
temp = np.arange(10 , 600,5)
thermal_capacity = []
theta_d = 438.0
for temp_i in temp:
    C = thermal_capacity.append((temp_i/theta_d)**3 *debye_integral(temp_i,theta_d ,'simpson' ,200))
thermal_capacity = np.array(thermal_capacity)


# plt.plot(temp ,thermal_capacity ,'r--' ,label = 'result of integral')
# plt.title('result of integral')
# plt.legend()
# plt.xlabel('T')
# plt.ylabel('I')

plt.plot(temp ,thermal_capacity ,'r--' ,label = 'thermal_capacity')
plt.title('result of integral')
plt.legend()
plt.xlabel('T')
plt.ylabel('thermal_capacity')