import numpy as np


def rate_3alpha(T: float) -> float:
    if T <= 0:
        raise ValueError("温度T必须大于0")
    T8 = T / 1.0e8
    return 5.09e11 * (T8 ** (-3.0)) * np.exp(-44.027 / T8)


def finite_diff_dq_dT(T0: float, h: float = 1e-8) -> float:
    # TODO A1: 使用前向差分实现 dq/dT
    delta_T=h*T0
    q0=rate_3alpha(T0)
    qh=rate_3alpha(T0+delta_T)
    return (qh-q0)/delta_T
    


def sensitivity_nu(T0: float, h: float = 1e-8) -> float:
    # TODO A2: 根据 nu = (T/q) * dq/dT 计算温度敏感性指数
    nu=(T0/rate_3alpha(T0))*finite_diff_dq_dT(T0,h)
    return nu


def nu_table(T_values, h: float = 1e-8):
    # TODO A3: 返回 [(T, nu(T)), ...]
    return [(T, sensitivity_nu(T, h)) for T in T_values]

Test_T=[1.0e8,1.5e8,2.0e8,2.5e8,3.0e8,5.0e8,1.0e9]
Table=nu_table(Test_T)
for T, nu in Table:
        print(f"{T:.4e}    | {nu:.2f}")