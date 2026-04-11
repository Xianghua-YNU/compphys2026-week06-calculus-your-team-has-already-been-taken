import numpy as np
from scipy.special import roots_legendre
import matplotlib.pyplot as plt

G = 6.674e-11

def roots_Legendre(n: int):
    """
    生成 n 点高斯-勒让德积分的节点 xi 和权重 wi
    这是你缺失的关键函数
    """
    xi, wi = np.polynomial.legendre.leggauss(n)
    return xi, wi

def gauss_legendre_2d(func, ax: float, bx: float, ay: float, by: float, n: int = 40) -> float:
    # TODO D1: 使用二维高斯-勒让德积分实现双重积分
    xi,wi=roots_Legendre(n)
    x = (bx - ax) / 2 * xi + (bx + ax) / 2
    y = (by - ay) / 2 * xi + (by + ay) / 2  # 二维积分x/y用同一组节点
    jacobian = (bx - ax) * (by - ay) / 4.0  # 雅可比行列式
    X, Y = np.meshgrid(x, x)  # 生成二维网格
    W = np.outer(wi, wi)      # 权重张量积
    integral = jacobian * np.sum(W * func(X, Y))
    
    return integral


def plate_force_z(z: float, L: float = 10.0, M_plate: float = 1.0e4, m_particle: float = 1.0, n: int = 40) -> float:
    # TODO D2: 计算方板中心正上方 z 位置的 Fz
    sigma = M_plate / (L * L)
    def f(x,y):
         r=x**2+y**2+z**2
         return z/r**(3/2)
    ax, bx = -L/2, L/2
    ay, by = -L/2, L/2
    integral_val = gauss_legendre_2d(f, ax, bx, ay, by, n)
    Fz = G * sigma * m_particle * integral_val
    return Fz
    


def force_curve(z_values, L: float = 10.0, M_plate: float = 1.0e4, m_particle: float = 1.0, n: int = 40):
    # TODO D3: 返回 z_values 对应的 Fz 数组
    Fz_array = np.array([
        plate_force_z(z, L, M_plate, m_particle, n)
        for z in z_values
    ])
    return Fz_array
    
if __name__ == "__main__":

    # 1. 输出要求的表格
    Test_Z = [0.2, 1.0, 5.0, 10.0]
    Fz_results = force_curve(Test_Z)

    print("z (m)      | Fz (N)")
    print("-" * 26)
    for z, fz in zip(Test_Z, Fz_results):
        print(f"{z:>5.1f}      | {fz:.6e}")

    # 2. 生成引力-距离曲线
    z_dense = np.linspace(0.2, 10.0, 100)  # 密集点，画平滑曲线
    fz_dense = force_curve(z_dense)

    plt.figure(figsize=(8, 5))
    plt.plot(z_dense, fz_dense, 'b-', linewidth=2.5, label='Gravitational Force $F_z$')
    plt.scatter(Test_Z, Fz_results, color='red', s=60, zorder=5, label='Key Points')
    plt.xlabel('z (m)', fontsize=12)
    plt.ylabel('$F_z$ (N)', fontsize=12)
    plt.title('Gravitational Force vs Distance above Square Plate', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()