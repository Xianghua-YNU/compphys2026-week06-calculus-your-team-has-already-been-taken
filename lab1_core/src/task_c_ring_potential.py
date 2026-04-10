import numpy as np


def ring_potential_point(x: float, y: float, z: float, a: float = 1.0, q: float = 1.0, n_phi: int = 720) -> float:
    # 用离散积分计算单点电势
    phi = np.linspace(0, 2 * np.pi, n_phi)
    dphi = 2 * np.pi / n_phi
    
    # 计算距离的平方，添加小量以提高数值稳定性
    r_squared = (x - a * np.cos(phi))**2 + (y - a * np.sin(phi))**2 + z**2
    epsilon = 1e-10  # 小量，避免除零
    r = np.sqrt(r_squared + epsilon)
    
    # 计算积分
    integral = np.sum(1 / r) * dphi
    potential = q / (2 * np.pi) * integral
    
    return potential


def ring_potential_grid(y_grid, z_grid, x0: float = 0.0, a: float = 1.0, q: float = 1.0, n_phi: int = 720):
    # 在 yz 网格上计算电势矩阵
    phi = np.linspace(0, 2 * np.pi, n_phi)
    dphi = 2 * np.pi / n_phi
    
    # 创建网格点
    Y, Z = np.meshgrid(y_grid, z_grid)
    
    # 初始化电势矩阵
    potential = np.zeros_like(Y)
    
    # 计算积分
    for phi_val in phi:
        # 计算每个 phi 值对应的距离
        r_squared = (x0 - a * np.cos(phi_val))**2 + (Y - a * np.sin(phi_val))**2 + Z**2
        epsilon = 1e-10  # 小量，避免除零
        r = np.sqrt(r_squared + epsilon)
        potential += 1 / r * dphi
    
    # 乘以系数
    potential *= q / (2 * np.pi)
    
    return potential


def axis_potential_analytic(z: float, a: float = 1.0, q: float = 1.0) -> float:
    return q / np.sqrt(a * a + z * z)


def calculate_electric_field(potential, y_grid, z_grid):
    # 计算电场：E = -∇V
    # 计算 y 方向的梯度
    dy = y_grid[1] - y_grid[0] if len(y_grid) > 1 else 1.0
    dVz_dy = np.gradient(potential, dy, axis=1)
    
    # 计算 z 方向的梯度
    dz = z_grid[1] - z_grid[0] if len(z_grid) > 1 else 1.0
    dVz_dz = np.gradient(potential, dz, axis=0)
    
    # 电场分量
    Ey = -dVz_dy
    Ez = -dVz_dz
    
    return Ey, Ez


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # 创建 yz 网格
    y = np.linspace(-3, 3, 100)
    z = np.linspace(-3, 3, 100)
    
    # 计算电势分布
    potential = ring_potential_grid(y, z, x0=0.0, a=1.0, q=1.0)
    
    # 计算电场分布
    Ey, Ez = calculate_electric_field(potential, y, z)
    
    # 创建网格点
    Y, Z = np.meshgrid(y, z)
    
    # 绘制等势线图和电场矢量
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # 第一张图：等势线图
    contour = ax1.contour(Y, Z, potential, levels=20, cmap='viridis')
    ax1.set_title('等势线图')
    ax1.set_xlabel('y')
    ax1.set_ylabel('z')
    ax1.set_aspect('equal')
    ax1.colorbar(contour, ax=ax1, label='电势')
    
    # 第二张图：电场矢量分布（箭头）
    # 为了避免箭头过于密集，使用步长采样
    stride = 5
    ax2.quiver(Y[::stride, ::stride], Z[::stride, ::stride], 
               Ey[::stride, ::stride], Ez[::stride, ::stride],
               scale=50, color='red')
    ax2.set_title('电场矢量分布')
    ax2.set_xlabel('y')
    ax2.set_ylabel('z')
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    
    # 保存图片
    plt.savefig('ring_potential.png', dpi=150)
    plt.show()
    
    print("计算完成，图片已保存为 ring_potential.png")
