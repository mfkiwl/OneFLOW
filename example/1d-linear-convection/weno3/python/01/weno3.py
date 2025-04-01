import numpy as np
import matplotlib.pyplot as plt

# 定义求解区域的网格和时间步长
nx = 100  # 网格数
L = 1     # 区域长度
dx = L / nx  # 网格步长
dt = 0.01  # 时间步长
nt = 10  # 时间步数
a = 1     # 对流速度

# 定义线性对流方程的初始条件和边界条件
x = np.linspace(0, L, nx)
u0 = np.sin(2 * np.pi * x)  # 初始条件
u = u0.copy()
u_bc = np.concatenate((u[-3:], u, u[:3]))  # 周期边界条件

# 定义WENO方法的参数
epsilon = 1e-6  # 避免分母为0
p = 2  # 权重参数

# 循环求解线性对流方程
for n in range(nt):
    # 计算每个子模板的光滑度指标
    beta_0 = (u_bc[2:nx+2] - u_bc[1:nx+1])**2
    beta_1 = (u_bc[3:nx+3] - u_bc[2:nx+2])**2
    beta_3 = (1/3) * ((u_bc[1:nx+1] - 2*u_bc[2:nx+2] + u_bc[3:nx+3]) / dx)**2 + \
             (1/4) * ((u_bc[3:nx+3] - u_bc[1:nx+1]) / dx)**2

    # 计算权重
    d_0 = 1/3
    d_1 = 2/3
    d_3 = 1/3
    alpha_0 = d_0 / (beta_0 + epsilon)**p
    alpha_1 = d_1 / (beta_1 + epsilon)**p
    alpha_3 = d_3 / (beta_3 + epsilon)**p
    w_0 = alpha_0 / (alpha_0 + alpha_1 + alpha_3)
    w_1 = alpha_1 / (alpha_0 + alpha_1 + alpha_3)
    w_3 = alpha_3 / (alpha_0 + alpha_1 + alpha_3)

    # 计算重构值
    f_0 = -0.5 * u_bc[1:nx+1] + 1.5 * u_bc[2:nx+2]
    f_1 = 0.5 * u_bc[2:nx+2] + 0.5 * u_bc[3:nx+3]
    f_3 = (-1/6) * u_bc[1:nx+1] + (5/6) * u_bc[2:nx+2] + (1/3) * u_bc[3:nx+3]
    f_weno = w_0 * f_0 + w_1 * f_1 + w_3 * f_3

    # 更新解
    u = u - dt / dx * (f_weno - np.roll(f_weno, 1))

    # 更新边界条件
    u_bc = np.concatenate((u[-3:], u, u[:3]))

# 计算理论解
t = nt * dt
u_exact = np.sin(2 * np.pi * (x - a * t))
print("at=",a * t)
print("u0=",u0)
print("u_exact=",u_exact)

# 绘制结果
#plt.plot(x, u, label='WENO Solution')
plt.plot(x, u_exact, label='Exact Solution', linestyle='--')
plt.plot(x, u0, label='init Solution', linestyle='--')
plt.legend()
plt.title('Comparison of WENO Solution and Exact Solution')
plt.xlabel('x')
plt.ylabel('u')
plt.show()