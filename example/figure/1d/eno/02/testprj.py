import numpy as np
import matplotlib.pyplot as plt
    
def plot_all_cell_center( xcc, yref ):
    plt.scatter(xcc, np.full_like(xcc, yref), s=20, facecolor='black', edgecolor='black', linewidth=1)
    return
    
def plot_cell_center( xcc, yref ):
    nx = xcc.size
    ii = nx // 2
    im = ii - 1
    ip = ii + 1
    xcc_new = []
    for i in range(0, nx):
        if i > 0 and i < im:
            continue
        if i > ip and i < nx-1:
            continue
        xcc_new.append( xcc[i] )
    plt.scatter(xcc_new, np.full_like(xcc_new, yref), s=20, facecolor='black', edgecolor='black', linewidth=1)
    return    

def plot_mesh( x, yref ):
    dx = x[1] - x[0]
    dy = 0.1 * dx
    for xm in x:
        plt.plot([xm, xm], [yref-dy, yref+dy], 'k-')  # 绘制垂直线
    #
    
    nxc = x.size - 1
    ii = nxc // 2
    im  = ii - 1
    im1 = ii - 2
    ip  = ii + 1
    ip1 = ii + 2
    
    for i in range(0, nxc):
        if i > 0 and i < im1:
            plt.plot([x[i], x[i+1]], [yref, yref], 'k--', linewidth=1)
        elif i > ip1 and i < nx-1:
            plt.plot([x[i], x[i+1]], [yref, yref], 'k--', linewidth=1)
        else :
            plt.plot([x[i], x[i+1]], [yref, yref], 'b-', linewidth=1)
    #plt.plot(x, np.full_like(x, yref), 'k--', linewidth=1)
    return
   
def plot_label(x, xcc, yref):
    dx = x[1] - x[0]
    dyb = 0.5 * dx
    dyt = dyb * 0.6
    yb = yref - dyb
    yt = yref + dyt
    ybc = yref - 0.5* dyb
    plt.text(x[0], yb, r'$x_{i-r-\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[1], yb, r'$x_{i-r+\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[2], yb, r'$x_{i-\frac{5}{2}}$', fontsize=12, ha='center')
    plt.text(x[3], yb, r'$x_{i-\frac{3}{2}}$', fontsize=12, ha='center')
    plt.text(x[4], yb, r'$x_{i-\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[5], yb, r'$x_{i+\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[6], yb, r'$x_{i+\frac{3}{2}}$', fontsize=12, ha='center')
    plt.text(x[7], yb, r'$x_{i+\frac{5}{2}}$', fontsize=12, ha='center')
    plt.text(x[8], yb, r'$x_{i+s-\frac{1}{2}}$', fontsize=12, ha='center')
    plt.text(x[9], yb, r'$x_{i+s+\frac{1}{2}}$', fontsize=12, ha='center')
    
    plt.text(xcc[0], ybc, r'$i-r$', fontsize=12, ha='center')
    plt.text(xcc[1], ybc, r'$\cdots\cdots$', fontsize=12, ha='center')
    plt.text(xcc[-2], ybc, r'$\cdots\cdots$', fontsize=12, ha='center')
    plt.text(xcc[-1], ybc, r'$i+s$', fontsize=12, ha='center')
    
    nx = xcc.size
    i = nx // 2
    print("i=",i)
    im  = i - 1
    im1 = i - 2
    ip  = i + 1
    ip1 = i + 2
    
    plt.text(xcc[im1], ybc, r'$i-2$', fontsize=12, ha='center')
    plt.text(xcc[im], ybc, r'$i-1$', fontsize=12, ha='center')
    plt.text(xcc[i], ybc, r'$i$', fontsize=12, ha='center')
    plt.text(xcc[ip], ybc, r'$i+1$', fontsize=12, ha='center')
    plt.text(xcc[ip1], ybc, r'$i+2$', fontsize=12, ha='center')
 
    return    

# 设置字体为 Times New Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

# 设置图形大小和样式
plt.figure(figsize=(12, 5))

nx = 9
L  = 1.0
x_l = 0.0
dx = L / nx

x   = np.zeros(nx+1, dtype=np.float64)
xcc = np.zeros(nx, dtype=np.float64)

for i in range(0, nx+1):
    x[i] = x_l + dx*(i)
    
for i in range(0, nx):
    xcc[i] = 0.5*(x[i]+x[i+1])
    
print("x=",x)
print("xcc=",xcc)

yref = 0.0

plot_cell_center( xcc, yref )
plot_mesh( x, yref )
plot_label(x, xcc, yref)
    
plt.axis('equal')
plt.axis('off')

plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
plt.show()