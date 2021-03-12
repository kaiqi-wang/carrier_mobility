from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

eps_inf = 4.7                       #需要输入的参数，高频介电常数
eps0 = 32.3                         #需要输入的参数，静介电常数
m0 = 9.1*10**(-31)
m = 0.117*m0                        #需要输入的参数，有效质量
T = 300                             #需要设定的参数，温度，单位K
e = 1.602*10**(-19)
h = 6.62606957*10**(-34)
c = 3*10**8
c_cm = 3*10**10                     #光速，单位cm/s
kB = 1.3806505*10**(-23)

def A(v, w):                        #表达式A
    if v == w:
        a = 3*(-log(2*pi*beta)/2)/beta
    else:
        a = 3 * (log(v / w) - log(2 * pi * beta) / 2 - log(sinh(v * beta / 2) / sinh(w * beta / 2))) / beta
    return a

def Y(x,v):                         #表达式Y
    y = (1 + exp(-1*v*beta) - exp(-1*v*x) - exp(v*(x-beta)))/(1-exp(-1*v*beta))
    return y

def X(v, w):                        #表达式C
    x = 3*(v**2 - w**2)*(cosh(v*beta/2)/sinh(v*beta/2)-2/(v*beta))/(4*v)
    return x

def B(v, w):                        #表达式B，其中数值n = 100是积分取点，可以视情况调整
    n = 100
    desc = beta / 2 / n
    y = 0
    for i in range(1, n):
        x = i * desc
        y = y + (exp(beta-x)+exp(x))/(sqrt(w**2*x*(1-x/beta)+Y(x,v)*(v**2-w**2)/v))*desc
    b = af*v/(sqrt(pi)*(exp(beta)-1))*y
    return b

def R(v,w):                         #表达式R
    r = (v**2-w**2)/(w**2*v)
    return r

def get_a(r,v):                     #表达式a
    a = sqrt((beta/2)**2+r*beta*cosh(beta*v/2)/sinh(beta*v/2))
    return a

def get_b(r,v):                     #表达式b，SI中的公式有误，sinh内无b
    b = r*beta/sinh(beta*v/2)
    return b

def K(a,b,v):                       #表达式K
    y = 0
    desc = 1/200000
    for i in range(1,200000):
        u = i*desc
        y1 = cos(u)/(u**2+a**2-b*cos(v*u))**(3/2)
        y2 = 1/(u**2)*cos(1/u)/(1/(u**2)+a**2-b*cos(v/u))**(3/2)
        y = y + (y1+y2)*desc
    return y

def Mob(k,v,w):                     #求解迁移率的函数
    u = 100*(3*sqrt(pi)*e/(omega*m*af*c*2*pi))*(sinh(beta/2)/beta**2.5)*(w**3/v**3)/k
    print('u = ', u)
    return u


#------------------------------------------------------------------#
#                          主函数部分
#------------------------------------------------------------------#

f = open('frequency.dat','r', encoding='utf-8')         #文件中第一列是WL0，第二列是WT0
wl0 = []
wt0 = []
WW = []
for line in f:
    s = line.strip().split('\t')
    wl0.append(float(s[0]))
    wt0.append(float(s[1]))
    WW.append((float(s[0])**2-float(s[1])**2)/eps_inf)
f.close()

sum = 0
ww = 0
for i in range(len(wl0)):
    sum = sum + WW[i]/wl0[i]/tanh(h*c_cm*wl0[i]/(2*kB*T))
    ww = ww + (wl0[i] ** 2 - wt0[i] ** 2) / eps_inf

i = -1
wl1 = 0.001
wl2 = 0.002
while i < 0:
    if (sum*tanh(h*c_cm*wl1/(2*kB*T))*wl1-ww)*(sum*tanh(h*c_cm*wl2/(2*kB*T))*wl2-ww) > 0:
        wl1 = wl2
        wl2 = wl2 + 0.001
    else:
        omega = wl1
        i = 1

beta = 4.8*10**(-3)*omega*300/T                                         #beta定义式
af = (1/eps_inf-1/eps0)*331.059*sqrt(m/m0)*sqrt(1/omega)                #af定义式，对应公式（3）

print("omega = ", omega)
print("af = ", af)
print("beta = ", beta)
print("eps_inf = ", eps_inf)
print("eps0 = ", eps0)
print("1/eps = ", 1/eps_inf-1/eps0)

x_mesh = []                                                             #求解w和v
y_mesh = []
z_mesh = []
for i in range(10,500):                                                 #第一次筛旋，w和v的精度0.1
    for j in range(10,500):
        v = i/10
        w = j/10
        F = -(A(v, w) + B(v, w) + X(v, w))
        x_mesh.append(v)
        y_mesh.append(w)
        z_mesh.append(F)

min_num = z_mesh[0]
min_index = 0
for i in range(len(z_mesh)):
    if z_mesh[i]<min_num:
        min_num = z_mesh[i]
        min_index = i

temp_v = int(x_mesh[min_index])
temp_w = int(y_mesh[min_index])

x_mesh = []
y_mesh = []
z_mesh = []
for i in range((temp_v-1)*100,(temp_v+1)*100):                          #第二次筛选，w和v的精度0.01
    for j in range((temp_w-1)*100,((temp_w+1)*100)):
        v = i/100
        w = j/100
        F = -(A(v, w) + B(v, w) + X(v, w))
        x_mesh.append(v)
        y_mesh.append(w)
        z_mesh.append(F)

min_num = z_mesh[0]
min_index = 0
for i in range(len(z_mesh)):
    if z_mesh[i]<min_num:
        min_num = z_mesh[i]
        min_index = i

print('v = ', x_mesh[min_index])
print('w = ', y_mesh[min_index])
v = x_mesh[min_index]
w = y_mesh[min_index]


#画图部分，可选
#fig = plt.figure()
#sub = fig.add_subplot(111, projection='3d')  # 3d表示三维图像
#sub.plot_trisurf(x_mesh, y_mesh, z_mesh, color='blue')

#plt.savefig("1.png")
#fig.show()

r = R(v,w)
a = get_a(r,v)
b = get_b(r,v)
k = K(a,b,v)
u = Mob(k,v,w)

