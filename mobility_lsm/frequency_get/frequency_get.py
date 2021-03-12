from math import *
import matplotlib.pylab as pyl

data_ion = open("diel_ion.dat", 'r',encoding='utf-8')
data_ele = open("diel_ele.dat", 'r',encoding='utf-8')

real_ion = []
real_ele = []

temp = 0
for line in data_ion:
    if line == ' \n':
        temp += 1
    else:
        if temp == 2:
            s = line.strip().split('\t')
            s1 = ' '.join(s) + '\n'
            s2 = s1.split()
            real_ion.append([float(s2[0])*5.30883742,float(s2[1])])
data_ion.close()

temp = 0
for line in data_ele:
    if  line == ' \n':
        temp += 1
    else:
        if temp == 2:
            s = line.strip().split('\t')
            s1 = ' '.join(s) + '\n'
            s2 = s1.split()
            real_ele.append([float(s2[0]),float(s2[1])])
data_ele.close()

def get_zeropoint(a1,a2):
    a = [0,0]
    if a1[1]*a2[1] <= 0:
        a[0] = a1[0]+abs(a1[1])/(abs(a1[1])+abs(a2[1]))*(a2[0]-a1[0])
        return a
    else:
        return 0

real_ion_zero = []
for i in range(len(real_ion)-1):
    a = (get_zeropoint(real_ion[i],real_ion[i+1]))
    if a != 0:
        real_ion_zero.append(a)
n_zero = len(real_ion_zero)

print('eps0_ion = ',real_ion[0][1])
print('eps0_ele = ',real_ele[0][1])

data_diel = open("diel.dat", 'w', encoding='utf-8')
x = []
y = []
x_zero = []
y_zero = []
for i in range(len(real_ion)):
    if i < n_zero:
        x_zero.append(real_ion_zero[i][0])
        y_zero.append(real_ion_zero[i][1])
        data_diel.write(str(real_ion[i][0])+'\t'+str(real_ion[i][1])+'\t'+str(real_ion_zero[i][0])+'\t'+str(real_ion_zero[i][1])+'\n')
    else:
        data_diel.write(str(real_ion[i][0]) + '\t' + str(real_ion[i][1]) + '\n')
    x.append(real_ion[i][0])
    y.append(real_ion[i][1])

pyl.plot(x,y)
pyl.plot(x_zero,y_zero,'o',markersize=1)
for a,b in zip(x_zero,y_zero):
    pyl.text(a,b+10,'%.0f'%a,ha = 'center',va = 'bottom',fontsize=7)
pyl.savefig('diel.png')
pyl.show()