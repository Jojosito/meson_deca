import convert
import model
import numpy as np
import matplotlib.pyplot as plt

N = 100000
#x = np.linspace(0.1,2,100)
#z = np.asarray([[convert.ComplexVectorForm(model.A_cv(5,[ _x,y,1,1,1])) for _x in x] for y in x])
#x = []#np.zeros([N,5])
#z = []#np.zeros(N,dtype=complex)

#for i in list(range(N)):
#    a = 0.2 + 1.6 * np.random.rand()
#    b = 0.2 + 1.6 * np.random.rand()
#    c = 0.2 + 1.6 * np.random.rand()
#    d = 0.2 + 1.6 * np.random.rand()
#    e = 0.2 + 1.6 * np.random.rand()
#    z[i] = convert.ComplexVectorForm(model.A_cv(5,[a,b,c,d,e]))

def f_list(num):
    # Return list of 5 random numbers between var_min and var_max
    var_min = 0
    var_max = 2
    return [var_min + (var_max - var_min) * np.random.rand() for i in list(range(num))]

x = np.linspace(0.0,8.0,100)
z = np.asarray([[convert.ComplexVectorForm(model.A_cv( 5, [_x, y] + f_list(3) )) for _x in x] for y in x])

z[np.isnan(z)] = 0.
var_names = ['F_P','F_R_1','F_R_2','Z_1','Z_2','T_R_1','T_R_2','Total']
#f, ax = plt.subplots(8, sharex=True)
plt.clf()
for i in list(range(8)):
    plt.subplot(3,3,i+1)
    plt.pcolor(x,x,np.abs(z[:,:,i]))
    plt.colorbar()
    plt.title(var_names[i])

plt.savefig('PWA_components_old.png')
plt.show()
