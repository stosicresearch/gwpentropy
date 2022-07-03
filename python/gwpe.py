import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import numpy as np
import gwpe

# Load data.
series = 'logisticmap'
x = np.loadtxt('../data/'+series+'.dat')
n = x.shape[0]

# Parameters.
w = 6           # word size
qdelta = 0.1
qmin = -10
qmax = 10

# Generalized Weighted Permutation Entropy.
Hlist = []
Clist = []
qlist = []
for q in np.arange(qmin,qmax+qdelta,qdelta):
    H, C = gwpe.calc(x, n, w, q)
    qlist.append(q)
    Hlist.append(H)
    Clist.append(C)
H = np.asarray(Hlist)
C = np.asarray(Clist)
q = np.asarray(qlist)

# Plot signature curve.
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(H, C, q)#, marker=m)
ax.set_xlabel('H')
ax.set_ylabel('C')
ax.set_zlabel('q')
plt.show()
fig.savefig(series + '.png')