title = "Worm Heart 10.4.13"
# requires VPython, MatPlotLib, easygui and numpy modules
from visual.graph import *
import matplotlib.pyplot as plt
from easygui import *
from random import *

c=ccbox(msg="Load parameters file?", title=' ', choices=('Yes','No'))
if c==1:
    f = open(fileopenbox(msg=None, title=None, default='', filetypes=None))
    l = f.readlines()
    fieldDefaults=[]
    for j in range(12):
        fieldDefaults.append(l[j])
    f.close()
else:
    fieldDefaults = [50, 2000, 0.01, 0.22, -65, 3.13, 5, 150, 0.01, 0.09, 1, 3]
fieldNames = ["n", "tstop", "a", "b", "c", "d", "i_gate [all units]",
              "i_syn [unit 1, 2, 3...]", "alpha", "beta", "p", "q"]
fieldValues = multenterbox("Parameters", title, fieldNames, fieldDefaults)
for j in range(11): fieldValues[j] = float(fieldValues[j])
n, tstop, a, b = fieldValues[0], fieldValues[1], fieldValues[2], fieldValues[3]
c, d, i_gate = fieldValues[4], fieldValues[5], fieldValues[6]
i_syn, alpha, beta = fieldValues[7], fieldValues[8], fieldValues[9]
p, q = fieldValues[10], fieldValues[11]
c=ccbox(msg="Save these parameters to file?", title=' ', choices=('Yes', 'No'))
if c==1:
    f = open(filesavebox(msg=None, title=None, default='', filetypes=None), "w")
    for j in range(12):
        print >>f, fieldValues[j]
    f.close()

yoffset = 110 # causes an offset along the y axis from one unit to the next.
tau = 0.02
f = 1
t, num = 0, 0
v, u, g, tstamp, DBV, s = [], [], [], [], [], []
v_out, d_out, syn = [], [], []
    
for j in range(n+1):
    v.append(-65)
    u.append(0)
    g.append(0)
    tstamp.append(0)
    DBV.append(0)
    v_out.append([-65])
    d_out.append([100])
    syn.append(0)

# for j in range(n+1):  # this is an array adaptation of how you had it work...
#    syn[j]=(n-j)/n

decay_coefficient = 1.75
for j in range(int(n/2)):
    syn[j]=(n-decay_coefficient*j)/n
    tmp = (n-decay_coefficient*j)/n
for j in range(int(n/2)+1):
    syn[int(n/2)+j-1] = tmp+decay_coefficient*j/n

for j in range(n): print syn[j]

#scene.autoscale=False
scene.autocenter=True
for j in range(n):
    DBV[j] = cylinder(pos=((j*50),0,1), axis=(50,0,0), radius=n)
    DBV[j].color, DBV[j].material = (1, 0.4, 0.4), materials.diffuse
    curve(x=(j*50)+25, y=-2*n+2*n*sin(arange(0,50,0.1)), z=n*2*cos(arange(0,50,0.1)),
             radius=n/5, color=(1, 0.4, 0.4))
body=cylinder(pos=(0,-2*n,0), axis=(50*n,0,0), radius=(n*3)+(0.2*n))
body.opacity, body.material, body.color = 0.5, materials.marble,(0.75, 0.5, 0.5)
gut=cylinder(pos=(0,-2*n,0), axis=(50*n,0,0), radius=n*2)
gut.material, gut.color = materials.diffuse, (0.9, 0.7, 0.5) 
#VBV=cylinder(pos=(1,-3.25*n,0), axis=(49*n,0,0), radius=n)
#VBV.material, VBV.color = materials.diffuse, (1, 0.4, 0.4)

while t<tstop:
    t += tau
    for j in range(n):
        if j == 0: s = 0
        s = s + i_gate*(n-j)/n  # variable depolarizing drive (stronger in post)
        v[j] += tau * (0.04*v[j]*v[j] + 5*v[j] + 140 - u[j] + s)
        #v[j] += 0.5*random() - 0.25 # +/- 0.5 mV membrane noise
        u[j] += tau * (a * (b*v[j] - u[j]))
        f = 1/(t-tstamp[j])
        if v[j] > 30:
            v[j] = c
            u[j] += d
            tstamp[j] = t
            s = i_syn*syn[j] # variable synaptic strength (stronger in post)
            v_out[j].append(30+yoffset*j)
        else:
            s = 0
            v_out[j].append(v[j]+yoffset*j)
        g[j] += tau * (alpha*f**p*(1-g[j]) - beta*g[j])
        tmp = (1-g[j]**p)
        d_out[j].append(100*tmp-65+yoffset*(j-1))
        DBV[j].radius = n*tmp
        #LBV[j].radius = n*tmp/5 # If they contract with the DBV segment
        #LBV[j].radius = n/10+n*(1-tmp)/5 # If they dilate when the DBV contracts

# plot
for j in range(n):
    plt.plot(v_out[j])
    plt.plot(d_out[j])
plt.show()
