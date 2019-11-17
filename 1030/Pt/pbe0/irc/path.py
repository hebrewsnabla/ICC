import numpy as np
import matplotlib.pyplot as plt

H2KCAL = 627.509

def get_xyz(raw_xyz):
    xyz = []
    for line in raw_xyz[:-1]:
        line = line.strip().split()
        x = float(line[-3])
        y = float(line[-2])
        z = float(line[-1])
        xyz.append(np.array([x,y,z]))
    return np.array(xyz)


with open('ircE.txt','r') as f1:
    ircdata = f1.readlines()
ircE = []
for line in ircdata:
    line = line.split('=')
    E = float(line[1].strip().split()[0])
    ircE.append(E)
ircE = ircE

with open('1Pt.xyz','r') as g1:
    Ptxyz = get_xyz(g1.readlines())
with open('28H.xyz','r') as g2:
    H28xyz = get_xyz(g2.readlines())
with open('29H.xyz','r') as g3:
    H29xyz = get_xyz(g3.readlines())

#print(Ptxyz)
godd_E = [-15.9,-10.09,-1.86,-0.56,-0.51,-0.13,0.7,1.4,2.34,2.05,1.4,0.52,0.]
godd_HF = [-16.21,-5.92,2.26,3.12,3.42,3.69,4.,4.04,3.47,2.37,1.26,0.15,0.]
godd_R = [1.166,1.417,1.542,1.667,1.792,1.917,2.042,2.167,2.417,2.667,2.917,3.167,4]
    
print(ircE)

PtH2xyz = list(np.linalg.norm(Ptxyz - 0.5*(H28xyz + H29xyz),axis=1))
print(PtH2xyz)
PtH2xyz = PtH2xyz[:20] + PtH2xyz[-50:]
print(PtH2xyz)
PtH2xyz2 = [1.206,4.0]
ircE2 = [-1042.434062,-1041.23404844-1.16787757399]
PtH2xyz += PtH2xyz2
ircE = ircE[:20] + ircE[-50:] + ircE2

"""
path1 = list(PtH2xyz[:51])
path1.reverse()
pathxyz = path1 + list(PtH2xyz[-51:])
print(pathxyz)
print(PtH2xyz) 
   """    
#print(ircE)
ircE = np.array(ircE)
ircE -= ircE[-1]
ircE *= H2KCAL
#print(ircE)
plt.figure()
plt.plot(PtH2xyz, ircE,'o')
plt.plot(godd_R, godd_E, '-*')
plt.plot(godd_R, godd_HF, '-*')
plt.legend(labels=['PBE0-D3/def2TZVP','GVB-RCI','HF'],fontsize=15)
plt.xlabel('R(Pt-H2)',fontsize=15)
plt.ylabel('E (kcal/mol)',fontsize=15)
plt.xticks(np.arange(1,5,0.5),('1.0','1.5', '2.0','2.5', '3.0','3.5', r'$\infty$'),fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
ax.invert_xaxis()
plt.savefig('irc.png',bbox_inches = 'tight')


