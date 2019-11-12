import os

solv = ['benzene','dichloromethane']
model = ['pcm','cpcm','smd']

for s in solv:
    for m in model:
        os.system("cp ts_def2tzvp_scrf.gjf ts_def2tzvp_scrf_%s_%s.gjf"%(m,s))
        os.system("sed -i '3 s/$/ scrf=(%s,solvent=%s)/' ts_def2tzvp_scrf_%s_%s.gjf"%(m,s,m,s))
        os.system("sed -i '1 s/ts/ts_%s_%s/' ts_def2tzvp_scrf_%s_%s.gjf"%(m,s,m,s))
        os.system("bsub -n 24 -q largemem g16 ts_def2tzvp_scrf_%s_%s.gjf"%(m,s))
