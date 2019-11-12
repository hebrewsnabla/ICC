"""
ATTENTION:
Add your path to Multiwfn and NBO6 to your .bashrc before running, e.g.
  export PATH=/share/home/wli/srwang/Multiwfn_3.5_bin_Linux:$PATH
  export Multiwfnpath=/share/home/wli/srwang/Multiwfn_3.5_bin_Linux
  export NBOHOME=/share/home/wli/srwang/nbo6
  export NBOBIN=$NBOHOME/bin
  export NBOEXE=$NBOBIN/nbo6.i8.exe
Modify 2 parameters in /.../Multiwfn_3.5_bin_Linux/setting.ini as
  isilent=1
  iloadasCart=1
if they are 0 in your Multiwfn version.

"""

import os

solv = ['benzene','dichloromethane']
model = ['pcm','cpcm','smd']
CleanMode = False
SingleNode = False

os.system("echo '100\n2\n8\noutput\n0\n-10\n' > nbo-M.in")
for s in solv:
    for m in model:
        title = "%s_%s" % (m,s)
        os.system("cp nbo-M.in %s/nbo_%s-M.in"%(m,title))
        os.system("sed -i 's/output/%s\/ts_%s.47/g' %s/nbo_%s-M.in"%(m,title,m,title))
        #os.system("sed -i '1 s/ts/ts_%s_%s/' ts_def2tzvp_scrf_%s_%s.gjf"%(m,s,m,s))
        os.system("Multiwfn %s/ts_%s_maugccpvtz.fchk < %s/nbo_%s-M.in > %s/nbo_%s-M.out"%(m,title,m,title,m,title))
        #os.system("echo '#!/bin/bash\nnbo6 %s/ts_%s.47 > %s/ts_%s-nbo.out' > %s/ts_%s-nbo.sh" % (m,title,m,title,m,title))
        os.system("echo '#!/bin/bash\n$NBOBIN/gennbo.i8.exe %s/ts_%s.47 > %s/ts_%s-nbo.out' > %s/ts_%s-nbo.sh" % (m,title,m,title,m,title))
        os.system("chmod 777 %s/ts_%s-nbo.sh" % (m,title))
        os.system("bsub -n 24 -q mpi ./%s/ts_%s-nbo.sh" % (m,title))
        #if CleanMode:
        #    os.system("rm %s/*-M.in %s/*-M.out %s/*-nbo.sh"%(m,m,m))
