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

import os,sys

#solv = ['benzene','dichloromethane']
#model = ['pcm','cpcm','smd']
#CleanMode = False
#SingleNode = False

jobname = sys.argv[1]

os.system("echo '100\n2\n8\noutput\n0\n-10\n' > nbo-M.in")
title = jobname[:-4] 
os.system("cp nbo-M.in nbo_%s-M.in"%(title))
os.system("sed -i 's/output/%s.47/g' %s/nbo_%s-M.in"%(title,title))
os.system("Multiwfn %s.fchk < nbo_%s-M.in > nbo_%s-M.out"%(title,title,title))
os.system("echo '#!/bin/bash\n$NBOBIN/gennbo.i8.exe %s.47 > %s-nbo.out' > %s-nbo.sh" % (title,title,title))
os.system("chmod 777 %s-nbo.sh" % (title))
os.system("bsub -n 24 -q mpi bash %s-nbo.sh" % (title))
#if CleanMode:
#    os.system("rm %s/*-M.in %s/*-M.out %s/*-nbo.sh"%(m,m,m))
