from seisflows.plugins.io import sem
import numpy as np
from seisflows.tools import array

PATH_MODEL_INIT = "/home/ql5/Tests/marmousi_SRVM/model_true/"

x = sem.read(PATH_MODEL_INIT, 'x', 0)
z = sem.read(PATH_MODEL_INIT, 'z', 0)

print x[:5],x[-5:]

mesh = array.stack(x, z)

vp = sem.read(PATH_MODEL_INIT, 'vp', 0)
vs = sem.read(PATH_MODEL_INIT, 'vs', 0)

vp = array.meshsmooth(1.0/vp, mesh, 3)
vs = array.meshsmooth(1.0/vs, mesh, 3)

vp = np.float32(1.0/vp)
vs = np.float32(1.0/vs)

#np.save("proc000000_vp.bin", vp)
#np.save("proc000000_vs.bin", vs)

fname1 = "proc000000_vp.bin"
fname2 = "proc000000_vs.bin"

sem.write(vp, "./", 'vp', 0)
sem.write(vs, "./", 'vs', 0)

#vp.tofile(fname1)
#vs.tofile(fname2)

#f1 = open(fname1,"rb")
#xx1 = np.fromfile(f1,dtype=np.float32)
#print xx1[:5],xx1[-5:]





