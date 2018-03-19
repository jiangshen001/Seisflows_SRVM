import sys
import numpy as np

from seisflows.tools import unix
from seisflows.tools.array import loadnpy, savenpy
from seisflows.tools.tools import savetxt, exists
from seisflows.tools.math import angle, dot

PATH = sys.modules['seisflows_paths']

class SRVM(object):
    """ Square root variable metric algorithm

        Includes optional safeguards: periodic restarting and descent
        conditions.

        To conserve memory, most vectors are read from disk rather than 
        passed from a calling routine.
    """

    def __init__(self, path='.', load=loadnpy, save=savenpy, memory=5, thresh=0., maxiter=np.inf, precond=None):
        assert exists(path)
        unix.cd(path)
        unix.mkdir('SRVM')

        self.path = path
        self.load = load
        self.save = save
        self.thresh = thresh
        self.maxiter = maxiter
        self.precond = precond
        self.memory = memory

        self.iter = 0
        self.memory_used = 0


    def __call__(self):
        """ Returns SRVM search direction
        """
        self.iter += 1

        unix.cd(self.path)
        g = self.load('g_new')
        if self.iter == 1:
	    ShatT_ghat = g
	    self.save('ShatT_ghat',g)
	    print 'max(ShatT_ghat)',max(ShatT_ghat)
            return -g, 0

        elif self.iter > self.maxiter:
            print 'restarting SRVM... [periodic restart]'
            self.restart()
            return -g, 1

	mu = self.loadtxt('alpha')
	g_old = self.load('g_old')


	print 'max(g_old-g)',max(g_old-g)

	dghat = g - g_old
	yhat = mu * g_old + dghat

	kk = self.iter - 1

	w = self.update(yhat,kk-1,1)

	self.save('w',w)
	unix.mv('w','w_%04d' % kk)
	print 'max(w)',max(w)
        ShatT_ghat = self.load('ShatT_ghat')
	#print 'max(ShatT_ghat)',max(ShatT_ghat)
	belta = w - mu * ShatT_ghat
	a = dot(w,belta)
	b = dot(w,w)
	nu = self.srvm_nu(a,b)

	self.savetxt('a',a)
	unix.mv('a','a_%04d' % kk)

	self.savetxt('nu',nu)
	unix.mv('nu','nu_%04d' % kk)
	
	ShatT_ghat = self.update(g,kk,1)
	q = self.update(ShatT_ghat,kk,0)
	#q = -q
	
	self.save('ShatT_ghat',ShatT_ghat)

	print 'max(g)',max(g)
	print 'max(q)',max(q)

        status = self.check_status(g,q)
	print 'status',status
	#status = 0
        if status != 0:
            self.restart()
            return -g, status
        else:
            return -q, status


    def update(self,chi,kk,tflag):
        """ Updates SRVM algorithm history
        """
        unix.cd(self.path)
		
	Shat_chi = chi

	if tflag == 1:
            for ii in range(kk):
		jj = ii + 1

		unix.cp('a_%04d' % jj,'A')
		A = self.loadtxt('A')

		unix.cp('nu_%04d' % jj,'Nu')
		nu = self.loadtxt('Nu')

		unix.cp('w_%04d' % jj,'W')
	        wtemp = self.load('W')

	        xtemp = dot(wtemp,Shat_chi)
		Shat_chi = Shat_chi - xtemp * nu / A *wtemp

	elif tflag == 0:
	    Shat_chi = chi
            for ii in range(kk):
		jj = kk - ii

		unix.cp('a_%04d' % jj,'A')
		A = self.loadtxt('A')

		unix.cp('nu_%04d' % jj,'Nu')
		nu = self.loadtxt('Nu')

		unix.cp('w_%04d' % jj,'W')
	        wtemp = self.load('W')

	        xtemp = dot(wtemp,Shat_chi)
		Shat_chi = Shat_chi - xtemp * nu / A *wtemp

	return	Shat_chi
			
	 		
    def srvm_nu(self,a,b):
        """ Determine nu
        """
	ratio = b / a

	if ratio < 1:
	    nu = (1 -np.sqrt(1 - ratio)) / ratio
	elif ratio == 1:
	    nu = 0.9
	else:
	    nu = 1.0

	return nu
	    
		

    def restart(self):
        """ Discards history and resets counters
        """



    def check_status(self, g, r):
        theta = 180.*np.pi**-1*angle(g,r)
	print 'theta =', theta
        if not 0. < theta < 90.:
            print 'restarting SRVM... [not a descent direction]'
            return 1
        elif theta > 90. - self.thresh:
            print 'restarting SRVM... [practical safeguard]'
            return 1
        else:
            return 0


    def dot(self,x,y):
        """ Computes inner product between vectors
        """
        return np.dot(
            np.squeeze(x),
            np.squeeze(y))

    def load(self, filename):
        return loadnpy(PATH.OPTIMIZE+'/'+filename)

    def save(self, filename, array):
        savenpy(PATH.OPTIMIZE+'/'+filename, array)


    def loadtxt(self, filename):
        return float(np.loadtxt(filename))

    def savetxt(self, filename, scalar):
        np.savetxt(PATH.OPTIMIZE+'/'+filename, [scalar], '%11.6e')


