
import sys
import numpy as np

from seisflows.config import custom_import, ParameterError
from seisflows.optimize.lib.SRVM import SRVM as lib
from seisflows.tools import unix

PAR = sys.modules['seisflows_parameters']
PATH = sys.modules['seisflows_paths']


class SRVM(custom_import('optimize', 'base')):
    """ Square-root variable metric algorithm
    """

    def check(self):
        """ Checks parameters, paths, and dependencies
        """
        # line search algorithm
        if 'LINESEARCH' not in PAR:
            setattr(PAR, 'LINESEARCH', 'Backtrack')

        # SRVM memory
        if 'SRVMMEM' not in PAR:
            setattr(PAR, 'SRVMMEM', 3)

        # SRVM periodic restart interval
        if 'SRVMMAX' not in PAR:
            setattr(PAR, 'SRVMMAX', np.inf)

        # SRVM angle restart threshold
        if 'SRVMTHRESH' not in PAR:
            setattr(PAR, 'SRVMTHRESH', 0)

        super(SRVM, self).check()


    def setup(self):
        super(SRVM, self).setup()

        self.SRVM = lib(
            path=PATH.OPTIMIZE,
            memory=PAR.SRVMMEM,
            maxiter=PAR.SRVMMAX,
            thresh=PAR.SRVMTHRESH,
            precond=self.precond())


    def compute_direction(self):
        g_new = self.load('g_new')
        p_new, self.restarted = self.SRVM()
        self.save('p_new', p_new)
        self.savetxt('s_new', self.dot(g_new, p_new))
        return p_new


    def update_SRVM(self):
        self.path=PATH.OPTIMIZE
        unix.cd(self.path)

	print 'here'
	if self.iter > 1:

	    g_old = self.load('g_old')
	    p_old = self.load('p_old')
	    g = self.load('g_new')

	    g *=1.0e5
            self.save('g_new',g)

	    mu = self.loadtxt('alpha')
	    #if self.iter > 2:
	    #    mu = -self.dot(g_old,p_old)/self.dot(p_old,p_old)

	    if mu > 2:
		mu = 2.0

	    if mu < 1.0e-4:
		mu = 1.0e-4

	    print 'mu =', mu
	    dghat = -g + g_old
	    yhat = mu * g_old - dghat

	    kk = self.iter - 1
	    w = self.update_w(yhat,kk-1)
	    self.save('w',w)
	    unix.mv('w','w_%04d' % kk)

            ShatT_ghat = self.load('ShatT_ghat')
	    belta = - w + mu * ShatT_ghat
	    a = self.dot(w,belta)
	    b = self.dot(w,w)
	    nu = self.srvm_nu(a,b)
	    print 'a,b,ratio,nu', a,b,b/a,nu

	    self.savetxt('a',a)
	    unix.mv('a','a_%04d' % kk)

	    self.savetxt('nu',nu)
	    unix.mv('nu','nu_%04d' % kk)


    #def update_coef(self):
	
	
    def srvm_nu(self,a,b):
        """ Determine nu
        """
	ratio = b / a
	#print 'ratio =',ratio

	if ratio > 1:
	    #nu = (1 -np.sqrt(abs(1 + ratio))) / ratio
	    nu = -1
	elif abs(ratio) < 1.0e-8:
	    nu = 1.0
	else:
	    nu = -(1 -np.sqrt(1 - ratio)) / ratio
	return nu


    def update_w(self,chi,kk):
        """ Updates SRVM algorithm history
        """

        self.path=PATH.OPTIMIZE
        unix.cd(self.path)
		
	mm = 5

	Shat_chi = chi

        for ii in range(mm):
	    jj = ii + 1 + kk - mm

	    if jj > 0 :

    	        unix.cp('a_%04d' % jj,'A')
	        a = self.loadtxt('A')

		unix.cp('nu_%04d' % jj,'Nu')
	    	nu = self.loadtxt('Nu')

		unix.cp('w_%04d' % jj,'W')
	        wtemp = self.load('W')

	        xtemp = self.dot(wtemp,Shat_chi)
		Shat_chi = Shat_chi - xtemp * nu / a * wtemp


	return	Shat_chi


    def restart(self):
        super(SRVM, self).restart()
        self.SRVM.restart()

