
import sys
import numpy as np

from seisflows.config import custom_import, ParameterError
from seisflows.optimize.lib.SRVM import SRVM as lib

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
            setattr(PAR, 'LINESEARCH', 'Bracket')

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


    def restart(self):
        super(SRVM, self).restart()
        self.SRVM.restart()

