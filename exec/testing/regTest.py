# Base class for regression tests
import sys
class regTest(object):

    def __init__(self, name):
        self.name = name
        self.modes = ['serial','mpi']
        self.compilers = ['gfortran']
        self.npes = [1,4]
        self.buildType = 'release'
        self.useBatch = 'no'
        self.modules = 'no'
        self.scratchDir = '.'
        self.resultsDir = '.'
        self.verification = 'compileOnly'
        self.endTime = 25
        self.unitTest = 'no'
        
    def getOpt(self, opt):
        if opt=='modes':
            return self.modes
        elif opt=='compilers':
            return self.compilers
        elif opt=='buildtype':
            return self.compile_only
        elif opt=='npes':
            nint = []
            for n in self.npes:
                nint.append(int(n))
            return nint
        elif opt=='usebatch':
            return self.useBatch
        elif opt=='modules':
            return self.modules
        elif opt=='resultsdir':
            return self.resultsDir
        elif opt=='scratchdir':
            return self.scratchDir
        elif opt=='verification':
            return self.verification
        elif opt=='endtime':
            return int(self.endTime)
        elif opt=='unittest':
            return self.unitTest

    def setOpts(self, deckconfig):
        for name,options in deckconfig.items():
            if self.name == name:
                for kk,vv in options.items():
                    if kk=='modes':
                        self.modes = vv
                    elif kk=='compilers':
                        self.compilers = vv
                    elif kk=='npes':
                        nint = []
                        for n in vv.split(','):
                            nint.append(int(n))
                        self.npes = nint 
                    elif kk=='buildtype':
                        self.compile_only = vv
                    elif kk=='usebatch':
                        self.useBatch = vv
                    elif kk=='modules':
                        self.modules = vv
                    elif kk=='scratchdir':
                        self.scratchDir = vv + '/scratch/'
                    elif kk=='resultsdir':
                        self.resultsDir = vv + '/results/'
                    elif kk=='verification':
                        self.verification = vv
                    elif kk=='endtime':
                        self.endTime = int(vv)
                    elif kk=='unittest':
                        self.unitTest = vv

    def dump(self):
        attrs = vars(self)
        print ', '.join("%s: %s" % item for item in attrs.items())

