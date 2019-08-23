import numpy as np

#coordination di m
_coorID_ = np.array([['x','1',0],['y','2',1],['z','3',2]],dtype='object')

class _element_:
    def __init__(self, id = '', index = 0, mass = 0):
        self.index = index
        self.id = id
        self.mass = mass

H = _element_(id='H',index=1,mass=1)
C = _element_(id='C',index=6,mass=12)
N = _element_(id='N',index=7,mass=14)
O = _element_(id='O',index=8,mass=16)
S = _element_(id='S',index=16,mass=32)
Au = _element_(id='Au',index=79,mass=198)

_ELEMENTS_ = [H,C,N,O,Au]

VOXESCALE = 0.5
THETA = 0.6
CUTOFF = 3
NUM_ELE = 4

_statusOutput_ = False #True
def __statusOutput__(status,localOutput = False):
    if _statusOutput_ or localOutput:
        print(status)


class _atom_:
    id = ''
    index = 0
    mass = 0

class _mol_:
    class basic:
        id = ''
        mass = ''
        aPosition = []
        aSpecies = []
        aMass = []
        more = {}
    class descriptor:
        moment = []
        gMoment = []
        voxel = []
        cVoxel = []
        smiles = ''
        pointNet = []
        more = {}
    class property:
        FDcurve = []
        stiffCurve = []
        stiff0 = float('-inf')
        stiff1 = float('-inf')
        stiff2 = float('-inf')
        diff = 0
        jump = []
        more = {}

