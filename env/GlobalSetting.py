import numpy as np

#Parameters and format used in the calculations



#Coordination (space dimension)
_coorID_ = np.array([['x','1',0],['y','2',1],['z','3',2]],dtype='object')

#Output
_statusOutput_ = False #True
def __statusOutput__(status,localOutput = _statusOutput_):
    if _statusOutput_ or localOutput:
        print(status)

#Saveing Format
class _element_:
    def __init__(self, id = '', index = 0, mass = 0, radius = 0):
        self.index = index
        self.id = id
        self.mass = mass
        self.radius = radius
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
    class descriptor:
        moment = []
        gMoment = []
        voxel = []
        cVoxel = []
        smiles = ''
        pointNet = []
    class property:
        FDcurve = []
        stiffCurve = []
        stiff0 = float('-inf')
        stiff1 = float('-inf')
        stiff2 = float('-inf')
        diff = 0
        jump = []
    class other:
        other = {}

#System elements
_ELEMENTS_ = {'H': _element_(id='H',index=1,mass=1,radius=1.2),
              'C': _element_(id='C',index=6,mass=12,radius=1.7),
              'N': _element_(id='N',index=7,mass=14,radius=1.55),
              'O': _element_(id='O',index=8,mass=16,radius=1.52),
              'S': _element_(id='S',index=16,mass=32,radius=1.8),
              'Au': _element_(id='Au', index=79, mass=198, radius=1.66)}

#voxelization parameters
VOXESCALE = 0.1
THETA = 0.6
CUTOFF = 3
NUM_ELE = 4
voxPlot = False
voxChannel = True

#FD-curve analysis parameters
_fdDelta_ = 0.0001
_fdPoints_ = 1000
_jumpDepthThreshold_ = 5
_jumpWidthThreshold_ = 1
_jumpSlopeThreshold_ = 2
_aveFracLimit_ = 1000
_aveThreshold_ = 40
_fdPlot_ = False



