import env
import StructuralDescriptor
import characterization
import numpy as np

# reading data file
linker0 = env.readData.all('linker0', lmpFile='linker0.lmpdat', xyzFile='linker0-stretched.xyz', readCoe=True)
linker0.read()


#position = linker0.atom['xyzPosition']
linker0.atom['xyzPosition'] = StructuralDescriptor.orient.centering(linker0.atom['xyzPosition'],weight=linker0.atom['Mass'])
linker0.atom['xyzPosition'] = StructuralDescriptor.orient.orienting(linker0.atom['xyzPosition'])

moments = StructuralDescriptor.moment.moment(linker0.atom['lmpPosition'],weight=linker0.atom['Mass'],momType='mom',upto=5)
voxel, vBox = StructuralDescriptor.voxel.channel(linker0.atom['xyzPosition'],linker0.atom['Type'])
plt = StructuralDescriptor.voxel.plot(voxel, vBox)

pointnet = StructuralDescriptor.pointNet.pointNet(linker0.atom['xyzPosition'],linker0.atom['Type'])
env.writeFile.addNpy('moments.npy',np.insert(moments,0,0))
env.writeFile.addNpy('voxel.npy',np.insert(voxel,0,0))
env.writeFile.addNpy('pointnet.npy',np.insert(pointnet,0,0))

plt.show()

properties = characterization.FDcurve.calcFileFeatures(open('linker283-ave-force.d','r'))
l = [0,len(properties['jumpC']),len(properties['jumpT']),properties['aveCB'],properties['aveCF'],properties['aveTB'],properties['aveTF']]

env.writeFile.writeFile('property.txt',l)