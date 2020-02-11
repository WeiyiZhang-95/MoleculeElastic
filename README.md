# MoleculeElastic

This project works on researching about molecules' elastic properties. 

The aim is to 1) classify molecules with different dynamic bahavior by computational simulation 
and 2) gain capability to predict/fast-characterize molecules kinetic features.

## Contents

- [Characterization](characterization/chara.md): => [:point_right::open_file_folder::point_left:](characterization)

    - [LAMMPS](https://lammps.sandia.gov) simulation script example
    - Stiffness/Deforming-related analysis from output files of LAMMPS simulation

- [env](env/env.md): => [:point_right::open_file_folder::point_left:](env)
    
    *_Change parameters from here to make it suitable for your system_*
    - Useful tools for reading required data from LAMMPS input files or .xyz formatted files
    - General Setting calculation setups/parameters
    - General info of used molecular system

- [StructuralDescriptor](StructuralDescriptor/descriptor.md): => [:point_right::open_file_folder::point_left:](StructuralDescriptor)

    *_Different ways to represent molecular structure numerically_*
    - SMILE-based Encoding
    - Voxelization
    - PointNet
    - High-order Moment
    - ... ...

## Manual/[Example](example)

- Modules:
  ```python
    import env
    import StructuralDescriptor
    import characterization
    import numpy as np
  ```
- Reading data file
  ```python
    linker0 = env.readData.all('linker0', lmpFile='linker0.lmpdat', xyzFile='linker0-stretched.xyz', readCoe=True)
    linker0.read()
    #position = linker0.atom['xyzPosition']
    linker0.atom['xyzPosition'] = StructuralDescriptor.orient.centering(linker0.atom['xyzPosition'],weight=linker0.atom['Mass'])
    linker0.atom['xyzPosition'] = StructuralDescriptor.orient.orienting(linker0.atom['xyzPosition'])
  ```
- Generating Descriptors
  ```python
    moments = StructuralDescriptor.moment.moment(linker0.atom['lmpPosition'],weight=linker0.atom['Mass'],momType='mom',upto=5)
    voxel, vBox = StructuralDescriptor.voxel.channel(linker0.atom['xyzPosition'],linker0.atom['Type'])
    plt = StructuralDescriptor.voxel.plot(voxel, vBox)
    plt.show()
    
    pointnet = StructuralDescriptor.pointNet.pointNet(linker0.atom['xyzPosition'],linker0.atom['Type'])
    env.writeFile.addNpy('moments.npy',np.insert(moments,0,0))
    env.writeFile.addNpy('voxel.npy',np.insert(voxel,0,0))
    env.writeFile.addNpy('pointnet.npy',np.insert(pointnet,0,0))
  ```

- Generating Property features
  ```python
    properties = characterization.FDcurve.calcFileFeatures(open('linker283-ave-force.d','r'))
    l = [0,len(properties['jumpC']),len(properties['jumpT']),properties['aveCB'],properties['aveCF'],properties['aveTB'],properties['aveTF']]
    env.writeFile.writeFile('property.txt',l)
  ```
