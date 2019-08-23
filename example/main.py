from env import readData
from StructuralDescriptor import moment
import numpy as np
import env
linker0 = readData.all('linker0', lmpFile='linker0.lmpdat', xyzFile='linker0-stretched.xyz', readCoe=True)
linker0.read()
print(linker0.coeff)

#print(moment.momentCalc(linker0.atom['lmpPosition'],weight=linker0.atom['Mass'],index=12312))
moments = moment.moment(linker0.atom['lmpPosition'],weight=linker0.atom['Mass'],index=['111',12212,'xyyxy'])
print(moments)