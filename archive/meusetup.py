from os import getcwd, makedirs, listdir, chdir, popen
from os import path as ospath
from shutil import rmtree, move
from sys import path, argv, version
from platform import system
from platform import version as platVersion

OS = system()
Version = platVersion().split()[0]
PyVersion = version.split()[0]
home = getcwd()

print('Runing Setup for MoleKing_util on {} -{}- with python {}:'.format(OS, Version, PyVersion))

if argv[1] == 'bin':
    pyPath = 'bin'
    try:
        rmtree('bin')
    except:
        pass
    makedirs('bin')
else:
    for i in path:
        if 'site-packages' in i:
            for j in listdir(i):
                if 'pybind' in j:
                    pyPath = i
                    break
            

chdir(pyPath)
try:
    rmtree('MoleKing_util')
except:
    pass
chdir(home)

d = '{}/src'.format(home)
CList = ['{}/main.cpp'.format(d)]
subdirs = [ospath.join(d, o) for o in listdir(d) if ospath.isdir(ospath.join(d,o))]
for directory in subdirs:
    for arq in listdir(directory):
        if arq.split('.')[1] == 'cpp':
            CList.append('{}/{}'.format(directory, arq))

makedirs('MoleKing_util')

if OS == 'Linux':
    if len(PyVersion.split('.')[1]) == 1:
        flags = '-O3 -Wall -shared -std=c++11 -fPIC `python{} -m pybind11 --includes`'.format(PyVersion[0:3])
        target = 'MoleKing_util`python{}-config --extension-suffix`'.format(PyVersion[0:3])
    else:
        flags = '-O3 -Wall -shared -std=c++11 -fPIC `python{} -m pybind11 --includes`'.format(PyVersion[0:4])
        target = 'MoleKing_util`python{}-config --extension-suffix`'.format(PyVersion[0:4])

        
elif OS == 'Darwin':
    flags = '-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes`'
    target = 'MoleKing_util`python3-config --extension-suffix`'
else:
    print('MoleKing_util does not work on DOS base systems to this date.')

objs = ' '.join(CList)
chdir('{}/MoleKing_util'.format(home))
obj = popen('c++ {0} {2} -o {1}'.format(flags, target, objs), 'r')
obj.read()

pyarq = open('__init__.py', 'w')
pyarq.write('############ MoleKing_util constructor file ############\n')
pyarq.write('from MoleKing_util.MoleKing_util import Atom\n')
pyarq.write('from MoleKing_util.MoleKing_util import ChargePoint\n')
pyarq.write('from MoleKing_util.MoleKing_util import Molecule\n')
pyarq.write('from MoleKing_util.MoleKing_util import Point\n')
pyarq.write('from MoleKing_util.MoleKing_util import SphericalCoords\n')
pyarq.write('from MoleKing_util.MoleKing_util import Vector3D\n')
pyarq.write('from MoleKing_util.MoleKing_util import Matrix\n')
pyarq.write('from MoleKing_util.MoleKing_util import PeriodicTable\n')
pyarq.write('from MoleKing_util.MoleKing_util import SupraMolecule\n')
pyarq.write('from MoleKing_util.MoleKing_util import G16LOGfile\n')
pyarq.write('from MoleKing_util.MoleKing_util import G16FCHKfile\n')
pyarq.close()

chdir(home)
move(home+"/MoleKing_util", pyPath)
print('Success, MoleKing_util was compiled on {} -{}- with python {}!'.format(OS, Version, PyVersion))


