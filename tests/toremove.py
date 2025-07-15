from MoleKing import G16LOGfile
from os import getcwd

def test_getAlpha():
    f = G16LOGfile('MK_polar.log', polarAsw=True).getFrequency()
    a1 = G16LOGfile('MK_polar.log', polarAsw=True).getAlpha(unit='esu', frequency=f[0])['iso']
    a2 = G16LOGfile('MK_polar.log', polarAsw=True).getAlpha(unit='SI', frequency=f[1])['xx']


def test_getBeta():
    f = G16LOGfile('MK_polar.log', polarAsw=True).getFrequency()
    print(f[1])
    #print(G16LOGfile('MK_polar.log', polarAsw=True).getBeta(unit='esu', frequency=f[1], orientation='input', BSHG=0))
    binp = G16LOGfile('MK_polar.log', polarAsw=True).getBeta(unit='esu', frequency=f[1], orientation='input', BSHG=0)['||(z)']
    bdip = G16LOGfile('MK_polar.log', polarAsw=True).getBeta(unit='esu', frequency=f[1], orientation='dipole', BSHG=0)['||(z)']
    b2 = G16LOGfile('MK_polar.log', polarAsw=True).getBeta(unit='au', frequency=f[1], BSHG=1)['_|_(z)']
    print(binp)
    print('-----------------------------------------')
    print(bdip)

def test_getGamma(sef):
    f = G16LOGfile('MK_polar.log', polarAsw=True).getFrequency()
    g1 = G16LOGfile('MK_polar.log', polarAsw=True).getGamma(unit='esu', frequency=f[0])['||']
    g2 = G16LOGfile('MK_polar.log', polarAsw=True).getGamma(unit='au', frequency=f[1], GSHG=1)['xxxx']


if __name__ == '__main__':
    test_getBeta()