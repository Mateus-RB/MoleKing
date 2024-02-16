from os import chdir, listdir

#Getting current version

chdir('../')

arq = open('CMakeLists.txt','r').read()

version = arq.split('MoleKing VERSION ')[1].split('\n')[0].split(')')[0]

chdir('src')
folders = [x for x in listdir() if '.' not in x]

for i in folders:
    chdir(i)
    files = [x for x in listdir() if '.cpp' in x or '.hpp' in x]
    for j in files:
        k = open(j,'r').readlines()
        
        for l in range(len(k)):
            if '//   Version:     ' in k[l]:
                k[l] = "//   Version:     ['{}']\n".format(version)
    
        m = open(j,'w')
        for l in k:
            m.write(l)
        m.close()

    chdir('..')

for i in ['main.cpp', 'main2.cpp']:
    k = open(i,'r').readlines()
    
    for l in range(len(k)):
        if '//   Version:     ' in k[l]:
            k[l] = "//   Version:     ['{}']\n".format(version)
    
    m = open(i,'w')
    for l in k:
        m.write(l)
    m.close()
