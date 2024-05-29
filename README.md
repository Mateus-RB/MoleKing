# This is a working in progress at LEEDMOL-UFG under manitence of Mateus Rodrigues Barbosa, Pedro Henrique F. Matias and Rafael F. Veríssimo. Feel free to use as a learning tool but don't trust the results.
# MoleKing: A Python Module for Theoretical Chemistry.
MoleKing is a Python module written in C++ with pybind11 Linkage under [LEEDMOL](leedmol.com) Research Group. This module contains several useful classes for those who program python scripts aimed at theoretical chemistry. This package's main goal is to introduce chemistry concepts, such as Molecules, Atoms, and Geometries, to python, making programming more intuitive and understandable to chemists. Additionally, MoleKing is capable of reading and writing inputs and outputs files for several theoretical chemistry programs.

---

# Requirements:

   <ul>
   <li> C++11 or greater.</li>
    <li>Cmake 3.15 or greater;</li>
    <li>libEingen;</li>
    <li>pybind11;</li>
    <li>PIP10 or greater</li>
    </ul>

---

# Installation:

Clone the repository 
                  	
For **python** installation 
```pip3 install ./MoleKing ``` under the main project directory;

For **Debugging** and **Testing**:

```
mkdir build
cd build
cmake ../ -DBuild_Python=OFF
cmake --build . 
```

For **python testing** and **debugging**

```
mkdir build
cd build
cmake ../ -DBuild_Python=ON
cmake --build . 
```

---

**Please feel free to report any issues or suggestions to mateus_barbosa@ufg.br or phfmatias@discente.ufg.br**

![MoleKing_Grad](https://user-images.githubusercontent.com/71854729/213286170-38170b42-8e1b-4bfb-9b9b-80aa8308444e.png)


