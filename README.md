# MoleKing
MoleKing is a Python module written in C++ with pybind11 Linkage under LEEDMOL Research Group. This module contains several useful classes for those who program python scripts aimed at theoretical chemistry. The main goal of this package is to introduce concepts of chemistry, such as Molecules, Atoms and Geometries, to python, making programing more intuitive and understandable to chemists. Additionally, MoleKing is capable of read and write inputs and outputs files for several theoretical chemistry programs·
#########################################################################################################################################################
#Requirements:
C++11 or greater;
Cmake 3.15 or greater;
libEingen;
libBoost;
pybind11;
PIP10 or greater.
#########################################################################################################################################################
#Installation:
Clone to repository ->
                  	For python installation pip3 install ./MoleKing under the main project directory
                  	For Debugging and Testing mkdir build -> cd build -> cmake ../ -DBuild_Python=OFF -> cmake --build .
                  	For python testing and debugging mkdir build -> cd build -> cmake ../ -DBuild_Python=ON -> cmake --build .
#########################################################################################################################################################
                          !! Please feel free to report any issues or suggestions to mateus_barbosa@ufg.discente.br !!
