from __future__ import annotations

import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# Map setuptools' Windows platform identifiers to the arguments CMake expects.
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}

ROOT_DIR = Path(__file__).parent.resolve()
CMAKELISTS = ROOT_DIR / "CMakeLists.txt"


class CMakeExtension(Extension):
    """A setuptools Extension that is built via an out-of-tree CMake project."""

    def __init__(self, name: str, sourcedir: str | Path = ROOT_DIR) -> None:
        super().__init__(name, sources=[])
        self.sourcedir = str(Path(sourcedir).resolve())


class CMakeBuild(build_ext):
    """Invoke CMake to build pybind11 extensions."""

    def build_extension(self, ext: CMakeExtension) -> None:
        ext_path = Path(self.get_ext_fullpath(ext.name))
        extdir = ext_path.parent.resolve()

        debug = bool(self.debug or int(os.environ.get("DEBUG", 0)))
        cfg = "Debug" if debug else "Release"

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DBuild_Python=ON",
        ]

        env_args = os.environ.get("CMAKE_ARGS")
        if env_args:
            cmake_args += shlex.split(env_args)

        build_args: list[str] = []
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        if self.compiler.compiler_type == "msvc":
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]
        else:
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja  # type: ignore
                except ImportError:
                    pass
                else:
                    ninja_path = Path(ninja.BIN_DIR) / "ninja"  # type: ignore[attr-defined]
                    cmake_args += [
                        "-GNinja",
                        f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_path}",
                    ]

        if sys.platform.startswith("darwin"):
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += [f"-DCMAKE_OSX_ARCHITECTURES={';'.join(archs)}"]

        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            if getattr(self, "parallel", None):
                build_args += [f"-j{self.parallel}"]

        build_temp = Path(self.build_temp) / ext.name
        build_temp.mkdir(parents=True, exist_ok=True)

        subprocess.check_call(
            ["cmake", ext.sourcedir, *cmake_args],
            cwd=build_temp,
        )
        subprocess.check_call(
            ["cmake", "--build", ".", *build_args],
            cwd=build_temp,
        )


def get_version_from_cmakelists(path: Path) -> str:
    pattern = re.compile(r"project\s*\(\s*MoleKing\s+VERSION\s+([^\s)]+)", re.IGNORECASE)
    with path.open(encoding="utf-8") as cmake_file:
        for line in cmake_file:
            match = pattern.search(line)
            if match:
                return match.group(1)
    raise RuntimeError("MoleKing VERSION not found in CMakeLists.txt")

setup(
    name="MoleKing",
    version = get_version_from_cmakelists('CMakeLists.txt'),
    author="LEEDMOL Research Group",
    author_email="mateus_barbosa@ufg.br",
    description="MoleKing is a python module for chemists aiming to add common principles to python. This module adds new types of python variables, MoleKing_Molecule; MoleKing_Atom; MoleKing_SupraMolecule, and MoleKing_Output, alongside many features considered common knowledge among chemists.",
    long_description="MoleKing is a Python module written in C++ with pybind11 Linkage under LEEDMOL Research Group. This module contains several useful classes for those who program python scripts aimed at theoretical chemistry. This package's main goal is to introduce chemistry concepts, such as Molecules, Atoms, and Geometries, to python, making programming more intuitive and understandable to chemists. Additionally, MoleKing is capable of reading and writing inputs and outputs files for several theoretical chemistry programs.",
    ext_modules=[CMakeExtension("MoleKing")],
    py_modules=["MoleKing"],
    include_package_data=True,
    package_data={"":["tests/*.log", "tests/*.out"]},
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.9",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
   )