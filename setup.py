from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension

extensions = [
    Extension("cpp_bind",
              ["cpp_bind.pyx"],
              language="c++",
              extra_compile_args=["-std=c++17"])
]

setup(
    name="pipe_model",
    ext_modules=cythonize(extensions)
)