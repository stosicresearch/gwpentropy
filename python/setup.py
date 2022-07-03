import torch
import setuptools
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

setup(
    name='gwpe',
    ext_modules=[
        CUDAExtension(
            name='gwpe',
            sources=['gwpentropy.cpp'],
            extra_compile_args={'cxx':['-O3'],
                                'nvcc':['-O3','--use_fast_math']})
    ],
    cmdclass={
        'build_ext': BuildExtension
})
