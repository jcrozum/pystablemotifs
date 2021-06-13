from setuptools import setup

setup(
    name="StableMotifs",
    version='1.0.0',
    author="Jordan Rozum",
    author_email="jcr52@psu.edu",
    description="Python package for analyzing Boolean Netowrk",
    url="https://github.com/jcrozum/StableMotifs",
    license='MIT',
    python_requires='>=3.5',
    packages=['StableMotifs'],
    install_requires=[
    'PyBoolNet @ git+https://github.com/hklarner/PyBoolNet@2.3.0',
    "networkx >= 2.4.0",
    "sympy >= 1.5.1",
    "pandas >= 1.0.0",
    "numpy >= 1.19.2",
    "matplotlib >= 3.2.1"
    ]
)
