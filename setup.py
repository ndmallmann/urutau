"""
    Install urutau package.
"""

from setuptools import setup, find_packages

setup(
    name="urutau",
    version="1.001",
    packages=find_packages(),
    author="Nicolas Dullius Mallmann & Rogerio Riffel",
    python_requires=">=3.10",
    install_requires=["astropy", "pandas", "scipy"]
)
