"""
    Install urutau package.
"""

from setuptools import setup, find_packages

setup(
    name="urutau",
    version="1.002",
    packages=find_packages(),
    author="Nicolas Dullius Mallmann",
    python_requires=">=3.10",
    install_requires=["astropy", "pandas", "scipy"]
)
