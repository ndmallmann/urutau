"""
    Install urutau package.
"""

from setuptools import setup, find_packages

setup(
    name="urutau",
    version="0.51",
    packages=find_packages(exclude=['main.py']),
    author="Nicolas Dullius Mallmann",
    python_requires=">=3.10",
)
