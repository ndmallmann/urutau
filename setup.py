"""
    Install urutau package.
"""

from setuptools import setup, find_packages

setup(
    name="urutau",
    version="1.2",
    packages=find_packages(),
    package_data={
        "urutau_gui": ["assets/*"],
    },
    include_package_data=True,
    author="Nicolas Dullius Mallmann & Rogerio Riffel",
    python_requires=">=3.10",
    install_requires=["astropy", "pandas", "scipy"],
    extras_require={
        "gui": ["PyQt5"],
    },
    entry_points={
        "console_scripts": [
            "fit_analyser=utils.fit_analyser:main",
            "urutau-gui=urutau_gui.main_gui:main",
        ],
    },
)
