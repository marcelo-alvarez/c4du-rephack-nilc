from setuptools import setup, find_packages

setup(
    name="act-nilc",
    version="0.1.0",
    description="NILC pipeline implementation for ACT cosmology data",
    author="ACT Collaboration",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.24.0",
        "scipy>=1.10.0",
        "matplotlib>=3.7.0",
        "pixell>=0.24.0",
        "healpy>=1.16.0",
    ],
    extras_require={
        "test": [
            "pytest>=7.4.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
