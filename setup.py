import os

from setuptools import find_packages, setup


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "hybracter",
            "hybracter.VERSION",
        )
    ) as f:
        return f.readline().strip()


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

# read long description
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="hybracter",
    packages=find_packages(),
    url="https://github.com/gbouras13/hybracter",
    python_requires=">=3.9",
    description="An automated long-read first bacterial genome assembly pipeline.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=get_version(),
    author="George Bouras",
    author_email="george.bouras@adelaide.edu.au",
    py_modules=["hybracter"],
    install_requires=[
        "snakemake>=8.0",
        "pyyaml>=6.0",
        "Click>=8.1.3",
        "biopython>=1.76",
        "pandas>=1.5.0",
    ],
    extras_require={
        # deps needed to run the fast unit tests (`pip install .[test]`)
        "test": ["pytest"],
    },
    entry_points={
        "console_scripts": [
            "hybracter=hybracter.__main__:main",
            "concat-gfas=hybracter.concat_GFAs_hybracter:main",
        ]
    },
    include_package_data=True,
)
