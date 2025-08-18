from setuptools import setup, find_packages

setup(
    name="simba",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pandas",
        "matplotlib",
        "scipy"
    ],
    entry_points={
        "console_scripts": [
            "simba=simba.cli:main"
        ]
    },
    author="Milton Simbarashe Kambarami",
    description="Codon-aware selection analysis and mutation hotspot mapping",
)

