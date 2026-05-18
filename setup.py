from setuptools import setup, find_packages

setup(
    name="mabsislandscanner",
    version="1.0.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "mabs-scan=mabsislandscanner.scanner:main",
        ],
    },
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.21",
    ],
)
