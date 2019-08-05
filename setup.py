import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rnaevol-lemerleau",
    version="0.0.3",
    author="Nono Saha Cyrille Merleau",
    author_email="nonosaha@mis.mpg.de",
    description="Evolving a population of rna sequence for an efficient rna design",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/strevol-mpi-mis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
