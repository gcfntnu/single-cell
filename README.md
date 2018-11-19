# single-cell
Single cell RNA seq analysis, 10Xgenomics platform

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/single-cell.svg)](https://hub.docker.com/r/single-cell)

### Introduction
single-cell: Single cell RNA seq analysis, 10Xgenomics platform

The pipeline is built using [Snakemake](https://bitbucket.org/snakemake/snakemake), a flexible pipeline tool to run tasks across multiple compute infrastructures. It comes with docker containers making installation trivial and results highly reproducible.


### Documentation
The single-cell pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline was written by Arnar Flatberg ([flatberg](https://github.com/flatberg)) at [Norwegian University of Science and Technology, NTNU](http://www.ntnu.no).
