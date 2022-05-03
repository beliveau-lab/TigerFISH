.. Tigerfish documentation master file, created by
   sphinx-quickstart on Thu Jan 20 15:33:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Getting Started
---------------

**Tigerfish** is a user-friendly and interactive software tool that designs oligonucleotide DNA probes for FISH experiments that target specific satellite DNA and repetitive genomic regions at the scale of genomes.

Citing
------
If you use Tigerfish in your work, please cite:

.. epigraph::
   Coming soon to biorxiv!

Introduction
------------

As emerging genomes are fully assembled, methods and resources must be developed to understand how previously unmappable and highly repetitive DNA sequences contribute to genome organiazation and stability. Here, **Tigerfish** provides a flexible framework that allows users to design oligo probes against highly repetitive DNA sequences in genomes for FISH experiments. With **Tigerfish**, only a handful or as little as one probe can be used to target highly repetitive DNA sequences with robust target specificity. **Tigerfish** includes several features that are accessible to users who are interested in planning repetitive DNA FISH experiments. 

Allows for de-novo repeat discovery and probe design |:microscope:|
        A pipeline that allows for repetitive probe discovery and design. Users may provide a particular region of interest for oligo probe design, or provide a genomics scaffold for repeat discovery and probe design.

Provides downstream oligo probe analyses |:bar_chart:|
        A post processing pipeline, where users may perform analysis on probes of interest to determine on and off target probe specificity in silico and generate karyoplots using ChromoMaps.

An upcoming web resource for model organisms |:tropical_fish:|
        An interactive web application for probe design that will include diverse model organism genomes containing information about relevant repetitive DNA probes for FISH experiments (FISHtank, forthcoming Q3 2022).

**Tigerfish** is intended to serve as a resource for exploring and visualizing abundantly repetitive DNA targets in genomes using oligo technology. You can find out more about features and functionality in these pages. Happy FISHing! |:fish:|


Installation
------------

1. Install `conda <https://docs.conda.io/en/latest/miniconda.html>`_ as needed for your system.

2. Proceed with installing Mamba to assist with snakemake installation, as recommended in the Snakemake installation tutorial.

3. Clone the **Tigerfish** repo and create and active the provided environment for implementing Snakemake for Tigerfish:
