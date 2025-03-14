# Release notes OSCAR v1.4.0

New OSCAR version available. It now accepts either tiff files and txt files. For the appropiate generation of txt files we have included a Fiji script, charge your own processed multichannel image in Fiji and obtain the needed txt for running OSCAR. This method greatly increases computational speed. On the other hand for single-channel tiff files OSCAR now includes an optional default processing pipeline.

###

Please find in ths repository the relevalt programs and scripts used in the mansucript entited "High-throughput three-dimensional characterization of morphogenetic signals during the formation of the vertebrate retina"
All code is written in julia 
- phantom.ipynb includes teh functions and scripts to generate the digital images uses as groudn thruth to comapre the performance of 3D segmentation solutions.
- src folder includes all necessary functions for OSCAR standalone application, Fiji script to obtain txt files and a test image.

# USAGE
1. Initialize Julia in your device
2. In terminal, run:
       include("your/path/OSCAR.jl")
4. In terminal, run:
    startOSCAR()
5. Follow the instructions and see the results in terminal
