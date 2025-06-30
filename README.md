# OSCAR v1.6 Release Notes

We are pleased to announce the release of **OSCAR v1.6**. This version brings several improvements and new features:

- **Enhanced File Support:**  
  OSCAR now accepts both TIFF and TXT files. A new Fiji script is included to generate the necessary TXT files. Simply load your processed multichannel image in Fiji, and the script will produce the required TXT file, significantly increasing computational speed.

- **Optional Processing Pipeline:**  
  For single-channel TIFF files, OSCAR now includes an optional default processing pipeline.

## Repository Contents

This repository contains the programs and scripts used in the manuscript titled **"Object Stitching by Clustering of Adjacent Regions designed for accurate quantification in dense three-dimensional tissues."** All code is written in Julia.

- **phantom.ipynb:**  
  Contains functions and scripts to generate digital images used as ground truth for comparing the performance of 3D segmentation solutions.

- **src folder:**  
  Includes all the necessary functions for the standalone OSCAR application, the Fiji script for generating TXT files, and a test image.

## Usage Instructions

1. **Initialize Julia:**  
   Ensure Julia is installed (https://julialang.org/downloads/) and initialized on your device.
   ```console
   julia -t auto
   ```

3. **Include the OSCAR Script:**  
   Open a terminal and run:
   
   ```julia
   include("your/path/OSCARvX.X.X.jl")
   ```

4. **Start OSCAR:**  
   In the terminal, run:

   ```julia
   startOSCAR()
   ```

5. **Follow On-Screen Instructions:**  
   Follow the prompts displayed in the terminal to view the results.
