# OSCAR v1.4.0 Release Notes

We are pleased to announce the release of **OSCAR v1.4.0**. This version brings several improvements and new features:

- **Enhanced File Support:**  
  OSCAR now accepts both TIFF and TXT files. A new Fiji script is included to generate the necessary TXT files. Simply load your processed multichannel image in Fiji, and the script will produce the required TXT file, significantly increasing computational speed.

- **Optional Processing Pipeline:**  
  For single-channel TIFF files, OSCAR now includes an optional default processing pipeline.

## Repository Contents

This repository contains the programs and scripts used in the manuscript titled **"High-throughput three-dimensional characterization of morphogenetic signals during the formation of the vertebrate retina."** All code is written in Julia.

- **phantom.ipynb:**  
  Contains functions and scripts to generate digital images used as ground truth for comparing the performance of 3D segmentation solutions.

- **src folder:**  
  Includes all the necessary functions for the standalone OSCAR application, the Fiji script for generating TXT files, and a test image.

## Usage Instructions

1. **Initialize Julia:**  
   Ensure Julia is installed and initialized on your device.

2. **Include the OSCAR Script:**  
   Open a terminal and run:
   
   ```julia
   include("your/path/OSCAR.jl")
   ```

3. **Start OSCAR:**
   
   In the terminal, run:

   ```julia
   startOSCAR()
   ```

5. **Follow On-Screen Instructions:**
   
   Follow the prompts displayed in the terminal to view the results.
     
   
