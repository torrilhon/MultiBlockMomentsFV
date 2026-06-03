This directory contains the following.

Files for the moment-FV-solver:
 Main.c               - main routines for moment solver
 Flux.c               - flux routines for moment solver
 Flux.h               - header for flux routines
 Input.c              - input data for moment simulation
 Input.h              - header for input data
 Boundary.c           - boundary procedures for moment solver
 database_range.c     - grid management 
 database_range.h     - header grid management 
 tools.c              - some usefull tools
 tools.c              - tool header
 defs.h               - fundamental definitions for moment solver
 primDef.h            - code-piece for primitive variable definition

Run-Time Files:
 blockfile.txt        - block file containing a square grid with 4 blocks
 blockCorner.txt      - block file containing an L-shape domain with 3 blocks

Reference Files:
 output.ref.txt               - text output of the standard simulation run
 Grad13Riemann2D.ref.dat      - tecplot file containing standard solution at 6 time levels
 TemperatureRiemann2D.ref.jpg - jpg-image of a contour plot of the standard solution

Rest: 
 CodeOveriew.pdf  - a pdf containing a printable collection of the most important code parts
 Readme.txt       - this file

------------
Instructions
------------

- Compile the program with the command given on top of Mainmoment.c
- Run the program with a.out (it uses blockfile.txt) 
- Compare the text output to the screen with output.ref.txt 
- Compare the data output in the file Grad13Riemann2D_nst40000_Kn0p001_tAll.dat
  with Grad13Riemann2D.ref.dat. This file is tecplot compatible.
- Compare a contour plot of temperature obtained from your simulation 
  with the image TemperatureRiemann2D.ref.jpg
- Try running a.out -bfile blockCorner.txt for another example






