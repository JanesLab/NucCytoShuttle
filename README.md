Introduction to NucCytoShuttle

NucCytoShuttle is a differential equations model of nucleocytoplasmic transport adapted from Riddick & Macara (2005, 2007) as described in Wang et al. (2023).  It is comprised of two .m files, a function that parses inputs and returns outputs (NucCytoShuttleModel.m), and a second function that encodes the system of differential equations (NucCytoShuttleODE.m).  The model relies on an ancillary file (ModelParametersDefault.xlsx), which contains the rate parameters used in the simulation.  The model requires a fixed “Spike” concentration (in µM) to be simulated for six categories of cargo:  1) SV40 T antigen monopartite NLS (“NLS”), 2) importin-beta binding domain of KPNA2 (“IBB”), 3) bipartite NLS of Cbp80 (“CBP80”), 4) NLS + nuclear export sequence (NES) of PKIA (“NNES”), 5) IBB + NES (“INES”), 6) CBP80 + NES (“CNES”).

Initial conditions are organized into four categories:  1) “CellType” — geometries that specify the nuclear and cytoplasmic volumes of the cell (in pL), 2) “CellConcs” — initial whole-cell concentrations of the ten bulk species in the model (Ran, RanBP1, RanGap, RCC1, NTF2, Impalpha, Impbeta, Cas, CRM1, RanBP3), 3) “LumpedReceptors” — bulk concentrations of four endogenous cargo types that compete with spike cargo (alpha-beta cargo, beta cargo, other transport receptors, NES cargo), 4)”NUPs” – nuclear pores classified by those that carry all cargo and those that only carry endogenous cargo (i.e., not spike cargo).  Default initial conditions for HeLa cells are used if arguments are not passed during the function call.

The function outputs a table of seven tables, each corresponding to the time evolution of species in the model for the six cargo, plus an initial cargo-free simulation that was used to establish the steady state before the cargo spike.  A second table outputs the calculated cytoplasmic volume and nuclear volume after accounting for the NPC compartment that creates a peri-nuclear annulus around a spherical nucleus.


System requirements and installation

NucCytoShuttle requires MATLAB (R2022a or later) v9.12 without any further dependencies.  It has been tested on MacOS v12.6 and Windows 10 Pro v22H2.  No non-standard hardware is required.  To install, download all files into a common folder and specify the common folder as the active folder or add the common folder to the MATLAB path.  Typical install time on a normal desktop computer should be instantaneous.


Example use cases:

%Case 1: Simulate 1 µM spike and return the results to output table “ot” and cell volume table “cvt”
[ot,cvt] = NucCytoShuttleModel(1);
(Run time on an Apple M1 Pro with 32 GB memory:  18 seconds.  Saved as NucCytoShuttleCase1.mat)

%Case 2: Replace with the cell geometries of B2B1 cells
[ot,cvt] = NucCytoShuttleModel(1,’CellType’,[0.52 1.45]);
(Run time on an Apple M1 Pro with 32 GB memory:  17 seconds.  Saved as NucCytoShuttleCase2.mat)

%Case 3: Replace with species concentrations of B2B1 cells
[ot,cvt] = NucCytoShuttleModel(1,’CellConcs’,[5.18; 2.69; 0.85; 0.45; 0.81; 1.52; 2.77; 5.26; 0.50; 0.15],’LumpedReceptors’,[10; 5; 6; 1]);
(Run time on an Apple M1 Pro with 32 GB memory:  20 seconds.  Saved as NucCytoShuttleCase3.mat)

%Case 4: Replace with NUPs of B2B1 cells
[ot,cvt] = NucCytoShuttleModel(1,’NUPs’,2600);
(Run time on an Apple M1 Pro with 32 GB memory:  16 seconds.  Saved as NucCytoShuttleCase4.mat)

%Case 5: Split the number of NUPs into separate classes
[ot,cvt] = NucCytoShuttleModel(1,’NUPs’,[1300; 1300]);
(Run time on an Apple M1 Pro with 32 GB memory:  29 seconds.  Saved as NucCytoShuttleCase5.mat)
