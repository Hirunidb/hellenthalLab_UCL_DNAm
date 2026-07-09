# hellenthalLab_UCL_DNAm
PhD Rotation II - A model for in-silico sorting of methylation reads from bulk-tissue by cell type

The core statistical model of the method follows 3 main steps:
  1. Comparing the similarity of the methylation pattern between a target read and the reference reads for a given region of the genome.
  2. Analysing the ability to distinguish different cell types within the reference reads of a given region
  3. Cumulation of step 1 and 2 to infer the probability of a cell type being the source of a given target read.

The simulated data used for the project originates from the methylation atlas published in 2023 by Loyfer et.at., where the blood samplese have been extracted. The analysis limits to the cell types with a minimum of 3 samples available, resulting in 27 samples belonging to 6 cell types. 

The supplimentary data provided by Loyfer et.al has been used to generate the colanno.csv file used throughout the project. Each row being a sample, and their metadata stored in the columns is as following, in order;
  1. File path for the .beta file
  2. Major cell type ('Refined group' as per supplimentary data)
  3. Sample name ('Sample name' as per supplimentary data)
  4. PatientID ('PatientID' as per supplimentary data)
  5. sampleID (last 9 character ID, extracted from sample name)
  6. stype ('Cell type' as per supplimentary data).
