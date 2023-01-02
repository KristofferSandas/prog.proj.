README


INTRODUCTION

visualizeKraken2 produces graphs and visualizations of the taxonomic composition of
metagenomic analyses performed in Kraken 2. visualizeKraken2 can produce a circular 
sunburst plot, an icicle chart or a ranked list of species/taxons. The plots approximate 
abundance of the taxons by displaying the percent count for each taxon, meaning the 
proportion of the total counts in the analysis belonging to each specific taxon. 
The count for a taxon is the number of sequences in the original multifasta file 
assigned to that specific taxon by Kraken 2. 


KRAKEN 2 OUTPUT FILES

Kraken 2 reads a mulitfasta file from a metagenomic sample and produces two types of files, 
a report variety with a .kreport extension, and a text file with a .kraken2 extension. 
The input for visualizeKraken2 is the text file variety. This file type contains columns 
with information about the analyzed sample: which sequences were successfully classified,
and which taxonomy IDs they map to. Kraken 2 aims to find the Lowest Common Ancestor (LCA) 
among all candidates for each sequence, meaning it tries to go as far down the taxonomy 
tree as possible. Sometimes, however, there is not enough information to pinpoint species 
or genus, and the algorithm simply adds the lowest taxonomy rank that the information in 
the sequence allows for. The result is that the ID assigned by Kraken is the Lowest Common 
Ancestor of all the candidate taxons mapped to a certain sequence. 


VISUALIZEKRAKEN2

visualizeKraken2 maps these IDs to the NCBI taxonomy database, extracting the name of 
the LCAs and their full taxonomic hierarchy or lineage. It also counts how many times a 
certain taxon, or LCA, was found by Kraken 2 and uses this information to produce plots 
where an approximation of taxon abundance can be seen. 

The root node of the NCBI taxonomy tree is called "root" and Kraken 2 often assigns this
ID if there is not enough information in a sequence to determine anything about the 
organism. Since a "root" designation tells us nothing more than that the sequence contains
nucleotides, and does not even differentiate between prokaryotes and eukaryotes, 
visualizeKraken2 removes these entries and displays a message to the user that a "root"
entry has been removed and the number of counts assigned to this entry.

If an entry has been assigned an ID that does not exist in the NCBI database, 
visualizeKraken2 outputs a message to the user that the ID was not found, after which it
is removed from the analysis.


The following python packages are used by the application and need to be installed
and functioning properly:

	numpy
	pandas
	plotly
	time 
	os
	shutil
	ftplib
	zipfile
	re

The files beloning to this application are:

	- This readme file
	- visualizeKraken2.ipynb: the main jupyter file that the program is run from
	- startupFunctions.py: contains the functions needed to download and build the 
	  taxonomy database
	- krakViewClass.py: the main class containing the functions needed to perform the 
	  visualizations
	- four example files called 'creek_water.kraken2', 'endotracheal.kraken2', 
	  'fecal.kraken2' and 'raw_sewage.kraken2'. These are output files from Kraken 2
	  analyses and can be used to familiarize yourself with the application.

These files should be placed in the same folder for the program to run properly.The program 
is designed to be run as a jupyter notebook. The main jupyter script is a template for the 
various functions and methods that are included in the program. Any number of analyses can 
be run in the same script, alternatively copies of the jupyter file can be made and modified
to make notebooks dedicated to a single analysis. Any file to be run in visualizeKraken2
should also be placed in this same directory, unless the full path is specified when 
loading the file.


The main functions are:

downloadDatabase(dbLocalDir, NCBIdbFile)

	This function is imported from startupFunctions.py and only needs to be run the first 
	time the program is used on a new machine, or if the NCBI database file is removed
	or needs to be updated or re-downloaded for some reason. 
	
	dbLocalDir is the directory where the NCBI database files will be downloaded. This 
	can be changed if another directory is desired. 
	
	NCBIdbFile is the name of the correct file in the NCBI server. This should always
	be "new_taxdump.zip". This file is updated hourly on the NCBI server. 
	
	The most straightforward way of running this function is to specify the directory and
	file name as variables and then run the function exactly as above.


buildDatabase(dbLocalDir)
	
	This function is also imported from startupFunctions.py and builds a database from 
	the downloaded NCBI file in the form of a python dictionary. This needs to be run
	at the beginning of every session to create the database. It should take about 20 
	seconds to run. Once the database is built, any number of .kraken2 files can be 
	analyzed in the session.

	dbLocalDir is the folder where the NCBI file used to build the database is located.
	It needs to be the same as in the downloadDatabase function.

	This function should be run saving the database as a variable:
		dataBase = buildDatabase(dbLocalDir)


analysisObject = krakView("filename.kraken2", pyDict = dataBase)

	This function is imported from krakViewClass.py and creates an object which can be 
	used to produce the various plots. The variable name is the name you wish to give
	the analysis, the filename is the .kraken2 file you want to display. 
	
	pyDict = dataBase specifies the name of the dictionary dataBase used for the analysis.
	This is the dataBase built using the buildDatabase function, so it needs to
	have the same name as given when built.
	
	
Methods that can be used on a krakView object:

.sourceFile
	
	Displays the name of the file used for the analysis object.

.data
	
	Displays the raw data in the file.

.showFormat(n = 5)

	Displays the first n lines of the data in the file. Default is n = 5.

.percentClassified(output = 1)

	Calculates the percentage of successfully classified entries in the kraken2 file.
	If this number is very low, there could have been problems with the sample or 
	something could have gone wrong in the kraken analysis process.
	Output=1 displays a descriptive message with the percentage. Output=0 disables 
	the message if the number is simply being stored in a variable for example.

.numberUniqueIDs(output = 1)

	Returns the total number of unique taxonomy IDs found in the kraken2 file. Output=1
	displays a descriptive message, output=0 disables it.

.countsDistribution()
	
	Outputs a data frame with each count and the number of taxons that have that count.
	This can be used to get a sense of how to filter the data when producing plots.
	
.showList(n = 25 , filterCounts = 5)

	This produces a list of the n taxons with the highest number of counts. The list 
	shows the number of counts for each taxon, and the percent counts, meaning the
	proportion of counts for the specific taxon compared to the total number of counts.

	filterCounts is the specified minimum number of counts needed for a taxon to be 
	included in the analysis. The distribution of counts can be seen using the 
	.countDistribution method. The default value for filtering is 5. 

	The messages "x unique taxons found" and "y entries in total after adding parents
	to all" mean that filtered by a specific count, x number of unique taxons were found
	and that after visualizeKraken2 added parents to all these, there are y entries in 
	the analysis list. 


.showSunburst(filterCounts= 5, width= 1100, height= 700, title = "Sunburst chart")

	Produces a sunburst chart of the resluts, categorized by percent counts, 
	meaning the number of counts corresponding to a certain taxon compared to the 
	total number of counts. 
	
	filterCounts is the specified minimum number of counts needed for a taxon to be 
	included in the analysis. The distribution of counts can be seen using the 
	.countDistribution method. The default value for filtering is 5. Some filtering 
	is recommended,	since the plots will be quite crowded if every taxon with a single 
	count is included. Plotly can also run into probelms if the sets to display are too
	large. 
	.showList can be used for an unfiltered overview.
	
	width, height can be used to change the size of the sunburst chart.
	
	title can be used to give the chart a specific title

	The messages "x unique taxons found" and "y entries in total after adding parents
	to all" mean that filtered by a specific count, x number of unique taxons were found
	and that after visualizeKraken2 added parents to all these, there are y entries in 
	the list used to make the plot. 

.showIcicle(filterCounts= 5, width= 1100, height= 700, title = "Icicle chart")

	Produces an icicle chart of the resluts, categorized by percent counts, 
	meaning the number of counts corresponding to a certain taxon compared to the 
	total number of counts. 
	
	filterCounts is the specified minimum number of counts needed for a taxon to be 
	included in the analysis. The distribution of counts can be seen using the 
	.countDistribution method. The default value for filtering is 5. Some filtering 
	is recommended,	since the plots will be quite crowded if every taxon with a single 
	count is included. Plotly can also run into probelms if the sets to display are too
	large. 
	.showList can be used for an unfiltered overview.
	
	width, height can be used to change the size of the chart.
	
	title can be used to give the chart a specific title
	
	The messages "x unique taxons found" and "y entries in total after adding parents
	to all" mean that filtered by a specific count, x number of unique taxons were found
	and that after visualizeKraken2 added parents to all these, there are y entries in 
	the list used to make the plot. 
	