"""
This is the main class containing all the functions needed to run the analysis. It is imported into the main script
and run from there.
"""

import numpy as np
import pandas as pd
import plotly.express as px

class krakView:
    """
    Defines a class for visualizing and displaying data from an imported kraken2 txt file.
    Input file must be in the same folder as the script, unless the entire path is provided.
    """
    
    def __init__ (self, krakenFile, pyDict): 
        """
        Creates an object from a kraken2 txt file which can be used to produce the various plots. 
        Automatically removes unclassified entries. 
        
        pyDict = dataBase specifies the name of the dictionary dataBase used for the analysis.
        This is the dataBase built using the buildDatabase function, so it needs to have the same name 
        as given when built.
        """
        self.sourceFile = krakenFile
        self.data = np.loadtxt(krakenFile, dtype=str, delimiter='\t')
        uniqueIDs, count = np.unique(self.data[:,2], return_counts=True)
        self.count = count[1:] 
        self.uniqueIDs= uniqueIDs[1:] # removes the ID 0, which represents unclassified entries
        self.dataBase = pyDict # for some reason, only this variable name works
        
    def showFormat(self, n = 5):
        """
        Displays the first n lines of the data in the file.
        """
        return self.data[0:n, :]
    
    def percentClassified(self, output = 1):
        """
        Calculates the percentage of successfully classified entries in the kraken2 file.
        If this number is very low, there could have been problems with the sample or something could 
        have gone wrong in the kraken analysis process.
        
        Output=1 displays a descriptive message with the percentage. Output=0 disables the message if 
        the number is simply being stored in a variable for example.
        """
        Us = np.sum(np.char.count(self.data[:, 0], 'U'))
        Cs = np.sum(np.char.count(self.data[:, 0], 'C'))
        percentageC = float(round((Cs/(Cs+Us))*100,1))
        if output == 1:
            print('The percentage of classified sequences in the fasta file: ' + str(percentageC) + '%')
        return percentageC
    
    def numberUniqueIDs(self, output = 1):
        """
        Returns the total number of unique taxonomy IDs found in the kraken2 file. Output=1	displays a 
        descriptive message, output=0 disables it.
        """
        if output == 1:
            print('Number of unique taxIDs found: ' + str(len(self.uniqueIDs)))
        return len(self.uniqueIDs)
    
    def countsDistribution(self):
        """
        Outputs a data frame with each count and the number of taxons that have that count.
        This can be used to get a sense of how to filter the data when producing plots.	
        """
        countsList, nrOfInstances = np.unique(self.count, return_counts=True)
        countsList = list(countsList)
        nrOfInstances = list(nrOfInstances)
        countDict = {}
        countDict['Count'] = countsList
        countDict['Nr of instances'] = nrOfInstances
        countDistribution = pd.DataFrame(countDict)
        countDistributionTransposed = countDistribution.transpose()
        pd.set_option('display.max_columns', 500)
        return countDistributionTransposed
    
    def mapIDs(self, filterCounts = 5):
        """
        Maps the IDs in the kraken file to the LCA names and ancestors in the NCBI database.
        
        filterCounts specifies filtering by number of counts found for a specific taxon,
        for example filterCounts = 5 maps only the taxons with more than 5 appearances in the kraken file.
        Too many entries may cause problems with plotly visualization, so filtering by 5 or 10 is 
        recommended. 
        
        "root" means that Kraken 2 could not place it more precisely than at the root of the NCBI taxonomy tree, 
        which basically tells us nothing except that the sequence contained nucleotides. These entries are removed.
        
        Commented out prints can be used for troubleshooting.   
        """
        speciesCounts = [] # number of counts for each ID
        speciesNames = [] # names of the LCAs
        speciesLineage = [] # full lineage in a list
        for i in range(len(self.uniqueIDs)): # takes each taxID at a time
            #print('\ni=', i, '\n')
            if self.count[i] > filterCounts: # Filter by number of matches.
                #print('i=', i)
                if self.uniqueIDs[i] == '1': # removes the 'root' entries
                    print('"' + str(self.dataBase[self.uniqueIDs[i]][0]) + '" ' + 'with ' + 
                          str(self.count[i]) + ' counts removed from list.')
                else:
                    try: # parse the dict database and add name and lineage list for each ID
                        #print(self.dataBase[self.uniqueIDs[i]][0])
                        speciesNames.append(self.dataBase[self.uniqueIDs[i]][0])
                        reversedAncestors = list(self.dataBase[self.uniqueIDs[i]][1])
                        reversedAncestors.reverse() # reversing to make parsing in ascending order easier
                        speciesLineage.append(reversedAncestors)
                        speciesCounts.append(self.count[i]) 
                    except: # if ID is not in database
                        #print('i=', i)
                        print('ID', self.uniqueIDs[i], 'not in database')
            
        if len(speciesCounts) == len(speciesNames) == len(speciesLineage): # plotly will not work if not even lenghts
            print(len(speciesLineage), 'unique taxons found.')
        else:
            print('Something went horribly wrong with database mapping. Visualization will probably not work :(')
        
        # This next part arranges the data in the format required by plotly to produce the visualizations.
        
        char = [] # creates a list of characters, each character is an entry in the plot
        parent = [] # the closest parent of the enrty in char. Every object in parent must also be in char.
        match = [] # the number of counts for the taxon in question.
        
        # when one of these root nodes is encountered, the plotly root parent '' is appended 
        listOfRootNodes = ['cellular organisms', 
                          'Viruses', 
                          'other entries', 
                          'unclassified entries']

        for i in range(len(speciesNames)): # take one LCA at a time
            #print('\n---------- entry', i,  '--------')
            #print('Name:', speciesNames[i], '\n')
            #print('Lineage list:', speciesLineage[i])
            if speciesNames[i] in listOfRootNodes: # if the LCA is a root node, add it and the '' parent
                char.append(speciesNames[i])
                match.append(speciesCounts[i])        
                parent.append('')
                #print('i=', i, '\nroot APPENDED as parent to', speciesNames[i], '\n')
            else: # add the LCA, its closest parent, and its count
                char.append(speciesNames[i])
                match.append(speciesCounts[i])
                parent.append(speciesLineage[i][0])
                #print('i=', i, speciesLineage[i][0], 'appended as PARENT to', speciesNames[i])
                for j in range(len(speciesLineage[i])): # check if the parent is already registered. If not, it's added w/ count=0.
                    #print('Current ancestor:', speciesLineage[i][j])
                    if speciesLineage[i][j] not in speciesNames and speciesLineage[i][j] not in char:
                        if speciesLineage[i][j] in listOfRootNodes: # this handles when the parents reach the root nodes
                            #print('i=', i)
                            char.append(speciesLineage[i][j])
                            #print(speciesLineage[i][j], 'appended to char')
                            match.append(0)
                            parent.append('')
                            #print('root appended as parent to', speciesLineage[i][j])
                        else: # this adds the current parent with a count of 0
                            char.append(speciesLineage[i][j])
                            #print(speciesLineage[i][j], 'appended to char')
                            match.append(0)
                            parent.append(speciesLineage[i][j+1])
        
        if len(char) == len(parent) == len(match):
            print(len(char), 'entries in total after adding parents to all.\n')
        else:
            print('Something went horribly wrong with parent mapping. Visualization will probably not work :(')
        
        return char, parent, match
    
    def organizeData(self, char, parent, match):
        """
        Organizes the data for plotly visualization.
        """
        percentMatches = []
        total = sum(self.count) 
        for i in match: # calculates relative abundance
            percent = round((i/total)*100, 2)
            percentMatches.append(percent)
            
        plotData = dict(LCA = char,
                        Parent = parent,
                        PercentCounts = percentMatches)
        return plotData
    
    def showList(self, n=25, filterCounts = 5):
        """
        This produces a list of the n taxons with the highest number of counts. The list shows 
        the number of counts for each taxon, and the percent counts, meaning the proportion of 
        counts for the specific taxon compared to the total number of counts.

    	filterCounts is the specified minimum number of counts needed for a taxon to be included 
        in the analysis. The distribution of counts can be seen using the .countDistribution method. 
        The default value for filtering is 5. 

    	The messages "x unique taxons found" and "y entries in total after adding parents to all" 
        mean that filtered by a specific count, x number of unique taxons were found and that after 
        visualizeKraken2 added parents to all these, there are y entries in the analysis list.
        """
        char, parent, match = self.mapIDs(filterCounts)
        percentMatches = []
        total = sum(self.count)
        for i in match:
            percent = round((i/total)*100, 2)
            percentMatches.append(percent)
            
        matchesAndLca = pd.DataFrame({'Percent counts' : percentMatches, 'Counts' : match, 'LCA' : char})
        matchesAndLcaSorted = matchesAndLca.sort_values('Percent counts', ascending=False)
        matchesAndLcaSortedTopResults = matchesAndLcaSorted.head(n=n)
        print('List of the top', n, 'taxons with highest number of counts:\n')
        print(matchesAndLcaSortedTopResults.to_string(index=False))
    
    def showSunburst(self, filterCounts = 5, width = 1100, height = 700, title = 'Sunburst chart'):
        """
        Produces a sunburst chart of the resluts, categorized by percent counts, meaning 
        the number of counts corresponding to a certain taxon compared to the total number of counts. 
	
       	filterCounts is the specified minimum number of counts needed for a taxon to be included 
        in the analysis. The distribution of counts can be seen using the .countDistribution method. 
        The default value for filtering is 5. Some filtering is recommended, since the plots will 
        be quite crowded if every taxon with a single count is included. Plotly can also run into 
        probelms if the sets to display are too large. 
        .showList can be used for an unfiltered overview.
	
    	width, height can be used to change the size of the sunburst chart.
	
    	title can be used to give the chart a specific title

    	The messages "x unique taxons found" and "y entries in total after adding parents to
        all" mean that filtered by a specific count, x number of unique taxons were found and 
        that after visualizeKraken2 added parents to all these, there are y entries in the 
        list used to make the plot. 
        """
        char, parent, match = self.mapIDs(filterCounts)
        plotData = self.organizeData(char, parent, match)
        fig = px.sunburst(plotData,
                          names = 'LCA',
                          parents = 'Parent',
                          values = 'PercentCounts',
                          color = 'PercentCounts',
                          title = title,
                          width = width,
                          height = height,
        )

        fig.show()
           
    def showIcicle(self, filterCounts = 5, width = 900, height = 2000, title = 'Icicle chart'):
        """
        Produces an icicle chart of the resluts, categorized by percent counts, 
        meaning the number of counts corresponding to a certain taxon compared to the 
        total number of counts. 
	
    	filterCounts is the specified minimum number of counts needed for a taxon to be 
		included in the analysis. The distribution of counts can be seen using the 
		.countDistribution method. The default value for filtering is 5. Some filtering is recommended,
		since the plots will be quite crowded if every taxon with a single count is 
		included. Plotly can also run into probelms if the sets to display are too large. 
        .showList can be used for an unfiltered overview.
	
    	width, height can be used to change the size of the chart.
	
    	title can be used to give the chart a specific title
	
    	The messages "x unique taxons found" and "y entries in total after adding parents to
        all" mean that filtered by a specific count, x number of unique taxons were found and
        that after visualizeKraken2 added parents to all these, there are y entries in the
        list used to make the plot. 	
        """
        char, parent, match = self.mapIDs(filterCounts)
        plotData = self.organizeData(char, parent, match)
        fig = px.icicle(plotData,
                        names ='LCA',
                        parents ='Parent',
                        values = 'PercentCounts',
                        color = 'PercentCounts',
                        title = title,
                        width = width,
                        height = height,
        )
        fig.update_traces(root_color="lightgrey")
        fig.update_layout(margin = dict(t=50, l=25, r=25, b=25)) # might need fine-tuning depending on machine
        fig.show()
