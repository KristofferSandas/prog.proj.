"""
This file contains the two startup functions needed for visualizeKraken2. The functions 
are imported into the main script and run from there. The first function downloads the 
file needed for the database from the NCBI server. The second function builds the data-
base as a python dictionary.
"""
import time
import os
import shutil
import ftplib
import zipfile
import re

version = 'beta0.8'

def downloadDatabase(dbLocalDir, NCBIdbFile):
    """
    If the correct NCBI file is not in its place locally, this function removes the whole database 
    folder and downloads everything from scratch from the NCBI server. The remove is done
    to avoid problems with old or incomplete files from previous downloads. If the correct file
    is in its place, the function will do nothing but display a message to the user that the needed
    file is already present.
    
    Malware and trackers on the computer might cause corruption of downloaded files. If 
    unexplainable errors occur in this step, try cleaning the computer of trackers and other
    things that might interfere with the download.
    """
    start = time.time()
    print('\nvisualizeKraken2 v.' + version)
    print('------------------------------\n')
    print('Checking for database source file...')
    
    if not os.path.isfile(str(dbLocalDir) + '/fullnamelineage.dmp'): # this is the file used to build the database
        if os.path.isfile(NCBIdbFile):
            os.remove(NCBIdbFile) # this removes any old zip file that might cause problems.
        if os.path.exists(dbLocalDir):
            shutil.rmtree(dbLocalDir) # os module stops you because of windows admin rights. shutil works.
        os.mkdir(dbLocalDir)
        print('\ndirectory "' + str(dbLocalDir) + '" created\n')
        with ftplib.FTP('ftp.ncbi.nlm.nih.gov') as ftp: # connecting to the NCBI server.
            print('Downloading taxonomy database from NCBI...\n')
            try:
                ftp.login()
                ftp.cwd('pub/taxonomy/new_taxdump/') # navigating to the correct directory in the server.
                with open(NCBIdbFile, "wb") as localFile:
                    ftp.retrbinary("RETR " + NCBIdbFile, 
                                   localFile.write, 8*1024) # download in binary to prevent file corruption
                    print('binary download successful')
            
            except ftplib.all_errors as error:
                print('FTP error: ',  error)
                print("This means something has gone wrong with the server. \
                      Make sure you're connected to the internet and try again.")
            
        with zipfile.ZipFile(NCBIdbFile, 'r') as zipFile: # unzips file
            print('unzipping')
            zipFile.extractall(dbLocalDir)
        
        os.remove(NCBIdbFile) # removes the zip file and any downloaded files which are not needed.
        for file in os.listdir(dbLocalDir):
            if file != 'fullnamelineage.dmp':
                #print(str(databaseFolder) + '/' + str(file) + ' removed') # can be used for troubleshooting.
                os.remove(str(dbLocalDir) + '/' + str(file))
        print('\nDatabase downloaded.\n')
        end = time.time()
        print('Time elapsed:', str(round(end-start, 1)), 'seconds\n')

    else:
        print('\nDatabase source file exists already. No download required.\n')
        
def buildDatabase(dbLocalDir):
    """
    Builds a python dictionary from the source file downloaded from NCBI.
    This is the most time-consuming step in the program and so is done before any analysis is started.
    Since the program is designed to be run as a jupyter notebook, this step only needs to be
    done once in each session. After the database is built, any number of files can be analysed
    quickly and efficiently. 
    
    Avoid doing other things on the computer while building the database
    as this process is prone to hanging itself if too many things are happening at once. It should
    only take about 20 seconds to build the database.
    """
    start = time.time()
    dataBase = {}
    print('Building database...\n')
    with open(str(dbLocalDir) + '/fullnamelineage.dmp', 
              errors='ignore') as inputF: # decoding the .dmp file = problematic. Errors are now ignored. Downstream problems?
        for line in inputF:
            tempList1 = re.split(r'\t\|\t', line)
            tempList2 = re.split(r'; ', tempList1[2])
            tempList3 = []
            tempList3.append(tempList1[1])
            tempList3.append(tempList2[:-1])
            dataBase[tempList1[0]] = tempList3
            #for key,item in dataBase.items(): # can be used for troubleshooting
                #print(key, item)
    
    print('Database built.\n')
    print('Length of dictionary "dataBase":',len(dataBase))
    end = time.time()
    print('Time to build database:', str(round(end-start, 1)), 'seconds\n')
    return dataBase