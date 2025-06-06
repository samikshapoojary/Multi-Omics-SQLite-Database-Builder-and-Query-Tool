#importing all necessary modules
import sys #to create command line arguements
import sqlite3 #to create and connect to database
import re #to find match in metabolite name and replace it

#to create dataframe and plotting
import pandas as pd 
from matplotlib import pyplot as plt 
import seaborn as sns

#obtaining all the required files from working directory
subject_file = 'Subject.csv' #file containing information of the subjects

#files containing entity ids, sample ids and their abundance values
transcriptome_file = 'HMP_transcriptome_abundance.tsv' 
proteome_file = 'HMP_proteome_abundance.tsv'
metabolome_file = 'HMP_metabolome_abundance.tsv'

metabolome_annot_file = 'HMP_metabolome_annotation.csv' #file containing annotation of each metabolome

#creating a function to create the database
def create_db(database):
    #writing out the sql commands to create the database 
    create_sql = '''--creating the Subject table
        CREATE TABLE Subject (
            SubjectID TEXT PRIMARY KEY,
            Race TEXT,
            Sex TEXT,
            Age REAL,
            BMI INT,
            SSPG_Level INT,
            Insulin_Status TEXT
        );
        
        --creating the Visits table
        CREATE TABLE Visits (
            VisitID TEXT,
            SubjectID TEXT,
            FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID)
        );
        
        --creating the Sample table
        CREATE TABLE Sample (
            SampleID TEXT PRIMARY KEY,
            SubjectID TEXT,
            VisitID TEXT,
            FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID),
            FOREIGN KEY (VisitID) REFERENCES Visits(VisitID)
        );
        
        --creating the Transcriptome table
        CREATE TABLE Transcriptome (
            SampleID TEXT,
            TranscriptomeID TEXT,
            Abundance REAL,
            FOREIGN KEY (SampleID) REFERENCES Sample(SampleID)
        );
        
        --creating the Proteome table
        CREATE TABLE Proteome (
            SampleID TEXT,
            ProteomeID TEXT,
            Abundance REAL,
            FOREIGN KEY (SampleID) REFERENCES Sample(SampleID)
        );
        
        --creating the Peak (Metabolome) table
        CREATE TABLE Peak (
            SampleID TEXT,
            PeakID TEXT,
            Abundance REAL,
            FOREIGN KEY (SampleID) REFERENCES Sample(SampleID)
        );
        
        --creating the Metabolome Annotation table
        CREATE TABLE MetabolomeAnnotation (
            PeakID TEXT,
            Metabolite_Name TEXT,
            KEGG_ID TEXT,
            HMDB_ID TEXT,
            Chemical_Class TEXT,
            Pathway TEXT,
            FOREIGN KEY (PeakID) REFERENCES Peak(PeakID)
        );'''

    try:
        conn = sqlite3.connect(database) #connecting to the database
        cur = conn.cursor()
        cur.executescript(create_sql) #executing above sql command to create database
        conn.commit() #committing the changes 
        conn.close() #closing the connection

    except Exception as error: #incase of an error, informing the user
        print(f'Database could not be created due to : {error}.')
        cur.close()
        conn.close() #closing the connection

#creating multiple functions to load data into each table

#creating a function to load Subject table
def load_subject_table(cur, file_path): #insert path of subjects data file
    subject_sql = '''INSERT INTO Subject (SubjectID, Race, Sex, Age, BMI, SSPG_Level, Insulin_Status)
    VALUES (?, ?, ?, ?, ?, ?, ?)''' #sql statement to load Subject table

    try:
        with open(file_path) as sub: #opening the csv file containing subjects data
            header = sub.readline() #skipping the header line
            for line in sub:
                SubjectID,Race,Sex,Age,BMI,SSPG_Level,Insulin_Status= line.strip().split(',') #splitting each row by ',' into its contents
                if Age == "NA": #incase of age being NA, storing it as None
                    Age = None
                else: #incase of age being a continuous value, storing as float
                    try:
                        Age = float(Age)  
                    except ValueError:
                        continue  #skipping if the value is not a number
                if BMI == "NA": #incase of bmi being NA, storing it as None
                    BMI = None
                else: #incase of bmi being a continuous value, storing as float
                    try:
                        BMI = float(BMI)  
                    except ValueError:
                        continue  #skipping if the value is not a number
                if SSPG_Level == "NA": #incase of sspg level being NA, storing it as None
                    SSPG_Level = None
                else: #incase of sspg level being a continuous value, storing as float
                    try:
                        SSPG_Level = float(SSPG_Level)  
                    except ValueError:
                        continue  #skipping if the value is not a number
                if Race == 'unknown': #incase of race being unknown, storing it as None
                    Race = None
                if Insulin_Status == 'Unknown':  #incase of insulin status being Unknown, storing it as None
                    Insulin_Status = None
                #inserting data into the Subject table using the sql statement above
                cur.execute(subject_sql, (SubjectID,Race,Sex,Age,BMI,SSPG_Level,Insulin_Status))

    except FileNotFoundError: #incase of an error in file name or path, informing the user
        print(f'{file_path} was not found.')
    except Exception as error: #incase of an error loading, informing the user
        print(f'Data could not be loaded into Subjects table due to : {error}')
        cur.close() #closing the connection

#creating a function to load Visits and Sample tables
def load_visits_and_sample_tables(cur, file_path): #insert omics tables containing the sample ids
    visits_sql = '''INSERT INTO Visits (VisitID, SubjectID) VALUES (?, ?)''' #sql statement to to load Visits table
    sample_sql = '''INSERT INTO Sample (SampleID, SubjectID, VisitID) VALUES (?, ?, ?)''' #sql statement to to load Sample table

    try:
        with open(file_path) as file: 
            header = file.readline()  #skipping the header
            for line in file:
                row = line.strip().split('\t')  #splitting each row by tab (as it is a tsv file) into its contents
                sample_id = row[0]  #extracting unique SampleID from the first column
                subject_id, visit_id = sample_id.split('-')  #splitting SubjectID and VisitID from SampleID

                #checking if sample id already exists in the sample table to avoid duplicates
                cur.execute('''SELECT COUNT(*) FROM Sample WHERE SampleID = ?''', (sample_id,))
                count = cur.fetchone()[0]
                
                if count == 0:  #if sample id does not exist
                    #inserting data into Visits and Sample tables using above sql statements
                    cur.execute(visits_sql, (visit_id, subject_id))
                    cur.execute(sample_sql, (sample_id, subject_id, visit_id))
                else: #if sample id already exists, skipping it
                    continue

    except FileNotFoundError: #incase of an error in file name or path, informing the user
        print(f'{file_path} was not found.')
    except Exception as error: #incase of an error loading, informing the user
        print(f'Data could not be loaded into table due to : {error}')
        cur.close() #closing the connection

#creating a function to load Transcriptome, Proteome and Peak (Metabolome) tables
def load_omics_table(cur, file_path, table_name): #insert the omics file path and table name
    try:
        with open(file_path) as file:
            header = file.readline().strip().split('\t')  #splitting header to obtain entity ids
            for line in file:
                row = line.strip().split('\t')  #splitting each row by tab (as it is a tsv file) into its contents
                sample_id = row[0]  #extracting unique SampleID from the first column
                for i in range(1, len(row)):  #obtaining abundance values starting from the second column till the last
                    entity_id = header[i]  #extracting entity ids 
                    try:
                        abundance = float(row[i])  #extracting abundance value
                    except ValueError:
                        continue  #skipping if the value is not a number

                    #inserting data into omics table 
                    cur.execute(f'''INSERT INTO {table_name} ({table_name}ID, SampleID, Abundance) 
                    VALUES (?, ?, ?)''', (entity_id, sample_id, abundance))

    except FileNotFoundError: #incase of an error in file name or path, informing the user
        print(f'{file_path} was not found.')
    except Exception as error: #incase of an error loading, informing the user
        print(f'Data could not be loaded into {table_name} table due to : {error}')
        cur.close() #closing the connection
        
#creating a function to load MetabolomeAnnotation table
def load_metabolome_annotations(cur, file_path): #insert file containing metabolome annotations
    met_annot_sql = '''INSERT INTO MetabolomeAnnotation (PeakID, Metabolite_Name, KEGG_ID, HMDB_ID, Chemical_Class, Pathway) 
    VALUES (?, ?, ?, ?, ?, ?)''' #sql statement to load MetabolomeAnnotation table

    try:
        with open(file_path) as ma: 
            header = ma.readline()  #skipping the header line
            for line in ma:
                PeakID, Metabolite, KEGG, HMDB, Chemical_Class, Pathway = line.strip().split(',') #splitting each row by ',' into its contents
                
                #splitting metabolite, kegg and hmdb columns if they have multiple values separated  by '|' 
                metabolites = Metabolite.split('|')
                kegg_ids = KEGG.split('|')
                hmdb_ids = HMDB.split('|')

                #removing suffixes like (1), (2), etc from metabolite names
                single_metabolites = [re.sub(r'\(\d+\)$', '', metabolite) for metabolite in metabolites]

                #iterating over each metabolite, its corresponding kegg and hmdb ids
                for i in range(0, len(single_metabolites)): 
                    metabolite = single_metabolites[i]  #getting each metabolite
                    kegg = kegg_ids[i] if i < len(kegg_ids) else ''  #extracting corresponding kegg 
                    hmdb = hmdb_ids[i] if i < len(hmdb_ids) else ''  #extracting corresponding hmdb 
                    
                    # Insert data into the MetabolomeAnnotation table for this metabolite and its corresponding KEGG and HMDB IDs
                    cur.execute(met_annot_sql, (PeakID, metabolite, kegg, hmdb, Chemical_Class, Pathway))

    except FileNotFoundError:  # If the file is not found
        print(f'{file_path} was not found.')
    except Exception as error:  # Handle any other errors
        print(f'Data could not be loaded into MetabolomeAnnotation table due to : {error}')
        cur.close()  # Close the connection

#creating a function to load all the data
def load_db(database):
    try:
        conn = sqlite3.connect(database) #connecting to database
        cur = conn.cursor()
        load_subject_table(cur, subject_file) #loading Subject table
        load_visits_and_sample_tables(cur, metabolome_file)  #loading Visits and Sample tables using sample id from metabolome file
        load_visits_and_sample_tables(cur, transcriptome_file) #loading Visits and Sample tables using sample id from transcriptome file
        load_visits_and_sample_tables(cur, proteome_file) #loading Visits and Sample tables using sample id from proteome file
        load_omics_table(cur, transcriptome_file, 'Transcriptome')  #loading Transcriptome table
        load_omics_table(cur, proteome_file, 'Proteome') #loading Proteome table
        load_omics_table(cur, metabolome_file, 'Peak') #loading Peak (Metabolome) table
        load_metabolome_annotations(cur, metabolome_annot_file)  #loading MetabolomeAnnotation data table
        conn.commit() #committing the changes
        conn.close() #closing the connection
    except Exception as error: #incase of an error loading the entire database, informing the user
        print(f'Data could not be loaded into database due to :{error}')
        cur.close()
        conn.close() #closing the connection

# Function to run queries on the database
def query_db(database, query_no):
    #writing out the sql query statements by number
    q1 = '''SELECT SubjectID, Age FROM Subject 
            WHERE Age > 70 AND Age != "NA" '''
    q2 = '''SELECT SubjectID FROM Subject 
            WHERE Sex = "F" AND BMI BETWEEN 18.5 AND 24.9 
            ORDER BY SubjectID DESC'''
    q3 = '''SELECT VisitID FROM Visits 
            WHERE SubjectID = "ZNQOVZV"'''
    q4 = '''SELECT DISTINCT Subject.SubjectID FROM Subject
            JOIN Sample ON Subject.SubjectID = Sample.SubjectID
            JOIN Peak ON Sample.SampleID = Peak.SampleID
            WHERE Subject.Insulin_Status = "IR"
            AND Sample.SampleID IS NOT NULL AND Peak.SampleID IS NOT NULL'''
    q5 = '''SELECT DISTINCT KEGG_ID FROM MetabolomeAnnotation 
            WHERE PeakID IN ('nHILIC_121.0505_3.5', 'nHILIC_130.0872_6.3', 'nHILIC_133.0506_2.3', 'nHILIC_133.0506_4.4')'''
    q6 = '''SELECT MIN(Age), MAX(Age), AVG(Age) FROM Subject'''
    q7 = '''SELECT Pathway, COUNT(*) AS AnnotationCount 
            FROM MetabolomeAnnotation WHERE Pathway != ''
            GROUP BY Pathway HAVING COUNT(*) >= 10
            ORDER BY AnnotationCount DESC'''
    q8 = '''SELECT MAX(Abundance) FROM Transcriptome
            JOIN Sample ON Transcriptome.SampleID = Sample.SampleID
            WHERE Transcriptome.TranscriptomeID = 'A1BG' AND Sample.SubjectID = 'ZOZOW1T' '''
    q9 = '''SELECT Age, BMI FROM Subject WHERE Age IS NOT NULL AND BMI IS NOT NULL'''

    #making a list of query statements 
    queries = [q1, q2, q3, q4, q5, q6, q7, q8, q9]

    try:
        conn = sqlite3.connect(database) #connecting to database
        cur = conn.cursor()

        if 1 <= query_no <= 8: #for queries 1 to 8
            results = cur.execute(queries[query_no - 1]).fetchall()  #to obtain the query statement from the list, subtracting 1 as it is 0 indexed
            for result in results:
                print("\t".join(str(value) for value in result)) #to obtain results separated by tab instead of tuples

        elif query_no == 9: #for query 9
            results = cur.execute(q9).fetchall()  #same logic as the previous queries
            for result in results: #selecting age and bmi values that are not empty
                print("\t".join(str(value) for value in result))
            #using pandas, seaborn and matplotlib to generate scatterplot
            age_vs_bmi = pd.DataFrame(results, columns=["Age", "BMI"]) #creating dataframe of selected values
            sns.scatterplot(data= age_vs_bmi, x="Age", y="BMI") #creating the scatterplot
            plt.title("Age vs BMI") #titling the plot
            plt.savefig('age_bmi_scatterplot.png') #saving it in the working directory
            plt.clf() #clearing the figure

        else: #incase the query number is not between 1 to 9
            print('Query number must be between 1-9.')
 
    except Exception as error: #incase of error, informing the user
        print(f'Query could not be obtained due to : {error}.')
        cur.close()
        conn.close() #closing the connection
        
#creating command line of the format - python mainfile.py function dbfile.db
operation = sys.argv[1]
database = sys.argv[2]

if operation == "--createdb": #if user inputs --createdb, creating the database
        create_db(database)

elif operation == "--loaddb": #if user inputs --loaddb, loading the database
        load_db(database)

elif operation.startswith("--querydb="): #if user inputs --querydb=n, parsing to get the query number
    try:
        query_no = int(operation.split("=")[1]) #splitting the function user input into queryno= and the digit and obtaining the digit
        if 1 <= query_no <= 9: #if query is between 1-9, executing the query statements
            query_db(database, query_no)
        else:
            print('Insert a number between 1-9.')
    except ValueError:
        print('Query number is invalid. Insert a number between 1-9.')
        sys.exit(1) #exiting the program

else: #incase of error in user input, informing the user
    print('''Use the following format -
          python mainfile.py (--createdb | --loaddb | --querydb=n where n = 1 to 9) databasefile.db''')
    sys.exit(1) #exiting the program
