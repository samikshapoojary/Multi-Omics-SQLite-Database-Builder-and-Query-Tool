# Multi-Omics-SQLite-Database-Builder-and-Query-Tool
Python scripts to create, populate, and query an SQLite database integrating transcriptome, proteome, and metabolome datasets. Developed as coursework for MSc Bioinformatics 2024/25, University of Glasgow.


This project creates and manages a SQLite database integrating multi-omics data (transcriptome, proteome, metabolome) along with subject and visit metadata. It allows creation of the database schema, loading of data from various files, and querying for specific biological insights. It also supports generating a scatter plot of Age vs BMI.

## Files used

* **Subject.csv**: Subject demographic and clinical information.
* **HMP\_transcriptome\_abundance.tsv**: Transcriptome abundance data.
* **HMP\_proteome\_abundance.tsv**: Proteome abundance data.
* **HMP\_metabolome\_abundance.tsv**: Metabolome abundance data.
* **HMP\_metabolome\_annotation.csv**: Metabolome peak annotations.

## Database schema

The SQLite database contains the following tables:

* **Subject**: Demographic and clinical info per subject.
* **Visits**: Visit records linked to subjects.
* **Sample**: Sample metadata linked to subjects and visits.
* **Transcriptome**, **Proteome**, **Peak**: Abundance data tables linked to SampleIDs.
* **MetabolomeAnnotation**: Annotations of metabolome peaks with metabolite names and IDs.

## Requirements
Python 3.x

Modules: sqlite3, pandas, matplotlib, seaborn, re

## Usage

Run the script `mainfile.py` with the following command-line arguments:

```
python mainfile.py (--createdb | --loaddb | --querydb=n) databasefile.db
```

* `--createdb`
  Creates the database schema in the specified database file.

* `--loaddb`
  Loads data from the input files into the database tables.

* `--querydb=n`
  Runs query number `n` (1 to 9) against the database and prints results.
  Query 9 additionally generates an Age vs BMI scatter plot (`age_bmi_scatterplot.png`).

## Queries

| Query No | Description                                                                             |
| -------- | --------------------------------------------------------------------------------------- |
| 1        | List SubjectIDs and Age for subjects older than 70.                                     |
| 2        | List female SubjectIDs with BMI between 18.5 and 24.9 in descending order.              |
| 3        | List VisitIDs for subject "ZNQOVZV".                                                    |
| 4        | List distinct SubjectIDs with Insulin\_Status "IR" having samples and metabolome peaks. |
| 5        | List distinct KEGG\_IDs for specific metabolome peaks.                                  |
| 6        | Get minimum, maximum, and average Age of subjects.                                      |
| 7        | List Pathways with at least 10 metabolome annotations, ordered by count.                |
| 8        | Get maximum abundance of TranscriptomeID "A1BG" for subject "ZOZOW1T".                  |
| 9        | List Age and BMI pairs (not null), and generate Age vs BMI scatter plot.                |

## Notes

* Data loading functions handle missing or malformed values.
* Metabolite names and IDs with multiple entries separated by `|` are parsed and stored individually.
* The script requires the input files to be in the working directory or provide paths accordingly.
