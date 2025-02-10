# EnsembleScoringAlgorithm-EnSCAN

This Python script processes SNP (Single Nucleotide Polymorphism) and LD (Linkage Disequilibrium) datasets to compute various scores and relationships between SNPs. It includes data preprocessing, filtering, and scoring mechanisms that assist in identifying significant SNPs and their relationships based on LD correlation coefficients. The final output is an Excel file containing SNV data along with their calculated scores, and a bar chart illustrating the distribution of the calculated scores. 


## Features 
-	**Data Preprocessing:** The script loads and preprocesses SNP and LD data from Excel files.
-	**LD Identification:** It computes the linkage disequilibrium (LD) between SNPs based on a user-defined correlation coefficient.
-	**Scoring Mechanism:** SNPs are scored based on various criteria, including their proximity to other SNPs and their relationship with the dataset.
-	**Output:** The processed data and scores are exported to Excel files, and a distribution plot of SNP scores is generated.

## Prerequisites
Python 3.6+
To run the script, you'll need the following Python libraries:
-	pandas
-	numpy
-	matplotlib
-	openpyxl (for Excel operations)
-	pip

You can install the required dependencies via pip by running the following:
```
bash
pip install pandas numpy matplotlib openpyxl
```

## Installation
1.	Clone the repository to your local machine:
2.	git clone https://github.com/erdogan-onur/EnsembleScoringAlgorithm-EnSCAN.git
3.	Navigate into the project directory:
4.	cd SNPData
5.	Install the necessary Python libraries (as mentioned in Prerequisites).

## Key Functions
- **findLinkageDisequilibrium:** Maps SNPs to their respective LDs based on the correlation coefficient.
- **calculateCenteredScore:** Computes a middle score for SNPs.
- **calculateFrameFirstScore:** Computes a score for SNPs based on their initial database representation.
- **calculateFrameLastScore:** Computes a score for SNPs based on their terminal database representation.
- **maxScore:** Calculates the maximum score for each SNP.

 ## File Structure
 ```
.
├── SNPDataEnsemble.py  # Main Python script for analysis
├── SNPDataset.xlsx     # Input SNP dataset (Excel file)
└── DeterminedLDOutput.csv # Output CSV for identified LD relationships

```

## How To Use
1. Place your SNP and LD datasets in an Excel file, e.g., SNPDataset.xlsx, with the following sheets:
      1. SNP Dataset: Contains SNP information with columns like Chromosome, Position, DATABASE, checkPRFRF, etc.
      2. Proxy Search Results: Contains LD data with columns like QRSID, RSID, R2, etc.
2. Update the file paths in the script:
     1. Replace C:/Users/<yourname>/Desktop/SNPDataset.xlsx with the path to your dataset.
3. Run the script:
     1. python script_name.py
4. Follow the on-screen prompt to input an LD correlation coefficient:
     1. A number between 0 and 1.
5. Processed files will be exported as:
     1. SNP data: path_to_file_SNP.xlsx
     2. Filtered LD data: path_to_file_LD_Reference_Table.xlsx
     3. SNP-LD relationships: path_to_file_SNP-LD_Relation.xlsx
     4. LD outputs: DeterminedLDOutput.csv
     5. Logs for scoring: MiddleScoreStartEnd.txt, FrameFirstStartEnd.txt, FrameLastStartEnd.txt

## Workflow
1.	**Data Loading:**
The script begins by loading the SNP and LD data from two sheets of an Excel file (SNPDataset.xlsx). SNPs are sorted based on their chromosome and position to create a well-organized dataset for further analysis.
2.	**LD Correlation Coefficient:**
The user is prompted to provide an LD correlation coefficient between 0 and 1. This value is used to filter out irrelevant LD pairs based on their R² value, which measures the correlation between SNPs.
3.	**Identifying Linkage Disequilibrium (LD):**
The findLinkageDisequilibrium function identifies SNPs that are in LD with each other, storing the relationships in a dictionary.
4. **Scoring:**
The script calculates three different scores for each SNP:
    1. **Score UP:** Based on its position relative to other SNPs.
    2. **Score Middle:** Incorporates LD information to assign a middle score.
    3. **Score BOTTOM:** A final score considering both neighboring SNPs and the SNP's LD relationships.

    These scores help assess the significance of each SNP within the dataset.
5. **Exporting Results:**
The results, including SNP data and LD relationships, are exported to separate Excel files for further analysis:
    1. SNPs and their scores
    2. LD relationships
    3. A distribution plot of the SNP scores
6. **Visualization:**
A bar chart of the distribution of SNP scores is generated and saved as a PNG file. The plot shows how the scores are distributed and highlights the most significant SNPs.

## Example Usage
To run the script, execute the following command in your terminal:
```
bash
python SNPDataEnsemble.py
```
You will be prompted to enter an LD correlation coefficient between 0 and 1. The script will proceed with the analysis and output the results to Excel files and a plot.

- LD correlation Coefficient Identification) Please enter a number between 1 and 0: 0.7
- Scoring Algorithm Has Been Started.....
- After filtering is implemented over LD Pairwise Data: 500
- After combining Significant SNPs from Different Databases: 1000
- RF importances have been validated by a Second-run RF: 800
- SNPs Successfully Exported to Excel
- Based on User Input, LDs Successfully Exported to Excel File
- SNPs are Successfully Mapped by LDs and Exported to Excel

## Output Files
The script generates multiple outputs such as SNP data enriched with scores and LD relationships, filtered LD datasets based on user input and logs for debugging and detailed scoring processes.
- SNP-LD Relationship File (path_to_file_SNP-LD_Relation.xlsx): Contains SNPs and their associated LD relationships.
- SNP File (path_to_file_SNP.xlsx): Contains SNP data and their calculated scores.
- LD File (path_to_file_LD_Reference_Table.xlsx): Contains LD pair data.
- Score Distribution Plot: A PNG file showing the distribution of SNP scores.

## Contributing
Contributions are welcome! If you have suggestions for improvements or additional features, feel free to fork the repository and create a pull request. Please ensure that all tests pass before submitting a pull request.

## License
You can adjust file paths, names, or additional details to fit your specific needs. Let me know if you'd like further customizations!

## Instructions on How to Use This “README”:
1. Replace the GitHub repository URL (“https://github.com/erdogan-onur/EnsembleScoringAlgorithm-EnSCAN.git”) with your own repository link if you have already created one on GitHub.
2. In the “How To Use” section, ensure that the paths (like “C:/Users/<yourname>/Desktop/SNPDataset.xlsx”) match the actual locations of your files.
3. If you need to modify the script (e.g., change file paths or other specifics), make sure to update the instructions accordingly.

---

This README follows GitHub's markdown guidelines and provides a clear, organized summary of the repository, installation steps, and usage instructions. Let me know if you need further assistance via onur.erdogan@metu.edu.tr and/or yesim@metu.edu.tr





   
