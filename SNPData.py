import pandas as pd
import numpy as np
from pandas import ExcelWriter
from pandas import ExcelFile
from pip._vendor.distlib.compat import raw_input
import csv
import matplotlib.pyplot as plot


class SNPDataEnsemble:
    LDCorrelationCoefficient=0.70 #default value but can be changed
    print("Reading SNPs Data (Combined):")
    df_SNPsData = pd.read_excel(io='C:/Users/Onur/Desktop/SNPDataset.xlsx',
                                sheet_name='SNPDataset (Combined) (3)') # SNPs are given in SNPDataset.xlsx
    df_SNPsData['ChrNumber'] = df_SNPsData['Chromosome'].str[3:]
    df_SNPsData['Label'] = df_SNPsData['Label'].astype(str)
    # Define the sorter
    sorter = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
              'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
              'chrY']
    # Create the dictionary that defines the order for sorting
    sorterIndex = dict(zip(sorter, range(len(sorter))))
    print("Sorter Index Created by Chromosome Number...")
    df_SNPsData['Chr_Rank'] = df_SNPsData['Chromosome'].map(sorterIndex)
    df_SNPsData.sort_values(by=['Chr_Rank','Position'], ascending=[True, True], inplace=False)
    df_SNPsData['AutoNumberAfterSort'] = df_SNPsData.reset_index().index
    print( df_SNPsData.info(verbose=True))
    df_SNPsData = df_SNPsData[
        ["Score UP", "Score Middle", "Score BOTTOM", "Label", "Order", "SNP", "Chromosome","Position", "Chr_Rank", "DATABASE", "checkPRFRF"]]
    print("Reading SNPs Data (Combined) Completed Successfully...")
    print("--------")
    print("\n Reading LDs Data (Proxy Data:)")
    df_LDsData = pd.read_excel(io='C:/Users/Onur/Desktop/SNPDataset.xlsx', sheet_name='proxySearch.results')
    df_LDsData = df_LDsData[["QRSID", "RSID", "R2", "POS2", "DIST"]]
    print("Reading LDs Data (Proxy Data:) Completed Successfully...")
    print("--------")
    while True:
        initial_input = raw_input(
            " \n (LD correlation Coefficient Identification) Please enter a number between 1 and 0: ").strip()
        if "." in initial_input or initial_input.isnumeric():
            initial_input = float(initial_input)
            if 0 <= initial_input <= 1:
                LDCorrelationCoefficient = initial_input
                break
        else:
            print("Please enter the value correctly, it must be a real number between 0 and 1")
    print("LDs Data is Being Preprocessed Based On Selected Coefficient: " ,LDCorrelationCoefficient)
    a = df_LDsData.sort_values(by=['DIST'], ascending=[False])
    df_LDsData = a[(a['QRSID'] != a['RSID']) & (a['R2'] > LDCorrelationCoefficient) ]



    def writeToExcelFile(x, processType):
        if processType == 'SNP':
            with ExcelWriter('path_to_file_SNP.xlsx') as writer:
                x.to_excel(writer)
                print('SNPs Succesfully Exported To Excel')
        if processType == 'LD':
            filename=('path_to_file_LD_Reference_Table.xlsx')
            with ExcelWriter(filename) as writer:
                x.to_excel(writer)
                print('Based on User Input, LDs Succesfully Exported To Excel File')
        if processType == 'LD-Relationship':
            with ExcelWriter('path_to_file_SNP-LD_Relation.xlsx') as writer:
                x.to_excel(writer)
                print('SNPs are Succesfully Mapped by LDs and Exported To Excel')



    def findLinkageDisequilibrium (df_SNPsData, df_LDsData):
        LDDict={}
        for index, row in df_SNPsData.iterrows():
            RelatedLD=[]
            if not pd.isna(row ["checkPRFRF"]):
                print("Processing for # of ", index, "->", row["checkPRFRF"])
                for indexLD,rowLD in df_LDsData.iterrows():
                    if (row ["checkPRFRF"]==rowLD["QRSID"]):
                        RelatedLD.append(rowLD['RSID'])
                print (RelatedLD)
                LDDict[row ["checkPRFRF"]]=RelatedLD
        w = csv.writer(open("DeterminedLDOutput.csv", "w"))
        for key, val in LDDict.items():
            w.writerow([key, val])
        return LDDict
    
	

    def calculateCenteredScore(SNPsData):
        for index, row in SNPsData.iterrows():
            if not pd.isna(row["checkPRFRF"]):
                if index>0:
                    prev=index-1
                    while ((str(index) + ".LD") in SNPsData.at [prev,'Label']):
                        prev = prev - 1
                    startframe=prev
                else:
                    startframe=prev
                if index<len(SNPsData)-1:
                    nxt = index + 1
                    while ((str(index) + ".LD") in SNPsData.at[nxt, 'Label']):
                        nxt = nxt + 1
                    endframe = nxt
                else:
                    endframe=None
                file1 = open("MiddleScoreStartEnd.txt", "a")
                st = ("For #" + str(index) + " index value [Start: " + str(startframe) + ", End:" + str(endframe) + "]\n")
                file1.writelines (st)
                file1.close()
                if(isinstance(startframe, int) & isinstance(endframe, int)):
                    SSforDB = SNPsData[int(startframe): int(endframe) + 1]
                    numberofDatabase = SSforDB['DATABASE'].nunique()
                if index>0:
                    preneighbourhood= index-1
                else: preneighbourhood=None
                if index<len(SNPsData)-1:
                    postneighbourhood=index+1
                else:
                    postneighbourhood=None
                if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int)):
                    SSforLD = SNPsData[int(preneighbourhood): int(postneighbourhood) + 1]
                    if SSforLD.at[index, 'Label'] != 'nan':
                        if ((str(index) + ".LD") in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=1
                        elif ((str(index) + ".LD") not in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") not in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=3
                        elif ((str(index) + ".LD") in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") not in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=2
                        elif ((str(index) + ".LD") not in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=2
                    elif SSforLD.at[index, 'Label'] == 'nan':
                        if (SSforLD.at[index, 'Label'] ==SSforLD.at[preneighbourhood, 'Label']) and (SSforLD.at[index, 'Label'] ==SSforLD.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=3
                        elif (SSforLD.at[index, 'Label'] !=SSforLD.at[preneighbourhood, 'Label']):
                            numberofIndependentLD=3
                        elif (SSforLD.at[index, 'Label'] !=SSforLD.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=3
                    else: numberofIndependentLD="*"
                    if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int) & isinstance(startframe, int) & isinstance(endframe, int)):
                        SNPsData.at[index, 'Score Middle'] = (1*numberofDatabase) + (1/10*numberofIndependentLD) + (1/100)
                    else:
                        SNPsData.at[index, 'Score Middle'] =None



    def calculateFrameFirstScore (SNPsData):
        for index, row in SNPsData.iterrows():
            if not pd.isna(row["checkPRFRF"]):
                if index==0:
                    startframe=index
                elif index>0:
                    prev=index-1
                    counter=0
                    while ((str(index) + ".LD") in SNPsData.at[prev, 'Label']):
                        counter = counter + 1
                        prev=prev-1
                        if prev==0:
                            break
                    if counter > 0:
                        startframe = prev+1
                    else:
                        startframe = index
                if index<len(SNPsData)-2:
                    nxt = index + 1
                    counter=0
                    while ((str(index) + ".LD") in SNPsData.at[nxt, 'Label']):
                        counter=counter+1
                        nxt = nxt + 1
                        if nxt>len(SNPsData)-2:
                            break
                    if counter==0 and index<len(SNPsData)-2 :
                        endframe = index + 2
                    elif    counter==0 and index>=len(SNPsData)-2:
                        endframe=None
                    elif counter>0 and nxt>len(SNPsData)-2:
                        endframe=None
                    else:
                        endframe=nxt+1
                else:
                    endframe=None
                file2 = open("FrameFirstStartEnd.txt", "a")
                st = ("For #" + str(index) + " index value [Start: " + str(startframe) + ", End:" + str(endframe) + "]\n")
                file2.writelines (st)
                file2.close()
                if (isinstance(startframe, int) & isinstance(endframe, int)):
                    SSforDB = SNPsData[int(startframe): int(endframe) + 1]
                    numberofDatabase = SSforDB['DATABASE'].nunique()
                if index>0:
                    preneighbourhood= index-1
                else: preneighbourhood=None
                if index<len(SNPsData)-1:
                    postneighbourhood=index+1
                else:
                    postneighbourhood=None
                if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int)):
                    SSforLD = SNPsData[int(preneighbourhood): int(postneighbourhood) + 1]
                    if SSforLD.at[index, 'Label'] != 'nan':
                        if ((str(index) + ".LD") in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=1
                        elif ((str(index) + ".LD") not in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") not in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=3
                        elif ((str(index) + ".LD") in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") not in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=2
                        elif ((str(index) + ".LD") not in SNPsData.at[preneighbourhood, 'Label']) and ((str(index) + ".LD") in SNPsData.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=2
                    elif SSforLD.at[index, 'Label'] == 'nan':
                        if (SSforLD.at[index, 'Label'] ==SSforLD.at[preneighbourhood, 'Label']) and (SSforLD.at[index, 'Label'] ==SSforLD.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=3
                        elif (SSforLD.at[index, 'Label'] !=SSforLD.at[preneighbourhood, 'Label']):
                            numberofIndependentLD=3
                        elif (SSforLD.at[index, 'Label'] !=SSforLD.at[postneighbourhood, 'Label']):
                            numberofIndependentLD=3
                    else: numberofIndependentLD="*"
                    if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int) & isinstance(startframe, int) & isinstance(endframe, int)):
                        SNPsData.at[index, 'Score UP'] = (1*numberofDatabase) + (1/10*numberofIndependentLD) + (1/100)
                    else:
                        SNPsData.at[index, 'Score UP'] =None



    def calculateFrameLastScore(SNPsData):
        for index, row in SNPsData.iterrows():
            if not pd.isna(row["checkPRFRF"]):
                if index<2:
                    startframe=None
                elif index>=2:
                    prev=index-1
                    counter=0
                    while ((str(index) + ".LD") in SNPsData.at[prev, 'Label']):
                        counter = counter + 1
                        prev=prev-1
                        if prev==0:
                            break
                    if counter == 0:
                        startframe = index - 2
                    if counter > 0:
                        if prev + 1 > 1:
                            startframe = prev-1
                        else:
                            startframe = None
                else:
                    startframe=None
                if index == len(SNPsData) - 1:
                    endframe=index
                elif index>1 and index < len(SNPsData) - 1 :
                    nxt = index + 1
                    counter = 0
                    while ((str(index) + ".LD") in SNPsData.at[nxt, 'Label']):
                        counter = counter + 1
                        nxt = nxt + 1
                        if nxt == len(SNPsData) - 1:
                            break
                    if counter > 0:
                        endframe = nxt - 1
                    elif counter==0:
                        endframe = index
                else:
                    endframe=None
                file3 = open("FrameLastStartEnd.txt", "a")
                st = ("For #" + str(index) + " index value [Start: " + str(startframe) + ", End:" + str(endframe) + "]\n")
                file3.writelines (st)
                file3.close()
                if (isinstance(startframe, int) & isinstance(endframe, int)):
                    SSforDB = SNPsData[int(startframe): int(endframe) + 1]
                    numberofDatabase = SSforDB['DATABASE'].nunique()
                if index>0:
                    preneighbourhood= index-1
                else: preneighbourhood=None
                if index<len(SNPsData)-1:
                    postneighbourhood=index+1
                else:
                    postneighbourhood=None
                if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int)):
                    SSforLD = SNPsData[int(preneighbourhood): int(postneighbourhood) + 1]
                    if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int)):
                        SSforLD = SNPsData[int(preneighbourhood): int(postneighbourhood) + 1]
                        if SSforLD.at[index, 'Label'] != 'nan':
                            if ((str(index) + ".LD") in SNPsData.at[preneighbourhood, 'Label']) and (
                                    (str(index) + ".LD") in SNPsData.at[postneighbourhood, 'Label']):
                                numberofIndependentLD = 1
                            elif ((str(index) + ".LD") not in SNPsData.at[preneighbourhood, 'Label']) and (
                                    (str(index) + ".LD") not in SNPsData.at[postneighbourhood, 'Label']):
                                numberofIndependentLD = 3
                            elif ((str(index) + ".LD") in SNPsData.at[preneighbourhood, 'Label']) and (
                                    (str(index) + ".LD") not in SNPsData.at[postneighbourhood, 'Label']):
                                numberofIndependentLD = 2
                            elif ((str(index) + ".LD") not in SNPsData.at[preneighbourhood, 'Label']) and (
                                    (str(index) + ".LD") in SNPsData.at[postneighbourhood, 'Label']):
                                numberofIndependentLD = 2
                        elif SSforLD.at[index, 'Label'] == 'nan':
                            if (SSforLD.at[index, 'Label'] == SSforLD.at[preneighbourhood, 'Label']) and (
                                    SSforLD.at[index, 'Label'] == SSforLD.at[postneighbourhood, 'Label']):
                                numberofIndependentLD = 3
                            elif (SSforLD.at[index, 'Label'] != SSforLD.at[preneighbourhood, 'Label']):
                                numberofIndependentLD = 3
                            elif (SSforLD.at[index, 'Label'] != SSforLD.at[postneighbourhood, 'Label']):
                                numberofIndependentLD = 3
                        else:
                            numberofIndependentLD = "*"
                    if (isinstance(preneighbourhood, int) & isinstance(postneighbourhood, int) & isinstance(startframe, int) & isinstance(endframe, int)):
                        SNPsData.at[index, 'Score BOTTOM'] = (1*numberofDatabase) + (1/10*numberofIndependentLD) + (1/100)
                    else:
                        SNPsData.at[index, 'Score BOTTOM'] = None



    def maxScore(df_SNPsData):
        maxScores=df_SNPsData[['Score UP', 'Score Middle', 'Score BOTTOM', 'SNP']].max(axis=1)
        df_SNPsData['Max Scores']=maxScores
        return df_SNPsData



    def assignLinkDiseqLabels (SNPsData, LDDict):
        for index, row in SNPsData.iterrows():
            if not pd.isna(row["checkPRFRF"]):
                key=row["checkPRFRF"]
                LDValues=LDDict[row["checkPRFRF"]]
                print (key, ": ",  LDValues, "->", len (LDValues) )
                print('\n')
                if len(LDValues)>0:
                    for i in LDValues:
                        ind=int(SNPsData[SNPsData["SNP"] == i].index.values)
                        if SNPsData.at [ind,'Label']=='nan':
                            st = str(index) + ".LD"
                            SNPsData.at[ind, 'Label'] = st
                            SNPsData.at[index, 'Label'] = st
                        else:
                            st = (SNPsData.at [ind,'Label']) + ", " +  str(index) + ".LD"
                            SNPsData.at[ind, 'Label'] = st
                            SNPsData.at[index, 'Label'] = st
        return SNPsData



    if __name__ == "__main__":
        print('\n',"Scoring Algorithm Has Been Started.....")
        print("******************************")
        print("After filtering is implemented over LD Pairwise Data: ", len(df_LDsData),'\n')
        print("After combining Significant SNPs from Different Databases: ",len(df_SNPsData),'\n')
        print ("RF importances have been validated by a Second-run RF: ", df_SNPsData.checkPRFRF.notnull().sum(), '\n' )
        writeToExcelFile(df_SNPsData, "SNP")
        writeToExcelFile(df_LDsData, "LD")
        LDDictionary=findLinkageDisequilibrium (df_SNPsData,df_LDsData)
        print(LDDictionary)
        df_SNPsData=assignLinkDiseqLabels(df_SNPsData, LDDictionary)
        writeToExcelFile(df_SNPsData, "LD-Relationship")
        calculateCenteredScore(df_SNPsData)
        writeToExcelFile(df_SNPsData, "LD-Relationship")
        calculateFrameFirstScore(df_SNPsData)
        writeToExcelFile(df_SNPsData, "LD-Relationship")
        calculateFrameLastScore(df_SNPsData)
        writeToExcelFile(df_SNPsData, "LD-Relationship")
        finalScores=maxScore(df_SNPsData)
        writeToExcelFile(df_SNPsData, "LD-Relationship")
        temp = df_SNPsData.groupby(['Max Scores'])['Max Scores'].count().to_frame('Count').reset_index()
        formattedPoints = ["%.2f" % member for member in temp.iloc[:, 0]]
        points = list(formattedPoints)
        print (points)
        count = list(temp.iloc[:, 1])
        print (count)
        temp['Max Scores'] = 'p' + temp['Max Scores'].astype(str)
        bars=plot.bar(points, count, color ='maroon', width = 0.3, align='center')
        plot.xlabel("Available Scores")
        plot.ylabel("Observation")
        plot.title("Distribution of Points")
        for bar in bars:
            yval = bar.get_height()
            plot.text(bar.get_x(), yval + .005, yval)
        filename = 'ScoreDistribution - ' + 'LDCoef' + str(initial_input) + '.png'
        plot.savefig(filename, dpi=200,bbox_inches='tight')
        plot.show()
