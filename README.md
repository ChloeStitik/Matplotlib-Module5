# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np
from scipy.stats import linregress

# Study data files
mouseMetadataPath = "data/Mouse_metadata.csv"
studyResultsPath = "data/Study_results.csv"

# Read the mouse data and the study results
mouseMetadata = pd.read_csv(mouse_metadata_path)
studyResults = pd.read_csv(study_results_path)

# Combine the above data into a single dataset
combinedDF = pd.merge(mouseMetadata, studyResults, how='outer', on="Mouse ID")
# Display the data table for preview
combinedDF.head()
Mouse ID	Drug Regimen	Sex	Age_months	Weight (g)	Timepoint	Tumor Volume (mm3)	Metastatic Sites
0	k403	Ramicane	Male	21	16	0	45.000000	0
1	k403	Ramicane	Male	21	16	5	38.825898	0
2	k403	Ramicane	Male	21	16	10	35.014271	1
3	k403	Ramicane	Male	21	16	15	34.223992	1
4	k403	Ramicane	Male	21	16	20	32.997729	1
# Checking the number of mice.
numMice = combinedDF["Mouse ID"].nunique()
numMice
249
# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint.
duplicateMiceID = combinedDF.loc[combinedDF.duplicated(subset=['Mouse ID', 'Timepoint']),'Mouse ID'].unique()
duplicateMiceID
array(['g989'], dtype=object)
# Optional: Get all the data for the duplicate mouse ID. 
duplicateMiceDF = combinedDF.loc[combinedDF["Mouse ID"] == "g989", :]
duplicateMiceDF
Mouse ID	Drug Regimen	Sex	Age_months	Weight (g)	Timepoint	Tumor Volume (mm3)	Metastatic Sites
908	g989	Propriva	Female	21	26	0	45.000000	0
909	g989	Propriva	Female	21	26	0	45.000000	0
910	g989	Propriva	Female	21	26	5	48.786801	0
911	g989	Propriva	Female	21	26	5	47.570392	0
912	g989	Propriva	Female	21	26	10	51.745156	0
913	g989	Propriva	Female	21	26	10	49.880528	0
914	g989	Propriva	Female	21	26	15	51.325852	1
915	g989	Propriva	Female	21	26	15	53.442020	0
916	g989	Propriva	Female	21	26	20	55.326122	1
917	g989	Propriva	Female	21	26	20	54.657650	1
918	g989	Propriva	Female	21	26	25	56.045564	1
919	g989	Propriva	Female	21	26	30	59.082294	1
920	g989	Propriva	Female	21	26	35	62.570880	2
# Create a clean DataFrame by dropping the duplicate mouse by its ID.
cleanDF = combinedDF[combinedDF['Mouse ID'].isin(duplicateMiceDF)==False]
cleanDF.head()
Mouse ID	Drug Regimen	Sex	Age_months	Weight (g)	Timepoint	Tumor Volume (mm3)	Metastatic Sites
0	k403	Ramicane	Male	21	16	0	45.000000	0
1	k403	Ramicane	Male	21	16	5	38.825898	0
2	k403	Ramicane	Male	21	16	10	35.014271	1
3	k403	Ramicane	Male	21	16	15	34.223992	1
4	k403	Ramicane	Male	21	16	20	32.997729	1
# Checking the number of mice in the clean DataFrame.
cleanMiceTot = cleanDF["Mouse ID"].nunique()

cleanMiceTot
249
Summary Statistics
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
mean = cleanDF['Tumor Volume (mm3)'].groupby(cleanDF['Drug Regimen']).mean()
median = cleanDF['Tumor Volume (mm3)'].groupby(cleanDF['Drug Regimen']).median()
var = cleanDF['Tumor Volume (mm3)'].groupby(cleanDF['Drug Regimen']).var()
std = cleanDF['Tumor Volume (mm3)'].groupby(cleanDF['Drug Regimen']).std()
sem = cleanDF['Tumor Volume (mm3)'].groupby(cleanDF['Drug Regimen']).sem()

summaryTable = pd.DataFrame({"Mean Tumor Volume":mean, 
                            "Median Tumor Volume":median, 
                           "Tumor Volume Variance":var, 
                           "Tumor Volume Std. Dev.":std, 
                           "Tumor Volume Std. Err.":sem})
# Display the Summary statistics table grouped by 'Drug Regimen' column
summaryTable
Mean Tumor Volume	Median Tumor Volume	Tumor Volume Variance	Tumor Volume Std. Dev.	Tumor Volume Std. Err.
Drug Regimen					
Capomulin	40.675741	41.557809	24.947764	4.994774	0.329346
Ceftamin	52.591172	51.776157	39.290177	6.268188	0.469821
Infubinol	52.884795	51.820584	43.128684	6.567243	0.492236
Ketapril	55.235638	53.698743	68.553577	8.279709	0.603860
Naftisol	54.331565	52.509285	66.173479	8.134708	0.596466
Placebo	54.033581	52.288934	61.168083	7.821003	0.581331
Propriva	52.322552	50.854632	42.351070	6.507770	0.512884
Ramicane	40.216745	40.673236	23.486704	4.846308	0.320955
Stelasyn	54.233149	52.431737	59.450562	7.710419	0.573111
Zoniferol	53.236507	51.818479	48.533355	6.966589	0.516398
# Using the aggregation method, produce the same summary statistics in a single line
summaryAg =  cleanDF.groupby(['Drug Regimen'])[['Tumor Volume (mm3)']].agg(['mean', 'median', 'var', 'std', 'sem'])
summaryAg
Tumor Volume (mm3)
mean	median	var	std	sem
Drug Regimen					
Capomulin	40.675741	41.557809	24.947764	4.994774	0.329346
Ceftamin	52.591172	51.776157	39.290177	6.268188	0.469821
Infubinol	52.884795	51.820584	43.128684	6.567243	0.492236
Ketapril	55.235638	53.698743	68.553577	8.279709	0.603860
Naftisol	54.331565	52.509285	66.173479	8.134708	0.596466
Placebo	54.033581	52.288934	61.168083	7.821003	0.581331
Propriva	52.322552	50.854632	42.351070	6.507770	0.512884
Ramicane	40.216745	40.673236	23.486704	4.846308	0.320955
Stelasyn	54.233149	52.431737	59.450562	7.710419	0.573111
Zoniferol	53.236507	51.818479	48.533355	6.966589	0.516398
Bar and Pie Charts
# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.
timePoint = cleanDF["Timepoint"].value_counts()
timePoint

pandasPlot= timePoint.plot.bar(color='b')  
plt.xlabel("Drug Regimen")
plt.ylabel("Timepoint")
plt.title("Number of Timepoints per Drug Regimen")
Text(0.5, 1.0, 'Number of Timepoints per Drug Regimen')

# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
xAxis = timePoint.index.values
yAxis = timePoint.values

plt.bar(xAxis, yAxis, color='b', align='center')

plt.title("Number of Timepoints per Drug Regimen")
plt.xlabel("Drug Regimen")
plt.ylabel("Timepoint")
plt.xticks(rotation="vertical")

plt.show()

# Generate a pie plot showing the distribution of female versus male mice using Pandas
genderDis = cleanDF["Sex"].value_counts()
plt.title("Female vs Male Mice Pop with Pandas")
genderDis.plot.pie(autopct= "%1.1f%%")
plt.show()

# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ['Female', 'Male']
sizes = [49.7999197, 50.200803]
plot = genderDis.plot.pie(y='Total Count', autopct="%1.1f%%")
plt.title('Female vs Male Mice Pop with PyPlot')
plt.ylabel('Sex')
plt.show()

Quartiles, Outliers and Boxplots
# Calculate the final tumor volume of each mouse across four of the treatment regimens: 

CapomulinDF = cleanDF.loc[cleanDF["Drug Regimen"] == "Capomulin",:]
RamicaneDF = cleanDF.loc[cleanDF["Drug Regimen"] == "Ramicane", :]
InfubinolDF = cleanDF.loc[cleanDF["Drug Regimen"] == "Infubinol", :]
CeftaminDF = cleanDF.loc[cleanDF["Drug Regimen"] == "Ceftamin", :]

# Start by getting the last (greatest) timepoint for each mouse

lastTime = drugs.groupby(["Drug Regimen", "Mouse ID"]).agg(tumorSize=("Tumor Volume (mm3)", lambda x: x.iloc[-1]))
lastTime = lastTime.stack(level=0).unstack(level=0)

for drug in drugList:
    print(drug)

# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
drugList = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
drugs = combinedDF[combinedDF["Drug Regimen"].isin(drugList)]
drugs.head()
Capomulin
Ramicane
Infubinol
Ceftamin
Mouse ID	Drug Regimen	Sex	Age_months	Weight (g)	Timepoint	Tumor Volume (mm3)	Metastatic Sites
0	k403	Ramicane	Male	21	16	0	45.000000	0
1	k403	Ramicane	Male	21	16	5	38.825898	0
2	k403	Ramicane	Male	21	16	10	35.014271	1
3	k403	Ramicane	Male	21	16	15	34.223992	1
4	k403	Ramicane	Male	21	16	20	32.997729	1
# Put treatments into a list for for loop (and later for plot labels)
treatments = 0

# Create empty list to fill with tumor vol data (for plotting)
# Calculate the IQR and quantitatively determine if there are any potential outliers. 
# Locate the rows which contain mice on each drug and get the tumor volumes
# add subset 
# Determine outliers using upper and lower bounds
for drug in drugList:
    quartiles = lastTime[drug].quantile([.25,.5,.75]).round(2)
    lowerq = quartiles[0.25].round(2)
    upperq = quartiles[0.75].round(2)
    iqr = round(upperq-lowerq,2)
    lowerBound = round(lowerq - (1.5*iqr),2)
    upperBound = round(upperq + (1.5*iqr),2)


    if treatments == 0:
        print(f"-------------------------------------------------")
    print(f"The lower quartile of {drug} treatments is: {lowerq}")
    print(f"The upper quartile of {drug} treatments is: {upperq}")
    print(f"The interquartile range of {drug} treatments is: {iqr}")
    print(f"Values below {lowerBound} could be {drug} outliers.")
    print(f"Values above {upperBound} could be {drug} outliers.")
    print(f"------------------------------------------------------")
    treatments+=1
    
-------------------------------------------------
The lower quartile of Capomulin treatments is: 32.38
The upper quartile of Capomulin treatments is: 40.16
The interquartile range of Capomulin treatments is: 7.78
Values below 20.71 could be Capomulin outliers.
Values above 51.83 could be Capomulin outliers.
------------------------------------------------------
The lower quartile of Ramicane treatments is: 31.56
The upper quartile of Ramicane treatments is: 40.66
The interquartile range of Ramicane treatments is: 9.1
Values below 17.91 could be Ramicane outliers.
Values above 54.31 could be Ramicane outliers.
------------------------------------------------------
The lower quartile of Infubinol treatments is: 54.05
The upper quartile of Infubinol treatments is: 65.53
The interquartile range of Infubinol treatments is: 11.48
Values below 36.83 could be Infubinol outliers.
Values above 82.75 could be Infubinol outliers.
------------------------------------------------------
The lower quartile of Ceftamin treatments is: 48.72
The upper quartile of Ceftamin treatments is: 64.3
The interquartile range of Ceftamin treatments is: 15.58
Values below 25.35 could be Ceftamin outliers.
Values above 87.67 could be Ceftamin outliers.
------------------------------------------------------
# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
boxplotList = []
for drug in drugList:
    boxplotList.append(list(lastTime[drug].dropna()))
    
fig = plt.figure()
plt.xlabel("Regimen")
plt.xticks([1,2,3,4], drugList,  rotation=45)
plt.ylabel("Tumor Volume")
plt.title("Tumor Volume by Drug Regimen")
plt.boxplot(boxplotList)
plt.show()   

Line and Scatter Plots
# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
forlineDF = CapomulinDF.loc[CapomulinDF["Mouse ID"] == "m957",:]
forlineDF.head()
xAxis = forlineDF["Timepoint"]
tumorSize = forlineDF["Tumor Volume (mm3)"]

fig1, ax1 = plt.subplots()
plt.title('Capomulin treatmeant of mouse 957')
plt.plot(xAxis, tumorSize,linewidth=1, markersize=10,marker="o",color="purple")
plt.xlabel('Timepoint (Days)')
plt.ylabel('Tumor Volume (mm3)')
Text(0, 0.5, 'Tumor Volume (mm3)')

# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
fig1, ax1 = plt.subplots()
avgVol =CapomulinDF.groupby(['Mouse ID']).mean()

markerSize=10
plt.scatter(avgVol['Weight (g)'],avgVol['Tumor Volume (mm3)'], color="green")
plt.title('Mouse Weight Versus Average Tumor Volume')
plt.xlabel('Weight (g)',fontsize =10)
plt.ylabel('Average Tumor Volume (mm3)')
Text(0, 0.5, 'Average Tumor Volume (mm3)')

Correlation and Regression
# Calculate the correlation coefficient and linear regression model for mouse weight and average tumor volume for the Capomulin regimen
correlation = st.pearsonr(avgVol['Weight (g)'],avgVol['Tumor Volume (mm3)'])
print(f"The correlation between the weight of the mice and the average tumor volume is {round(correlation[0],2)}.")
The correlation between the weight of the mice and the average tumor volume is 0.84.
#add the linear regression line to the scatter plot from above
(slope, intercept,rvalue, pvalue, stderr)= linregress(avgVol["Weight (g)"],avgVol["Tumor Volume (mm3)"])
regressValues=avgVol["Weight (g)"]* slope + intercept
lineEq= f"y = {round(slope, 2)} x + {round(intercept, 2)}"


markersize=10
plt.scatter(avgVol['Weight (g)'],avgVol['Tumor Volume (mm3)'], color="green")
plt.plot(avgVol["Weight (g)"], regressValues, color='orange')
plt.annotate(lineEq,(20,36), fontsize=8)
plt.xlabel('Weight (g)',fontsize =10)
plt.ylabel('Tumor Volume (mm3)')
plt.title('Mouse Weight Vs. Average Tumor Volume for Capomulin')
print(f"The r-squared is: {round(rvalue**2,3)}")
plt.show()
The r-squared is: 0.709

 
