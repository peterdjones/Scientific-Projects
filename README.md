# portfolio
A collection of my R and python scripts used in my scientific research. 

Contents:

ApoptosisFlowAnalysis
  This script extracts complex data from all data files in a folder, normalises and transforms it, then ultimately presents several graphs summarising the data and exports a PDF of graphs as well as a .csv file containing the numerical results. I used this to drastically speed up my analysis of experiemental data, analysis the data manually would take 5-10 minutes for each sample, this script allows processing of a large number of samples with minimal setup. 

FlowRNAFlowAnaysisDESeq
  A script applying k-means clustering to some gene expression from artery cells exposed to different blood flow conditions. Using k-means clustering, 5 clusters were identified. This allowed me to examine which other genes our gene of interest clustered with and whether this told me anything meaningful. 

Heparin motif search
  A Python script to search a protein sequence for a Heparin binding motif, a variable short protein sequence, using a regular expression. Using regex is the only efficient way of carrying out this search. 
