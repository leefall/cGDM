#!/usr/bin/bash
import glob, slate, os

lPDFlist=glob.glob("*.pdf")

fp=open(lPDFlist[0],"r")
doc=slate.PDF(fp)

Linelist=doc[0].split("\n")

nTotalNumberofReads = Linelist[-31]
nTotalnumberofBases = Linelist[-29]
nTotalNumberofBasesAQ20=Linelist[-27]
nMeanCoverageDepth=Linelist[-25]
nCoveragewithTargetregion=  Linelist[-23]
nMeanReadLengthAQ20=Linelist[-21]
#nMeanReadLength(AQ30)=Linelist[-19]

lSecondlist=doc[1].split("\n")
TotalMappedFusionPanelReads =lSecondlist[-13]

sCurrentDir=os.getcwd()
sID=sCurrentDir.split("/")[-2]



fout=open("../../QC_Report_DNA_"+sID+".txt","w")

fout.write("Total Reads\tTotal Aligned Reads\t% Reads Aligned\tTotal Bases\tTotal Mapped Bases\tAverage on target depth\tStandard deviation on target depth \tOn Target Bases\n")

fout.write(str(nTotalNumberofReads)+"\t")
fout.write("\t")
fout.write("\t")
fout.write(str(float(nTotalnumberofBases)*1000000)+"\t")
fout.write(str(nMeanCoverageDepth)+"\t")
fout.write("\t")
#print nTotalNumberofBasesAQ20
#print nMeanCoverageDepth
fout.write(str(nCoveragewithTargetregion)+"\n")
fout.close()
fout=open("../../QC_Report_RNA_"+sID+".txt","w")
fout.write("Total Reads\tTotal Aligned Reads\t% Reads Aligned\tTotal Bases\tTotal Mapped Bases\tAverage on target depth\tStandard deviation on target depth \tOn Target Bases\n")
fout.write(TotalMappedFusionPanelReads)



