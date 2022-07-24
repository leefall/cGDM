#!/usr/bin/env python
import os, glob, json, unicodedata, sys



def HGVS(sChr, nPosition, sRef, sAlt):

	
	sChr=sChr.replace("chr","NC_0000").replace("X","23").replace("Y","24").replace("MT","M")
	
	
	#print (sChr, nPosition, sRef, sAlt)
	
	if len(sRef)<len(sAlt):
		return sChr+".13:g."+str(nPosition)+"_"+str(int(nPosition)+len(sAlt))+"ins"+sAlt
	elif len(sRef)>len(sAlt):
		return sChr+".13:g."+str(nPosition)+"_"+str(int(nPosition)+len(sAlt))+"del"+sRef
	else:
		return sChr+".13:g."+str(nPosition)+sRef+">"+sAlt



def HGVS_CNV(sChr, nStartPosition, nEndPosition, nRatio):
	sChr=sChr.replace("chr","NC_0000").replace("X","23").replace("Y","24").replace("MT","M")
	
	
	#print (sChr, nPosition, sRef, sAlt)
	
	if float(nRatio)>1:
		return sChr+".13:g."+str(nStartPosition)+"_"+str(int(nEndPosition))+"dup"
	elif len(nRatio)<1:
		return sChr+".13:g."+str(nStartPosition)+"_"+str(int(nEndPosition))+"del"
		

def HGVS_Fusion(sLocus1, sLocus2):
	
	(sChr1,Position1)=sLocus1.split(":")
	(sChr2,Position2)=sLocus2.split(":")
	
	
	
	sChr1=sChr1.replace("chr","NC_0000").replace("X","23").replace("Y","24").replace("MT","M")
	
	sChr2=sChr2.replace("chr","NC_0000").replace("X","23").replace("Y","24").replace("MT","M")
	
	
	
	
	return sChr1+".13:g."+str(nPosition1)+"::"+sChr2+"13.g."+str(nPosition2)

def GetNormalized(sChr,nPosition,sRef):
	rawVcf=glob.glob("*Non-Filtered*.vcf")
	fVCF=open(rawVcf[0],"r")
	
	sKey=sChr+"_"+str(nPosition)+"_"+sRef
	
	for sLine in fVCF.xreadlines():
		if sLine[0]=="#":
			pass
		else:
			t=sLine.split("\t")
			sLineKey=t[0]+"_"+t[1]+"_"+t[3]
			if sLineKey==sKey:
				#print sChr,nPosition,sRef
				sFUNC=t[-3]
				sFUNC=sFUNC.split("=[")[1]
				sFUNC=sFUNC.replace("]","")
				sFUNC=sFUNC.replace("'",'"')
				#print sFUNC
				sjsonDict=json.loads(sFUNC)
#				print sLine
				try:
					(sNormalPos, sNormalRef, sNormalAlt)=(sjsonDict["normalizedPos"],sjsonDict['normalizedRef'],sjsonDict['normalizedAlt'])
					nNormalPos=unicodedata.normalize('NFKD', sNormalPos).encode('ascii','ignore')
					sNormalRef=unicodedata.normalize('NFKD', sNormalRef).encode('ascii','ignore')
					sNormalAlt=unicodedata.normalize('NFKD', sNormalAlt).encode('ascii','ignore')
					sGenotype=t[-1]
					
					sGenotype=sGenotype.split(":")[0]
					(sAlt1, sAlt2)=sGenotype.split("/")
				
					if sAlt1==sAlt2:
						sAllelicState="Homozygous"
					else:
						sAllelicState="Heterozygous"
				
					return (nNormalPos, sNormalRef, sNormalAlt,sAllelicState)
				except:
					sGenotype=t[-1]
					sGenotype=sGenotype.split(":")[0]
					(sAlt1, sAlt2)=sGenotype.split("/")
					if sAlt1==sAlt2:
						sAllelicState="Homozygous"
					else:
						sAllelicState="Heterozygous"
					
					return (t[1], t[3], t[4], sAllelicState)


def GetTranscriptID(sChr,nPosition,sRef):
	rawVcf=glob.glob("*Non-Filtered*.vcf")
	fVCF=open(rawVcf[0],"r")
	
	sKey=sChr+"_"+str(nPosition)+"_"+sRef
	
	for sLine in fVCF.xreadlines():
		if sLine[0]=="#":
			pass
		else:
			t=sLine.split("\t")
			sLineKey=t[0]+"_"+t[1]+"_"+t[3]
			if sLineKey==sKey:
				#print sChr,nPosition,sRef
				
				sFUNC=t[-3]
				print sFUNC
				sFUNC=sFUNC.split("=[")[1]
				#if 
				sFUNC=sFUNC.replace("]","")
				sFUNC=sFUNC.replace("'",'"')
				#print sFUNC
				sjsonDict=json.loads(sFUNC)
				(sTranscript, sGene,location)=(sjsonDict["transcript"],sjsonDict['gene'],sjsonDict['location'])
				sTranscript=unicodedata.normalize('NFKD', sTranscript).encode('ascii','ignore')
				sGene=unicodedata.normalize('NFKD', sGene).encode('ascii','ignore')
				location=unicodedata.normalize('NFKD', location).encode('ascii','ignore')
				
				if location=="exonic":
					nExonicPosition=sjsonDict["exon"]
				else:
					nExonicPosition=''
				
				return (sTranscript, sGene, location,nExonicPosition)



def SNV_INDEL_Parser(sLine,fSNV,sGenomebuild):
	sLine=sLine.strip()
	t=sLine.split("\t")
	(sLocus, sRef, nLength, sGenotype, sFilter, sGene, sNMID, sLocation, sExon, sProtein, sCoding, sClinvar, sCosmic, sDbsnp, sType,sFunction)=\
	(t[0], t[2], t[3], t[4], t[5], t[17], t[18], t[19], t[22], t[23], t[24], t[30], t[31], t[32], t[1],t[20])
	
	
	
	
	
	(sChr, nPosition)=sLocus.split(":")
	(sAlt1, sAlt2)=(sGenotype.split("/")[0], sGenotype.split("/")[1])
	
	
	
	#(sNMID, sGene, sLocation, sExon)=GetTranscriptID(sChr,nPosition,sRef)
	
	if sAlt1==sAlt2:
		sAllelicState="Homozygous"
	else:
		sAllelicState="Heterozygous"
	
	
	if sAlt1!=sRef:
		sAlt=sAlt1
	elif sAlt2!=sRef:
		sAlt=sAlt2
	else:
		print "Exception Error"
		print sLine
		sys.exit()
	
	if ((sType=="SNV") and (len(sRef)!=1) and (len(sAlt1)!=1 or len(sAlt2)!=1)):
		(nPosition, sRef, sAlt, sAllelicState)=GetNormalized(sChr,nPosition,sRef)
	
	
	if "|" in sNMID:
		sNMID=sNMID.split("|")[1]
	
	
	sHGVS=HGVS(sChr, nPosition, sRef, sAlt)
	
	if ((len(sRef)==1) and (len(sAlt)==1)):
		sVariantType="Substitution"
	elif (len(sRef)>len(sAlt)):
		sVariantType="Deletion"
	elif (len(sRef)<len(sAlt)):
		sVariantType="Insertion"
	else:
		print "Error"
		print sLine
		sys.exit()
	
	fSNV.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\n".format(\
	sHGVS, sNMID+":"+sCoding,sProtein.replace("|",""),"HGVS version 15.11",sGenomebuild,"Somatic",sGene, '','',sChr.replace("chr",""),nPosition, sRef, sAlt, "","","",sExon.replace("|",""), sFunction.replace("|",""), sVariantType,\
	sAllelicState, sDbsnp, sClinvar, sCosmic,"","","",""))
	
	


def CNV_Parser(sLine, fCNV):
	
	(sLocus, sRef, nLength, sGenotype, sFilter, sGene, sNMID, sLocation, sExon, sProtein, sCoding, sClinvar, sCosmic, sDbsnp, sType,sFunction,sIscn)=\
	(t[0], t[2], t[3], t[4], t[5], t[17], t[18], t[19], t[22], t[23], t[24], t[30], t[31], t[32], t[1],t[20],t[11])
	
	sChr=sLocus.split(":")[0]
	
	(sCytoband, nRatio)=(sIscn.split(")x")[0], sIscn.split("x")[1])
	
	(sCytoband, nPosition)=sCytoband.split("(")
	
	(nStartPosition, nEndPosition)=nPosition.split("-")
	
	sHGVS=HGVS_CNV(sChr, nStartPosition, nEndPosition, nRatio)
	
	if float(nRatio)>2.0:
		sCopynumberType="GN"
	else:
		sCopynumberType="LS"
	
	
	fCNV.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\n".format(\
	sHGVS,"","","HGVS version 15.11",sChr.replace("chr",""), nStartPosition, nEndPosition, "","",nRatio,sCopynumberType,sGene,"","",sGenomebuild,"Somatic","","","","",sCosmic,"","","","","","","","",""))
	
	



def Fusion_Parser(sLine, fFusion):
	t=sLine.split("\t")
	
	
	(sLocus, sGene, sExon)=(t[0], t[17], t[22])
	
	sFusion_Locus=sLocus.split("_")[0]
	
	(sLocus1, sLocus2)=sFusion_Locus.split("-")[0], sFusion_Locus.split("-")[1]
	(sGene1, sGene2)=sGene.split("|")
	(nExon1, nExon2)=sExon.split("|")
	
	sHGVS=HGVS_Fusion(sLocus1, sLocus2)
	
	fFusion.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\n".format(\
	sHGVS, "", "", "HGVS version 15.11", sGene1, "", "", sGene2, "","","","Upstream","Downstream", sGenomebuild,"","","","","","","","","","","","","","","",""))
	
	
	
	
	

	
	
lFile=glob.glob("*RNA-full.tsv")

fp=open(lFile[0],"r")

sReference=fp.readline()
sReference=sReference.strip()
sGenomebuild=sReference.split("=")[1]

fp.readline()
fp.readline()

fSNV=open("SNVINDEL_"+lFile[0],"w")
fCNV=open("CNV_"+lFile[0],"w")
fFusion=open("Translocation_"+lFile[0],"w")

fCNV.write("HGVS genomic change\tHGVS coding change\tHGVS protein change\tHGVS version\tChromosome\tGenomic Start Position\tGenomic Stop Position\tBreakpoint exon, From\tBreakpoint exon, To\tCopy number\tCopy number type\tHGNC gene symbol(s)\tEntrez gene ID(s)\tEnsembl gene ID(s)\tGenome Build\tGenomic source\tCopy number ratio type\tCopy number ratio\tFunctional Domain\tdbVar Identification Number\tCOSMIC Identification Number\tCoding Size\tBiomarker name\tBiomarker state\tWildtype biomarker YN\n")
fSNV.write("HGVS genomic change\tHGVS coding change\tHGVS protein change\tHGVS version\tGenome build\tGenomic source\tHGNC gene symbol\tEntrez gene ID\tEnsembl gene ID\tChromosome\tPosition\tReference allele\tAlternative allele\tCytogenetic location\tStrand orientation\tCodon\tExon\tMolecular Effects\tVariant type\tGenotype\tdbSNP Identification Number\tclinVar Variation Identification Number\tCOSMIC Identification Number\tBiomarker name\tBiomarker state\tWildtype biomarker YN\tFunctional Domain\n")






for sLine in fp.xreadlines():
	t=sLine.split("\t")
	if t[5]=="PASS":
		if (t[1]=="SNV") or (t[1]=="INDEL"):
			SNV_INDEL_Parser(sLine,fSNV,sGenomebuild)
		elif (t[1]=="CNV"):
			CNV_Parser(sLine, fCNV)
		elif (t[1]=="FUSION"):
			Fusion_Parser(sLine, fFusion)
		else:
			pass
		
	
	



	
	
	
	
	









