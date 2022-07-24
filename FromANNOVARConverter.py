#!/usr/bin/env python
import os, subprocess, glob, time, tabix
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager


logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


def Clinvar(sFile):
	#print sFile
	#fp=open(sFile,"r")
	dClinvar=dict()
	url="/storage/home/leefall2/clara/Personal/"
	
	url+="clinvar_20171231.vcf.gz"

	tb=tabix.open(url)
	fp=open(sFile,"r")
	
	for sLine in fp.xreadlines():
		t=sLine.split("\t")
		sChr=t[0]
		sStartPosition=t[1]
		sEndPosition=t[2]
		sReference_Position=t[3]
		sAlternativeOne=t[4]
			#sAlternativeTwo=t[8]
	
		sChr=sChr.replace("chr","")
	
		tb=tabix.open(url)
	
		records=tb.query(sChr,int(sStartPosition)-1,int(sEndPosition)+1)
	
	
		for record in records:
	#		sClinvarID=record[2]
			sRef=record[3]
			sAlt=record[4]
			sInfo=record[7]
			lInfo=sInfo.split(";")
			sVariants=sAlt+","+sRef
			lVariants=sVariants.split(",")
			sClinvarID=record[2]
			if ((sRef==sReference_Position) and ((sAlt==sAlternativeOne)) ):
				
				sKey="|".join([sChr, sStartPosition, sReference_Position, sAlternativeOne])
				dClinvar[sKey]=sClinvarID
		
	
	return dClinvar






def Intron_Dict(sFile):
	sIntronDict=dict()
	fp=open(sFile,'r')
	dGeneDict=dict()
	for sLine in fp.xreadlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#print sLine
		(sType, sGene, sChr, sStartPosition, sEndPosition, sRef, sAlt, sInfo)=\
		(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[-1])
		#print lFrequency
		

		
		if sType=="exonic":
			pass
		else:
			sKey="|".join([sChr,sStartPosition,sRef,sAlt])
			if sType=="intergenic":
				sGene=''
			sIntronDict[sKey]=[sType,sGene,sInfo]
	
			
	return sIntronDict


def ExonicDict(sFile):
	dexonicANNOVAR=dict()
	fGeneANNOVAR=open(sFile,"r")

	for sLine in fGeneANNOVAR.xreadlines():
	#	print sLine
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition, sRef, sAlt)=(t[3], t[4], t[6], t[7])
		sKey="|".join([sChr,nPosition,sRef,sAlt])
		(sHGVS)=t[2]
		lHGVS=sHGVS.split(",")
		Allele=t[-1]
		
		
		lNMIDs=[]
		nExons=[]
		sCodons=[]
		
		for aHGVS in lHGVS:
			#print aHGVS
			ag=aHGVS.split(":")
		#	print ag
			try:
				(sNMID, sExonNumber, sCodonChange, sAAChange)=(ag[1], ag[2], ag[3], ag[4])
				(sGene)=ag[0]
				#lNMIDs.append(sNMID)
				nExonNumber=sExonNumber.replace("exon","")
				lNMIDs.append(sNMID+":"+sCodonChange)
				nExons.append(nExonNumber)
				sCodons.append(sAAChange)
			except IndexError:
				pass
			
		
		dexonicANNOVAR[sKey]=[lNMIDs,nExons,sCodons,Allele,sGene]
		
	return dexonicANNOVAR
	
	
def CosmicCoding(sFile):
	dCosmicCoding=dict()
	fGeneANNOVAR=open(sFile,"r")

	for sLine in fGeneANNOVAR.xreadlines():
	#	print sLine
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition, sRef, sAlt)=(t[2], t[3], t[5], t[6])
		(sCosmicID)=t[1]
		sCosmicID=sCosmicID.split(";")[0]
		sCosmicID=sCosmicID.replace("ID=","")
		sKey="|".join([sChr,nPosition,sRef,sAlt])
		dCosmicCoding[sKey]=sCosmicID
		
		
	return dCosmicCoding
	






def ENSG_Exonic(sFile):
	dCosmicCoding=dict()
	fGeneANNOVAR=open(sFile,"r")

	for sLine in fGeneANNOVAR.xreadlines():
	#	print sLine
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition, sRef, sAlt)=(t[3], t[4], t[6], t[7])
		(sCosmicID)=t[2]
		sCosmicID=sCosmicID.split(":")[0]
		
		sKey="|".join([sChr,nPosition,sRef,sAlt])
		dCosmicCoding[sKey]=sCosmicID
		
		
	return dCosmicCoding

def ENSG_dict(sFile):
	dCosmicCoding=dict()
	fGeneANNOVAR=open(sFile,"r")

	for sLine in fGeneANNOVAR.xreadlines():
	#	print sLine
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sType,sChr, nPosition, sRef, sAlt)=(t[0],t[2], t[3], t[5], t[6])
		(sENSG)=t[1]
		
		if "(" in sENSG:
			sENSG=sENSG.split("(")[0]
		
		if sType=="intergenic":
			sENSG=''
		else:
			pass
		
		
		sKey="|".join([sChr,nPosition,sRef,sAlt])
		dCosmicCoding[sKey]=sENSG
		
		
	return dCosmicCoding




	
	
def Cytoband(sFile):
	dCosmicCoding=dict()
	fGeneANNOVAR=open(sFile,"r")

	for sLine in fGeneANNOVAR.xreadlines():
	#	print sLine
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr, nPosition, sRef, sAlt)=(t[2], t[3], t[5], t[6])
		(sCosmicID)=t[1]
		sKey="|".join([sChr,nPosition,sRef,sAlt])
		dCosmicCoding[sKey]=sCosmicID
		
		
	return dCosmicCoding
	


def HGVS(sKey):
	
	(sChr, nPosition, sRef, sAlt)=sKey.split("|")
	
	
	sChr=sChr.replace("chr","NC_0000").replace("X","23").replace("Y","24").replace("MT","M")
	
	
	print (sChr, nPosition, sRef, sAlt)
	
	if "-" in sRef:
		return sChr+".10:g."+str(nPosition)+"_"+str(int(nPosition)+len(sAlt))+"ins"+sAlt
	elif "-" in sAlt:
		return sChr+".10:g."+str(nPosition)+"_"+str(int(nPosition)+len(sAlt))+"del"+sRef
	else:
		return sChr+".10:g."+str(nPosition)+sRef+">"+sAlt
	
	


	




def Converter(sFile):
	
	
	n=1
	
	dIntron=Intron_Dict("./ANNOVAR_Temporaly/Result_Singleton_Intersected_Normalize_Decompose_"+sFile+".variant_function")
	dExon=ExonicDict("./ANNOVAR_Temporaly/Result_Singleton_Intersected_Normalize_Decompose_"+sFile+".exonic_variant_function")
	dCosmicCoding=CosmicCoding("./ANNOVAR_Temporaly/Result_Cosmic_coding_Intersected_Normalize_Decompose_"+sFile+".hg19_cosmic83_coding_dropped")
	dCosmicNoncoding=CosmicCoding("./ANNOVAR_Temporaly/Result_Cosmic_noncoding_Intersected_Normalize_Decompose_"+sFile+".hg19_cosmic83_noncoding_dropped")
	dCytoband=Cytoband("./ANNOVAR_Temporaly/Result_CytoBand_Intersected_Normalize_Decompose_"+sFile+".hg19_cytoBand")
	dClinvar=Clinvar("./ANNOVAR_Input/For_ANNOVARIntersected_Normalize_Decompose_"+sFile)
	dENSG=ENSG_dict("./ANNOVAR_Temporaly/Result_ENSG_Intersected_Normalize_Decompose_"+sFile+".variant_function")
	
	dDBsnp=Cytoband("./ANNOVAR_Temporaly/Result_dbSNP_Intersected_Normalize_Decompose_"+sFile+".hg19_avsnp150_dropped")
	
	
	fout=open("SNV_INDEL_"+sFile,"w")
	
	fout.write("HGVS genomic change\tHGVS coding change\tHGVS protein change\tHGVS version\tGenome build\tGenomic source\tHGNC gene symbol\tEntrez gene ID\tEnsembl gene ID\thromosome\tPosition\tReference allele\tAlternative allele\tCytogenetic location\tStrand orientation\tCodon\tExon\tMolecular Effects\tVariant type\tGenotype\tdbSNP Identification Number\tclinVar Variation Identification Number\tCOSMIC Identification Number\tBiomarker name\tBiomarker state\tWildtype biomarker YN\tFunctional Domain\n")
	for sKey in dExon.keys():
		
		
		sChr,nPosition,sRef,sAlt=sKey.split("|")
		
		if "-" in sRef:
			sVariantType="Insertion"
		elif "-" in sAlt:
			sVariantType="Deletion"
		else:
			sVariantType="Substitution"
		
		
		lExon=dExon[sKey]
		lNMIDs=lExon[0]
		nExons=lExon[1]
		sCodons=lExon[2]
		Allele=lExon[3]
		sGene=lExon[4]
		
		
		lCodonNumber=[]
		
		if not sCodons==[]:
			for sCodon in sCodons:
				sCodon=sCodon.replace("p.","")
				sCodon=sCodon[1:-1]
				lCodonNumber.append(sCodon)
		
		
		
		
		if Allele=="1/1":
			sAlleleStatus="Homozygous"
		else:
			sAlleleStatus="Heterozygous"
		
		
		sCosmic=''
		sCytoband=''
		sClinvar=''
		sDBsnp=''
		
		
		if sKey in dCosmicCoding.keys():
			sCosmic=dCosmicCoding[sKey]
		
		if sKey in dCosmicNoncoding.keys():
			sCosmic=dCosmicNoncoding[sKey]
			
		if sKey in dClinvar.keys():
			sClinvar=dClinvar[sKey]
		
		if sKey in dCytoband.keys():
			sCytoband=dCytoband[sKey]
		
		if sKey in dDBsnp.keys():
			sDBsnp=dDBsnp[sKey]
		
		
		
		fout.write("{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{27}\t{28}\t{29}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\n".\
		format("SI_"+str(n),"",HGVS(sKey),",".join(lNMIDs), ",".join(sCodons), "15.11","GRCh37/hg19","U",sGene, '',dENSG[sKey],sChr ,sCytoband, '0',",".join(lCodonNumber) ,",".join(nExons),"Missense", sVariantType, sAlleleStatus,sDBsnp,"",sClinvar,sCosmic,'','','','',nPosition,sRef,sAlt))
		
		n+=1
	
	
	
	
	for sKey in dIntron.keys():
		
		
		sChr,nPosition,sRef,sAlt=sKey.split("|")
		
		if "-" in sRef:
			sVariantType="Insertion"
		elif "-" in sAlt:
			sVariantType="Deletion"
		else:
			sVariantType="Substitution"
		
		
		lIntron=dIntron[sKey]
		#[sType,sGene,sInfo]
		
		sType=lIntron[0]
		sGene=lIntron[1]
		Allele=lIntron[2]
		
		if Allele=="1/1":
			sAlleleStatus="Homozygous"
		else:
			sAlleleStatus="Heterozygous"
		
		
		
		
		
		
		sCosmic=''
		sCytoband=''
		sClinvar=''
		sDBsnp=''
		
		if sKey in dCosmicCoding.keys():
			sCosmic=dCosmicCoding[sKey]
		
		if sKey in dCosmicNoncoding.keys():
			sCosmic=dCosmicNoncoding[sKey]
			
		if sKey in dClinvar.keys():
			sClinvar=dClinvar[sKey]
		
		if sKey in dCytoband.keys():
			sCytoband=dCytoband[sKey]
		
		if sKey in dDBsnp.keys():
			sDBsnp=dDBsnp[sKey]
		
		
		
		fout.write("{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{27}\t{28}\t{29}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\n".\
		format("SI_"+str(n),"",HGVS(sKey),"", "", "15.11","GRCh37/hg19","U",sGene, '',dENSG[sKey],sChr ,sCytoband, '0',"" ,"","", sVariantType, sAlleleStatus,sDBsnp,"",sClinvar,sCosmic,'','','','',nPosition,sRef,sAlt))
		n+=1
	
	
	
	
	






def producer_task(q, cosmic_dict):
	sTarget001list=glob.glob("Filter*.vcf")
	sFilelist=sTarget001list
	for i in sFilelist:
		value=i
		cosmic_dict[value]=None
	
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)


def consumer_task(q, cosmic_dict):
	sTarget="/storage/home/SNUH/3068381_Covered.bed"
	while not q.empty():
		value=q.get(True, 0.05)
		
		Converter(value)
		
		
		
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
#	number_of_cpus=cpu_count()-2
	number_of_cpus=8
	manager=Manager()
	fibo_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, fibo_dict))
	producer.start()
	producer.join()
	consumer_list=[]
	for i in range(number_of_cpus):
		consumer=Process(target=consumer_task, args=(data_queue,fibo_dict))
		consumer.start()
		consumer_list.append(consumer)

	[consumer.join() for consumer in consumer_list]

	logger.info(fibo_dict)
	
	
	
	
	
	print "Start Time"
	print StartTime
	print "End Time"
	print(time.ctime())





