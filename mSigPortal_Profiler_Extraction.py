#!/usr/bin/python
#encoding=utf8
import re,os,argparse,sys,time
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import sigProfilerPlotting as sigPlt

'''
Name:       mSigPortal_Format_Convert
Function:	Generate Input File for mSigPortal
Version:    8.0
Date:       June-17-2020
Update:     (1) Add Error 233: A indicator for format Error
			(2) Add sigProfilerPlotting to generate PDF and SVG
			(3) Update sigProfilerPlotting
			(4) SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGeneratorFunc.py:
				Comment the line of 312, or a exception will happen
				#log_out.write("SigProfilerPlotting version: "+sigPlt.__version__+"\n")
			(5) Catelog results should be tsv, no matter if input is csv or tsv
			(6) Default Output:'mSigPortal_Project_%s' % time.strftime('%Y%m%d%H%M%S',time.localtime(time.time())) (+ ProjectID will be updated later)
			(7) Add -b option for Bed file in SigProfilerMatrixGenerator function
'''

########################################################################
###################### 0 Define Basic Function #########################
########################################################################


####### 01-1 Get Options
def Parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--Input_Format', required=True, type=str, help = "Define the formats of input data. Only 'vcf', 'csv', 'tsv', 'catalog_tsv', 'catalog_csv' are supported.")
	parser.add_argument('-i', '--Input_Path', required=True, type=str, help = "Define the absolute path for input file")
	parser.add_argument('-p', '--Project_ID', required=True, type=str, help = "Define the ID of Project")
	parser.add_argument('-o', '--Output_Dir', required=True, nargs='?', const='mSigPortal_%s' % time.strftime('%Y%m%d%H%M%S',time.localtime(time.time())),type=str, help = "Define the absolute path for Output Dir")
	parser.add_argument('-g', '--Genome_Building', required=True, nargs='?', const='GRCh37', type=str, help = "Define the version of Genome_Building, default:GRCh37 ")
	parser.add_argument('-t', '--Data_Type', required=True, type=str, help = "Define the Data_Type: 'WGS' or 'WES' ")
	parser.add_argument('-F', '--Filter', required=False, type=str, help="Define the Terms for filtring. Different terms should be seperated by \'@\' ")
	parser.add_argument('-c', '--Collapse', required=False, type=str, nargs='?', const='Sample_Collpase', help="Whether collapses multiple samples into one group, If yes: add '-c Group_Name'; Default:'Sample_Collpase'")
	parser.add_argument('-z', '--gzip', required=False, type=str, nargs='?', help="Whether gzip results? samples into one group, If yes: add '-z True'")
	parser.add_argument('-s', '--vcf_split_all_filter', required=False, type=str, nargs='?', help="Whether split all vcf filter? If yes: add '-s True'")
	parser.add_argument('-b', '--Bed', required=False, type=str, nargs='?', help="Bed File for SigProfilerMatrixGenerator")


	args = parser.parse_args()
	return args.Input_Format,args.Input_Path,args.Project_ID,args.Output_Dir,args.Genome_Building,args.Data_Type,args.Filter,args.Collapse,args.gzip,args.vcf_split_all_filter,args.Bed

####### 01-14 If input file is compressed?
def If_Compressed():
	### 000 Parse Options
	Input_Format,Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse,gzip,vcf_split_all_filter,Bed = Parser()
	
	Input_Path_New_Name = Input_Path
	
	
	### 001 if Output_Dir exists, delete it!
	if os.path.exists(Output_Dir):
		os.system("rm -rf %s" % Output_Dir)
	GenerateDir(Output_Dir) 
	
	### 002 if in 3 types compressed files?
	if re.search(r'zip$',Input_Path):
		String = "unzip %s -d %s" % (Input_Path,Output_Dir) 
		print(String)
		os.system(String)
		name = Input_Path.split("/")[1].split(".zip")[0]
		Input_Path_New_Name = "%s/%s" % (Output_Dir,name)
		print(Input_Path_New_Name)

	if re.search(r'tar$',Input_Path):
		String = "tar -zxvf %s -C %s" % (Input_Path,Output_Dir) 
		print(String)
		os.system(String)
		name = Input_Path.split("/")[1].split(".tar")[0]
		Input_Path_New_Name = "%s/%s" % (Output_Dir,name)
		print(Input_Path_New_Name)
		
	if re.search(r'gz$',Input_Path):
		if "tar.gz" not in Input_Path:
			name = Input_Path.split("/")[1].split(".gz")[0]
			String = "gunzip -c %s > %s/%s" % (Input_Path,Output_Dir,name) 
			print(String)
			os.system(String)

			Input_Path_New_Name = "%s/%s" % (Output_Dir,name)
			print(Input_Path_New_Name)
		else:
			String = "tar -zxvf %s -C %s" % (Input_Path,Output_Dir) 
			print(String)
			os.system(String)
			name = Input_Path.split("/")[1].split(".tar")[0]
			Input_Path_New_Name = "%s/%s" % (Output_Dir,name)
			print(Input_Path_New_Name)
	#print(Input_Path)
	return Input_Path_New_Name




####### 01-2 Parse Options
####### Choose Sub Function according to Options
def Parse_Options():
	Input_Format,Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse,gzip,vcf_split_all_filter,Bed = Parser()
	
	###### Parse Options 001 Check Compressed.Generate New Input_Path, still named as Input_Path
	Input_Path = If_Compressed()

	###### Parse Options 002 Choose Transforming Function based on Input_Format
	if Input_Format == "vcf_Single":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports vcf format")
			sys.exit()

		if Collapse == None:
			pass
		else:
			String = "Error, \"Collapse\" option doese not support \"%s\" format " % (Input_Format)
			print(String)
			sys.exit()

		if Filter == None:
			vcf_Single_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
		else:
			vcf_Single_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter)


	elif Input_Format == "vcf_Multiple":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports vcf format")
			sys.exit()

		if Filter == None:
			vcf_Multiple_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
		else:
			vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)

		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == None:
			pass
		else:
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)



	elif Input_Format == "csv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports vcf format")
			sys.exit()

		###### Note Default of an option is None,not ''
		if Filter == None:
			csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
		else:
			csv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)

		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == None:
			pass
		else:
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)


	elif Input_Format == "tsv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports vcf format")
			sys.exit()

		###### Note Default of an option is None,not ''
		if Filter == None:
			tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
		else:
			tsv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)

		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == None:
			pass
		else:
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)



	elif Input_Format == "catalog_tsv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports vcf format")
			sys.exit()

		###### Note Default of an option is None,not ''
		if Collapse == None:
			pass
		else:
			String = "Error, \"-c/--Collapse\" option doese not support \"%s\" format " % (Input_Format)
			print(String)
			sys.exit()

		if Filter == None:
			catalog_tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
		else:
			String = "Error, \"%s\" format do not support option \"-F/--Filter\"" % (Input_Format)
			print(String)



	elif Input_Format == "catalog_csv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports vcf format")
			sys.exit()
		
		###### Note Default of an option is None,not ''
		if Collapse == None:
			pass
		else:
			String = "Error, \"Collapse\" option doese not support \"%s\" format " % (Input_Format)
			print(String)
			sys.exit()

		if Filter == None:
			catalog_csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
		else:
			String = "Error, \"%s\" format do not support option \"-F/--Filter\"" % (Input_Format)
			print(String)
			sys.exit()


	elif Input_Format == "vcf":
		if vcf_split_all_filter == "True":
			if Filter == None:
				vcf_Multiple_Convert_Split_All_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
			else:
				print("Error, You have spcified '-s True', in this case the option -F is not supported!")
				sys.exit()

		else:
			if Filter == None:
				vcf_Multiple_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
			else:
				vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)

		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == None:
			pass
		else:
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)

	else:
		String = "Error in Format: Only 'vcf', 'csv', 'tsv', 'catalog_csv', 'catalog_tsv' are validated formats!"
		print(String)
		sys.exit()

	###### Parse Options 003 Print Basic Output Statistic.
	Print_Statistic(Output_Dir)


	###### Parse Options 004 Print Basic Output Statistic.
	sigProfilerPlotting(Input_Format,Output_Dir,Project_ID,Genome_Building,Bed)


	###### Parse Options 005 Checke if gzip is on.
	## 004 Check gzip:
	if gzip == None:
		pass
	elif gzip == "True":
		gzip_Output(Output_Dir)
	elif gzip == "False":
		pass
	else:
		String = "Error: Cannot Compress Output File. -z/--gzip option can only follow 'True' or 'False'"
		print(String)
		sys.exit()


####### 01-3 vcf_Single_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def vcf_Single_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in VCF format with single sample. ******* ')
	
	####### 01-3-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Sample_ID = "Sample"
	
	for line in Input_File:
		if re.match(r'##',line):
			pass
		elif re.match(r'#CHROM',line):
			ss = line.strip().split("	")
			Sample_ID = ss[-1]
		else:
			ss = line.strip().split("	")
			REF = ss[3]
			ALT = ss[4]
			Chr = ss[0]
			Start = ss[1]
			
			if len(REF) == len(ALT):
				End = ss[1]
				Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
				mSigPortal_Format_SNV_File.write(Output_String)
				#print(Output_String)
			else:
				if len(REF) == 1:
					Start = ss[1]
					End = ss[1]
				else:
					Start = ss[1]
					End = int(Start) + len(REF) - 1
				Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
				mSigPortal_Format_INDEL_File.write(Output_String)


	Input_File.close()
	
	#print(Sample_ID)
	mSigPortal_Format_INDEL_File.close()
	mSigPortal_Format_SNV_File.close()


####### 01-4 GenerateDir
def GenerateDir(Dir):
	if not os.path.exists(Dir):
		os.system("mkdir %s" % Dir)

####### 01-5 csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in CSV format. ******* ')


	####### 01-3-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Header = "Sample_ID,Chrom,Start,End,REF,ALT,Filter"
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split(",")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The CSV option requires each line in input file should include 7 Item: Sample,Chrom,Start,End,Ref,Alt,Filter.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"Sample_ID,Chrom,Start,End,REF,ALT,Filter\" is required!")
		sys.exit()
		
	####### 01-3-2 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'Sample_ID',f):
			pass
		else:
			ss = f.split(",")
			if len(ss) == 7:
				Sample_ID = ss[0]
				Chr = ss[1]
				Start = ss[2]
				End = ss[3]
				REF = ss[4]
				ALT = ss[5]
				if len(REF) == len(ALT):
					Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_SNV_File.write(Output_String)
				else:
					Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_INDEL_File.write(Output_String)
	Input_File.close()


####### 01-6 csv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)
def csv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in CSV format, with filtration: %s ******* ' % Filter)

	####### 01-6-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-6-1 Parse File 
	Input_File = open(Input_Path)
	Header = "Sample_ID,Chrom,Start,End,REF,ALT,Filter"
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split(",")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The CSV option requires each line in input file should include 7 Item: Sample,Chrom,Start,End,Ref,Alt,Filter.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"Sample_ID,Chrom,Start,End,REF,ALT,Filter\" is required!")
		sys.exit()

	####### 01-6-2 Parse Filter:
	option_Filter_Arr = Filter.split("@")

	####### 01-6-3 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'Sample_ID',f):
			pass
		else:
			ss = f.split(",")
			if len(ss) == 7:
				Sample_ID = ss[0]
				Chr = ss[1]
				Start = ss[2]
				End = ss[3]
				REF = ss[4]
				ALT = ss[5]
				Filter = ss[6].strip()
				if Filter in option_Filter_Arr:
					Sample_ID = "%s@%s" % (Sample_ID,Filter)
					if len(REF) == len(ALT):
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)
	Input_File.close()


####### 01-7 tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in TSV format. ******* ')


	####### 01-3-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Header = "Sample_ID	Chrom	Start	End	Ref	Alt	Filter"
	#Header = ""
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split("	")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The TSV format requires each line in input file should include 7 Item: Sample_ID	Chrom	Start	End	Ref	Alt	Filter.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"Sample_ID	Chrom	Start	End	Ref	Alt	Filter\" is required!")
		sys.exit()
	#print(String_File)
	####### 01-3-2 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'Sample_ID',f):
			pass
		else:
			ss = f.split("	")
			if len(ss) == 7:
				Sample_ID = ss[0]
				Chr = ss[1]
				Start = ss[2]
				End = ss[3]
				REF = ss[4]
				ALT = ss[5].strip()
				if len(REF) == len(ALT):
					Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_SNV_File.write(Output_String)
				else:
					Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_INDEL_File.write(Output_String)

	Input_File.close()



####### 01-8 tsv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)
def tsv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in TSV format, with filtration: %s ******* ' % Filter)

	####### 01-6-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-6-1 Parse File 
	Input_File = open(Input_Path)
	Header = "Sample_ID	Chrom	Start	End	Ref	Alt	Filter"
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split("	")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The CSV optformat requires each line in input file should include 7 Item: Sample_ID	Chrom	Start	End	Ref	Alt	Filter.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"Sample_ID	Chrom	Start	End	Ref	Alt	Filter\" is required!")
		sys.exit()

	####### 01-6-2 Parse Filter:
	option_Filter_Arr = Filter.split("@")

	####### 01-6-3 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'Sample_ID',f):
			pass
		else:
			ss = f.split("	")
			if len(ss) == 7:
				Sample_ID = ss[0]
				Chr = ss[1]
				Start = ss[2]
				End = ss[3]
				REF = ss[4]
				ALT = ss[5]
				Filter = ss[6].strip()
				if Filter in option_Filter_Arr:
					Sample_ID = "%s@%s" % (Sample_ID,Filter)
					if len(REF) == len(ALT):
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)
	Input_File.close()



####### 01-9 vcf_Multple_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def vcf_Multiple_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in VCF format. ******* ')
	
	####### 01-3-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	
	####### 01-3-20 Check VCF fromat
	header = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
	String_File = ""
	
	Input_File = open(Input_Path)
	for line in Input_File:
		String_File += line
	Input_File.close()

	if header not in String_File:
		print("Error 233: A header line: \"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1...\" is required!")
		sys.exit()


	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Sample_ID = "Sample"
	
	for line in Input_File:
		if re.match(r'##',line):
			pass
		elif re.match(r'#CHROM',line):
			#######vcf_Multiple Sample id will be stored in a arr
			ss = line.strip().split("	")
			Sample_ID_arr = ss[9:len(ss)]
		else:
			ss = line.strip().split("	")
			REF = ss[3]
			ALT = ss[4]
			Chr = ss[0]
			Start = ss[1]
			
			####### 对每个样本 Sample_ID_arr中的成员进行便利，判断是否有基因型信息
			for i in range(0,len(Sample_ID_arr)):
				Genotype_String = ss[9+i]
				gg = Genotype_String.split(":")[0]
				if re.search(r'1',gg):
					#print i,Sample_ID_arr[i]
					Sample_ID = Sample_ID_arr[i]
					
					if len(REF) == len(ALT):
						End = ss[1]
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
						#print(Output_String)
					else:
						if len(REF) == 1:
							Start = ss[1]
							End = ss[1]
						else:
							Start = ss[1]
							End = int(Start) + len(REF) - 1
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)
	
	print(Sample_ID_arr)
	mSigPortal_Format_INDEL_File.close()
	mSigPortal_Format_SNV_File.close()


####### 01-10 vcf_Multiple_Convert_Split_All_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def vcf_Multiple_Convert_Split_All_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in VCF format. ******* ')
	
	####### 01-3-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Sample_ID = "Sample"
	
	for line in Input_File:
		if re.match(r'##',line):
			pass
		elif re.match(r'#CHROM',line):
			#######vcf_Multiple Sample id will be stored in a arr
			ss = line.strip().split("	")
			Sample_ID_arr = ss[9:len(ss)]
		else:
			ss = line.strip().split("	")
			REF = ss[3]
			ALT = ss[4]
			Chr = ss[0]
			Start = ss[1]
			Filter_arr = ss[6].split(";")
			
			####### 对每个样本 Sample_ID_arr中的成员进行便利，判断是否有基因型信息
			for i in range(0,len(Sample_ID_arr)):
				Genotype_String = ss[9+i]
				gg = Genotype_String.split(":")[0]
				if re.search(r'1',gg):
					#print i,Sample_ID_arr[i]
					Sample_ID = Sample_ID_arr[i]
					
					for i in Filter_arr:
						Sample_ID_Filter = "%s@%s" % (Sample_ID,i)
						
						if len(REF) == len(ALT):
							End = ss[1]
							Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Filter,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
							mSigPortal_Format_SNV_File.write(Output_String)
							#print(Output_String)
						else:
							if len(REF) == 1:
								Start = ss[1]
								End = ss[1]
							else:
								Start = ss[1]
								End = int(Start) + len(REF) - 1
							Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Filter,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
							mSigPortal_Format_INDEL_File.write(Output_String)
	
	print(Sample_ID_arr)
	mSigPortal_Format_INDEL_File.close()
	mSigPortal_Format_SNV_File.close()





####### 01-11 catalog_tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def catalog_tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in CATALOG tsv format. ******* ')

	####### 01-3-0 Output_Path:

	mSigPortal_Format_catalog_Path = "%s/%s_mSigPortal_catalog_tsv.txt" % (Output_Dir,Project_ID)
	mSigPortal_Format_catalog_File = open(mSigPortal_Format_catalog_Path,'w')

	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Header = "MutationType	"
	String_File = ""

	for line in Input_File:
		if re.match(r'\n',line):
			pass
		else:
			String_File += line
		
	if Header not in String_File:
		print("Error 233: A header line like \"MutationType	Sample1	Sample2	Sample3...\" is required!")
		sys.exit()
		
	####### 01-3-2 Generate Result
	ff = String_File.split("\n")
	Header = ""
	for f in ff:
		if re.match(r'MutationType	',f):
			Header = f
			mSigPortal_Format_catalog_File.write("%s	Count\n" % Header)

		else:
			ss = f.split("	")
			if len(ss) > 1:
				Count_Total = 0
				Count_Arr = ss[1:len(ss)]
				for i in Count_Arr:
					Count_Total = Count_Total + int(i)
				
				Output_String = "%s	%s\n" % (f,Count_Total)
				mSigPortal_Format_catalog_File.write(Output_String)

	Input_File.close()
	mSigPortal_Format_catalog_File.close()


####### 01-12 catalog_csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def catalog_csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in CATALOG csv format. ******* ')

	####### 01-3-0 Output_Path:

	mSigPortal_Format_catalog_Path = "%s/%s_mSigPortal_catalog_csv.txt" % (Output_Dir,Project_ID)
	mSigPortal_Format_catalog_File = open(mSigPortal_Format_catalog_Path,'w')

	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Header = "MutationType,"
	String_File = ""
	#Count = 1
	for line in Input_File:
		if re.match(r'\n',line):
			pass
		else:
			String_File += line
		#print(line)
	if Header not in String_File:
		print("Error 233: A header line like \"MutationType,Sample1,Sample2,Sample3...\" is required!")
		sys.exit()
	#mSigPortal_Format_catalog_File.write("ss")
	####### 01-3-2 Generate Result
	ff = String_File.split("\n")
	Header = ""
	for f in ff:
		if re.match(r'MutationType,',f):
			Header = f
			mSigPortal_Format_catalog_File.write("%s	Count\n" % Header.replace(",","\t"))
			#print(f)
		elif re.match(r'\n',f):
			pass
		else:
			ss = f.split(",")
			Count_Total = 0
			if len(ss) > 1:
				Count_Arr = ss[1:len(ss)]
				for i in Count_Arr:
					Count_Total = Count_Total + int(i)
				
				Output_String = "%s	%s\n" % (f.replace(",","\t"),Count_Total)
				mSigPortal_Format_catalog_File.write(Output_String)

	Input_File.close()
	mSigPortal_Format_catalog_File.close()

	#print("ssss")


####### 01-13 vcf_Single_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter)
def vcf_Single_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter):

	GenerateDir(Output_Dir)
	print("******* Your Input File is in VCF format. Your Filter term is: %s ******* " % Filter)
	
	####### 01-10-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-10-1 Parse Filter:
	option_Filter_Arr = Filter.split("@")

	
	####### 01-10-2 Parse File 
	Input_File = open(Input_Path)
	Sample_ID = "Sample"
	
	for line in Input_File:
		if re.match(r'##',line):
			pass
		elif re.match(r'#CHROM',line):
			ss = line.strip().split("	")
			Sample_ID = ss[-1]
		else:
			ss = line.strip().split("	")
			REF = ss[3]
			ALT = ss[4]
			Chr = ss[0]
			Start = ss[1]
			vcf_Filter_arr = ss[6].split(";")
			
			if len(REF) == len(ALT):
				End = ss[1]
				
				### Pay attention to option "-"
				if "-" in option_Filter_Arr:
					Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_SNV_File.write(Output_String)
				for x in vcf_Filter_arr:
					if x in option_Filter_Arr:
						Sample_ID_Final = "%s@%s" % (Sample_ID,x)
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
				#print(Output_String)
			else:
				if len(REF) == 1:
					Start = ss[1]
					End = ss[1]
				else:
					Start = ss[1]
					End = int(Start) + len(REF) - 1

				### Pay attention to option "-"
				if "-" in option_Filter_Arr:
					Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_INDEL_File.write(Output_String)

				for x in vcf_Filter_arr:
					if x in option_Filter_Arr:
						Sample_ID_Final = "%s@%s" % (Sample_ID,x)
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)

	Input_File.close()
	
	print(Sample_ID)
	mSigPortal_Format_INDEL_File.close()
	mSigPortal_Format_SNV_File.close()


####### 01-14 vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter)
def vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in VCF format with Multple sample. Your Filter term is: %s ******* ' % (Filter))
	
	####### 01-3-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)


	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-10-1 Parse Filter:
	option_Filter_Arr = Filter.split("@")
	
	
	####### 01-3-1 Parse File 
	Input_File = open(Input_Path)
	Sample_ID = "Sample"
	
	for line in Input_File:
		if re.match(r'##',line):
			pass
		elif re.match(r'#CHROM',line):
			#######vcf_Multiple Sample id will be stored in a arr
			ss = line.strip().split("	")
			Sample_ID_arr = ss[9:len(ss)]
		else:
			ss = line.strip().split("	")
			REF = ss[3]
			ALT = ss[4]
			Chr = ss[0]
			Start = ss[1]
			
			vcf_Filter_arr = ss[6].split(";")


			### Pay attention to option "-"
			if "-" in option_Filter_Arr:
				if len(REF) == len(ALT):
					End = ss[1]
					Output_String = "%s	Sample	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_SNV_File.write(Output_String)
				else:
					if len(REF) == 1:
						Start = ss[1]
						End = ss[1]
					else:
						Start = ss[1]
						End = int(Start) + len(REF) - 1
					Output_String = "%s	Sample	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
					mSigPortal_Format_INDEL_File.write(Output_String)


			####### 对每个样本 Sample_ID_arr中的成员进行便利，判断是否有基因型信息
			for i in range(0,len(Sample_ID_arr)):
				Genotype_String = ss[9+i]
				gg = Genotype_String.split(":")[0]
				if re.search(r'1',gg):
					#print i,Sample_ID_arr[i]
					Sample_ID = Sample_ID_arr[i]
					
					if len(REF) == len(ALT):
						End = ss[1]
						
						for x in vcf_Filter_arr:
							if x in option_Filter_Arr:
								Sample_ID_Final = "%s@%s" % (Sample_ID,x)
								Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
								mSigPortal_Format_SNV_File.write(Output_String)

					else:
						if len(REF) == 1:
							Start = ss[1]
							End = ss[1]
						else:
							Start = ss[1]
							End = int(Start) + len(REF) - 1

						for x in vcf_Filter_arr:
							if x in option_Filter_Arr:
								Sample_ID_Final = "%s@%s" % (Sample_ID,x)
								Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
								mSigPortal_Format_INDEL_File.write(Output_String)
	
	print(Sample_ID_arr)
	mSigPortal_Format_INDEL_File.close()
	mSigPortal_Format_SNV_File.close()



####### 01-15 Convert_Collapse(Output_Dir)
def Convert_Collapse(Output_Dir,Collapse,Project_ID):
	
	####### 01-12-0 Input_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	####### 01-12-1 Output_Path:
	mSigPortal_Format_SNV_Collapse_Path = "%s/%s_mSigPortal_SNV_Collapse.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Collapse_Path = "%s/%s_mSigPortal_INDEL_Collapse.txt" % (Output_Dir,Project_ID)

	####### 01-12-2 Count the Sample in Input_Path: mSigPortal_Format_SNV_File
	SNV_Sample_Count = []
	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path)
	for line in mSigPortal_Format_SNV_File:
		ss = line.strip().split("	")
		Sample_Name = ss[1].split("@")[0]
		SNV_Sample_Count.append(Sample_Name)
	mSigPortal_Format_SNV_File.close()
	
	count = len(SNV_Sample_Count)
	String_Collpase = "There are %d samples in the original Output File:%s" % (count,mSigPortal_Format_SNV_Path)
	
	####### Only the Input file with more than 1 sample can support “Collapse Option” !

	####### 01-12-3 Parse Input_Path: mSigPortal_Format_SNV_File
	if count > 1:
		mSigPortal_Format_SNV_Collapse_File = open(mSigPortal_Format_SNV_Collapse_Path,'w')
		mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path)
		for line in mSigPortal_Format_SNV_File:
			ss = line.strip().split("	")
			Sample_Name = ss[1].split("@")[0]
			String_1 = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			String_2 = "%s	%s_Collapse	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],Sample_Name,ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			mSigPortal_Format_SNV_Collapse_File.write(String_1)
			mSigPortal_Format_SNV_Collapse_File.write(String_2)
		mSigPortal_Format_SNV_File.close()
		mSigPortal_Format_SNV_Collapse_File.close()

	else:
		print("Error: Only the Input file with more than 1 sample can support “Collapse Option”")


	####### 01-12-4 Count the Sample in Input_Path: mSigPortal_Format_INDEL_File
	INDEL_Sample_Count = []

	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path)
	for line in mSigPortal_Format_INDEL_File:
		ss = line.strip().split("	")
		Sample_Name = ss[1].split("@")[0]
		INDEL_Sample_Count.append(Sample_Name)
	count = len(INDEL_Sample_Count)
	String_Collpase = "There are %d samples in the original Output File:%s" % (count,mSigPortal_Format_SNV_Path)
	mSigPortal_Format_INDEL_File.close()


	####### 01-12-5 Parse Input_Path: mSigPortal_Format_INDEL_File
	if count > 1:
		mSigPortal_Format_INDEL_Collapse_File = open(mSigPortal_Format_INDEL_Collapse_Path,'w')
		mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path)
		for line in mSigPortal_Format_INDEL_File:
			ss = line.strip().split("	")
			Sample_Name = ss[1].split("@")[0]
			String_1 = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			String_2 = "%s	%s_Collapse	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],Sample_Name,ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			mSigPortal_Format_INDEL_Collapse_File.write(String_1)
			mSigPortal_Format_INDEL_Collapse_File.write(String_2)
		mSigPortal_Format_INDEL_Collapse_File.close()
		mSigPortal_Format_INDEL_File.close()
	else:
		print("Error: Only the Input file with more than 1 sample can support “Collapse Option”")


	####### 01-6-0 RM mSigPortal_Format_SNV_Path and mSigPortal_Format_INDEL_Path
	cmd_1 = "rm %s" % (mSigPortal_Format_SNV_Path)
	cmd_2 = "rm %s" % (mSigPortal_Format_INDEL_Path)

	os.system(cmd_1)
	os.system(cmd_2)



####### 01-16 gzip_Output(Output_Dir)
def gzip_Output(Output_Dir):
	print("Compress Output File")
	####### 01-16-0 Input_Path_arr:
	
	Input_Path_arr = ["mSigPortal_SNV.txt","mSigPortal_INDEL.txt","mSigPortal_SNV_Collapse.txt","mSigPortal_INDEL_Collapse.txt","mSigPortal_catalog_csv.txt","mSigPortal_catalog_tsv.txt"]

	for OutputFile in os.listdir(Output_Dir):
		for i in Input_Path_arr:
			if re.search(r'%s$' % i,OutputFile):
				Output_Path = "%s/%s" % (Output_Dir,OutputFile)
				cmd_String = "gzip %s" % (Output_Path)
				os.system(cmd_String)

####### 01-17 sigProfilerPlotting
def sigProfilerPlotting(Input_Format,Output_Dir,Project_ID,Genome_Building,Bed):
	#print(Bed)
	
	Input_Format_arr_1 = ['vcf', 'csv', 'tsv']
	Input_Format_arr_2 = ['catalog_csv', 'catalog_tsv']
	
	SBS_Arr = [6,24,96,384,1536,6144]
	ID_Arr = [28,83,415,8268]
	DBS_Arr = [78,186,312,1248,2976]
	
	# ####### Which format is the input file
	if Input_Format in Input_Format_arr_1:
		matrices = matGen.SigProfilerMatrixGeneratorFunc(Project_ID, Genome_Building, Output_Dir, exome=False, bed_file=Bed, chrom_based=False, plot=True, tsb_stat=False, seqInfo=False)

	elif Input_Format in Input_Format_arr_2:
		for matrix_name in os.listdir(Output_Dir):
			count = 0
			matrix_path = "%s/%s" % (Output_Dir,matrix_name)
			matrix_File = open(matrix_path)
			for line in matrix_File:
				ss = line.split("	")
				if len(ss) > 1:
					count += 1
			Type = count-1
			matrix_File.close()
			
			
			###### Plotting the Matrix Based on Type
			Final_output_Dir = "%s/" % (Output_Dir)
			Final_Type = "%d" % Type
			if Type in SBS_Arr:
				sigPlt.plotSBS(matrix_path, Final_output_Dir, Project_ID, Final_Type, percentage=False)

			elif Type in DBS_Arr:
				sigPlt.plotDBS(matrix_path, Final_output_Dir, Project_ID, Final_Type, percentage=False)

			elif Type in ID_Arr:
				sigPlt.plotID(matrix_path, Final_output_Dir, Project_ID, Final_Type, percentage=False)

			else:
				print("Error 233: Your input type is not supported yet!")
				sys.exit()


####### 01-18 Print_Statistic(Output_Dir)
def Print_Statistic(Output_Dir):
	Input_Path_arr = ["mSigPortal_SNV.txt","mSigPortal_INDEL.txt","mSigPortal_SNV_Collapse.txt","mSigPortal_INDEL_Collapse.txt","mSigPortal_catalog_csv.txt","mSigPortal_catalog_tsv.txt"]

	for OutputFile in os.listdir(Output_Dir):
		for i in Input_Path_arr:
			if re.search(r'%s$' % i,OutputFile):
				Output_Path = "%s/%s" % (Output_Dir,OutputFile)
				arr = []
				Output_File = open(Output_Path)
				for line in Output_File:
					arr.append(line)
				Output_File.close()
				String = "There are %d items in the OutputFile: %s" % (len(arr),Output_Path)
				print(String)
				Output_File.close()


if __name__ == "__main__":
	#If_Compressed()
	
	Parse_Options()


########################################################################
###########################  Usage Example #############################
########################################################################


### Usage for csv
# python mSigPortal_Format_Convert.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for tsv
# python mSigPortal_Format_Convert.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for catalog_csv
# python mSigPortal_Format_Convert.py -f catalog_csv -i Demo_input/demo_input_catalog.csv -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for catalog_tsv
# python mSigPortal_Format_Convert.py -f catalog_tsv -i Demo_input/demo_input_catalog.tsv -p Project -o Test_Output -g GRCh37 -t WGS

### Collpase ###
# python mSigPortal_Format_Convert.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS -c True
# python mSigPortal_Format_Convert.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output -g GRCh37 -t WGS -c True
# python mSigPortal_Format_Convert.py -f vcf_Multiple -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS -c True

### Usage for vcf
# python mSigPortal_Format_Convert.py -f vcf -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for Filter
# python mSigPortal_Format_Convert.py -f vcf -F PASS -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for Compressed File
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_single.vcf.gz -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_single.vcf.tar.gz -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_single.vcf.zip -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Format_Convert.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_single.vcf.tar -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for vcf_split_all_filter File
# python mSigPortal_Format_Convert.py -f vcf -s True -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS

