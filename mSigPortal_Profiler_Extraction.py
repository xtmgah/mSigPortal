#!/usr/bin/python
#encoding=utf8
import re,os,argparse,sys,time
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import sigProfilerPlotting as sigPlt
from zipfile import ZipFile
import pandas as pd

'''
Name:		mSigPortal_Profiler_Extraction
Function:	Generate Input File for mSigPortal
Version:	1.36
Date:		Sep-20-2022
Update:		
			(17) tar compressing without directory structure, This is very complicated better with zip not gzip, command is following:
				 cmd = "zip -jr %s/File_Dir_Name.zip %s/File_Dir_Name" % (zip_Dir,Original_Dir)
			(18) Support MAF format ["Tumor_Sample_Barcode", "Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
			(19) Enable sigPlt to support percentage
			(20) Support R32 and CNV48 for catalog_TSV and catalog_CSV
			(21) Add Cluster Function     # 06-16-2022
			(22) Fix the bug of VAF=-1.5  # 07-30-2022
			(23) Add Plotting and SeqInfo parameters # 09-20-2022
			(24) Fix the bug: If Empline file with 0 line was generated as input, delete them. # 09-20-2022
 '''


########################################################################
###################### 0 Define Basic Function #########################
########################################################################
####### 01-0 Rewrite arg_Parse help Document
Help_String = '''
Program:	mSigPortal_Profiler_Extraction
Version:	Alpha v 1.34
Function:	Extract and generate standard format for mSigPortal analysis from multiple type of input file.
Updated:	June-16-2022

Usage:	 	python3 mSigPortal_Profiler_Extraction.py [options]

Required Options:
  -f		(--Input_Format) Define the format of Input files. 
		   Only 'vcf', 'csv', 'tsv', 'maf', 'catalog_tsv', 'catalog_csv' are supported.
  -i		(--Input_Path) Absolute path for input file.
  -p		(--Project_ID) Project ID.
  -o		(--Output_Dir) Absolute path for output Directory.
		   A unique time stamp will be added automatically.
  -g		(--Genome_Building) Define the version of Genome_Building.
		   Default:GRCh37 .
  -t		(--Data_Type) Define the project source type: 
		   'WGS' for 'Whole Genome Sequencing Porject'; 'WES' for "Whole Exome Sequencing Project".
  -c        (--Cluster) Whether perform SigProfilerClusters
  -P        (--Cluster) Whether Generate Plotting
  -S        (--Cluster) Whether Generate SeqInfo

Optional Options:
  -F		(--Filter) Define the terms for variations filtration.
		   Different terms should be seperated by '@'.
  -c		(--Collapse) Do you want to collapse multiple samples into one group?
		   If yes, add '-c/--Collpase True'.
  -z		(--gzip) Do you want to compress your results?
		   If yes, add '-z/--gzip True'.
  -s		(--vcf_split_all_filter) Do you want to split all filtration terms for input variants?.
		   If yes, add '-s/--vcf_split_all_filter True'.
  -b		(--Bed) Provide a absolute path of bed file for regional-based filtration result.
		   
'''


####### 01-1 Get Options
def Parser():
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument('-f', '--Input_Format', required=True, type=str, help = "Define the formats of input data. Only 'vcf', 'csv', 'tsv', 'maf', 'catalog_tsv' and 'catalog_csv' are supported.")
	parser.add_argument('-i', '--Input_Path', required=True, type=str, help = "Define the absolute path for input file")
	parser.add_argument('-p', '--Project_ID', required=True, type=str, help = "Define the ID of Project")
	parser.add_argument('-o', '--Output_Dir', required=True, nargs='?', const='mSigPortal_%s' % time.strftime('%Y%m%d%H%M%S',time.localtime(time.time())),type=str, help = "Define the absolute path for Output Dir")
	parser.add_argument('-g', '--Genome_Building', required=True, nargs='?', const='GRCh37', type=str, help = "Define the version of Genome_Building, default:GRCh37 ")
	parser.add_argument('-t', '--Data_Type', required=True, type=str, help = "Define the Data_Type: 'WGS' or 'WES' ")
	parser.add_argument('-F', '--Filter', required=False, type=str, help="Define the Terms for filtring. Different terms should be seperated by \'@\' ")
	parser.add_argument('-c', '--Collapse', required=False, type=str, nargs='?', const='Sample_Collpase', help="Whether collapses multiple samples into one group, If yes: add '-c True'")
	parser.add_argument('-z', '--gzip', required=False, type=str, nargs='?', help="Whether gzip results? samples into one group, If yes: add '-z True'")
	parser.add_argument('-s', '--vcf_split_all_filter', required=False, type=str, nargs='?', help="Whether split all vcf filter? If yes: add '-s True'")
	parser.add_argument('-b', '--Bed', required=False, type=str, nargs='?', help="Bed File for SigProfilerMatrixGenerator")
	parser.add_argument('-C', '--Cluster', required=False, type=str, nargs='?', help="Whether perform SigProfilerClusters, If yes: add '-C True'")
	parser.add_argument('-P', '--Plotting', required=False, type=str, nargs='?', default="False", help="Whether Generate Plotting, If yes: add '-P True'")
	parser.add_argument('-S', '--SeqInfo', required=False, type=str, nargs='?', default="False", help="Whether Generate seqInfo, If yes: add '-S True'")

	
	if len(sys.argv) < 2:
		print(Help_String)
		sys.exit(1)
	args = parser.parse_args()

	return args.Input_Format, args.Input_Path, args.Project_ID, args.Output_Dir, args.Genome_Building, args.Data_Type, args.Filter, args.Collapse, args.gzip, args.vcf_split_all_filter, args.Bed, args.Cluster, args.Plotting, args.SeqInfo

####### 01-14 If input file is compressed?
def If_Compressed():
	### 000 Parse Options
	Input_Format,Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse,gzip,vcf_split_all_filter,Bed,Cluster,Plotting,SeqInfo = Parser()
	
	Input_Path_New_Name = Input_Path
	
	####### 获取输入文件所在的当前路径 ####### 
	Input_Dir = os.path.dirname(Input_Path)
	if Input_Dir == "":
		Input_Dir = "."
	Input_Dir_tmp = "%s/tmp" % (Input_Dir)
	
	if os.path.exists(Input_Dir_tmp):
		os.system("rm -rf %s" % Input_Dir_tmp)
	GenerateDir(Input_Dir_tmp)


	if re.search(r'zip$',Input_Path):
		String = "unzip %s -d %s" % (Input_Path,Input_Dir_tmp) 
		print(String)
		os.system(String)
		name = Input_Path.split("/")[-1].split(".zip")[0]
		Input_Path_New_Name = "%s/%s" % (Input_Dir_tmp,name)
		print(Input_Path_New_Name)


	if re.search(r'tar$',Input_Path):
		String = "tar -zxvf %s -C %s" % (Input_Path,Input_Dir_tmp) 
		print(String)
		os.system(String)
		name = Input_Path.split("/")[-1].split(".tar")[0]
		Input_Path_New_Name = "%s/%s" % (Input_Dir_tmp,name)
		print(Input_Path_New_Name)
		
	if re.search(r'gz$',Input_Path):
		if "tar.gz" not in Input_Path:
			name = Input_Path.split("/")[-1].split(".gz")[0]
			String = "gunzip -c %s > %s/%s" % (Input_Path,Input_Dir_tmp,name) 
			print(String)
			os.system(String)

			Input_Path_New_Name = "%s/%s" % (Input_Dir_tmp,name)
			print(Input_Path_New_Name)
		else:
			String = "tar -zxvf %s -C %s" % (Input_Path,Input_Dir_tmp) 
			print(String)
			os.system(String)
			name = Input_Path.split("/")[-1].split(".tar")[0]
			Input_Path_New_Name = "%s/%s" % (Input_Dir_tmp,name)
			print(Input_Path_New_Name)
	print(Input_Path_New_Name)
	return Input_Path_New_Name




####### 01-2 Parse Options
####### Choose Sub Function according to Options
def Parse_Options():
	Input_Format,Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse,gzip,vcf_split_all_filter,Bed,Cluster,Plotting,SeqInfo = Parser()


	### 000 if Output_Dir exists, delete it! #######
	if os.path.exists(Output_Dir):
		os.system("rm -rf %s" % Output_Dir)
	GenerateDir(Output_Dir) 
	
	
	###### Parse Options 001 Check Compressed.Generate New Input_Path, still named as Input_Path
	Input_Path = If_Compressed()

	###### Parse Options 002 Choose Transforming Function based on Input_Format


	if Input_Format == "csv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter == "True":
			if Filter == None:
				csv_Convert_Split(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
			else:
				print("Error 233, You have spcified '-s True', in this case the option -F is not supported!")
				sys.exit()
		else:
			###### Note Default of an option is None,not ''
			if Filter == None:
				csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
			else:
				csv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)


		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == "True":
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)


	elif Input_Format == "tsv":
		
		if vcf_split_all_filter == "True":
			if Filter == None:
				tsv_Convert_Split(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
			else:
				print("Error 233, You have spcified '-s True', in this case the option -F is not supported!")
				sys.exit()

		else:
			###### Note Default of an option is None,not ''
			if Filter == None:
				tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
			else:
				tsv_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)

		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == "True":
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)


	elif Input_Format == "catalog_tsv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports csv, tsv and vcf format")
			sys.exit()

		if Filter == None:
			pass
		else:
			String = "Error, \"%s\" format do not support option \"-F/--Filter\"" % (Input_Format)
			print(String)
			sys.exit()

		###### Note Default of an option is None,not ''
		if Collapse == "True":
			catalog_tsv_Convert_Collapse(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
		else:
			catalog_tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)



	elif Input_Format == "catalog_csv":
		###### if vcf_split_all_filter option is on
		if vcf_split_all_filter != None:
			print("Error: -s option only supports csv, tsv and vcf format")
			sys.exit()

		if Filter == None:
			pass
		else:
			String = "Error, \"%s\" format do not support option \"-F/--Filter\"" % (Input_Format)
			print(String)
			sys.exit()


		###### Note Default of an option is None,not ''
		if Collapse == "True":
			catalog_csv_Convert_Collapse(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
		else:
			catalog_csv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)



	elif Input_Format == "vcf":
		if vcf_split_all_filter == "True":
			if Filter == None:
				vcf_Multiple_Convert_Split_All_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
			else:
				print("Error 233, You have spcified '-s True', in this case the option -F is not supported!")
				sys.exit()

		else:
			if Filter == None:
				vcf_Multiple_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
			else:
				vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse)

		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == "True":
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)



	elif Input_Format == "maf":
		
		if vcf_split_all_filter != None:
			print("Error: -s option only supports csv, tsv and vcf format")
			sys.exit()

		if Filter == None:
			maf_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
		else:
			String = "Error, \"%s\" format do not support option \"-F/--Filter\"" % (Input_Format)
			print(String)
			sys.exit()


		###### Only if Collapse option is given, tsv_Convert_collpase will run
		if Collapse == "True":
			Convert_Collapse(Output_Dir,Collapse,Project_ID)
			String = "Finish Running Collapse Step"
			print(String)


	else:
		String = "Error in Format: Only 'vcf', 'csv', 'tsv', 'maf','catalog_csv', 'catalog_tsv' are validated formats!"
		print(String)
		sys.exit()

	###### Parse Options 003 Print Basic Output Statistic.
	Print_Statistic(Output_Dir)


	###### Parse Options 004 sigProfilerPlotting. Need to test Plotting and SeqInfo
	if Plotting not in ["True", "False"]:
		print("Error in Format, only \'True\' and \'False\' are supported for option: -p/--Plotting")
		sys.exit()
	if SeqInfo not in ["True", "False"]:
		print("Error in Format, only \'True\' and \'False\' are supported for option: -S/--SeqInfo_Para")
		sys.exit()
	else:
		print(Plotting, SeqInfo)
		sigProfilerPlotting(Input_Format,Output_Dir,Project_ID,Genome_Building,Bed,Plotting,SeqInfo)



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


	###### Parse Options 006 Run_SigProfilerClusters
	###### Only if Cluster option is given, Run_SigProfilerClusters will run
	if Input_Format == "vcf":
		if Cluster == "True":
			Run_SigProfilerClusters(Input_Path, Output_Dir, Project_ID, Genome_Building)
			String = "Finish Run_SigProfilerClusters Step\n"
			print(String)

	# if Input_Format != "vcf":
		# if Cluster == "True":
			# String = "Error: Option -c/Cluster only support VCF file as input \n"
			# print(String)
			# sys.exit()

####### 01-4 GenerateDir
def GenerateDir(Dir):
	if not os.path.exists(Dir):
		os.system("mkdir %s" % Dir)


####### 01-5 Delete_File
def Delete_File(File_Path):
	if os.path.exists(File_Path):
		os.system("rm %s" % (File_Path))


####### 01-6 Delete_Dir
def Delete_Dir(Dir_Path):
	if os.path.exists(Dir_Path):
		os.system("rm -rf %s" % (Dir_Path))


####### 31-1 maf_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def maf_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in MAF format. ******* ')

	####### 01-3-0 Output_Path:
	mSigPortal_Format_Tem_Path = "%s/%s_mSigPortal_Tem.txt" %  (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')


	####### 01-3-1 Check headers:
	Check_headers_Arr = ["Tumor_Sample_Barcode", "Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
	
	df = pd.read_table(Input_Path)
	
	
	for cc in Check_headers_Arr:
		if cc not in df.columns:
			print("Error 233: The column of %s can not be found from your MAF file!" % (cc))
			sys.exit()

	####### 01-3-2 Build Clean df:
	Clean_df = df[["Tumor_Sample_Barcode", "Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]]
	print(Clean_df)

	Clean_df.to_csv(mSigPortal_Format_Tem_Path, sep="\t", header=True, quoting=None, index=None)


	####### 01-3-3 Parse File 
	Input_File = open(mSigPortal_Format_Tem_Path)
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split(",")
		String_File += line
	Input_File.close()

	####### 01-3-4 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'Tumor_Sample_Barcode',f):
			pass
		else:
			ss = f.strip().split("\t")
			if len(ss) > 6:
				Sample_ID = ss[0]
				Chr = ss[1]
				Start = ss[2]
				End = ss[3]
				REF = ss[4]
				ALT_1 = ss[5]
				ALT_2 = ss[6]
				ALT_Final = ""
				if ALT_1 == REF:
					ALT_Final = ALT_2
				else:
					ALT_Final = ALT_1

				if "," in REF:
					pass
				else:
					if "-" in REF or "-" in ALT_1 or "-" in ALT_2:
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT_Final)
						mSigPortal_Format_INDEL_File.write(Output_String)
					elif len(REF) != len(ALT_2):
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT_Final)
						mSigPortal_Format_INDEL_File.write(Output_String)
					else:
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT_Final)
						mSigPortal_Format_SNV_File.write(Output_String)


	mSigPortal_Format_SNV_File.close()
	mSigPortal_Format_INDEL_File.close()


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
	Header = "SAMPLE,CHROM,START,END,REF,ALT,FILTER"
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split(",")
		if len(ss) == 7:
			String_File += line
		else:
			print(ss)
			print("Error 233: The CSV option requires each line in input file should include 7 Item: SAMPLE,CHROM,START,END,REF,ALT,FILTER.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"SAMPLE,CHROM,START,END,REF,ALT,FILTER\" is required!")
		sys.exit()
		
	####### 01-3-2 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'SAMPLE',f):
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
				
				if "," in ALT:
					pass
				else:
					if len(REF) == len(ALT):
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)
	Input_File.close()
	mSigPortal_Format_SNV_File.close()
	mSigPortal_Format_INDEL_File.close()


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
	Header = "SAMPLE,CHROM,START,END,REF,ALT,FILTER"
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
		print("Error 233: A header line: \"SAMPLE,CHROM,START,END,REF,ALT,FILTER\" is required!")
		sys.exit()

	####### 01-6-2 Parse Filter:
	option_Filter_Arr = Filter.split("@")

	####### 01-6-3 Generate Result
	ff = String_File.split("\n")
	
	####### 如果filter中有“-”用csv_Filter_tmp_arr_SNV和csv_Filter_tmp_arr_INDEL来存储这句话，并去重复
	csv_Filter_tmp_arr_SNV = []
	csv_Filter_tmp_arr_INDEL = []
	
	for f in ff:
		if re.match(r'SAMPLE',f):
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
				
				Filter_arr = ss[6].strip().split(";")

				###### 首先，不管有木有filter，如果filter中有“-”，都要先输出不filter时的样本名称
				if "-" in option_Filter_Arr:
					
					
					if len(REF) == len(ALT):
						Output_String = "%s	All_Samples	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						csv_Filter_tmp_arr_SNV.append(Output_String)
						#mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	All_Samples	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						csv_Filter_tmp_arr_INDEL.append(Output_String)
						#mSigPortal_Format_INDEL_File.write(Output_String)


				###### 其次，进行filter
				if "," in ALT:
					pass
				else:
					for Filter in option_Filter_Arr:
						if Filter in Filter_arr:
							Sample_ID = "%s@%s" % (Sample_ID,Filter)
							if len(REF) == len(ALT):
								Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
								mSigPortal_Format_SNV_File.write(Output_String)
							else:
								Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
								mSigPortal_Format_INDEL_File.write(Output_String)

	for i in set(csv_Filter_tmp_arr_SNV):
		mSigPortal_Format_SNV_File.write(i)
	for i in set(csv_Filter_tmp_arr_INDEL):
		mSigPortal_Format_INDEL_File.write(i)
	
	mSigPortal_Format_SNV_File.close()
	mSigPortal_Format_INDEL_File.close()
	
	Input_File.close()



####### 01-16 csv_Convert_Split(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
def csv_Convert_Split(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in TSV format, with filtration: -s True ******* ')

	####### 01-6-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-6-1 Parse File 
	Input_File = open(Input_Path)
	Header = "SAMPLE,CHROM,START,END,REF,ALT,FILTER"
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
		print("Error 233: A header line: \"SAMPLE,CHROM,START,END,REF,ALT,FILTER\" is required!")
		sys.exit()

	####### 01-6-2 Parse Filter:
	#option_Filter_Arr = Filter.split("@")

	####### 01-6-3 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'SAMPLE',f):
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
				Filter_arr = ss[6].strip().split(";")

				if "," in ALT:
					pass
				else:
					
					###### 首先，不管有木有filter，都要先输出不filter时的样本名称
					if len(REF) == len(ALT):
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)

					###### 其次，进行Filter
					for Filter in Filter_arr:
						Sample_ID_Final = "%s@%s" % (Sample_ID,Filter)
						if len(REF) == len(ALT):
							Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
							mSigPortal_Format_SNV_File.write(Output_String)
						else:
							Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
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
	Header = "SAMPLE	CHROM	START	END	REF	ALT	FILTER"
	#Header = ""
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split("	")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The TSV format requires each line in input file should include 7 Item: SAMPLE	CHROM	START	END	REF	ALT	FILTER.\nHowever, the following line contains %d items!" % (len(ss)))
			print(ss)
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"SAMPLE	CHROM	START	END	REF	ALT	FILTER\" is required!")
		sys.exit()
	#print(String_File)
	####### 01-3-2 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'SAMPLE',f):
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
				if "," in ALT:
					pass
				else:
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
	Header = "SAMPLE	CHROM	START	END	REF	ALT	FILTER"
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split("	")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The TSV optformat requires each line in input file should include 7 Item: Sample_ID	Chrom	Start	End	Ref	Alt	Filter.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"SAMPLE	CHROM	START	END	REF	ALT	FILTER\" is required!")
		sys.exit()

	####### 01-6-2 Parse Filter:
	option_Filter_Arr = Filter.split("@")

	####### 01-6-3 Generate Result
	ff = String_File.split("\n")
	
	####### 如果filter中有“-”用csv_Filter_tmp_arr_SNV和csv_Filter_tmp_arr_INDEL来存储这句话，并去重复
	csv_Filter_tmp_arr_SNV = []
	csv_Filter_tmp_arr_INDEL = []

	for f in ff:
		if re.match(r'SAMPLE',f):
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
				Filter_arr = ss[6].strip().split(";")
				
				###### 首先，不管有木有filter，如果filter中有“-”，都要先输出不filter时的样本名称
				if "-" in option_Filter_Arr:
					if len(REF) == len(ALT):
						Output_String = "%s	All_Samples	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						csv_Filter_tmp_arr_SNV.append(Output_String)
						#mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	All_Samples	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						csv_Filter_tmp_arr_INDEL.append(Output_String)
						#mSigPortal_Format_INDEL_File.write(Output_String)

				###### 其次，filter
				if "," in ALT:
					pass
				else:
					for Filter in option_Filter_Arr:
						if Filter in Filter_arr:
							Sample_ID = "%s@%s" % (Sample_ID,Filter)
							if len(REF) == len(ALT):
								Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
								mSigPortal_Format_SNV_File.write(Output_String)
							else:
								Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
								mSigPortal_Format_INDEL_File.write(Output_String)

	for i in set(csv_Filter_tmp_arr_SNV):
		mSigPortal_Format_SNV_File.write(i)
	for i in set(csv_Filter_tmp_arr_INDEL):
		mSigPortal_Format_INDEL_File.write(i)
	
	mSigPortal_Format_SNV_File.close()
	mSigPortal_Format_INDEL_File.close()

	Input_File.close()



####### 01-18 tsv_Convert_Split(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse)
def tsv_Convert_Split(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in TSV format, with filtration: -s True ******* ')

	####### 01-6-0 Output_Path:
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)

	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path,'w')
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path,'w')

	####### 01-6-1 Parse File 
	Input_File = open(Input_Path)
	Header = "SAMPLE	CHROM	START	END	REF	ALT	FILTER"
	String_File = ""
	Count = 1
	for line in Input_File:
		ss = line.split("	")
		if len(ss) == 7:
			String_File += line
		else:
			print("Error 233: The TSV optformat requires each line in input file should include 7 Item: Sample_ID	Chrom	Start	End	Ref	Alt	Filter.\nHowever, the %d line contains %d items!" % (Count,len(ss)))
			sys.exit()
			Count += 1
	if Header not in String_File:
		print("Error 233: A header line: \"SAMPLE	CHROM	START	END	REF	ALT	FILTER\" is required!")
		sys.exit()

	####### 01-6-3 Generate Result
	ff = String_File.split("\n")
	for f in ff:
		if re.match(r'SAMPLE',f):
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
				Filter_arr = ss[6].strip().split(";")

				if "," in ALT:
					pass
				else:
					###### 首先，不管有木有filter，都要先输出不filter时的样本名称
					if len(REF) == len(ALT):
						Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
					else:
						Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_INDEL_File.write(Output_String)


					for Filter in Filter_arr:
						Sample_ID_Final = "%s@%s" % (Sample_ID,Filter)
						#print(Sample_ID)
						if len(REF) == len(ALT):
							Output_String = "%s	%s	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
							#print(Output_String)
							mSigPortal_Format_SNV_File.write(Output_String)
						else:
							Output_String = "%s	%s	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Sample_ID_Final,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
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
			
			
			####### Sometimes, if there is a "," in ALT, then pass
			if "," in ALT:
				pass
			else:
				
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
			Filter_arr = ss[6].split(";")
			
			####### Sometimes, if there is a "," in ALT, then pass
			if "," in ALT:
				pass
			else:
			
				####### 对每个样本 Sample_ID_arr中的成员进行便利，判断是否有基因型信息
				for i in range(0,len(Sample_ID_arr)):
					Genotype_String = ss[9+i]
					gg = Genotype_String.split(":")[0]
					if re.search(r'1',gg):
						#print i,Sample_ID_arr[i]
						Sample_ID = Sample_ID_arr[i]
						
						
						####### 在 Filter之前先打印所有样本未filter的信息 #######
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
						
						
						
						####### 对每个sample 进行filter #######
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




####### 01-11 catalog_tsv_Convert_Collapse(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def catalog_tsv_Convert_Collapse(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
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
			mSigPortal_Format_catalog_File.write("%s	All_Samples\n" % (Header))

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


####### 01-11-2 catalog_tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
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
	Input_File.close()

	####### 01-3-2 Generate Result
	Input_File = open(Input_Path)
	for New_line in Input_File:
		if re.match(r'\n',line):
			pass
		else:
			mSigPortal_Format_catalog_File.write(New_line)

	Input_File.close()
	mSigPortal_Format_catalog_File.close()


####### 01-12-2 catalog_tsv_Convert(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
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

	for line in Input_File:
		if re.match(r'\n',line):
			pass
		else:
			String_File += line
		
	if Header not in String_File:
		print("Error 233: A header line like \"MutationType,Sample1,Sample2,Sample3...\" is required!")
		sys.exit()
	Input_File.close()

	####### 01-3-2 Generate Result
	Input_File = open(Input_Path)
	for New_line in Input_File:
		if re.match(r'\n',line):
			pass
		else:
			mSigPortal_Format_catalog_File.write(New_line.replace(",","	"))

	Input_File.close()
	mSigPortal_Format_catalog_File.close()


####### 01-12 catalog_csv_Convert_Collapse(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type)
def catalog_csv_Convert_Collapse(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Collapse):
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
			mSigPortal_Format_catalog_File.write("%s	All_Samples\n" % Header.replace(",","\t"))
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


####### 01-14 vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter)
def vcf_Multiple_Convert_Filter(Input_Path,Project_ID,Output_Dir,Genome_Building,Data_Type,Filter,Collapse):
	GenerateDir(Output_Dir)
	print('******* Your Input File is in VCF format with Multple sample. Your Filter term is: %s ******* ' % (Filter))
	
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

			####### Sometimes, if there is a "," in ALT, then pass
			if "," in ALT:
				pass
			else:

				### Pay attention to option "-"
				if "-" in option_Filter_Arr:
					if len(REF) == len(ALT):
						End = ss[1]
						Output_String = "%s	All_Samples	%s	%s	SNV	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
						mSigPortal_Format_SNV_File.write(Output_String)
					else:
						if len(REF) == 1:
							Start = ss[1]
							End = ss[1]
						else:
							Start = ss[1]
							End = int(Start) + len(REF) - 1
						Output_String = "%s	All_Samples	%s	%s	INDEL	%s	%s	%s	%s	%s	SOMATIC\n" % (Project_ID,Data_Type,Genome_Building,Chr,Start,End,REF,ALT)
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
	
	count = len(set(SNV_Sample_Count))
	String_Collpase = "There are %d samples in the original Output File:%s" % (count,mSigPortal_Format_SNV_Path)
	
	####### Only the Input file with more than 1 sample can support “Collapse Option” !
	####### 还有一个问题，collpase之后，会出现，重复字段。所以这里面要设置一个 collapse_String_arr
	collapse_String_arr = []

	####### 01-12-3 Parse Input_Path: mSigPortal_Format_SNV_File
	if count > 1:
		mSigPortal_Format_SNV_Collapse_File = open(mSigPortal_Format_SNV_Collapse_Path,'w')
		mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path)
		for line in mSigPortal_Format_SNV_File:
			ss = line.strip().split("	")
			
			###### 要判断 sample后面是否有filter
			Sample_and_Filter_arr = ss[1].split("@")
			Sample_Name = ss[1].split("@")[0]
			
			Collapse_Sample_Name = "All_Samples"
			if len(Sample_and_Filter_arr) == 2:
				####### 这就说明是带了filter的比如 SC00101@PASS
				Collpase_Filter = Sample_and_Filter_arr[1]
				Collapse_Sample_Name = "All_Samples@%s" % Collpase_Filter
				
			String_1 = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			String_2 = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],Collapse_Sample_Name,ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			mSigPortal_Format_SNV_Collapse_File.write(String_1)
			
			collapse_String_arr.append(String_2)

		for line_2 in set(collapse_String_arr):
			mSigPortal_Format_SNV_Collapse_File.write(line_2)
		
		mSigPortal_Format_SNV_File.close()
		mSigPortal_Format_SNV_Collapse_File.close()

		####### 01-6-0 RM mSigPortal_Format_SNV_Path
		cmd_1 = "rm %s" % (mSigPortal_Format_SNV_Path)
		os.system(cmd_1)



	####### 01-12-4 Count the Sample in Input_Path: mSigPortal_Format_INDEL_File
	INDEL_Sample_Count = []

	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path)
	for line in mSigPortal_Format_INDEL_File:
		ss = line.strip().split("	")
		Sample_Name = ss[1].split("@")[0]
		INDEL_Sample_Count.append(Sample_Name)
	count = len(set(INDEL_Sample_Count))
	String_Collpase = "There are %d samples in the original Output File:%s" % (count,mSigPortal_Format_SNV_Path)
	mSigPortal_Format_INDEL_File.close()

	####### Only the Input file with more than 1 sample can support “Collapse Option” !
	####### 还有一个问题，collpase之后，会出现，重复字段。所以这里面要设置一个 collapse_String_arr
	collapse_String_arr = []


	####### 01-12-5 Parse Input_Path: mSigPortal_Format_INDEL_File
	if count > 1:
		mSigPortal_Format_INDEL_Collapse_File = open(mSigPortal_Format_INDEL_Collapse_Path,'w')
		mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path)
		for line in mSigPortal_Format_INDEL_File:
			ss = line.strip().split("	")
			
			###### 要判断 sample后面是否有filter
			Sample_and_Filter_arr = ss[1].split("@")
			Sample_Name = ss[1].split("@")[0]
			
			Collapse_Sample_Name = "All_Samples"
			if len(Sample_and_Filter_arr) == 2:
				####### 这就说明是带了filter的比如 SC00101@PASS
				Collpase_Filter = Sample_and_Filter_arr[1]
				Collapse_Sample_Name = "All_Samples@%s" % Collpase_Filter

			Sample_Name = ss[1].split("@")[0]
			String_1 = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],ss[1],ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			String_2 = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (ss[0],Collapse_Sample_Name,ss[2],ss[3],ss[4],ss[5],ss[6],ss[7],ss[8],ss[9],ss[10])
			mSigPortal_Format_INDEL_Collapse_File.write(String_1)

			collapse_String_arr.append(String_2)

		for line_2 in set(collapse_String_arr):
			mSigPortal_Format_INDEL_Collapse_File.write(line_2)

		mSigPortal_Format_INDEL_Collapse_File.close()
		mSigPortal_Format_INDEL_File.close()

		####### 01-6-0 RM  mSigPortal_Format_INDEL_Path
		cmd_2 = "rm %s" % (mSigPortal_Format_INDEL_Path)
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
#def sigProfilerPlotting(Input_Format, Output_Dir, Project_ID, Genome_Building, Bed):
def sigProfilerPlotting(Input_Format, Output_Dir, Project_ID, Genome_Building,Bed,Plotting,SeqInfo):
	####### 01-17-0 If there is 0 line in the intput file, delete it.
	mSigPortal_Format_SNV_Path = "%s/%s_mSigPortal_SNV.txt" %  (Output_Dir,Project_ID)
	n_SNV = 0
	mSigPortal_Format_SNV_File = open(mSigPortal_Format_SNV_Path)
	for line in mSigPortal_Format_SNV_File:
		if re.match(r'\n', line):
			pass
		else:
			n_SNV += 1
	mSigPortal_Format_SNV_File.close()
	if n_SNV == 0:
		Note = "\nThere is 0 line in input file: %s, will delete it!\n" % mSigPortal_Format_SNV_Path
		print(Note)
		os.system("rm %s" % mSigPortal_Format_SNV_Path)


	mSigPortal_Format_INDEL_Path = "%s/%s_mSigPortal_INDEL.txt" % (Output_Dir,Project_ID)
	n_INDEL = 0
	mSigPortal_Format_INDEL_File = open(mSigPortal_Format_INDEL_Path)
	for line in mSigPortal_Format_INDEL_File:
		if re.match(r'\n', line):
			pass
		else:
			n_INDEL += 1
	mSigPortal_Format_INDEL_File.close()

	if n_INDEL == 0:
		Note = "\nThere is 0 line in input file: %s, will delete it!\n" % mSigPortal_Format_INDEL_Path
		print(Note)
		os.system("rm %s" % mSigPortal_Format_INDEL_Path)



	####### 01-17-1 Which format is the input file
	Input_Format_arr_1 = ['vcf', 'csv', 'tsv', 'maf']
	Input_Format_arr_2 = ['catalog_csv', 'catalog_tsv']
	
	SBS_Arr = [6, 24, 96, 384, 1536, 6144]
	ID_Arr = [28, 83, 415, 8268]
	DBS_Arr = [78, 186, 312, 1248, 2976]
	CNV_Arr = [48]
	R32_Arr = [32]

	Input_Path_arr = ["mSigPortal_SNV.txt", "mSigPortal_INDEL.txt", "mSigPortal_SNV_Collapse.txt", "mSigPortal_INDEL_Collapse.txt","mSigPortal_catalog_csv.txt", "mSigPortal_catalog_tsv.txt"]


	####### Which format is the input file
	if Plotting == "False":
		Plotting = False
	if Plotting == "True":
		plotting = True
	if SeqInfo == "False":
		SeqInfo = False
	if SeqInfo == "True":
		SeqInfo = True


	####### Which format is the input file
	if Input_Format in Input_Format_arr_1:
		matrices = matGen.SigProfilerMatrixGeneratorFunc(Project_ID, Genome_Building, Output_Dir, exome=False, bed_file=Bed, chrom_based=False, plot=Plotting, tsb_stat=False, seqInfo=SeqInfo)
		
		####### Generate Summary File
		summary_Path = "%s/svg_files_list.txt" % (Output_Dir)
		summary_File = open(summary_Path, 'w')
		Header = "Sample_Name,Profile_Type,Matrix_Size,Filter,Path\n"
		summary_File.write(Header)
		SVG_Ouput_Dir = "%s/output/plots/svg" % (Output_Dir)
		#print(SVG_Ouput_Dir)
		SVG_New_Output_Dir = "%s/output/svg" % (Output_Dir)
		if os.path.exists(SVG_New_Output_Dir):
			os.system("mv %s %s" % (SVG_Ouput_Dir, SVG_New_Output_Dir))

		####### Generate Download File and Matrix_List_File #######
		Matrix_List_Path = "%s/matrix_files_list.txt" % (Output_Dir)
		
		DBS_Path = "%s/output/DBS" % (Output_Dir)
		ID_Path = "%s/output/ID" % (Output_Dir)
		PDF_Path = "%s/output/plots" % (Output_Dir)
		SBS_Path = "%s/output/SBS" % (Output_Dir)
		Matrix_Path = "%s/output/vcf_files" % (Output_Dir)
		Matrix_ZiP_Path = "%s/output/vcf_files_zip" % (Output_Dir) # tar directory without whole structure // 04-01-2021
		GenerateDir(Matrix_ZiP_Path)


		Matrix_List_File = open(Matrix_List_Path, 'w')
		Header = "Profile_Type,Matrix_Size,Path\n"
		Matrix_List_File.write(Header)

		if os.path.exists(DBS_Path):
			os.system("zip -jr %s.zip %s" % (DBS_Path,DBS_Path))

			Catalog = "DBS"
			for ii in os.listdir(DBS_Path):
				Type = ii.split(".")[1].split("DBS")[1]
				Path = "%s/%s" % (DBS_Path,ii)
				Path = os.path.abspath(Path)
				Final_String = "%s,%s,%s\n" % (Catalog,Type,Path)
				Matrix_List_File.write(Final_String)


		if os.path.exists(ID_Path):
			#os.system("tar -zcvf %s.tar.gz %s" % (ID_Path,ID_Path))
			os.system("zip -jr %s.zip %s" % (ID_Path,ID_Path))

			Catalog = "ID"
			for ii in os.listdir(ID_Path):
				Type = ii.split(".")[1].split("ID")[1]
				Path = "%s/%s" % (ID_Path,ii)
				Path = os.path.abspath(Path)
				Final_String = "%s,%s,%s\n" % (Catalog,Type,Path)
				Matrix_List_File.write(Final_String)


		if os.path.exists(PDF_Path):
			#os.system("tar -zcvf %s.tar.gz %s" % (PDF_Path,PDF_Path))
			os.system("zip -jr %s.zip %s" % (PDF_Path,PDF_Path))


		if os.path.exists(SBS_Path):
			os.system("zip -jr %s.zip %s" % (SBS_Path,SBS_Path))

			Catalog = "SBS"
			for ii in os.listdir(SBS_Path):
				Type = ii.split(".")[1].split("SBS")[1]
				Path = "%s/%s" % (SBS_Path,ii)
				Path = os.path.abspath(Path)
				Final_String = "%s,%s,%s\n" % (Catalog,Type,Path)
				Matrix_List_File.write(Final_String)

		if os.path.exists(Matrix_Path):
			#os.system("tar -zcvf %s.tar.gz %s" % (Matrix_Path,Matrix_Path))
			for iii in os.listdir(Matrix_Path):
				Sub_iii = "%s/%s" % (Matrix_Path,iii)
				zip_Cmd = "zip -jr %s/%s_txt.zip %s" % (Matrix_ZiP_Path,iii,Sub_iii)
				os.system(zip_Cmd)

		Matrix_List_File.close()


		####### zip #######
		if os.path.exists(SVG_New_Output_Dir):
			for svg in os.listdir(SVG_New_Output_Dir):
				if "_plots_" in svg:
					#print(svg)
					Type = svg.split("_plots_")[0]
					Profile_Type = Type.split("_")[0]
					Matrix = "%s" % (Type.split("_")[1])

					Tag = "NA"
					sample_Name = ""
					sample_Name_Tag = svg.split("%s_" % Project_ID)[1].strip(".svg")
					if "@" in sample_Name_Tag:
						Tag = sample_Name_Tag.split("@")[1]
						sample_Name = sample_Name_Tag.split("@")[0]
					else:
						sample_Name = sample_Name_Tag
					if sample_Name == "filter":
						pass
					else:
						svg_Location = "%s/%s" % (SVG_New_Output_Dir,svg)
						abs_path = os.path.abspath(svg_Location)
						String = "%s,%s,%s,%s,%s\n" % (sample_Name,Profile_Type,Matrix,Tag,abs_path)
						summary_File.write(String)
		summary_File.close()

	elif Input_Format in Input_Format_arr_2:
		####### Generate Download File and Matrix_List_File #######
		Matrix_List_Path = "%s/matrix_files_list.txt" % (Output_Dir)
		Matrix_List_File = open(Matrix_List_Path, 'w')
		Header = "Profile_Type,Matrix_Size,Path\n"
		Matrix_List_File.write(Header)

		for matrix_name in os.listdir(Output_Dir):
			for i in Input_Path_arr:
				if re.search(i, matrix_name):
					count = 0
					matrix_path = "%s/%s" % (Output_Dir, matrix_name)
					matrix_File = open(matrix_path)
					for line in matrix_File:
						ss = line.split("	")
						if len(ss) > 1:
							count += 1
					Type = count-1
					matrix_File.close()
					
					#print(Type)
					###### Plotting the Matrix Based on Type
					Final_output_Dir = "%s/" % (Output_Dir)
					Final_Type = "%d" % Type

					FF_Dir = "%s/output" % Output_Dir
					GenerateDir(FF_Dir)


					#print(Type)
					if Type in SBS_Arr:
						sigPlt.plotSBS(matrix_path, Final_output_Dir, Project_ID, Final_Type, percentage=True)
						FF_temp_Dir = "%s/SBS" % FF_Dir
						FF_temp_cmd_1 = "mkdir %s" % FF_temp_Dir
						FF_temp_cmd_2 = "cp %s %s/%s.SBS%s.all" % (matrix_path, FF_temp_Dir, Project_ID, Final_Type)

						os.system(FF_temp_cmd_1)
						os.system(FF_temp_cmd_2)

						Path = "%s/%s.SBS%s.all" % (FF_temp_Dir, Project_ID, Final_Type)
						Path = os.path.abspath(Path)
						Catalog = "SBS"
						Final_String = "%s,%s,%s\n" % (Catalog, Type, Path)
						Matrix_List_File.write(Final_String)
						#os.system("tar -zcvf %s.tar.gz %s" % (FF_temp_Dir, FF_temp_Dir))
						os.system("zip -jr %s.zip %s" % (FF_temp_Dir, FF_temp_Dir))


					elif Type in DBS_Arr:
						sigPlt.plotDBS(matrix_path, Final_output_Dir, Project_ID, Final_Type, percentage=True)
						FF_temp_Dir = "%s/DBS" % FF_Dir
						FF_temp_cmd_1 = "mkdir %s" % FF_temp_Dir
						FF_temp_cmd_2 = "cp %s %s/%s.DBS%s.all" % (matrix_path, FF_temp_Dir, Project_ID, Final_Type)

						os.system(FF_temp_cmd_1)
						os.system(FF_temp_cmd_2)

						Path = "%s/%s.DBS%s.all" % (FF_temp_Dir, Project_ID, Final_Type)
						Path = os.path.abspath(Path)
						Catalog = "DBS"
						Final_String = "%s,%s,%s\n" % (Catalog, Type, Path)
						Matrix_List_File.write(Final_String)

						#os.system("tar -zcvf %s.tar.gz %s" % (FF_temp_Dir, FF_temp_Dir))
						os.system("zip -jr %s.zip %s" % (FF_temp_Dir, FF_temp_Dir))


					elif Type in ID_Arr:
						sigPlt.plotID(matrix_path, Final_output_Dir, Project_ID, Final_Type, percentage=True)

						FF_temp_Dir = "%s/ID" % FF_Dir
						FF_temp_cmd_1 = "mkdir %s" % FF_temp_Dir
						FF_temp_cmd_2 = "cp %s %s/%s.ID%s.all" % (matrix_path, FF_temp_Dir, Project_ID, Final_Type)

						os.system(FF_temp_cmd_1)
						os.system(FF_temp_cmd_2)

						Path = "%s/%s.ID%s.all" % (FF_temp_Dir, Project_ID, Final_Type)
						Path = os.path.abspath(Path)
						Catalog = "ID"
						Final_String = "%s,%s,%s\n" % (Catalog, Type, Path)
						Matrix_List_File.write(Final_String)

						#os.system("tar -zcvf %s.tar.gz %s" % (FF_temp_Dir, FF_temp_Dir))
						os.system("zip -jr %s.zip %s" % (FF_temp_Dir, FF_temp_Dir))


					elif Type in CNV_Arr:
						sigPlt.plotCNV(matrix_path, Final_output_Dir, Project_ID, percentage=True, aggregate=False)

						#sigPlt.plotCNV(matrix_path, output_path, project, "pdf", percentage=True, aggregate=False)

						FF_temp_Dir = "%s/CNV" % FF_Dir
						FF_temp_cmd_1 = "mkdir %s" % FF_temp_Dir
						FF_temp_cmd_2 = "cp %s %s/%s.CNV%s.all" % (matrix_path, FF_temp_Dir, Project_ID, Final_Type)

						os.system(FF_temp_cmd_1)
						os.system(FF_temp_cmd_2)

						Path = "%s/%s.CNV%s.all" % (FF_temp_Dir, Project_ID, Final_Type)
						Path = os.path.abspath(Path)
						Catalog = "CNV"
						Final_String = "%s,%s,%s\n" % (Catalog, Type, Path)
						Matrix_List_File.write(Final_String)

						#os.system("tar -zcvf %s.tar.gz %s" % (FF_temp_Dir, FF_temp_Dir))
						os.system("zip -jr %s.zip %s" % (FF_temp_Dir, FF_temp_Dir))

					elif Type in R32_Arr:
						sigPlt.plotSV(matrix_path, Final_output_Dir, Project_ID, percentage=True, aggregate=False)

						#sigPlt.plotCNV(matrix_path, output_path, project, "pdf", percentage=True, aggregate=False)

						FF_temp_Dir = "%s/RS" % FF_Dir
						FF_temp_cmd_1 = "mkdir %s" % FF_temp_Dir
						FF_temp_cmd_2 = "cp %s %s/%s.RS%s.all" % (matrix_path, FF_temp_Dir, Project_ID, Final_Type)

						os.system(FF_temp_cmd_1)
						os.system(FF_temp_cmd_2)

						Path = "%s/%s.RS%s.all" % (FF_temp_Dir, Project_ID, Final_Type)
						Path = os.path.abspath(Path)
						Catalog = "RS"
						Final_String = "%s,%s,%s\n" % (Catalog, Type, Path)
						Matrix_List_File.write(Final_String)

						#os.system("tar -zcvf %s.tar.gz %s" % (FF_temp_Dir, FF_temp_Dir))
						os.system("zip -jr %s.zip %s" % (FF_temp_Dir, FF_temp_Dir))



					else:
						print("Error 233: Your input type in the file is not supported yet!" )
						sys.exit()

					###### Try to Generate same output Structure as non catalog format
					FF_Dir_pdf = "%s/plots" % (FF_Dir)
					GenerateDir(FF_Dir_pdf)
					FF_command_1 = "mv %s*pdf %s/" % (Final_output_Dir, FF_Dir_pdf)
					os.system(FF_command_1)
					FF_Dir_svg = "%s/svg" % (FF_Dir)
					
					if os.path.exists(FF_Dir_svg):
						FF_command_2 = "mv %s/svg %s/" % (Output_Dir, FF_Dir)
						print(FF_command_2)
						os.system(FF_command_2)
					#os.system("tar -zcvf %s.tar.gz %s" % (FF_Dir_pdf, FF_Dir_pdf))
					#os.system("tar -zcvf %s.tar.gz %s" % (FF_Dir_svg, FF_Dir_svg))

					os.system("zip -jr %s.zip %s" % (FF_Dir_pdf, FF_Dir_pdf))
					os.system("zip -jr %s.zip %s" % (FF_Dir_svg, FF_Dir_svg))


		Matrix_List_File.close()
		print("Finisheh !!!!!!!!!!!!!!!!!!!!!!!!!!")
		print(Final_Type)




		####### Generate Summary File
		summary_Path = "%s/svg_files_list.txt" % (Output_Dir)
		summary_File = open(summary_Path,'w')
		Header = "Sample_Name,Profile_Type,Matrix_Size,Filter,Path\n"
		summary_File.write(Header)
		SVG_Ouput_Dir = "%s/output/svg" % (Output_Dir)

		for svg in os.listdir(SVG_Ouput_Dir):
			if "_plots_" in svg:
				#print(svg)
				Type = svg.split("_plots_")[0]
				Profile_Type = Type.split("_")[0]
				Matrix = "%s" % (Type.split("_")[1])
				
				Tag = "NA"
				sample_Name = ""
				sample_Name_Tag = svg.split("%s_" % Project_ID)[1].strip(".svg")
				if "@" in sample_Name_Tag:
					Tag = sample_Name_Tag.split("@")[1]
					sample_Name = sample_Name_Tag.split("@")[0]
				else:
					sample_Name = sample_Name_Tag
				if sample_Name == "filter":
					pass
				else:
					svg_Location = "%s/%s" % (SVG_Ouput_Dir,svg)
					String = "%s,%s,%s,%s,%s\n" % (sample_Name,Profile_Type,Matrix,Tag,svg_Location)
					summary_File.write(String)
		summary_File.close()


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


	Input_Path_arr_No_Catlog = ["mSigPortal_SNV.txt","mSigPortal_INDEL.txt","mSigPortal_SNV_Collapse.txt","mSigPortal_INDEL_Collapse.txt"]
	for OutputFile in os.listdir(Output_Dir):
		for i in Input_Path_arr:
			if re.search(r'%s$' % i,OutputFile):
				Output_Path = "%s/%s" % (Output_Dir,OutputFile)
				output_Dict = {}
				Output_File = open(Output_Path)
				for line in Output_File:
					ss = line.split("	")
					Sample_ID = ss[1]
					output_Dict.setdefault(Sample_ID,[]).append(line)
				Output_File.close()
				
				String = "The counts for each samples@filtration in OutputFile: %s is as following:" % (Output_Path)
				print(String)
				for key in output_Dict:
					Value = len(output_Dict[key])
					print("%s	%s" % (key,Value))


####### 01-19 Run_SigProfilerClusters(Input_Path, Output_Dir)
####### This function only support VCF format

def Run_SigProfilerClusters(Input_Path, Output_Dir, Project_ID, Genome_Building):
	# need bcftools to be installed
	from SigProfilerSimulator import SigProfilerSimulator as sigSim
	from SigProfilerClusters import SigProfilerClusters as hp
	import subprocess
	
	print("\n\n Run_SigProfilerClusters \n\n")

	Run_SigProfilerClusters_Output = "%s/Cluster/" % (Output_Dir)
	GenerateDir(Run_SigProfilerClusters_Output)

	### 001 First cp and split Input vcf file into Dir
	cmd = "bcftools query -l %s" % (Input_Path)
	result = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	m_stdout, m_stderr = result.communicate()
	sample_ID_str = m_stdout.decode("utf8")
	sample_ID_arr = sample_ID_str.strip().split("\n")
	
	for sample in sample_ID_arr:
		
		sub_cmd_Str = "bcftools view -c1 -s %s -o %s/%s.vcf %s" % (sample, Run_SigProfilerClusters_Output, sample, Input_Path)
		os.system(sub_cmd_Str)
	
	#cmd_001_cp_VCF = "cp %s %s" % (Input_Path, Run_SigProfilerClusters_Output)
	#os.system(cmd_001_cp_VCF)

	### 002 SigProfilerSimulator
	sigSim.SigProfilerSimulator(Project_ID, Run_SigProfilerClusters_Output, Genome_Building, contexts = ['96'], chrom_based=True, simulations=100)

	### 003 SigProfilerClusters
	hp.analysis(Project_ID, Genome_Building, "96", ["96"], Run_SigProfilerClusters_Output, analysis="all", sortSims=True, subClassify=True, correction=True, calculateIMD=True, max_cpu=6, standardVC=True, TCGA=False, sanger=False)




	### 004 Integrate SigProfilerClusters results of text

	print("\n####### Arrange Final Results: Clusters File\n")
	
	SigProfilerClusters_Result = "%sResult" % (Run_SigProfilerClusters_Output)

	GenerateDir(SigProfilerClusters_Result)
	
	

	c1a_Path = "%soutput/vcf_files_corrected/%s_clustered/subclasses/class1a/%s_clustered_class1a.txt" % (Run_SigProfilerClusters_Output,Project_ID, Project_ID)
	c1b_Path = "%soutput/vcf_files_corrected/%s_clustered/subclasses/class1b/%s_clustered_class1b.txt" % (Run_SigProfilerClusters_Output,Project_ID, Project_ID)
	c1c_Path = "%soutput/vcf_files_corrected/%s_clustered/subclasses/class1c/%s_clustered_class1c.txt" % (Run_SigProfilerClusters_Output,Project_ID, Project_ID)
	c2_Path = "%soutput/vcf_files_corrected/%s_clustered/subclasses/class2/%s_clustered_class2.txt" % (Run_SigProfilerClusters_Output,Project_ID, Project_ID)
	c3_Path = "%soutput/vcf_files_corrected/%s_clustered/subclasses/class3/%s_clustered_class3.txt" % (Run_SigProfilerClusters_Output,Project_ID, Project_ID)

	clustered_class_All_Output_Path = "%s/%s_clustered_class_All.txt" % (SigProfilerClusters_Result, Project_ID)
	clustered_class_All_Output_File = open(clustered_class_All_Output_Path, 'w')
	clustered_class_All_Output_Header = "project	samples	ID	genome	mutType	chr	start	end	ref	alt	mutClass	IMDplot	group	IMD	VAF/CCF	subclass\n"
	clustered_class_All_Output_File.write(clustered_class_All_Output_Header)

	cN_Path = "%soutput/vcf_files_corrected/%s_nonClustered/SNV/%s_nonClustered_Vaf.txt" % (Run_SigProfilerClusters_Output,Project_ID, Project_ID)


	if os.path.exists(c1a_Path):
		File = open(c1a_Path)
		for line in File:
			if re.match(r'project	samples', line):
				pass
			else:
				clustered_class_All_Output_File.write(line)
		File.close()
	else:
		print("Error:  There is no file: %s" % c1a_Path)

	if os.path.exists(c1b_Path):
		File = open(c1b_Path)
		for line in File:
			if re.match(r'project	samples', line):
				pass
			else:
				clustered_class_All_Output_File.write(line)
		File.close()
	else:
		print("Error:  There is no file: %s" % c1b_Path)

	if os.path.exists(c1c_Path):
		File = open(c1c_Path)
		for line in File:
			if re.match(r'project	samples', line):
				pass
			else:
				clustered_class_All_Output_File.write(line)
		File.close()
	else:
		print("Error:  There is no file: %s" % c1c_Path)

	if os.path.exists(c2_Path):
		File = open(c2_Path)
		for line in File:
			if re.match(r'project	samples', line):
				pass
			else:
				clustered_class_All_Output_File.write(line)
		File.close()
	else:
		print("Error:  There is no file: %s" % c2_Path)

	if os.path.exists(c3_Path):
		File = open(c3_Path)
		for line in File:
			if re.match(r'project	samples', line):
				pass
			else:
				clustered_class_All_Output_File.write(line)
		File.close()
	else:
		print("Error:  There is no file: %s" % c3_Path)

	if os.path.exists(cN_Path):
		File = open(cN_Path)
		for line in File:
			if re.match(r'project	samples', line):
				pass
			elif re.match(r"\n", line):
				pass
			else:
				ss = line.strip().split("	")
				
				Project = ss[0]
				Samples = ss[1]
				ID = ss[2]
				genome = ss[3]
				mutType = ss[4]
				Chr = ss[5]
				start = ss[6]
				end = ss[7]
				ref = ss[8]
				alt = ss[9]
				mutClass = ss[10]
				IMDplot = ss[11]
				IMD = ss[12]
				group = "N"
				VAF = ss[13]
				subclass = "Non-clust"
				Final_Str = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (Project, Samples, ID, genome, mutType, Chr, start, end, ref, alt, mutClass, IMDplot, group, IMD, VAF, subclass)
				clustered_class_All_Output_File.write(Final_Str)
		File.close()
	else:
		print("Error:  There is no file: %s" % c3_Path)


	clustered_class_All_Output_File.close()



	### 005 Integrate SigProfilerClusters results of image
	print(" \n####### Arrange Final Results: Clusters image\n")
	#print("\n\nHello World 1")

	clustered_image_Output_Dir = "%soutput/plots" % (Run_SigProfilerClusters_Output)
	image_CP_cmd = "cp %s/* %s" % (clustered_image_Output_Dir, SigProfilerClusters_Result)
	print(image_CP_cmd)
	os.system(image_CP_cmd)
	#print("\n\nHello World 2")


	### 006 Clean temp files.
	print(" \n####### 006 Clean temp files.\n")
	Cluster_Input_Dir = "%sinput" % (Run_SigProfilerClusters_Output)
	Cluster_Output_Dir = "%soutput" % (Run_SigProfilerClusters_Output)

	Cluster_logs_Dir = "%slogs" % (Run_SigProfilerClusters_Output)
	Cluster_Statistics_File = "%sStatistics.txt" % (Run_SigProfilerClusters_Output)
	Cluster_VCF_File = "%s*vcf" % (Run_SigProfilerClusters_Output)
	
	Delete_Dir(Cluster_Input_Dir)
	Delete_Dir(Cluster_Output_Dir)
	Delete_Dir(Cluster_logs_Dir)
	Delete_File(Cluster_Statistics_File)
	Delete_File(Cluster_VCF_File)


if __name__ == "__main__":
	#If_Compressed()
	
	Parse_Options()


########################################################################
###########################  Usage Example #############################
########################################################################


### Usage for csv
# python mSigPortal_Profiler_Extraction.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output -g GRCh37 -t WGS

### Usage for tsv
# python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS

### Collpase ###
# python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS -c True
# python mSigPortal_Profiler_Extraction.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output -g GRCh37 -t WGS -c True
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS -c True

### Usage for vcf
# python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS


### Usage for Filter
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal -i Demo_input/demo_input_single.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f tsv -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS



### Usage for Compressed File
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal@- -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_single.vcf.gz -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal@- -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_single.vcf.tar.gz -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal@- -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_single.vcf.tar -p Project -o Test_Output-6-22 -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f catalog_tsv -i Demo_input/demo_input_catalog.tsv.zip -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f vcf -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_single.vcf.tar -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f catalog_csv -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_catalog.csv.zip -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f catalog_tsv -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_catalog.tsv.zip -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv.zip -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f csv -i Demo_input/demo_input_multi.csv.zip -p Project -o Test_Output -g GRCh37 -t WGS


### Usage for vcf_split_all_filter File
# python mSigPortal_Profiler_Extraction.py -f vcf -s True -i /Users/sangj2/z-0-Projects/2-mSigPortal/Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output -g GRCh37 -t WGS -s True

### Final 
# time python mSigPortal_Profiler_Extraction_v28.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS -c True
# time python mSigPortal_Profiler_Extraction_v28.py -f vcf -F PASS@alt_allele_in_normal -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS -c True

# python mSigPortal_Profiler_Extraction_v32.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_multi.vcf.gz -p Project -o Test_Output -g GRCh37 -t WGS -c True

### Usage for catalog_csv
# python mSigPortal_Profiler_Extraction_v32.py -f catalog_csv -i Demo_input/demo_input_catalog.csv -p Project -o Test_Output_Catlog_CSV -g GRCh37 -t WGS

### Usage for catalog_tsv
# python mSigPortal_Profiler_Extraction_V32.py -f catalog_tsv -i Demo_input/demo_input_catalog.tsv -p Project -o z-9-Test_Output_Catlog_TSV -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction_V32.py -f catalog_csv -i Demo_input_2/tmp2.csv -p Project -o z-9-Test_Output_Catlog_TSV -g GRCh37 -t WGS


### Usage for catalog_maf
# python mSigPortal_Profiler_Extraction_v32.py -f maf -i Demo_input/demo_input_multi_MAF.txt -p Project -o z-9-Test_Output_Catlog_TSV -g GRCh37 -t WGS


### Usage for catalog_tsv RS32
# python mSigPortal_Profiler_Extraction_v32.py -f catalog_tsv -i Demo_input_2/breast_cancer_samples_example.RS32.all -p Project -o z-9-Test_Catalog_tsv_R32 -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction_v32.py -f catalog_tsv -i Demo_input_3/Cancer_Reference_Signatures_GRCh37_RS32 -p Project -o z-9-Test_Catalog_tsv_R32_TW1 -g GRCh37 -t WGS
# python mSigPortal_Profiler_Extraction_v32.py -f catalog_tsv -i Demo_input_3/Cancer_Reference_Signatures_GRCh37_RS32 -p Project -o Cancer_Reference_Signatures_GRCh37_RS32 -g GRCh37 -t WGS



### Usage for catalog_tsv CNV48
# python mSigPortal_Profiler_Extraction_v35.py -f catalog_tsv -i Demo_input_2/breast_cancer_samples_example.CNV48.all -p Project -o z-9-Test_Catalog_tsv_CNV48 -g GRCh37 -t WGS


# python mSigPortal_Profiler_Extraction_V35.py -f catalog_tsv -i /Users/sangj2/0-Project/2-1-New-Package/z-8-Reports/R-3-2021-01-19/Test-02-Plot-Signatures/1-Raw-To-Catalog.txt -p Project -o /Users/sangj2/0-Project/2-1-New-Package/z-8-Reports/R-3-2021-01-19/Test-02-Plot-Signatures/z-9-Test_Output_Catlog_TSV -g GRCh37 -t WGS

# python mSigPortal_Profiler_Extraction_V35.py -f catalog_tsv -i /Users/sangj2/0-Project/3-Tongwu/0-mSigPortal/0-mSigPortal_Profiler_Extraction/Demo_input_4_2022_0222/SBS96.txt -p Project -o /Users/sangj2/0-Project/3-Tongwu/0-mSigPortal/0-mSigPortal_Profiler_Extraction/z-10-Demo_input_4_2022_0222 -g GRCh37 -t WGS


# python mSigPortal_Profiler_Extraction_V35.py -f catalog_tsv -i Demo_input_Test/Demo_input_4_2022_0222/DBS78.txt -p Project -o test-14 -g GRCh37 -t WGS



### Usage for Run_SigProfilerClusters
# python3 mSigPortal_Profiler_Extraction_V36.py -f vcf -i Demo_input_Test/demo_input_multi.vcf -p Project -o test-19 -g GRCh37 -t WGS -C True

# python3 mSigPortal_Profiler_Extraction_V36.py -f vcf -i Demo_input_Test/demo_input_multi.vcf.gz -p Project -o test-20 -g GRCh37 -t WGS -C True


# python3 mSigPortal_Profiler_Extraction_V36.py -f tsv -i Demo_input_Test/demo_input_multi.tsv -p Project -o test-17-tsv -g GRCh37 -t WGS





