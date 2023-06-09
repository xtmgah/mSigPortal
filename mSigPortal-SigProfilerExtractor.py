#!/usr/bin/python
#encoding=utf8
'''
Name: 		mSigPortal-SigProfilerExtractor (pipeline for running SigProfilerExtractor)
Version: 	V6
Data: 		Mar-29-2023

'''
import re,os,argparse,sys,time
from SigProfilerExtractor import sigpro as sig

def Parser():
	parser = argparse.ArgumentParser(add_help=False)
	
	### 001 Input Data
	parser.add_argument('--input_type', required=True, type=str, help = "Define the formats of input data. Only 'vcf' and 'matirx' are supported.")
	parser.add_argument('--output', required=True, nargs='?', const='mSigPortal_%s' % time.strftime('%Y%m%d%H%M%S',time.localtime(time.time())),type=str, help = "Define the absolute path for Output Dir.")
	parser.add_argument('--input_data', required=True, type=str, help = "Define the absolute path for input file")
	parser.add_argument('--reference_genome', required=True, default='GRCh37', type=str, help = "Define the version of Genome_Building, default:GRCh37.")
	parser.add_argument('--opportunity_genome', required=False, nargs='?', const='GRCh37', type=str, help = "Define the version of Genome_Building, default:GRCh37.")
	parser.add_argument('--context_type', required=False, default='SBS96', type=str, help = "Define the context_type: 'SBS96', 'DINUC', 'ID'. Notes: selected context_type should be matched with selected Signature_Database.")
	parser.add_argument('--exome', required=False, type=str, default="False", help = "Defines if the exomes will be extracted. The default value is \"False\".")
	parser.add_argument('--signature_database', required=True, nargs='?', type=str, help = "Define the absolute path for Signature_Database. Notes: selected Signature_Database should be matched with selected context_type.")

	### 002 Execution
	parser.add_argument('--cpu', required=False, default='-1', type=int, help = "The number of processors to be used to extract the signatures.")
	parser.add_argument('--gpu', required=False, default="False", type=str, help = "Defines if the GPU resource will used if available. ")
	parser.add_argument('--batch_size', required=False, default="1", type=str, help = "Will be effective only if the GPU is used.")


	### 003 NMF Replicates	
	parser.add_argument('--minimum_signatures', required=True, type=str, help = "Number of minimum_signatures.")
	parser.add_argument('--maximum_signatures', required=True, type=str, help = "Number of maximum_signatures.")
	parser.add_argument('--nmf_replicates', required=True, type=str, help = "Number of nmf_replicates")
	parser.add_argument('--resample', required=False, type=str, default="True", help = "Default is True. If True, add poisson noise to samples by resampling.")
	parser.add_argument('--seeds', required=False, type=str, default="random", help = "It can be used to get reproducible resamples for the NMF replicates. ")


	### 004 NMF Engines	
	parser.add_argument('--matrix_normalization', required=False, type=str, default="gmm", help = "Method of normalizing the genome matrix before it is analyzed by NMF.")
	parser.add_argument('--nmf_init', required=False, type=str, default="random", help = "The initialization algorithm for W and H matrix of NMF.")
	parser.add_argument('--precision', required=False, type=str, default="single", help = "Values should be single or double.")
	parser.add_argument('--nmf_tolerance', required=False, type=float, default=1e-15, help = "Value defines the tolerance to achieve to converge.")
	parser.add_argument('--min_nmf_iterations', required=True, type=str, help = "Number of min_nmf_iterations.")
	parser.add_argument('--max_nmf_iterations', required=True, type=str, help = "Number of max_nmf_iterations.")
	parser.add_argument('--nmf_test_conv', required=True, type=str, help = "Number of nmf_test_conv")


	### 005 Solution Estimation Thresholds
	parser.add_argument('--stability', required=False, type=float, default=0.8, help = "The cutoff thresh-hold of the average stability. ")
	parser.add_argument('--min_stability', required=False, type=float, default=0.2, help = "The cutoff thresh-hold of the minimum stability. ")
	parser.add_argument('--combined_stability', required=False, type=float, default=1.0, help = "The cutoff thresh-hold of the combined stability. ")
	parser.add_argument('--allow_stability_drop', required=False, type=str, default="False", help = "Defines if solutions with a drop in stability with respect to the highest stable number of signatures will be considered. ")

	
	### 006 Decomposition
	parser.add_argument('--make_decomposition_plots', required=False, type=str, default="True", help = "If True, Denovo to Cosmic sigantures decompostion plots will be created as a part the results.")
	parser.add_argument('--collapse_to_SBS96', required=False, type=str, default="True", help = "If True, SBS288 and SBS1536 Denovo signatures will be mapped to SBS96 reference signatures.")

	### 007 Others
	parser.add_argument('--get_all_signature_matrices', required=False, type=str, default="True", help = "If True, the Ws and Hs from all the NMF iterations are generated in the output.")
	parser.add_argument('--export_probabilities', required=False, type=str, default="True", help = "If False, then doesn't create the probability matrix.")

	
	
	
	
	

	args = parser.parse_args()
	return args.input_type, args.output, args.input_data, args.reference_genome, args.opportunity_genome, args.context_type, args.exome, args.signature_database, args.cpu, args.gpu, args.minimum_signatures, args.maximum_signatures, args.min_nmf_iterations, args.max_nmf_iterations, args.nmf_test_conv, args.nmf_replicates, args.resample, args.seeds, args.matrix_normalization, args.nmf_init, args.precision, args.nmf_tolerance, args.stability, args.min_stability, args.combined_stability, args.allow_stability_drop, args.make_decomposition_plots, args.collapse_to_SBS96, args.get_all_signature_matrices, args.export_probabilities



def main_function():    
	input_type, Output_Dir, Input_Path, reference_genome, i_opportunity_genome, i_context_type, i_exome, i_Signature_Database, i_CPU, i_GPU, i_minimum_signatures, i_maximum_signatures, i_min_nmf_iterations, i_max_nmf_iterations, i_nmf_test_conv, i_nmf_replicates, i_resample, i_seeds, i_matrix_normalization, i_nmf_init, i_precision, i_nmf_tolerance, i_stability, i_min_stability, i_combined_stability, i_allow_stability_drop, i_make_decomposition_plots, i_collapse_to_SBS96, i_get_all_signature_matrices, i_export_probabilities = Parser()



	print("\n 001 Parse Parameter:\n")

	#print(i_minimum_signatures)
	i_minimum_signatures = int(i_minimum_signatures)
	#print(i_maximum_signatures)
	i_maximum_signatures = int(i_maximum_signatures)
	# print(i_min_nmf_iterations)
	i_min_nmf_iterations = int(i_min_nmf_iterations)
	# print(i_max_nmf_iterations)
	i_max_nmf_iterations = int(i_max_nmf_iterations)
	i_nmf_test_conv = int(i_nmf_test_conv)
	
	i_nmf_replicates = int(i_nmf_replicates)
	
	print("input_type:", input_type)
	print("Output_Dir:",Output_Dir)
	print("Input_Path:", Input_Path)
	print("reference_genome:", reference_genome)
	print("Signature_Database:", i_Signature_Database)
	print("context_type:", i_context_type)
	print("minimum_signatures:", i_minimum_signatures)
	print("min_nmf_iterations:", i_min_nmf_iterations)
	print("max_nmf_iterations", i_max_nmf_iterations)
	print("nmf_test_conv", i_nmf_test_conv)
	
	if i_exome == "False":
		i_exome = False
	else:
		i_exome = True
	print("exome:", i_exome)
	
	i_CPU = int(i_CPU)
	print("CPU:", i_CPU)

	if i_GPU == "False":
		i_GPU = False
	else:
		i_GPU = True
	print("GPU:", i_GPU)


	if i_resample == "False":
		i_resample = False
	else:
		i_resample = True
	print("i_resample:", i_resample)

	print("matrix_normalization:", i_matrix_normalization)

	print("nmf_init:", i_nmf_init)

	print("nmf_init:", i_nmf_init)

	print("precision:", i_precision)

	print("nmf_tolerance:", i_nmf_tolerance)

	print("stability:", i_stability)

	print("min_stability:", i_min_stability)

	print("combined_stability:", i_combined_stability)

	print("allow_stability_drop:", i_allow_stability_drop)

	if i_make_decomposition_plots == "False":
		i_make_decomposition_plots = False
	else:
		i_make_decomposition_plots = True
	print("make_decomposition_plots:", i_make_decomposition_plots)

	if i_collapse_to_SBS96 == "False":
		i_collapse_to_SBS96 = False
	else:
		i_collapse_to_SBS96 = True
	print("collapse_to_SBS96:", i_collapse_to_SBS96)


	if i_get_all_signature_matrices == "False":
		i_get_all_signature_matrices = False
	else:
		i_get_all_signature_matrices = True
	print("get_all_signature_matrices:", i_get_all_signature_matrices)


	if i_export_probabilities == "False":
		i_export_probabilities = False
	else:
		i_export_probabilities = True
	print("export_probabilities:", i_export_probabilities)


	print("\n 002 Running sigProfilerExtractor:\n")
	
	sig.sigProfilerExtractor(input_type, Output_Dir, Input_Path, reference_genome=reference_genome, 
	opportunity_genome=i_opportunity_genome, 
	Signature_Database=i_Signature_Database, 
	context_type=i_context_type,
	exome = i_exome,
	cpu = i_CPU,
	gpu = i_GPU,
	minimum_signatures = i_minimum_signatures, 
	maximum_signatures = i_maximum_signatures, 
	min_nmf_iterations = i_min_nmf_iterations, 
	max_nmf_iterations = i_max_nmf_iterations, 
	nmf_test_conv = i_nmf_test_conv, 
	nmf_replicates=i_nmf_replicates,
	resample = i_resample,
	seeds = i_seeds,
	matrix_normalization = i_matrix_normalization,
	nmf_init = i_nmf_init,
	precision = i_precision,
	nmf_tolerance = i_nmf_tolerance,
	stability = i_stability,
	min_stability = i_min_stability, 
	combined_stability = i_combined_stability,
	allow_stability_drop = i_allow_stability_drop,
	make_decomposition_plots = i_make_decomposition_plots,
	collapse_to_SBS96 = i_collapse_to_SBS96,
	get_all_signature_matrices = i_get_all_signature_matrices,
	export_probabilities = i_export_probabilities,
	)


if __name__=="__main__":
	main_function()








'''
python mSigPortal-SigProfilerExtractor.py --input_type vcf --input_data  vcftest --output z-Test-Matrix-6-SBS --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.1_SBS_GRCh37_exome.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10
python mSigPortal-SigProfilerExtractor-V5.py --input_type vcf --input_data  vcftest --output z-Test-Matrix-6-DINUC --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.2_DBS_GRCh37_exome.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type DINUC --nmf_replicates 10
python mSigPortal-SigProfilerExtractor-V5.py --input_type vcf --input_data  vcftest --output z-Test-Matrix-6-ID --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.2_ID_GRCh37.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type ID --nmf_replicates 10

python mSigPortal-SigProfilerExtractor.py --input_type matrix --input_data Demo_input_Test/vcftest.SBS96.all --output z-Test-Matrix-6-SBS --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.1_SBS_GRCh37_exome.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10


python mSigPortal-SigProfilerExtractor.py --input_type matrix --input_data Workshop-Session5/output/SBS/S5.SBS96.all --output Workshop-Session5-SigProfilerExtractor --reference_genome GRCh38 --signature_database Reference_Signatures/GRCh38/COSMIC_v3.3_SBS_GRCh38.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10

# not working :
# python mSigPortal-SigProfilerExtractor.py --input_type matrix --input_data Workshop-Session5/output/SBS/S5.SBS96.all --output Workshop-Session5-SigProfilerExtractor --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.3_SBS_GRCh37.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10

python mSigPortal-SigProfilerExtractor.py --input_type matrix --input_data Workshop-Session5-37/output/SBS/S5.SBS96.all --output Workshop-Session5-SigProfilerExtractor --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.1_SBS_GRCh37.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10



# The following code is working:

python mSigPortal-SigProfilerExtractor.py --input_type matrix --input_data Workshop-Session5-37/output/SBS/S5.SBS96.all --output Workshop-Session5-SigProfilerExtractor --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.1_SBS_GRCh37.txt --minimum_signatures 1 --maximum_signatures 2 --min_nmf_iterations 2 --max_nmf_iterations 4 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10
python mSigPortal-SigProfilerExtractor.py --input_type matrix --input_data Workshop-Session5-37/output/SBS/S5.SBS96.all --output Workshop-Session5-SigProfilerExtractor-Test-2 --reference_genome GRCh37 --signature_database Reference_Signatures/GRCh37/COSMIC_v3.1_SBS_GRCh37.txt --minimum_signatures 1 --maximum_signatures 6 --min_nmf_iterations 2 --max_nmf_iterations 10 --nmf_test_conv 2 --context_type SBS96 --nmf_replicates 10



'''



