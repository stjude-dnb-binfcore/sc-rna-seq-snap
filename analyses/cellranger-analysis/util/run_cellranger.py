#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author:
	Name: 			Cody Alexander Ramirez
	Email: 			cody.ramirez@stjude.org
	Affiliation: 	St. Jude Children's Research Hospital, Memphis, TN
	Date: 			June 22nd, 2022 - Updated on July 29th, 2024 to be compatible with CellRanger v8.0.1
"""


import os, sys, argparse

#Beginning the argparse module to extract command-line arguments
parser = argparse.ArgumentParser(description="This is a wrapper script that will submit cellranger commands to the HPC cluster with all required parameters. It accepts either a	single sample's info input (e.g. id, sample and fastqs) or a file containing the info for multiple samples")
#Creating the following optional arguments; some have default values
parser.add_argument('--id', help='(Required). A unique run ID string: e.g. DYE3755', type=str)
parser.add_argument('--output_dir', help='Optional. Path to a custom output directory for storing the results of the run. When this argument is specified, the outputs are placed in the user designated directory, following the format: /path/to/<custom name>/outs/. If excluded, the outputs are directed to the default path: /path/to/ID/outs/.', type=str)
#--description Optional. Sample description to embed in output files. We are not going to add this parameter just yet as we don't request description information from collaborators at this time
parser.add_argument('--transcriptome', help='(Required). Path to the Cell Ranger compatible transcriptome reference or you can use pre-existing reference genomes such as \"GRCh38\", \"GRCm39\" or \"GRCh38ANDGRCm39\".', type=str)
parser.add_argument('--fastqs', help='(Required). The full path to the corresponding fafstq file.', type=str)
#--project Optional. Name of the project folder within a mkfastq, bcl2fastq, or bcl-convert-generated folder from which to pick FASTQs. We are not going to add this parameter just yet as we don't request project information from collaborators at this time.
parser.add_argument('--sample', help='(Required). Allowable characters in sample names are letters, numbers, hyphens and underscores.', type=str)
parser.add_argument('--create_bam', help='(Required). <true|false> Enable or disable BAM file generation. Setting --create-bam=false reduces the total computation time and the size of the output directory (BAM file not generated). We recommend setting --create-bam=true if unsure. See https://10xgen.com/create-bam for additional guidance [possible values: true, false]', type=str)
parser.add_argument('--lanes', help='Optional. Only use FASTQs from selected lanes.', type=int)
parser.add_argument('--libraries', help='Optional. Path to a Libraries CSV file declaring FASTQ paths and library types of input libraries. Required for Gene Expression + Feature Barcode analysis. When using this argument, neither --fastqs nor --sample must be used. Additionally, this argument is inappropriate for gene expression-only analysis; for such cases, use --fastqs instead.', type=str)
parser.add_argument('--feature_ref', help='Required for Feature Barcode analysis. Path to a Feature Reference CSV file declaring the Feature Barcode reagents used in the experiment.', type=str)
parser.add_argument('--expect_cells', help='Optional. Override the pipeline’s auto-estimate. See cell calling algorithm for details on how this parameter is used. If used, enter the expected number of recovered cells.', type=int)
parser.add_argument('--force_cells', help='Optional. Force pipeline to use this number of cells, bypassing the cell detection algorithm. Use this if the number of cells estimated by Cell Ranger is not consistent with the barcode rank plot.', type=int)
parser.add_argument('--r1_length', help='Optional. Limit the length of the input Read 1 sequence of Gene Expression (and any Feature Barcode) library to the first N bases, where N is a user-supplied value. Note that the length includes the 10x Barcode and UMI sequences so do not set this below 26 for Single Cell 3′ v2 or Single Cell 5′v1/v2. This and --r2-length are useful options for determining the optimal read length for sequencing.', type=int)
parser.add_argument('--r2_length', help='Optional. Limit the length of the input R2 sequence to the first N bases, where N is a user-supplied value. Trimming occurs before sequencing metrics are computed and therefore, limiting R2 read length may affect Q30 scores (true, false).', type=int)
parser.add_argument('--include_introns', help='Optional. Set to false to exclude intronic reads in count. Including introns in analysis is recommended to maximize sensitivity. Refer to the guide on including introns for gene expression analysis for comprehensive recommendations. Default: true (in Cell Ranger v7.0 and later, intronic reads are counted by default)', action="store_true")
parser.add_argument('--chemistry', help='Optional. Assay configuration. By default the assay configuration is detected automatically, which is the recommended mode. You should only specify chemistry if there is an error in automatic detection or for some special cases.', type=str)
parser.add_argument('--check_library_compatibility', help='Optional. This option allows users to disable the check that evaluates 10x Barcode overlap between libraries when multiple libraries are specified (e.g., Gene Expression + Antibody Capture). Setting this option to false will disable the check across all library combinations. We recommend running this check (default), however if the pipeline errors out, users can bypass the check to generate outputs for troubleshooting (true, false; default: true).', type=str)
parser.add_argument('--nosecondary', help='Optional. Adding --nosecondary to the command disables the generation of secondary analysis results, such as clustering and dimensionality reduction (tSNE/UMAP) and the .cloupe output file. This flag does not require a true or false value to be specified.', action="store_true")
parser.add_argument('--min_crispr_umi', help='Optional. Set the minimum number of CRISPR guide RNA UMIs required for protospacer detection. Default: 3 UMIs (in Cell Ranger v8.0 and later). If a lower or higher sensitivity is desired for detection, this value can be customized according to specific experimental needs. Applicable only to datasets that include a CRISPR Guide Capture library.', type=int)
parser.add_argument('--file', help='Optional. Tsv file containing the following	header and corresponding info (#ID SAMPLE FASTQ) for batch submission.', type=argparse.FileType('r'))


#Converts argument strings to objects and assigns them as attributes of the namespace; e.g. --id -> args.id
args = parser.parse_args()

#If all info for a single sample command are enter, then submit a single cellranger count command to the HPC
if args.id != None and args.sample != None and args.fastqs != None and args.transcriptome != None:
	
	SingleCommand = 'bsub -P CellRangerCount -J CellRangerCount_' + args.id + ' -q standard -o ./results/01_logs/' + args.id + '.out -e ./results/01_logs/' + args.id + '.err -n 4 -R "rusage[mem=4000] span[hosts=1]" ' + 'cellranger count --id=' + args.id + ' --fastqs=' + args.fastqs + ' --sample=' + args.sample + ' --transcriptome=' + args.transcriptome + ' --create-bam=' + args.create_bam + ' --jobmode=lsf '
	
	if args.lanes != None:
		SingleCommand = SingleCommand + '--lanes ' + str(args.lanes) + ' '
	if args.libraries != None:
		SingleCommand = SingleCommand + '--libraries ' + args.libraries + ' '
	if args.feature_ref != None:
		SingleCommand = SingleCommand + '--feature_ref ' + args.feature_ref + ' '
	if args.expect_cells != None:
		SingleCommand = SingleCommand + '--expect-cells ' + str(args.expect_cells) + ' '
	if args.force_cells != None:
		SingleCommand = SingleCommand + '--force-cells ' + str(args.force_cells) + ' '
	if args.r1_length != None:
		SingleCommand = SingleCommand + '--r1-length ' + str(args.r1_length) + ' '
	if args.r2_length != None:
		SingleCommand = SingleCommand + '--r2-length ' + str(args.r2_length) + ' '
	if args.include_introns == True:
		SingleCommand = SingleCommand + '--include-introns '
	if args.chemistry != None:
		SingleCommand = SingleCommand + '--chemistry ' + args.chemistry + ' '
	if args.check_library_compatibility != None:
		SingleCommand = SingleCommand + '--check-library-compatibility ' + args.check_library_compatibility + ' '
	if args.nosecondary == True:
		SingleCommand = SingleCommand + '--nosecondary '
	if args.min_crispr_umi != None:
		SingleCommand = SingleCommand + '--min-crispr-umi ' + str(args.min_crispr_umi) + ' '
	if args.output_dir != None:
		SingleCommand = SingleCommand + '--output-dir ' + str(args.output_dir) + args.id + ' '
		
		#Submit command to the terminal
		# print(SingleCommand)
		os.system(SingleCommand)
		

#If a file is submitted, then submit a cellranger count command for every sample in the file to the HPC
elif args.file != None:
	#Reading through each line in the file
	for line in args.file:
		#Accessing the header of the file
		if line.startswith('ID'):
			#Identifying the index of the columns
			IdIndex = line.strip().split('\t').index('ID')
			SampleIndex = line.strip().split().index('SAMPLE')
			FastqIndex = line.strip().split().index('FASTQ')
		else:
			#Collecting the current sample's required infos
			Id = line.strip().split('\t')[IdIndex]
			Sample = line.strip().split('\t')[SampleIndex]
			Fastq = line.strip().split('\t')[FastqIndex]

			#If cellranger output already exist for the file, then skip
			if os.path.isfile(Id + '/outs/filtered_feature_bc_matrix.h5'):
				print('CellRanger has already been run for ' + Id + '. Skipping this sample and proceeding to next sample.')
			#Otherwise, submit a cellranger count command for the current sample
			else:
				print('Still need to run CellRanger for ' + Id + ' Submitting command now.')

				SingleCommand = 'bsub -P CellRangerCount -J CellRangerCount_' + Id + ' -q standard -o ./results/01_logs/' + Id + '.out -e ./results/01_logs/' + Id + '.err -n 4 -R "rusage[mem=4000] span[hosts=1]" ' + 'cellranger count --id=' + Id + ' --fastqs=' + Fastq + ' --sample=' + Sample + ' --transcriptome=' + args.transcriptome + ' --create-bam=' + args.create_bam + ' --jobmode=lsf '
				
				if args.lanes != None:
					SingleCommand = SingleCommand + '--lanes ' + str(args.lanes) + ' '
				if args.libraries != None:
					SingleCommand = SingleCommand + '--libraries ' + args.libraries + ' '
				if args.feature_ref != None:
					SingleCommand = SingleCommand + '--feature_ref ' + args.feature_ref + ' '
				if args.expect_cells != None:
					SingleCommand = SingleCommand + '--expect-cells ' + str(args.expect_cells) + ' '
				if args.force_cells != None:
					SingleCommand = SingleCommand + '--force-cells ' + str(args.force_cells) + ' '
				if args.r1_length != None:
					SingleCommand = SingleCommand + '--r1-length ' + str(args.r1_length) + ' '
				if args.r2_length != None:
					SingleCommand = SingleCommand + '--r2-length ' + str(args.r2_length) + ' '
				if args.include_introns == True:
					SingleCommand = SingleCommand + '--include-introns '
				if args.chemistry != None:
					SingleCommand = SingleCommand + '--chemistry ' + args.chemistry + ' '
				if args.check_library_compatibility != None:
					SingleCommand = SingleCommand + '--check-library-compatibility ' + args.check_library_compatibility + ' '
				if args.nosecondary == True:
					SingleCommand = SingleCommand + '--nosecondary '
				if args.min_crispr_umi != None:
					SingleCommand = SingleCommand + '--min-crispr-umi ' + str(args.min_crispr_umi) + ' '
				if args.output_dir != None:
					SingleCommand = SingleCommand + '--output-dir ' + str(args.output_dir) + Id + ' '
					
					#Submit command to the terminal
					# print(SingleCommand)
					os.system(SingleCommand)
					
					#If improper arguments are used, then print the help statement and exits the program
					
				else:
					parser.print_help()
					sys.exit()
						
						
