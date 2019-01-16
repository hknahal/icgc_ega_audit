#!/usr/bin/python
# Hardeep Nahal
# Jan 2016
# parse_xml_seqstrat.py
# Parses EGA XML files to generate an EGA-DCC Audit


import os
import string
import xmltodict
import sys, getopt
import fnmatch
import re
import collections
import bz2
import gzip
import json
import psycopg2
from configparser import ConfigParser
#import cStringIO
import pycurl
import json
from io import BytesIO

import parse_ega_xmls as ega
import parse_dcc_data as dcc
#from EGA import DCCSamples, Dataset
import EGA
import get_file_info

ega_dataset = ''
ega_study = ''
project_code = ''
project_name = ''
mappingFileName = ''
output_dir = ''
ega_datasets = []
ega_file_info = {}
ega_file_sizes = {}
rootDir = ''
datasetFileInfo = {}
test = open("Testing.tsv", "w")
donor_summary = open("Donor_Summary_3Dec2018.tsv", "a")

try:
   opts, args = getopt.getopt(sys.argv[1:], "hm:o:r:", ["mapping_file","output_dir", "root_dir"])
except getopt.GetoptError:
      print('parse_xmls.py -m <mapping_file> -o <output_dir> -r <root ega metadata directory>')
      sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print('parse_xmls.py -m <mapping_file> -o <output_dir>')
      sys.exit()
   elif opt in ("-m", "--mapping_file"):
      mappingFileName = arg
   elif opt in ("-o", "--output_dir"):
      output_dir = arg
   elif opt in ("-r", "--root__dir"):
      rootDir = arg

rootDir="/Users/hnahal/EGA_Metadata/ftp-private.ebi.ac.uk/ICGC_metadata"
#rootDir="/Users/hnahal/EGA_Metadata/cleaned_code/using_egaf_accession/mela-au_ega_datasets"
#rootDir="/Users/hnahal/EGA_Metadata/ftp-private.ebi.ac.uk/dcc_api"
submissions_dir = "/Users/hnahal/submission_files/ICGC28"
#submissions_dir = "/Users/hnahal/Release27"
#submissions_dir = "/Users/hnahal/Release26/test_submissions"
#submissions_dir = "/Users/hnahal/submission_files/test"
#submissions_dir = "ICGC24"
#fileInfo = "EGA_File_Sizes_20Sept2018.tsv"

specimen_types = {}
specimen_type_list = { '101': 'Normal - solid tissue', '102': 'Normal - blood derived', '103': 'Normal - bone marrow', '104': 'Normal - tissue adjacent to primary', '105': 'Normal - buccal cell', '106': 'Normal - EBV immortalized', '107': 'Normal - lymph node', '108': 'Normal - other', '109': 'Primary tumour - solid tissue', '110': 'Primary tumour - blood derived (peripheral blood)', '111': 'Primary tumour - blood derived (bone marrow)', '112': ' Primary tumour - additional new primary', '113': 'Primary tumour - other', '114': 'Recurrent tumour - solid tissue', '115': 'Recurrent tumour - blood derived (peripheral blood)', '116': 'Recurrent tumour - blood derived (bone marrow)', '117': 'Recurrent tumour - other', '118': 'Metastatic tumour - NOS', '119': 'Metastatic tumour - lymph node', '120': 'Metastatic tumour - metastasis local to lymph node', '121': 'Metastatic tumour - metastasis to distant location', '122': 'Metastatic tumour - additional metastatic', '123': 'Xenograft - derived from primary tumour', '124': 'Xenograft - derived from tumour cell line', '125': 'Cell line - derived from tumour', '126': 'Primary tumour - lymph node' , '127': 'Metastatic tumour - other', '128': 'Cell line - derived from xenograft tumour'}

project_info = collections.defaultdict(dict)
summary_file = open("EGA_DCC_Summary_3Dec2018.tsv", "a")
#summary_file = open("EGA_DCC_Summary_MELA-AU.txt", "a")
#summary_file.write("Project\tCountry\tProject Name\tSequencing Strategy\tTotal Number of analyzed samples\tTotal found in EGA\tTotal Not Found in EGA\tDonors with EGA Data\tDonors without EGA Data\tTotal Donors\tSamples with Files but missing EGAF IDs\tExclude PCAWG samples\tPCAWG samples not found in non-PCAWG EGA datasets, but exist in PCAWG EGA datasets\tTotal PCAWG samples\t\n")
logfile = open("EGA_DCC_Audit.log", "w")

# some projects use the value in specimen_types. Store keys and values
for code in specimen_type_list:
   value = specimen_type_list[code]
   specimen_types[code] = value
   specimen_types[value] = value

#   if code.isnumeric():
#      logfile.write("code %s is numeric\n"%code)
#      specimen_types[code] = value
#   else:
#      logfile.write("value = %s\n"%value)
#      specimen_types[code] = code

for key in specimen_types:
   logfile.write("key = %s\tvalue = %s\n"%(key, specimen_types[key]))

# Connect to Postgres DB
def config(db_config='db_config.ini', section='postgresql'):
   parser = ConfigParser()
   #print "in config"
   parser.read(db_config)
   db = {}
   if parser.has_section(section):
      params = parser.items(section)
      for param in params:
         db[param[0]] = param[1]
      #   print "%s: %s"%(db[param[0]], param[1])
   else:
     raise Exception('Section {0} not found in the {1} file'.format(section, db_config))
   return db


def connect():
   conn = None
   cur = None
   try:
      params = config()
      #print "Connecting to the PostgreSQL database..."
      #print "params = %s"%params
      conn = psycopg2.connect(**params)
      cur = conn.cursor()
      #print "PostgreSQL database verson:"
      #cur.execute("SELECT version()")
      #db_version = cur.fetchone()
      #print(db_version)
   except (Exception, psycopg2.DatabaseError) as error:
      print(error)
   return cur

# DB query to get EGA dataset ID for a given EGA File ID
def get_egaDatasetID(ega_file_id):
   cur = None
   try:
      params = config()
      cur = connect()
      cur.execute("SELECT dataset_id FROM ega.view_ega_sample_mapping WHERE file_id = '%s'"%ega_file_id)
      #cur.execute("SELECT dataset_id FROM ega.ega_sample_mapping_1536023150 WHERE file_id = '%s'"%ega_file_id)
      rows = cur.fetchall()
   #   for row in rows:
   #      print "dataset_id from query = %s"%row
      cur.close()
      return rows[0][0]
   except (Exception, psycopg2.DatabaseError) as error:
      print(error)

# DB query to get EGA File IDs for a given sample ID
def get_egaFileID(ega_sample_id):
   cur = None
   try:
      params = config()
      cur = connect()
      cur.execute("SELECT file_id FROM ega.view_ega_sample_mapping WHERE sample_id = '%s'"%ega_sample_id)
      #cur.execute("SELECT file_id FROM ega.ega_sample_mapping_1536023150 WHERE sample_id = '%s'"%ega_sample_id)
      rows = cur.fetchall()
      #for row in rows:
      #   print "EGA File ID = %s"%row
      cur.close()
   except (Exception, psycopg2.DatabaseError) as error:
      print(error)
   return rows
 


def getFileSizes(ega_dataset_id):
   egaFiles = {}
   logfile.write("Getting File info for EGA dataset %s\n"%ega_dataset_id)
   api_endpoint = "https://ega-archive.org/metadata/v2/files?queryBy=dataset&queryId=%s&skip=0&limit=0"%ega_dataset_id
   buf = BytesIO()
   c = pycurl.Curl()
   c.setopt(c.URL, api_endpoint)
   c.setopt(c.WRITEFUNCTION, buf.write)
   c.setopt(c.TIMEOUT, 120)
   try:
      c.perform()
      body = buf.getvalue()
      #output = body.decode('iso-8859-1')
      output = body.decode('utf-8')
      file_sizes = json.loads(output)
     # transform output for easier indexing later (to search by EGA File ID)
      for item in file_sizes['response']['result']:
         ega_file_id = item['egaStableId']
         egaFiles[ega_file_id] = {}
         #egaFiles[item['egaStableId']['md5']] = item['md5']
         egaFiles[ega_file_id]['fileSize'] = item['fileSize']
   except pycurl.error as error:
      ret = error.args[0]
      logfile.write("An error occurred trying to retrieve file sizes from EGA for EGA dataset %s:\n%s\n"%(ega_dataset_id, ret))
   return egaFiles
  


# Get list of EGA datasets in directory
def get_ega_datasets(metadata_dir):
   ega_datasets = []
   for ega_dir in os.listdir(metadata_dir):
      if os.path.exists("%s/%s/xmls/samples"%(metadata_dir, ega_dir)):
         ega_datasets.append(ega_dir)
   return ega_datasets


# Parse EGA Study-Dataset mapping file
def parse_ega_dcc_mapping(mappingFileName):
   mappingFile = open(mappingFileName, "r")
   file_contents = mappingFile.readlines()
   file_contents.pop(0)
   for data in file_contents:
      if re.match("^#", data):
         continue
      data = data.rstrip("\n")
      line = data.split("\t")
      project = line[0]
      project_info[project] = {}
      project_info[project]['name'] = line[1]
      project_info[project]['country'] = line[2]
   mappingFile.close()
   return project_info   

def parse_ega_file_info(fileInfo):
   egaFileInfo = {}
   egaFileSizes = open(fileInfo, "r")
   file_contents = egaFileSizes.readlines()
   for line in file_contents:
      line = line.rstrip("\n")
      data = line.split("\t")
      egaFileInfo[data[0]] = {}
      egaFileInfo[data[0]]['ega_dataset'] = data[1]
      egaFileInfo[data[0]]['fileName'] = data[2]
      egaFileInfo[data[0]]['fileSize'] = data[3]
   egaFileSizes.close()
   return egaFileInfo

def parse_ega_mapping(egaDataset):
   ega_file_info = {}
   sample_to_egaf_mapping = {}
   egaf_to_sample_mapping = {}

   egaFile = open("~/EGA_Metadata/ftp-private.ebi.ac.uk/ICGC_metadata/%s/delimited_maps/Sample_File.map", "r")
   file_contents = egaFile.readlines()
   for line in file_contents:
      line = line.rstrip("\n")
      data = line.split("\t")
      ega_sample_name = data[0]
      ega_sample_acc = data[1]
      raw_file_name = data[2]
      ega_file_acc = data[3]
      ega_file_info[raw_file_name] = {}
      ega_file_info[raw_file_name]['ega_sample_name'] = ega_sample_name
      ega_file_info[raw_file_name]['ega_sample_acc'] = ega_sample_acc
      ega_file_info[raw_file_name]['ega_file_acc'] = ega_file_acc

      if ega_sample_name not in sample_to_egaf_mapping:
         sample_to_egaf_mapping[ega_sample_name] = []
      sample_to_egaf_mapping[ega_sample_name].append(ega_file_acc)

      if egaf_file_acc not in egaf_to_sample_mapping:
         egaf_to_sample_mapping[egaf_file_acc] = ega_sample_name
   return ega_file_info, sample_to_egaf_mapping, egaf_to_sample_mapping


def process_metadata(ega_samples, dcc_sample, seqstrat, ega_sample_info, exp_info, ega_files, dcc_sample_info, submitted_ega_files, ega_dataset, matched_sample_id):
   global sample_foundEGA
   global donor_withEGAData
   global sample_notFoundEGA
   global donor_withoutEGAData
   global fileExists_notEGAF
   global sample_mismatchedEGAF
   global donor_mismatchedEGAF
   logfile.write("Looking at dcc_sample %s [%s]\n"%(dcc_sample, seqstrat))
   logfile.write("ega_dataset = %s\n"%ega_dataset)
   # for MELA-AU (EGA dataset EGAD00001003388)
   mod_dcc_sample = dcc_sample.replace("-", "_")
   if dcc_sample in ega_samples:
      # print out matching DCC and EGA records
      # check if dcc seqstrat matches EGA libstrat
      libstrat = seqstrat
      logfile.write("library strategies for dcc sample %s = %s\n"%(dcc_sample, ega_samples[dcc_sample]))
      if ( (len(ega_samples[dcc_sample]) == 1) and ('none' in ega_samples[dcc_sample]) ):
         libstrat = "none"
      if seqstrat == 'RNA-Seq' and seqstrat not in ega_samples[dcc_sample]:
         libstrat = "none"
      elif seqstrat in ('non-NGS', 'miRNA-Seq', 'AMPLICON'):
         libstrat = "none"
      logfile.write("Checking to see if libstrat %s is in dcc_sample %s\n"%(libstrat, dcc_sample))
      if libstrat in ega_samples[dcc_sample]:
         for ega_info in ega_samples[dcc_sample][libstrat]: 
            ega_dataset_acc = ega_info['ega_dataset_acc']
            ega_study = ega.get_study(rootDir, ega_dataset_acc)
            ega_exp_acc = ega_info['ega_exp_acc']
            exp_erx_acc = ega_info['exp_erx']
            sample_ers_acc = ega_info['sample_ers']
            ega_analysis_acc = ega_info['ega_analysis_acc']
            pairedEnd = ega_info['pairedEnd']
            insertSize = ega_info['insertSize']
            referenceGenome = ega_info['referenceGenome']
            aligned = ega_info['aligned']
            if ega_dataset_acc not in datasetFileInfo:
               datasetFileInfo[ega_dataset_acc] = getFileSizes(ega_dataset_acc)
            if sample_ers_acc not in ega_sample_info:
               ega_sample_acc = ''
             #  logfile.write("[%s] No EGA sample accession (EGAN*) found for DCC sample = %s"%(project, dcc_sample))
             #  logfile.write("[%s] \tega_dataset_acc = %s"%(project, ega_dataset_acc))
             #  logfile.write("[%s]\tega_exp_acc = %s"%(project, ega_exp_acc))
             #  logfile.write("[%s]\texp_erx_acc = %s"%(project, exp_erx_acc))
             #  logfile.write("[%s]\tsample_ers_acc = %s\n"%(project, sample_ers_acc))
            else: 
               ega_sample_acc = ega_sample_info[sample_ers_acc]['ega_sample_acc']
            ega_libstrat  = libstrat
            ega_sample_data = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(ega_study, ega_dataset_acc, ega_libstrat, ega_exp_acc, exp_erx_acc, sample_ers_acc, ega_sample_acc)
            #print "getting egaf mapping for %s"%ega_dataset_acc
            egaf_mapping, egaf_accs = ega.parse_egaf_mapping(ega_dataset_acc)
            # If DCC has submitted EGAF IDs, then use these first to find raw file. Else, search whole dataset.
            if (not exp_erx_acc == '') and (exp_erx_acc in exp_info):
              # logfile.write("exp_erx_acc for dcc_sample %s = %s\n"%(dcc_sample, exp_erx_acc))
              # logfile.write("ega_dataset_acc for dcc_sample %s = %s\n"%(dcc_sample, ega_dataset_acc))
               for ega_run_acc in exp_info[exp_erx_acc]:
                  for run_file in ega_files[ega_run_acc]:
                     egaf_acc = ''
                     fileSize = ''
                     egaf_mismatch = 0
                     md5_checksum = ega_files[ega_run_acc][run_file]['checksum']
                     unencrypted_checksum = ega_files[ega_run_acc][run_file]['unencrypted_checksum']
                     fileType = ega_files[ega_run_acc][run_file]['fileType']
                     for egaFileName in egaf_mapping:
                        if run_file in egaFileName or egaFileName in run_file or run_file.rsplit('.',1)[0] in egaFileName or egaFileName.rsplit('.',1)[0] in run_file:
                           logfile.write("ega_run_acc = %s\texp_erx_acc = %s\n"%(ega_run_acc, exp_erx_acc))
                           logfile.write("run_file = %s\tegaFileName = %s\n"%(run_file, egaFileName))
                           egaf_acc = egaf_mapping[egaFileName]['egaf_acc']
                           logfile.write("egaf_acc in process_metadata = %s\n"%egaf_acc)
                           #print "egaf_acc in process_metadata = %s"%egaf_acc
                           if len(submitted_ega_files) != 0:
                              if egaf_acc not in submitted_ega_files:
                                 logfile.write("EGAF acc %s for sample %s does not exist in list of submitted EGAF IDs %s\n"%(egaf_acc, dcc_sample, submitted_ega_files))
                                 project_logFile.write("EGAF acc %s for sample %s does not exist in list of submitted EGAF IDs %s\n"%(egaf_acc, dcc_sample, submitted_ega_files))
                              # Commented out for some projects (does not apply for pre-Release 22 projects which may not have EGAF IDs and only have EGAD IDs)
                                 egaf_mismatch = 1
                                 ##outFile.write("%s"%dcc_sample_info)
                                 ##outFile.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\n"%donor_sex)
                                 sample_mismatchedEGAF[seqstrat][dcc_sample] = 1
                                 donor_mismatchedEGAF[seqstrat][donor_id] = 1
                              else:
                                 if egaf_acc in datasetFileInfo[ega_dataset_acc]:
                                    fileSize = datasetFileInfo[ega_dataset_acc][egaf_acc]['fileSize']
                                 if fileSize == '':
                                    logfile.write("Cannot get file size for %s\n"%egaf_acc)
                     if egaf_acc == '':
                        logfile.write("[%s] Cannot find EGAF ID for run file %s in EGA Dataset %s\n"%(project, run_file, ega_dataset_acc))
                        fileExists_notEGAF[seqstrat][dcc_sample] = 1
                     if egaf_mismatch == 0 or egaf_mismatch == 1:
                        #logfile.write("in process_metadata, printing to outFile\n")
                        outFile.write("%s"%dcc_sample_info)
                        outFile.write("%s"%ega_sample_data)
                        outFile.write("%s\t%s\t"%(ega_run_acc, ega_analysis_acc))
                        outFile.write("%s\t"%run_file)
                        outFile.write("%s\t"%egaf_acc)
                        outFile.write("%s\t"%md5_checksum)
                        outFile.write("%s\t"%unencrypted_checksum)
                        outFile.write("%s\t"%fileSize)
                        outFile.write("%s\t"%fileType)
                        outFile.write("%s\t"%matched_sample_id)
                        outFile.write("%s\t"%pairedEnd)
                        outFile.write("%s\t"%insertSize)
                        outFile.write("%s\t"%aligned)
                        outFile.write("%s\t"%referenceGenome)
                        outFile.write("%s\n"%donor_sex)
                     if egaf_acc != '':
                        sample_foundEGA[seqstrat][dcc_sample] = 1
                        donor_withEGAData[seqstrat][donor_id] = 1
            elif not ega_analysis_acc == '':
               for analysis_file in ega_files[ega_analysis_acc]:
                  egaf_acc = ''
                  fileSize = ''
                  md5_checksum = ega_files[ega_analysis_acc][analysis_file]['checksum']
                  unencrypted_checksum = ega_files[ega_analysis_acc][analysis_file]['unencrypted_checksum']
                  fileType = ega_files[ega_analysis_acc][analysis_file]['fileType']
                  for egaFileName in egaf_mapping:
                     if analysis_file in egaFileName or egaFileName in analysis_file or analysis_file.rsplit('.',1)[0] in egaFileName or egaFileName.rsplit('.',1)[0] in analysis_file:
                        egaf_acc = egaf_mapping[egaFileName]['egaf_acc']
                        if egaf_acc in datasetFileInfo[ega_dataset_acc]:
                           fileSize = datasetFileInfo[ega_dataset_acc][egaf_acc]['fileSize']
                        if fileSize == '':
                           logfile.write("Cannot get file size for %s\n"%egaf_acc)
                  if egaf_acc == '':
                     logfile.write("[%s] Cannot find EGAF ID for analysis file %s in EGA Dataset %s\n"%(project, analysis_file, ega_dataset_acc))
                     fileExists_notEGAF[seqstrat][dcc_sample] = 1
                  outFile.write("%s"%dcc_sample_info)
                  outFile.write("%s"%ega_sample_data)
                  outFile.write("%s\t%s\t"%('', ega_analysis_acc))
                  outFile.write("%s\t"%analysis_file)
                  outFile.write("%s\t"%egaf_acc)
                  outFile.write("%s\t"%md5_checksum)
                  outFile.write("%s\t"%unencrypted_checksum)
                  outFile.write("%s\t"%fileSize)
                  outFile.write("%s\t"%fileType)
                  outFile.write("%s\t"%matched_sample_id)
                  outFile.write("%s\t"%pairedEnd)
                  outFile.write("%s\t"%insertSize)
                  outFile.write("%s\t"%aligned)
                  outFile.write("%s\t"%referenceGenome)
                  outFile.write("%s\n"%donor_sex)
               sample_foundEGA[seqstrat][dcc_sample] = 1
               donor_withEGAData[seqstrat][donor_id] = 1
            #else:
               #logfile.write("exp_erx_acc is not found in exp_info for dcc_sample %s = %s\n"%(dcc_sample, exp_erx_acc))
      else:
         logfile.write("libstrat %s was not found in dcc_sample %s\n"%(libstrat, dcc_sample))
         outFile.write("%s"%dcc_sample_info)
         outFile.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\t\t\t%s\n"%(matched_sample_id, donor_sex))
         sample_notFoundEGA[seqstrat][dcc_sample] = 1
         donor_withoutEGAData[seqstrat][donor_id] = 1
   else:
      logfile.write("Did not find dcc sample %s anywhere in EGA\n"%dcc_sample)
      outFile.write("%s"%dcc_sample_info)
      outFile.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\t\t\t%s\n"%(matched_sample_id, donor_sex))
      sample_notFoundEGA[seqstrat][dcc_sample] = 1
      donor_withoutEGAData[seqstrat][donor_id] = 1
   return sample_foundEGA, donor_withEGAData, sample_notFoundEGA, donor_withoutEGAData, fileExists_notEGAF, sample_mismatchedEGAF


#connect()
#ega_dataset_test = get_egaDatasetID("EGAF00001589952")
#print "ega_dataset_test = %s"%ega_dataset_test
#exit()

#egaFileInfo = parse_ega_file_info(fileInfo)
#logfile.write("size of egaFileInfo = %s"%len(egaFileInfo))
project_info = parse_ega_dcc_mapping(mappingFileName)
ega_datasets = get_ega_datasets(rootDir)
logfile.write("ega_datasets = %s"%ega_datasets)
logfile.write("EGA samples in ega_samples\n")
dataset_info = EGA.Dataset()
all_ega_runs, all_ega_analyses, all_ega_files, all_exp_info, all_ega_samples, all_ega_sample_info, dataset_info = ega.parse_datasetFile(rootDir, ega_datasets, {}, {}, {}, dataset_info, parseAll=True)
logfile.write("size of all_ega_samples = %s\n"%len(all_ega_samples.get_samples()))
# Preform audit for each ICGC Project
for project in project_info:
   print("PROJECT = %s"%project)
   unique_samples = EGA.DCCSamples()        # dictionary all of unique analyzed samples in project (sequencing strategy keyed by analyzed sample ID)
   accessions = {}                # accession subnmitted in "raw_data_accession" field keyed by analyzed sample ID
   meta_donors = {}                # sequencing strategy keyed by donor

   donors = {}        
   specimens = {}
   samples = {}

   specimen_to_donor = {}
   sample_to_donor = {}
   sample_to_specimen = {}
   specimen_to_type = {}
   pcawg_samples = {}
 
   analyses = {}

   auditFile = "%s/%s_Audit_ICGC28.tsv"%(output_dir, project)
   if os.path.exists(auditFile):
      logfile.write("auditFile %s already exists\n"%auditFile)
      continue
   logfile.write("\nPROJECT: %s\n"%project)
   outFile = open(auditFile, "w")
   project_logFile = open("%s/%s.log"%(output_dir, project), "w")
   projectDir = "%s/%s"%(submissions_dir, project)
   #ega_study = project_info[project]['ega_study']
   #ega_datasets = project_info[project]['ega_datasets']
   
   sample_foundEGA = {}
   donor_withEGAData = {}
   sample_notFoundEGA = {}
   donor_withoutEGAData = {}
   fileExists_notEGAF = {}
   sample_mismatchedEGAF = {}
   donor_mismatchedEGAF = {}
    
   donorComplete = {}
   donorIncomplete = {}
   donorIssue = {}
   
   # Get all data from project's ICGC submission
   donorFiles = dcc.getFiles('^donor(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   specimenFiles = dcc.getFiles('^specimen(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   sampleFiles = dcc.getFiles('^sample(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   ssmFiles = dcc.getFiles('^ssm_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   cnsmFiles = dcc.getFiles('^cnsm_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   stsmFiles = dcc.getFiles('^stsm_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   sgvFiles = dcc.getFiles('^sgv_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   mirnaSeqFiles = dcc.getFiles('^mirna_seq_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   expSeqFiles = dcc.getFiles('^exp_seq_m(\.[a-zA-Z0-9]+)?\.txt(?:\\.gz|\.bz2)?$', projectDir)
   #expArrayFiles = dcc.getFiles('^exp_array_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   jcnFiles = dcc.getFiles('^jcn_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   #pexpFiles = dcc.getFiles('^pexp_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   methSeqFiles = dcc.getFiles('^meth_seq_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   #methArrayFiles = dcc.getFiles('^meth_array_m(\.[a-zA-Z0-9]+)?\.txt(?:\.gz|\.bz2)?$', projectDir)
   project_logFile.write("sgv files = %s\n"%sgvFiles)
   # Parse clinical files and get donor/specimen/sample information 
   donors = dcc.parse_donors(donorFiles, projectDir)
   specimen_to_donor, specimen_to_type = dcc.parse_specimens(specimenFiles, projectDir)
   sample_to_specimen, sample_to_donor, pcawg_samples = dcc.parse_samples(sampleFiles, specimen_to_donor, projectDir)
   #print "PCAWG samples = %s"%(len(pcawg_samples))
   project_logFile.write("[%s]\tTotal number of donors = %s\n"%(project, len(donors)))
   project_logFile.write("[%s]\tTotal number of specimens = %s\n"%(project, len(specimen_to_donor)))
   project_logFile.write("[%s]\tTotal number of samples = %s\n"%(project, len(sample_to_specimen)))
   
   # Get unique sample_ids across all the different analysis types
   unique_samples, analyses = dcc.get_uniqueSamples(ssmFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(cnsmFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(stsmFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(sgvFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(jcnFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(mirnaSeqFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(expSeqFiles, unique_samples, projectDir, sample_to_donor, analyses)
   #unique_samples, analyses = dcc.get_uniqueSamples(expArrayFiles, unique_samples, projectDir, sample_to_donor, analyses)
   #unique_samples, analyses = dcc.get_uniqueSamples(pexpFiles, unique_samples, projectDir, sample_to_donor, analyses)
   unique_samples, analyses = dcc.get_uniqueSamples(methSeqFiles, unique_samples, projectDir, sample_to_donor, analyses)
   #unique_samples, analyses = dcc.get_uniqueSamples(methArrayFiles, unique_samples, projectDir, sample_to_donor, analyses)
   

   # Print audit to output file
   outFile.write("ICGC DCC Project Code\tICGC Submitted Donor ID\tICGC Submitted Specimen ID\tICGC Submitted Specimen Type\tICGC Submitted Sample ID\tICGC Submitted Sequencing Strategy\tICGC Raw Data Accession\tDCC Submission File\tPCAWG sample\tEGA Study Accession\tEGA Dataset Accession\tEGA Library Strategy\tEGA Experiment Accession\tEGA ERX Accession\tEGA ERS Accession\tEGA Sample Accession\tEGA Run Accession\tEGA Analysis Accession\tEGA Raw Sequence Filename\tEGA File Accession\tMD5 Checksum\tUnencrypted Checksum\tFile Size\tFile Type\tMatched DCC Sample ID\tPaired-End\tInsert Size\tAligned\tReference Genome\tDonor Gender\n") 

   project_logFile.write("[%s]\tTotal number of samples in %s = %s\n"%(project, project, len(unique_samples.unique_dcc_samples)))
   # if analyzed sample ID exists in EGA data, print out corresponding EGA data
   for seqstrat in unique_samples:
      sample_foundEGA[seqstrat] = {}
      sample_notFoundEGA[seqstrat] = {}
      donor_withEGAData[seqstrat] = {}
      donor_withoutEGAData[seqstrat] = {}
      fileExists_notEGAF[seqstrat] = {}
      sample_mismatchedEGAF[seqstrat] = {}
      donor_mismatchedEGAF[seqstrat] = {}


      total_donors = {}
      project_logFile.write("[%s]\tSequencing Strategy = %s\n"%(project, seqstrat))
      project_logFile.write("[%s]\tTotal number of unique samples = %s\n"%(project, len(unique_samples[seqstrat])))
      for dcc_sample in unique_samples[seqstrat]:
         logfile.write("DCC SAMPLE %s\n"%dcc_sample)
         print("DCC SAMPLE %s"%dcc_sample)
         matched_sample_id = ""
         if seqstrat not in total_donors:
            total_donors[seqstrat] = []
         donor_id = sample_to_donor[dcc_sample] # donor_id
         donor_sex = donors[donor_id]
         if donor_id not in total_donors[seqstrat]:
            total_donors[seqstrat].append(donor_id)
         specimen_id = sample_to_specimen[dcc_sample] #specimen_id
         specimen_type = specimen_types[specimen_to_type[sample_to_specimen[dcc_sample]]] #specimen_type
         logfile.write("specimen_type for %s = %s\n"%(dcc_sample, specimen_type))
         # matched_normal_sample
         if dcc_sample in analyses[seqstrat][donor_id]:
            matched_sample_id = analyses[seqstrat][donor_id][dcc_sample]
         #ega_study = project_info[project]['ega_study'] # ega study acc
         raw_data_accession = unique_samples[seqstrat][dcc_sample]['submitted_accession'] #raw data accession
         submitted_ega_datasets = unique_samples[seqstrat][dcc_sample]['ega_accessions']['ega_datasets']
         #logfile.write("submitted_ega_datasets for dcc_sample %s = %s\n"%(dcc_sample, submitted_ega_datasets))
         submitted_ega_runs = unique_samples[seqstrat][dcc_sample]['ega_accessions']['ega_runs']
         submitted_ega_analyses = unique_samples[seqstrat][dcc_sample]['ega_accessions']['ega_analyses']
         submitted_ega_experiments = unique_samples[seqstrat][dcc_sample]['ega_accessions']['ega_experiments']
         submitted_ega_samples = unique_samples[seqstrat][dcc_sample]['ega_accessions']['ega_samples']
         submitted_ega_files = unique_samples[seqstrat][dcc_sample]['ega_accessions']['ega_files']
         dcc_fileName = unique_samples[seqstrat][dcc_sample]['fileName']
         is_pcawg_sample = ""
         if (dcc_sample in pcawg_samples):
            is_pcawg_sample = 'PCAWG'
         dcc_sample_info = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(project, donor_id, specimen_id, specimen_type, dcc_sample, seqstrat, raw_data_accession, dcc_fileName, is_pcawg_sample)
         if len(submitted_ega_files) != 0:
            project_logFile.write("DCC Sample = %s has multiple EGAFs %s\n"%(dcc_sample, submitted_ega_files))
            processed_datasets = {}
            sample_ega_datasets = {}
            for egaf in submitted_ega_files:
               #print "egaf = %s"%egaf
               # get EGA dataset ID from DB
               egad = get_egaDatasetID(egaf)
               print("EGA dataset from DB = %s"%egad)
               if egad not in datasetFileInfo:
                  datasetFileInfo[egad] = getFileSizes(egad)
               if egad not in sample_ega_datasets:
                  sample_ega_datasets[egad] = egad
               elif egad is None and len(submitted_ega_datasets) != 0:
                  project_logFile.write("egad for EGAF %s is none\n"%egaf)
                  for submitted_ega_dataset in submitted_ega_datasets:
                     sample_ega_datasets[submitted_ega_dataset] = submitted_ega_dataset
            for ega_dataset in sample_ega_datasets: 
               if ega_dataset is None or ega_dataset not in ega_datasets:
                  project_logFile.write("ega_dataset is None\n")
                  outFile.write("%s"%dcc_sample_info)
                  outFile.write("\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\t\t\t%s\n"%(matched_sample_id, donor_sex))
                  sample_notFoundEGA[seqstrat][dcc_sample] = 1
                  donor_withoutEGAData[seqstrat][donor_id] = 1
                  continue
               if ega_dataset in processed_datasets:
                  sample_foundEGA, donor_withEGAData, sample_notFoundEGA, donor_withoutEGAData, fileExists_notEGAF, sample_mismatchedEGAF = process_metadata(dataset_info[ega_dataset]['ega_samples'], dcc_sample, seqstrat, dataset_info[ega_dataset]['ega_sample_info'], dataset_info[ega_dataset]['ega_experiments'], dataset_info[ega_dataset]['ega_files'], dcc_sample_info, submitted_ega_files, ega_dataset, matched_sample_id)
                  continue
               else:
                  # get mapping information from Sample_File.map in EGA dataset 
                  egaf_mapping, egaf_accs = ega.parse_egaf_mapping(ega_dataset)
                  processed_datasets[ega_dataset] = dcc_sample
                  if ega_dataset in dataset_info.get_datasets():
                     project_logFile.write("ega_dataset in dataset_info = %s\n"%ega_dataset)
                     project_logFile.write("samples in %s: %s\n"%(ega_dataset, dataset_info[ega_dataset]['ega_samples']))
                     if dcc_sample in dataset_info[ega_dataset]['ega_samples']:
                        project_logFile.write("dcc_sample %s in dataset_info ega_samples\n"%dcc_sample)
                        # get rest of EGA info
                        sample_foundEGA, donor_withEGAData, sample_notFoundEGA, donor_withoutEGAData, fileExists_notEGAF, sample_mismatchedEGAF = process_metadata(dataset_info[ega_dataset]['ega_samples'], dcc_sample, seqstrat, dataset_info[ega_dataset]['ega_sample_info'], dataset_info[ega_dataset]['ega_experiments'], dataset_info[ega_dataset]['ega_files'], dcc_sample_info, submitted_ega_files, ega_dataset, matched_sample_id)
                  else:
                     project_logFile.write("ega_dataset '%s' does not exist in dataset list\n"%(ega_dataset))
         elif len(submitted_ega_datasets) != 0:
            logfile.write("submitted_ega_datasets = %s\n"%submitted_ega_datasets)
            submitted_dataset_info = EGA.Dataset()
            ega_runs, ega_analyses, ega_files, exp_info, ega_samples, ega_sample_info, submitted_dataset_info = ega.parse_datasetFile(rootDir, submitted_ega_datasets, submitted_ega_runs, submitted_ega_experiments, submitted_ega_samples, submitted_dataset_info)
            #logfile.write("EGA Dataset in dataset_info = %s\n"%submitted_dataset_info.get_datasets())
            for ega_dataset in submitted_dataset_info.get_datasets():
               if ega_dataset not in datasetFileInfo:
                  datasetFileInfo[ega_dataset] = getFileSizes(ega_dataset)
               if dcc_sample not in sample_foundEGA[seqstrat]:
                  sample_foundEGA, donor_withEGAData, sample_notFoundEGA, donor_withoutEGAData, fileExists_notEGAF, sample_mismatchedEGAF = process_metadata(submitted_dataset_info[ega_dataset]['ega_samples'], dcc_sample, seqstrat, submitted_dataset_info[ega_dataset]['ega_sample_info'], submitted_dataset_info[ega_dataset]['ega_experiments'], submitted_dataset_info[ega_dataset]['ega_files'], dcc_sample_info, submitted_ega_files, ega_dataset, matched_sample_id)
               #sample_foundEGA, donor_withEGAData, sample_notFoundEGA, donor_withoutEGAData, fileExists_notEGAF = process_metadata(ega_samples, dcc_sample, seqstrat, ega_sample_info, exp_info, ega_files, dcc_sample_info, submitted_ega_files, ega_dataset)
         else:
            # find EGA info in other EGA datasets
            project_logFile.write("Finding sample %s in other EGA datasets\n"%dcc_sample)
            sample_foundEGA, donor_withEGAData, sample_notFoundEGA, donor_withoutEGAData, fileExists_notEGAF, sample_mismatchedEGAF = process_metadata(all_ega_samples, dcc_sample, seqstrat, all_ega_sample_info, all_exp_info, all_ega_files, dcc_sample_info, submitted_ega_files, ega_dataset, matched_sample_id)
      #print "dcc_sample = %s\tseqstrat = %s\tEGAN Sample ID = %s"%(dcc_sample, seqstrat, submitted_ega_samples)
   
      project_logFile.write("Project Summary for %s samples in project %s\n"%(seqstrat, project))
      project_logFile.write("Total # of analyzed samples = %s\n"%len(unique_samples[seqstrat]))
      project_logFile.write("Number of %s samples that exist in EGA = %s/%s\n"%(seqstrat, len(sample_foundEGA[seqstrat]), len(unique_samples[seqstrat])))
      project_logFile.write("Number of %s samples not found in EGA = %s/%s\n"%(seqstrat, len(sample_notFoundEGA[seqstrat]), len(unique_samples[seqstrat])))
      project_logFile.write("Number of %s samples with existing raw file but no EGAF ID = %s/%s\n"%(seqstrat, len(fileExists_notEGAF[seqstrat]), len(unique_samples[seqstrat])))
      project_logFile.write("Number of %s samples with mismatched EGAF IDs = %s/%s\n"%(seqstrat, len(sample_mismatchedEGAF[seqstrat]), len(unique_samples[seqstrat])))
      
      logfile.write("Project Summary for %s samples in project %s\n"%(seqstrat, project))
      logfile.write("Total # of analyzed samples = %s\n"%len(unique_samples[seqstrat]))
      logfile.write("Number of %s samples that exist in EGA = %s/%s\n"%(seqstrat, len(sample_foundEGA[seqstrat]), len(unique_samples[seqstrat])))
      logfile.write("Number of %s samples not found in EGA = %s/%s\n"%(seqstrat, len(sample_notFoundEGA[seqstrat]), len(unique_samples[seqstrat])))
      logfile.write("Number of %s samples with existing raw file but no EGAF ID = %s/%s\n"%(seqstrat, len(fileExists_notEGAF[seqstrat]), len(unique_samples[seqstrat])))
      logfile.write("Number of %s samples with mismatched EGAF IDs = %s/%s\n"%(seqstrat, len(sample_mismatchedEGAF[seqstrat]), len(unique_samples[seqstrat])))
      
      print("Project Summary for %s samples in project %s"%(seqstrat, project))
      print("Total # of analyzed samples = %s"%len(unique_samples[seqstrat]))
      print("Number of %s samples that exist in EGA = %s/%s"%(seqstrat, len(sample_foundEGA[seqstrat]), len(unique_samples[seqstrat])))
      print("Number of %s samples not found in EGA = %s/%s"%(seqstrat, len(sample_notFoundEGA[seqstrat]), len(unique_samples[seqstrat])))
      print("Number of %s samples with existing raw file but no EGAF ID = %s/%s"%(seqstrat, len(fileExists_notEGAF[seqstrat]), len(unique_samples[seqstrat])))
      print("Number of %s samples with mismatched EGAF IDs = %s/%s"%(seqstrat, len(sample_mismatchedEGAF[seqstrat]), len(unique_samples[seqstrat])))
      samplefound_exists = 0
      exclude_pcawg_samples = {}
      total_pcawg_samples = {}
      existsInPCAWGDataset = {}

      project_logFile.write("PCAWG samples:\n")
      for sample in pcawg_samples:
         project_logFile.write("%s\n"%sample)
      project_logFile.write("\n")

  
      test.write("Analyzed %s samples:\n"%seqstrat)
      for sample in unique_samples[seqstrat]:
         test.write("%s\n"%sample)
         if sample in pcawg_samples:
            total_pcawg_samples[sample] = sample
      test.write("\n\n")
      print("Total number of PCAWG samples = %s"%(len(total_pcawg_samples)))
 
      test.write("Not found %s samples:\n"%seqstrat)
      for sample in sample_notFoundEGA[seqstrat]:
         test.write("%s\n"%sample)
         if seqstrat == 'WGS' or seqstrat == 'RNA-Seq':
            if sample in pcawg_samples:
               existsInPCAWGDataset[sample] = sample
      test.write("\n\n")
      print("Total number of PCAWG samples that only exist in PCAWG EGA datasets = %s"%(len(existsInPCAWGDataset)))
      
      print("Number of PCAWG samples = %s\n"%(len(pcawg_samples)))
      test.write("Found %s samples:\n"%seqstrat)
      for sample in sample_foundEGA[seqstrat]:
         test.write("%s\n"%sample)
         if sample in pcawg_samples:
            exclude_pcawg_samples[sample] = sample
      test.write("\n\n")
  
      test.write("%s Samples with raw file but no EGAF ID:\n"%seqstrat)
      for sample in fileExists_notEGAF[seqstrat]:
         if sample in sample_foundEGA[seqstrat]:
            samplefound_exists += 1
            test.write("%s: Has other EGAF IDs\n"%sample)
            if sample in pcawg_samples:
               exclude_pcawg_samples[sample] = sample
         else:
            test.write("%s\n"%sample)
      test.write("\n\n")
      print("Total number of PCAWG samples which are found in non-PCAWG EGA datasets = %s"%(len(exclude_pcawg_samples)))

      
      project_logFile.write("%s"%project)
      project_logFile.write("\t%s"%project_info[project]['country'])
      project_logFile.write("\t%s"%project_info[project]['name'])
      project_logFile.write("\t%s"%seqstrat)
      project_logFile.write("\t%s"%len(unique_samples[seqstrat]))
      project_logFile.write("\t%s"%(len(sample_foundEGA[seqstrat])))
      project_logFile.write("\t%s"%(len(sample_notFoundEGA[seqstrat])))
      project_logFile.write("\t%s"%(len(donor_withEGAData[seqstrat])))
      project_logFile.write("\t%s"%(len(donor_withoutEGAData[seqstrat])))
      project_logFile.write("\t%s"%(len(total_donors[seqstrat])))
      project_logFile.write("\t%s"%(len(fileExists_notEGAF[seqstrat])))
      project_logFile.write("\t%s"%(len(exclude_pcawg_samples)))
      project_logFile.write("\t%s"%(len(existsInPCAWGDataset)))
      project_logFile.write("\t%s"%(len(total_pcawg_samples)))
      project_logFile.write("\t%s\n"%(samplefound_exists))
      
      summary_file.write("%s"%project)
      summary_file.write("\t%s"%project_info[project]['country'])
      summary_file.write("\t%s"%project_info[project]['name'])
      summary_file.write("\t%s"%seqstrat)
      summary_file.write("\t%s"%len(unique_samples[seqstrat]))
      summary_file.write("\t%s"%(len(sample_foundEGA[seqstrat])))
      summary_file.write("\t%s"%(len(sample_notFoundEGA[seqstrat])))
      summary_file.write("\t%s"%(len(donor_withEGAData[seqstrat])))
      summary_file.write("\t%s"%(len(donor_withoutEGAData[seqstrat])))
      summary_file.write("\t%s"%(len(total_donors[seqstrat])))
      summary_file.write("\t%s"%(len(fileExists_notEGAF[seqstrat])))
      summary_file.write("\t%s"%(len(exclude_pcawg_samples)))
      summary_file.write("\t%s"%(len(existsInPCAWGDataset)))
      summary_file.write("\t%s"%(len(total_pcawg_samples)))
      summary_file.write("\t%s\n"%(samplefound_exists))



      for seqstrat in analyses:
         print("Number of %s analyses = %s"%(seqstrat, len(analyses[seqstrat])))
         if seqstrat not in donorComplete:
            donorComplete[seqstrat] = {}
         if seqstrat not in donorIncomplete:
            donorIncomplete[seqstrat] = {}
         if seqstrat not in donorIssue:
            donorIssue[seqstrat] = {}
         for donor_id in analyses[seqstrat]:
            for analyzed_sample in analyses[seqstrat][donor_id]:
               matched_sample = analyses[seqstrat][donor_id][analyzed_sample]
               if seqstrat in sample_foundEGA:
                  if matched_sample is not None:
                     if analyzed_sample in sample_foundEGA[seqstrat] and matched_sample in sample_foundEGA[seqstrat]:
                        donorComplete[seqstrat][donor_id] = 1
                     else:
                        donorIncomplete[seqstrat][donor_id] = 1
                  else:
                     if analyzed_sample in sample_foundEGA[seqstrat]:
                        donorComplete[seqstrat][donor_id] = 1
                     else:
                        donorIncomplete[seqstrat][donor_id] = 1
               else:
                  if seqstrat in fileExists_notEGAF:
                     if analyses_sample in fileExists_notEGAF[seqstrat] or matched_sample in fileExists_notEGAF[seqstrat]:
                        donorIssue[seqstrat][donor_id] = 1
                  elif seqstrat in sample_notFoundEGA:
                     if analyzed_sample_id in sample_notFoundEGA[seqstrat] or matched_sample in sample_notFoundEGA[seqstrat]:
                        donorIncomplete[seqstrat][donor_id] = 1 


         
      for seqstrat in donorComplete:
         project_logFile.write("Total %s donors complete = %s\n"%(seqstrat, len(donorComplete[seqstrat])))
         print("Total %s donors complete = %s"%(seqstrat, len(donorComplete[seqstrat])))
      for seqstrat in donorIncomplete:
         project_logFile.write("Total %s donors incomplete = %s\n"%(seqstrat, len(donorIncomplete[seqstrat])))
         print("Total %s donors incomplete = %s"%(seqstrat, len(donorIncomplete[seqstrat])))
      for seqstrat in donorIssue:
         project_logFile.write("Total %s donors issue = %s\n"%(seqstrat, len(donorIssue[seqstrat])))
         print("Total %s donors issue = %s"%(seqstrat, len(donorIssue[seqstrat])))


   seqstrat_list = ['WGS', 'WXS', 'RNA-Seq']
   ##donor_summary.write("Project\tCompleted WGS Donors\tIncomplete WGS Donors\tPartially Complete WGS Donors\tCompleted WXS Donors\tIncomplete WXS Donors\tPartially Complete WXS Donors\tCompleted RNA-Seq Donors\tIncomplete RNA-Seq Donors\tPartially Complete RNA-Seq Donors\n")
   donor_summary.write("%s"%(project))
   for ss in seqstrat_list:
      if ss in donorComplete:
         donor_summary.write("\t%s"%len(donorComplete[ss]))
      else:
         donor_summary.write("\tNone")
      if ss in donorIncomplete:
         donor_summary.write("\t%s"%len(donorIncomplete[ss]))
      else:
         donor_summary.write("\tNone")
      if ss in donorIssue:
         donor_summary.write("\t%s"%len(donorIssue[ss]))
      else:
         donor_summary.write("\tNone")
   donor_summary.write("\n")
          

summary_file.close()
 


