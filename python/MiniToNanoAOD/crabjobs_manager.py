import datetime
import os
date = str(datetime.date.today())
submit = True

submit_MC = True
submit_data = False

storage_location = "/store/user/daebi/HHbbWW"	 				#Location of stored output files from CRAB
eos_location = "/eos/user/d/daebi/HHbbWW_crabfiles_March30_2022"		#Location of crab working files (crab status files/logs)
storage_site = "T3_US_FNALLPC" 							#T3_CH_CERNBOX or T3_US_FNALLPC
local_dir_for_cfg_files = "NanoAODproduction_" 					#"txt_"+{year} format
publish = "True"								#T/F on publishing to DAS

radion_and_graviton_2016 = ["GluGluToRadionToHHTo2B2VTo2L2Nu", "GluGluToBulkGravitonToHHTo2B2VTo2L2Nu", "GluGluToRadionToHHTo2B2Tau", "GluGluToBulkGravitonToHHTo2B2Tau", "GluGluToRadionToHHTo2B2WToLNu2J", "GluGluToRadionToHHTo2B2VToLNu2J", "GluGluToBulkGravitonToHHTo2B2WToLNu2J"]
ST_2016 = ["ST"] #["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"] some of these but this is weird ask tao, may need to submit individually
TT_2016 = ["ttHJetToNonbb", "ttHJetTobb", "THQ", "THW", "TTZToLLNuNu", "TTZToLL", "TTZToQQ", "TTWJetsToLNu", "TTWJetsToQQ", "TTWW", "TTWH", "TTZH"]
TT_2016_SPECIAL = ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic"] #DON'T SUBMIT THESE!!! Many jobs related to uncertanties, submit by hand. Checking jobs with this list is okay
DY_2016 = ["DYJetsToLL", "DYToLL"]

to_do_list_2016 = TT_2016_SPECIAL + ST_2016
done_list_2016 = DY_2016

to_do_list = to_do_list_2016
done_list = done_list_2016

def make_crab_files(yearlist):
  cwd = os.getcwd()
  for year in yearlist:
    folder_name = "{local_dir}_{year}/".format(year = year, local_dir = local_dir_for_cfg_files)
    if not os.path.exists(folder_name): os.mkdir(folder_name)

    os.chdir(folder_name)
    if year == 2016:  os.system("cmsDriver.py NanoAODproduction_2016_cfg -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --no_exec --conditions 102X_mcRun2_asymptotic_v8 --era Run2_2016,run2_nanoAOD_94X2016 --nThreads=8")
    if year == 2017:  os.system("cmsDriver.py NanoAODproduction_2017_cfg -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --no_exec --conditions 102X_mc2017_realistic_v8 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --nThreads=8")
    if year == 2018:  os.system("cmsDriver.py NanoAODproduction_2018_cfg -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --no_exec --conditions 102X_upgrade2018_realistic_v21 --era Run2_2018,run2_nanoAOD_102Xv1 --nThreads=8")
    os.chdir("..")

    if submit_MC:
      f = open('my_dataset_lists/{year}/MC_{year}_datasets.txt'.format(year = year), 'r')
      lines = f.readlines()
      for line in lines[:-1]:
        process_name = line.split("_")[0][1:]
        file_title = line.split("/")[1]
        folder = folder_name +"MC/" + process_name
        if not os.path.exists(folder): os.makedirs(folder)
        requestName = "{file_title}".format(year = year, file_title = file_title)
        fileName_crab = folder+'/crab3_'+requestName+'.py'
        fileName_cfg = '../../NanoAODproduction_{year}_cfg_NANO.py'.format(year = year)

        crab_script = open(fileName_crab, "write")
        crab_script.write("from CRABClient.UserUtilities import config\n")
        crab_script.write("config = config()\n")
        crab_script.write("config.General.requestName = '{requestName}'\n".format(requestName = requestName))
        crab_script.write("config.General.workArea = '{eos_location}/{year}/{process_name}'\n".format(year = year, process_name = process_name, eos_location = eos_location))
        crab_script.write("config.General.transferOutputs = True\n")
        crab_script.write("config.General.transferLogs = True\n")
        crab_script.write("config.JobType.pluginName = 'Analysis'\n")
        crab_script.write("config.JobType.psetName = '{fileName_cfg}'\n".format(fileName_cfg = fileName_cfg))
        crab_script.write("config.JobType.maxMemoryMB = 8000\n")
        crab_script.write("config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours\n")
        crab_script.write("config.JobType.numCores = 8\n")
        crab_script.write("config.JobType.allowUndistributedCMSSW = True\n")
        crab_script.write("config.Data.inputDataset = '{line}'\n".format(line = line[:-1]))
        crab_script.write("config.Data.inputDBS = 'global'\n")
        crab_script.write("config.Data.splitting = 'FileBased'\n")
        #crab_script.write("config.Data.splitting = 'Automatic'\n")
        crab_script.write("config.Data.unitsPerJob = 1\n")
        crab_script.write("config.Data.outLFNDirBase = '{storage_location}/HHbbWW_datasets/{year}/{process_name}/'\n".format(year = year, process_name = process_name, storage_location = storage_location))
        crab_script.write("config.Data.publication = {publish}\n".format(publish = publish))
        crab_script.write("config.Data.outputDatasetTag = config.General.requestName\n")
        crab_script.write("config.Site.storageSite = '{storage_site}'\n".format(storage_site = storage_site))
        crab_script.write("config.Site.ignoreGlobalBlacklist = True\n")
    os.chdir(cwd)
        
def submit_jobs(yearlist):
  print "Only submitting those in to do list ", to_do_list
  print "Skipping those in done list ", done_list
  cwd = os.getcwd()
  for year in yearlist:
    os.chdir(cwd+"/{local_dir}_{year}/".format(year = year, local_dir = local_dir_for_cfg_files))
    for datatype in os.listdir("."):
      if ".py" in datatype: continue
      os.chdir(datatype)
      for subdir in os.listdir("."):
        if subdir not in to_do_list:
          print subdir, " not in to do list, skipping for now"
          continue
        os.chdir(subdir)
        for filename in os.listdir("."):
          if 'crab3' in filename:
            print "Submitting ", filename
            os.system("crab submit {filename}".format(filename = filename))
        os.chdir("..")
        print "Finished all ", subdir
      os.chdir("..")
      print "Finished all ", datatype
    os.chdir("..")  
    print "Finished all ", year

def check_jobs(yearlist):
  print "Only checking those in to do list ", to_do_list
  print "Skipping those in done list ", done_list
  cwd = os.getcwd()
  for year in yearlist:
    os.chdir("{eos_location}/{year}/".format(year = year, eos_location = eos_location))
    for datatype in os.listdir("."):
      if datatype not in to_do_list:
        print datatype, " not in to do list, skipping for now"
        continue
      os.chdir(datatype)
      for filename in os.listdir("."):
        print "Checking ", filename
        os.system("crab status {filename}".format(filename = filename))
      os.chdir("..")
      print "Finished all ", datatype
    os.chdir("..")
    print "Finished all ", year
  os.chdir(cwd)



def kill_jobs(yearlist):
  print "Only killing those in to do list ", to_do_list
  print "Skipping those in done list ", done_list
  cwd = os.getcwd()
  for year in yearlist:
    os.chdir("{eos_location}/{year}/".format(year = year, eos_location = eos_location))
    for datatype in os.listdir("."):
      if datatype not in to_do_list:
        print datatype, " not in to do list, skipping for now"
        continue
      os.chdir(datatype)
      for filename in os.listdir("."):
        print "Checking ", filename
        os.system("crab kill {filename}".format(filename = filename))
      os.chdir("..")
      print "Finished all ", datatype
    os.chdir("..")
    print "Finished all ", year
  os.chdir(cwd)


def resubmit_jobs(yearlist):
  print "Only resubmitting those in to do list ", to_do_list
  print "Skipping those in done list ", done_list
  cwd = os.getcwd()
  for year in yearlist:
    os.chdir("{eos_location}/{year}/".format(year = year, eos_location = eos_location))
    for datatype in os.listdir("."):
      if datatype not in to_do_list:
        print datatype, " not in to do list, skipping for now"
      os.chdir(datatype)
      for filename in os.listdir("."):
        print "Reading ", filename
        readout = os.popen("crab status {filename} --long | grep 'Jobs status:                    failed'".format(filename = filename)).read()
        if readout:
          print "Resubmitting failed jobs ", filename
          os.system("crab resubmit {filename}/{filename_two}".format(filename = filename, filename_two = filename_two))
        readout = os.popen("crab status {filename}/{filename_two} --long | grep 'SUBMITFAILED'".format(filename = filename, filename_two = filename_two)).read()
        if readout:
          print "Resubmitting submitfailed ", filename
          os.system("crab resubmit {filename}/{filename_two}".format(filename = filename, filename_two = filename_two))

      os.chdir("../..")
      print "Finished all ", datatype
    os.chdir("..")
    print "Finished all ", year
  os.chdir(cwd)




def help():
  print "To make the crab cfg files: make_crab_files(yearlist)"
  print "To submit the crab jobs   : submit_jobs(yearlist)"
  print "To check the jobs         : check_jobs(yearlist)"
  print "To kill the jobs          : kill_jobs(yearlist)"
  print "To resubmit failed jobs   : resubmit_jobs(yearlist)"
  print "yearlist = [2016, 2017, 2018]"
  print "Current to do list = ", to_do_list
  print "Current done list = ", done_list
  print "help() to repeat"

help()
