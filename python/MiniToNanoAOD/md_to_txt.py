

f = open('datasets.md', 'r')


data_2016_file = open('my_dataset_lists/2016/data_2016_datasets.txt', 'w')
MC_2016_file = open('my_dataset_lists/2016/MC_2016_datasets.txt', 'w')

data_2017_file = open('my_dataset_lists/2017/data_2017_datasets.txt', 'w')
MC_2017_file = open('my_dataset_lists/2017/MC_2017_datasets.txt', 'w')

data_2018_file = open('my_dataset_lists/2018/data_2018_datasets.txt', 'w')
MC_2018_file = open('my_dataset_lists/2018/MC_2018_datasets.txt', 'w')

data_files = [data_2016_file, data_2017_file, data_2018_file]
MC_files = [MC_2016_file, MC_2017_file, MC_2018_file]

read = f.readlines()

for line in read:
  data_strings = ["Run2016", "Run2017", "Run2018"]
  MC_strings = ["RunIISummer16MiniAODv3", "RunIIFall17MiniAODv2", "RunIIAutumn18MiniAOD"]

  MC_string = "RunIISummer16MiniAODv3" #2016 MC
  data_2016_string = "Run2016" #2016 Data

  for year in [2016, 2017, 2018]:
    dataset_string = [x for x in line.split("|") if data_strings[year-2016] in x and "https" not in x]
    MCset_string = [x for x in line.split("|") if MC_strings[year-2016] in x and "https" not in x]
    if dataset_string:
      if data_strings[year-2016] == dataset_string[0].strip(): continue
      data_files[year-2016].write(dataset_string[0].strip()+"\n")
    if MCset_string:
      if "<sup>" in MCset_string[0]:
        MCset_string = [x for x in MCset_string[0].split("<") if MC_strings[year-2016] in x]
        MCset_string = [x for x in MCset_string[0].split(">") if MC_strings[year-2016] in x]
      if MC_strings[year-2016] == MCset_string[0].strip(): continue
      MC_files[year-2016].write(MCset_string[0].strip()+"\n")

