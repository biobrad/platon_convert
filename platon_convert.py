#!/bin/python
import os
import pandas as pd
import json
pd.options.display.width = 0
dirname=os.path.dirname(os.path.abspath(__file__))
#
#Gets all filenames in platon folder that has returned a result
#
ext=('.json')
filelist = []
filelist2 = []
for file in os.listdir(dirname):
 if file.endswith(ext):
  filelist.append(file)
for i in filelist:
 filelist2.append(i.split('.')[0])
#
#Gets results for plasmids that have returned a result with AMR presence
#
for i in filelist2:
 tsv=pd.read_table(i + '.tsv')
 frame=pd.DataFrame(data=tsv)
 amrframe=frame[frame['# AMRs']!=0]
 amrframe2=amrframe[amrframe['# Plasmid Hits']!=0]
 dftsv=pd.DataFrame(amrframe2, columns=['Genome', 'ID', '# Plasmid Hits', '# AMRs'])
 dftsv['Genome'] = i
#
#creates individual csv for these results which could be appended into a total ie: Genome: xyz Plasmids: 8 AMRs: 7
#
 dftsv.to_csv(i + "_amrPlasmids.csv", index=False)
#
#This loop uses the above .csv files and contigs to interogate the .json files for the plasmid details
#
full=pd.DataFrame({'Genome' : '', 'id' : '', 'contig_start' : '', 'contig_end' : '', 'AMRs' : '', 'plasmid.id' : '', 'plasmid_start' : '', 'plasmid_end' : '', 'coverage' : '', 'identity' : '', 'plasmid.length' : '','AMRs' : '' }, index=[0])
for x in filelist2:
 f=open(x + ".json")
 d=json.load(f)
 df=pd.read_csv(x + "_amrPlasmids.csv")
 tmpdf=pd.DataFrame({'Genome' : '', 'id' : '', 'contig_start' : '', 'contig_end' : '', 'plasmid.id' : '', 'plasmid_start' : '', 'plasmid_end' : '', 'coverage' : '', 'identity' : '', 'plasmid.length' : ''}, index=[0])
 for y in range(len(df)):
  data=[]
			#extract contig details
  getcontig=pd.json_normalize(data=d[df.loc[y, "ID"]])
			#extract plasmid details
  getplasmid=pd.json_normalize(data=d[df.loc[y, "ID"]], record_path='plasmid_hits', meta=['plasmid'], errors='ignore')
			#extract AMR details and add to list
  getamr=pd.json_normalize(data=d[df.loc[y, "ID"]], record_path='amr_hits')
			#adds the genome name to each line
  data.extend(getamr['gene'].tolist())
			#concatentates the contig ID with the plasmid details
  ccconplas=pd.concat([getcontig['id'], getplasmid], axis=1)
			#Concatenates the contig details with the plasmid details
  ccplasinfo=ccconplas[['id', 'contig_start', 'contig_end', 'plasmid.id', 'plasmid_start', 'plasmid_end', 'coverage', 'identity', 'plasmid.length']]
			#creates the final concatenated data from from json extraction (good print() point for testing purposes)
  finframe=pd.DataFrame(ccplasinfo)
			#next two lines insert the sequence identifier into each row of the dataframe and adds the AMR information from the list
  finframe.insert(0, 'Genome', x)
  finframe.insert(4, 'AMRs', [data])
			#probably unnecessary rerender of the dataframe
  ccfin=pd.DataFrame(finframe)
			#appending concatenated details to the dataframe outside of the second loop
  tmpdf=tmpdf.append(ccfin, ignore_index=True)
		#dropping the axis details
 tmpdf=tmpdf.drop([0])
		# appending the looped data to the frame outside of the first loop
 full=full.append(tmpdf, ignore_index=True)
	#dropping the index again
full=full.drop([0])
#
#writes the complete .csv to file
#
full.to_csv('AllGenomesPlasmidsWithAMRs.csv', index=False)
