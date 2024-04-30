#!/usr/bin/env python

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import os
import pandas as pd
import argparse



# define the necessary argument (directory name)
parser=argparse.ArgumentParser(description="DESCRIPTION",formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-d", "--dir", help="Folder where all the BED files are stored", default='')
args=parser.parse_args()
# define species name and the list of clusters

species = args.dir[args.dir.find('5.2_')+4:args.dir.find('_results')]
path = os.path.join(os.getcwd(), args.dir)
parent_dir=os.path.abspath(os.path.join(os.path.join(path, os.pardir), os.pardir))
output_dir=os.path.join(parent_dir,'6_'+species+'_plots')
clusters = []

for file in os.listdir(path):
	clusters.append(file[:file.find('.')])

clusters=list(set(clusters))
# column names
genomecov_col = ["chr","chrPos","depth"]
bed_col = ["chrom","chromStart","chromEnd","name","score","strand"]

# empty lists, to be filled and transformed as a single DF
chr=[]
chrPos=[]
depth_total=[]
depth_5=[]
depth_3=[]
centroid=[]

# read each file and fill lists
for cluster in clusters:
	total_genomecov=pd.read_csv(os.path.join(path,cluster)+'.genomecov',sep='\t',names=genomecov_col)
	total_genomecov_5=pd.read_csv(os.path.join(path,cluster)+'.genomecov.5',sep='\t',names=genomecov_col)
	total_genomecov_3=pd.read_csv(os.path.join(path,cluster)+'.genomecov.3',sep='\t',names=genomecov_col)

	total_genomecov = total_genomecov[(total_genomecov[['depth']] != 0).any(axis=1)]
	total_genomecov_5 = total_genomecov_5[(total_genomecov_5[['depth']] != 0).any(axis=1)]
	total_genomecov_3 = total_genomecov_3[(total_genomecov_3[['depth']] != 0).any(axis=1)]

	df_partial = total_genomecov.merge(total_genomecov_5, on = ['chr', 'chrPos'], how = 'left')

	df_temp = df_partial.merge(total_genomecov_3, on = ['chr', 'chrPos'], how = 'left')

	df_temp.columns=['chr','chrPos','depth_total','depth_5', 'depth_3']

	df_temp['centroid']=cluster

	chr += df_temp['chr'].to_list()
	chrPos += df_temp['chrPos'].to_list()
	depth_total += df_temp['depth_total'].to_list()
	depth_5 += df_temp['depth_5'].to_list()
	depth_3 += df_temp['depth_3'].to_list()
	centroid += df_temp['centroid'].to_list()


	bed = pd.read_csv(os.path.join(path,cluster)+'.bed', sep='\t', names=bed_col)
	bed_reads= bed[['name']]
	

# create the final DF
df=pd.DataFrame({'chr': chr, 'chrPos': chrPos, 'depth_total' : depth_total, 'depth_5': depth_5, 'depth_3' : depth_3, 'centroid' : centroid})


# groupby the id
grps = list(df.groupby('centroid'))
# create a start and end step because we have three graphs for each id
start=0
end=3

buttons = []
figs = go.Figure()
# iterate over the groupby object
for idx,v in enumerate(grps):
	# use the centroid as the name of each graph
	name = v[0]
	# set the visibility of each graph on selection
	visible = [False] * len(grps)*3
	visible[start:end] = [True]*3
	# update the start/end values
	start+=3
	end+=3
	# create your graph for each id
	fig = px.line(v[1], x='chrPos', y='depth_total',color_discrete_sequence=["#333131"])
	fig.add_bar(x=v[1]['chrPos'], y=v[1]['depth_5'], name="5' clusters",marker_color='forestgreen')
	fig.add_bar(x=v[1]['chrPos'], y=v[1]['depth_3'], name="3' clusters",marker_color='firebrick')
  
	
	
	# set the visibility for the first three traces as the default
	
	if idx == 0:
		fig.data[0].visible = True
		fig.data[1].visible = True
		fig.data[2].visible = True
	else:
		fig.data[0].visible = False
		fig.data[1].visible = False
		fig.data[2].visible = False
	# add all your traces to the empty go.Figure
	figs.add_traces(fig.data)
	figs.update_xaxes(tick0=0, dtick=2, ticks='inside')
	# append your dropdowns
	buttons.append(dict(label=f'Cluster {name}',
						method="update",
						args=[{"visible": visible},
							  {"title":f"Cluster n°: {name}"}]))


updatemenus = [{'active': 0, "buttons": buttons, 'showactive': True}]
figs.update_layout(title=f'Cluster n°: {grps[0][0]}', title_x=0.5, updatemenus=updatemenus,template='plotly_white')
figs.write_html(os.path.join(output_dir, 'Plot.html'))
print('HTML plot file saved in:', os.path.join(output_dir, 'Plot.html'))

