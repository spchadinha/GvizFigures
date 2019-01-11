
import sys
from collections import defaultdict

# need code file and dataframe
template = sys.argv[1]
data = sys.argv[2]
output = template.replace(".R","_params.R")
# map each biotype to a color, represented by a position in an array of colors
biotype2color = defaultdict(int)
# map each transcript to a biotype
transcript2biotype = defaultdict(str)
# nonredundant set of biotypes, needed to count total number 
biotypes = set()

with open(data,'r') as dat:
	# skip the header
	a = dat.readline()
	print a
	# map transcript ID to biotype and count biotypes
	for line in dat.readlines():
		d = line.rstrip().split("\t")
		transcript = d[6]
		biotype = d[8]
		transcript2biotype[transcript] = biotype
		biotypes.add(biotype)

# map a color to a biotype
color = 1
for bio in biotypes:
	biotype2color[bio] = color
	color += 1

cmd = ""
for transcript in transcript2biotype.keys():
	# color corresponding to the biotype of the transcript
	col = biotype2color[transcript2biotype[transcript]]
	cmd += "\t{} = transcript_cols[{}],\n".format(transcript, col)
print transcript2biotype
with open(template,'r') as temp, open(output,'w') as out:
	t = temp.read()
	o = t.replace("#insert#\n",cmd)
	out.write(o)