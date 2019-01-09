
import sys

# need code file and dataframe
template = sys.argv[1]
data = sys.argv[2]
output = template.replace(".R","_params.R")

transcripts = set()
with open(data,'r') as dat:
	# skip the header
	dat.readline()
	# add each ENST to the set of transcripts
	for line in dat.readlines():
		d = line.split("\t")
		transcripts.add(d[6])

transcripts = list(transcripts)
cmd = ""
for i in range(len(transcripts)):
	cmd += "\t{} = transcript_cols[{}],\n".format(transcripts[i], i+1)

with open(template,'r') as temp, open(output,'w') as out:
	t = temp.read()
	o = t.replace("#insert#\n",cmd)
	out.write(o)