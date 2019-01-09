
import sys
import re

with open(sys.argv[1], 'r') as inf, open(sys.argv[2], 'w') as out:
	keep = []
	count = 0
	transcript = 'ENST00000000001'
	for line in inf:
		count += 1
		start = line.split('\t')[1]
		end = line.split('\t')[2]

		query = start + "\t" + end
		novel = True 

		for k in keep:
			if query in k:
				novel = False

		if novel:
			# print re.sub(r'ENST[0-9]{11}', transcript, line)
			keep.append(re.sub(r'ENST[0-9]{11}', transcript, line))

	for l in keep:
		out.write(l)


	
