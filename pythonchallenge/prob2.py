fo = open("prob2text","r")
instr = fo.read()
outstr = ""

for c in instr:
	if c.isalpha():
		outstr += c 

fo.close()

print outstr