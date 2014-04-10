#regex? 

import re 

fo = open("prob3text.html","r")

text = fo.read()


match = re.findall("[^A-Z][A-Z]{3}[a-z][A-Z]{3}[^A-Z]",text)

print match

"""
if match:
	print match.group(0)
else:
	print "no match found"
"""
fo.close()