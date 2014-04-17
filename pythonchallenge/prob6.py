import zipfile 

file = zipfile.ZipFile("./channel.zip", "r")

name = "90052.txt"


commentstring = ""



while True: 
	info = file.getinfo(name)
	text = file.read(name)
	text.strip()
	#print text

	thing = text.split()	 

	if unicode(thing[len(thing) - 1]).isnumeric() == False or thing[0] != 'Next':
		print "breaking"
		print thing
		break

	name = thing[len(thing) - 1]
	name += ".txt"

	commentstring += info.comment

print commentstring

