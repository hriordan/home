import urllib 

defurl = "http://www.pythonchallenge.com/pc/def/linkedlist.php?nothing=" #default url
url = defurl + "8022" #add midway node

x = 1 

while x <= 200:
	site = urllib.urlopen(url)
	text = site.read()
	text.strip()
	thing = text.split()
	print text
	url = defurl + thing[len(thing)-1]
	#print url
	x += 1


