import pprint
import cPickle as pickle 

pkl_file = open('banner.pkl', 'rb')

data1 = pickle.load(pkl_file)



for line in data1:
	linestr = ""
	for tupl in line:
		(char, num) = tupl
		for x in range(num):
			linestr += char 
	
	print linestr


#pprint.pprint(data1)

pkl_file.close()

