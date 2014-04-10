#For problem "1" of the python challenge series

def rotc(char, n):
	ascii = ord(char) + n
	newchar = chr(((ascii - 97) % 26) + 97)
	return newchar   

str = raw_input("Enter text to decode:")
str = list(str)

for x in range(len(str)):
	if str[x].isalpha():
		str[x] = rotc(str[x],2)


str = ''.join(str)
print str