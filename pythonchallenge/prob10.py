

def tokenize(nums):
	tokens = []
	temptoken = ""
	curNum = nums[0]

	for i in range(len(nums)):

		if nums[i] != curNum:
			tokens.append(temptoken)
			curNum = nums[i]
			temptoken = ""
			temptoken += nums[i]
		else:
			temptoken += nums[i]

	tokens.append(temptoken)
	return tokens 


def readtokens(tokens):
	new = ""
	for string in tokens:
		new += string[0]
		new += str(len(string))
	return new

def lasseq(seed, n):
 	tokens = tokenize(seed)	#tokenize '1'
 	lis = []
 	lis.append(tokens)
 	x = 0
 	while x <= n:
 		lis.append(readtokens(lis[x]))
 		tokenize(lis[x+1])
 		x = x+1

 	return lis

SEQ = lasseq("1",30)

print SEQ[30]






