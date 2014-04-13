import Image 

im = Image.open("oxygen.png")

print im.format, im.size, im.mode

answer = ""
for x in range (0,629,7):
	pix = im.getpixel((x,47))
	(r,g,b,a) = pix
	print pix 
	answer += chr(r)

	if r != g and g != b and r != b:
		print "stops being gray at" + str(x)
		break 

print "answer is: " + answer

anslist = [105, 110, 116, 101, 103, 114, 105, 116, 121]
final = ""
for num in anslist: 
	final += chr(num)

print final
