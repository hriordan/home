import Image, ImageDraw

def coordsList(file):
	f = open(file, "r")
	coords = f.read()
	coords = coords.replace("\n","")
	coords = coords.replace("\r","")
	coords = coords.split(",")
	coords = map(int, coords)
	return coords

im = Image.open("dots.jpg")
draw = ImageDraw.Draw(im)

coords1 = coordsList("coords.txt")
coords2 = coordsList("coords2.txt")

print coords1
print coords2



draw.line(coords1,fill=128)
draw.line(coords2,fill=10)

im.save("dots2.jpg")



#im = Image.new("RGB",(512,512),"white")

