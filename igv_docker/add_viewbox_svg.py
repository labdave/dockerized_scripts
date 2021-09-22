from PIL import Image
import os
import sys

img_dir = sys.argv[1]
out_dir = sys.argv[2]

svgs = []
for filename in os.listdir(img_dir):
	if filename.endswith(".svg"):
		# print(os.path.abspath(os.path.join("outputs/", filename)))
		prefix = filename.rstrip(".svg")
		if not os.path.isfile(os.path.join(img_dir, prefix+".png")):
			print(prefix+".png not found!")
		else:
			svgs.append(filename)
			# print(filename)

for svg in svgs:
	png = svg.replace("svg", "png")
	img = Image.open(os.path.join(img_dir, png))
	height = img.height
	width = img.width

	lines = []
	with open(os.path.join(img_dir, svg), "r") as f:
		for line in f:
			if "viewBox" in line:
				continue
			line = line.replace("xmlns=\"http://www.w3.org/2000/svg\"", f"xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width} {height}\"")
			line = line.replace("xmlns=\"http://www.w3.org/2000/svg\"", f"xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\"")
			lines.append(line)
	with open(os.path.join(out_dir, svg), "w+") as f:
		f.write("".join(lines))