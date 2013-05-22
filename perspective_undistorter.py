#!/usr/bin/env python
#Copyright (c) 2013 Johannes Reinhardt <jreinhardt@ist-dein-freund.de>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
#of the Software, and to permit persons to whom the Software is furnished to do
#so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import numpy as np
import sys
import argparse
import numpy as np
from numpy.linalg import norm
import xml.dom.minidom as dom
from os.path import splitext
from random import choice

#STL units are assumed to be mm
#SVG units are px, conversion:
px_to_mm = 0.282222225
mm_to_px = 1./px_to_mm

#An inputfile reads a object, and provides information about the bounding box
#of the object and about points and lines between these points
class InputFile:
	def __init__(self):
		self.points = []
		self.lines = []

	def get_max_coords(self):
		raise NotImplementedError

	def get_points(self):
		return np.array(self.points)

	def get_lines(self):
		return np.array(self.lines)

class InputSVG(InputFile):
	def __init__(self,fid,distance,epsilon,h_o):
		InputFile.__init__(self)

		self.distance = distance
		self.epsilon = epsilon
		self.h_o = h_o

		self.dom = dom.parse(fid)
		self.svg = self.dom.getElementsByTagName('svg')[0]

		#determine dimensions
		self.width, self.height = None, None
		if self.svg.getAttribute('width')[-2:] == "mm":
			self.width = float(self.svg.getAttribute('width')[:-2])
		else:
			self.width = float(self.svg.getAttribute('width'))*px_to_mm
		if self.svg.getAttribute('height')[-2:] == "mm":
			self.height = float(self.svg.getAttribute('height')[:-2])
		else:
			self.height = float(self.svg.getAttribute('height'))*px_to_mm

		#get points and lines
		groups = self.svg.getElementsByTagName('g')
		for g in groups:
			paths = self.svg.getElementsByTagName('path')
			for path in paths:
				relative = None
				first_point_idx = len(self.points)
				prev_point_idx = None

				#The current variables are in svg coordinates
				r_hat = np.zeros((2,))
				coordinates = path.getAttribute('d').split()
				#handle closed paths
				for c in coordinates:
					if c == 'm'  or c == 'l':
						relative = True
					elif c == 'M' or c == 'L':
						relative = False
					elif c == 'z':
						relative = True
						self.lines.append([prev_point_idx,first_point_idx])
						break
					elif  c == 'Z':
						relative = False
						self.lines.append([prev_point_idx,first_point_idx])
						break
					else:
						if relative:
							step = np.array(map(float,c.split(',')))
							r_hat += step
						else:
							r_hat = np.array(map(float,c.split(',')))

						#transform to 3D coordinates
						r = np.zeros((3,))
						r[0] = self.distance
						r[1] = 0.5*self.width -r_hat[0]
						r[2] = self.h_o - self.epsilon - r_hat[1]
						self.points.append(r)

						if not prev_point_idx is None:
							#not first point in path
							self.lines.append([prev_point_idx,len(self.points) - 1])
						prev_point_idx = len(self.points) - 1

	def get_max_coords(self):
		return np.array([self.distance,self.width/2,self.h_o - self.epsilon])

class InputSTL(InputFile):
	def __init__(self,fid):
		InputFile.__init__(self)
		self.name = fid.readline().split()[1]
		self.bbox = np.zeros((3,))

		triangles = []

		line = fid.readline()
		while line.split()[0] != "endsolid":
			#skip normal and loop start
			fid.readline()
			idx = []
			for i in range(3):
				parts = fid.readline().split()
				point = map(float,parts[1:])
				fabs_point = np.fabs(np.array(point))
				self.bbox[0] = max(self.bbox[0],point[0])
				self.bbox[1] = max(self.bbox[1],np.fabs(point[1]))
				self.bbox[2] = max(self.bbox[2],point[2])
				if not point in self.points:
					idx.append(len(self.points))
					self.points.append(point)
				else:
					idx.append(self.points.index(point))
			triangles.append(idx)
			#skip loop and facet end
			fid.readline()
			fid.readline()
			line = fid.readline()

		for triangle in triangles:
			for i,j in [(0,1),(0,2),(1,2)]:
				line = [triangle[i],triangle[j]]
				rline = [triangle[j],triangle[i]]

				if line in self.lines or rline in self.lines:
					continue
				self.lines.append(line)
	def get_max_coords(self):
		return self.bbox

def parse_filter_file(fid):
	line = fid.readline()
	lines = {'filtered':[], 'partial':[], 'extra':[]}
	mode = None
	while not line == '':
		line = line.strip()
		if len(line) == 0 or line[0] == '#':
			#empty line or comment
			pass
		elif line[0] == '[':
			#section header line
			mode = line[1:-1]
			if not mode in lines.keys():
				print "Unknown section '%s' encountered when parsing filter file"
		else:
			if mode == 'extra':
				points = map(int,line.split())
				lines[mode].append(points)
			else:
				line = int(line)
				lines[mode].append(line)
		line = fid.readline()
	return lines['filtered'],lines['partial'], lines['extra']
			
			

	return filtered_lines, partial_lines, extra_lines

def filter_and_dictify(points,lines,filtered_lines, extra_lines):
	if len(extra_lines) > 0:
		lines = np.vstack((input.get_lines(), np.array(extra)))
	new_points = {}
	new_lines = {}
	#filter lines
	for i in range(len(lines)):
		if i in filtered_lines:
			continue
		new_lines[i] = lines[i]
		for idx in lines[i]:
			new_points[idx] = points[idx]
	return new_points, new_lines

#An output file allows to output lines and points to a useful representation
class OutputFile:
	def __init__(self):
		self.points = {}
		self.lines = {}
		self.lines_partial = []
	#point is a 2 component np.array with coordinates
	def add_point(self,point,label):
		self.points[label] = point
	#point is a 2 component list with labels referencing points
	def add_line(self,points,label,partial=False):
		self.lines[label] = points
		if partial:
			self.lines_partial.append(label)
	def to_file(self,fid):
		pass

class OutputTable(OutputFile):
	def __init__(self,origin, reference):
		OutputFile.__init__(self)
		self.origin = origin
		self.reference = reference

	def to_file(self,fid):
		#sort points in nearest neighbor order for efficient painting
		l_nn = []
		current = self.origin
		for i in range(len(self.points)):
			cand = None
			for l,p in self.points.iteritems():
				if l in l_nn:
					continue
				if cand is None or norm(current - p) < norm(current - self.points[cand]):
					cand = l
			l_nn.append(cand)
			current = self.points[cand]

		fid.write("Point List\n")
		fid.write("==========\n\n")
		fid.write("ID\td_o \td_r\n")
		d_or = norm(self.reference - self.origin)/1e3
		fid.write("%s\t%2.2f\t%2.2f\n" % ("Or",0,d_or))
		fid.write("%s\t%2.2f\t%2.2f\n" % ("Rf",d_or,0))
		for l in l_nn:
			#in m
			d_o = norm(self.origin - self.points[l])/1e3
			d_r = norm(self.reference - self.points[l])/1e3
			fid.write("%s\t%2.2f\t%2.2f\n" % (l,d_o,d_r))

		fid.write("\nLine List\n")
		fid.write("=========\n\n")
		fid.write("ID\tsrc\tdst\n")
		for l, points in self.lines.iteritems():
			if not l in self.lines_partial:
				fid.write("%s\t%s\t%s\n" % (l, points[0], points[1]))

		fid.write("\nPartial Line List\n")
		fid.write("=========\n\n")
		fid.write("ID\tsrc\tdst\n")
		for l, points in self.lines.iteritems():
			if l in self.lines_partial:
				fid.write("%s\t%s\t%s\n" % (l, points[0], points[1]))


class OutputSVG(OutputFile):
	def __init__(self,max_coords, colorful):
		OutputFile.__init__(self)

		if colorful:
			self.colors = ['#0000ff','#ff0000','#ff00ff','#ffff00','#00ffff']
		else:
			self.colors = ['#000000']
		self.partial_color = '#888888'

		self.max_coords = max_coords

		#the scale factor is heuristically determined to give a good
		#font size on A4 paper
		self.fontsize = (max_coords).max()/30

		#set up svg DOM
		self.dom = dom.Document()
		self.svg = self.dom.createElement('svg')
		self.svg.setAttribute('xmlns:dc',"http://purl.org/dc/elements/1.1/")
		self.svg.setAttribute('xmlns:cc',"http://creativecommons.org/ns#")
		self.svg.setAttribute('xmlns:rdf',"http://www.w3.org/1999/02/22-rdf-syntax-ns#")
		self.svg.setAttribute('xmlns:svg',"http://www.w3.org/2000/svg")
		self.svg.setAttribute('xmlns',"http://www.w3.org/2000/svg")
		self.svg.setAttribute('width',"%fmm" % (2*max_coords[1],))
		self.svg.setAttribute('height',"%fmm" % (max_coords[0],))
		self.svg.setAttribute('version',"1.1")
		self.svg.setAttribute('id',"svg2")
		self.dom.appendChild(self.svg)

		defs = self.dom.createElement('defs')
		defs.setAttribute('id',"defs4")
		self.svg.appendChild(defs)

	def to_svg(self,point):
		res = np.zeros((2,))
		res[0] = self.max_coords[1] - point[1]
		res[1] = self.max_coords[0] - point[0]
		return mm_to_px*res

	def add_point(self,point,label):
		self.points[label] = point
		self._add_point_svg(self.to_svg(point),label)

	def _add_point_svg(self,point,label):
		text = self.dom.createElement('text')
		text.setAttribute('xml:space',"preserve")
		text.setAttribute('style',"font-size:%dpx;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:#00ff00;fill-opacity:1;stroke:none;font-family:Arial;-inkscape-font-specification:Arial" % self.fontsize)
		text.setAttribute('x',str(point[0]))
		text.setAttribute('y',str(point[1]))
		
		tspan = self.dom.createElement('tspan')
		tspan.setAttribute('x',str(point[0]))
		tspan.setAttribute('y',str(point[1]))
		tspan.appendChild(self.dom.createTextNode(label))
		text.appendChild(tspan)

		self.svg.appendChild(text)

	def add_line(self,p_refs,label,partial=True):
		points = np.array([self.points[str(p_ref)] for p_ref in p_refs])

		color = None
		if partial:
			color = self.partial_color
		else:
			color = choice(self.colors)
		self._add_line_svg(np.array([self.to_svg(p) for p in points]),label, color)

	def _add_line_svg(self,points,label,color):
		g = self.dom.createElement('g')

		#line
		path = self.dom.createElement('path')
		path_string = 'M '
		for i in range(points.shape[0]):
			path_string += "%.6f,%.6f " % tuple(points[i])
		path.setAttribute('d',path_string)
		path.setAttribute('style',"fill:none;stroke:%s;stroke-width:%dpx;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" % (color,self.fontsize/10))
		g.appendChild(path)

		#label
		text = self.dom.createElement('text')
		text.setAttribute('xml:space',"preserve")
		text.setAttribute('style',"font-size:%dpx;font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:%s;fill-opacity:1;stroke:none;font-family:Arial;-inkscape-font-specification:Arial" % (self.fontsize,color))
		text.setAttribute('x',str(points[:,0].sum()*0.5))
		text.setAttribute('y',str(points[:,1].sum()*0.5))
		
		tspan = self.dom.createElement('tspan')
		tspan.setAttribute('x',str(points[:,0].sum()*0.5))
		tspan.setAttribute('y',str(points[:,1].sum()*0.5))
		tspan.appendChild(self.dom.createTextNode('%s' % label))
		text.appendChild(tspan)

		g.appendChild(text)
		self.svg.appendChild(g)

	def to_file(self,fid):
		self.dom.writexml(fid)


if __name__ == "__main__":
	#command line argument parsing
	parser = argparse.ArgumentParser(description="Project a 2D drawing or a 3D model to the floor plane")
	parser.add_argument("--height", type=int, help="Height of the observer in mm", default=1700)
	parser.add_argument("--distance", type=int, help="Distance of the drawing to the observer, only relevant for projection of a 2D drawing", default=100)
	parser.add_argument("--epsilon", type=int, help="Distance of the upper edge of the drawing to the height of the observer. The smaller, the large the projection will get, only applies to 2D drawing", default=50)
	parser.add_argument("--colorful", type=bool, nargs='?', help="Color lines randomly in the output SVG. This is useful while creating a filter file", default=False, const=True)
	parser.add_argument("--filter", type=int, nargs='+', help="Line indices that should be suppressed in the output", default=[])
	parser.add_argument("--partial", type=int, nargs='+', help="Line indices of partially hidden lines, that should be marked in the output", default=[])
	parser.add_argument("--extra", type=int, nargs='+', help="Pairs of point indices for which extra lines should be added", default=[])
	parser.add_argument("--filter-file", type=str, help="File that contains the indices of the edges that should be suppressed in the output", default=None)
	parser.add_argument("inputfile", type=str, help="A .svg or a .stl file")

	params = vars(parser.parse_args())

	#get filtered, partial and extra edges
	filtered = None
	partial = None
	extra = None
	if not params['filter_file'] is None:
		filtered,partial,extra = parse_filter_file(open(params['filter_file']))
	else:
		filtered = params['filter']
		partial = params['partial']
		if not len(params['extra']) % 2 == 0:
			print "Invalid parameters for --extra. Pairs of point indices are required"
		extra = map(list,zip(params['extra'][0::2],params['extra'][1::2]))

	basename = splitext(params["inputfile"])[0]

	input = None
	fid = open(params["inputfile"])
	if splitext(params["inputfile"])[1] == ".stl":
		input =  InputSTL(fid)
	elif splitext(params["inputfile"])[1] == ".svg":
		input =  InputSVG(fid,params["distance"],params["epsilon"],params["height"])
	fid.close()

	#Projection
	r_obj = input.get_points()
	r_proj = np.zeros((r_obj.shape[0],2))
	r_proj[:,0] = params["height"]*r_obj[:,0]/(params["height"] - r_obj[:,2])
	r_proj[:,1] = params["height"]*r_obj[:,1]/(params["height"] - r_obj[:,2])

	points,lines = filter_and_dictify(r_proj,input.get_lines(),filtered, np.array(extra))

	#Maximum proj extensions
	r_max = input.get_max_coords()
	r_pmax = np.zeros((2,))
	for coord in points.values():
		r_pmax = np.maximum(r_pmax,coord)

	origin = np.array([0,0])
	reference = np.array([0.5*r_pmax[0],r_pmax[1]])

	#output svg and table
	svg_out = OutputSVG(r_pmax,params['colorful'])
	tbl_out = OutputTable(origin, reference)
	for i,point in points.iteritems():
		svg_out.add_point(r_proj[i,:],str(i))
		tbl_out.add_point(r_proj[i,:],str(i))
	svg_out.add_point(origin,"Origin")
	svg_out.add_point(reference,"Reference")

	for i,line in lines.iteritems():
		svg_out.add_line(line,str(i),i in partial)
		tbl_out.add_line(line,str(i),i in partial)

	fid = open(basename + "_perspective.svg","w")
	svg_out.to_file(fid)
	fid.close()

	fid = open(basename + "_perspective.table","w")
	tbl_out.to_file(fid)
	fid.close()
