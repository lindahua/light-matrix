# Python module for analyzing the codes in Light Matrix library

import os.path
from glob import glob


class LibConfig:
	"""The class to represent the configuration of the library"""

	def __init__(self):

		tools_path = os.path.dirname(__file__)
		rpath = os.path.dirname(tools_path)

		self.rootpath = rpath
		self.incpath = os.path.join(rpath, 'light_mat')

		if not os.path.isdir(self.incpath):
			raise RuntimeError("The incpath {0} does not exist".format(self.incpath))

		module_list = os.path.join(tools_path, "modules.lst")
		with open(module_list, 'r') as f:
			mnames = f.readlines()

		mnames = [x.strip() for x in mnames if len(x.strip()) > 0]

		for md in mnames:
			mpath = os.path.join(self.incpath, md)
			if not os.path.isdir(mpath):
				raise RuntimeError("The module path {0} does not exist".format(mpath))

		self.module_names = mnames


	def get_module_path(self, mname):
		"""Get the full path of a module"""

		if mname not in self.module_names:
			raise RuntimeError("The module {0} does not exist".format(mname))

		return os.path.join(self.incpath, mname)



class LibModule:
	"""The class to represent a library module"""

	def __init__(self, path, name):
		self.path = path
		self.name = name
		self.files = []
		self.stats = CodeStats()


	def get_filepath(self, fname):
		return os.path.join(self.path, fname)


class CodeStats:
	"""The class to capture statistics of codes"""

	def __init__(self):
		self.num_files = 0
		self.total_lines = 0
		self.total_nonempty_lines = 0

	def add_stats(self, s):

		self.num_files = self.num_files + s.num_files
		self.total_lines = self.total_lines + s.total_lines
		self.total_nonempty_lines = self.total_nonempty_lines + s.total_nonempty_lines


class LibCodeFile:
	"""The class to represent a single library file (usually .h)""" 

	def __init__(self, module, filename, title):

		self.module = module
		self.filename = filename
		self.title = title
		self.head = None
		self.preamble = None
		self.mainbody = None
		self.stats = None



# The codes to parse a C++ library code file


def parse_file_head(title, lines, sec):
	pass

def parse_file_preamble(title, lines, sec):
	pass


def parse_file_mainbody(title, lines, sec):
	pass


def parse_libcodes(title, lines):
	"""Parse a section of library codes, given in a list of lines"""

	[hs, ps, ms] = divide_code_sections(lines)

	head = parse_file_head(title, lines, hs)
	preamble = parse_file_preamble(title, lines, ps)
	mainbody = parse_file_mainbody(title, lines, ms)

	return (head, preamble, mainbody)


def parse_libcode_file(module, filename, verbose=0):
	"""Parse a library code file"""

	filepath = module.get_filepath(filename)
	title = os.path.splitext(filename)[0]

	if verbose > 1:
		print "parsing file {0}.h ...".format(title)

	with open(filepath, 'r') as f:
		lines = f.readlines()

	if verbose > 1:
		print "stats: total {0} lines".format(len(lines))

	stats = CodeStats()
	stats.total_lines = len(lines)

	doc = LibCodeFile(module, filename, title)
	doc.stats = stats

	return doc




# The codes to scan a module

def scan_module(libcfg, module_name, verbose=0):

	if verbose > 0:
		print "scanning module {0} ...".format(module_name)

	mpath = libcfg.get_module_path(module_name)

	md = LibModule(mpath, module_name)

	files = glob(os.path.join(mpath, '*.h'))
	for fpath in files:
		fname = os.path.basename(fpath)
		doc = parse_libcode_file(md, fname, verbose)
		md.stats.add_stats(doc.stats)

	internal_files = glob(os.path.join(mpath, 'internal', '*.h'))
	for fpath in internal_files:
		fname = os.path.join('internal', os.path.basename(fpath))
		doc = parse_libcode_file(md, fname, verbose)
		md.stats.add_stats(doc.stats)

	if verbose > 0:
		print "module stats: total {0} lines".format(md.stats.total_lines)
		print " "

	return md
