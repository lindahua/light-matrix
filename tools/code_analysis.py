# Python module for analyzing the codes in Light Matrix library

import sys
import os.path
from glob import glob

import re


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


def count_lines(lines):
	"""Get line counting stats"""

	s = CodeStats()
	s.total_lines = len(lines)

	ne = 0
	for line in lines:
		line = line.strip()
		if len(line) > 0:
			ne = ne + 1

	s.total_nonempty_lines = ne
	return s


class ParseError:

	def __init__(self, title, lineno, cause):
		self.title = title
		self.lineno = lineno
		self.cause = cause

	def message(self):
		return "[ParseError]: {0}({1}): {2}".format(self.title, self.lineno, self.cause)



def divide_code_sections(title, lines):
	"""Divide a code file into three sections: head, preamble, and mainbody"""

	hs0 = None
	hs1 = None
	ps0 = None
	ps1 = None
	ms0 = None
	ms1 = None

	status = 'start'
	pragma_guard_open = False
	macro_guard_open = False
	guarded = False

	for i in range(len(lines)):
		lineno = i + 1

		line = lines[i].strip()

		if len(line) == 0:  # empty line
			continue

		if status == 'start':

			if line.startswith('/*'):
				status = 'head'
				hs0 = i
			else:
				raise ParseError(title, lineno, 'Invalid starting line')

		elif status == 'head':
			if line.startswith('*/'):
				hs1 = i + 1
				status = 'preamble'

			elif line.startswith('*'):
				pass 
			else:
				raise ParseError(title, lineno, 'Invalid line in head section')

		elif status == 'preamble':

			if ps0 is None: ps0 = i

			in_preamble = False

			if line == '#ifdef _MSC_VER':
				pragma_guard_open = True
				in_preamble = True

			elif line.startswith('#pragma'):
				in_preamble = True

			elif line.startswith('#ifndef LIGHTMAT_'):
				in_preamble = True
				macro_guard_open = True

			elif line.startswith('#define LIGHTMAT_'):
				if macro_guard_open:
					in_preamble = True
					guarded = True
				macro_guard_open = False

			elif line.startswith('#include'):
				in_preamble = True

			elif line == '#endif':
				in_preamble = pragma_guard_open
				pragma_guard_open = False

			if not in_preamble:
				ps1 = i
				ms0 = i
				ms1 = len(lines)
				break

	if not guarded:
		raise ParseError(title, lineno, "The file is not guarded")

	return (hs0, hs1), (ps0, ps1), (ms0, ms1)


file_decl_pat = re.compile('\*\s+@file (\w+).h')

def parse_file_head(title, lines, sec):

	first = lines[sec[0]].strip()
	last = lines[sec[1]-1].strip()

	assert first.startswith('/*')
	assert last.startswith('*/')

	decl_filename = None

	expect_title = os.path.basename(title)

	for i in range(sec[0], sec[1]):
		lineno = i + 1
		line = lines[i].strip()
		r = file_decl_pat.match(line)
		if r is not None:
			decl_filename = r.group(1)
			if decl_filename != expect_title:
				raise ParseError(title, lineno, 
					"Incorrect file declaration: {0} (expect {1})".format(
						decl_filename, expect_title))

	if decl_filename is None:
		raise ParseError(title, sec[1], "file declaration is not found.")

	return None


def parse_file_preamble(title, lines, sec):

	pragma_guard = None
	macro_guard = None

	gmacro = 'LIGHTMAT_' + os.path.basename(title).upper() + '_H_'

	for i in range(sec[0], sec[1]):
		lineno = i + 1
		line = lines[i].strip()

		if line == '#ifdef _MSC_VER':
			pragma_guard = 'open'

		elif line == '#pragma once':
			if pragma_guard is not 'open':
				raise ParseError(title, lineno, '#pragma once found in invalid place')

			pragma_guard = 'found'

		elif line == '#endif':
			if pragma_guard == 'found':
				pragma_guard = 'closed'

		elif line.startswith('#ifndef LIGHTMAT_'):
			if pragma_guard is not 'closed':
				raise ParseError(title, sec[1], "The pragma guard is not correctly put.")

			cur_macro = line[8:]
			if cur_macro != gmacro:
				raise ParseError(title, lineno, 
					"Incorrect guard macro: {0} (expect {1})".format(cur_macro, gmacro))

			macro_guard = 'open'

		elif line.startswith('#define LIGHTMAT_'):
			if macro_guard is not 'open':
				raise ParseError(title, lineno, "macro guard found in wrong place.")

			cur_macro = line[8:]
			if cur_macro != gmacro:
				raise ParseError(title, lineno, 
					"Incorrect guard macro: {0} (expect {1})".format(cur_macro, gmacro))

	pass


def parse_file_mainbody(title, lines, sec):
	pass



def parse_libcodes(module, filename, lines, verbose=0):
	"""Parse a section of library codes, given in a list of lines"""

	title = os.path.splitext(filename)[0]
	if verbose > 1:
		print "parsing file {0}.h".format(title)

	doc = LibCodeFile(module, filename, title)
	doc.stats = count_lines(lines)

	try:
		hs, ps, ms = divide_code_sections(title, lines)
		doc.head = parse_file_head(title, lines, hs)
		doc.preamble = parse_file_preamble(title, lines, ps)

	except ParseError as err:
		print err.message()
		sys.exit(1)

	if verbose > 1:
		tlines = doc.stats.total_lines
		tlines1 = doc.stats.total_nonempty_lines
		print "stats: total {0} lines | {1} non-empty".format(tlines, tlines1)

	return doc


def parse_libcode_file(module, filename, verbose=0):
	"""Parse a library code file"""

	filepath = module.get_filepath(filename)
	title = os.path.splitext(filename)[0]

	with open(filepath, 'r') as f:
		lines = f.readlines()

	return parse_libcodes(module, filename, lines, verbose)



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
		print "module stats: total {0} lines | {1} non-empty".format(
			md.stats.total_lines, md.stats.total_nonempty_lines)
		print " "

	return md
