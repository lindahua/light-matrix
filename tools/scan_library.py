# A python script to scan the lbrary

from code_analysis import *


def run_library_scan():

	# initialize configuration

	print "Initializing configuration ..."

	cfg = LibConfig()
	print "rootpath =", cfg.rootpath
	print "incpath  =", cfg.incpath
	print "the library has {0} modules:".format(len(cfg.module_names))
	for md in cfg.module_names:
		print "    {0}".format(md)

	stats = CodeStats()

	verbose = 3
	for mname in cfg.module_names:
		md = scan_module(cfg, mname, verbose)
		stats.add_stats(md.stats)

	print "overall stats: total {0} lines".format(stats.total_lines)
	print "Done!\n"


if __name__ == '__main__':

	run_library_scan()

