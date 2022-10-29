import sys
from importlib.metadata import version 
from .create_windows import main as create_windows
from .rename_fasta import main as rename_fasta

def main()->None:
	"""Entry point to isocomp"""

	isocomp_version = version('isocomp')

	version_stmnt = f'Isocomp version: {isocomp_version}'  
	# Stuff to print before offering a list of tools
	preamble = [
		version_stmnt,
		'Description: A set of tools to discover novel isoforms',
		'Available tools are:',
	]

    # list of available tools
	# preface each with \t
	tool_list = [
		'\tcreate_windows',
		'\trename_fasta'
	]

	# Anything to print after the list of tools
	epigraph = [
		'For usage instructions, enter isocomp <tool_name> --help'
	]

	help_str = "\n".join(
		["\n".join(preamble),"\n".join(tool_list),"\n".join(epigraph)])
	
	try:
		tool = sys.argv[1]
	except IndexError:
		tool = 'default'

	if tool == '--version': 
		print(version_stmnt)
	elif tool == '--help': 
		print(help_str)
	elif tool == 'create_windows': 
		create_windows(sys.argv[2:])
	elif tool == 'rename_fasta':
		rename_fasta(sys.argv[2:])
	else:
		print(help_str)

if __name__ == "__main__":
    sys.exit(main())
