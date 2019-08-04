from Bio.Application import _Option, _Switch, AbstractCommandline
from Bio.Emboss.Applications import _EmbossCommandLine, _EmbossMinimalCommandLine

class EmmaCommandline(_EmbossMinimalCommandLine):
	"""
	Commandline object for the emma program from EMBOSS. 
	"""

	def __init__(self, cmd="emma", **kwargs):

		self.parameters = [ 
		 _Option(["-sequence", "sequence"], 
				 "seqall", 
				 filename=True, 
				 is_required=True),
		 _Option(["-outseq", "outseq"], 
				 "seqoutset", 
				 filename=True, 
				 is_required=True), 
		 _Option(["-dendoutfile", "dendoutfile"], 
				 "outfile", 
				 filename=True, 
				 is_required=True)]
		 # _Option(["-gapextend", "gapextend"], 
			# 	 "Gap extension penalty", 
			# 	 is_required=True), 
		 # _Option(["-datafile", "datafile"], 
			# 	 "Matrix file", 
			# 	 filename=True), 

		 # _Switch(["-nobrief", "nobrief"], 
			# 	 "Display extended identity and similarity"), 
		 # _Switch(["-brief", "brief"], 
			# 	 "Display brief identity and similarity"),
		 # _Option(["-similarity", "similarity"], 
			# 	 "Display percent identity and similarity"), 
		 # _Option(["-snucleotide", "snucleotide"], 
			# 	 "Sequences are nucleotide (boolean)"), 
		 # _Option(["-sprotein", "sprotein"], 
			# 	 "Sequences are protein (boolean)"), 
		 # _Option(["-aformat", "aformat"], 
			# 	 "Display output in a different specified output format")] 

		_EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

class SuperMatcherCommandline(_EmbossCommandLine):
	"""
	Commandline object for the supermatcher program from EMBOSS. 
	"""

	def __init__(self, cmd="supermatcher", **kwargs):
		
		self.parameters = [ 
		 _Option(["-asequence", "asequence"], 
				 "First sequence to align", 
				 filename=True, 
				 is_required=True),
		 _Option(["-bsequence", "bsequence"], 
				 "Second sequence to align", 
				 filename=True, 
				 is_required=True), 
		 _Option(["-gapopen", "gapopen"], 
				 "Gap open penalty", 
				 is_required=True), 
		 _Option(["-gapextend", "gapextend"], 
				 "Gap extension penalty", 
				 is_required=True), 
		 _Option(["-datafile", "datafile"], 
				 "Matrix file", 
				 filename=True),

		 _Switch(["-nobrief", "nobrief"], 
				 "Display extended identity and similarity"), 
		 _Switch(["-brief", "brief"], 
				 "Display brief identity and similarity"),
		 _Option(["-similarity", "similarity"], 
				 "Display percent identity and similarity"), 
		 _Option(["-snucleotide", "snucleotide"], 
				 "Sequences are nucleotide (boolean)"), 
		 _Option(["-sprotein", "sprotein"], 
				 "Sequences are protein (boolean)"), 
		 _Option(["-aformat", "aformat"], 
				 "Display output in a different specified output format")] 
		_EmbossCommandLine.__init__(self, cmd, **kwargs)
