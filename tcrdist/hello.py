# hello.py
def say_hello():
	"""
	say hello

	An example function loosely using the NumPy doc string conventions

	Parameters
	---------

	Returns
	-------
	greeting : string

	Notes
	-----
	https://numpydoc.readthedocs.io/en/latest/format.html#documenting-modules

	Example
	--------
	"""

	greeting = "Hello: 'By recombination, random insertion, deletion and substitution, " \
			   "the small set of genes that encode the T-cell receptor has the potential to create " \
			   "between 10^15 and 10^20 TCR clonotypes ... However, the actual diversity of a persons " \
			   "TCR repertoire cannot possibly lie in this range. There are only an estimated 10^13 cells " \
			   "in the human body [3]' -- Laydon et al. 2015. PMC4528489"
	print(greeting)


if __name__ == '__main__':
	say_hello()