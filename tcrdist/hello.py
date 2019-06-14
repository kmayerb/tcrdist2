# hello.py
def say_hello(x=None):
	"""
	say hello

	An example function loosely using the NumPy doc string conventions

	Parameters
	---------
	x : string (defaults to None)

	Returns
	-------
	greeting : string

	Notes
	-----
	https://numpydoc.readthedocs.io/en/latest/format.html#documenting-modules

	Example
	--------
	>>> say_hello("Mike")
	Hello Mike, the VIDD is part of Fred Hutch
	"""
	if x is None:
		greeting = "Hello User, tcrdist is being improved at Fred Hutch"
	else:
		greeting = "Hello {}, tcrdist is being improved at Fred Hutch".format(x)
	return greeting


if __name__ == '__main__':
	say_hello()