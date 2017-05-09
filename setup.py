from distutils.core import setup

setup(name='trm.observing',
      version='0.5',
      packages = ['trm', 'trm.observing'],
      scripts=['scripts/eplanner.py', ],

      # metadata
      author='Tom Marsh',
      author_email='t.r.marsh@warwick.ac.uk',
      description="Python observing utility module",
      url='http://www.astro.warwick.ac.uk/',
      long_description="""
observing mainly provides a script to help schedule eclipses along
with a few classes to help other build other scripts that need to use
the position / ephemeris files.
""",

      )

