from distutils.core import setup

setup(
      name='PyClone',
      version='0.3.0',
      description='Python tools for analysing clonal evolution using NGS data.',
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',
      url='http://compbio.bccrc.ca',
      package_dir = {'': 'lib'},    
      packages=[ 
                'pyclone'
                ],
      scripts=['PyClone']
     )
