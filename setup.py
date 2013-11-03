from distutils.core import setup

setup(
      name='PyClone-SCN',
      version='0.1.0',
      description='Python tools for analysing clonal evolution using NGS data which accounts for sub-clonal copy number information.',
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',
      url='http://compbio.bccrc.ca',
      package_dir = {'': 'lib'},    
      packages=[ 
                'pyclone',
                'pyclone.post_process',
                'pyclone.post_process.plot'
                ],
      scripts=['PyClone-SCN']
     )
