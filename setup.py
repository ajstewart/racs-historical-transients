from setuptools import setup, find_packages

setup(name='askap-image-diagnostic',
      version='0.2',
      description='Perform analysis tasks on ASKAP images.',
      url='http://github.com/ajstewart/askap-image-diagnostic',
      author='Adam Stewart',
      author_email='adam.stewart@sydney.edu.au',
      license='',
      packages=find_packages(),
      install_requires=['numpy',
                        'astropy',
                        'matplotlib',
                        'aplpy',
                        'AegeanTools',
                        'pandas',
                        'astroquery',
                        'colorlog',
                        # 'jupyter'
                        ],
      scripts=["bin/processASKAPimage.py"],
      include_package_data=True)
      