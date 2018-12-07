from setuptools import setup, find_packages

setup(name='askap-image-diagnostic',
      version='0.1a',
      description='Perform analysis tasks on ASKAP images.',
      url='http://github.com/ajstewart/askap-image-diagnostic',
      author='Adam Stewart',
      author_email='adam.stewart@sydney.edu.au',
      license='',
      packages=find_packages(),
      # install_requires=['astropy',
                        # 'numpy'],
      scripts=["bin/processASKAPimage.py"],
      include_package_data=True)
      