from setuptools import setup

setup(name='spekpy',
      version='1.0.0',
      description='A Python software toolkit for modelling x-ray spectra from tungsten anode x-ray tubes',
      url='https://bitbucket.org/spekpy/spekpy_release',
      author='Robert Bujila & Gavin Poludniowski',
      author_email='gpoludniowski@gmail.com',
      license='MIT License',
      install_requires=['scipy'],
      packages=['spekpy'],
      zip_safe=False,
      package_data={'spekpy':[r'data/tables/*.dat', r'data/matl_def/*.comp', r'data/matl_usr/*.comp', r'data/state_usr/*.state', r'data/state_def/*.state']},
      include_package_data=True)
