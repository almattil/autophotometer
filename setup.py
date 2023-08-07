from setuptools import setup, find_packages
import os

setup(
    name='autophotometer',
    version='1.0.0',
    url='https://github.com/almattil/autophotometer.git',
    author='Aleksi Mattila, Vitaly Neustroev',
    author_email='aleksi.a.mattila@gmail.com',
    description='Description of my package',
    python_requires='>=3.7',
    #install_requires=['numpy >= 1.20.2', 'matplotlib >= 3.4.1'],
    install_requires=[
        'numpy >= 1.20.2',      #I have 1.21.4
        'astropy >= 4.2.1',
        'ccdproc >= 2.2.0',
        'photutils >= 1.2.0',   #1.3.0
        'plotille >= 3.8.0',    #4.0.2
        'requests >= 2.25.1',   #2.28.1
        'scipy >= 1.4.1',
        'matplotlib >= 3.4.1'], #3.5.0   
    packages=['autophotometer'],
#    package_dir={'autophoto': 'src/mypkg'},
#    package_data={'mypkg': ['data/*.dat']},

    package_data={'autophotometer': ['conf/*']},
    data_files=[(os.path.expanduser('~/.autophotometer'), ['conf/conf_autophot.ini','conf/default.nnw','conf/default.psf',
                                   'conf/default.sex','conf/ph2conf.sex','conf/run1.param','conf/run2.param','conf/scamp.conf'])],
    include_package_data=True,
    entry_points = {'console_scripts': ['AUTOPHOTOMETER=autophotometer.autophot:cli_main']},
)
