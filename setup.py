from setuptools import setup, find_packages
import os

setup(
    name='autophot',
    version='0.2.0a1',
    url='https://github.com/almattil/autophot.git',
    author='Aleksi Mattila, Vitaly Neustroev',
    author_email='author@gmail.com',
    description='Description of my package',
    python_requires='>=3.7',
    install_requires=['numpy >= 1.20.2', 'matplotlib >= 3.4.1'],
    packages=['autophot'],
#    package_dir={'autophoto': 'src/mypkg'},
#    package_data={'mypkg': ['data/*.dat']},

    package_data={'autophot': ['conf/*']},
    data_files=[(os.path.expanduser('~/.autophot'), ['conf/conf_autophot.ini','conf/default.nnw','conf/default.psf',
                                   'conf/default.sex','conf/ph2conf.sex','conf/run1.param','conf/run2.param','conf/scamp.conf'])],
    include_package_data=True,
    entry_points = {'console_scripts': ['AUTOPHOT=autophot.autophot:cli_main']},
)
