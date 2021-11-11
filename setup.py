from setuptools import setup, find_packages

setup(
    name='AutoPhot',
    version='0.1.0-alpha',
    url='https://github.com/almattil/autophot.git',
    author='Author Name',
    author_email='author@gmail.com',
    description='Description of my package',
    scripts=['bin/AUTOPHOT'],
    data_files = [('',["conf/*"])],
    packages=find_packages(),    
    install_requires=['numpy >= 1.20.2', 'matplotlib >= 3.4.1'],
)
