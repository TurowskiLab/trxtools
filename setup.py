#!/usr/bin/env python

import setuptools

# TODO requirements

setuptools.setup(
    name='trxtools',
    version='0.2.1',
    author="Tomasz W. Turowski, Jan Mikołajczyk",
    author_email="twturowski@gmail.com",
    description='Python tools facilitating bioinformatic analysis of nascent transcripts and transcriptomic data.',
    long_description=open('README.md').read(),
    license='LICENSE.txt',
    keywords="bioinformatics RNA RNAfolding transcription BigWig",
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License', 
        'Operating System :: OS Independent',        
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Software Development :: Libraries :: Python Modules'],
    scripts=['scripts/SAM2profiles.py',
        'scripts/SAM2profilesGenomic.py',
        'scripts/mergeSalmon.py',
        'scripts/csv2pickle.py',
        'scripts/genome2WindowsGTF.py'],
    install_requires=[
        'scikit-learn','pandas', 'numpy', 
        'pyBigWig', 'adjustText', 'scipy',
        'seaborn', 'matplotlib'],
    python_requires='>=3.6'
)
