empty - so far

#python packages
python 3.7
pandas
numpy
biopython: conda install -c conda-forge biopython
deeptools

#folding software
ViennaRNA downloaded from:
https://www.tbi.univie.ac.at/RNA/index.html#download


### how to generate documentation (just in case)
1. run Sphinx Quickstart
2. in conf.py file:
- uncomment `import os, sys` and `sys.path.insert(0, os.path.abspath('.'))`
- add 'sphinx.ext.autodoc' in `extensions`
- change `html_theme` = 'sphinx_rtd_theme' (or other preffered)
3. install sphinx_rtd_theme to python (pip or conda install)
4. add code.rts to source folder
5. add each module by typing within code.rts
    `..  automodule:: TTools.methods` 
       `:members:'`
    (this can by divided into sections etc)
6. type make html in terminal to generate documentation
