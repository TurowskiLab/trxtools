import subprocess
import os

################################################
#############        bash

def listPaths(folder=".", suffix=str(), nameElem=None):
    '''List all file paths in a folder with a specific suffix and optional name element.

    :param folder: Folder to list files from, defaults to "."
    :type folder: str, optional
    :param suffix: Suffix of the files to list, defaults to an empty string
    :type suffix: str, optional
    :param nameElem: Optional element that must be in the file name, defaults to None
    :type nameElem: str, optional

    :return: List of file paths that match the criteria
    :rtype: list

    :example:

    >>> listPaths(folder=".", suffix=".txt", nameElem="file")
    ['file1.txt', 'file2.txt']
    '''
    if nameElem:
        return [folder + file for file in os.listdir(folder) if file.endswith(suffix) and nameElem in file]
    else:
        return [folder + file for file in os.listdir(folder) if file.endswith(suffix)]

def bashCommand(bashCommand=str()):
    '''Run command in bash using subprocess.call()
    
    :param bashCommand: str with command

    :return: None
    '''
    pipes = subprocess.Popen(bashCommand,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
            err_msg = "%s. Code: %s" % (std_err.strip(), pipes.returncode)
            raise Exception(err_msg)
