# pypopgen3

This repository provides tools for population genomic analyses. It mainly contains python3 libraries and scripts, but also some shell and awk scripts. 
## Installation
You can either just download a .zip and decompress, or clone the respository
```
git clone https://github.com/feilchenfeldt/pypopgen3.git
```
In order for the python modules to be found by your system python, you need to add the local path to the repository to your python path. To achieve this permanently, add  a `.pth` file containing this path to your python site-packages folder. For example, you to once:
```
#locate your python site package folder
python -m site --user-site
#outputs <your_site_package_folder>
cd <your_site_package_folder>
echo <path_to_parent_folder_of_pypopgen3> >> mypackages.pth 
```
Alternatively, you can just add the folder to the python path in each python session in which you want to use pypopgen3. 
```
import sys
sys.path.append('<path_to_parent_folder_of_pypopgen3>')
```

## Python modules
### Dependencies
The python modules depend on some python packages, including `numpy`, `pandas`, `ete3`, `matplotlib`. You will need to install these packages manually, for example with `conda` or `pip`. 
### Usage
Start python(3) and load modules from pypopgen3 like
```
from pypopgen3.modules import treetools
```
