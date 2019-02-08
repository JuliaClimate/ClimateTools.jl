# Installation

## Approach no. 1 Use main system python distribution

ClimateTools need some Python dependencies for mapping purpose. To ensure that ClimateTools works properly, it is recommended to use a Python distribution that can properly load the following python modules and build `PyCall` with the same python distribution.

**1.1 Python dependencies**

* matplotlib (tested with version 2.0.1)
* basemap (tested with version 1.0.7)
* scipy (tested with version 1.0.1)
* cmocean

*Note. Installing Basemap for python 3.6+ seems problematic.*

**1.2 Building PyCall**
After the confirmation that the Python dependencies can be loaded in Python, the user needs to build PyCall with the same Python version. Alternatively, if PyCall is already built, it may be only a matter of installing the Python dependencies with the PyCall's Python version by using `pip`.

```julia
ENV["PYTHON"]="path_to_python_distribution"
pkg> build PyCall
```

## Approach no. 2. Build a python virtual environment and link PyCall.jl to it

One approach to ensure that the right python dependencies are installed is to use a virtual environment. More information can be found in [PyCall](https://github.com/JuliaPy/PyCall.jl) documentation.

**2.1 Create a virtual environment with Python 2.7.x.**

```bash
$ virtualenv --python=/usr/bin/python2 /path/to/venv
$ /path/to/venv/bin/python -m pip install numpy
$ /path/to/venv/bin/python -m pip install scipy
$ /path/to/venv/bin/python -m pip install matplotlib
$ /path/to/venv/bin/python -m pip install https://github.com/matplotlib/basemap/archive/v1.0.7rel.tar.gz
$ /path/to/venv/bin/python -m pip install git+https://github.com/matplotlib/cmocean
```

**2.2 Testing Python installation**

Launch newly created virtual env python

```python
#bash
$ /path/to/venv/bin/python # launch virtual env python
```
Ensure that you can load the appropriate python packages inside the python interpreter.

```python
#python
>>> import mpl_toolkits.basemap as basemap
>>> import matplotlib.pyplot as plt
>>> import cmocean as cm
>>> import scipy as sc
```

**2.3 Build PyCall with the new venv python**

Once the virtual environment python is tested, it's a matter of telling PyCall to use this distribution.

```julia
# julia
julia> ENV["PYTHON"] = "/path/to/venv/bin/python"
julia> using Pkg;Pkg.build("PyCall")
julia> exit()
# re-enter julia
julia> using ClimateTools
julia> using Pkg; Pkg.test("ClimateTools")
```

## Installing ClimateTools.jl

```julia
pkg> add ClimateTools # Tagged release
```
