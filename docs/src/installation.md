# Installation

ClimateTools.jl requires some python packages for mapping purpose. it is currently using Basemap under the hood. To ensure that ClimateTools works properly, it is recommended to use a Python distribution that can properly load the required python packages and build `PyCall` with the same python distribution.

## Approach no. 1 Use main system python distribution

The first approach is to use the system python distribution and ensure that you can import the following packages.

**1.1 Python dependencies**

* matplotlib (tested with version 2.0.1)
* basemap (tested with version 1.0.7)
* scipy (tested with version 1.0.1)
* cmocean

*Note. Installing Basemap with `pip` for python 3.6+ seems problematic. See approach no 2a below for an installation with Anaconda*

**1.2 Building PyCall**
After the confirmation that the Python dependencies can be loaded in Python, the user needs to build PyCall with the same Python version. Alternatively, if PyCall is already built, it may be only a matter of installing the Python dependencies with the PyCall's Python version by using `pip`.

```julia
ENV["PYTHON"]="path_to_python_distribution"
pkg> build PyCall
```

## Approach no. 2a. Use a Conda virtual environment and link PyCall.jl to it

If the main system python distribution is not set-up properly, the recommended approach is to use a Conda environment. Basemap is much easier to install through this approach. The is simply a matter of creating a virtual python environment with conda and installing the required packages and then linking Julia's PyCall.jl package to this virtual conda environment. More information can be found in [PyCall](https://github.com/JuliaPy/PyCall.jl) documentation for custom Conda environment.

**2.1a Create a virtual environment with Conda**

```bash
$ conda create -name ClimateTools pytyhon=3.6
$ conda activate ClimateTools
$ pip install numpy scipy matplotlib
$ conda install -c anaconda basemap
$ pip install cmocean
$ which python # gives you the path of Conda virtual environment to use in the next steps.
```

Once those packages are installed, you need to tell PyCall.jl (in Julia) to use this conda environement in Julia.

```julia
julia> using PyCall
julia> ENV["PYTHON"]="...PATH TO CONDA PYTHON..." # find the path with "which python" at previous step
julia> using Pkg; Pkg.build("PyCall")
```

## Approach no. 2b. Build a python virtual environment with `virtualenv` and link PyCall.jl to it

This is another way to install the required python package. This approach is less recommended than `Approach no. 2a`. This approach consist of using `virtualenv` command and install the required packages inside this virtual env and then link PyCall.jl to this virtual python environment. More information can be found in [PyCall](https://github.com/JuliaPy/PyCall.jl) documentation.

**2.1b Create a virtual environment with Python 2.7.x.**

```bash
$ virtualenv --python=/usr/bin/python2 /path/to/venv
$ /path/to/venv/bin/python -m pip install numpy
$ /path/to/venv/bin/python -m pip install scipy
$ /path/to/venv/bin/python -m pip install matplotlib
$ /path/to/venv/bin/python -m pip install https://github.com/matplotlib/basemap/archive/v1.0.7rel.tar.gz
$ /path/to/venv/bin/python -m pip install git+https://github.com/matplotlib/cmocean
```

**2.2 Testing Python installation**

no matter the approach used, it is recommended to test the python environment and see if you can import the required python packages in the newly created virtual env python.

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

once the python dependencies are properly installed you can then install ClimateTools in Julia.

```julia
pkg> add ClimateTools # Tagged release
```
