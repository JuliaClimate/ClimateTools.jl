# Installation

## Installing ClimateTools.jl

Once the Python dependencies are properly installed you can then install ClimateTools in Julia.

```julia
pkg> add ClimateTools # Tagged release
```

## Installing ClimatePlots.jl


Plotting data and results is done through the [ClimatePlots.jl](https://github.com/JuliaClimate/ClimatePlots.jl) package.

Under the hood, ClimatePlots.jl requires some Python packages for mapping purpose and is currently using Cartopy. To ensure that ClimatePlots works properly, it is recommended to use Julia's built-in Python distribution which can properly load the required Python packages.

Hence, the best approach is to simply configure PyCall as follows:

```julia
julia> using Pkg; Pkg.add("PyCall")
julia> ENV["PYTHON"]=""
julia> Pkg.build("PyCall")
```
Then, once PyCall configured, installation is the standard command.

```julia
pkg> add ClimatePlots # Tagged release
```

## Alternatives Python distribution

It is also possible to use your own Python distribution instead of Julia's Python.

### Approach no. 1 Use Anaconda or main system Python distribution

The first approach is to use Anaconda or the system Python distribution and ensure that you can import the following packages. By far, the easiest solution is to use Anaconda, which can easily install the Basemap dependency.

**1.1 Python dependencies**

* matplotlib
* cartopy
* cmocean

**1.2 Testing Python installation**

No matter the approach used, it is recommended to test the Python environment and see if you can import the required Python packages.

```python
#python
>>> import cartopy
>>> import matplotlib.pyplot as plt
>>> import cmocean as cm
```

**1.2 Building PyCall**
After the confirmation that the Python dependencies can be loaded in Python, the user needs to build PyCall with the same Python version. Alternatively, if PyCall is already built, it may be only a matter of installing the Python dependencies with the PyCall's Python version by using `pip`.

```julia
ENV["PYTHON"]="path_to_python_distribution"
pkg> build PyCall
```

## Approach no. 2. Use a Conda virtual environment and link PyCall.jl to it

If the main system Python distribution is not set-up properly, the recommended approach is to use a Conda environment. Cartopy is much easier to install through this approach. It is simply a matter of creating a virtual Python environment with Conda and installing the required packages and then linking Julia's PyCall.jl package to this virtual Conda environment. More information can be found in [PyCall](https://github.com/JuliaPy/PyCall.jl) documentation for custom Conda environment.

**Create a virtual environment with Conda**

```bash
$ conda create -name ClimateTools python=3.6
$ conda activate ClimateTools
$ conda install -c conda-forge matplotlib --strict-channel-priority
$ conda install -c conda-forge cartopy --strict-channel-priority
$ conda install -c conda-forge cmocean --strict-channel-priority
$ which python # gives you the path of Conda virtual environment to use in the next steps.
```

Once those packages are installed, you need to tell PyCall.jl (in Julia) to use this Conda environement in Julia.

```julia
julia> using PyCall
julia> ENV["PYTHON"]="...PATH TO CONDA PYTHON..." # find the path with "which python" at previous step
julia> using Pkg; Pkg.build("PyCall")
```
