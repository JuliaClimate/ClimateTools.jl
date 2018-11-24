using Pkg

using BinDeps
import CondaBinDeps

if lowercase(get(ENV, "CI", "false")) == "true"

    let basepython = get(ENV, "PYTHON", "python2")
        envpath = joinpath(@__DIR__, "env")
        run(`pip install --user virtualenv`)
        run(`virtualenv --python=$basepython $envpath`)

        if Sys.iswindows()
            python = joinpath(envpath, "Scripts", "python.exe")
        else
            python = joinpath(envpath, "bin", "python2")
            # if Sys.islinux()
            #     run(`sudo apt get install libnetcdf-dev`)
            # elseif Sys.isapple()
            #     run(`sudo apt get install libnetcdf-dev`)
            # end
        end
        # run(`sudo apt-get install libnetcdf-dev`)
        run(`$python -m pip --user install numpy`)
        run(`$python -m pip --user install scipy`)
        run(`$python -m pip --user install matplotlib`)
        run(`$python -m pip --user install https://github.com/matplotlib/basemap/archive/v1.0.7rel.tar.gz`)
        run(`$python -m pip --user install git+https://github.com/matplotlib/cmocean`)
        run(`$python -m pip --user install netcdf4`)

        ENV["PYTHON"] = python
        Pkg.build("PyCall")

        if VERSION >= v"0.7.0"
            using Libdl
        end

        function validate_netcdf_version(name,handle)
            f = Libdl.dlsym_e(handle, "nc_inq_libvers")
            #
            # Example
            # banner = "4.5.6.7 of Apr 1 2000 00:00:00"
            banner = unsafe_string(ccall(f,Ptr{UInt8},()))

            # remove the date
            # verstr = "4.5.6.7"
            verstr = split(banner)[1]

            # vernumbers = ["4","5","6","7"]
            vernumbers = split(verstr,r"[.-]")

            # major_minor_patch = ["4","5","6"]
            major_minor_patch = vernumbers[1:min(3,length(vernumbers))]

            # ver = v"4.5.6"
            ver = VersionNumber([parse(Int,s) for s in major_minor_patch]...)

            return ver > v"4.2"
        end

        @BinDeps.setup
        libnetcdf = library_dependency("libnetcdf", aliases = ["libnetcdf4","libnetcdf-7","netcdf"], validate = validate_netcdf_version)

        CondaBinDeps.Conda.add_channel("anaconda")
        provides(CondaBinDeps.Manager, "libnetcdf", libnetcdf)
        provides(AptGet, "libnetcdf-dev", libnetcdf, os = :Linux)
        # provides(Yum, "netcdf-devel", libnetcdf, os = :Linux)

        @BinDeps.install Dict(:libnetcdf => :libnetcdf)

    end
end
