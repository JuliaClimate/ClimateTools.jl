using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"
    let basepython = get(ENV, "PYTHON", "python")
        envpath = joinpath(@__DIR__, "env")
        run(`virtualenv --python=$basepython $envpath`)

        if Sys.iswindows()
            python = joinpath(envpath, "Scripts", "python.exe")
        else
            python = joinpath(envpath, "bin", "python")
        end
        run(`$python -m pip install basemap`)
        run(`$python -m pip install scipy`)
        run(`$python -m pip install matplotlib`)
        run(`$python -m pip install cmocean`)

        ENV["PYTHON"] = python
        Pkg.build("PyCall")
    end
end
