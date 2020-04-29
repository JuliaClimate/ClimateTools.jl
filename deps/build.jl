using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"    

        ENV["PYTHON"] = ""
        Pkg.build("PyCall")
      
end
