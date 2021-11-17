@ECHO OFF

IF "%1"=="install_dev" (
    julia -e "using Pkg; Pkg.add(path=\"../OpenQuantumSystems.jl\")"
)
IF "%1"=="install_dev2" (
    julia -e "using Pkg; Pkg.develop(path=\"../OpenQuantumSystems.jl\")"
)
IF "%1"=="remove" (
    julia -e "using Pkg; Pkg.rm(\"OpenQuantumSystems\")"
)
IF "%1"=="tests" (
    julia --project -e "using Pkg; Pkg.build(); Pkg.test()"
)
IF "%1"=="test" (
    julia --project -e "using Pkg; Pkg.build(); Pkg.test(test_args = [\""%2"\"])"
)
IF "%1"=="test_dev" (
    julia test/runtests_dev.jl
)
IF "%1"=="update_devel" (
	git checkout master
	git pull
	git checkout devel
	git merge master
)
IF "%1"=="git_update_master" (
	git fetch upstream
	git checkout master
	git merge upstream/master
)
IF "%1"=="git_add_upstream" (
	git remote add upstream https://github.com/detrin/OpenQuantumSystems.jl
)
IF "%1"=="format" (
	julia -e "using JuliaFormatter; format(\".\")"
)
IF "%1"=="docs_generate" (
	julia --project=docs docs/make.jl local
)
IF "%1"=="benchmark" (
	julia -e "using Pkg; Pkg.add(path=\"../OpenQuantumSystems.jl\")"
	julia .\benchmark\benchmark.jl
	python .\benchmark\benchmark_plot.py
)
:: Set-Alias -Name julia1.7 -Value C:\Users\daniel.herman\AppData\Local\Programs\Julia-1.7.0-rc3\bin\julia.exe