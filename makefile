
REPOSITORY=https://github.com/detrin/OpenQuantumSystems.jl

install_from_repo:
	julia -e 'using Pkg; Pkg.add(path="../OpenQuantumSystems.jl")'

remove:
	julia -e 'using Pkg; Pkg.rm("OpenQuantumSystems")'

tests:
	julia --project -e 'using Pkg; Pkg.build(); Pkg.test()'

tests_dev:
	julia test/runtests_dev.jl

update_devel:
	git checkout master
	git pull
	git checkout devel
	git merge master

git_update_master:
	git fetch upstream
	git checkout master
	git merge upstream/master

git_add_upstream:
	git remote add upstream ${REPOSITORY}

benchmark:
	julia test/benchmark.jl | tee Benchmarks.md

