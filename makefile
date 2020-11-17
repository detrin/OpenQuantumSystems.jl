
REPOSITORY=https://github.com/detrin/OpenQuantumSystems.jl

tests:
	julia --project -e 'using Pkg; Pkg.build(); Pkg.test()'

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
