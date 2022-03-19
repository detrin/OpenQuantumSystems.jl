empty!(Base.DEPOT_PATH)
pushfirst!(Base.DEPOT_PATH, "/opt/simplecontainergenerator_containers/depot_backup_simplecontainergenerator")
import Pkg
import Pkg; Pkg.add(Pkg.PackageSpec(name = "SimpleContainerGenerator"));

