empty!(Base.DEPOT_PATH)
pushfirst!(Base.DEPOT_PATH, "/opt/simplecontainergenerator_containers/packagecompiler_depot")
import Pkg
import Pkg; Pkg.add(Pkg.PackageSpec(name = "PackageCompiler", version = "1.2.2 - 1"));

