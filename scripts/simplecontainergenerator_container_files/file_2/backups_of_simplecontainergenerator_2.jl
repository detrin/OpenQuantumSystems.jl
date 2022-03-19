empty!(Base.DEPOT_PATH)
pushfirst!(Base.DEPOT_PATH, "/opt/simplecontainergenerator_containers/depot_backup_packagecompiler")
import Pkg
import Pkg; Pkg.add(Pkg.PackageSpec(name = "PackageCompiler", version = "1.2.2 - 1"));

