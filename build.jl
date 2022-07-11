# load standard packages
using Pkg
# install Downloads packages
Pkg.add("Downloads")
using Downloads

# change into directory of this file
cd(@__DIR__)

# install custom registry
Pkg.Registry.add(Pkg.RegistrySpec(url = "https://github.com/LudwigBoess/LMBRegistry.git"))

# activate environment and install packages
Pkg.activate(".")
Pkg.instantiate()

# download data
@info "downloading test data..."
Downloads.download("http://www.usm.uni-muenchen.de/~lboess/CRESCEND/data.tar.gz", "./data.tar.gz")
@info "done!"

# unpack data 
@info "unpacking data..."
run(`tar –xvzf data.tar.gz`)
@info "done!"