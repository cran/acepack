.First.lib <- function(lib, pkg) {
  library.dynam("acepack", pkg, lib)
}

if (version$minor < "0.62")
  library.dynam("acepack")

