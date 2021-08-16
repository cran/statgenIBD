.onUnload <- function(libpath) {
  library.dynam.unload("statgenIBD", libpath)
}
