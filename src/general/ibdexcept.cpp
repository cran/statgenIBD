#include "convert.h"
#include "ibdexcept.h"

ibd::ibd_file_error::ibd_file_error(const std::string& filename,
                                    int line_nr,
                                    const std::string& what_arg)
  : ibd_error("file: " + filename + ",line " + stringify(line_nr) + ": " + what_arg) {}

