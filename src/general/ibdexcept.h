#ifndef IBD_EXCEPTION_HEADER
#define IBD_EXCEPTION_HEADER

#include <string>
#include <stdexcept>
#include <iostream>

namespace ibd
{

class ibd_error : public std::runtime_error
{
public:
  ibd_error(const std::string& what_arg) : std::runtime_error(what_arg) {}
  virtual ~ibd_error() throw() {;}
};

class ibd_file_error : public ibd_error
{
public:
  ibd_file_error(const std::string& filename,
                 int line_nr,
                 const std::string& what_arg);
  virtual ~ibd_file_error() throw() {;}
};

typedef ibd_error MQMlib_error;
typedef ibd_file_error MQMlib_file_error;

}

#endif

