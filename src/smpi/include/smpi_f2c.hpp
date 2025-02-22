/* Handle Fortan - C conversion for MPI Types*/

/* Copyright (c) 2010-2019. The SimGrid Team.
 * All rights reserved.                                                     */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SMPI_F2C_HPP_INCLUDED
#define SMPI_F2C_HPP_INCLUDED

#include <unordered_map>
#include <string>

#define KEY_SIZE (sizeof(int) * 2 + 1)

namespace simgrid{
namespace smpi{

class F2C {
  private:
    // We use a single lookup table for every type.
    // Beware of collisions if id in mpif.h is not unique
    static std::unordered_map<std::string, F2C*>* f2c_lookup_;
    static int f2c_id_;
    int my_f2c_id_;

  protected:
    static std::unordered_map<std::string, F2C*>* f2c_lookup();
    static void set_f2c_lookup(std::unordered_map<std::string, F2C*>* map);
    static int f2c_id();
    static void f2c_id_increment();

  public:
    char* get_my_key(char* key);
    static char* get_key(char* key, int id);
    static void delete_lookup();
    static std::unordered_map<std::string, F2C*>* lookup();
    F2C() : my_f2c_id_(-1){}
    virtual ~F2C() = default;

    //Override these to handle specific values.
    virtual int add_f();
    static void free_f(int id);
    virtual int c2f();
    static void print_f2c_lookup();
    // This method should be overridden in all subclasses to avoid casting the result when calling it.
    // For the default one, the MPI_*_NULL returned is assumed to be NULL.
    static F2C* f2c(int id);
};

}
}

#endif
