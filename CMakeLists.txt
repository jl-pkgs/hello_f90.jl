cmake_minimum_required (VERSION 3.16.3)

project(hello VERSION 1.0.1 DESCRIPTION "mylib description" LANGUAGES C Fortran)


add_library(hello SHARED
  mod_diag_functions.f90
)

add_library(hello_static STATIC
  mod_diag_functions.f90
)
