/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <iostream>
#include <Utils.h>
#include <time.h>
#include <stdio.h>
#include <cstring>
using std::cout;
using std::endl;
//--------------------------------------------------
// MPI Rank 0 will print to stdout
void print(const char format[],...)
{
  int rank = 0;
  if(!rank) {
    va_list Argp;
    va_start(Argp, format);
    vprintf(format, Argp);
    va_end(Argp);
  }
  return;
}

//--------------------------------------------------
// MPI Rank 0 will print to stdout in red color
void print_error(const char format[],...)
{
  int rank = 0;
  if(!rank) {

    char format_colored[strlen(format)+40] = "";
    strcat(format_colored, "\033[0;31m");
    strcat(format_colored, format);
    strcat(format_colored, "\033[0m");

    va_list Argp;
    va_start(Argp, format);
    vprintf(format_colored, Argp);
    va_end(Argp);
  }
  return;
}

//--------------------------------------------------
// MPI Rank 0 will print to a file
void print(FILE* fd, const char format[],...)
{
  int rank = 0;
  if(!rank) {
    va_list Argp;
    va_start(Argp, format);
    vfprintf(fd, format, Argp);
    va_end(Argp);
  }

  return;
}

//--------------------------------------------------
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string getCurrentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X ", &tstruct);
    if(strlen(tzname[1]) != 0)
      strcat(buf, tzname[1]); //daylight saving time
    else
      strcat(buf, tzname[0]); //standard time


    return buf;
}

//--------------------------------------------------

// Terminate program properly
void exit_mpi()
{
  exit(-1);
}

//--------------------------------------------------










