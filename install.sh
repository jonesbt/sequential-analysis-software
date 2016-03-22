#!/bin/bash

echo "Installing the dependencies..."
# Determine the operating system.
if [ `uname`='Darwin' ];
then
    # If OS X, use Homebrew as a package manager.
    # Check if Homebrew is installed. If not, notify and install it.
    hash brew 2>/dev/null || {
	read -n 1 -s -p "Installing the Homebrew package manager, press any key to continue or Ctrl-C to quit..."
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    }
    # Install BOOST and GSL.
    read -n 1 -s -p "Installing the BOOST libraries, press any key to continue or Ctrl-C to quit..."
    brew install boost
    read -n 1 -s -p "Installing GSL, press any key to continue or Ctrl-C to quit..."
    brew install gsl
elif [ `uname`="Linux" ];
then
     echo "Support for this operating system is not yet available...aborting"
     exit -1;     
else
    echo "Unsupported operating system. Aborting."
    exit -1;
fi

echo "All of the dependencies have been installed, moving on to the library..."

# Download the source code.
read -n 1 -s -p "Downloading the source code, press any key to continue or Ctrl-C to quit..."
git clone ssh://git@github.com/btjones16/sequential-analysis-software
# Build the library.
read -n 1 -s -p "Building the library, press any key to continue or Ctrl-C to quit..."
cd sequential-analysis-software/src
make lib
read -n 1 -s -p "Building the Python package, press any key to continue or Ctrl-C to quit..."
make python

