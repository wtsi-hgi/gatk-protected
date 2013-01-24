#!/bin/zsh

# This script will update all files in the GATK with the license information
# in the files living in the license directory. This script is meant to be
# run manually (and hopefully only once). To run it, you must be in the
# $GATK root directory. If you want to run it from a different directory
# you will have to update the relative paths on this file
#
# author: Mauricio Carneiro
# date: 1/9/13

ls private/**/*.java  | private/python/UpdateLicense.py
ls private/**/*.scala | private/python/UpdateLicense.py
