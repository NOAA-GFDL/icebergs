#!/bin/bash

#download and extract data for the A68a experiment from the GFDL FTP server:
wget ftp://ftp.gfdl.noaa.gov/perm/Alexander.Huth/a68a_data.tar.gz
tar -xvf a68a_data.tar.gz
rm a68a_data.tar.gz
