To create list of files to compile:
```
../../mkmf/bin/list_paths ../../FMS/{mpp,diag_manager,time_manager,include,memutils,constants,platform,fms,random_numbers,mosaic,exchange} ../src/ ../driver/
```
To create a `Makefile`:
```
../../mkmf/bin/mkmf -t ../../mkmf/templates/ncrc-gnu.mk -p bergs.x path_names 
```
To compile a "debug" executable:
```
make DEBUG=1 -j
```
