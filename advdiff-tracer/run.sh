#!/bin/sh

# arguments:
# - movie: create movies (say ./run.sh movie)
# - concerning the script call, there are some arguments, cf. *.f90

# TODO:
# - scheme and bc as run.sh arguments, shared with the fortran script and plotit.m

# clear previous computations
rm -f matlab.output
rm -f nohup.out
rm -f advect1d
rm -f distributions.dat
rm -rf data
rm -rf output
rm -f *.avi

# new data subdirectories
mkdir data
mkdir data/b
mkdir data/f
mkdir data/s
mkdir data/wind
mkdir data/wind/w
mkdir output
mkdir output/3D
mkdir output/map

# compile
echo 'compiling...'
f95 td4.f90 -o advectodiffusion2d

# run
echo 'running...'
./advectodiffusion2d b n

# ok then, either...

# 1. plot?

#xmgrace data/.../d_*

# 2. run matlab, output plots?

# plotit
echo 'matlab processing...'
#/usr/bin/nohup nice -10 en préfixe de la commande suivante, éventuellement
matlab -nodesktop -nosplash -r plotit
echo 'matlab processing done.'

# 3. merge them into a good movie?

for param in $*; do
  if [ "$param" = "movie" ]; then
    echo 'creating movies'
    #echo 'converting from eps to png...'
    ## pseudo-3D plots
    #for f in $(ls output/3D/*.eps); do
      #echo -n $f;
      #echo -n " ";
      #convert -density 100 $f -flatten ${f%.*}.png;
    #done
    ## maps
    #for f in $(ls output/map/*.eps); do
      #echo -n $f;
      #echo -n " ";
      #convert -density 100 $f -flatten ${f%.*}.png;
    #done
    ### ou bien :
    ##ls -1 output/ | while read filename; do
      ##echo $filename;
      ##convert $f -flatten output/${filename%.*}.png;
    ##done
    #echo 'done'

    #echo 'creating movies...'
    #mencoder mf://output/3D/*.tif -mf w=478:h=380:fps=24:type=tif -ovc copy -o animation-3D.mpg
    mencoder mf://output/map/*.png -mf w=478:h=380:fps=24:type=png -ovc copy -o animation-map.avi
    #echo
    #echo 'movies created.'
  fi
done

echo
