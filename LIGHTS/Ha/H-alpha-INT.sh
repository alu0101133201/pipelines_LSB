#!/bin/bash
#
# Pipeline to reduce LBT images.
#
# Copyright (C) 2021 Golini ggolini <ggolini.golini@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.






###########################################################################

############      Brief explanation about this pipeline      ##############

###########################################################################
# This pipeline is created to reduce data from LBT
# cameras (LBC-Blue or LBC-Red) for LIGHTS survey. The usual setup for this
# observations are Sloan-g in the LBC-Blue and Sloan-r filter in the
# LBC-Red. The observation is splitted into 30 pointings of of 180s
# . The FoV of the
# camera is 23' * 25'. The dimension of the entire mosaic depends on the
# dithering pattern that was performed for the specific source.

# The reduction is done ccd after ccd (4 in total) once they are in the
# final corrected step they are taken together to build the final mosaic.

#                The main steps of the pipeline are:
#
# Prepare the images to be processed.
#     and save the dimension of the entire mosaic.
# Images should be float32, bad columns masked, renamed.
#
# Creation of the masterbias with resistant mean (remove INTEGER)
#
# Subtract the masterbias to science images. Bad-saturated pixels masked.
#
# Crop to avoid overscan from the masterbias.
#
# Creation of the ring profile to normalize the images.
#
# Normlization of the science images.
#
# Combination with sigma clipping median to get the first iteration
#     flat.
#
# Divide the cropped images (step 3) by the iteration 1 flat.
#
# Delete temporary images #1.
#
# Run noisechisel over the images from the step 7 to detect sources.
#
# Mask the step 3 images with the step 9 masks.
#
# Normalization of these masked images with the ring.
#
# Combination with sigma clipping median to get the second iteration
#     flat.
#
# Divide the cropped images (step 3) by the iteration 2 flat.
#
# Deleting temporary images #2.
#
# Mask vignetting corners. Remove pixel that have < 0.9 S/N in the flat.
#
# Interpolate saturated pixels to improve astrometry
#
# Downloading the GAIA (edr3,check for new catalogs) catalog for astrometry.
#
# Making the indexes from the GAIA catalog.
#
# Solving the astrometry in the images with SOLVE-FIELD.
#
# Compute distortion by using SExtractor and Scamp and correct
#       of distortion
#
# Mask back aturated pixels
#
# Resample frames into a common grid of the dimension of the final mosaic.
#   here with swarp or warp
#
# Run noisechisel to estimate the sky value.
#
# The sky value is subtracted.
#
# Deleting temporary images #3.
#
# Conversion to physical units (alpha factor = WHT_adu/SDSS_adu)
#      To have every pix value in nanomaggies.
#      This is done using aperture photometry with Gnuastro.
#       2" radius, range of stars between 19 and 22 (g band)
#
# Weight estimate (skystd_best / skystd_i)
#
# Remove outliers out of 3sigma + - mean
#
# First co-add with weighted mean
#
# Create bigger Mask (+ Manual mask around bright objects) and combine with
#      noisechisel mask.
#
# Mask frames with the bigger mask.
#
# Sky subtraction again on resampled individual frames
#
# Second Conversion to physical units (alpha factor = LBT_adu/SDSS_adu)
#
# Second Weight estimate (related to the 2nd sky std)
#
# Remove outliers out of 3sigma + - mean
#
# Final co-add with weighted mean
#




# NECESSARY CONFIG FILES TO WORK
# This script needs some config files to run:
#default.conv
#default.param
#scamp.cfg
#swarp.cfg

# To change the path
#default.sex


# data should be located in the folder DATA
# bias images should be located in the folder bias-"filter"
# the sdss mosaic is also needed: /sdss-mosaic/sdss-"filter".fits
# You must download it from https://dr12.sdss.org/mosaics/
# all the reduction is saved in the folder data-reduction



###########################################################################

###############                    Usage                  #################

###########################################################################
# Although LBC-Blue
# has the Sloan-g filter, the images from this telescope are labeled as
# "lbcb". For this reason, to run the Sloan-g reduction you have to put "b"
# instead of "g". No problems with the r band.
#
# To run this script you have to specify the object and the band to reduce:
#
# ./LBT_red_script.sh OBJ BAND
# ./LBT_red_script.sh ngc1042 r
#





# exit script in case of failure
set -e





# lang to read correctly the awk command
export LANG=C






# Import needed modules to work
# Import needed modules to work
module load scamp/2.10.0
module load gnuastro/0.19.45-0551
module load sextractor/2.25.0





###########################################################################

###############              Previous things              #################

###########################################################################
# This variables are to skip some steps like noisechisel, remove temporary
# files, etc. 1 is for do the steps and 0 for not do they. Note that some
# steps are necessary to continue with the reduction (for instance
# noisechisel steps).

nc_iter1_flat=1       # Noichisel step for iteration 1 flat
nc_iter2_flat=1       # Noichisel step for iteration 2 flat
remove_tmp_1=1
remove_tmp_2=1
remove_tmp_3=1
remove_tmp_4=1
remove_tmp_5=1
remove_tmp_6=1
remove_tmp_7=1
remove_tmp_8=1
gaia_query=1          # Download GAIA catalog for astrometry
mask_corner=1         # Perform corner correction
coadd_it2=1           # Coadd iteration 2.







###########################################################################

###############                 Variables                 #################

###########################################################################
# Each of the 30 pointing *.fits.gz is a cube frame composed of 4 extensions,
# one for each CCD. To rename the images we use
# as first argument the name of the galaxy and also the
# filter of the images that we want to reduce.
# Name of the object (ngc1042, ngc3486, ...)
obj=$1
# Filter (Sloan g, r, i... But just only the band "r" or "g")
c=$2
# ra
ra_gal=$3
# dec
dec_gal=$4

# affecting star
ra_star=40.2705
dec_star=-8.2556





# The telescopes are LBTR and LBTB, and so the images of each telescope
# have a diferent names (lbcr and lbcb). For this reason, this if creates a
# variable which depend on if the image is from LBTR or LBTB and set the
# filter properly.
if test "$c" = "b"; then
    k=g
else
    k=halpha
fi





# SET THE PATH
#path=/Users/giuliagolini/Desktop/int-halpha
path=/scratch/ggolini/1042-"$k"





# DEFINE INPUT DIRECTORIES
INDIR=$path/DATA
DIR=$path/





# Making the project build directory if it doesn't exist
BDIR=$path/data-reduction
if ! [ -d $BDIR ]; then mkdir $BDIR; fi





# DEFINE NUMBER OF POINTINGS
# initial index
start=1
# number of files for the final index (n pointings)
stop=$(ls $INDIR/*.fit | wc -l)






########################################################################
###########                       FUNCTIONS             ################
########################################################################


#1
#########################################################################

#################            Mosaic dimension             ###############

#########################################################################
# Compute the dimension of the mosaic.
# The camera is composed of 4 ccds, 3 vertical (3 2 1) and a horizontal
# one on the top of the three.
# The gap between the ccds is 18" = 80 pixels
# Each ccd is 2304 * 4608 pix (remember we cut the overscan region 52:2096,
# but better keep the entire field of view).
# The camera FoV is 23' * 25'. This dimension is computed in pix.
# The center of the mosaic is the center of the first pointing.
# The coordinates of the pointings are read and saved into a file.
# Since they are in hh:mm:ss then they are converted into degrees.
# Then, compute the offset in x and y of the pointings from the
# center and
# save the pointing with larger distance.

# The final dimensions of the mosaic are:
# X = 2 * (1/2 dim_camera +  (xp-xc)  )
# Y = 2 * (1/2 dim_camera +  (yp-yc)  )


# Center
# The first pointing should be at the center of the coadd
# so the x and y coord of the ccd1 pointing 1 are the reference
#ref=$INDIR/"$obj"_Sloan-"$k"_1_ccd1_all.txt
#xc=$(awk 'NR=='3'{print $1}' $ref)          # degrees
#yc=$(awk 'NR=='3'{print $2}' $ref)          # degrees
# NO use directly the ccords of the galaxy that is easier and more precise.

funct_mosaicdim() {
    xc=$ra_gal
    yc=$dec_gal

    # Now look at all the offsets of the pointings

    # RA offset with to respect to the galaxy center
    echo "# Column 1: X [X ,float,] Column from arithmetic operation " > deltax.txt
    for a in $(seq $start $stop); do
        i=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd1_all.txt
        xp=$(asttable $i  --txtf32format=fixed --txtf32precision=32 | awk '{print $1}' );
        deltax=$(astarithmetic $xp $xc - --quiet)       # delta x
        abs_deltax=$(astarithmetic $deltax abs --quiet) # abs value
        echo $abs_deltax >> deltax.txt                  # save into a file
    done

    # DEC offset with to respect to the galaxy center
    echo "# Column 1: Y [Y ,float,] Column from arithmetic operation " > deltay.txt
    for a in $(seq $start $stop); do
        i=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd1_all.txt
        yp=$(asttable $i  --txtf32format=fixed --txtf32precision=32 | awk '{print $2}');
        deltay=$(astarithmetic $yp $yc - --quiet)       # delta y
        abs_deltay=$(astarithmetic $deltay abs --quiet) # abs value
        echo $abs_deltay >> deltay.txt                  # save into a file
    done
        
    # Camera foV
    # the camera is 23' x 25' then you have to add the offset in x and y to have
    # the entire field. Take 12' * 22.5 and 13'.* 60
    # dimension camera on X
    pix_scale=0.224
    to_arcsec=$(astarithmetic 12. 60. x --quiet)                   # arcmin to arcsec
    to_pix_camera=$(astarithmetic $to_arcsec $pix_scale / --quiet) # arcsec to pix
    # Save the further pointing
    max_deltara=$(asttable  deltax.txt  --sort=X --descending --head=1 --txtf32format=fixed)
    to_arcsec=$(astarithmetic $max_deltara 3600. x --quiet)
    to_pix_offset=$(astarithmetic $to_arcsec $pix_scale / --quiet)
    # sum camera fov and offset
    halfxxdim=$(astarithmetic $to_pix_camera $to_pix_offset + --quiet)
    xxdim=$(astarithmetic $halfxxdim 2. x --quiet)

    # dimension camera on Y
    to_arcsec=$(astarithmetic 13. 60. x --quiet)
    to_pix_camera=$(astarithmetic $to_arcsec $pix_scale / --quiet)
    max_deltadec=$(asttable  deltay.txt  --sort=Y --descending --head=1 --txtf32format=fixed)
    to_arcsec=$(astarithmetic $max_deltadec 3600. x --quiet)
    to_pix_offset=$(astarithmetic $to_arcsec $pix_scale / --quiet)
    # sum camera fov and offset
    halfyydim=$(astarithmetic $to_pix_camera $to_pix_offset + --quiet)
    yydim=$(astarithmetic $halfyydim 2. x --quiet)

    # Save dimensions into a file
    echo $xxdim $yydim > $DIR/mosaic_dim_camera.txt

}



# 2
#########################################################################

#################             Star distances             ###############

#########################################################################
# distance from stars that can affect the flat

# assuming for now that we know the mag of the bright stars that
# could affect the flat and at which distance from the star (and galaxy) we
# have to remove frames to make the flat, how I can compute the real center of
# each single ccd? Because in the header I have only the coords of the pointing
#(center of the camera)
# I can use CRPIX1 and 2 in the header
# CCD2 IS In the center at(1035,2924). (10 pix from the gal)
# CCD1 IS (3157,2924) = 2048 ccd pix + 70 pix of gap from the center
# ccd3 is (-1087,2924) = 2048 ccd pix + 70 pix of gap from the center
# ccd4 is horizontal in the camera! at (-1709,..) so at 1709 + 1035 from the center

funct_diststar(){
    # coords star that couls affect the flat estimation
    ra_star=$1
    dec_star=$2


    # Work for every single ccd
    # ----------------------------------------------------------------------------
    # 1   # frame : pointing + decrease ra of 2048 pix decrease DEC of 10 pix = 2.24"
    # Save the distance from the center of each frame to the star and to the galaxy
    dDIR=$BDIR/distances
    if ! [ -d $dDIR ]; then mkdir $dDIR; fi
    h=1
    distdone=$dDIR/done_ccd"$h".txt
    if [ -f $distdone ]; then
        echo " "
        echo " Distances computed ccd1"
        echo " "
    else
        for a in $(seq 1 $stop); do
            #   READ coords of the frame
            # frame : pointing + decrease DEC of 10 pix = 2.24" = 0.00062 deg
            # decrease ra of 2048 pmt=fixed --txtf32precision=32 x = 0.127 deg
            ra_=$(asttable $INDIR/"$a"_ra_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 );
            ra=$(astarithmetic $ra_ 0.127 - --quiet  ) # floating
            dec_=$(asttable $INDIR/"$a"_dec_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 );
            dec=$(astarithmetic $dec_ 0.00062 - --quiet ) # floating

            # distance from central galaxy
            dist_g_ra=$(astarithmetic $ra $ra_gal - set-d d d x -q)
            dist_g_dec=$(astarithmetic $dec $dec_gal - set-d d d x -q)
            dist_g_arcmin=$(astarithmetic $dist_g_ra $dist_g_dec + sqrt 60 x -q)
            echo $dist_g_arcmin > $dDIR/"$a"_dist-gal_ccd"$h".txt # once
            
            # distance from bright star
            dist_s_ra=$(astarithmetic $ra $ra_star - set-d d d x -q)
            dist_s_dec=$(astarithmetic $dec $dec_star - set-d d d x -q)
            dist_s_arcmin=$(astarithmetic $dist_s_ra $dist_s_dec + sqrt 60 x -q)
            echo $dist_s_arcmin >> $dDIR/"$a"_dist-star_ccd"$h".txt #for every star
             
        done
        echo done > $distdone
    fi




    # ----------------------------------------------------------------------------
    # 2   # frame : pointing + only decrease DEC of 10 pix = 2.24"
    h=2
    distdone=$dDIR/done_ccd"$h".txt
    if [ -f $distdone ]; then
        echo " "
        echo " Distances computed ccd2"
        echo " "
    else
        
        for a in $(seq 1 $stop); do
            # frame : pointing + decrease DEC of 10 pix = 2.24" = 0.00062 deg
            ra=$(asttable $INDIR/"$a"_ra_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 );
            dec_=$(asttable $INDIR/"$a"_dec_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 )
            dec=$(astarithmetic $dec_ 0.00062 - --quiet) # floating
        
            # distance from central galaxy
            dist_g_ra=$(astarithmetic $ra $ra_gal - set-d d d x -q)
            dist_g_dec=$(astarithmetic $dec $dec_gal - set-d d d x -q)
            dist_g_arcmin=$(astarithmetic $dist_g_ra $dist_g_dec + sqrt 60 x -q)
            echo $dist_g_arcmin > $dDIR/"$a"_dist-gal_ccd"$h".txt # once
            
            # distance from bright star
            dist_s_ra=$(astarithmetic $ra $ra_star - set-d d d x -q)
            dist_s_dec=$(astarithmetic $dec $dec_star - set-d d d x -q)
            dist_s_arcmin=$(astarithmetic $dist_s_ra $dist_s_dec + sqrt 60 x -q)
            echo $dist_s_arcmin >> $dDIR/"$a"_dist-star_ccd"$h".txt #for every star
                
        done
        echo done > $distdone
    fi




    # -------------------------------------
    # 3   # frame : pointing + increase ra of 2048 pix decrease DEC of 10 pix = 2.24"

    h=3
    distdone=$dDIR/done_ccd"$h".txt
    if [ -f $distdone ]; then
        echo " "
        echo " Distances computed ccd3"
        echo " "
    else
        for a in $(seq 1 $stop); do
            # frame : pointing + decrease DEC of 10 pix = 2.24" = 0.00062 deg
            # decrease ra of 2048 pix = 0.127 deg
            ra_=$(asttable $INDIR/"$a"_ra_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 );
            ra=$(astarithmetic $ra_ 0.127 + --quiet  ) # floating
        
            dec_=$(asttable $INDIR/"$a"_dec_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 );
            dec=$(astarithmetic $dec_ 0.00062 - --quiet ) # floating
            
            # distance from central galaxy
            dist_g_ra=$(astarithmetic $ra $ra_gal - set-d d d x -q)
            dist_g_dec=$(astarithmetic $dec $dec_gal - set-d d d x -q)
            dist_g_arcmin=$(astarithmetic $dist_g_ra $dist_g_dec + sqrt 60 x -q)
            echo $dist_g_arcmin > $dDIR/"$a"_dist-gal_ccd"$h".txt # once
            
            # distance from bright star
            dist_s_ra=$(astarithmetic $ra $ra_star - set-d d d x -q)
            dist_s_dec=$(astarithmetic $dec $dec_star - set-d d d x -q)
            dist_s_arcmin=$(astarithmetic $dist_s_ra $dist_s_dec + sqrt 60 x -q)
            echo $dist_s_arcmin >> $dDIR/"$a"_dist-star_ccd"$h".txt #for every star

        done
        echo done > $distdone
    fi


    # -------------------------------------
    # 4   # frame : pointing + increase DEC of 1709 + 1024 pix = 2733

    h=4

    distdone=$dDIR/done_ccd"$h".txt

    if [ -f $distdone ]; then
        echo " "
        echo " Distances computed ccd4"
        echo " "
    else
        for a in $(seq 1 $stop); do
            # frame : pointing + decrease DEC of 10 pix = 2.24" = 0.00062 deg
            # increase dec of 2733 pix = 0.170 deg

            ra=$(asttable $INDIR/"$a"_ra_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 );
            dec_=$(asttable $INDIR/"$a"_dec_deg_ccd"$h".txt  --txtf32format=fixed --txtf32precision=32 )
            dec=$(astarithmetic $dec_ 0.170 + --quiet) # floating
            
            # distance from central galaxy
            dist_g_ra=$(astarithmetic $ra $ra_gal - set-d d d x -q)
            dist_g_dec=$(astarithmetic $dec $dec_gal - set-d d d x -q)
            dist_g_arcmin=$(astarithmetic $dist_g_ra $dist_g_dec + sqrt 60 x -q)
            echo $dist_g_arcmin > $dDIR/"$a"_dist-gal_ccd"$h".txt # once
            
            # distance from bright star
            dist_s_ra=$(astarithmetic $ra $ra_star - set-d d d x -q)
            dist_s_dec=$(astarithmetic $dec $dec_star - set-d d d x -q)
            dist_s_arcmin=$(astarithmetic $dist_s_ra $dist_s_dec + sqrt 60 x -q)
            echo $dist_s_arcmin >> $dDIR/"$a"_dist-star_ccd"$h".txt #for every star
                
            
        done
        echo done > $distdone
    fi
    
    echo done > distances_done.txt

}


#3
#########################################################################

#################             Select flat images          ###############

#########################################################################

# Now choose the images that can be used for the flats

funct_select_framesforflat(){
    min_d_star=5
    min_d_gal=5
    # Move into folder for flats only the frames further than a certain distance
    # from the bright sources and the galaxy
    imforflatsdir=$BDIR/images_for_flats
    if ! [ -d $imforflatsdir ]; then mkdir $imforflatsdir; fi
    imforflatdone=$imforflatsdir/done_"$k"_ccd"$h".txt.txt
    if [ -f $imforflatdone ]; then
        echo " "
        echo "Images for flat selected"
        echo " "
    else
        for h in 1 2 3 4; do
            for a in $(seq 1 $stop); do
                g=$dDIR/"$a"_dist-gal_ccd"$h".txt
                d_gal=$(asttable $g --txtf32format=fixed --txtf32precision=32 | awk '{printf("%.16f\n", $1)}')
                if (( $(echo "$d_gal > $min_d_gal" |bc -l) ));then
                    for i in $(seq 1 $n_stars); do
                        # read distances for each star
                        f=$dDIR/"$a"_dist-star_ccd"$h".txt
                        d_star=$(asttable $f --txtf32format=fixed --txtf32precision=32 | awk '{printf("%.16f\n", $1)}')
                            if (( $(echo "$d_star > $min_d_star" |bc -l) )); then
                                cp $INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits $imforflatsdir/
                            fi
                    done
                fi
            done
            # save how many flat images we have for each ccd
            n=$(ls $imforflatsdir/*ccd"$h".fits | wc -l)
            echo $n > $imforflatsdir/n_imforflat_ccd"$h".txt
        done
        echo done > $imforflatdone
    fi

}

# 4
###########################################################################

############                Resampling                      ##############

###########################################################################

funct_resampling_with_SWARP() {
    
    for h in $(seq 1 4);do
    #for h in 1;do
      # swarp into bigger grid and correct of distorsion to improve photometry
      # The grid should have the same dimension of the mask/ final coadd
      cp $sexdir/*Sloan-"$k"_*_ccd$h.head $astroimadir
      # Running swarp over both images and bleeding masks
      swarpcfg=$CDIR/swarp.cfg
      swarpdir=$BDIR/swarp-single-frames-SWARP
      
      swarpdone=$swarpdir/done_"$k"_ccd"$h".txt
      if ! [ -d $swarpdir ]; then mkdir $swarpdir; fi
        cp $sexdir/*Sloan-"$k"_*_ccd$h.head $swarpdir
        if [ -f $swarpdone ]; then
          echo " "
          echo "All the images are already in the same grid for extension $h using SWARP"
          echo " "
        else
            for a in $(seq 1 $n_exp); do
                base="$obj"_Sloan-"$k"_"$a"_ccd"$h"
                i=$astroimadir/"$base".fits
                SWarp -c $swarpcfg $i -CENTER $ra,$dec \
			-IMAGE_SIZE 10508,10508 \
			-IMAGEOUT_NAME $swarpdir/"$base"_m.fits \
			-WEIGHTOUT_NAME $swarpdir/"$base"_m_w.fits \
			-SUBTRACT_BACK N -PIXEL_SCALE 0.224 -PIXELSCALE_TYPE MANUAL
                
                # Mask bad pixels
                astarithmetic $swarpdir/"$base"_m_w.fits -h0 \
			$swarpdir/"$base"_m_w.fits -h0 0 lt nan where -otemp"$base".fits
                # set zero to nan
                # and save as fits  to save space
                astarithmetic $swarpdir/"$base"_m.fits -h0 \
			temp"$base".fits -h1 0 eq nan where -o$swarpdir/"$base".fits
               
                rm -f $swarpdir/"$base"_m.fits
                rm -f $swarpdir/"$base"_m_w.fits
                rm -f temp"$base".fits
               
                
            done
            echo done > $swarpdone
      fi
    done
}





# PIPELINE STARTS HERE
#########################################################################

#################            Prepare data                ###############

#########################################################################

# Preparation of the images:
# -- Transformation to float32 fits format (avoid bad arithmetic)
# -- Save each single ccd frame extension with index name (help find outliers)
# -- Save ra and dec for ech pointing.
# -- Save file with the regions of the pointings. (dithering pattern)
# -- Compute the dimension of the final mosaic





            

# Since I have 4 extension I need 2 loops; one for each ccd
# one for each pointing.
renamedone=$INDIR/renamedone_.txt
if [ -f $renamedone ]; then
  echo " "
  echo "Science images are already renamed"
  echo " "
else
    for h in $(seq 1 4); do
        a=$start
        for i in $INDIR/*.fit; do
            tmp=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_integer.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            astfits $i --copy="$h" -o$tmp            # copy extension
            astarithmetic $tmp -h1 float32 -o$out    # Convert to floating point32
            a=$((a+1));                               # increase index
        done
    done
    echo done > $renamedone
fi




#########################################################################

#################            Bad Columns                  ###############

#########################################################################

# -----------------------------------------------------------------------------
# Mask bad columns CCD1
if [ "$k" == g ]; then
    h=1
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            # Vertical bad columns CCD1
            astcrop $input --mode=img --polygon=632,0:632,4605:634,4605:634,0 --polygonout -otmp.fits
            astcrop tmp.fits --mode=img --polygon=752,0:752,4605:754,4605:754,0 --polygonout -otmp2.fits
            astcrop tmp2.fits --mode=img --polygon=1100,0:1100,4605:1102,4605:1102,0 --polygonout -otmp3.fits
            astcrop tmp3.fits --mode=img --polygon=1117,0:1117,4605:1120,4605:1120,0 --polygonout -otmp4.fits
            astcrop tmp4.fits --mode=img --polygon=1434,0:1434,4605:1436,4605:1436,0 --polygonout -otmp5.fits
            astcrop tmp5.fits --mode=img --polygon=1698,0:1698,4605:1700,4605:1700,0 --polygonout -otmp6.fits
            astcrop tmp6.fits --mode=img --polygon=111,0:111,4605:114,4605:114,0 --polygonout -o$out
            # Remove temporary files
            rm tmp.fits tmp2.fits tmp3.fits tmp4.fits tmp5.fits tmp6.fits
        done
        echo done > $badcolumnsdone
    fi
else
    h=1
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            mv $input $out
        done
        echo done > $badcolumnsdone
    fi
fi


# Mask bad columns CCD2
if [ "$k" == g ]; then
    h=2
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            astcrop $input --mode=img --polygon=1879,0:1879,4605:1883,4605:1883,0 --polygonout -o$out
        done
        echo done > $badcolumnsdone
    fi
else
    h=2
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            mv $input $out
        done
        echo done > $badcolumnsdone
    fi
fi
    
    
    
    
    
# Mask bad columns CCD3
if [ "$k" == g ]; then
    h=3
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            astcrop $input --mode=img --polygon=1670,0:1670,4605:1675,4605:1675,0 --polygonout -o$out
        done
        echo done > $badcolumnsdone
    fi
else
    h=3
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            mv $input $out
        done
        echo done > $badcolumnsdone
    fi
fi





# Mask bad columns CCD4
if [ "$k" == g ]; then
    h=4
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            astcrop $input --mode=img --polygon=80,0:80,4605:83,4605:83,0 --polygonout -otmp.fits
            astcrop tmp.fits --mode=img --polygon=94,0:94,4605:96,4605:96,0 --polygonout -o$out
            rm tmp.fits
        done
        echo done > $badcolumnsdone
    fi
else
    h=4
    badcolumnsdone=$INDIR/done_badcolumns_"$h".txt
    if [ -f $badcolumnsdone ]; then
      echo " "
      echo "Bad columns masked for ccd $h"
      echo " "
    else
        for a in $(seq 1 $stop); do
            input=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h"_e.fits
            out=$INDIR/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            mv $input $out
        done
        echo done > $badcolumnsdone
    fi
fi


# coord pointing in file raw e keywords ra and dec
# sarebbe da fare il loop





#########################################################################

#################         Compute dimension Mosaic        ###############

#########################################################################
# compute dimensions of the mosaic
#if ! [ -f $DIR/mosaic_dim_camera.txt ]; then
#    funct_mosaicdim
#fi
#echo " "
#echo " Mosaic dimension computed"
#echo " "
mos_dim_x=11001
mos_dim_y=11001
#########################################################################

#################         Compute distances               ###############

#########################################################################
# compute distances from bright star and galaxy

#funct_diststar $ra_star $dec_star
#echo " "
#echo " Distances computed for each pointing"
#echo " "


#########################################################################

#################         Images for flat                ###############

#########################################################################

# select frames for flat
#funct_select_framesforflat
#echo " "
#echo "Images for flat selected"
#echo " "





# -------------------------------------------------------
# Now the data are prepared, the data reduction can start.
# These variables are the folders where the data and the results were being
# stored. Set, following your preferences.
# Set the build and the data directory

INDIR=$path
SDIR=$path
CDIR=$SDIR/config

# Get the coordinates of the first point
ra=$ra_gal
dec=$dec_gal

# Get the total number of exposures per ccd for the given filter
n_exp=$(ls $INDIR/DATA/*.fit | wc -l)
sat_level=50000         # saturation level

# Basic parameters for the photometry
SDSS_pix=0.396      # Pixel scales
lbt_pix=0.333       # Pixel scales
zp_sdss=22.5        # Reference zeropoint

#mos=$INDIR/sdss-mosaic/sdss-"$k".fits   # SDSS mosaic image of the galaxy


m_min=17    # Range of magnitude for the photometry
m_max=22

# Why are you choosing these values?
# 2" aperture photometry for both telescopes. Why aperture with sextractor and not noisechisel? Because noisechisel catalog is affected by the noise of the sky and we want to be sure to take always the same stars and the same area to compute the flux. This is the range of sdss magnitude that has to be chosen in order to not be affected of neither faint stars (too faint for sdss catalog to be correct) nor very bright stars (saturated in LBT).
# The value 2"x 2" is because it is big enough for the seeing and little enough to avoid contamination from other sources nearby.
# min max values for magnitudes of stars for photometry

# Select the dimension of the aperture for the aperture photometry(RADIUS)
d_aperture=2.5                # arcsec


# Covert into pixels the diameter
r_sdss_pix_=$(astarithmetic $d_aperture $SDSS_pix / --quiet |awk '{print $1}')
r_sdss_pix=${r_sdss_pix_/\.*}       # to convert in integer
r_lbt_pix_=$(astarithmetic $d_aperture $lbt_pix / --quiet |awk '{print $1}' )
r_lbt_pix=${r_lbt_pix_/\.*}
echo "number of pixels for the aperture photometry"
echo $r_lbt_pix $r_sdss_pix         # to convert in integer


# Read the dimension of the mosaic
#mos_dim_x=$(awk 'NR=='1'{printf("%.16f\n", $1)}' $INDIR/mosaic_dim_camera.txt)
#mos_dim_y=$(awk 'NR=='1'{printf("%.16f\n", $1)}' $INDIR/mosaic_dim_camera.txt)





echo " Start with the pipeline"

# Here the loop for all the ccds begins.
for h in $(seq 1 4); do

  #########################################################################

  #################         Creating the master bias        ###############

  #########################################################################
  mbiasdir=$BDIR/masterbias
  mbiasdone=$mbiasdir/mbias_"$k"_ccd"$h".fits
  if ! [ -d $mbiasdir ]; then mkdir $mbiasdir; fi
  if [ -f $mbiasdone ]; then
    echo " "
    echo "Masterbias is already done for extension $h"
    echo " "
  else
    astarithmetic $(ls $INDIR/bias-"$k"/* ) \
                  $(ls $INDIR/bias-"$k"/* | wc -l) \
                  3 0.2 sigclip-mean -g$h \
                  -o $mbiasdir/mbias_"$k"_ccd$h.fits
  fi
  
  
  # Save airmass
  skydir=$BDIR/airmass-analysis
  skydone=$skydir/done_.txt
  if ! [ -d $skydir ]; then mkdir $skydir; fi
  if [ -f $skydone ]; then
    echo " "
    echo "Airmass already saved"
    echo " "
  else
    for i in $(ls $INDIR/DATA/*.fit ); do
      air=$(astfits $i -h0 --keyvalue=AIRMASS | awk '{print $2}')
      UT=$(astfits $i -h0 --keyvalue=UTC_OBS | awk '{print $2}')
      AZ=$(astfits $i -h0 --keyvalue=TELAZ | awk '{print $2}')
      echo $air $UT $AZ >> $skydir/airmass.txt
    done
    echo done > $skydone
  fi
  
  
  #########################################################################

  #################        Subtract master bias            ###############

  #########################################################################
  
  # Now using indexes that could be better
  # Substracting mbias to science images.
  # and also from images for flat.
  # Also a counter variable is
  # created to rename the images. Bad and saturated pixels are masked.
  mbiascorrdir=$BDIR/bias-corrected
  mbiascorrdone=$mbiascorrdir/done_"$k"_ccd"$h".txt
  if ! [ -d $mbiascorrdir ]; then mkdir $mbiascorrdir; fi
  if [ -f $mbiascorrdone ]; then
    echo " "
    echo "Science images are already bias corrected for extension $h"
    echo " "
  else
    for a in $(seq 1 $n_exp); do
      base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
      i=$INDIR/DATA/$base
      astarithmetic $i -h1 set-i $mbiasdir/mbias_"$k"_ccd"$h".fits -h1 set-m \
                i i $sat_level gt i isblank or 2 dilate 2 dilate nan where m -    \
                --wcsfile=none -o $mbiascorrdir/$base
      
     done
     echo done > $mbiascorrdone
  fi
  
  
  
  

  #########################################################################

  #################        Remove overscan                  ###############

  #########################################################################

  # The images must be cropped due to and a kind of overscan zone in master
  # bias image.
  cropdir=$BDIR/cropped
  cropdone=$cropdir/done_"$k"_ccd"$h".txt
  if ! [ -d $cropdir ]; then mkdir $cropdir; fi
  
  if [ -f $cropdone ]; then
    echo " "
    echo "Science images overscan corrected for extension $h"
    echo " "
  else
    for a in $(seq 1 $n_exp); do
      base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
      i=$mbiascorrdir/$base
      astcrop $i -h1 --mode=img --section=52:2000,0:4001 -o $cropdir/$base
      
    done
    echo done > $cropdone
  fi

  
  

  #########################################################################

  #################         Creating the ring mask          ###############

  #########################################################################
  ringdir=$BDIR/ring
  ringdone=$ringdir/ring_ccd"$h".fits
  if ! [ -d $ringdir ]; then mkdir $ringdir; fi
  if [ -f $ringdone ]; then
    echo " "
    echo "Ring mask are already created for extension $h"
    echo " "
  else
    astmkprof --background=$cropdir/"$obj"_Sloan-"$k"_1_ccd"$h".fits -h1 --mforflatpix --mode=img --type=uint8 --circumwidth=200 --clearcanvas -o $ringdir/ring_ccd"$h".fits $CDIR/flat_ring_ccd"$h".txt
  fi
  
  # NOW SINCE I WANT 2 FLATS I DIVIDE THIS STEP INTO 2
  
  # FIRST 10
  #########################################################################

  ##########         Creating the it1 master flat image         ###########

  #########################################################################
  # Creating iteration 1 flat_it1. First we need to normalise the science
  # images.
  normit1dir=$BDIR/norm-it1-images-firstflat
  normit1done=$normit1dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit1dir ]; then mkdir $normit1dir; fi
  if [ -f $normit1done ]; then
    echo " "
    echo "Science images are already normalized for extension $h"
    echo " "
  else
    for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit1dir/$base;
    done
    echo done > $normit1done
  fi

  


  # Then we combine the normalize images with a sigma clipping median
  flatit1dir=$BDIR/flat-it1-first-flat
  flatit1done=$flatit1dir/flat-it1_"$k"_ccd"$h".fits
  if ! [ -d $flatit1dir ]; then mkdir $flatit1dir; fi
  if [ -f $flatit1done ]; then
    echo " "
    echo "Science images are stacked for it1 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits -1) \
                  $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                  -o $flatit1dir/flat-it1_"$k"_ccd"$h".fits
  
  fi

  


  # Dividing the science image by the it1 flat
  flatit1imadir=$BDIR/flat-it1-ima-firstflat
  flatit1imadone=$flatit1imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
  if [ -f $flatit1imadone ]; then
    echo " "
    echo "Science images are divided by the it1 flat for extension $h"
    echo " "
  else
    for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit1dir/flat-it1_"$k"_ccd"$h".fits -h1 / -o $flatit1imadir/$base
    done
    echo done > $flatit1imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_1 = 1 ]; then
    #rm -f $mbiascorrdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $normit1dir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files $h removed"
    echo " "
  fi



  
  #########################################################################

  ##########         Creating the it2 master flat image         ###########

  #########################################################################
  # Now let's make the second iteration of the flat. First we need to apply
  # noisechisel to flat corrected images. We found that pseudo detection
  # S/N in all the images are around 5.2 (you can see it in the output of
  # noisechisel). With this conditional you can skip this step simply
  # changing true to false.
  noisechisel_param="--tilesize=50,50 \
                     --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  if [ $nc_iter1_flat = 1 ]; then
    noiseit2dir=$BDIR/noise-it2-firstflat
    noiseit2done=$noiseit2dir/done_"$k"_ccd"$h".txt
    if ! [ -d $noiseit2dir ]; then mkdir $noiseit2dir; fi
    if [ -f $noiseit2done ]; then
      echo " "
      echo "Science images are 'noisechiseled' for it2 flat for extension $h"
      echo " "
    else
      for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit1imadir/$base
        astnoisechisel $i $noisechisel_param -o $noiseit2dir/$base
      done
      echo done > $noiseit2done
    fi
  fi


 
 
  # Substituting ones in flat corrected images image for nan values. WHY?MASKING!
  maskedit2dir=$BDIR/masked-it2-firstflat
  maskedit2done=$maskedit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $maskedit2dir ]; then mkdir $maskedit2dir; fi
  if [ -f $maskedit2done ]; then
    echo " "
    echo "Science images are masked for extension $h"
    echo " "
  else
    for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $noiseit2dir/$base -h1 1 eq nan where -o $maskedit2dir/$base;
    done
    echo done > $maskedit2done
  fi




  # Normalizing masked images
  normit2dir=$BDIR/norm-it2-images-firstflat
  normit2done=$normit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit2dir ]; then mkdir $normit2dir; fi
  if [ -f $normit2done ]; then
    echo " "
    echo "Masked science images are normalized for extension $h"
    echo " "
  else
    for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedit2dir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit2dir/$base;
    done
    echo done > $normit2done
  fi




  # Combining masked normalized images to make it2 flat
  flatit2dir=$BDIR/flat-it2-firstflat
  flatit2done=$flatit2dir/flat-it2_"$k"_ccd"$h".fits
  if ! [ -d $flatit2dir ]; then mkdir $flatit2dir; fi
  if [ -f $flatit2done ]; then
    echo " "
    echo "Science images are stacked for it2 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits) \
                  $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                 -o $flatit2dir/flat-it2_"$k"_ccd"$h".fits
  fi


  # Dividing the science image by the it2 flat
  flatit2imadir=$BDIR/flat-it2-ima-firstflat
  flatit2imadone=$flatit2imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
  if [ -f $flatit2imadone ]; then
    echo " "
    echo "Science images are divided by the it2 flat for extension $h"
    echo " "
  else
    for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 / -o $flatit2imadir/$base
    done
    echo done > $flatit2imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_2 = 0 ]; then
    #rm -f $cropdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $noiseit2dir/*.fits
    #rm -f $maskedit2dir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $normit2dir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $flatit1imadir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files #2 removed"
    echo " "
  fi




  #########################################################################

  ##########            Masking the vignetting zones            ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdir=$BDIR/masked-corner-firstflat
    maskedcornerdone=$maskedcornerdir/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
    if [ -f $maskedcornerdone ]; then
      echo " "
      echo "Corners are already masked for extension $h"
      echo " "
    else
      for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit2imadir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 0.9 lt nan where -o $maskedcornerdir/$base;
      done
      echo done > $maskedcornerdone
    fi
  fi
  
  #########################################################################

  ##########            Masking remnant bad columns           ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdirbad=$BDIR/masked-corner-firstflat-badcolumns
    maskedcornerdonebad=$maskedcornerdirbad/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdirbad ]; then mkdir $maskedcornerdirbad; fi
    if [ -f $maskedcornerdonebad ]; then
      echo " "
      echo "Corners are already masked for extension $h bad"
      echo " "
    else
      for a in $(seq 1 29); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedcornerdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 2. gt nan where -o $maskedcornerdirbad/$base;
      done
      echo done > $maskedcornerdonebad
    fi
  fi
  
  
  
  
  # SECOND 10
  #########################################################################

  ##########         Creating the it1 master flat image         ###########

  #########################################################################
  # Creating iteration 1 flat_it1. First we need to normalise the science
  # images.
  normit1dir=$BDIR/norm-it1-images-secondflat
  normit1done=$normit1dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit1dir ]; then mkdir $normit1dir; fi
  if [ -f $normit1done ]; then
    echo " "
    echo "Science images are already normalized for extension $h"
    echo " "
  else
    for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit1dir/$base;
    done
    echo done > $normit1done
  fi

  


  # Then we combine the normalize images with a sigma clipping median
  flatit1dir=$BDIR/flat-it1-secondflat
  flatit1done=$flatit1dir/flat-it1_"$k"_ccd"$h".fits
  if ! [ -d $flatit1dir ]; then mkdir $flatit1dir; fi
  if [ -f $flatit1done ]; then
    echo " "
    echo "Science images are stacked for it1 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits -1) \
                  $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                  -o $flatit1dir/flat-it1_"$k"_ccd"$h".fits
  
  fi

  


  # Dividing the science image by the it1 flat
  flatit1imadir=$BDIR/flat-it1-ima-secondflat
  flatit1imadone=$flatit1imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
  if [ -f $flatit1imadone ]; then
    echo " "
    echo "Science images are divided by the it1 flat for extension $h"
    echo " "
  else
    for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit1dir/flat-it1_"$k"_ccd"$h".fits -h1 / -o $flatit1imadir/$base
    done
    echo done > $flatit1imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_1 = 1 ]; then
    #rm -f $mbiascorrdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $normit1dir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files $h removed"
    echo " "
  fi



  
  #########################################################################

  ##########         Creating the it2 master flat image         ###########

  #########################################################################
  # Now let's make the second iteration of the flat. First we need to apply
  # noisechisel to flat corrected images. We found that pseudo detection
  # S/N in all the images are around 5.2 (you can see it in the output of
  # noisechisel). With this conditional you can skip this step simply
  # changing true to false.
  noisechisel_param="--tilesize=50,50 \
                     --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  if [ $nc_iter1_flat = 1 ]; then
    noiseit2dir=$BDIR/noise-it2-secondflat
    noiseit2done=$noiseit2dir/done_"$k"_ccd"$h".txt
    if ! [ -d $noiseit2dir ]; then mkdir $noiseit2dir; fi
    if [ -f $noiseit2done ]; then
      echo " "
      echo "Science images are 'noisechiseled' for it2 flat for extension $h"
      echo " "
    else
      for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit1imadir/$base
        astnoisechisel $i $noisechisel_param -o $noiseit2dir/$base
      done
      echo done > $noiseit2done
    fi
  fi


 
 
  # Substituting ones in flat corrected images image for nan values. WHY?MASKING!
  maskedit2dir=$BDIR/masked-it2-secondflat
  maskedit2done=$maskedit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $maskedit2dir ]; then mkdir $maskedit2dir; fi
  if [ -f $maskedit2done ]; then
    echo " "
    echo "Science images are masked for extension $h"
    echo " "
  else
    for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $noiseit2dir/$base -h1 1 eq nan where -o $maskedit2dir/$base;
    done
    echo done > $maskedit2done
  fi




  # Normalizing masked images
  normit2dir=$BDIR/norm-it2-images-secondflat
  normit2done=$normit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit2dir ]; then mkdir $normit2dir; fi
  if [ -f $normit2done ]; then
    echo " "
    echo "Masked science images are normalized for extension $h"
    echo " "
  else
    for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedit2dir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit2dir/$base;
    done
    echo done > $normit2done
  fi




  # Combining masked normalized images to make it2 flat
  flatit2dir=$BDIR/flat-it2-secondflat
  flatit2done=$flatit2dir/flat-it2_"$k"_ccd"$h".fits
  if ! [ -d $flatit2dir ]; then mkdir $flatit2dir; fi
  if [ -f $flatit2done ]; then
    echo " "
    echo "Science images are stacked for it2 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits) \
                  $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                 -o $flatit2dir/flat-it2_"$k"_ccd"$h".fits
  fi


  # Dividing the science image by the it2 flat
  flatit2imadir=$BDIR/flat-it2-ima-secondflat
  flatit2imadone=$flatit2imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
  if [ -f $flatit2imadone ]; then
    echo " "
    echo "Science images are divided by the it2 flat for extension $h"
    echo " "
  else
    for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 / -o $flatit2imadir/$base
    done
    echo done > $flatit2imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_2 = 1 ]; then
    #rm -f $cropdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $noiseit2dir/*.fits
    #rm -f $maskedit2dir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $normit2dir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $flatit1imadir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files #2 removed"
    echo " "
  fi




  #########################################################################

  ##########            Masking the vignetting zones            ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdir=$BDIR/masked-corner-secondflat
    maskedcornerdone=$maskedcornerdir/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
    if [ -f $maskedcornerdone ]; then
      echo " "
      echo "Corners are already masked for extension $h"
      echo " "
    else
      for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit2imadir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 0.9 lt nan where -o $maskedcornerdir/$base;
      done
      echo done > $maskedcornerdone
    fi
  fi
  
  
  
  #########################################################################

  ##########            Masking remnant bad columns           ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdirbad=$BDIR/masked-corner-secondflat-badcolumns
    maskedcornerdonebad=$maskedcornerdirbad/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdirbad ]; then mkdir $maskedcornerdirbad; fi
    if [ -f $maskedcornerdonebad ]; then
      echo " "
      echo "Corners are already masked for extension $h bad"
      echo " "
    else
      for a in $(seq 30 49); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedcornerdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 2. gt nan where -o $maskedcornerdirbad/$base;
      done
      echo done > $maskedcornerdonebad
    fi
  fi
  
  # ---------third night

  #########################################################################

  ##########         Creating the it1 master flat image         ###########

  #########################################################################
  # Creating iteration 1 flat_it1. First we need to normalise the science
  # images.
  normit1dir=$BDIR/norm-it1-images-thirdflat
  normit1done=$normit1dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit1dir ]; then mkdir $normit1dir; fi
  if [ -f $normit1done ]; then
    echo " "
    echo "Science images are already normalized for extension $h"
    echo " "
  else
    for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit1dir/$base;
    done
    echo done > $normit1done
  fi

  


  # Then we combine the normalize images with a sigma clipping median
  flatit1dir=$BDIR/flat-it1-thirdflat
  flatit1done=$flatit1dir/flat-it1_"$k"_ccd"$h".fits
  if ! [ -d $flatit1dir ]; then mkdir $flatit1dir; fi
  if [ -f $flatit1done ]; then
    echo " "
    echo "Science images are stacked for it1 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits -1) \
                  $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                  -o $flatit1dir/flat-it1_"$k"_ccd"$h".fits
  
  fi

  


  # Dividing the science image by the it1 flat
  flatit1imadir=$BDIR/flat-it1-ima-thirdflat
  flatit1imadone=$flatit1imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
  if [ -f $flatit1imadone ]; then
    echo " "
    echo "Science images are divided by the it1 flat for extension $h"
    echo " "
  else
    for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit1dir/flat-it1_"$k"_ccd"$h".fits -h1 / -o $flatit1imadir/$base
    done
    echo done > $flatit1imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_1 = 1 ]; then
    #rm -f $mbiascorrdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $normit1dir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files $h removed"
    echo " "
  fi



  
  #########################################################################

  ##########         Creating the it2 master flat image         ###########

  #########################################################################
  # Now let's make the second iteration of the flat. First we need to apply
  # noisechisel to flat corrected images. We found that pseudo detection
  # S/N in all the images are around 5.2 (you can see it in the output of
  # noisechisel). With this conditional you can skip this step simply
  # changing true to false.
  noisechisel_param="--tilesize=50,50 \
                     --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  if [ $nc_iter1_flat = 1 ]; then
    noiseit2dir=$BDIR/noise-it2-thirdflat
    noiseit2done=$noiseit2dir/done_"$k"_ccd"$h".txt
    if ! [ -d $noiseit2dir ]; then mkdir $noiseit2dir; fi
    if [ -f $noiseit2done ]; then
      echo " "
      echo "Science images are 'noisechiseled' for it2 flat for extension $h"
      echo " "
    else
      for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit1imadir/$base
        astnoisechisel $i $noisechisel_param -o $noiseit2dir/$base
      done
      echo done > $noiseit2done
    fi
  fi


 
 
  # Substituting ones in flat corrected images image for nan values. WHY?MASKING!
  maskedit2dir=$BDIR/masked-it2-thirdflat
  maskedit2done=$maskedit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $maskedit2dir ]; then mkdir $maskedit2dir; fi
  if [ -f $maskedit2done ]; then
    echo " "
    echo "Science images are masked for extension $h"
    echo " "
  else
    for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $noiseit2dir/$base -h1 1 eq nan where -o $maskedit2dir/$base;
    done
    echo done > $maskedit2done
  fi




  # Normalizing masked images
  normit2dir=$BDIR/norm-it2-images-thirdflat
  normit2done=$normit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit2dir ]; then mkdir $normit2dir; fi
  if [ -f $normit2done ]; then
    echo " "
    echo "Masked science images are normalized for extension $h"
    echo " "
  else
    for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedit2dir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit2dir/$base;
    done
    echo done > $normit2done
  fi




  # Combining masked normalized images to make it2 flat
  flatit2dir=$BDIR/flat-it2-thirdflat
  flatit2done=$flatit2dir/flat-it2_"$k"_ccd"$h".fits
  if ! [ -d $flatit2dir ]; then mkdir $flatit2dir; fi
  if [ -f $flatit2done ]; then
    echo " "
    echo "Science images are stacked for it2 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits) \
                  $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                 -o $flatit2dir/flat-it2_"$k"_ccd"$h".fits
  fi


  # Dividing the science image by the it2 flat
  flatit2imadir=$BDIR/flat-it2-ima-thirdflat
  flatit2imadone=$flatit2imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
  if [ -f $flatit2imadone ]; then
    echo " "
    echo "Science images are divided by the it2 flat for extension $h"
    echo " "
  else
    for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 / -o $flatit2imadir/$base
    done
    echo done > $flatit2imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_2 = 1 ]; then
    #rm -f $cropdir/*Sloan-"$k"*_ccd"$h".fits
    rm -f $noiseit2dir/*.fits
    rm -f $maskedit2dir/*Sloan-"$k"*_ccd"$h".fits
    rm -f $normit2dir/*Sloan-"$k"*_ccd"$h".fits
    rm -f $flatit1imadir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files #2 removed"
    echo " "
  fi




  #########################################################################

  ##########            Masking the vignetting zones            ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdir=$BDIR/masked-corner-thirdflat
    maskedcornerdone=$maskedcornerdir/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
    if [ -f $maskedcornerdone ]; then
      echo " "
      echo "Corners are already masked for extension $h"
      echo " "
    else
      for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit2imadir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 0.9 lt nan where -o $maskedcornerdir/$base;
      done
      echo done > $maskedcornerdone
    fi
  fi
  
  
  
  #########################################################################

  ##########            Masking remnant bad columns           ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdirbad=$BDIR/masked-corner-thirdflat-badcolumns
    maskedcornerdonebad=$maskedcornerdirbad/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdirbad ]; then mkdir $maskedcornerdirbad; fi
    if [ -f $maskedcornerdonebad ]; then
      echo " "
      echo "Corners are already masked for extension $h bad"
      echo " "
    else
      for a in $(seq 50 82); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedcornerdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 2. gt nan where -o $maskedcornerdirbad/$base;
      done
      echo done > $maskedcornerdonebad
    fi
  fi
      
    # ---------four night

  #########################################################################

  ##########         Creating the it1 master flat image         ###########

  #########################################################################
  # Creating iteration 1 flat_it1. First we need to normalise the science
  # images.
  normit1dir=$BDIR/norm-it1-images-fourflat
  normit1done=$normit1dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit1dir ]; then mkdir $normit1dir; fi
  if [ -f $normit1done ]; then
    echo " "
    echo "Science images are already normalized for extension $h"
    echo " "
  else
    for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit1dir/$base;
    done
    echo done > $normit1done
  fi

  


  # Then we combine the normalize images with a sigma clipping median
  flatit1dir=$BDIR/flat-it1-fourflat
  flatit1done=$flatit1dir/flat-it1_"$k"_ccd"$h".fits
  if ! [ -d $flatit1dir ]; then mkdir $flatit1dir; fi
  if [ -f $flatit1done ]; then
    echo " "
    echo "Science images are stacked for it1 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits -1) \
                  $(ls $normit1dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                  -o $flatit1dir/flat-it1_"$k"_ccd"$h".fits
  
  fi

  


  # Dividing the science image by the it1 flat
  flatit1imadir=$BDIR/flat-it1-ima-fourflat
  flatit1imadone=$flatit1imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit1imadir ]; then mkdir $flatit1imadir; fi
  if [ -f $flatit1imadone ]; then
    echo " "
    echo "Science images are divided by the it1 flat for extension $h"
    echo " "
  else
    for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit1dir/flat-it1_"$k"_ccd"$h".fits -h1 / -o $flatit1imadir/$base
    done
    echo done > $flatit1imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_1 = 1 ]; then
    #rm -f $mbiascorrdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $normit1dir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files $h removed"
    echo " "
  fi



  
  #########################################################################

  ##########         Creating the it2 master flat image         ###########

  #########################################################################
  # Now let's make the second iteration of the flat. First we need to apply
  # noisechisel to flat corrected images. We found that pseudo detection
  # S/N in all the images are around 5.2 (you can see it in the output of
  # noisechisel). With this conditional you can skip this step simply
  # changing true to false.
  noisechisel_param="--tilesize=50,50 \
                     --minskyfrac=0.9 \
                     --meanmedqdiff=0.01 \
                     --snthresh=5.2 \
                     --detgrowquant=0.7 \
                     --detgrowmaxholesize=1000 \
                     --rawoutput"
  if [ $nc_iter1_flat = 1 ]; then
    noiseit2dir=$BDIR/noise-it2-fourflat
    noiseit2done=$noiseit2dir/done_"$k"_ccd"$h".txt
    if ! [ -d $noiseit2dir ]; then mkdir $noiseit2dir; fi
    if [ -f $noiseit2done ]; then
      echo " "
      echo "Science images are 'noisechiseled' for it2 flat for extension $h"
      echo " "
    else
      for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit1imadir/$base
        astnoisechisel $i $noisechisel_param -o $noiseit2dir/$base
      done
      echo done > $noiseit2done
    fi
  fi


 
 
  # Substituting ones in flat corrected images image for nan values. WHY?MASKING!
  maskedit2dir=$BDIR/masked-it2-fourflat
  maskedit2done=$maskedit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $maskedit2dir ]; then mkdir $maskedit2dir; fi
  if [ -f $maskedit2done ]; then
    echo " "
    echo "Science images are masked for extension $h"
    echo " "
  else
    for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $noiseit2dir/$base -h1 1 eq nan where -o $maskedit2dir/$base;
    done
    echo done > $maskedit2done
  fi




  # Normalizing masked images
  normit2dir=$BDIR/norm-it2-images-fourflat
  normit2done=$normit2dir/done_"$k"_ccd"$h".txt
  if ! [ -d $normit2dir ]; then mkdir $normit2dir; fi
  if [ -f $normit2done ]; then
    echo " "
    echo "Masked science images are normalized for extension $h"
    echo " "
  else
    for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedit2dir/$base
        me=$(astarithmetic $i -h1 $ringdir/ring_ccd"$h".fits -h1 0 eq nan where medianvalue --quiet)
        astarithmetic $i -h1 $me / -o $normit2dir/$base;
    done
    echo done > $normit2done
  fi




  # Combining masked normalized images to make it2 flat
  flatit2dir=$BDIR/flat-it2-fourflat
  flatit2done=$flatit2dir/flat-it2_"$k"_ccd"$h".fits
  if ! [ -d $flatit2dir ]; then mkdir $flatit2dir; fi
  if [ -f $flatit2done ]; then
    echo " "
    echo "Science images are stacked for it2 flat for extension $h"
    echo " "
  else
    astarithmetic $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits) \
                  $(ls $normit2dir/*Sloan-"$k"*_ccd"$h".fits | wc -l) \
                  3 10 sigclip-median -g1 \
                 -o $flatit2dir/flat-it2_"$k"_ccd"$h".fits
  fi


  # Dividing the science image by the it2 flat
  flatit2imadir=$BDIR/flat-it2-ima-fourflat
  flatit2imadone=$flatit2imadir/done_"$k"_ccd"$h".txt
  if ! [ -d $flatit2imadir ]; then mkdir $flatit2imadir; fi
  if [ -f $flatit2imadone ]; then
    echo " "
    echo "Science images are divided by the it2 flat for extension $h"
    echo " "
  else
    for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$cropdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 / -o $flatit2imadir/$base
    done
    echo done > $flatit2imadone
  fi




  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_2 = 1 ]; then
    #rm -f $cropdir/*Sloan-"$k"*_ccd"$h".fits
    rm -f $noiseit2dir/*.fits
    rm -f $maskedit2dir/*Sloan-"$k"*_ccd"$h".fits
    rm -f $normit2dir/*Sloan-"$k"*_ccd"$h".fits
    rm -f $flatit1imadir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files #2 removed"
    echo " "
  fi




  #########################################################################

  ##########            Masking the vignetting zones            ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdir=$BDIR/masked-corner-fourflat
    maskedcornerdone=$maskedcornerdir/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
    if [ -f $maskedcornerdone ]; then
      echo " "
      echo "Corners are already masked for extension $h"
      echo " "
    else
      for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$flatit2imadir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 0.9 lt nan where -o $maskedcornerdir/$base;
      done
      echo done > $maskedcornerdone
    fi
  fi
  
  
  
  #########################################################################

  ##########            Masking remnant bad columns           ###########

  #########################################################################
  # Enmascarando las esquinas
  if [ $mask_corner = 1 ]; then
    maskedcornerdirbad=$BDIR/masked-corner-fourflat-badcolumns
    maskedcornerdonebad=$maskedcornerdirbad/done_"$k"_ccd"$h".txt
    if ! [ -d $maskedcornerdirbad ]; then mkdir $maskedcornerdirbad; fi
    if [ -f $maskedcornerdonebad ]; then
      echo " "
      echo "Corners are already masked for extension $h bad"
      echo " "
    else
      for a in $(seq 83 106); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedcornerdir/$base
        astarithmetic $i -h1 $flatit2dir/flat-it2_"$k"_ccd"$h".fits -h1 2. gt nan where -o $maskedcornerdirbad/$base;
      done
      echo done > $maskedcornerdonebad
    fi
  fi
      
  # Now put all the images together
  maskedcornerdir=$BDIR/masked-corner
  maskedcornerdone=$maskedcornerdir/done_"$k"_ccd"$h".txt
  if ! [ -d $maskedcornerdir ]; then mkdir $maskedcornerdir; fi
  if [ -f $maskedcornerdone ]; then
    echo " "
    echo "Masked images moved"
    echo " "
  else
    mv $BDIR/masked-corner-secondflat-badcolumns/*_ccd"$h".fits $maskedcornerdir
    mv $BDIR/masked-corner-firstflat-badcolumns/*_ccd"$h".fits $maskedcornerdir
    mv $BDIR/masked-corner-fourflat-badcolumns/*_ccd"$h".fits $maskedcornerdir
    mv $BDIR/masked-corner-thirdflat-badcolumns/*_ccd"$h".fits $maskedcornerdir
    echo done > $maskedcornerdone
  fi
  
done



#########################################################################
#                       PREPARE FOR ASTROMETRY                          #
#########################################################################
# First of all, be sure to have Astrometry.net installed, scamp,swarp and
# sextractor.
# First you have to download the Gaia catalog to build the indexes that
# solve-field needs to compute the astrometric solutions for each frame.
# Second, you have to create these indexes. Third, run solve-field and
# compute the astrometric solutions.
#At this point the images have a first astrometric solution but there
# could be distortions of the camera plane so that the solutions computed
# are not the real position in the sky ( solve field assumes no distortions).
# To compute the distortions coefficients we first have to run Sextractor
# to create the catalog of the stars in that frame and then we run scamp
# to compute the distortions.
# Scamp will create files.head containing such coefficients.
# To run Sextractor, scamp and solve field you will need some configuration
# files



#########################################################################
#                              ASTROMETRY                               #
#########################################################################

for h in $(seq 1 4);do
#for h in 1;do
  # Download Gaia catalog to build the indexes needed by solve-field
  if [ $gaia_query = 1 ]; then
    query_param="gaia \
                 --dataset=edr3 \
                 --center=$ra_gal,$dec_gal \
                 --radius=0.9 \
                 --column=ra,dec,phot_g_mean_mag"
    catdir=$BDIR/catalogs
    catdone=$catdir/"$obj"_Gaia_eDR3.fits
    if ! [ -d $catdir ]; then mkdir $catdir; fi
    if [ -f $catdone ]; then
      echo " "
      echo "Catalog is already downloaded "
      echo " "
    else
      astquery $query_param -o $catdir/"$obj"_Gaia_eDR3.fits
    fi
  fi


  # Making the indexes
  indexdir=$BDIR/indexes
  indexdone=$indexdir/done_"$k".txt
  if ! [ -d $indexdir ]; then mkdir $indexdir; fi
  if [ -f $indexdone ]; then
    echo " "
    echo "Gaia eDR3 indexes are already created"
    echo " "
  else
    for re in $(seq -5 6); do
      build-astrometry-index -i $catdir/"$obj"_Gaia_eDR3.fits -e1 \
                             -P $re \
                             -S phot_g_mean_mag\
                             -E -A RA -D DEC\
                             -o $indexdir/index_$re.fits;
    done
    echo done > $indexdone
  fi


  
  # Solving the images
  astrocfg=$CDIR/astrometry_$obj.cfg
  if [ -f $astrocfg ]; then
    echo " "
    echo "Astrometry config file are already created "
    echo " "
  else
    echo inparallel > $astrocfg
    echo cpulimit 300 >> $astrocfg
    echo "add_path $indexdir" >> $astrocfg
    echo autoindex >> $astrocfg
  fi
  astroimadir=$BDIR/astro-ima
  astroimadone=$astroimadir/done_"$k"_ccd"$h".txt
  if ! [ -d $astroimadir ]; then mkdir $astroimadir; fi
  if [ -f $astroimadone ]; then
    echo " "
    echo "Images are already astrometrized for extension $h"
    echo " "
  else
    for a in $(seq 1 $n_exp); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$maskedcornerdir/$base
        solve-field $i --no-plots                    \
		-L 0.32 -H 0.34 -u arcsecperpix              \
		--overwrite --extension 1 --config $astrocfg \
        --use-sextractor --sextractor-path=/usr/local/bin/sextractor \
		-Unone --temp-axy -Snone -Mnone -Rnone -Bnone -N$astroimadir/$base --no-remove-lines --uniformize 0;
        #--use-sextractor --sextractor-path=/usr/local/bin/sextractor
        #--use-source-extractor                 \
        #--source-extractor-path=/Users/giuliagolini/Downloads/sextractor-2.25.0/src/sex
    done
    rm -f $flatit2imadir/*.wcs
    echo done > $astroimadone
  fi
  
  #########################################################################

  #################       Distorsion correction        ####################

  #########################################################################
  #
  # Making sex catalogs
  sexcfg=$CDIR/default.sex
  sexparam=$CDIR/default.param
  sexconv=$CDIR/default.conv
  sexdir=$BDIR/sex-it1
  sexdone=$sexdir/done_"$k"_ccd"$h".txt
  if ! [ -d $sexdir ]; then mkdir $sexdir; fi
  if [ -f $sexdone ]; then
    echo " "
    echo "Sex catalogs are already done for extension $h"
    echo " "
  else
    for a in $(seq 1 $n_exp); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$astroimadir/$base
        sex $i -c $sexcfg -PARAMETERS_NAME $sexparam -FILTER_NAME $sexconv -CATALOG_NAME $sexdir/"$obj"_Sloan-"$k"_"$a"_ccd"$h".cat;
        
    done
    echo done > $sexdone
  fi

  


  # Making scamp headers
  scampcfg=$CDIR/scamp.cfg
  scampdir=$BDIR/scamp-it1
  scampres=$scampdir/results_Sloan-"$k"_ccd"$h"
  scampdone=$scampdir/done_"$k"_ccd"$h".txt
  if ! [ -d $scampdir ]; then mkdir $scampdir; fi
  if ! [ -d $scampres ]; then mkdir $scampres; fi
  if [ -f $scampdone ]; then
    echo " "
    echo "Scamp headers are already done for extension $h"
    echo " "
  else
    scamp -c $scampcfg $(ls $sexdir/*Sloan-"$k"*ccd"$h".cat) -ASTREF_WEIGHT 10 -DISTORT_DEGREES 1
    #scamp -c $scampcfg  $sexdir/NGC1042_Sloan-r_3_ccd1.cat -ASTREF_WEIGHT 10 -STABILITY_TYPE INSTRUMENT
    #mv *.pdf $scampres/
    
    echo done > $scampdone
  fi
  
  


  # At this point, is possible to remove the temporary files
  if [ $remove_tmp_4 = 1 ]; then
    #rm -f $maskedcornerdir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $subskydir/*Sloan-"$k"*_ccd"$h".fits
    #rm -f $flatit2imadir/*Sloan-"$k"*_ccd"$h".fits
    echo " "
    echo "Temporary files #4 removed"
    echo " "
  fi
  
  
done






# Here I can choose to swarp with SWARP or ASTWARP
# swarp into bigger grid and correct of distorsion to improve photometry
###########################################################################

############      Resampling each ccd into common grid       ##############

###########################################################################

#resampling_with_SWARP
funct_resampling_with_SWARP
#resampling_with_ASTWARP
#funct_resampling_with_ASTWARP




#########################################################################

#################           Sky Subtraction           ####################

#########################################################################



for h in $(seq 1 4);do
#for h in 1;do
  # SAVE SKY VALUE
  # The sky substraction is done by using the --checksky option in
  # noisechisel. Then, the sigma clipping median of the sky is computed
  # and substracted of the science images. Apart of this, all the sky 
  # medians are stored in a separated file, all together, for a late
  # clasification of the quality of the images.
  if [ $nc_iter2_flat = 1 ]; then
    noiseskydir=$BDIR/noise-sky
    skyfile=$skydir/sky_Sloan-"$k"_ccd"$h".txt
    noiseskydone=$noiseskydir/done_"$k"_ccd"$h".txt
    if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
    if [ -f $noiseskydone ]; then
      echo " "
      echo "Science images are 'noisechiseled' for sky substraction for extension $h"
      echo " "
    else
      for a in $(seq 1 $n_exp); do
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        i=$swarpdir/$base
        sky=$(echo $base | sed 's/.fits/_sky.fits/')
        out=$(echo $base | sed 's/.fits/.txt/')
        astnoisechisel $i $noisechisel_param --checksky -o $noiseskydir/$base
        sky_row=$(wc -l $skyfile | awk '{print $1}')
        if [[ "$sky_row" == "$n_exp" ]]; then
          echo " "
          echo "Sky file for ccd"$h" and "$k" band already exist"
          echo " "
        else
          mean=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
          std=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
          echo $base $mean $std >> $skyfile
        fi
        # save both mean and std
        m=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
        s=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
        echo "$m $s" > $noiseskydir/$out
        
        #astfits $noiseskydir/$sky --copy=DETECTED -o$noiseskydir/$base
        rm -f $noiseskydir/$sky
      done
      echo done > $noiseskydone
    fi
  fi

  

  
  # Substract the sky to the images
  subskydir=$BDIR/sub-sky-it1
  subskydone=$subskydir/done_"$k"_ccd"$h".txt
  if ! [ -d $subskydir ]; then mkdir $subskydir; fi
  if [ -f $subskydone ]; then
    echo " "
    echo "Sky substraction is done for extension $h"
    echo " "
  else
    for a in $(seq 1 $n_exp); do
        # read file with sky values
        base="$obj"_Sloan-"$k"_"$a"_ccd"$h".txt
        i=$noiseskydir/$base
        me=$(awk 'NR=='1'{print $1}' $i)
        input=$swarpdir/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        output=$subskydir/"$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
        astarithmetic $input -h1 $me - -o$output;
        # remove input to save space (we cannot! we need it for the second iteration)
    done
    echo done > $subskydone
  fi
  
done



# ora sigma clip mean per fare una prima prova

#astarithmetic $(ls $subskydir/*.fits) $(ls $subskydir/*.fits | wc -l) 3 0.2 sigclip-median -g1 -o 3198_halpha.fits
#exit 1


# Now I have to put the ccd together
entiredir=$BDIR/pointings
if ! [ -d $entiredir ]; then mkdir $entiredir; fi

entiredone=$entiredir/done_pointings.txt
if [ -f $entiredone ]; then
  echo " "
  echo "Pointings images done (after sky subtraction)"
  echo " "
else
    for a in $(seq 1 $n_exp); do
        ccd1=$subskydir/"$obj"_Sloan-"$k"_"$a"_ccd1.fits
        ccd2=$subskydir/"$obj"_Sloan-"$k"_"$a"_ccd2.fits
        ccd3=$subskydir/"$obj"_Sloan-"$k"_"$a"_ccd3.fits
        ccd4=$subskydir/"$obj"_Sloan-"$k"_"$a"_ccd4.fits
        astarithmetic $ccd1 $ccd2 $ccd3 $ccd4 -g1 4 sum -g1 -o$entiredir/"$obj"_Sloan-"$k"_"$a".fits
    done
    echo done > $entiredone
fi





# Choose the best frame ( the one that has lower rms sky)
# in file rms_best_ccd' + h + '.txt'
# only ccd4
python3.6 find_best_frame.py "$k" 1 $n_exp 4 $INDIR $obj


index=$(awk 'NR=='1'{print $1}' $BDIR/rms_best_ccd4.txt)

echo $index


# MOS IN THE LBT PIPELINE IS NOW THE POINTING WITH the ccd with BEST RMS
mos=$entiredir/"$obj"_Sloan-"$k"_"$index".fits

###########################################################################

############              Photometry correction                ############

###########################################################################


###########################################################################

############     Create the catalog for reference image        ############

###########################################################################
# Build a catalog with sextractor with an aperture of 2" on SDSS image
# Extracting the sources of the Sloan image.

retcatdir=retrived-catalog
retcatdone=$retcatdir/done_"$n".txt
if ! [ -d $retcatdir ]; then mkdir $retcatdir; fi
if [ -f $retcatdone ]; then
    echo " "
    echo " gaia dr3 catalog retrived"
    echo " "
else
    # TXT TO MASK
    astquery gaia --dataset=edr3  \
	    --center=$ra_gal,$dec_gal -r 0.9 \
	    -csource_id,ra,dec,parallax -o $retcatdir/retrived_cat.cat
    echo done > $retcatdone
fi

sdssdir=$BDIR/sdss-catalogs
sdssdone=$sdssdir/done_"$k".txt
if ! [ -d $sdssdir ]; then mkdir $sdssdir; fi
if [ -f sdss.cat ]; then
  echo " "
  echo "Sources extracted for $k band in best image"
  echo " "
else
    asttable $retcatdir/retrived_cat.cat -hSOURCE_ID -cRA,DEC  \
	    | awk '{print NR, $1, $2, 5, '$r_sdss_pix', 0, 0, 1, NR, 1}' > apertures_sdss.txt
    
    astmkprof apertures_sdss.txt --background=$mos --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=apertures_sdss.fits

    astmkcatalog apertures_sdss.fits -h1 --zeropoint=22.5 \
            --valuesfile=$mos --valueshdu=1 \
            --ids --ra --dec --magnitude --sum \
            --output=sdss.cat
    echo done > $sdssdone
fi




###########################################################################

############       Create the catalog for your frame         ############

###########################################################################
## Extracting the sources of the lbt image.
lbtdir=$BDIR/lbt-catalogs
lbtdone=$lbtdir/done_"$k".txt
if ! [ -d $lbtdir ]; then mkdir $lbtdir; fi
if [ -f $lbtdone ]; then
  echo " "
  echo "Sources  extracted for $k band in lbt frame "
  echo " "
else
  asttable $retcatdir/retrived_cat.cat -hSOURCE_ID -cRA,DEC \
	  | awk '!/^#/{print NR, $1, $2, 5, '$r_lbt_pix', 0, 0, 1, NR, 1}' > apertures_lbt.txt
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    i=$entiredir/$base
    astmkprof apertures_lbt.txt --background=$i --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$lbtdir/"$obj"_Sloan-"$k"_"$a"_apertures.fits
        
    astmkcatalog $lbtdir/"$obj"_Sloan-"$k"_"$a"_apertures.fits -h1 --zeropoint=0 \
            --valuesfile=$i --valueshdu=1 \
            --ids --ra --dec --magnitude --sum \
            --output=$lbtdir/$base
  done
  echo done > $lbtdone
fi





###########################################################################

######       Match between Sloan correct and frame catalog        ########

###########################################################################
# Match with the Sloan catalog
matchdir2=$BDIR/match-sdss-frames
matchdir2done=$matchdir2/done_"$k".txt
if ! [ -d $matchdir2 ]; then mkdir $matchdir2; fi
if [ -f $matchdir2done ]; then
  echo " "
  echo "Match between SDSS catalog and LBT catalogs done"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    i=$lbtdir/$base
    out=$matchdir2/$base
    astmatch sdss.cat --hdu=1 $i --hdu=1 --ccol1=RA,DEC \
	    --ccol2=RA,DEC --aperture=1/3600 \
	    --outcols=aRA,aDEC,aMAGNITUDE,aSUM,bMAGNITUDE,bSUM -o$out

  done
  echo done > $matchdir2done
fi



###########################################################################

######           Save only those alpha we can trust                ########

###########################################################################
# so only if Sloan_MAG is between 18 and 22.

alphatruedir=$BDIR/alpha-stars-true
alphatruedone=$alphatruedir/done_"$k".txt
if ! [ -d $alphatruedir ]; then mkdir $alphatruedir; fi
if [ -f $alphatruedone ]; then
  echo " "
  echo "true catalogs filtered computed for extension $h"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    f=$matchdir2/$base
    alphatruet=$alphatruedir/"$obj"_Sloan-"$k"_"$a".txt
    asttable $f -h1 --range=MAGNITUDE,12,22 -o$alphatruet
    asttable $alphatruet -h1 -c1,2,3,'arith $4 $6 /' -o$alphatruedir/$base
    mean=$(asttable $alphatruedir/$base -c'ARITH_1'| aststatistics --sclipparams=3,0.1 --sigclip-median)
    std=$(asttable $alphatruedir/$base -c'ARITH_1'| aststatistics --sigclip-std)
    echo "$mean $std" > $alphatruedir/alpha_"$obj"_Sloan-"$k"_"$a"_ccd"$h".txt;
    rm -f $alphatruedir/$base
  done
  echo done > $alphatruedone
fi




###########################################################################

############                Multiply for alpha              ##############

###########################################################################
# Multiply each frame for the corresponding alpha value
muldir=$BDIR/phot-corr-dir
muldone=$muldir/done_"$k"_ccd"$h".txt
if ! [ -d $muldir ]; then mkdir $muldir; fi
if [ -f $muldone ]; then
    echo " "
    echo "Multiplication for alpha is done for extension $h"
    echo " "
else
    for a in $(seq 1 $n_exp); do
        base="$obj"_Sloan-"$k"_"$a".fits
        f=$entiredir/$base
        alpha_cat=$alphatruedir/alpha_"$obj"_Sloan-"$k"_"$a"_ccd4.txt
        alpha=$(awk 'NR=='1'{print $1}' $alpha_cat)
        # check if the values are read in floating points
        astarithmetic $f -h1 $alpha x -o $muldir/$base
        # remove sky sub to save space
        #rm -f $subskydir/$base;
    done
    echo done > $muldone
fi

rm -r $sdssdir



echo "here we have all the pointings at the same level of a reference image"
# Photometry part done
# ----------------------------------------------------------



# UNTIL here we have all the pointings is physical units
###########################################################################

############                ASSIGN WEIGHT                   ##############

###########################################################################
# Compute rms and of the photometrized frames

noiseskydir=$BDIR/noise-sky-after-photometry1
noiseskydone=$noiseskydir/done_"$k".txt
if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
if [ -f $noiseskydone ]; then
  echo " "
  echo "rms frame computed after photometry"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    i=$muldir/$base
    sky=$(echo $base | sed 's/.fits/_sky.fits/')
    out=$(echo $base | sed 's/.fits/.txt/')
    astnoisechisel $i --checksky -o $noiseskydir/$base $noisechisel_param
    #--tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2
    # save both mean and std
    m=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
    s=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
    echo "$m $s" > $noiseskydir/$out
    rm -f $noiseskydir/$sky
  done
  echo done > $noiseskydone
fi

###########################################################################

############        The best frame after photometry          ##############
########## it should be the same but maybe the others have less rms #######

###########################################################################
# as the std of the sky necessary to compute the weights
python3.6 find_rms_entire.py "$k" 1 $n_exp $h $INDIR $obj

#exit 1
###########################################################################

############                Compute weights                  ##############

###########################################################################

# Reading std values and save weight for each frame

wdir=$BDIR/weight-dir
wdone=$wdir/done_"$k".txt
if ! [ -d $wdir ]; then mkdir $wdir; fi

wonlydir=$BDIR/only-w-dir
wonlydone=$wonlydir/done_"$k".txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
if [ -f $wdone ]; then
  echo " "
  echo " weights computation done for extension $h"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    f=$muldir/$base
    rms_min=$(awk 'NR=='1'{print $1}' $BDIR/rms_min_val-1.txt)
    echo $rms_min
    rms_f=$(awk 'NR=='1'{print $2}' $noiseskydir/"$obj"_Sloan-"$k"_"$a".txt)
    echo $rms_f
    weight=$(astarithmetic $rms_min $rms_f / --quiet)
    echo "$weight" > $wdir/"$obj"_Sloan-"$k"_"$a".txt    #  saving into file
    # multiply each image for its weight
    wixi_im=$wdir/$base           # frame x weight
    w_im=$wonlydir/$base           # only weight
    astarithmetic $f -h1 $weight x -o$wixi_im
    astarithmetic $wixi_im -h1 $f -h1 / -o$w_im
    # now you can remove files in the photometry corrected dir
    rm -f $f
  done
  echo done > $wdone
  echo done > $wonlydone
fi



###########################################################################

############                Outliers treatment               ##############

###########################################################################
# Common part
# Remove outliers before the final coadd by using sigclip-median and sigclip
# -std
# This is particularly important to remove cosmic rays:
# more here https://www.bnl.gov/isd/documents/89281.pdf
# and https://www.astronomy.ohio-state.edu/MDM/OSMOS/CCD_CosmicRays_groom.pdf
# to build a mask of the ouliers pixels
# I reject above and belov 3 std

###########################################################################

############             save up and low limits              ##############

###########################################################################

clippingdir=$BDIR/clipping-outliers
clippingdone=$clippingdir/done_"$k".txt
if ! [ -d $clippingdir ]; then mkdir $clippingdir; fi
if [ -f $clippingdone ]; then
    echo " "
    echo "Upper and lower limits computed"
    echo " "
else
    # Compute clipped median and std
    med_im=$clippingdir/median_image.fits
    std_im=$clippingdir/std_image.fits
    astarithmetic $(ls $wdir/*.fits) $(ls $wdir/*.fits | wc -l) 3 0.2 sigclip-median -g1 -o$med_im
    astarithmetic $(ls $wdir/*.fits) $(ls $wdir/*.fits | wc -l) 3 0.2 sigclip-std -g1 -o$std_im
    # Compute "borders" images
    up_lim=$clippingdir/upperlim.fits
    lo_lim=$clippingdir/lowerlim.fits
    astarithmetic 1.5 $std_im x -o thresh.fits
    
    astarithmetic $med_im thresh.fits + -g1 float32 -o $up_lim
    astarithmetic $med_im thresh.fits - -g1 float32 -o $lo_lim
    
    rm -f $med_im $std_im
    echo done > $clippingdone
fi



###########################################################################

############             Masks and mask                      ##############

###########################################################################
# Now I need a loop for each ccd to compute the masks and mask the rejected
# pixels
# w*I
mowdir=$BDIR/weight-dir-no-outliers
mowdone=$mowdir/done_"$k".txt
if ! [ -d $mowdir ]; then mkdir $mowdir; fi

# only weight
moonwdir=$BDIR/only-weight-dir-no-outliers
moonwdone=$moonwdir/done_"$k".txt
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi


mowdone=$mowdir/done_"$k".txt
if [ -f $mowdone ]; then
    echo " "
    echo "numerator and denominater masked"
    echo " "
else
    for a in $(seq 1 $n_exp); do
        base="$obj"_Sloan-"$k"_"$a".fits
        tmp_ab=$mowdir/"$obj"_Sloan-"$k"_"$a"_maskabove.fits
        wom=$mowdir/$base
        astarithmetic $wdir/$base -h1 set-i i i $clippingdir/upperlim.fits -h1 gt nan where float32 -o $tmp_ab
        # here equal to the wi * Ii image masked
        astarithmetic $tmp_ab -h1 set-i i i $clippingdir/lowerlim.fits -h1 lt nan where float32 -o$wom
        # save the new mask
        mask=$mowdir/"$obj"_Sloan-"$k"_"$a"_mask.fits
        astarithmetic $wom -h1 isblank float32 -o $mask
        # mask the onlyweight image
        owom=$moonwdir/$base
        astarithmetic $wonlydir/$base $mask -g1 1 eq nan where float32  -o $owom
        
        # Remove temporary files
        rm -f $tmp_ab
        rm -f $mask
        
        # Remove images with outliers
        rm -f $wdir/$base
        rm -f $wonlydir/$base
    echo done > $mowdone
    done
fi





###########################################################################

############                Final coadd                      ##############

###########################################################################

# Doing the final coadd and weight image
coaddir=$BDIR/final-coadd-1
coaddone=$coaddir/done_"$k".txt
if ! [ -d $coaddir ]; then mkdir $coaddir; fi
if [ -f $coaddone ]; then
    echo " "
    echo "sum of weight*frame coadd is done for 1st iteration"
    echo " "
else
    astarithmetic $(ls $mowdir/*.fits) $(ls $mowdir/*.fits | wc -l) sum -g1 -o$coaddir/"$k"_wx.fits
    astarithmetic $(ls $moonwdir/*.fits ) $(ls $moonwdir/*.fits  | wc -l) sum -g1 -o$coaddir/"$k"_w.fits
    astarithmetic $coaddir/"$k"_wx.fits -h1 $coaddir/"$k"_w.fits -h1 / -o$coaddir/"$obj"_Sloan-"$k".fits
    echo done > $coaddone
fi

exit 1

# remove temporary files
#rm -f $mowdir/*.fits
#rm -f $moonwdir/*.fits


###########################################################################

#############              Create bigger mask                 #############

###########################################################################
# Now we have to run noisechisel to the final coadd trying to get only
# those pixels with sky. With this mask, more deep than the individial
# ones, we're going to mask the individual images after the sky
# substraction. And, then, redo the final coadd.
noisechisel_param_coadd="--tilesize=80,80 \
                       --minskyfrac=0.9 \
                       --meanmedqdiff=0.1 \
                       --snthresh=5.2 \
                       --detgrowquant=0.7 \
                       --detgrowmaxholesize=50000 \
                       --interpnumngb=40"
                       
noisechisel_param_coadd=" "
coaddmaskdir=$BDIR/coadd-mask
coaddmaskdone=$coaddmaskdir/done_"$k".txt
if ! [ -d $coaddmaskdir ]; then mkdir $coaddmaskdir; fi
if [ -f $coaddmaskdone ]; then
    echo " "
    echo "First iteration coadd is 'noisechiseled'"
    echo " "
else
    astnoisechisel $coaddir/"$obj"_Sloan-"$k".fits $noisechisel_param_coadd -o $coaddmaskdir/coadd_mask_full.fits
    #Save only detections to save space
    astfits $coaddmaskdir/coadd_mask_full.fits --copy=2 float32 -o$coaddmaskdir/coadd_mask.fits
    rm -f $coaddmaskdir/coadd_mask_full.fits
    echo done > $coaddmaskdone
fi







# Loop for 4 ccd begins again
for h in $(seq 1 4); do
    ###########################################################################

    ############     Apply the new mask and subtract sky         ##############

    ###########################################################################
    
    subskydir2=$BDIR/sub-sky-ima2               # sky sub ima
    if ! [ -d $subskydir2 ]; then mkdir $subskydir2; fi

    
    skydir=$BDIR/sky-sub2                       # sky values
    if ! [ -d $skydir ]; then mkdir $skydir; fi
    
    cropcoadmaskdir=$BDIR/crop-coadd-mask       # mask image
    if ! [ -d $cropcoadmaskdir ]; then mkdir $cropcoadmaskdir; fi
    
    skyimadone=$subskydir2/done_"$k"_ccd"$h".txt
    skyfile2=$skydir/sky_Sloan-"$k"_ccd"$h".txt
    if [ -f $skyimadone ]; then
        echo " "
        echo "sky subtracted it2 images resampled done for extension $h"
        echo " "
    else
        for a in $(seq 1 $n_exp); do
            base="$obj"_Sloan-"$k"_"$a"_ccd"$h".fits
            i=$swarpdir/$base
            j=$cropcoadmaskdir/$base

            # CREATE PREVIOUS NOISECHISEL MASK
            sky=$(echo $base | sed 's/.fits/_sky.fits/')
            out=$(echo $base | sed 's/.fits/.txt/')

            astnoisechisel $i $noisechisel_param -o $cropcoadmaskdir/temp2.fits
            # Add togheter the 2 masks
            astarithmetic $cropcoadmaskdir/temp2.fits $coaddmaskdir/coadd_mask.fits -g1 + float32 -o$cropcoadmaskdir/bigmask.fits
            # Mask the frame
            astarithmetic $i -h1 $cropcoadmaskdir/bigmask.fits -h1 0 gt nan where float32 -o$j
            # save both mean and std
            mean=$(aststatistics $j -h1 -sky --quiet | awk 'NR=='6'{print $2}')
            std=$(aststatistics $j -h1 -sky --quiet | awk 'NR=='7'{print $3}')
            echo "$mean $std" > $skydir/"$obj"_Sloan-"$k"_"$a"_ccd"$h".txt #  saving into file
            echo "$mean $std" >> $skyfile2
            # sky subtraction
            astarithmetic $i -h1 $mean - -o$subskydir2/$base
         done
         echo done > $skyimadone
        
    fi


done





# Now I have to put the ccd together
entiredir=$BDIR/pointings2
if ! [ -d $entiredir ]; then mkdir $entiredir; fi

entiredone=$entiredir/done_pointings.txt
if [ -f $entiredone ]; then
  echo " "
  echo "Pointings images done (after sky subtraction)"
  echo " "
else
    for a in $(seq 1 $n_exp); do
        ccd1=$subskydir2/"$obj"_Sloan-"$k"_"$a"_ccd1.fits
        ccd2=$subskydir2/"$obj"_Sloan-"$k"_"$a"_ccd2.fits
        ccd3=$subskydir2/"$obj"_Sloan-"$k"_"$a"_ccd3.fits
        ccd4=$subskydir2/"$obj"_Sloan-"$k"_"$a"_ccd4.fits
        astarithmetic $ccd1 $ccd2 $ccd3 $ccd4 -g1 4 sum -g1 -o$entiredir/"$obj"_Sloan-"$k"_"$a".fits
    done
    echo done > $entiredone
fi





# Choose the best frame ( the one that has lower rms sky)
# in file rms_best_ccd' + h + '.txt'
# only ccd4
python3.6 find_best_frame2.py "$k" 1 $n_exp 4 $INDIR $obj


index=$(awk 'NR=='1'{print $1}' $BDIR/rms_best_ccd4-2.txt)

echo $index


# MOS IN THE LBT PIPELINE IS NOW THE POINTING WITH the ccd with BEST RMS
mos=$entiredir/"$obj"_Sloan-"$k"_"$index".fits

###########################################################################

############              Photometry correction                ############

###########################################################################


###########################################################################

############     Create the catalog for reference image        ############

###########################################################################
# Build a catalog with sextractor with an aperture of 2" on SDSS image
# Extracting the sources of the Sloan image.

sdssdir=$BDIR/sdss-catalogs2
sdssdone=$sdssdir/done_"$k".txt
if ! [ -d $sdssdir ]; then mkdir $sdssdir; fi
if [ -f sdss.cat ]; then
  echo " "
  echo "Sources extracted for $k band in best image"
  echo " "
else
    asttable $retcatdir/retrived_cat.cat -hSOURCE_ID -cRA,DEC  \
        | awk '{print NR, $1, $2, 5, '$r_sdss_pix', 0, 0, 1, NR, 1}' > apertures_sdss.txt
    
    astmkprof apertures_sdss.txt --background=$mos --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=apertures_sdss.fits

    astmkcatalog apertures_sdss.fits -h1 --zeropoint=22.5 \
            --valuesfile=$mos --valueshdu=1 \
            --ids --ra --dec --magnitude --sum \
            --output=sdss.cat
    echo done > $sdssdone
fi




###########################################################################

############       Create the catalog for your frame         ############

###########################################################################
## Extracting the sources of the lbt image.
lbtdir=$BDIR/lbt-catalogs2
lbtdone=$lbtdir/done_"$k".txt
if ! [ -d $lbtdir ]; then mkdir $lbtdir; fi
if [ -f $lbtdone ]; then
  echo " "
  echo "Sources  extracted for $k band in lbt frame "
  echo " "
else
  asttable $retcatdir/retrived_cat.cat -hSOURCE_ID -cRA,DEC \
      | awk '!/^#/{print NR, $1, $2, 5, '$r_lbt_pix', 0, 0, 1, NR, 1}' > apertures_lbt.txt
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    i=$entiredir/$base
    astmkprof apertures_lbt.txt --background=$i --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$lbtdir/"$obj"_Sloan-"$k"_"$a"_apertures.fits
        
    astmkcatalog $lbtdir/"$obj"_Sloan-"$k"_"$a"_apertures.fits -h1 --zeropoint=0 \
            --valuesfile=$i --valueshdu=1 \
            --ids --ra --dec --magnitude --sum \
            --output=$lbtdir/$base
  done
  echo done > $lbtdone
fi





###########################################################################

######       Match between Sloan correct and frame catalog        ########

###########################################################################
# Match with the Sloan catalog
matchdir2=$BDIR/match-sdss-frames2
matchdir2done=$matchdir2/done_"$k".txt
if ! [ -d $matchdir2 ]; then mkdir $matchdir2; fi
if [ -f $matchdir2done ]; then
  echo " "
  echo "Match between SDSS catalog and LBT catalogs done"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    i=$lbtdir/$base
    out=$matchdir2/$base
    astmatch sdss.cat --hdu=1 $i --hdu=1 --ccol1=RA,DEC \
        --ccol2=RA,DEC --aperture=1/3600 \
        --outcols=aRA,aDEC,aMAGNITUDE,aSUM,bMAGNITUDE,bSUM -o$out

  done
  echo done > $matchdir2done
fi



###########################################################################

######           Save only those alpha we can trust                ########

###########################################################################
# so only if Sloan_MAG is between 12 and 22.

alphatruedir=$BDIR/alpha-stars-true2
alphatruedone=$alphatruedir/done_"$k".txt
if ! [ -d $alphatruedir ]; then mkdir $alphatruedir; fi
if [ -f $alphatruedone ]; then
  echo " "
  echo "true catalogs filtered computed for extension $h"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    f=$matchdir2/$base
    alphatruet=$alphatruedir/"$obj"_Sloan-"$k"_"$a".txt
    asttable $f -h1 --range=MAGNITUDE,12,22 -o$alphatruet
    asttable $alphatruet -h1 -c1,2,3,'arith $4 $6 /' -o$alphatruedir/$base
    mean=$(asttable $alphatruedir/$base -c'ARITH_1'| aststatistics --sclipparams=3,0.1 --sigclip-median)
    std=$(asttable $alphatruedir/$base -c'ARITH_1'| aststatistics --sigclip-std)
    echo "$mean $std" > $alphatruedir/alpha_"$obj"_Sloan-"$k"_"$a"_ccd"$h".txt;
    rm -f $alphatruedir/$base
  done
  echo done > $alphatruedone
fi




###########################################################################

############                Multiply for alpha              ##############

###########################################################################
# Multiply each frame for the corresponding alpha value
muldir=$BDIR/phot-corr-dir2
muldone=$muldir/done_"$k"_ccd"$h".txt
if ! [ -d $muldir ]; then mkdir $muldir; fi
if [ -f $muldone ]; then
    echo " "
    echo "Multiplication for alpha is done for extension $h"
    echo " "
else
    for a in $(seq 1 $n_exp); do
        base="$obj"_Sloan-"$k"_"$a".fits
        f=$entiredir/$base
        alpha_cat=$alphatruedir/alpha_"$obj"_Sloan-"$k"_"$a"_ccd4.txt
        alpha=$(awk 'NR=='1'{print $1}' $alpha_cat)
        # check if the values are read in floating points
        astarithmetic $f -h1 $alpha x -o $muldir/$base
        # remove sky sub to save space
        #rm -f $subskydir/$base;
    done
    echo done > $muldone
fi

rm -r $sdssdir




echo "here we have all the pointings is physical units of nanomaggies"
echo "again"


# UNTIL here we have all the pointings is physical units
###########################################################################

############                ASSIGN WEIGHT                   ##############

###########################################################################
# Compute rms and of the photometrized frames

noiseskydir=$BDIR/noise-sky-after-photometry2
noiseskydone=$noiseskydir/done_"$k".txt
if ! [ -d $noiseskydir ]; then mkdir $noiseskydir; fi
if [ -f $noiseskydone ]; then
  echo " "
  echo "rms frame computed after photometry"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    i=$muldir/$base
    sky=$(echo $base | sed 's/.fits/_sky.fits/')
    out=$(echo $base | sed 's/.fits/.txt/')
    astnoisechisel $i --checksky -o $noiseskydir/$base $noisechisel_param
    #--tilesize=20,20 --interpnumngb=5 --dthresh=0.1 --snminarea=2
    # save both mean and std
    m=$(aststatistics $noiseskydir/$sky -hSKY --sigclip-mean)
    s=$(aststatistics $noiseskydir/$sky -hSTD --sigclip-mean)
    echo "$m $s" > $noiseskydir/$out
    rm -f $noiseskydir/$sky
  done
  echo done > $noiseskydone
fi

###########################################################################

############        The best frame after photometry          ##############
########## it should be the same but maybe the others have less rms #######

###########################################################################
# as the std of the sky necessary to compute the weights
python3.6 find_rms_entire2.py "$k" 1 $n_exp $h $INDIR $obj


###########################################################################

############                Compute weights                  ##############

###########################################################################

# Reading std values and save weight for each frame

wdir=$BDIR/weight-dir2
wdone=$wdir/done_"$k".txt
if ! [ -d $wdir ]; then mkdir $wdir; fi

wonlydir=$BDIR/only-w-dir2
wonlydone=$wonlydir/done_"$k".txt
if ! [ -d $wonlydir ]; then mkdir $wonlydir; fi
if [ -f $wdone ]; then
  echo " "
  echo " weights computation done for extension $h"
  echo " "
else
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a".fits
    f=$muldir/$base
    rms_min=$(awk 'NR=='1'{print $1}' $BDIR/rms_min_val-2.txt)
    echo $rms_min
    rms_f=$(awk 'NR=='1'{print $2}' $noiseskydir/"$obj"_Sloan-"$k"_"$a".txt)
    echo $rms_f
    weight=$(astarithmetic $rms_min $rms_f / --quiet)
    echo "$weight" > $wdir/"$obj"_Sloan-"$k"_"$a".txt    #  saving into file
    # multiply each image for its weight
    wixi_im=$wdir/$base           # frame x weight
    w_im=$wonlydir/$base           # only weight
    astarithmetic $f -h1 $weight x -o$wixi_im
    astarithmetic $wixi_im -h1 $f -h1 / -o$w_im
    # now you can remove files in the photometry corrected dir
    rm -f $f
  done
  echo done > $wdone
  echo done > $wonlydone
fi




###########################################################################

############                Improve distortion               ##############

###########################################################################


    


#########################################################################

#################       Distortion correction        ####################

#########################################################################
#
# Making sex catalogs
sexcfg=$CDIR/default.sex
sexparam=$CDIR/default.param
sexconv=$CDIR/default.conv
sexdir=$BDIR/sex-it-final
sexdone=$sexdir/done_"$k".txt
if ! [ -d $sexdir ]; then mkdir $sexdir; fi
if [ -f $sexdone ]; then
echo " "
echo "Sex catalogs are already done for pointings"
echo " "
else
  for i in $muldir/"$obj"_Sloan-"$k"_"$a".fits; do
    base=$(basename $i)
    sex $i -c $sexcfg -PARAMETERS_NAME $sexparam -FILTER_NAME $sexconv -CATALOG_NAME $sexdir/"$base".cat;
  done
  echo done > $sexdone
fi




# Making scamp headers
scampcfg=$CDIR/scamp.cfg
scampdir=$BDIR/scamp-it-final
scampres=$scampdir/results_Sloan-"$k"_ccd"$h"
scampdone=$scampdir/done_"$k".txt
if ! [ -d $scampdir ]; then mkdir $scampdir; fi
if ! [ -d $scampres ]; then mkdir $scampres; fi
if [ -f $scampdone ]; then
echo " "
echo "Scamp headers are already done for pointings"
echo " "
else
scamp -c $scampcfg $(ls $sexdir/"$obj"_Sloan-"$k"_*.cat)
#mv *.pdf $scampres/
echo done > $scampdone
fi

 
 
# NOW i have the distorsion headers for the pointings
# swarp into bigger grid and correct of distorsion
###########################################################################

############      Resampling each POINTING into common grid  ##############

###########################################################################

# Running swarp over both images and bleeding masks
swarpcfg=$CDIR/swarp.cfg
swarpdir=$BDIR/improved_wi

wswarpdir=$BDIR/improved_w
if ! [ -d $wswarpdir ]; then mkdir $wswarpdir; fi


swarpdone=$swarpdir/done_"$k".txt
if ! [ -d $swarpdir ]; then mkdir $swarpdir; fi
if [ -f $swarpdone ]; then
  echo " "
  echo "All the pointings are already distortion corrected"
  echo " "
else
  cp $sexdir/*Sloan-"$k"_*.head $wdir
  cp $sexdir/*Sloan-"$k"_*.head $wonlydir
  for a in $(seq 1 $n_exp); do
    base="$obj"_Sloan-"$k"_"$a"
    i=$wdir/"$base".fits
    SWarp -c $swarpcfg $i -CENTER $ra_gal,$dec_gal -IMAGE_SIZE $mos_dim_x,$mos_dim_y -IMAGEOUT_NAME $swarpdir/"$base"_m.fits -WEIGHTOUT_NAME $swarpdir/"$base"_m_w.fits -SUBTRACT_BACK N -PIXEL_SCALE 0.333 -PIXELSCALE_TYPE MANUAL
    
    # Mask bad pixels
    astarithmetic $swarpdir/"$base"_m_w.fits -h0 $swarpdir/"$base"_m_w.fits -h0 0 lt nan where -otemp"$base".fits
    # set zero to nan
    astarithmetic $swarpdir/"$base"_m.fits -h0 temp"$base".fits -h1 0 eq nan where -o$swarpdir/"$base".fits
    
    rm -f $swarpdir/"$base"_m.fits
    rm -f $swarpdir/"$base"_m_w.fits
    rm -f temp"$base".fits
    
    
    i=$wonlydir/"$base".fits
    SWarp -c $swarpcfg $i -CENTER $ra_gal,$dec_gal -IMAGE_SIZE $mos_dim_x,$mos_dim_y -IMAGEOUT_NAME $wswarpdir/"$base"_m.fits -WEIGHTOUT_NAME $wswarpdir/"$base"_m_w.fits -SUBTRACT_BACK N -PIXEL_SCALE 0.333 -PIXELSCALE_TYPE MANUAL
    
    # Mask bad pixels
    astarithmetic $wswarpdir/"$base"_m_w.fits -h0 $wswarpdir/"$base"_m_w.fits -h0 0 lt nan where -otemp"$base".fits
    # set zero to nan
    astarithmetic $wswarpdir/"$base"_m.fits -h0 temp"$base".fits -h1 0 eq nan where -o$wswarpdir/"$base".fits
    
    rm -f $wswarpdir/"$base"_m.fits
    rm -f $wswarpdir/"$base"_m_w.fits
    rm -f temp"$base".fits
    
 done
 echo done > $swarpdone
fi




# Now I have to remove outliers before stacking

CDIR=$SDIR/config
# Now again remove outliers
# Common part
# Remove outliers before the final coadd by using sigclip-median and sigclip-std
# to build a mask of the ouliers pixels
# I reject above and belov 3 std
# Doing the sigmaclip median image
clippingdir=$BDIR/clipping-outliers-f
clippingdone=$clippingdir/done_"$k".txt
if ! [ -d $clippingdir ]; then mkdir $clippingdir; fi
if [ -f $clippingdone ]; then
    echo " "
    echo "Images to remove outliers pointings computed"
    echo " "
else
    # Compute clipped median and std
    med_im=$clippingdir/median_image.fits
    std_im=$clippingdir/std_image.fits
    astarithmetic $(ls $swarpdir/*.fits) $(ls $swarpdir/*.fits | wc -l) 3 0.2 sigclip-median -g1 -o$med_im
    astarithmetic $(ls $swarpdir/*.fits) $(ls $swarpdir/*.fits | wc -l) 3 0.2 sigclip-std -g1 -o$std_im
    # Compute "borders" images
    up_lim=$clippingdir/upperlim.fits
    lo_lim=$clippingdir/lowerlim.fits
    thresh_up=$clippingdir/threesigmalim_up.fits
    thresh_d=$clippingdir/threesigmalim_d.fits
    astarithmetic 3 $std_im x -o$thresh_up
    astarithmetic 3 $std_im x -o$thresh_d
    astarithmetic $med_im $thresh_up + -g1 -o $up_lim
    astarithmetic $med_im $thresh_d - -g1 -o $lo_lim
    rm -f $thresh
    echo done > $clippingdone
fi



# Now I need a loop for each ccd to compute the masks and mask the rejected pixels
# w*I
mowdir=$BDIR/weight-dir-no-outliers-f
mowdone=$mowdir/done_"$k"_ccd"$h".txt
if ! [ -d $mowdir ]; then mkdir $mowdir; fi

# only weight
moonwdir=$BDIR/only-weight-dir-no-outliers-f
moonwdone=$moonwdir/done_"$k"_ccd"$h".txt
if ! [ -d $moonwdir ]; then mkdir $moonwdir; fi



mowdone=$mowdir/done_"$k".txt
if [ -f $mowdone ]; then
    echo " "
    echo "w_i*I_i with new mask"
    echo " "
else
      for i in $swarpdir/"$obj"_Sloan-"$k"_*.fits; do
            base=$(basename $i)
            tmp_ab=$mowdir/"$obj"_Sloan-"$k"_"$a"_maskabove.fits
            wom=$mowdir/$base
            astarithmetic $swarpdir/$base -h1 set-i i i $clippingdir/upperlim.fits -h1 gt nan where float32 -o $tmp_ab
            # here equal to the wi * Ii image masked
            astarithmetic $tmp_ab -h1 set-i i i $clippingdir/lowerlim.fits -h1 lt nan where float32 -o$wom
    
            # save the new mask
            mask=$mowdir/"$obj"_Sloan-"$k"_"$a"_mask.fits
            astarithmetic $wom -h1 isblank float32 -o $mask
            # mask the onlyweight image
            owom=$moonwdir/$base
            astarithmetic $wswarpdir/$base $mask -g1 1 eq nan where float32 -o $owom
        
            # Remove temporary files
            rm -f $tmp_ab
            rm -f $mask
            
            # Remove images with outliers
            rm -f $wdir/$base
            rm -f $wonlydir/$base
            echo done > $mowdone
        done
fi






###########################################################################

############                Final coadd 2                     ##############

###########################################################################

# Doing the final coadd and weight image
coaddir=$BDIR/final-coadd-2
coaddone=$coaddir/done_"$k".txt
if ! [ -d $coaddir ]; then mkdir $coaddir; fi
if [ -f $coaddone ]; then
    echo " "
    echo "sum of weight*frame coadd is done for 2nd iteration"
    echo " "
else
    astarithmetic $(ls $mowdir/*.fits) $(ls $mowdir/*.fits | wc -l) sum -g1 -o$coaddir/"$k"_wx.fits
    astarithmetic $(ls $moonwdir/*.fits ) $(ls $moonwdir/*.fits  | wc -l) sum -g1 -o$coaddir/"$k"_w.fits
    astarithmetic $coaddir/"$k"_wx.fits -h1 $coaddir/"$k"_w.fits -h1 / -o$coaddir/"$obj"_Sloan-"$k".fits
    echo done > $coaddone
fi

