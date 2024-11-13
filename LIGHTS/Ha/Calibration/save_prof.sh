#/bin/bash 

# This script is used to obtain the only Ha image.
# To do so I need a color image ( -2.5 log10 (fg/fr))
# and mask those pixels where the s/n is below 5 sigma.
# Copyright (C) 2021 Golini Giulia <giulia.golini@gmail.com>
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





# Stop the execution if something crashes
set -e





# Save the current system language, and then change it to English to avoid
# problems in some parts of the code (AWK with `,' instead of `.' for
# decimal separator).
export LANG=C





gal=$1
ra=$2
dec=$3
# input directory
indir=inputs-"$gal"




g=$indir/"$gal"-g_nanomaggy.fits        # g data
r=$indir/"$gal"-r_nanomaggy.fits        # r data
ha=build-"$gal"/Ha_calibrated-"$gal".fits     # Ha data
m=$indir/mask_"$gal".fits               # mask
#astarithmetic $m -h1 isblank float32 -omask_ok.fits



# Tmp Directory
tmpdir=tmpdir-"$gal"
if ! [ -d $tmpdir ]; then mkdir $tmpdir; fi






###########################################################################

#                           Create g-r image                              #

###########################################################################


gc=$indir/g_toha.fits
rc=$indir/r_toha.fits


mc=$tmpdir/maskcorrected.fits
sgc=$tmpdir/sg_toha.fits
src=$tmpdir/sr_toha.fits
suc=$tmpdir/su_toha.fits



# masked
gw=$tmpdir/g_m.fits
rw=$tmpdir/r_m.fits
haw=$tmpdir/ha_m.fits
haw_sub=$tmpdir/ha_sub_m.fits




# Rescale the mask
astwarp $m -h1 --gridfile=$ha -o $mc

# mask all of them

#astarithmetic $rc $mc -g1 0 ne nan where -o $rw
#astarithmetic $gc $mc -g1 0 ne nan where -o $gw
#astarithmetic $ha $mc -g1 0 ne nan where -o $haw
astarithmetic $rc $mc -g1 2000 eq nan where -o $rw
astarithmetic $gc $mc -g1 2000 eq nan where -o $gw
astarithmetic $ha $mc -g1 2000 eq nan where -o $haw





# Parmeters that have to be read from the table for each galaxy
ar=$4
pa=$5
Ag=$6
Ar=$7
mdim=$8         #apparent size of the galaxy
pixscale=0.333
telescope=lights




# -------------------------------------------------------------------------
#                           DEFINE FUNCTIONS
# -------------------------------------------------------------------------
# FUNCTION TO COMPUTE SKY
funct_computesky_ell() {
    # variables
    start_sky_as=$1
    end_sky_as=$2
    pix_scale=$3
    f=$4
    undersample=$5
    telescope=$6
    image=$7
    hdu=$8
    
    # convert sky as into pix
    start_sky_pix_g=$(astarithmetic $start_sky_as $pix_scale / --quiet |awk '{print $1}' )
    end_sky_pix_g=$(astarithmetic $end_sky_as $pix_scale / --quiet |awk '{print $1}' )
    
    # save sky value
    astscript-radial-profile $image -h$hdu --positionangle=$pa --axisratio=$ar --mode=wcs --center=$ra,$dec --rmax=$stop  --measure=sigclip-median --measure=sigclip-std --measure=sigclip-number --undersample=$undersample -o$bdir/to_prof_$f.cat 


    asttable $bdir/to_prof_"$f".cat  --range=1,$start_sky_pix_g,$end_sky_pix_g -o $bdir/skyregion.cat
    skyvalue=$(aststatistics $bdir/skyregion.cat --column=2 --median --quiet)
    skystd=$(aststatistics $bdir/skyregion.cat --column=3 --median --quiet)

    echo $skyvalue $skystd > $bdir/sky_"$f"_"$telescope".txt
    rm $bdir/to_prof_"$f".cat
}


# FUNCTION TO COMPUTE SB PROFILE
funct_computesbprofile() {
    # variables
    skystd=$1
    f=$2
    undersample=$3
    telescope=$4
    image=$5
    outputcat=$6
    hdu=$7
        
    
    astscript-radial-profile $image -h$hdu --positionangle=$pa --axisratio=$ar --mode=wcs --center=$ra,$dec --rmax=$stop  --measure=sigclip-median --measure=sigclip-std --measure=sigclip-number --undersample=$undersample -oprof_"$f".cat 

    #add errorbar!!
    asttable prof_"$f".cat -c1,2,3,4 -c'arith  SIGCLIP-NUMBER sqrt '  -ot0.txt
    asttable t0.txt -c1,2,3,4 -c'arith 1 $5 /' -ot.txt
    asttable t.txt -c1,2,3,4 -c'arith $5 '$skystd' x' -c'arith $5 SIGCLIP-STD x' -ott.txt
    asttable tt.txt -c1,2,3,4,5,6 -c'arith $5 $5 x $6 $6 x +' -o $outputcat
    rm t.txt tt.txt t0.txt prof_"$f".cat
    
}





# FUNCTION TO COMPUTE SB PROFILE WEDGES
funct_wedges() {
    # variables
    skystd=$1
    f=$2
    undersample=$3
    telescope=$4
    image=$5
    outputcat=$6
    hdu=$7
    alpha_i=$8
    alpha_f=$9
    
    astscript-radial-profile $image -h$hdu --positionangle=$pa --axisratio=$ar --mode=wcs --center=$ra,$dec --rmax=$stop  --measure=sigclip-median --measure=sigclip-std --measure=sigclip-number --undersample=$undersample -oprof_"$f".cat --azimuth=$alpha_i,$alpha_f

    #add errorbar!!
    asttable prof_"$f".cat -c1,2,3,4 -c'arith SIGCLIP-NUMBER sqrt '  -ot0.txt
    asttable t0.txt -c1,2,3,4 -c'arith 1 $5 /' -ot.txt
    asttable t.txt -c1,2,3,4 -c'arith $5 '$skystd' x' -c'arith $5 SIGCLIP-STD x' -ott.txt
    asttable tt.txt -c1,2,3,4,5,6 -c'arith $5 $5 x $6 $6 x +' -o $outputcat
    rm t.txt tt.txt t0.txt prof_"$f".cat
    
}

# FUNCTION TO COMPUTE SKY WEDGES
funct_computesky_we() {
    # variables
    start_sky_as=$1
    end_sky_as=$2
    pix_scale=$3
    f=$4
    undersample=$5
    telescope=$6
    image=$7
    hdu=$8
    alpha_i=$9
    alpha_f=${10}
    
    # convert sky as into pix
    start_sky_pix_g=$(astarithmetic $start_sky_as $pix_scale / --quiet |awk '{print $1}' )
    end_sky_pix_g=$(astarithmetic $end_sky_as $pix_scale / --quiet |awk '{print $1}' )
    
    # save sky value
    astscript-radial-profile $image -h$hdu --positionangle=$pa                      \
                    --axisratio=$ar --mode=wcs --center=$ra,$dec                    \
                    --rmax=$stop  --measure=sigclip-median --measure=sigclip-std    \
                    --measure=sigclip-number --azimuth=$alpha_i,$alpha_f    --keeptmp        \
                    --undersample=$undersample -o$bdir/to_prof_$f.cat


    asttable $bdir/to_prof_$f.cat  --range=1,$start_sky_pix_g,$end_sky_pix_g -o $bdir/skyregion.cat
    skyvalue=$(aststatistics $bdir/skyregion.cat --column=2 --median --quiet)
    skystd=$(aststatistics $bdir/skyregion.cat --column=3 --median --quiet)

    echo $skyvalue $skystd > $bdir/sky_"$f"_"$telescope".txt
    rm $bdir/to_prof_"$f".cat
}





####################################################################
##################                  LBT
# Create build directory for that galaxy
bmdir=build-"$gal"
if ! [ -d $bmdir ]; then mkdir $bmdir; fi

# Create build directory LBT
bdir=$bmdir/LBT
if ! [ -d $bdir ]; then mkdir $bdir; fi



# where should I stop to compute the profile?
stop=$(astarithmetic $mdim 7 x --quiet)
# where the sky starts and ends? in arcseconds
startsky=$(astarithmetic $mdim 1.1 x --quiet)            # 330"
endsky=$(astarithmetic $mdim 1.2 x --quiet)             # 360"
#step in pix?
step=10
pixscale=0.333
width_wedge=5




for f in "r" "g" "ha"; do
    
    # ###################################################################
    # Ellipses
    im=$tmpdir/"$f"_m.fits
    # compute sky,subtract and extract profile
    funct_computesky_ell $startsky $endsky $pixscale $f $step $telescope $im 1

    # Subtract the sky
    skysub=$bdir/"$telescope"_"$gal"_Sloan-"$f"_masked_skysub.fits
    astarithmetic $im -h1 $skyvalue - -o $skysub

    # Compute profile
    prof=$bdir/prof_"$telescope"_"$gal"_Sloan-"$f".txt
    funct_computesbprofile $skystd $f $step $telescope $skysub $prof 1
    
    # ###################################################################

    # Wedges
    # Compute sky
    #alpha_i=0
    #alpha_f=10
    #startsky=$(astarithmetic $mdim 1.1 x --quiet)            # 330"
    #endsky=$(astarithmetic $mdim 1.2 x --quiet)             # 360"

    #funct_computesky_we $startsky $endsky $pixscale $f $step $telescope $im 1 $alpha_i $alpha_f
    
    # Subtract the sky
    #skysub=$bdir/"$telescope"_"$gal"_Sloan-"$f"_masked_skysub_we.fits
    #astarithmetic $im -h1 $skyvalue - -o $skysub
    
    # Compute profile
    #prof_w=$bdir/prof_"$telescope"_"$gal"_Sloan-"$f"_upmayor.txt
    #funct_wedges $skystd $f $step $telescope $skysub $prof_w 1 $alpha_i $alpha_f

    #prof_w=$bdir/prof_"$telescope"_"$gal"_Sloan-"$f"_downmayor.txt
    #alpha_i=170
    #alpha_f=180
    #funct_wedges $skystd $f $step $telescope $skysub $prof_w 1 $alpha_i $alpha_f

    # ###################################################################
done



#python3.11 im_up.py "$gal"
#python3.11 im_down.py "$gal"
python3.11 image_ell.py "$gal"
