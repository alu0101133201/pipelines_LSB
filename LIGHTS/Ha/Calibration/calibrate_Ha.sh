#bin/bash


# first go to Skyviever, then you find the correct plate
# around the galaxy ( sdss photometry )
# I need to save a list with ra and dec of all the stars
# i the field for with I can retrieve the spectra
# https://dr16.sdss.org/optical/spectrum/search
# at this website search the corresponding plate
# select only ra and dec and copy and paste into a file.txt
# I want to download all the spectra of the stars in the SDSS
# and find the relation between F_Ha and F r_Sloan
# I download stars with 0<z<0.1 --> 230 000 stars
# try first with the stars in theplate of the glaxy (455 + 456)
# I used the spectra of stars in various field of sdss
# to understand the relation between the flux in h alpha and the flux in
# r Sloan of LBT.
# For each spectra I compute the color g - r
# using the convolution of the filters.
# to select the correct range of color to find the relation factor
# I study the distribution of the colors of my stars in the field
# once I remove the red ones ( degeneracy in the Fh/Fr factor)
# and select the range of values where I have the majority of stars
# I save the factor relative to that range of colors.
# Once I have the factor I have to calibrate the INT images knowing
# which is the expected value and the real Flux values of the stars.
# for each star I compute the factor and then I take the sigma clip mean od these
# values as a conversion factor.
# The fluxes in the data (mages) are computed using aperture photometry
# with diameter of 3 as = dimension of the diameter of the fiber to get the spectra.
#
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




# fits files are stored in the folder inputs    inputs
# you need a folder with the filters            filters
# a list of of stars to retrieve                list.txt
# This code need the files plot_lbtcolors.py and save_factor.py





# Stop the execution if something crashes
set -e





# Save the current system language, and then change it to English to avoid
# problems in some parts of the code (AWK with `,' instead of `.' for
# decimal separator).
export LANG=C





###########################################################################

#                           DOwnload spectra                            #

###########################################################################
# for each ra and dec in the list I run the python code to download the table
# directory to save spectras
spectra=spectra
if ! [ -d $spectra ]; then mkdir $spectra; fi
# code to download spectra in the list
#python3.11 ret_spectra.py
#exit 1


###########################################################################

#                            Save    SED                                 #

###########################################################################
# directory to save SEDs
seds=seds
if ! [ -d $seds ]; then mkdir $seds; fi
# save table with only the flux and the wavelength
sedsaved=$seds/done_sed.txt
if [ -f $sedsaved ]; then
    echo " "
    echo "SED already saved"
    echo " "
else
    for i in $(ls $spectra/*.fits); do
        base=$(basename $i)
        asttable $i -h1 -c1,2 -o $seds/$base.txt
    done
    echo done > $sedsaved
fi




###########################################################################
#                       Working on specific field
###########################################################################
# It can be that each field has its own "better range of stars
# LBT IMAGE choose a range of color and compute the factor to calibrate
# the range of color with higher number of star possibly
# target RA and DEC
#NGC6015 237.856 62.3109 0.36 120  0.041 0.029 300
gal=$1
ra=$2
dec=$3
ar=$4
pa=$5
Ag=$6
Ar=$7
mdim=$8         #apparent size of the galaxy
grmax=$9        # maximum color for stars
mulsigma=${10}  # divide the sigma around 2nd peak by this number

# input directory
indir=inputs-"$gal"



# Converto to nanomaggy ( if not already calibrated)
if [ -f $indir/"$gal"-g_nanomaggy.fits ]; then
    echo " "
    echo "Inputs in nanomaggy"
    echo " "
else
    for f in "r" "g"; do
        g=$indir/"$gal"_Sloan-"$f".fits
        g_n=$indir/"$gal"-"$f"_nanomaggy.fits
        zp=$(astfits $g -h1 | grep ZP | awk '{print $3}')
        astarithmetic $g -h1 $zp counts-to-nanomaggy -o $g_n
    done
fi



g=$indir/"$gal"-g_nanomaggy.fits
r=$indir/"$gal"-r_nanomaggy.fits
ha=$indir/"$gal"_Sloan-halpha.fits

# Prepare data
# resample lbt at int pixel scale (larger area so gain s/n)
lbtp=0.224
intp=0.333
###########################################################################

#                      Prepare data same dimension                        #

###########################################################################
#  create color image
gc=$indir/g_toha.fits
rc=$indir/r_toha.fits
    
g_r=$indir/g-r.fits
if [ -f $g_r ]; then
    echo " "
    echo "Inputs ready"
    echo " "
else

    # Rescale LBT images to Ha image ( So then they are all with the same size)
    astwarp $r -h1 --gridfile=$ha -o $rc
    astwarp $g -h1 --gridfile=$ha -o $gc
    
    # color image
    astarithmetic $gc -h1 set-i $rc -h1 set-j i j / log10 -2.5 x -o $g_r
fi



# Output Directories
# General directory for that galaxy
bdir=build-"$gal"
if ! [ -d $bdir ]; then mkdir $bdir; fi

plots=$bdir/plots
if ! [ -d $plots ]; then mkdir $plots; fi


outcat=$bdir/output-catalogs
if ! [ -d $outcat ]; then mkdir $outcat; fi





###########################################################################

########### GAIA OK FOR A FIRST TEST OF THE STARS THAT ARE OK ############

###########################################################################
# Remove sources on top of the galaxy
# GAIA CATALOG TO SELECT STARS WITH PARALLAX
# and All stars OUTSIDE THE GALAXY
ra_max=$(astarithmetic $ra 0.01 - --quiet)
ra_min=$(astarithmetic $ra 0.01 + --quiet)
de_max=$(astarithmetic $dec 0.01 + --quiet)
de_min=$(astarithmetic $dec 0.01 - --quiet)
retcatdir=$bdir/retrived-catalog
retcatdone=$retcatdir/done_"$n".txt
if ! [ -d $retcatdir ]; then mkdir $retcatdir; fi
if [ -f $retcatdone ]; then
    echo " "
    echo "Gaia dr3 catalog retrived"
    echo " "
else
    # TXT TO MASK
    astquery gaia --dataset=edr3  --center=$ra,$dec -r 0.3  -csource_id,ra,dec,parallax,phot_g_mean_mag  -o $retcatdir/r.cat
    asttable $retcatdir/r.cat -c1,2,3,4,5 --outpolygon=ra,dec --polygon="$ra_max,$de_min:$ra_min,$de_min:$ra_min,$de_max:$ra_max,$de_max" --range=parallax,0:4 --range=phot_g_mean_mag,19,25 -o $retcatdir/retrived.cat
    # NB:halpha_c.fits is the image in Ha
    astscript-psf-select-stars $ha -h1 --catalog=$retcatdir/retrived.cat --magnituderange=19,25  --mindistdeg=0.0001 -o $retcatdir/retrived_cat.cat
    rm $retcatdir/r.cat $retcatdir/retrived.cat
    echo done > $retcatdone
fi





###########################################################################

#                      Compute colors From Data                          #

###########################################################################
# Now that list should be read by gnuastro and I have to build
# apertures and compute the flux
# Dimension of the aperture?
# This should be the same as the dimension of the Fiber Aperture
# https://www.sdss4.org/dr12/spectro/spectro_basics/
# Select the dimension of the aperture for the aperture photometry(RADIUS)
d_aperture=2.                 # arcsec
intp=0.333                    # pix scale of the image
# Covert into pixels
r_pix_=$(astarithmetic $d_aperture $intp /  --quiet | awk '{printf("%.16f\n", $1)}'  )
r_pix=${r_pix_/\.*}
echo " "
echo "Number of pixels for the aperture photometry"
echo $r_pix         # to convert in integer


# Color

intdir=$bdir/INT-apertures
intdone=$intdir/done_"$k".txt
if ! [ -d $intdir ]; then mkdir $intdir; fi
if [ -f $outcat/hc.cat ]; then
  echo " "
  echo "Sources extracted for g-r image"
  echo " "
else
    asttable $retcatdir/retrived_cat.cat -h1 -cRA,DEC  | awk '{print NR, $1, $2, 5, '$r_pix', 0, 0, 1, NR, 1}' > $intdir/apertures_INT.txt
    
    astmkprof $intdir/apertures_INT.txt --background=$g_r --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$intdir/apertures_hc.fits

    astmkcatalog $intdir/apertures_hc.fits -h1 --zeropoint=0 \
            --valuesfile=$g_r --valueshdu=1 \
            --ids --ra --dec --mean \
            --output=$outcat/hc_.cat
    asttable $outcat/hc_.cat --range=mean,0.,1.5 -o$outcat/hc.cat     # basically all I have
    #rm hc_.cat
    echo done > $intdone
fi



# plot LBT stars with color
if [ -f $outcat/range.txt ]; then
  echo " "
  echo "Sources extracted for g-r image"
  echo " "
else
    python3.11 plot_lbtcolors.py "$gal" $grmax $mulsigma
fi



###########################################################################

#                      Compute magnitudes From Spectra                    #

###########################################################################
###########################################################################

#                            Compute Factor                               #

###########################################################################
# Now for every star
# I have to save into a table the f halpha, f r and color
# then plot, depending on teh color the point
# divide by groups of color
if [ -f $bdir/table_mag_reference_"$gal".txt ]; then
    echo " "
    echo "Magnitudes already saved"
    echo " "
else
    # Plot the SEDs
    python3.11 save_factor.py "$gal"
fi



# read correct range
ref=$outcat/range.txt
min=$(awk 'NR=='1'{print $1}' $ref)
max=$(awk 'NR=='1'{print $2}' $ref)

echo $min $max
asttable $outcat/hc.cat --range=mean,$min:$max -o$outcat/correct_range.cat



###########################################################################

#                      Compute FLUX  From HALPHA                   #

###########################################################################

intdir=$bdir/INT-apertures
intdone=$intdir/done_"$k".txt
if ! [ -d $intdir ]; then mkdir $intdir; fi
if [ -f $outcat/halpha.cat ]; then
  echo " "
  echo "Sources extracted for H ALPHA image"
  echo " "
else
    asttable $outcat/correct_range.cat -h1 -cRA,DEC  | \
        awk '{print NR, $1, $2, 5, '$r_pix', 0, 0, 1, NR, 1}' > $intdir/apertures_INT.txt
    
    astmkprof $intdir/apertures_INT.txt --background=$ha --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$intdir/apertures_hc.fits

    astmkcatalog $intdir/apertures_hc.fits -h1 --zeropoint=0 \
            --valuesfile=$ha --valueshdu=1 \
            --ids --ra --dec --sum \
            --output=$outcat/halpha.cat
    echo done > $intdone
fi



###########################################################################

#                      Compute FLUX  From LBT R                           #

###########################################################################

intdir=$bdir/LBT-apertures
intdone=$intdir/done_"$k".txt
if ! [ -d $intdir ]; then mkdir $intdir; fi
if [ -f $outcat/rc.cat ]; then
  echo " "
  echo "Sources extracted for LBT r Sloan image"
  echo " "
else
    asttable $outcat/correct_range.cat -h1 -cRA,DEC  \
        | awk '{print NR, $1, $2, 5, '$r_pix', 0, 0, 1, NR, 1}' > $intdir/apertures_INT.txt
    astmkprof $intdir/apertures_INT.txt --background=$rc --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$intdir/apertures_rc.fits

    astmkcatalog $intdir/apertures_rc.fits -h1 --zeropoint=0 \
            --valuesfile=$rc --valueshdu=1 \
            --ids --ra --dec --sum \
            --output=$outcat/rc.cat
    echo done > $intdone
fi


###########################################################################

#                      Compute FLUX  From LBT G                           #

###########################################################################

intdir=$bdir/LBT-apertures-g
intdone=$intdir/done_"$k".txt
if ! [ -d $intdir ]; then mkdir $intdir; fi
if [ -f $outcat/gc.cat ]; then
  echo " "
  echo "Sources extracted for LBT g Sloan image"
  echo " "
else
    asttable $outcat/correct_range.cat -h1 -cRA,DEC \
        | awk '{print NR, $1, $2, 5, '$r_pix', 0, 0, 1, NR, 1}' > $intdir/apertures_INT.txt
    astmkprof $intdir/apertures_INT.txt --background=$gc --backhdu=1 \
        --clearcanvas --replace --type=int16 --mforflatpix \
        --mode=wcs --output=$intdir/apertures_rc.fits

    astmkcatalog $intdir/apertures_rc.fits -h1 --zeropoint=0 \
            --valuesfile=$gc --valueshdu=1 \
            --ids --ra --dec --sum \
            --output=$outcat/gc.cat
    echo done > $intdone
fi







# plot LBT stars with color
if [ -f build-"$gal"/plots/f_from_data.png ]; then
  echo " "
  echo "Plot done"
  echo " "
else
    python3.11 plot_lbt_fluxes.py "$gal"
fi





###########################################################################

#                           Save factor                                   #

###########################################################################
# read the theoretical correction factor
# Read the table and remove zeros, save the factor
# table subgruop is of sdss spectra stars

f=$(asttable $outcat/table_subgroup.txt --noblankend=4 | aststatistics --mean -c4)
std_f=$(asttable $outcat/table_subgroup.txt --noblankend=4 | aststatistics --std -c4)
echo $f $std_f> $outcat/factor_relationR-Ha.txt


ref=$bdir/pol_coefficients.txt
p0=$(awk 'NR=='1'{print $1}' $ref)
p1=$(awk 'NR=='1'{print $2}' $ref)
    
echo $p0 $p1
# Now here I want to save a table with F_LBT, F_Ha_real, F_expected
# Now the problem is that rc and ha cat are not selected on the base of the color
# but I am taking all the stars in the retrieved cat
# So I should first select the stars with g-r in the same range subgroup.







asttable $outcat/rc.cat -c1,2,3,4 --catcolumnfile=$outcat/halpha.cat -ot.txt
asttable t.txt -c1,2,3,4,8 -o t1.txt
asttable t1.txt -c1,2,3,4,5 -c'arith $4 '$p0' x '$p1' +' -ot2.txt
# for the istogram ( here stars with g < 0.5)
asttable t2.txt -c'arith $6 $5 /' -o conv_factors.txt
rm t.txt t1.txt

images_factor=$(aststatistics conv_factors.txt -c1 --sigclip-median)
im_std=$(aststatistics conv_factors.txt -c1 --sigclip-std)
echo $images_factor $im_std > $outcat/factor_relationR-Ha_myimages.txt

rm conv_factors.txt t2.txt

# did you check that all the stras ae also present int? They may be too faint!
###########################################################################

#                           Calibrate image                               #

###########################################################################
# Last, I have to multiply my Halpha image to have it at the same zeropoint
# of LBT_r ( 22.5)

echo $images_factor
echo $f
if [ -f $bdir/Ha_calibrated-"$gal".fits ]; then
  echo " "
  echo "H image calibrated"
  echo " "
else
    astarithmetic $ha -h1 $images_factor x -o $bdir/Ha_calibrated-"$gal".fits
    # Multiply lbt to have the r at the same level of h alpha + f
    # test of consistency
    astarithmetic $rc -h1 $f x -o $bdir/r_atHa.fits
    astarithmetic $bdir/Ha_calibrated-"$gal".fits -h1 \
                $bdir/r_atHa.fits -h1 - -o $bdir/Ha_sub-"$gal".fits
fi






/bin/bash save_prof.sh "$gal" $ra $dec $ar $pa $Ag $Ar $mdim


rm -r $bdir/INT-apertures
rm -r $bdir/LBT-apertures
rm -r $bdir/LBT-apertures-g
