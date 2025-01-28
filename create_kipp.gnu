# Script for creating Kippenhahn diagrams from a MESA history file

# Owned by: Ryosuke Hirai

# Prerequisites: pdflatex, ImageMagick, gnuplot, pdftoppm

#  Most prerequisites can be installed via standard repositories
#  Linux example) $ sudo dnf install ImageMagick gnuplot poppler-utils
#  MacOS example) $ brew install imagemagick gnuplot poppler

# Required columns in the history file
# 1. 'mass' for mass Kippenhahn diagram
#    'log_R' or 'radius' or 'radius_cm' for radius Kippenhahn diagram
# 2. 'burning_regions <integer>' or 'burn_relr_regions <integer>'
# 3. 'mixing_regions <integer>' or 'mix_relr_regions <integer>'
# 4. (optional) 'he_core_mass' or 'he_core_radius' and similar for co_core, one_core, fe_core

# Settings ===================================================================
yaxis=1 # 1 for mass coordinate, 2 for log radius, anything else for linear radius
xaxis=2 # 1 for Model Number, 2 for Linear time, 3 for log time, else for Time to end (log)
filename='history.data' # Input history file name
outfile='kippdiagram.pdf' # Name of output pdf file

age_units=2 # 1 for yr, 2 for Myr (only used when xaxis=2)
minrad=1e-3 # Minimum radius in units of Rsun for radius Kipp diagram (only used when yaxis=2)
marginfrac=1.03 # How much margin to leave above the stellar mass/radius
magick='convert' # 'convert' for ImageMagick v6, 'magick' for v7

# Plot ranges (negative for default)
leftbound=-1  # left bound for xaxis
rightbound=-1 # right bound for xaxis

# set colour palette
convcolor='green'  # convection zone
overcolor='yellow' # overshoot
semicolor='purple' # semiconvection zone
thercolor='grey'   # thermohaline mixing

# ===========================================================================
reset

# Preamble
hisfile=sprintf('<tail -n +6 %s',filename)
print "Reading file: ",filename

rsun=6.96e10
colexist(col,file) = sprintf("if grep -q %s %s; then echo '1';else  echo '2';fi",col,file) # command to tell if a column exists

set style fill solid
set tics front
set linetype 1 lc rgb convcolor # convection zone
set linetype 2 lc rgb overcolor # overshoot
set linetype 3 lc rgb semicolor # semiconvection zone
set linetype 4 lc rgb thercolor # thermohaline mixing

# Set horizontal axis
if (xaxis==1){
    stats hisfile u (column('model_number')):($0) nooutput
    if( leftbound<0 ) { leftbound =STATS_min_x }
    if(rightbound<0 ) { rightbound=STATS_max_x }
    set xr [leftbound:rightbound]
    set xl 'Model Number'
    xlabel='model_number'
    xfunc(t)=t

    print "Using 'model_number' as x-axis"

} else {
    stats hisfile u (column('star_age')):($0) nooutput
    if( leftbound<0 ) { leftbound =STATS_min_x }
    if(rightbound<0 ) { rightbound=STATS_max_x }

    # Set age units
    if (age_units==1) {
	ageunit=1
	ageunitlabel='yr'
    } else {
	ageunit=1e6
	ageunitlabel='Myr'
    }

    if (xaxis==2||xaxis==3){
	set xr [leftbound/ageunit:rightbound/ageunit]
	set xl sprintf('Age [%s]',ageunitlabel)
	xlabel='star_age'
	xfunc(t)=t/1e6

	if(xaxis==3){set log x;set form x '$10^{%g}$'}

	print "Using 'star_age' as x-axis"

    } else {
	tailhisfile=sprintf("<awk 'NR==6; {buffer[NR%%3]=$0} END{print buffer[(NR+1)%%3]; print buffer[(NR+2)%%3]; print buffer[(NR+3)%%3]}' %s",filename)
	STATS_min_y=-99
	stats tailhisfile u (column('star_age')):($0) nooutput
	mindt=log10((STATS_max_x-STATS_min_x)/2.)

	set xr [log10(rightbound-leftbound):mindt]
	set xl 'Time to end [yr]'
	set form x '$10^{%g}$'
	xlabel='star_age'
	xfunc(t)=(t<rightbound?log10(rightbound-t):1/0)

	print "Using 'Time to end' as x-axis"
	
    }
}

# Set vertical axis
shorthisfile=sprintf("<awk 'NR>=6&&NR<=10' %s",filename)

if (yaxis==1) {
    ### Mass Kippenhahn diagram ###
    stats [0:1e99] hisfile u ($0):(column('star_mass')) nooutput
    mass=STATS_max_y
    set yr [0:mass*marginfrac]
    set yl 'Mass [M$_\odot$]'

    convqtop(i)=sprintf('mix_qtop_%d',i)
    convtype(i)=sprintf('mix_type_%d',i)
    burnqtop(i)=sprintf('burn_qtop_%d',i)
    burntype(i)=sprintf('burn_type_%d',i)

    he_core='he_core_mass'
    co_core='co_core_mass'
    one_core='one_core_mass'
    fe_core='fe_core_mass'

    upper(x)=x
    upperlbl='star_mass'

    print "Using 'mass coordinate' as y-axis"

} else {
    ### Radius Kippenhahn diagram ###
#    STATS_min_y=-99
    #    stats [0:1e99] shorthisfile u ($0):(column('log_R')) nooutput
    logr_exist=system(colexist('log_R',filename))
    if (logr_exist==1) {
	radius(x)=10**x
	rad='log_R'
	print "Using 'radius' as y-axis"
    } else {
	rad_exist=system(colexist('radius',filename))
	if (rad_exist) {
	    radius(x)=x
	    rad='radius'
	} else {
	    rcm_exist=system(colexist('radius',filename))
	    if (rcm_exist) {
		radius(x) = x/rsun
		rad='radius_cm'
	    } else {
		print 'Radius coordinate not found'
		pause -1
	    }
	}
    }
    print "Using 'radius' as y-axis"
    stats [0:1e99] hisfile u ($0):(radius(column(rad))) nooutput
    maxrad=STATS_max_y

    set yr [minrad:maxrad*marginfrac]
    set yl 'Radius [R$_\odot$]'
    if( yaxis==2 ) {
	set log y
	set form y '$10^{%L}$'
    }

    convqtop(i)=sprintf('mix_relr_top_%d',i)
    convtype(i)=sprintf('mix_relr_type_%d',i)
    burnqtop(i)=sprintf('burn_relr_top_%d',i)
    burntype(i)=sprintf('burn_relr_type_%d',i)

    he_core='he_core_radius'
    co_core='co_core_radius'
    one_core='one_core_radius'
    fe_core='fe_core_radius'

    upper(x)=radius(x)
    upperlbl=rad

}

# Count number of convective zones in file

iconvmax=0
do for [i=1:1000]{
    iconv=system(colexist(convqtop(i),filename))
    if (iconv>1){
	iconvmax=i-1
	break
    }
}

print "Maximum number of convective zones is iconvmax=",iconvmax

# Count number of burning zones in file
iburnmax=0
do for [i=1:1000]{
    iburn=system(colexist(burnqtop(i),filename))
    if (iburn>1){
	iburnmax=i-1
	break
    }
}

print "Maximum number of burning zones is iburnmax=",iburnmax

# Check if core masses/radii are outputted
ihecore =system(colexist(he_core,filename))
icocore =system(colexist(co_core,filename))
ionecore=system(colexist(one_core,filename))
ifecore =system(colexist(fe_core,filename))

# Start plotting --------------------------------------------------------------

print "Plotting burning/cooling zones..."
# Burning
unset colorbox

# Find maximum burning rate
maxburnrate=0
do for [i=1:iburnmax] {
    stats [*:*][*:*] hisfile u 0:(column(burntype(i))>-100?abs(column(burntype(i))):0) nooutput
    maxburnrate=(maxburnrate>STATS_max_y?maxburnrate:STATS_max_y)
}

set term cairolatex standalone pdf size 8.in,5 font 'Times-Roman,12'
set out 'burn.tex'

set pal defined (-maxburnrate 'blue',0 'white',maxburnrate 'red')
set cbr [-maxburnrate:maxburnrate]
pl for [i=iburnmax:1:-1] hisfile \
   u (xfunc(column(xlabel))):\
   (column(burnqtop(i))<=1.?column(burnqtop(i))*upper(column(upperlbl)):1/0):\
   (column(burntype(i))>-100?column(burntype(i)):0) w boxes fc pal not,\
   '' u (xfunc(column(xlabel))):(upper(column(upperlbl))) w l lw 3 lc rgb 'black' not,\
   '' u (xfunc(column(xlabel))):(ihecore==1?column(he_core):NaN) w l lw 4 lc rgb 'forest-green' not,\
   '' u (xfunc(column(xlabel))):(icocore==1?column(co_core):NaN) w l lw 4 lc rgb 'blue' not,\
   '' u (xfunc(column(xlabel))):(ionecore==1?column(one_core):NaN) w l lw 4 lc rgb 'red' not,\
   '' u (xfunc(column(xlabel))):(ifecore==1?column(fe_core):NaN) w l lw 4 lc rgb 'brown' not

print "Plotting mixing zones..."
# Convection zones
    set out 'kipp.tex'
    set pal defined (0 'white',1 convcolor,2 overcolor,3 semicolor,4 thercolor)
    set cbr [0:4]
    pl  for [i=iconvmax:1:-1]hisfile \
    u (xfunc(column(xlabel))):\
    (column(convqtop(i))<=1?column(convqtop(i))*upper(column(upperlbl)):1/0):\
    (column(convtype(i)) >0?column(convtype(i)):1/0) w boxes fc pal not


# compile tex file
set out
system('pdftoppm burn-inc.pdf burnpng -png')
system('pdftoppm kipp-inc.pdf kipppng -png')
system(sprintf('%s kipppng-1.png -transparent white -channel Alpha -evaluate Divide 4 kipptrans.png',magick))
system(sprintf('%s burnpng-1.png kipptrans.png -composite -format png hoge.png',magick))
system(sprintf('%s hoge.png -transparent white burn-inc.png',magick))
system(sprintf('%s burn-inc.png -format pdf burn-inc.pdf',magick))
print "Adding axis labels with latex..."
system('pdflatex -interaction=nonstopmode burn.tex > /dev/null && mv burn.pdf tempfile.pdf')
system('rm -f burn* *-inc* kipp*png hoge* *tex')
system(sprintf('mv tempfile.pdf %s 2>/dev/null',outfile))

print "Outputted ",outfile
