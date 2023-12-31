reset
yaxis=1 # 1 for mass coordinate, anything else for radius
xaxis=3 # 1 for Model Number, 2 for Linear time, else for Time to CC (log)
filename='history.data'

# set colour palette
convcolor='green'
overcolor='yellow'
semicolor='purple'
thercolor='grey'

# Preamble
hisfile=sprintf('<tail -n +6 %s',filename)

rsun=6.96e10

set term cairolatex standalone pdf size 8.in,5 font 'Times-Roman,12'

set style fill solid
set tics front
set linetype 1 lc rgb convcolor # convection zone
set linetype 2 lc rgb overcolor# overshoot
set linetype 3 lc rgb semicolor # semiconvection zone
set linetype 4 lc rgb thercolor # thermohaline mixing

# Set horizontal axis
if (xaxis==1){
    stats hisfile u (column('model_number')):($0) nooutput
    lastmodel=STATS_max_x
    firstmodel=STATS_min_x
    set xr [firstmodel:lastmodel]
    set xl 'Model Number'
    set form x '{%g}'
    xlabel='model_number'
    xfunc(t)=t

} else {
    stats hisfile u (column('star_age')):($0) nooutput
    age=STATS_max_x
    start=STATS_min_x
    if (xaxis==2){
	set xr [start/1e6:age/1e6]
	set xl 'Age [Myr]'
	set form x '{%g}'
	xlabel='star_age'
	xfunc(t)=t/1e6

    } else {
	tailhisfile=sprintf("<awk 'NR==6; {buffer[NR%%3]=$0} END{print buffer[(NR+1)%%3]; print buffer[(NR+2)%%3]}' %s",filename)
	STATS_min_y=-99
	stats tailhisfile u (column('star_age')):($0) nooutput
	mindt=log10(STATS_max_x-STATS_min_x)

	set xr [log10(age-start):mindt]
	set xl 'Time to CC [yr]'
	set form x '$10^{%g}$'
	xlabel='star_age'
	xfunc(t)=(t<age?log10(age-t):1/0)
    }
}

# Set vertical axis
shorthisfile=sprintf("<awk 'NR>=6&&NR<=10' %s",filename)

if (yaxis==1) {
    ### Mass Kippenhahn diagram ###
    stats [0:1e99] hisfile u ($0):(column('star_mass')) nooutput
    mass=STATS_max_y
    set yr [0:mass]
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

} else {
    ### Radius Kippenhahn diagram ###
    STATS_min_y=-99
    stats [0:1e99] shorthisfile u ($0):(column('log_R')) nooutput
    if (STATS_min_y>-50) {
	radius(x)=10**x
	rad='log_R'
    } else {
	stats [0:1e99] shorthisfile u ($0):(column('radius')) nooutput
	if (STATS_min_y>-50) {
	    radius(x)=x
	    rad='radius'
	} else {
	    stats [0:1e99] shorthisfile u ($0):(column('radius_cm')) nooutput
	    if (STATS_min_y>-50) {
		radius(x) = x/rsun
		rad='radius_cm'
	    } else {
		print 'Radius coordinate not found'
		pause -1
	    }
	}
    }
    stats [0:1e99] hisfile u ($0):(radius(column(rad))) nooutput
    maxrad=STATS_max_y

    set yr [1e-3:maxrad]
    set yl 'Radius [R$_\odot$]'
    set form y '$10^{%L}$'
    set log y

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
    STATS_max_y=-99
    cmd=sprintf("if grep -q %s %s; then echo '1';else  echo '2';fi",convqtop(i),filename)
    iconv=system(cmd)
    if (iconv>1){
	iconvmax=i-1
	break
    }
}

# Count number of burning zones in file
iburnmax=0
do for [i=1:1000]{
    STATS_max_y=-99
    cmd=sprintf("if grep -q %s %s; then echo '1';else  echo '2';fi",burnqtop(i),filename)
    iburn=system(cmd)
    if (iburn>1){
	iburnmax=i-1
	break
    }
}

# Start plotting --------------------------------------------------------------

    
# Burning
    unset colorbox
    set out 'burn.tex'
    set pal defined (-15 'blue',0 'white',15 'red')
    set cbr [-20:20]
    pl for [i=iburnmax:1:-1] hisfile \
       u (xfunc(column(xlabel))):\
       (column(burnqtop(i))<=1.?column(burnqtop(i))*upper(column(upperlbl)):1/0):\
       (column(burntype(i))>-100?column(burntype(i)):0) w boxes fc pal not,\
       '' u (xfunc(column(xlabel))):(upper(column(upperlbl))) w l lw 3 lc rgb 'black' not,\
       '' u (xfunc(column(xlabel))):(column(he_core)) w l lw 4 lc rgb 'forest-green' not,\
       '' u (xfunc(column(xlabel))):(column(co_core)) w l lw 4 lc rgb 'blue' not,\
       '' u (xfunc(column(xlabel))):(column(one_core)) w l lw 4 lc rgb 'red' not,\
       '' u (xfunc(column(xlabel))):(column(fe_core)) w l lw 4 lc rgb 'brown' not
    
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
system('convert kipppng-1.png -transparent white -channel Alpha -evaluate Divide 4 kipptrans.png')
system('convert burnpng-1.png kipptrans.png -composite -format png hoge.png')
system('convert hoge.png -transparent white burn-inc.png')
system('convert burn-inc.png -format pdf burn-inc.pdf')
system('pdflatex burn.tex && mv burn.pdf kippdiagram.pdf')
system('rm burn* *-inc* kipp*png hoge* *tex')
