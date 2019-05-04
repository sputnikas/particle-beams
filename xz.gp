if (exists("FILENAME")) {
	if (!exists("FONTSIZE")) { FONTSIZE = '24' }
	if (!exists("WX")) { WX = 944 }
	if (!exists("WY")) { WY = 944 }
	
	if (!exists("MULT")) { MULT = 1000 }
	if (!exists("XMULT")) { XMULT = MULT }
	if (!exists("YMULT")) { YMULT = MULT }
	
	if (!exists("XLABEL") && (XMULT == 1000)) { XLABEL = '{/Times:Italic z}, мм' }
	if (!exists("YLABEL") && (YMULT == 1000)) { YLABEL = '{/Times:Italic x}, мм' }
	if (!exists("XLABEL") && (XMULT == 100)) { XLABEL = '{/Times:Italic z}, см' }
	if (!exists("YLABEL") && (YMULT == 100)) { YLABEL = '{/Times:Italic x}, см' }
	if (!exists("XLABEL") && (XMULT == 10)) { XLABEL = '{/Times:Italic z}, дм' }
	if (!exists("YLABEL") && (YMULT == 10)) { YLABEL = '{/Times:Italic x}, дм' }
	if (!exists("XLABEL") && (XMULT == 1)) { XLABEL = '{/Times:Italic z}, м' }
	if (!exists("YLABEL") && (YMULT == 1)) { YLABEL = '{/Times:Italic x}, м' }
	
	if (!exists("XMIN")) { XMIN = -5 }
	if (!exists("XMAX")) { XMAX = 5  }
	if (!exists("YMIN")) { YMIN = -5 }
	if (!exists("YMAX")) { YMAX = 5  }
	if (!exists("IW")) { IW = 0.1}
	if (!exists("IC")) { IC = 1.5}

	if (!exists("TERM")) { TERM = 'png'}

	print(TERM)
	if (!exists("TERM") || (TERM eq "wxt")) {
		set term wxt 1 size WX,WY enhanced font 'Times New Roman,'.FONTSIZE
	}
	if (TERM eq "png") {
		set term pngcairo size 944,944 enhanced font 'Times,'.FONTSIZE
		set output FILENAME.'.png'
	}
	if (TERM eq "pdf") {
		set term cairolatex pdf size 24cm,12cm font 'ptm,bx'
		set output FILENAME.'.pdf'
	}
	if (TERM eq "gif") {
		set term gif font 'Times,'.FONTSIZE animate transparent opt delay 3 size WX,WY 
		set output FILENAME.'.gif'
		stats FILENAME.'.gpdat' nooutput
	}

	unset key
	set grid
	set xrange [XMIN: XMAX]
	set yrange [YMIN: YMAX]
	set xlabel XLABEL
	set ylabel YLABEL
	set parametric
	#set xtics format ""
	#set xtics add ("-2" -2, "-1,5" -1.5, "-1" -1, "-0,5" -0.5, "0" 0, "0,5" 0.5, "1" 1, "1,5" 1.5, "2" 2)
	#set ytics format ""
	#set ytics add ("-2" -2, "-1,5" -1.5, "-1" -1, "-0,5" -0.5, "0" 0, "0,5" 0.5, "1" 1, "1,5" 1.5, "2" 2)
	
	if (TERM eq "gif") {
		print 'STATS_blocks = '.STATS_blocks
		print 'STATS_blank = '.STATS_blank
		do for [i=1:int(STATS_blocks)] {
			plot FILENAME.'.gpdat' index (i-1) using ($1*XMULT):($2*YMULT) w p lt 1 lc 'black',\
				 [t=0:2*pi] 2*cos(t), 2*sin(t) w l dt 1 lw 4 lc 'black',\
				 [t=0:2*pi] cos(t), sin(t) w l dt 1 lw 4 lc 'black',\
				 [t=0:2*pi] IC + IW*t/2/pi, 0 w l dt 1 lw 4 lc 'red'
		}
	} else {
		plot FILENAME.'.dat' using ($4*XMULT):($2*YMULT) w p lt 1 lc 'black'
			 
	}
} else {
	print 'You not send FILENAME!'
}

