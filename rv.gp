if (exists("FILENAME")) {

	if (!exists("FONTSIZE")) { FONTSIZE = '24' }
	if (!exists("WX")) { WX = 944 }
	if (!exists("WY")) { WY = 944 }
	
	if (!exists("MULT")) { MULT = 100 }
	if (!exists("XMULT")) { XMULT = MULT }
	if (!exists("YMULT")) { YMULT = 1 }
	
	if (!exists("YLABEL")) { YLABEL = '{/Times:Italic v/c}'}
	
	if (!exists("XLABEL") && (XMULT == 1000)) { XLABEL = '{/Times:Italic r}, мм' }
	if (!exists("XLABEL") && (XMULT == 100)) { XLABEL = '{/Times:Italic r}, см' }
	if (!exists("XLABEL") && (XMULT == 10)) { XLABEL = '{/Times:Italic r}, дм' }
	if (!exists("XLABEL") && (XMULT == 1)) { XLABEL = '{/Times:Italic r}, м' }
	
	if (!exists("XMIN")) { XMIN = 0 }
	if (!exists("IW")) { IW = 0.1}
	if (!exists("IC")) { IC = 1.5}
	if (!exists("RA0")) {RA0 = 2}
	if (!exists("RC0")) {RC0 = 1}
	if (!exists("RAL")) {RAL = 1.6}
	if (!exists("RCL")) {RCL = 1.4}
	if (!exists("L")) {L = 10}

	if (!exists("TERM")) { TERM = 'png'}
	
	CL = 2.99792458E8

	print(TERM)
	if (!exists("TERM") || (TERM eq "wxt")) {
		set term wxt 1 size WX,WY enhanced font 'Times New Roman,'.FONTSIZE
	}
	if (TERM eq "png") {
		set term pngcairo size WX,WY enhanced font 'Times,'.FONTSIZE
		set output FILENAME.'_zvz.png'
	}
	if (TERM eq "pdf") {
		set term cairolatex pdf size 24cm,12cm font 'ptm,bx'
		set output FILENAME.'_zvz.pdf'
	}
	if (TERM eq "gif") {
		set term gif font 'Times,'.FONTSIZE animate transparent opt delay 3 size WX,WY 
		set output FILENAME.'_zvz.gif'
		stats FILENAME.'.dat' nooutput
	}

	stats FILENAME.'.dat' using 4 name 'x' nooutput
	if (!exists("XMAX")) { XMAX = 1.01*x_max*XMULT  }
	
	stats FILENAME.'.dat' using 7 name 'y' nooutput
	if (!exists("YMAX")) {YMAX = 1.1*y_max*YMULT  }
	
	if (!exists("YMIN")) { YMIN = - YMAX }
	#print XMAX
	#print floor(XMAX)
	#XMAX = (floor(XMAX) > 100) ? (floor(XMAX/10) + 1)*10 : ((floor(XMAX) > 10) ? (floor(XMAX/5) + 1)*5 : floor(XMAX) + 1)
	#XMAX = (floor(XMAX) > 100) ? (floor(XMAX/20) + 1)*20 : ((floor(XMAX) > 50) ? (floor(XMAX/10) + 1)*10 :((floor(XMAX) > 25) ? (floor(XMAX/5) + 1)*5 : (floor(XMAX) > 10) ? (floor(XMAX/2) + 1)*2 : floor(XMAX)+1))
	#XTICS = (floor(XMAX) > 100) ? 20 : ((floor(XMAX) > 50) ? 10 :((floor(XMAX) > 25) ? 5 : (floor(XMAX) > 10) ? 2 : 1))
	unset key
	set style line 12 lc rgb 'black' dt 2 lw 2
	set style line 13 lc rgb 'black' dt 3 lw 1
	set grid xtics ytics mxtics mytics ls 12, ls 13
	set mxtics 2
	set mytics 2
	set xrange [XMIN: XMAX]
	set yrange [YMIN: YMAX]
	#set xtics XTICS
	#set ytics 0.2
	set xlabel XLABEL
	set ylabel YLABEL
	set parametric
	#set xtics format ""
	#set xtics add ("-2" -2, "-1,5" -1.5, "-1" -1, "-0,5" -0.5, "0" 0, "0,5" 0.5, "1" 1, "1,5" 1.5, "2" 2)
	#set ytics format ""
	#set ytics add ("-2" -2, "-1,5" -1.5, "-1" -1, "-0,5" -0.5, "0" 0, "0,5" 0.5, "1" 1, "1,5" 1.5, "2" 2)
	
	#print FILENAME2
	
	if (exists("FILENAME2")) {
		plot FILENAME.'.dat' using (sqrt($3**2 + $2**2 + $1**2)*XMULT):(sqrt($4**2 + $6**2 + $5**2)*YMULT) w p lt 1 lw 3 lc 'black',\
			FILENAME2.'.gpdat' using ($1*XMULT):2 w l lt 1 lw 4 lc 'red'
	} else {
		plot FILENAME.'.dat' using (sqrt($3**2 + $2**2 + $1**2)*XMULT):(sqrt($4**2 + $6**2 + $5**2)*YMULT) w p lt 1 lw 3 lc 'black'
	}
} else {
	print 'You not send FILENAME!'
}

