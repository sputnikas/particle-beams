REM gnuplot -e "FILENAME='test2';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;XMAX=250;YMIN=-20;YMAX=20;WX=1416" xz.gp

REM for %%i in (test1.600.dat) do (
	REM gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='600';YMIN=0;WX=1416" rv.gp
REM )
REM for %%i in (test1.500.dat) do (
	REM gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='500';YMIN=0;WX=1416" rv.gp
REM )
REM for %%i in (test1.400.dat) do (
	REM gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='400';YMIN=0;WX=1416" rv.gp
REM )
for %%i in (test2.300.dat) do (
	gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='300';YMIN=0;XMAX=10;YMAX=0.15;WX=1416" rv.gp
)
for %%i in (test2.200.dat) do (
	gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='200';YMIN=0;XMAX=7;YMAX=0.15;WX=1416" rv.gp
)
REM for %%i in (test2.100.dat) do (
	REM gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='100';YMIN=0;XMAX=4;YMAX=0.15;WX=1416" rv.gp
REM )
REM for %%i in (test2.050.dat) do (
	REM gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='50';YMIN=0;XMAX=2;YMAX=0.1;WX=1416" rv.gp
REM )
REM for %%i in (test2.025.dat) do (
	REM gnuplot -e "FILENAME='%%~ni';TERM='png';FONTSIZE=37;XMULT=1000;XMIN=0;FILENAME2='25';YMIN=0;YMAX=0.05;XMAX=1.2;WX=1416" rv.gp
REM )

pause

