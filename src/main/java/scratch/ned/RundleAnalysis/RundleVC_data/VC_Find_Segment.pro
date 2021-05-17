;
;	PROGRAM VC_Find_Segment.PRO
;
;	This program allows the user to find a particular segment and shows
;		it on a map .....
;	
;	Some Admonitions:
;
;	1.  Make sure the dimensions here match up with those in the
;		parameter file EQPARM_4.FOR
;	2.  You need to have the IDL package of software visualization
;		codes for these scripts to work.
;	3.  The data files you will need are:
;
;		i.   The output from VC simulator
;		ii.  The fault data file that describes the fault
;			system topology
;
;	*****************************************************
;
;

	fname1 		=''
	fn1		=''
	fname2		=''
	fname3		=''
	fname5		=''
	fn5 		=''
;
	rmenu		=''
	resp		=''
	respfb		=''
;
;	******************************************************
;
;	Set up color printer
;	
 	set_plot,'ps'
        device,/color,bits=8
	
;
;     Query for output file from EQSYN
;
      Print, ' '
      print, ' Enter the name of the fault file'
      print, ' for VC Model'

      read,fname2

      openr, 10, fname2
 ;

      readf, 10, nfault, mend
;
     
	NFT=11000
	NTSTP=11000
	NF = nfault
	NF1= nfault+1
	NF2= 800
	NFF = 640000
	NMF = 9600000
	NPT = 40000

	db		= fltarr(NF)
	dt		= fltarr(NF)
	xfe		= fltarr(NF)
	yfe		= fltarr(NF)
	xfw		= fltarr(NF)
	yfw		= fltarr(NF)
	slpvel		= fltarr(NF)
	nseg		= intarr(NF)
	data		= fltarr(7)

		for i=0L,nfault-1 do begin
		readf,10,data
		db(i) = data(0)
		dt(i) = data(1)
		xfe(i)= data(2)
		yfe(i)= data(3)
		xfw(i)= data(4)
		yfw(i)= data(5)
		slpvel(i) = data(6)
		endfor

	readf, 10, vplatx, vplaty

	close, 10

;
      rmenu = 'A'

      while (rmenu ne 'Z') do begin
;
      print, ' '
      print, ' A. FIND A PARTICULAR FAULT SEGMENT & PLOT ON FAULT MAP'
      print, ' '

      print, ' Z. EXIT PROGRAM'
      print, ' '
;
      read, format ='(a1)', rmenu
;

        if (rmenu eq 'A') then begin
        find_segment, nfault, vplatx, vplaty, xfe, xfw, yfe, yfw
        endif


;
;	End of the Main "While" Loop
;

        endwhile
;
	end



