package scratch.aftershockStatisticsETAS;

import java.awt.Font;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.ListIterator;
import java.util.TimeZone;

import javax.swing.JTable;

import org.apache.commons.io.FileUtils;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

// template SVG file

public class GraphicalForecast{

	private static boolean D = false;
	
	private String location = "somewhere";
	private double lat0 = -33.333;
	private double lon0 = -123.333;
	private double mag0 = 9.9;
	private double depth0 = 0.0;
	private double tPredStart = 0.67;
	private GregorianCalendar forecastStartDate = new GregorianCalendar();
	private GregorianCalendar eventDate = new GregorianCalendar();
	private String msName = "MS_NAME";
	private String versionNumber = "";
	private GregorianCalendar versionDate = new GregorianCalendar();
	private ETAS_AftershockModel aftershockModel;
	private File outFile;
	private double b;
	private String shakeMapURL = null;
	private String eventURL = "#";
	
	public GraphicalForecast(){
		// run with dummy data. This is all made up.
		eventDate.set(1980, 3, 4, 11, 59, 59);
		
		number = new double[predictionMagnitudes.length][predictionIntervals.length][3];
		probability = new double[predictionMagnitudes.length][predictionIntervals.length];
		
		double[] fractiles = new double[3];
		
		double baseRate = 1000;
		for (int i = 0; i < predictionMagnitudes.length; i++){
			for (int j = 0; j < predictionIntervals.length; j++){
				double rate = baseRate*(j+1d)/Math.pow(10d, i) ;
				probability[i][j] = 1 - Math.exp( -rate );
				
				fractiles[0] = rate;
				fractiles[1] = rate - Math.sqrt(rate);
				fractiles[2] = rate*10d;
				
				number[i][j][0] = fractiles[0];
				number[i][j][1] = fractiles[1];
				number[i][j][2] = fractiles[2];
			}
		}
		
		// do weekly probabilities of felt and damaging quakes
		double rate = baseRate*(1+1d)/Math.pow(10d, 3);
		pDamage_wk = 1 - Math.exp( -rate );
		
		rate = baseRate*(1+1d)/Math.pow(10d, 0);
		fractiles[0] = rate;
		fractiles[1] = rate - Math.sqrt(rate);
		fractiles[2] = rate*10d;
		
		nfelt_wk[0] = Math.max(fractiles[0],0);
		nfelt_wk[1] = Math.max(fractiles[1],0);
		nfelt_wk[2] = Math.max(fractiles[2],1);
	}
	
	public GraphicalForecast(File outFile, ETAS_AftershockModel aftershockModel, GregorianCalendar eventDate,
			GregorianCalendar startDate) {
		
		startDate.setTimeZone(TimeZone.getTimeZone("UTC"));
		eventDate.setTimeZone(TimeZone.getTimeZone("UTC"));

		this.forecastStartDate = startDate;
		this.eventDate = eventDate;
		this.aftershockModel = aftershockModel;
		this.outFile = outFile;
		this.msName = aftershockModel.mainShock.getEventId(); //get USGS name
		this.eventURL = "https://earthquake.usgs.gov/earthquakes/eventpage/" + msName + "#executive";
		versionDate.setTimeInMillis(System.currentTimeMillis());
		this.location = getStringParameter(aftershockModel.mainShock, "description");
		this.lat0 = aftershockModel.mainShock.getHypocenterLocation().getLatitude();
		this.lon0 = aftershockModel.mainShock.getHypocenterLocation().getLongitude();
		this.mag0 = aftershockModel.mainShock.getMag();
		this.depth0 = aftershockModel.mainShock.getHypocenterLocation().getDepth();
		this.tPredStart = getDateDelta(eventDate, forecastStartDate);
		this.b = aftershockModel.get_b();
		
	}

	private String getStringParameter(ObsEqkRupture rup, String paramName){
		ListIterator<?> iter = rup.getAddedParametersIterator();
		if (iter != null) {
			while (iter.hasNext()){
				Parameter<?> param = (Parameter<?>) iter.next();
				if (param.getName().equals(paramName));
				return param.getValue().toString();
			}
		} 
		return "";
	}
	
	private final static int DAY = 1;
	private final static int WEEK = 7;
	private final static int MONTH = 30;
	private final static int YEAR = 365;
		
	
	//	parameters for word wrap algorithm and text summary
	private boolean smartRoundPercentages = true;   //give rounded range of probabilities?
	private double feltMag = 3;  //% magnitude to use for calculating number of 'felt' earthquakes
	private double damageMag = 6.0;   // magnitude to use for damaging earthquake
	private double pDamage_wk = 0;
	private double[] nfelt_wk = new double[3];
	
	
	private double[] predictionMagnitudes = new double[]{3,4,5,6,7,8,9};
	private double[] predictionIntervals = new double[]{DAY,WEEK,MONTH,YEAR}; //day,week,month,year
//	private String[] predictionIntervalStrings = new String[]{"day","week","month","year"}; //day,week,month,year
	private double[][][] number;  //dimensions: [predMag][predInterval][range: exp, lower, upper];
	private String[][] numberString;
	private double[][] probability;
	private String[][] probString;
	private String forecastHorizon;
	private String spatialForecastInterval = "week";
	
	
	private HashMap<String, String> tags = new HashMap<String, String>();

	public void constructForecast(){
		number = new double[predictionMagnitudes.length][predictionIntervals.length][3];
		probability = new double[predictionMagnitudes.length][predictionIntervals.length];
		
		double tMinDays = getDateDelta(eventDate, forecastStartDate);
		double[] calcFractiles = new double[]{0.5,0.025,0.975};
		double[] fractiles = new double[3];
		double tMaxDays;
		
		for (int i = 0; i < predictionMagnitudes.length; i++){
			for (int j = 0; j < predictionIntervals.length; j++){
				tMaxDays = tMinDays + predictionIntervals[j];
				
				fractiles = aftershockModel.getCumNumFractileWithAleatory(calcFractiles, predictionMagnitudes[i], tMinDays, tMaxDays);
				number[i][j][0] = fractiles[0];
				number[i][j][1] = fractiles[1];
				number[i][j][2] = fractiles[2];
				probability[i][j] = aftershockModel.getProbabilityWithAleatory(predictionMagnitudes[i], tMinDays, tMaxDays);
			}
		}
		
		// do weekly probabilities of felt and damaging quakes
		tMaxDays = tMinDays + WEEK;
		pDamage_wk = aftershockModel.getProbabilityWithAleatory(damageMag, tMinDays, tMaxDays);
		fractiles = aftershockModel.getCumNumFractileWithAleatory(calcFractiles, feltMag, tMinDays, tMaxDays);
		nfelt_wk[0] = Math.max(fractiles[0],0);
		nfelt_wk[1] = Math.max(fractiles[1],0);
		nfelt_wk[2] = Math.max(fractiles[2],0);
		
		processForecastStrings();
		
		setForecastHorizon();
		
		assignForecastStrings();
		
//		writeHTML(outFile); //call it explicitly
	}
	
	private void processForecastStrings(){
		
		// numbers and probabilities
		// extract numbers and probabilities from ETAS_forecast
		numberString = new String[predictionMagnitudes.length][predictionIntervals.length];
		probString = new String[predictionMagnitudes.length][predictionIntervals.length];
		
		for (int i = 0; i < predictionMagnitudes.length; i++){
			for (int j = 0; j < predictionIntervals.length; j++){
				if (number[i][j][1] == 0 && number[i][j][2] == 0)
					numberString[i][j] = "0*";
				else
					numberString[i][j] = numberRange(number[i][j][1], number[i][j][2]);

				if (probability[i][j] < 0.001)
					probString[i][j] = "";
				else if (probability[i][j] <= 0.005)
					probString[i][j] = "<1%";
//				else if (probability[i][j] < 0.01)
//					probString[i][j] = String.format("%2.1f%%", 100*probability[i][j]);
				else if (probability[i][j] > 0.99)
					probString[i][j] = ">99%";
				else
					probString[i][j] = String.format("%1.0f%%", 100*probability[i][j]);
				
				if(D){System.out.println("M" + predictionMagnitudes[i] + " " + predictionIntervals[j] 
						+ " " + String.format("%5.4f", probability[i][j]) +" " + probString[i][j] + " " + numberString[i][j]);}
			}
		}
	}

	private void setForecastHorizon(){
		if (number[0][3][0] >= number[0][2][0] + 1 || //number[predMagIndex][predIntervalIndex][median, lower,  upper]
				number[0][3][2] >= number[0][2][2] + 1)
		    forecastHorizon = "year";
		else if (number[0][2][0] >= number[0][2][0] + 1 ||
				number[0][3][2] >= number[0][2][2] + 1)
			forecastHorizon = "month";
		else
			forecastHorizon = "week";
	}
	
	// Set all the variable bits of text 
	private void assignForecastStrings(){
		tags.put("N1_DA", numberString[0][0]);
		tags.put("P1_DA", probString[0][0]);
		tags.put("N1_WK", numberString[0][1]);
		tags.put("P1_WK", probString[0][1]);
		tags.put("N1_MO", numberString[0][2]);
		tags.put("P1_MO", probString[0][2]);
		tags.put("N1_YR", numberString[0][3]);
		tags.put("P1_YR", probString[0][3]);
		
		tags.put("N2_DA", numberString[1][0]);
		tags.put("P2_DA", probString[1][0]);
		tags.put("N2_WK", numberString[1][1]);
		tags.put("P2_WK", probString[1][1]);
		tags.put("N2_MO", numberString[1][2]);
		tags.put("P2_MO", probString[1][2]);
		tags.put("N2_YR", numberString[1][3]);
		tags.put("P2_YR", probString[1][3]);
		
		tags.put("N3_DA", numberString[2][0]);
		tags.put("P3_DA", probString[2][0]);
		tags.put("N3_WK", numberString[2][1]);
		tags.put("P3_WK", probString[2][1]);
		tags.put("N3_MO", numberString[2][2]);
		tags.put("P3_MO", probString[2][2]);
		tags.put("N3_YR", numberString[2][3]);
		tags.put("P3_YR", probString[2][3]);
		
		tags.put("N4_DA", numberString[3][0]);
		tags.put("P4_DA", probString[3][0]);
		tags.put("N4_WK", numberString[3][1]);
		tags.put("P4_WK", probString[3][1]);
		tags.put("N4_MO", numberString[3][2]);
		tags.put("P4_MO", probString[3][2]);
		tags.put("N4_YR", numberString[3][3]);
		tags.put("P4_YR", probString[3][3]);
		
		tags.put("N5_DA", numberString[4][0]);
		tags.put("P5_DA", probString[4][0]);
		tags.put("N5_WK", numberString[4][1]);
		tags.put("P5_WK", probString[4][1]);
		tags.put("N5_MO", numberString[4][2]);
		tags.put("P5_MO", probString[4][2]);
		tags.put("N5_YR", numberString[4][3]);
		tags.put("P5_YR", probString[4][3]);
		
		tags.put("N6_DA", numberString[5][0]);
		tags.put("P6_DA", probString[5][0]);
		tags.put("N6_WK", numberString[5][1]);
		tags.put("P6_WK", probString[5][1]);
		tags.put("N6_MO", numberString[5][2]);
		tags.put("P6_MO", probString[5][2]);
		tags.put("N6_YR", numberString[5][3]);
		tags.put("P6_YR", probString[5][3]);
		
		tags.put("N7_DA", numberString[6][0]);
		tags.put("P7_DA", probString[6][0]);
		tags.put("N7_WK", numberString[6][1]);
		tags.put("P7_WK", probString[6][1]);
		tags.put("N7_MO", numberString[6][2]);
		tags.put("P7_MO", probString[6][2]);
		tags.put("N7_YR", numberString[6][3]);
		tags.put("P7_YR", probString[6][3]);
		
//		tags.put("M1_R", String.format("%2.1f", predictionMagnitudes[0]));
//		tags.put("M2_R", String.format("%2.1f", predictionMagnitudes[1]));
//		tags.put("M3_R", String.format("%2.1f", predictionMagnitudes[2]));
		
		tags.put("NFELT_WK", numberRange(nfelt_wk[1], nfelt_wk[2]));
		
		if (smartRoundPercentages) {
			if (pDamage_wk < 0.001)
				tags.put("PDAMAGE_WK", "much less than 1");
			else if (pDamage_wk < 0.01)
				tags.put("PDAMAGE_WK", "less than 1");
			else if (pDamage_wk < 0.05)
				tags.put("PDAMAGE_WK", "1-5");
			else if (pDamage_wk < 0.95)
				tags.put("PDAMAGE_WK", String.format("%d", smartRound(Math.floor(pDamage_wk*100/5)*5)) + " - " +
						String.format("%d", (int) Math.ceil(pDamage_wk*100/5)*5));
			else if (pDamage_wk < 0.99)
				tags.put("PDAMAGE_WK", "greater than 95");
			else if (pDamage_wk <= 1.0) 
				tags.put("PDAMAGE_WK", "greater than 99");
			else
				tags.put("PDAMAGE_WK", "NaN");
		} else {
			String formatStr;
			if (pDamage_wk < 0.01)
				formatStr = "%2.1f";
			else
				formatStr = "%1.0f";
			tags.put("PDAMAGE_WK", String.format(formatStr, pDamage_wk*100));
		}

		// spatial forecast time interval
		if (spatialForecastInterval.equals("week"))
			tags.put("F_PLOT_T", "7");
		else if (spatialForecastInterval.equals("month"))
			tags.put("F_PLOT_T", "30");
		else if (spatialForecastInterval.equals("year"))
			tags.put("F_PLOT_T", "365");
	
		//	lat lon of mainshock
		String degSym = "&deg;";
		String tag;
		if (lat0 >= 0)
			tag = degSym + "N";
		else
			tag = degSym + "S";
		tags.put("MS_LAT", String.format("%4.3f" + tag, Math.abs(lat0)));
		
		lon0 = unwrap(lon0);
		
		if (lon0 >= 0) 
			tag = degSym + "E";
		else
			tag = degSym + "W";
		tags.put("MS_LON", String.format("%4.3f" + tag, Math.abs(lon0)));
	
		//		construct descriptive text
		SimpleDateFormat formatter=new SimpleDateFormat("d MMM yyyy, HH:mm:ss");  
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
		
		tags.put("MS_MAG", String.format("M%2.1f", mag0));
		tags.put("MS_NAME", msName);
		tags.put("MS_DATETIME", formatter.format(eventDate.getTime()));
		tags.put("MS_DEPTH", String.format("%2.1f km", depth0));
		tags.put("V_NUM", versionNumber);
		tags.put("V_DATE", formatter.format(versionDate.getTime()));
		
		tags.put("F_START_ABS", formatter.format(forecastStartDate.getTime()));
		
		String MS_LOC = location;
		String F_START_REL_DAYS = String.format("%2.1f", tPredStart);
		String DAMAGE_MAG = String.format("%.0f",  damageMag);
		String F_HORIZON = forecastHorizon;
		String DAYS;
		
		if (tPredStart == 1)
		    DAYS = "day";
		else
		    DAYS = "days";
		
		String DESCRIPTIVE_TEXT = 
		    "An earthquake of magnitude " + tags.get("MS_MAG") + " occurred " + F_START_REL_DAYS + " " + DAYS + " ago " + MS_LOC + ". " +
		    "More earthquakes than usual will continue to occur in the area, decreasing " +
		    "in frequency over the following " + F_HORIZON + " or longer. During the next week " +
		    "there are likely to be " + tags.get("NFELT_WK") + " aftershocks large enough to be felt locally, " +
		    "and there is a " + tags.get("PDAMAGE_WK") + "% chance of at least one damaging M" + DAMAGE_MAG + " (or larger) aftershock." +
		    " The earthquake rate may be re-invigorated in response to large aftershocks, should they occur. ";
		
		tags.put("DESCRIPTIVE_TEXT", DESCRIPTIVE_TEXT);
	}
	
	public void setFeltMag(double mag){
		this.feltMag = mag;
	}

	public double getFeltMag(){
		return feltMag;
	}

	public void setDamageMag(int mag){
		this.damageMag = mag;
	}

	public double getDamageMag(){
		return damageMag;
	}

//	This is disabled until the rest of the code can accomodate a flexible magnitude range
//	public void setPredictionMagnitudes(double[] predictionMagnitudes){
//		this.predictionMagnitudes = predictionMagnitudes;
//	}

	//	returns a string token of form "N1 - N2" where N1 and N2 represent a range of values
	private String numberRange(double number1, double number2){
		if (smartRound(number2) == 0)
			return "0*";
		else
			return smartRound(number1) + " - " + smartRound(number2);
		
	}
	
	// rounds the number to an conversational level of accuracy
	private int smartRound(double number){
		long roundNumber;
		
		if (number < 20)
			roundNumber = Math.round(number);
		else if (number < 100)
			roundNumber = Math.round(number/5d)*5;
		else if (number < 200)
			roundNumber = Math.round(number/10d)*10;
		else if (number < 1000)
			roundNumber = Math.round(number/50d)*50;
		else 
			roundNumber = Math.round(number/100d)*100;

		return (int) roundNumber;
	}

	// return the time difference between two GregorianCalendar objects. why is this not a GregCal.method()?
	private static double getDateDelta(GregorianCalendar start, GregorianCalendar end) {
		long diff = end.getTimeInMillis() - start.getTimeInMillis();
		return (double)diff/(double)ProbabilityModelsCalc.MILLISEC_PER_DAY;
	}
	
	public void writeHTMLTable(File outputFile){
		StringBuilder tableString = new StringBuilder();

		//header
		tableString.append(""
				+"<html>\n"
				+"	<head>\n"
				+"		<style>\n"
		        +"			body {font-family:helvetica; background-color: white; page-break-inside:avoid}\n"
		        +"	        .tableForecast {font-size:14px; border:0px solid gray; border-collapse:collapse; margin:0px; text-align:center}\n"
		        +"			.tfElem2 {font-size:14px; border:0px solid gray; border-collapse:collapse; padding:1px; margin:0; text-align:center}\n"
		        +"			.tfElem1 {font-size:14px; background-color:#eeeeee; border-collapse:collapse; padding:1px; margin:0; text-align:center}\n"
		        +" 			.tableFootnote {font-size:12px;text-align:right;}\n"
		        +"		</style>\n"
		        +"	</head>\n"
		        +"	<body>\n");

		//generate forecast table
		tableString.append(""
				+" 	<table class=\"tableForecast\" style=\"width:525px\">\n"
				+"		<tr class=\"forecastHeader\"><th style=\"font-weight:bold; padding:4px\">Forecast Interval</th><th style=\"font-weight:bold\">Magnitude</th><th style=\"font-weight:bold\">Number</th><th style=\"font-weight:bold\">Number Range</th><th style=\"font-weight:bold\">Probability</th></tr>\n"
				+"		<tr><td colspan=\"5\"></td></tr>\n");
			
		String DATE_START = tags.get("F_START_ABS");
		String DATE_END;
		SimpleDateFormat formatter=new SimpleDateFormat("d MMM yyyy, HH:mm:ss");  
		formatter.setTimeZone(TimeZone.getTimeZone("UTC"));
		GregorianCalendar forecastEndDate = new GregorianCalendar();
		
		
		// find the largest magnitude to plot
		// assign variables
		int minMag = 3;
		double maxObsMag = mag0;
		for (int i = 0; i < aftershockModel.magAftershocks.length; i++) {
			if (aftershockModel.magAftershocks[i] > maxObsMag)
				maxObsMag = aftershockModel.magAftershocks[i];
		}
		int maxMag = Math.max(minMag + 3, Math.min(9, (int) Math.ceil(maxObsMag + 0.5)));
		if(D) System.out.println("maxMag: " + maxMag + " largestShockMag: " + maxObsMag);
		

		String[] durString = new String[]{"Day","Week","Month","Year"};
		for (int j = 0; j<4; j++){
			forecastEndDate.setTimeInMillis(forecastStartDate.getTimeInMillis() + (long) (predictionIntervals[j]*ETAS_StatsCalc.MILLISEC_PER_DAY));
			DATE_END = formatter.format(forecastEndDate.getTime());
			tableString.append(""
					+"		<tr class=\"tfElem2\"><td rowspan=\"5\"><div style=\"font-weight:bold;display:inline\">1 " + durString[j]
			 		+"		</div><br>"+ DATE_START +"<br>through<br>"+ DATE_END +"</td>\n");
			 		
			int n = 0;
			for (int i = 0; i < predictionMagnitudes.length; i++) {
				String classStr = (Math.floorMod(n, 2) == 0)?"tfElem1":"tfElem2";
				int mag = (int) predictionMagnitudes[i];
//				if( mag == 3 || i >= predictionMagnitudes.length-4) {//always plot the M3s and then the last four
				if (mag == 3 || (mag > 3 && mag > maxMag - 4 && mag <= maxMag)) {
					tableString.append(""
							+" 		<td class=\""+ classStr +"\">"
							+ "M ≥ " + mag + "</td><td class=\""+ classStr +"\">"
							+ Math.round(number[i][j][0]) + "</td><td class=\""+ classStr +"\">"
							+ numberRange(number[i][j][1], number[i][j][2]) +"</td><td class=\""+ classStr +"\">"
							+ ((probability[i][j] >= 0.001)?probString[i][j]:"<0.1%") +"</td></tr>\n");
					
					n++;
				}
			}
			while(n++ <= 5)
				if (n == 6 && j == 3)
					tableString.append("<tr><td colspan=\"5\" class=\"tableFootnote\">*Earthquake possible but with low probability</td><tr>\n");
				else
					tableString.append("	<tr><td colspan=\"5\"><br></td></tr>\n");
		}
		tableString.append(""
	            +"		</table>\n"
	            +" 	</body>\n"
	            +"</html>\n");
		
		//write to a file
		FileWriter fw;
		try {
			fw = new FileWriter(outputFile, false);
		} catch (IOException e1) {
			//				e1.printStackTrace();
			System.err.println("Couldn't save to file " + outputFile.getAbsolutePath());
			return;
		}

		try {
			fw.append(tableString);
		} catch (IOException e) {
			//					e.printStackTrace();
			System.err.println("Couldn't save to file " + outputFile.getAbsolutePath());
		}

		try {
			fw.close();
		} catch (IOException e) {
			//				e.printStackTrace();
			System.err.println("Problem closing file.");
		}
	}

	
	// Build an html document for displaying the advisory
	public void writeHTML(File outputFile){
	
		StringBuilder outputString = new StringBuilder();
	
		StringBuilder headString = new StringBuilder();
		StringBuilder infoString = new StringBuilder();
		StringBuilder probTableString = new StringBuilder();
		StringBuilder keyString = new StringBuilder();
		StringBuilder imgString = new StringBuilder(); 
		
		//parameters for the svg table
		double barHeight = 40;
		double barWidth = 50;
		double barPadding = 2;
		String[] colors = new String[]{"#abe0ff","#81fff3","#aaff63","#ffd100","#ff5700","#800000","#800000"};
		String[][] tableTags = new String[][]{{"P1_WK","P2_WK","P3_WK","P4_WK","P5_WK","P6_WK","P7_WK"},
			{"P1_MO","P2_MO","P3_MO","P4_MO","P5_MO","P6_MO","P7_MO"},
			{"P1_YR","P2_YR","P3_YR","P4_YR","P5_YR","P6_YR","P7_YR"}};
		
		headString.append(""
				+"    <!DOCTYPE html>\n"
				+"	  <html>\n"
				+"    <head>\n"
				+" 		<meta http-equiv=\"Content-Type\" content=\"text/html;charset=utf-8\">\n"
				+" 		<meta name=\"source\" content=\"USGS AftershockForecaster\">\n"
				+"		<meta name=\"software-version\" content=\"Beta\">\n"
				+"		<meta name=\"software-version-date\" content=\"2018-05-01\">\n"
				+"		<meta name=\"advisory-generated-date\" content=\"" + tags.get("V_DATE") + "\">\n"
				+"		<meta name=\"disclaimer\" content=\"This advisory was generated using softare developed on the OpenSHA platform by the\n"
				+"			US Geological Survey and the US-AID Office of Foreign Disaster Assistance. The information provided in this document\n"
				+"			is not an official statement of the US Geological Survey or any other US Government Entity.\">\n"
				+"		<meta name=\"OpenSHA\" content=\"www.opensha.org/apps\">\n"
				+"        <style>\n"
				+"            body {font-family:helvetica; background-color: white; page-break-inside:avoid}\n"
				+"            h1 {color: black; font-family:helvetica; font-size:20pt; margin:0}\n"
				+"            h2 {color: black; font-family:helvetica; font-size:14pt; margin:0; text-align:center}\n"
				+"            h3 {color: black; font-family:helvetica; font-size:12pt; margin:0;}\n"
				+"            th {font-weight:normal}\n"
				+"			  tr {padding:0px}\n"
				+"            td {padding:0px}\n"
				+"            \n"
				+" 			  .imageButton {color:#000000;background-color:#eeeeee; border:1px solid black; width:150px; padding:2px; font-size:10pt;}\n"
	            +"			  .durButton {color:#000000; background-color:#eeeeee; border:1px solid black; width:100px; padding:2px; font-size:10pt;}\n"
	            +"			  .durButtonDisabled {color:#cccccc; background-color:#eeeeee; border:1px solid black; width:100px; padding:2px; font-size:10pt;}\n"
	            +" 			  .activeImageButton {color:#dd0000; background-color:#ffeeee; border:1px solid black; width:150px; padding:2px; font-size:10pt;}\n"
	            +"			  .activeDurButton {color:#dd0000; background-color:#ffeeee; border:1px solid black; width:100px; padding:2px; font-size:10pt;}\n"
	            +"			  .activeDurButtonDisabled {color:#cccccc; background-color:#eeeeee; border:1px solid black; width:100px; padding:2px; font-size:10pt;}\n"
				+"            .forecast { font-size:14px; border:0px solid gray; border-collapse:collapse; margin:0; text-align:center}\n"
				+"            .forecastHeader { font-size:16px; border:0px solid gray; border-collapse:collapse; text-align:center; vertical-align:center; font-weight:normal; height:20px;}\n"
				+"            .forecastValue { font-size:12px; border:0px solid gray; border-collapse:collapse; text-align:center; margin:0; vertical-align:bottom; color:#666666;}\n"
				+"            .forecastKey { font-size:12px; border:0px solid gray; border-collapse:collapse; text-align:center; margin:0; vertical-align:bottom; padding:0px}\n"
				+"            .forecastBar { height:"+(int)barHeight+"px; width:" + (int)barWidth + "px; padding-top:1px;}\n"
				+"            .forecastRow {border:1px solid #dddddd; border-collapse:collapse; padding-top:1px;}\n"
				+"\n"				
				+"            .forecastBox {stroke-width:0px; x:"+(int)barPadding+"px; width:"+(int)(barWidth-2*barPadding)+"px}\n"
				+" 			  .forecastBoxText {text-anchor:middle; fill:#666666;}\n"
				+" 			  .key {width:30px;height:12px}\n"
				+"            div {min-width:800px; max-width:800px; border-collapse:collapse}\n"
				+"            \n"
				+"        </style>\n"
				+"\n"
				+"		<script>\n"
				+" 			// Script to crop Shakemap constant amount from bottom, accounting for different Shakemap heights\n"
				+" 			function cropShakemap(){\n"
				+ "				var imageData = new Image();\n"
				+" 				imageData.src = document.getElementById('shakemap').src;\n"
		        +" 				var shakemapHeight = imageData.height;\n"
		        +" 				var shakemapWidth = imageData.width;\n"
		        +" 				var shakemapCropHeight = shakemapHeight * (250/shakemapWidth) - 63;\n"
		        +" 				document.getElementById('shakemapCrop').style.height = shakemapCropHeight + \"px\";\n"
		        +" 			}\n"
		        +"\n"
		        +"		</script>\n"
				+"\n"
				+"	    <title>Aftershock Advisory</title>\n"
				+"    </head>\n\n");

		infoString.append("    <body>\n"
				+"        <div>\n"
				+"            <table style=\"text-align:center\">\n"
				+"                <tr>\n"
				+"                    <td style=\"width:150px\">\n"  
				+"						<img id=\"logo\" src=\"Logo.png\" alt=\"\" style=\"max-height:60px;max-width:150px\">\n"
				+"  				  </td>\n"
				+"                    <td style=\"width:500px;height:60px\"><h1 >Aftershock Advisory and Forecast</h1></td>\n"
				+"					  </td>\n"
				+"                    <td style=\"width:150px\">\n"  
				+"						<img id=\"logo\" src=\"Logo.png\" alt=\"\" style=\"max-height:60px;max-width:150px\">\n"
				+"  				  </td>\n"
				+"                </tr>\n"
				+"            </table>\n"
				+"            <br>\n"
				+"            <table>\n"
				+"                <tr>\n"
				+"					  <!-- Mainshock and forecast information -->\n"	
				+"                    <td style=\"width:400px\"><h2 style=\"text-align:left\">" + tags.get("MS_MAG") + " eventID:" + tags.get("MS_NAME") + "</h2></td>\n"
				+"                    <td style=\"width:400px\"><h3 style=\"text-align:right\">Forecast Generated: " + tags.get("V_DATE") + " UTC</h3></td>\n"
				+"                </tr>\n"
				+"            </table>\n"
				+"            <table>\n"
				+"                <tr>\n"
				+"					  <!-- Mainshock origin information -->\n"
				+"                    <td style=\"margin-left:0px;width:300px\"><h3>Origin: " + tags.get("MS_DATETIME") + " UTC</h3></td>\n"
				+"                    <td style=\"text-align:center;margin:auto;width:300px\"><h3>Location: " + tags.get("MS_LAT") + " " + tags.get("MS_LON") + "</h3></td>\n"
				+"                    <td style=\"text-align:right;margin-right:0px;width:200px\"><h3>Depth: " + tags.get("MS_DEPTH") + "</h3></td>\n"
				+"                </tr>\n"
				+"            </table>\n"
				+"        </div>\n"
				+"        <div>\n"
				+"			  <!-- Descriptive forecast text -->\n"
				+"            <p style=\"text-align:justify\">" + tags.get("DESCRIPTIVE_TEXT") + "</p>\n"
				+"        </div>\n\n");

		probTableString.append("        <div style=\"width:800px;text-align:center;\">\n"
				+"            <p><h2>Anticipated aftershock activity</h2>\n"
				+"            <p style=\"margin-top:0px\">Forecast start date: " + tags.get("F_START_ABS") + " UTC</p>\n"
				+"        </div>\n"
				+"\n"
				+"		  <!-- Forecast table with probabilities of different magnitudes -->\n"
				+"        <table style=\"width:650px;margin-left:60px;margin-right:100px\"><!-- Table of Probabilities -->\n"
				+"            <tr>\n"
				+"                <td>\n"
				+"                    <table>\n"
				+"                        <tr>\n"
				+"                            <th></th>\n"
				+"                            <th class=\"forecast\">Probability of at least one aftershock larger than:</th>\n"
				+"                        </tr>\n"
				+"                        <tr>\n"
				+"                            <td>\n"
				+"                                <table>\n"
				+"                                    <tr><th class=\"forecastHeader\" style=\"width:60px\"></th> </tr>\n"
				+"                                    <tr><th class=\"forecastBar\">Week</th></tr>\n"
				+"                                    <tr><th class=\"forecastBar\">Month</th></tr>\n"
				+"                                    <tr><th class=\"forecastBar\">Year</th></tr>\n"
				+"                                </table>\n"
				+"                            </td>\n"
				+"                            <td>\n"
				+"                                <table class=\"forecast\">\n"
				+"                                    <tr>\n"
				+"                                        <th class=\"forecastHeader\">M3</th>\n"
				+"                                        <th class=\"forecastHeader\">M4</th>\n"
				+"                                        <th class=\"forecastHeader\">M5</th>\n"
				+"                                        <th class=\"forecastHeader\">M6</th>\n"
				+"                                        <th class=\"forecastHeader\">M7</th>\n"
				+"                                        <th class=\"forecastHeader\">M8</th>\n"
				+"                                        <th class=\"forecastHeader\">M9</th>\n"
				+"                                    </tr>\n"
				+"                                    \n");

		
		
		//generate forecast table
		for (int j = 0; j<3; j++){
			probTableString.append(""
					+"                                    <tr class = \"forecastRow\">\n");

			for (int i = 0; i<colors.length; i++){
				double probVal = probability[i][j+1];
				double height = barHeight*probVal;
				String probStr = tags.get(tableTags[j][i]); 
				double yVal;
				if (probVal > 0.50) yVal = 11 + barHeight*(1 - probVal);
				else yVal = barHeight*(1 - probVal) - 3;

				probTableString.append(""
						+"                                        <td class=\"forecastValue\">\n"
						+"												<svg class=\"forecastBar\">\n"
						+"													<rect class = \"forecastBox\" y=\"" + (int) (barHeight - (int) height) + "px\" height=\"" + ((int) height) + "px\" width=\"" + barWidth + "px\" fill=\""+colors[i]+"\" />\n"
						+"													<text class = \"forecastBoxText\" x=\"25px\" y=" + String.format("\"%.0fpx\"", yVal) + ">"+probStr+"</text>\n"
						+"												</svg>\n"
						+"                                        </td>\n"
						);
			}
			probTableString.append(""
					+"                                    </tr>\n");
		}
		
		probTableString.append(""
				+"                                </table>\n"
				+"                            </td>\n"
				+"                        </tr>\n"
				+"                    </table>\n"
				+"                </td>\n\n");
		
		keyString.append("               <td style=\"vertical-align:top\">\n"
				+"                    <!-- MMI key-->\n"
				+"                    <table>\n"
				+"                        <tr><td style=\"height:23px\"></td></tr>\n"
				+"                        <tr><td class=\"forecast\">Key to colors</td></tr>\n"
				+"                        <tr><td>\n"
				+"                            <table class=\"forecastKey\" style=\"border:1px solid #dddddd;\">\n"
				+"                                <tr style=\"font-weight:bold\">\n"
				+"                                    <th style=\"width:30px; font-weight:bold\"></th>\n"
				+"                                    <th style=\"width:70px; font-weight:bold\">peak MMI</th>\n"
				+"                                    <th style=\"width:70px; font-weight:bold\">Potential Shaking</th>\n"
				+"                                    <th style=\"width:70px; font-weight:bold\">Potential Damage*</th>\n"
				+"                                </tr>\n"
				+"                                <tr>\n"
				+"                                    <td><svg class=\"key\"><rect class=\"key\" width=\"30px\" height=\"12px\" fill=\"#abe0ff\"/></svg></td>\n"
				+"                                    <td>II-III</td>\n"
				+"                                    <td>weak</td>\n"
				+"                                    <td>none</td>\n"
				+"                                </tr>\n"
				+"                                <tr>\n"
				+"                                    <td><svg class=\"key\"><rect class=\"key\" width=\"30px\" height=\"12px\" fill=\"#81fff3\"/></svg></td>\n"
				+"                                    <td>III-IV</td>\n"
				+"                                    <td>light</td>\n"
				+"                                    <td>very light</td>\n"
				+"                                </tr>\n"
				+"                                <tr>\n"
				+"                                    <td><svg class=\"key\"><rect class=\"key\" width=\"30px\" height=\"12px\" fill=\"#aaff63\"/></svg></td>\n"
				+"                                    <td>IV-V</td>\n"
				+"                                    <td>moderate</td>\n"
				+"                                    <td>light</td>\n"
				+"                                </tr>\n"
				+"                                <tr>\n"
				+"                                    <td><svg class=\"key\"><rect class=\"key\" width=\"30px\" height=\"12px\" fill=\"#ffd100\"/></svg></td>\n"
				+"                                    <td>V-VII</td>\n"
				+"                                    <td>strong</td>\n"
				+"                                    <td>moderate</td>\n"
				+"                                </tr>\n"
				+"                                <tr>\n"
				+"                                    <td><svg class=\"key\"><rect class=\"key\" width=\"30px\" height=\"12px\" fill=\"#ff5700\"/></svg></td>\n"
				+"                                    <td>VIII-IX</td>\n"
				+"                                    <td>severe</td>\n"
				+"                                    <td>heavy</td>\n"
				+"                                </tr>\n"
				+"                                <tr>\n"
				+"                                    <td><svg class=\"key\"><rect class=\"key\" width=\"30px\" height=\"12px\" fill=\"#800000\"/></svg></td>\n"
				+"                                    <td>IX-X</td>\n"
				+"                                    <td>violent</td>\n"
				+"                                    <td>very heavy</td>\n"
				+"                                </tr>\n"
				+"                            </table>\n"
				+"                        </td></tr>\n"
				+"                    </table>\n"
				+"						<p class=\"forecastKey\" style=\"font-size:10px\">*Damage may be higher in vulnerable structures</p>\n"
				+"                </td>\n"
				+"            </tr>\n"
				+"        </table>\n\n");

		// set up javascript to select between different image products. Durations and styles. 
		imgString.append(""
				+" 		<br>\n"
				+"		<div style=\"width:800px;height:500px;margin-left:0px;\">\n"
		);
		
		
		imgString.append(" "
				+"  	<!-- Mainshock shakemap -->\n"
				+"		<table style=\"width:800px;vertical-align:top;text-align:center;\">\n"
				+"	    	<tr>\n"
				+"	        	<td style=\"width:250px;display:inline\">\n"
				+"			    	<br>\n"
				+" 					<p class=\"forecast\" style=\"white-space:pre\">Mainshock ShakeMap\n(previous shaking)</p>\n"
				+" 					<div style=\"height:230px;overflow:hidden;margin: 0 -275px 0px -275px;\" id=\"shakemapCrop\">\n"
				+" 						<!-- Change the link URL to point to your local event summary if preferred. Default is to go to USGS summary -->\n"
				+" 						<a href=\"" + eventURL + "\">\n");
		if (shakeMapURL != null)
			imgString.append(" "
					+"	        			<img style=\"width:250px\" src=\"" + shakeMapURL + "\" alt=\"Mainshock shakemap\" id=\"shakemap\" onload=\"cropShakemap()\">\n"
					+" 						</a>\n");
		else 	
			imgString.append(" "
					+"	        			<p>EventPage"
					+" 						</a>\n");
		imgString.append(" "
				+" 					</div>\n"
				+"		 	   </td>\n"
				+"			    <td style=\"width:550px\" id=\"imageBox\">\n"
				+"	 				<!-- To manually specifiy the image you want to see displayed, replace src=\"...\" with the desired filename. -->\n"
				+"	          	<img style=\"margin:auto;width:550px;max-height:480px\" src=\"ratemap.png\" alt=\"Graphical Forecast\" id=\"theimage\">\n"
				+"	      	</td>\n"
				+"		    </tr>\n"
				+"		</table>\n"
				);


		imgString.append(""
			
		);
		
		imgString.append(""
				+"		</div>\n"
				+" 		<div class=\"forecast\">\n"
				+"			<!-- This would be a good place for a link to an online source of the forecast, if applicable -->\n"
				+"			This forecast will be updated as new information becomes available.\n"
				+"		</div>\n"
				+" 		<br>\n"
				+"\n"
				+"  <div style=\"margin-left:0px;width:800px;text-align:center\">\n"
				+"		<!-- Set up buttons for changing which image type is displayed -->\n"
				+" 		<input type=\"button\" class=\"imageButton\" value=\"Table only\" id=\"imageButton0\" onClick=\"showTable();\">\n"
				+"		<input type=\"button\" class=\"imageButton\" value=\"Magnitude Distribution\" id=\"imageButton1\" onClick=\"changeImage('1');\">\n"
				+"		<input type=\"button\" class=\"imageButton\" value=\"Number with time\" id=\"imageButton2\" onClick=\"changeImage('2');\">\n"
				+"		<input type=\"button\" class=\"activeImageButton\" value=\"Rate map\" id=\"imageButton3\" onClick=\"changeImage('3');\">\n"
				+"		<input type=\"button\" class=\"imageButton\" value=\"Shaking map\" id=\"imageButton4\" onClick=\"changeImage('4');\">\n"
				+"  </div>\n"
				+"\n"
				+"  <div style=\"margin-left:0px;width:800px;text-align:center\">\n"
				+"  	<!-- Set up buttons for changing the image duration -->\n"
				+"		<input type=\"button\" class=\"durButtonDisabled\" value=\"Day\" id=\"durationButton1\" onClick=\"changeDuration('1');\">\n"
				+"		<input type=\"button\" class=\"durButtonDisabled\" value=\"Week\" id=\"durationButton2\" onClick=\"changeDuration('2');\">\n"
				+"		<input type=\"button\" class=\"activeDurButtonDisabled\" value=\"Month\" id=\"durationButton3\" onClick=\"changeDuration('3');\">\n"
				+"		<input type=\"button\" class=\"durButtonDisabled\" value=\"Year\" id=\"durationButton4\" onClick=\"changeDuration('4');\">\n"
				+"  </div>\n"
				+"\n"
				+"	<script>\n"
				+"		function showTable(){\n"
				+"			document.getElementById('imageBox').innerHTML = '<iframe style=\"width:550px;height:480px;border:none\" src=\"Table.html\"></iframe>';\n"
				+"			for (i = 1; i < 5; i++){\n"
				+"				document.getElementById('imageButton' + i).className = 'imageButton';\n"
				+"			}\n"
				+"			document.getElementById('imageButton0').className = 'activeImageButton';\n"
				+"\n"
				+"			for (i = 1; i < 5; i++){\n"
                +"			    if (document.getElementById('durationButton' + i).className == 'durButton')\n"
                +"			        document.getElementById('durationButton' + i).className = 'durButtonDisabled';\n"
                +"			    if (document.getElementById('durationButton' + i).className == 'activeDurButton')\n"
                +"			        document.getElementById('durationButton' + i).className = 'activeDurButtonDisabled';\n"
                +"				}\n"
				+"		}\n"
				+"\n"
				+" 		// Script for changing the image in response to button click\n"
				+"		function changeImage(buttonNumber){\n"
				+" 		// first make sure we've got an image\n"
		        +"    	document.getElementById('imageBox').innerHTML = '<img style=\"margin:auto;width:550px;max-height:480px\" src=\"ratemap.png\" alt=\"Graphical Forecast\" id=\"theimage\">';\n"		            
		        +"\n"
		        +" 		// now decide which image\n"
				+"			durationIndex = 1;\n"
				+"			for (i = 1; i < 5; i++){\n"
				+"				if (document.getElementById('durationButton' + i).className.includes('active')){\n"
				+"	 				var durationIndex = i;\n"
				+"		   		}\n"
				+"			}\n"
				+"			console.log(durationIndex);\n"
				+"\n"
				+"			var dur;\n"
				+"			switch (durationIndex){\n"
				+"				case 1:\n"
				+"			    	dur = 'day';\n"
				+"					break;\n"
				+"				case 2:\n"
				+"					dur = 'week';\n"
				+"					break;\n"
				+"				case 3:\n"
				+"					dur = 'month';\n"
				+"					break;\n"
				+"				case 4:\n"
				+"					dur = 'year';\n"
				+"					break;\n"
				+"				default:\n"
				+"					dur = '';\n"
				+"					break;\n"
				+"			}\n"
				+"\n"
				+"			var imgName = document.getElementById('theimage').src;\n"
				+"\n"
				+"			switch (buttonNumber){\n"
				+"				case '1':\n"
				+"					var imgName = 'number';\n"
				+"				    var durNeeded = true;\n"
				+"					break;\n"
				+"				case '2':\n"
				+"					var imgName = 'forecastCmlNum';\n"
				+"		            var durNeeded = false;\n"
				+"	 				break;\n"
				+"				case '3':\n"
				+"					var imgName = 'ratemap';\n"
				+"	        	    var durNeeded = false;\n"
				+"					break;\n"
				+"				case '4':\n"
				+"					var imgName = 'shaking';\n"
				+"	        	    var durNeeded = true;\n"
				+"					break;\n"
				+"			}\n"
				+"\n"
				+"			if (durNeeded) imgName = imgName + dur;\n"
				+"\n"
				+"			var image = document.getElementById('theimage');\n"
				+"			image.src = imgName + '.png';\n"
				+"\n"
				+"			for (i = 0; i < 5; i++){\n"
				+"				document.getElementById('imageButton' + i).className = 'imageButton';\n"
				+"			}\n"
				+"			document.getElementById('imageButton' + buttonNumber).className = 'activeImageButton';\n"
				+"\n"
				+"			if (durNeeded){\n"
				+"				for (i = 1; i < 5; i++){\n"
                +"			    	if (document.getElementById('durationButton' + i).className == 'durButtonDisabled')\n"
                +"			    	    document.getElementById('durationButton' + i).className = 'durButton';\n"
                +"			    	if (document.getElementById('durationButton' + i).className == 'activeDurButtonDisabled')\n"
                +"			    	    document.getElementById('durationButton' + i).className = 'activeDurButton';\n"
                +"				}\n"
                +"			} else {\n"
                +"				for (i = 1; i < 5; i++){\n"
                +"				    if (document.getElementById('durationButton' + i).className == 'durButton')\n"
                +"				        document.getElementById('durationButton' + i).className = 'durButtonDisabled';\n"
                +"				    if (document.getElementById('durationButton' + i).className == 'activeDurButton')\n"
                +"				        document.getElementById('durationButton' + i).className = 'activeDurButtonDisabled';\n"
                +"				}\n"
                +"		 	}\n"
				+"		}\n"
				+"\n"
				+"		function changeDuration(buttonNumber){\n"
				+"			var imgName = document.getElementById('theimage').src;\n"
				+"\n"
				+"	        if (imgName.includes('number')) {\n"
				+"	    		var imgtype = 'number';\n"
				+"	          	durNeeded = true;\n"
				+"			} else if (imgName.includes('forecastCmlNum')){\n"
				+"				var imgtype = 'forecastCmlNum';\n"
				+"				durNeeded = false;\n"
				+"			} else if (imgName.includes('rate')){\n"
				+"				var imgtype = 'ratemap';\n"
				+"				durNeeded = false;\n"
				+"			} else if (imgName.includes('shaking')){\n"
				+"				var imgtype = 'shaking';\n"
				+"				durNeeded = true;\n"
				+"			} else {\n"
				+"				var imgtype = 'ImageNotRecognized';\n"
				+"				durNeeded = false;\n"
				+"			}\n"
				+"\n"
				+"			if (durNeeded) {\n"
				+"				switch (buttonNumber){\n"
				+"				    case '1':\n"
				+"	    				imgName = imgtype+'day';\n"
				+"	 					break;\n"
				+"					case '2':\n"
				+"						imgName = imgtype+'week';\n"
				+"						break;\n"
				+"					case '3':\n"
				+"						imgName = imgtype+'month';\n"
				+"						break;\n"
				+"					case '4':\n"
				+"						imgName = imgtype+'year';\n"
				+"						break;\n"
				+"			    }\n"
				+"			} else {\n"
				+"				imgName = imgtype;\n"
				+"			}\n"
				+"\n"
				+"			var image = document.getElementById('theimage');\n"
				+"			image.src = imgName + '.png';\n"
				+"\n"
				+"			if (durNeeded){\n"
				+"        		for (i = 1; i < 5; i++){\n"
		        +"            		document.getElementById('durationButton' + i).className = 'durButton';\n"
		        +"        		}\n"
		        +"        		document.getElementById('durationButton' + buttonNumber).className = 'activeDurButton';\n"
		        +"    		} else {\n"
		        +"       		for (i = 1; i < 5; i++){\n"
		        +"            		document.getElementById('durationButton' + i).className = 'durButtonDisabled';\n"
		        +"        		}\n"
		        +"        		document.getElementById('durationButton' + buttonNumber).className = 'activeDurButtonDisabled';\n"
		        +"    		}\n"
				+"		}\n"
				+"	</script>\n"
				);

		imgString.append(""
				+" 	</body>\n"
				+"</html>\n");

		outputString.append(headString);
		outputString.append(infoString);
		outputString.append(probTableString);
		outputString.append(keyString);
		outputString.append(imgString);
		
		// write file
		FileWriter fw;
		try {
			fw = new FileWriter(outputFile, false);
		} catch (IOException e1) {
			//				e1.printStackTrace();
			System.err.println("Couldn't save to file " + outputFile.getAbsolutePath());
			return;
		}

		try {
			fw.append(outputString);
		} catch (IOException e) {
			//					e.printStackTrace();
			System.err.println("Couldn't save to file " + outputFile.getAbsolutePath());
		}
		
		try {
			fw.close();
		} catch (IOException e) {
			//				e.printStackTrace();
			System.err.println("Problem closing file.");
		}
	}

	public void setShakeMapURL(String shakeMapURL) {
		this.shakeMapURL = shakeMapURL;
	}

	private static double unwrap(double lon){
		if (lon > 180)
			lon -= 360;
		return lon;
	}
	
	public static void main(String[] args){
		GraphicalForecast gf = new GraphicalForecast();
		gf.processForecastStrings();
		gf.assignForecastStrings();
		gf.setShakeMapURL("https://earthquake.usgs.gov/archive/product/shakemap/atlas20100404224043/atlas/1520888708106/download/intensity.jpg");
		try{
			gf.writeHTML(new File(System.getenv("HOME") + "/example_forecast.html"));
			gf.writeHTMLTable(new File(System.getenv("HOME") + "/Table.html"));
		}catch(Exception e){
			e.printStackTrace();
		}
	}

//	private static String logoUSGSxml(){
//		return new String(""
//				+ "<svg xmlns=\"http://www.w3.org/2000/svg\" xml:space=\"preserve\" height=\"60\" width=\"150\" version=\"1.1\" overflow=\"visible\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" viewBox=\"0 0 500 200\" enable-background=\"new 0 0 500 200\">\n"
//				+ "<path id=\"USGS\" d=\"m234.95 15.44v85.037c0 17.938-10.132 36.871-40.691 36.871-27.569 0-40.859-14.281-40.859-36.871v-85.04h25.08v83.377c0 14.783 6.311 20.593 15.447 20.593 10.959 0 15.943-7.307 15.943-20.593v-83.377h25.08m40.79 121.91c-31.058 0-36.871-18.27-35.542-39.03h25.078c0 11.462 0.5 21.092 14.283 21.092 8.471 0 12.619-5.482 12.619-13.618 0-21.592-50.486-22.922-50.486-58.631 0-18.769 8.968-33.715 39.525-33.715 24.42 0 36.543 10.963 34.883 36.043h-24.419c0-8.974-1.491-18.106-11.627-18.106-8.136 0-12.952 4.486-12.952 12.787 0 22.757 50.492 20.763 50.492 58.465 0 31.06-22.75 34.72-41.85 34.72m168.6 0c-31.06 0-36.871-18.27-35.539-39.03h25.075c0 11.462 0.502 21.092 14.285 21.092 8.475 0 12.625-5.482 12.625-13.618 0-21.592-50.494-22.922-50.494-58.631 0-18.769 8.969-33.715 39.531-33.715 24.412 0 36.536 10.963 34.875 36.043h-24.412c0-8.974-1.494-18.106-11.625-18.106-8.145 0-12.955 4.486-12.955 12.787 0 22.757 50.486 20.763 50.486 58.465 0 31.06-22.75 34.72-41.85 34.72m-79.89-46.684h14.76v26.461l-1.229 0.454c-3.815 1.332-8.301 2.327-12.452 2.327-14.287 0-17.943-6.645-17.943-44.177 0-23.256 0-44.348 15.615-44.348 12.146 0 14.711 8.198 14.933 18.107h24.981c0.197-23.271-14.789-36.043-38.42-36.043-41.021 0-42.521 30.724-42.521 60.954 0 45.507 4.938 63.167 47.12 63.167 9.784 0 25.36-2.211 32.554-4.18 0.437-0.115 1.212-0.597 1.212-1.217v-59.598h-38.611v18.09\" fill=\"#007150\"/>\n"
//				+ "<path id=\"waves\" d=\"m48.736 55.595l0.419 0.403c11.752 9.844 24.431 8.886 34.092 2.464 6.088-4.049 33.633-22.367 49.202-32.718v-10.344h-116.03v27.309c7.071-1.224 18.47-0.022 32.316 12.886m43.651 45.425l-13.705-13.142c-1.926-1.753-3.571-3.04-3.927-3.313-11.204-7.867-21.646-5.476-26.149-3.802-1.362 0.544-2.665 1.287-3.586 1.869l-28.602 19.13v34.666h116.03v-24.95c-2.55 1.62-18.27 10.12-40.063-10.46m-44.677-42.322c-0.619-0.578-1.304-1.194-1.915-1.698-13.702-10.6-26.646-5.409-29.376-4.116v11.931l6.714-4.523s10.346-7.674 26.446 0.195l-1.869-1.789m16.028 15.409c-0.603-0.534-1.214-1.083-1.823-1.664-12.157-10.285-23.908-7.67-28.781-5.864-1.382 0.554-2.7 1.303-3.629 1.887l-13.086 8.754v12.288l21.888-14.748s10.228-7.589 26.166 0.054l-0.735-0.707m68.722 12.865c-4.563 3.078-9.203 6.203-11.048 7.441-4.128 2.765-13.678 9.614-29.577 2.015l1.869 1.797c0.699 0.63 1.554 1.362 2.481 2.077 11.418 8.53 23.62 7.304 32.769 1.243 1.267-0.838 2.424-1.609 3.507-2.334v-12.234m0-24.61c-10.02 6.738-23.546 15.833-26.085 17.536-4.127 2.765-13.82 9.708-29.379 2.273l1.804 1.729c0.205 0.19 0.409 0.375 0.612 0.571l-0.01 0.01 0.01-0.01c12.079 10.22 25.379 8.657 34.501 2.563 5.146-3.436 12.461-8.38 18.548-12.507l-0.01-12.165m0-24.481c-14.452 9.682-38.162 25.568-41.031 27.493-4.162 2.789-13.974 9.836-29.335 2.5l1.864 1.796c1.111 1.004 2.605 2.259 4.192 3.295 10.632 6.792 21.759 5.591 30.817-0.455 6.512-4.351 22.528-14.998 33.493-22.285v-12.344\" fill=\"#007150\"/>\n"
//				+ "<path id=\"tagline\" d=\"m22.329 172.13c-0.247 0.962-0.401 1.888-0.251 2.554 0.195 0.68 0.749 1.011 1.923 1.011 1.171 0 2.341-0.757 2.642-2.183 0.954-4.479-9.653-3.479-8.218-10.224 0.972-4.567 5.792-5.954 9.607-5.954 4.022 0 7.257 1.928 5.951 6.495h-5.783c0.312-1.466 0.33-2.347-0.007-2.722-0.298-0.379-0.783-0.463-1.413-0.463-1.297 0-2.188 0.841-2.492 2.264-0.714 3.354 9.718 3.189 8.271 9.975-0.781 3.688-4.388 6.457-9.29 6.457-5.157 0-8.316-1.306-6.724-7.21h5.784m25.284-6.85c0.667-3.141 0.093-4.188-1.75-4.188-2.513 0-3.193 2.22-4.13 6.619-1.373 6.455-1.124 7.838 1.057 7.838 1.844 0 3.08-1.676 3.667-4.439h5.909c-1.218 5.741-4.847 8.215-10.382 8.215-7.627 0-7.645-4.654-6.234-11.273 1.229-5.784 3.119-10.729 10.915-10.729 5.447 0 8.033 2.433 6.856 7.964h-5.908m18.389-15.69l-0.989 4.651h-5.909l0.989-4.651h5.909m-6.233 29.31h-5.909l4.501-21.165h5.909l-4.501 21.16zm282.77-29.31l-0.991 4.651h-5.911l0.99-4.651h5.91m-6.23 29.31h-5.906l4.496-21.165h5.911l-4.5 21.16zm-259.03-12.95c0.438-2.052 1.144-4.984-1.664-4.984-2.727 0-3.36 3.186-3.743 4.984h5.407m-6.111 3.31c-0.533 2.516-1.251 6.284 1.345 6.284 2.097 0 2.945-2.012 3.318-3.771h5.992c-0.574 2.306-1.728 4.192-3.429 5.489-1.66 1.298-3.916 2.055-6.681 2.055-7.63 0-7.645-4.654-6.239-11.273 1.229-5.784 3.12-10.729 10.915-10.729 7.965 0 7.75 5.152 6.097 11.944h-11.318zm22.462-9.38h0.083c1.575-1.886 3.31-2.557 5.534-2.557 2.808 0 4.923 1.676 4.3 4.607l-3.608 16.978h-5.909l3.099-14.584c0.401-1.887 0.38-3.353-1.507-3.353-1.886 0-2.536 1.468-2.936 3.353l-3.098 14.584h-5.909l4.5-21.165h5.905l-0.452 2.14m23.465 5.4c0.667-3.141 0.093-4.188-1.751-4.188-2.512 0-3.194 2.22-4.131 6.619-1.373 6.455-1.122 7.838 1.058 7.838 1.843 0 3.079-1.676 3.668-4.439h5.909c-1.222 5.741-4.846 8.215-10.382 8.215-7.627 0-7.644-4.654-6.235-11.273 1.229-5.784 3.116-10.729 10.912-10.729 5.45 0 8.037 2.433 6.86 7.964h-5.908m19.61 0.67c0.434-2.052 1.145-4.984-1.664-4.984-2.725 0-3.36 3.186-3.743 4.984h5.4m-6.11 3.31c-0.54 2.516-1.255 6.284 1.344 6.284 2.095 0 2.94-2.012 3.316-3.771h5.992c-0.574 2.306-1.728 4.192-3.432 5.489-1.656 1.298-3.912 2.055-6.68 2.055-7.627 0-7.647-4.654-6.237-11.273 1.231-5.784 3.12-10.729 10.915-10.729 7.961 0 7.747 5.152 6.093 11.944h-11.31zm36.12-15.9c-2.352-0.168-3.051 0.758-3.507 2.896l-0.42 1.481h2.77l-0.775 3.646h-2.768l-3.723 17.521h-5.909l3.722-17.521h-2.638l0.774-3.646h2.682c1.188-5.292 2.251-8.231 8.516-8.231 0.713 0 1.376 0.041 2.08 0.082l-0.79 3.77m9.56 14.35c0.937-4.399 1.198-6.619-1.317-6.619-2.512 0-3.196 2.22-4.13 6.619-1.373 6.455-1.122 7.838 1.057 7.838 2.17 0 3.01-1.38 4.39-7.84m-11.43 0.34c1.229-5.784 3.117-10.729 10.912-10.729s7.586 4.945 6.355 10.729c-1.409 6.619-3.403 11.274-11.032 11.274-7.63 0-7.65-4.65-6.24-11.27zm27.86-10.31l-0.577 2.723h0.082c1.607-2.431 3.77-3.143 6.162-3.143l-1.122 5.279c-5.129-0.335-5.854 2.682-6.298 4.779l-2.448 11.524h-5.909l4.496-21.165h5.62m17.83 14.58c-0.32 1.51-0.464 3.354 1.465 3.354 3.479 0 3.935-4.693 4.421-7-2.95 0.13-5.08-0.12-5.88 3.65m10.34 2.64c-0.281 1.305-0.395 2.645-0.546 3.941h-5.491l0.345-2.808h-0.082c-1.721 2.18-3.664 3.229-6.223 3.229-4.105 0-4.961-3.06-4.181-6.748 1.488-7 6.958-7.295 12.43-7.206l0.347-1.638c0.385-1.804 0.41-3.104-1.729-3.104-2.054 0-2.549 1.551-2.908 3.229h-5.784c0.545-2.558 1.69-4.19 3.278-5.151 1.553-1.011 3.561-1.388 5.827-1.388 7.5 0 7.777 3.229 6.959 7.084l-2.25 10.57zm23.39-9.68c0.667-3.141 0.093-4.188-1.749-4.188-2.515 0-3.196 2.22-4.132 6.619-1.373 6.455-1.122 7.838 1.059 7.838 1.842 0 3.08-1.676 3.668-4.439h5.909c-1.221 5.741-4.848 8.215-10.382 8.215-7.627 0-7.642-4.654-6.237-11.273 1.232-5.784 3.121-10.729 10.916-10.729 5.447 0 8.034 2.433 6.857 7.964h-5.909m32.45 7.04c-0.322 1.51-0.463 3.354 1.465 3.354 3.479 0 3.936-4.693 4.42-7-2.96 0.13-5.09-0.12-5.89 3.65m10.34 2.64c-0.281 1.305-0.396 2.645-0.547 3.941h-5.49l0.344-2.808h-0.081c-1.72 2.18-3.659 3.229-6.226 3.229-4.104 0-4.961-3.06-4.18-6.748 1.487-7 6.957-7.295 12.432-7.206l0.351-1.638c0.378-1.804 0.405-3.104-1.733-3.104-2.057 0-2.549 1.551-2.909 3.229h-5.781c0.543-2.558 1.69-4.19 3.276-5.151 1.561-1.011 3.563-1.388 5.828-1.388 7.5 0 7.776 3.229 6.957 7.084l-2.24 10.57zm12.97-15.08h0.08c1.58-1.886 3.311-2.557 5.533-2.557 2.805 0 4.924 1.676 4.297 4.607l-3.608 16.978h-5.905l3.101-14.584c0.397-1.887 0.377-3.353-1.507-3.353-1.889 0-2.534 1.468-2.936 3.353l-3.104 14.584h-5.912l4.505-21.165h5.911l-0.45 2.14m18.63 15.08c2.141 0 2.903-2.22 3.854-6.703 0.986-4.652 1.342-7.291-0.843-7.291-2.22 0-2.926 1.549-4.3 8.004-0.42 1.97-1.56 5.99 1.29 5.99m12-17.22l-4.686 22.045c-0.313 1.465-1.461 7.25-9.59 7.25-4.398 0-7.937-1.129-6.965-6.285h5.785c-0.187 0.883-0.222 1.637 0.048 2.137 0.265 0.547 0.87 0.841 1.794 0.841 1.469 0 2.473-1.388 2.93-3.521l0.862-4.063h-0.084c-1.229 1.632-3.082 2.47-5.008 2.47-6.499 0-4.946-5.949-3.928-10.729 0.992-4.65 2.33-10.563 8.49-10.563 2.096 0 3.704 0.92 4.12 2.893h0.085l0.525-2.473h5.616v-0.002h0.03zm19.78 2.14h0.092c1.572-1.886 3.306-2.557 5.525-2.557 2.809 0 4.924 1.676 4.301 4.607l-3.606 16.978h-5.912l3.104-14.584c0.397-1.887 0.377-3.353-1.51-3.353-1.889 0-2.531 1.468-2.932 3.353l-3.104 14.584h-5.91l4.5-21.165h5.91l-0.45 2.14m18.63 15.08c2.137 0 2.901-2.22 3.854-6.703 0.992-4.652 1.344-7.291-0.836-7.291-2.222 0-2.928 1.549-4.301 8.004-0.41 1.97-1.56 5.99 1.29 5.99m12-17.22l-4.687 22.045c-0.313 1.465-1.454 7.25-9.593 7.25-4.398 0-7.931-1.129-6.957-6.285h5.785c-0.191 0.883-0.226 1.637 0.043 2.137 0.264 0.547 0.874 0.841 1.791 0.841 1.471 0 2.474-1.388 2.929-3.521l0.863-4.063h-0.079c-1.23 1.632-3.086 2.47-5.011 2.47-6.497 0-4.939-5.949-3.928-10.729 0.993-4.65 2.33-10.563 8.496-10.563 2.094 0 3.695 0.92 4.117 2.893h0.085l0.525-2.473h5.615v-0.002h0.02zm9.23 0h5.869l-0.482 15.926h0.079l7.044-15.926h6.278l0.08 15.926h0.083l6.438-15.926h5.666l-10.043 21.165h-6.195l-0.608-14.039h-0.081l-7.083 14.039h-6.281l-0.78-21.16m41.29 9.97c0.937-4.399 1.196-6.619-1.313-6.619-2.521 0-3.203 2.22-4.133 6.619-1.373 6.455-1.121 7.838 1.059 7.838 2.17 0 3.01-1.38 4.38-7.84m-11.43 0.34c1.226-5.784 3.114-10.729 10.911-10.729 7.796 0 7.586 4.945 6.354 10.729-1.409 6.619-3.401 11.274-11.028 11.274s-7.64-4.65-6.23-11.27zm28.19-10.31l-0.582 2.723h0.086c1.608-2.431 3.771-3.143 6.16-3.143l-1.123 5.279c-5.125-0.335-5.851 2.682-6.297 4.779l-2.449 11.524h-5.906l4.496-21.165h5.61m-182.4-0.42c-2.219 0-3.955 0.671-5.53 2.557h-0.086l2.188-10.285h-5.91l-6.231 29.314h5.91l3.103-14.585c0.4-1.885 1.047-3.353 2.935-3.353s1.906 1.468 1.51 3.353l-3.102 14.585h5.907l3.605-16.976c0.62-2.93-1.49-4.6-4.3-4.6m192.26-7.73l-6.231 29.313h5.912l6.23-29.313h-5.91m15.8 18.54c-1.133 5.324-1.979 7.545-3.992 7.545-2.135 0-2.043-2.221-0.912-7.545 0.9-4.231 1.48-7.166 4.039-7.166 2.43 0 1.76 2.93 0.86 7.16m3.95-18.54l-2.158 10.16h-0.084c-0.834-1.804-2.165-2.433-4.301-2.433-5.953 0-7.187 6.582-8.097 10.856-0.924 4.355-2.578 11.146 3.541 11.146 2.268 0 4.051-0.715 5.574-2.768h0.087l-0.504 2.349h5.62l6.229-29.315h-5.909v-0.01z\" fill=\"#007150\"/>\n"
//				+ "</svg>\n"
//				);
//	}
	
	
}
