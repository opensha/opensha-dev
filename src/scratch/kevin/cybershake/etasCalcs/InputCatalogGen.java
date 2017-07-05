package scratch.kevin.cybershake.etasCalcs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import scratch.UCERF3.erf.utils.ProbabilityModelsCalc;

import com.google.common.collect.Lists;

public class InputCatalogGen {

	public static void main(String[] args) throws IOException {
		// this will parse the input catalog, and change the times so that the last one occurs right at 2014
		
		List<String> inputCat = Lists.newArrayList();
		inputCat.add("14432456	3.1	3/21/2009	13:12:18	33	18.8	N	115	43.8	W	4.7	4");
		inputCat.add("14432496	3.3	3/21/2009	13:17:03	33	18.9	N	115	44	W	4.7	4");
		inputCat.add("14432576	3.1	3/21/2009	13:40:06	33	18.9	N	115	44	W	4.9	4");
		inputCat.add("14433456	4.8	3/24/2009	4:55:43	33	19.1	N	115	43.7	W	5.8	4");
		inputCat.add("14433696	3.1	3/24/2009	6:52:51	33	18.8	N	115	44.2	W	4.2	4");
		inputCat.add("14434264	3.1	3/24/2009	16:43:15	33	18.4	N	115	44.6	W	5	5");
		inputCat.add("14434504	3	3/24/2009	22:50:55	33	17.4	N	115	43.5	W	6	7");
		inputCat.add("14434688	3.6	3/25/2009	0:51:23	33	17.4	N	115	43.3	W	4.4	7");
		inputCat.add("14435296	3.7	3/25/2009	12:59:44	33	17.5	N	115	43.2	W	6.5	7");
		inputCat.add("14435776	4	3/25/2009	20:25:21	33	17.6	N	115	43.3	W	7.4	6");
		inputCat.add("10111842	3.4	3/25/2009	20:25:50	33	18.1	N	115	42.4	W	5	6");
		inputCat.add("14436024	3	3/25/2009	22:12:33	33	17.5	N	115	43.1	W	11	7");
		
		long targetOT = Math.round((2014.0-1970.0)*ProbabilityModelsCalc.MILLISEC_PER_YEAR)
				+(8*ProbabilityModelsCalc.MILLISEC_PER_DAY)/24; // occurs at 2014;
		
		List<ObsEqkRupture> origRups = Lists.newArrayList();
		
		long maxTime = 0;
		
		SimpleDateFormat df = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss z");
		int id = 0;
		for (String line : inputCat) {
			String[] split = line.trim().split("\t");
			double mag = Double.parseDouble(split[1]);
			String dateStr = split[2]+" "+split[3]+" PDT";
			Date date = null;
			try {
				date = df.parse(dateStr);
			} catch (ParseException e) {
				e.printStackTrace();
				System.exit(0);
			}
			
			long millis = date.getTime();
			if (millis > maxTime)
				maxTime = millis;
			
			double lat = Double.parseDouble(split[4]) + Double.parseDouble(split[5])/60d;
			double lon = -(Double.parseDouble(split[7]) + Double.parseDouble(split[8])/60d);
			double dep = Double.parseDouble(split[10]);
			
			ObsEqkRupture rup = new ObsEqkRupture((id++)+"", millis, new Location(lat, lon, dep), mag);
			origRups.add(rup);
		}
		
		long timeDelta = targetOT - maxTime;
		
		FileWriter fw = new FileWriter(new File("/home/kevin/OpenSHA/UCERF3/cybershake_etas/bombay_catalog.txt"));
		DateFormat outDF = new SimpleDateFormat("yyyy MM dd HH mm ss.SSS");
		for (ObsEqkRupture rup : origRups) {
			long origTime = rup.getOriginTime();
			long modTime = origTime + timeDelta;
			
			// FORMAT:
			// 0-Year 1-month 2-day 3-hour 4-minute 5-second 6-lat 7-long 8-depth 9-mag
			// 10-magType 11-magSource 12-magErr 13-magRounding 14-EarthquakeID
			
			Date date = new Date(modTime);
			Location loc = rup.getHypocenterLocation();
			
			String line = outDF.format(date)+" "+loc.getLatitude()+" "+loc.getLongitude()+" "+loc.getDepth()+" "+rup.getMag()
					+" 0.00 0.00 NaN NaN "+rup.getEventId();
			fw.write(line+"\n");
		}
		fw.close();
	}

}
