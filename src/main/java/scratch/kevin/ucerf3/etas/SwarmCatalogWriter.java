package scratch.kevin.ucerf3.etas;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collections;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.TimeZone;

import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.event.EventQuery;
import gov.usgs.earthquake.event.EventWebService;
import gov.usgs.earthquake.event.Format;
import gov.usgs.earthquake.event.JsonEvent;

public class SwarmCatalogWriter {
	
	private static final boolean D = false;
	
	private static ObsEqkRupList fetchEvents(Region region, long startTime, long endTime, double minMag)
			throws MalformedURLException {
		EventWebService service = new EventWebService(new URL("https://earthquake.usgs.gov/fdsnws/event/1/"));
		
		EventQuery query = new EventQuery();
		
//		Preconditions.checkState(minDepth < maxDepth, "Min depth must be less than max depth");
//		query.setMinDepth(new BigDecimal(minDepth));
//		query.setMaxDepth(new BigDecimal(maxDepth));
		
//		Preconditions.checkState(minDays < maxDays, "Min days must be less than max days");
		// time zones shouldn't be an issue since we're just adding to the original catalog time, whatever
		// time zone that is in.
//		long startTime = eventTime + (long)(minDays*day_millis);
//		long endTime = eventTime + (long)(maxDays*day_millis);
		query.setStartTime(new Date(startTime));
		if (endTime > Long.MIN_VALUE)
			query.setEndTime(new Date(endTime));
		
		query.setMinLatitude(new BigDecimal(region.getMinLat()));
		query.setMaxLatitude(new BigDecimal(region.getMaxLat()));
		query.setMinLongitude(new BigDecimal(region.getMinLon()));
		query.setMaxLongitude(new BigDecimal(region.getMaxLon()));
		
		if (minMag > 0d)
			query.setMinMagnitude(new BigDecimal(minMag));
		
		if (D)
			try {
				System.out.println(service.getUrl(query, Format.GEOJSON));
			} catch (MalformedURLException e) {
				e.printStackTrace();
			}
		List<JsonEvent> events;
		try {
			events = service.getEvents(query);
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}
		
		ObsEqkRupList rups = new ObsEqkRupList();
		for (JsonEvent event : events) {
			ObsEqkRupture rup = eventToObsRup(event);
			rups.add(rup);
		}
		
		if (!region.isRectangular()) {
			if (D) System.out.println("Fetched "+rups.size()+" events before region filtering");
			for (int i=rups.size(); --i>=0;)
				if (!region.contains(rups.get(i).getHypocenterLocation()))
					rups.remove(i);
		}
		
		if (D) System.out.println("Returning "+rups.size()+" aftershocks");
		
		return rups;
	}
	
	public static ObsEqkRupture eventToObsRup(JsonEvent event) {
		double lat = event.getLatitude().doubleValue();
		double lon = event.getLongitude().doubleValue();
		double dep = event.getDepth().doubleValue();
		Location hypo = new Location(lat, lon, dep);
		double mag = event.getMag().doubleValue();
		
		ObsEqkRupture rup = new ObsEqkRupture(event.getEventId().toString(),
				event.getTime().getTime(), hypo, mag);
		
		return rup;
	}

	public static void main(String[] args) throws IOException {
//		Region region = new Region(new Location(33.298, -115.713), 30d);
//		long startTime = 1474889900000l;
//		long endTime = 1474990200000l;
////		long endTime = Long.MIN_VALUE;
//		File outputFile = new File("/tmp/2016_bombay_swarm.txt");
		
//		Region region = new Region(new Location(33.298, -115.713), 30d);
//		long startTime = 1237593600000l;
//		long endTime = 1238630400000l;
////		long endTime = Long.MIN_VALUE;
//		File outputFile = new File("/tmp/2009_bombay_swarm.txt");
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		long startTime = 	1466665200000l;
		long endTime = 		1498201200000l;
//		long endTime = Long.MIN_VALUE;
		File outputFile = new File("/tmp/csep_bench_inputs.txt");
		
		ObsEqkRupList events = fetchEvents(region, startTime, endTime, 2.5d);
		events.sortByOriginTime();
		
		System.out.println(events.size()+" events");
		FileWriter fw = new FileWriter(outputFile);
		for (ObsEqkRupture rup : events) {
			// FORMAT:
			// 0-Year 1-month 2-day 3-hour 4-minute 5-second 6-lat 7-long 8-depth 9-mag
			// 10-magType 11-magSource 12-magErr 13-magRounding 14-EarthquakeID
			
//			GregorianCalendar cal = rup.getOriginTimeCal();
////			cal.setTimeZone(TimeZone.getTimeZone("GMT-0:00"));
//			cal.setTimeZone(TimeZone.getTimeZone("UTC-0:00"));
//			System.out.println("Millis: "+rup.getOriginTime());
			GregorianCalendar cal = new GregorianCalendar(TimeZone.getTimeZone("UTC-0:00"));
			cal.setTimeInMillis(rup.getOriginTime());
			
			int year = cal.get(GregorianCalendar.YEAR);
			int month = cal.get(GregorianCalendar.MONTH)+1; // make it 1-based
			int day = cal.get(GregorianCalendar.DAY_OF_MONTH);
			int hour = cal.get(GregorianCalendar.HOUR_OF_DAY);
			int minute = cal.get(GregorianCalendar.MINUTE);
			int second = cal.get(GregorianCalendar.SECOND);
			
			Location hypo = rup.getHypocenterLocation();
			
			String line = year+"\t"+month+"\t"+day+"\t"+hour+"\t"+minute+"\t"+second+"\t";
			line += hypo.getLatitude()+"\t"+hypo.getLongitude()+"\t"+hypo.getDepth()+"\t";
			line += rup.getMag()+"\t0\t0\t0\t0\t"+rup.getEventId().replaceAll("[^\\d.]", ""); // strip out non numeric
			
			System.out.println("M"+(float)rup.getMag()+", "+rup.getOriginTime()+":\t"+line);
			fw.write(line+"\n");
		}
		
		fw.close();
	}

}
