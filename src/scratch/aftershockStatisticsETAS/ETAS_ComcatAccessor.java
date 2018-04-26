package scratch.aftershockStatisticsETAS;

import gov.usgs.earthquake.event.EventQuery;
import gov.usgs.earthquake.event.EventWebService;
import gov.usgs.earthquake.event.Format;
import gov.usgs.earthquake.event.JsonEvent;
import scratch.aftershockStatistics.ComcatAccessor;
import scratch.aftershockStatistics.cmu.EqkEventIdComparator;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.json.simple.JSONArray;
import org.apache.commons.io.IOUtils;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.opensha.commons.geo.GeoTools;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.impl.StringParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.faultSurface.FaultTrace;

import com.google.common.base.Preconditions;

public class ETAS_ComcatAccessor extends ComcatAccessor{
	
	private static final boolean D = false;
	
	public FaultTrace fetchShakemapSource(String eventID){
//		ETASEqkRupture shakeMapSource = new ETASEqkRupture();
		FaultTrace faultTrace = null;
		
		EventQuery query = new EventQuery();
		query.setEventId(eventID);
		List<JsonEvent> events;
		try {
			events = service.getEvents(query);
		} catch (FileNotFoundException e) {
			return null;
		} catch (IOException e){
			return null;
		} catch (Exception e) {
//			throw ExceptionUtils.asRuntimeException(e);
			Exception e2 = new FileNotFoundException("Could not retrieve event '"+ eventID +"' from Comcat.");
			e2.initCause(e);
			throw ExceptionUtils.asRuntimeException(e2);
		}
		if (events.isEmpty())
			return null;
		
		Preconditions.checkState(events.size() == 1, "More that 1 match? "+events.size());
		
//		JsonEvent event = events.get(0);
		JSONObject obj = events.get(0);
//		JSONParser jsonParser = new JSONParser();s
		
//		printJSON(obj);
		
		
		JSONObject prop = (JSONObject) obj.get("properties");
		JSONObject prods = (JSONObject) prop.get("products");
		JSONArray shakemaps = (JSONArray) prods.get("shakemap");
		JSONObject shakemap;
		
		boolean faultFound = false;

		if (shakemaps != null) {

			shakemap = (JSONObject) shakemaps.get(0);

			JSONObject contents = (JSONObject) shakemap.get("contents");
			JSONObject fault = new JSONObject();

			if (contents != null){
				Set<?> keys = contents.keySet();
				Iterator<?> i = keys.iterator();
				while (i.hasNext()) {
					String str = i.next().toString();
					if (str.endsWith("fault.txt")){
						if(D) System.out.println(str);
						fault = (JSONObject) contents.get(str);
						faultFound = true;
						if(D) System.out.println(fault.get("url"));
						break;
					}
				}

				if (faultFound){
					InputStream webin;
					String results = new String();
					try {
						webin = new URL(fault.get("url").toString()).openStream();
						results = IOUtils.toString(webin, "UTF-8");
						IOUtils.closeQuietly(webin);

					} catch (MalformedURLException e1) {
						// TODO Auto-generated catch block
						//				e1.printStackTrace();
						if(D) System.out.println("No Shakemap");
						return null;
					} catch (IOException e1) {
						//				 TODO Auto-generated catch block
						//				e1.printStackTrace();
						if(D) System.out.println("No Shakemap");
						return null;
					} 

					//			System.out.println(results);
					String[] lines = results.split("\n");

					faultTrace = new FaultTrace("Shakemap Source");

					for (String line : lines){
						if (line.charAt(0) != '#' && line.charAt(0) != '>'){
							String[] splitAll = line.split(" ");
							String[] split = new String[3];
							
							int n = 0;
							for(String str : splitAll){
								if (!str.isEmpty()) split[n++]=str;
							}
								
							Location loc;
							try{
								loc = new Location(new Double(split[0]), new Double(split[1]));
							} catch (Exception e) {
								System.err.println("Problem parsing finite source at line:") ;
								StringBuilder outString = new StringBuilder();
								for (String elem : split)
									outString.append(elem + ",");
								System.err.println(outString);	
								
								return null;
							}
							faultTrace.add(loc);
						}
					}
				} else {
					if(D) System.out.println("No Shakemap");
					return null;
				}
			}
		} else {
			if(D) System.out.println("No Shakemap");
			return null;
		}
		
		return faultTrace;
	}
	
	public String fetchShakemapURL(String eventID){
		String shakeMapURL = null;
		
		EventQuery query = new EventQuery();
		query.setEventId(eventID);
		List<JsonEvent> events;
		try {
			events = service.getEvents(query);
		} catch (FileNotFoundException e) {
			return null;
		} catch (IOException e){
			return null;
		} catch (Exception e) {
//			throw ExceptionUtils.asRuntimeException(e);
			Exception e2 = new FileNotFoundException("Could not retrieve event '"+ eventID +"' from Comcat.");
			e2.initCause(e);
			throw ExceptionUtils.asRuntimeException(e2);
		}
		if (events.isEmpty())
			return null;
		
		Preconditions.checkState(events.size() == 1, "More that 1 match? "+events.size());
		
		JSONObject obj = events.get(0);
		JSONObject prop = (JSONObject) obj.get("properties");
		JSONObject prods = (JSONObject) prop.get("products");
		JSONArray shakemaps = (JSONArray) prods.get("shakemap");
		JSONObject shakemap;
		
		

		if (shakemaps != null) {

			shakemap = (JSONObject) shakemaps.get(0);

			JSONObject contents = (JSONObject) shakemap.get("contents");
			JSONObject intensity = new JSONObject();

			if (contents != null){
				Set<?> keys = contents.keySet();
				Iterator<?> i = keys.iterator();
				while (i.hasNext()) {
					String str = i.next().toString();
					if (str.endsWith("intensity.jpg")){
						if(D) System.out.println(str);
						intensity = (JSONObject) contents.get(str);
						if(D) System.out.println(intensity.get("url"));
						break;
					}
				}

				shakeMapURL = intensity.get("url").toString();
			}
		} else {
			if(D) System.out.println("No Shakemap");
			return null;
		}
		
		return shakeMapURL;
	}
	
}
