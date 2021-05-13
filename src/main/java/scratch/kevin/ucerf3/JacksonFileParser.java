package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.faultSurface.ApproxEvenlyGriddedSurface;
import org.opensha.sha.faultSurface.EvenlyGriddedSurface;
import org.opensha.sha.faultSurface.FaultTrace;

public class JacksonFileParser {
	
	public static void parseJacksonFile(File file) throws IOException {
		ArrayList<String> lines = FileUtils.loadFile(file.getAbsolutePath(), false);
		
		HashMap<String, FaultTrace> uppersMap = new HashMap<String, FaultTrace>();
		HashMap<String, FaultTrace> lowersMap = new HashMap<String, FaultTrace>();
		
		for (String line : lines) {
			line = line.trim();
			if (line.isEmpty())
				continue;
			
			StringTokenizer tok = new StringTokenizer(line);
			
			String id = tok.nextToken()+" "+tok.nextToken()+" "+tok.nextToken()+" "+tok.nextToken()+" "+tok.nextToken();
			
			if (!uppersMap.containsKey(id)) {
				uppersMap.put(id, new FaultTrace(id+" (upper)"));
				lowersMap.put(id, new FaultTrace(id+" (lower)"));
			}
			FaultTrace upper = uppersMap.get(id);
			FaultTrace lower = lowersMap.get(id);
			
			// degrees
			double lat = Double.parseDouble(tok.nextToken());
			// minutes
			lat += Double.parseDouble(tok.nextToken()) / 60d;
			
			// degrees
			double lon = Double.parseDouble(tok.nextToken());
			// minutes
			lon += Double.parseDouble(tok.nextToken()) / 60d;
			// make it negative
			lon = -lon;
			
			tok.nextToken(); // skip
			tok.nextToken(); // skip
			tok.nextToken(); // skip
			double strike = Double.parseDouble(tok.nextToken());
			double dip = Double.parseDouble(tok.nextToken());
			
//          instead of lat(i),lon(i), use lat(i)+/- (6 km *sin(strike)*tan(90-dip)) and lon(i)+/- (6 km *cos(strike)*tan(90-dip))
//   		(if the two points don't have the same strike, compute an average strike
//			or use the azimuth between (lat(i),lon(i)) and (lat(i+1),lon(i+1)))
			
			System.out.println("mid:\tlat:\t"+lat+"\tlon:\t"+lon);
			double latDelta = 6d * Math.sin(Math.toRadians(strike)) * Math.tan(Math.toRadians(90-dip));
			double lonDelta = 6d * Math.cos(Math.toRadians(strike)) * Math.tan(Math.toRadians(90-dip));
			
			Location upperLoc = new Location(lat, lon, 0d);
			Location lowerLoc = new Location(lat, lon, 12d);
			
			upperLoc = LocationUtils.location(upperLoc, 0, latDelta);
			lowerLoc = LocationUtils.location(lowerLoc, 0, -latDelta);
			
			upperLoc = LocationUtils.location(upperLoc, Math.PI*0.5, -lonDelta);
			lowerLoc = LocationUtils.location(lowerLoc, Math.PI*0.5, lonDelta);
			
//			upperLoc = LocationUtils.location(upperLoc, 0, -latDelta);
//			lowerLoc = LocationUtils.location(lowerLoc, 0, latDelta);
//			
//			upperLoc = LocationUtils.location(upperLoc, Math.PI*0.5, lonDelta);
//			lowerLoc = LocationUtils.location(lowerLoc, Math.PI*0.5, -lonDelta);
			
//			double upperLat = lat + ;
//			double lowerLat = lat - (6d * Math.sin(Math.toRadians(strike)) * Math.tan(Math.toRadians(90-dip)));
//			
//			double upperLon = lon  - (6d * Math.cos(Math.toRadians(strike)) * Math.tan(Math.toRadians(90-dip)));
//			double lowerLon = lon  + (6d * Math.cos(Math.toRadians(strike)) * Math.tan(Math.toRadians(90-dip)));
//			System.out.println("up:\tlat:"+upperLat+"\tlon:\t"+upperLon);
//			System.out.println("dwn:\tlat:"+lowerLat+"\tlon:\t"+lowerLon);
//			
//			Location upperPt = new Location(upperLat, upperLon, 0);
//			Location lowerPt = new Location(lowerLat, lowerLon, 12);
			
			upper.add(upperLoc);
			lower.add(lowerLoc);
		}
		
		int cnt = 0;
		for (String id : uppersMap.keySet()) {
			FaultTrace upper = uppersMap.get(id);
			FaultTrace lower = lowersMap.get(id);
			
			EvenlyGriddedSurface surf = new ApproxEvenlyGriddedSurface(upper, lower, 1d);
			
			CatalogWriter.writeFiniteSurfaceFile(surf, new File("/tmp/surf_"+id.replaceAll(" ", "_")+".txt"), 0d);
		}
	}
	
	public static void parseFourCornersFile(File file) throws IOException {
		ArrayList<String> lines = FileUtils.loadFile(file.getAbsolutePath());
		String prevComment = null;
		String curID = null;
		int curCNT = 0;
		
		Location[] locs = null;
		for (String line : lines) {
			if (line.startsWith("#")) {
				prevComment = line;
				continue;
			}
			if (curCNT == 0) {
				StringTokenizer idTok = new StringTokenizer(prevComment);
				idTok.nextToken();
				curID = idTok.nextToken()+"_"+idTok.nextToken()+"_"+idTok.nextToken()+"_"+idTok.nextToken()+"_"+idTok.nextToken();
				locs = new Location[4];
			}
			
			String[] split = line.split(" ");
			Location loc = new Location(Double.parseDouble(split[1]), Double.parseDouble(split[0]), Double.parseDouble(split[2]));
			
			locs[curCNT] = loc;
			
			curCNT++;
			
			if (curCNT == 4) {
				FaultTrace upper = new FaultTrace(curID);
				FaultTrace lower = new FaultTrace(curID);
				
				lower.add(locs[0]);
				lower.add(locs[1]);
				upper.add(locs[3]);
				upper.add(locs[2]);
				
				ApproxEvenlyGriddedSurface surf = new ApproxEvenlyGriddedSurface(upper, lower, 1d);
				
				CatalogWriter.writeFiniteSurfaceFile(surf, new File("/tmp/surf_"+curID+".txt"), 0d);
				
				curID = null;
				curCNT = 0;
				locs = null;
			}
		}
	}
	
	public static void main(String[] args) throws IOException {
		parseJacksonFile(new File("/home/kevin/OpenSHA/UCERF3/cal_extended.dat"));
//		parseFourCornersFile(new File("/tmp/Finite_Faults.txt"));
	}

}
