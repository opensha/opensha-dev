package scratch.ned.ETAS_Tests.obsLanders;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.faultSurface.FaultTrace;
import org.opensha.sha.faultSurface.StirlingGriddedSurface;
import org.opensha.commons.gui.plot.GraphWindow;


public class LandersAnalysis {
	
	final static String path = "/Users/field/workspace/OpenSHA/dev/scratch/ned/ETAS_Tests/obsLanders/";
	StirlingGriddedSurface landersSurf1, landersSurf2, landersSurf3;
	
	int[] year, month, day, hour, minute;
	double[] second, mag, distFromMainShock;
	Location[] hypocenter;
	
	int numQks;

	/**
	 * Constructor
	 */
	public LandersAnalysis() {
		mkLandersSurfaces();
		readQuakeData(2.5);
		computeDistances();
		plotDistanceCDF();
	}
	
	private void readQuakeData(double minMag) {
		StringTokenizer st;
		String line;
		ArrayList<String> lines = new ArrayList<String>();
	    try {
			BufferedReader reader = new BufferedReader(new FileReader(path + "LandersCatalog.txt"));
			reader.readLine(); // skip header
	        while ((line = reader.readLine()) != null) {
	        	st = new StringTokenizer(line);
		        for(int r=0; r<9; r++) st.nextToken();
		        double mag = Double.valueOf(st.nextToken());
		        System.out.println(mag);
		        if(mag >= minMag)
		        	lines.add(line);	        	
	        }
       } catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        numQks = lines.size();
        
        year = new int[numQks];
        month = new int[numQks];
        day = new int[numQks];
        hour = new int[numQks];
        minute = new int[numQks];
        second = new double[numQks];
        mag = new double[numQks];
        hypocenter = new Location[numQks];
        
        for(int l=0; l<numQks; l++) {
        	st = new StringTokenizer(lines.get(l));
        	year[l] = Integer.valueOf(st.nextToken());
        	month[l] = Integer.valueOf(st.nextToken());
        	day[l] = Integer.valueOf(st.nextToken());
        	hour[l] = Integer.valueOf(st.nextToken());
        	minute[l] = Integer.valueOf(st.nextToken());
        	second[l] = Double.valueOf(st.nextToken());
        	hypocenter[l] = new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken()));
        	mag[l] = Double.valueOf(st.nextToken());
        }
//		System.out.println(hypocenter[0]);
//		System.out.println(hypocenter[numQks-1]);
        System.out.println("numQks="+numQks);

	}
	
	private void mkLandersSurfaces() {
		StringTokenizer st;
		String line;
        try {
			BufferedReader reader = new BufferedReader(new FileReader(path + "shakeMapLandersTrace1.txt"));
			reader.readLine();	// skip header line
			FaultTrace trace1 = new FaultTrace("Landers Surf 1");
			st = new StringTokenizer(reader.readLine());
			trace1.add(new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken())));
			st = new StringTokenizer(reader.readLine());
			trace1.add(new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken())));
			System.out.println(trace1);
			landersSurf1 = new StirlingGriddedSurface(trace1, 90.0, 0.0, 15.0, 1.0);

			reader = new BufferedReader(new FileReader(path + "shakeMapLandersTrace2.txt"));
			reader.readLine();	// skip header line
			FaultTrace trace2 = new FaultTrace("Landers Surf 2");
			st = new StringTokenizer(reader.readLine());
			trace2.add(new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken())));
			st = new StringTokenizer(reader.readLine());
			trace2.add(new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken())));
			System.out.println(trace2);
			landersSurf2 = new StirlingGriddedSurface(trace2, 90.0, 0.0, 15.0, 1.0);
		

			reader = new BufferedReader(new FileReader(path + "shakeMapLandersTrace3.txt"));
			reader.readLine();	// skip header line
			FaultTrace trace3 = new FaultTrace("Landers Surf 3");
			st = new StringTokenizer(reader.readLine());
			trace3.add(new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken())));
			st = new StringTokenizer(reader.readLine());
			trace3.add(new Location(Double.valueOf(st.nextToken()),Double.valueOf(st.nextToken())));
			System.out.println(trace3);
			landersSurf3 = new StirlingGriddedSurface(trace3, 90.0, 0.0, 15.0, 1.0);

        } catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public void computeDistances() {
		distFromMainShock = new double[numQks];
		double minDist;
		for(int qk=0; qk<numQks; qk++) {
			Location loc = hypocenter[qk];
			minDist = LocationUtils.distanceToSurf(loc, landersSurf1);
			double dist = LocationUtils.distanceToSurf(loc, landersSurf2);
			if(dist < minDist) minDist = dist;
			dist = LocationUtils.distanceToSurf(loc, landersSurf3);
			if(dist < minDist) minDist = dist;
			distFromMainShock[qk] = minDist;
//			System.out.println(minDist);
		}
		Arrays.sort(distFromMainShock);
	}
	
	public void plotDistanceCDF() {
		ArbitrarilyDiscretizedFunc distCDF = new ArbitrarilyDiscretizedFunc();
		distCDF.setName("Distance to Landers CDF");
		
		for(int qk=0; qk<numQks; qk++) {
			distCDF.set(distFromMainShock[qk], (double)(qk+1)/(double)numQks);
		}
		
		ArrayList funcs = new ArrayList();
		funcs.add(distCDF);
		
		GraphWindow graph = new GraphWindow(funcs, " "); 
		graph.setX_AxisLabel("Distance (km)");
		graph.setY_AxisLabel(" ");
//		graph.setX_AxisRange(0.4, 360);
//		graph.setY_AxisRange(0.1, graph.getY_AxisMax());
//		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
//		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLUE));
//		graph.setPlottingFeatures(plotChars);
//		graph.setYLog(true);
//		graph.setXLog(true);
/*		try {
			graph.saveAsPDF(dirToSaveData+"numAshocksVsTime.pdf");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
*/
	}

	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		LandersAnalysis doit = new LandersAnalysis();
		
	}


}
