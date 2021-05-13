package scratch.kevin.magDepth;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;

import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupture;
import org.opensha.sha.earthquake.observedEarthquake.parsers.UCERF3_CatalogParser;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

public class DepthPlotMaker {
	
	private static class Mapping {
		private double depth, mag;
		private Mapping(double depth, double mag) {
			this.depth = depth;
			this.mag = mag;
		}
	}
	
	private static List<Mapping> loadSCFormat(File file) throws IOException {
		ArrayList<Mapping> list = new ArrayList<Mapping>();
		for (String line : FileUtils.loadFile(file.getAbsolutePath())) {
			StringTokenizer tok = new StringTokenizer(line.trim());
			if (tok.countTokens() < 5)
				continue;
			tok.nextToken();
			tok.nextToken();
			tok.nextToken();
			tok.nextToken();
			tok.nextToken();
			tok.nextToken();
			tok.nextToken();
			double lon = Double.parseDouble(tok.nextToken());
			double lat = Double.parseDouble(tok.nextToken());
			double depth = Double.parseDouble(tok.nextToken());
			double mag = Double.parseDouble(tok.nextToken());
			if (depth != Math.floor(depth))
				list.add(new Mapping(depth, mag));
		}
		
		return list;
	}
	
	private static List<Mapping> loadHauksson(double maxError, File file) throws IOException {
		ArrayList<Mapping> list = new ArrayList<Mapping>();
		int skipped = 0;
		int tot = 0;
		for (String line : FileUtils.loadFile(file.getAbsolutePath())) {
			ArrayList<String> vals = Lists.newArrayList(Splitter.on(' ').omitEmptyStrings().split(line));
			try {
				double depth = Double.parseDouble(vals.get(9));
				double mag = Double.parseDouble(vals.get(10));
				if (mag > 7)
					System.out.println(line);
				double depthError = Double.parseDouble(vals.get(20));
				if (vals.get(23).equals("qb"))
					continue;
				if (vals.get(24).equals("xx"))
					continue;
				if (depthError > maxError)
					skipped++;
				else
					list.add(new Mapping(depth, mag));
				tot++;
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.out.println(line);
				System.out.println(vals.get(9));
				System.out.println(vals.get(10));
				System.out.println(vals.get(20));
				System.exit(0);
			}
		}
		
		System.out.println("Skipped "+skipped+"/"+tot+" ("+(float)(100d * (double)skipped / (double)tot)+" %)");
		
		return list;
	}
	
	private static List<Mapping> loadUCERF3(int lastYear, File file) throws IOException {
		GregorianCalendar cal = new GregorianCalendar(lastYear, 0, 0);
		ObsEqkRupList cat = UCERF3_CatalogParser.loadCatalog(file);
		
		ArrayList<Mapping> list = new ArrayList<Mapping>();
		for (ObsEqkRupture rup : cat.getRupsAfter(cal.getTimeInMillis())) {
			double depth = rup.getHypocenterLocation().getDepth();
//			if (depth != 6)
			if (depth != Math.floor(depth))
				list.add(new Mapping(depth, rup.getMag()));
			if (depth < 6 && rup.getMag() > 6)
				System.out.println("UCERF3: "+rup.getInfo());
		}
		
		return list;
	}
	
	private static class DepthStatTracker implements Comparable<DepthStatTracker>{
		int num;
		double depth;
		public DepthStatTracker(double depth) {
			this.depth = depth;
			this.num = 1;
		}
		
		public void increment() {
			this.num++;
		}

		@Override
		public int compareTo(DepthStatTracker o) {
			return -Double.compare(num, o.num);
		}
	}
	
	private static List<Mapping> loadANSS(int minStations, File file) throws NumberFormatException, FileNotFoundException, IOException {
		ArrayList<String> lines = FileUtils.loadFile(file.getAbsolutePath());
		ArrayList<Mapping> list = new ArrayList<Mapping>();
		
		HashMap<Double, DepthStatTracker> trackers = new HashMap<Double, DepthStatTracker>();
		
		for (String line : lines.subList(2, lines.size())) {
			StringTokenizer tok = new StringTokenizer(line.trim());
			if (tok.countTokens() < 5)
				continue;
			tok.nextToken(); // date
			tok.nextToken(); // time
			tok.nextToken(); // lat 
			tok.nextToken(); // lon
			double depth = Double.parseDouble(tok.nextToken());
			double mag = Double.parseDouble(tok.nextToken());
			if (mag > 6.8)
				System.out.println(line);
			DepthStatTracker track = trackers.get(depth);
			if (track == null)
				trackers.put(depth, new DepthStatTracker(depth));
			else
				track.increment();
			
			if (depth == Math.floor(depth))
				continue;
//			if (depth == 10)
//				continue;
//			if (depth == 5)
//				continue;
//			if (depth == 6)
//				continue;
//			if (mag > 7 && depth <= 1) {
//				System.out.println(line);
//				System.exit(0);
//			}
			tok.nextToken(); // Magt
			try {
				int nst = Integer.parseInt(tok.nextToken());
				if (nst < minStations)
					continue;
			} catch (Exception e) {
				continue;
			}
			
			list.add(new Mapping(depth, mag));
		}
		
//		ArrayList<DepthStatTracker> tlist = Lists.newArrayList();
//		tlist.addAll(trackers.values());
//		Collections.sort(tlist);
//		for (int i=0; i<20; i++) {
//			DepthStatTracker d = tlist.get(i);
//			System.out.println("DEPTH: "+d.depth+"\tNUM: "+d.num);
//		}
		
		return list;
	}
	
	private static void plot(List<Mapping> list, double maxDepth) {
		if (maxDepth <= 0)
			maxDepth = Double.POSITIVE_INFINITY;
		System.out.println("Plotting "+list.size()+" mappings");
		DefaultXY_DataSet func = new DefaultXY_DataSet();
		for (Mapping val : list) {
			if (val.depth < maxDepth) {
				func.set(val.mag, -val.depth);
			}
		}
		System.out.println("Done setting func, displaying.");
		ArrayList<XY_DataSet> funcs = Lists.newArrayList();
		funcs.add(func);
		ArrayList<PlotCurveCharacterstics> plotChars = Lists.newArrayList(new PlotCurveCharacterstics(PlotSymbol.X, 2f, Color.BLACK));
		GraphWindow gw = new GraphWindow(funcs, "Mag Vs. Depth", plotChars);
		gw.setX_AxisLabel("Magnitude");
		gw.setY_AxisLabel("Altitude");
	}
	
	private static void writeCatalog(List<Mapping> cat, File file) throws IOException {
		FileWriter fw = new FileWriter(file);
		for (Mapping eq : cat) {
			fw.write(eq.mag+"\t"+eq.depth+"\n");
		}
		fw.close();
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File landersFile = new File("/home/kevin/workspace/scec_vdo_ucerf3/data/Catalogs/landers_etas_new.sc");
		File ucerf3File = new File("/home/kevin/OpenSHA/UCERF3/Felzer_UCERF3_Catalog4_0.txt");
		File anssFile = new File("/home/kevin/Documents/geol_440_eq_depth/anss_cat_1986");
		File haukFile = new File("/home/kevin/Documents/geol_440_eq_depth/hs_1981_2011_06_comb_K2_A.cat_so_SCSN_v01");
		File outDir = new File("/home/kevin/Documents/geol_440_eq_depth/");
		List<Mapping> hauk = loadHauksson(1d, haukFile);
		writeCatalog(hauk, new File(outDir, "hauk_output.txt"));
		List<Mapping> anss = loadANSS(10, anssFile);
		writeCatalog(anss, new File(outDir, "anss_output.txt"));
//		plot(anss, 15);
		List<Mapping> landersList = loadSCFormat(landersFile);
		writeCatalog(landersList, new File(outDir, "landers_output.txt"));
//		plot(landersList, 15);
		List<Mapping> ucerf3List = loadUCERF3(1986, ucerf3File);
		writeCatalog(ucerf3List, new File(outDir, "ucerf3_output.txt"));
//		plot(ucerf3List, 15);
	}

}
