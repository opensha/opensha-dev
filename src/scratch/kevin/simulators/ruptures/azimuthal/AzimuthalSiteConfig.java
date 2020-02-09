package scratch.kevin.simulators.ruptures.azimuthal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;

import scratch.kevin.bbp.BBP_Site;
import scratch.kevin.simulators.ruptures.BBP_PartBValidationConfig.Scenario;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;
import scratch.kevin.bbp.BBP_Module.VelocityModel;

public abstract class AzimuthalSiteConfig<E> {
	
	private Scenario scenario;
	private List<E> ruptures;
	private EvenlyDiscrXYZ_DataSet gc2xyz;
	
	protected AzimuthalSiteConfig(Scenario scenario, List<E> ruptures) {
		this.scenario = scenario;
		this.ruptures = ruptures;
	}
	
	protected void init(double spacing, double buffer, double length, boolean positiveX) {
		int bufferPoints = (int)Math.round(buffer/spacing);
		System.out.println(bufferPoints+" buffer points");

		int ny = 2*bufferPoints + (int)Math.round(length/spacing) + 1;
		double minY = -bufferPoints*spacing;
		
		int nx;
		double minX;
		if (positiveX) {
			minX = 0d;
			nx = bufferPoints+1;
		} else {
			minX = -bufferPoints*spacing;
			nx = 2*bufferPoints+1;
		}
		
		System.out.println("minX="+(float)minX+", minY="+(float)minX+", nx="+nx+", ny="+ny);
		
		init(new EvenlyDiscrXYZ_DataSet(nx, ny, minX, minY, spacing));
	}
	
	protected void init(EvenlyDiscrXYZ_DataSet gc2xyz) {
		this.gc2xyz = gc2xyz;
	}
	
	public Scenario getScenario() {
		return scenario;
	}
	
	public List<E> getRuptures() {
		return ruptures;
	}
	
	public final EvenlyDiscrXYZ_DataSet getGC2XYZ() {
		return gc2xyz;
	}
	
	public abstract List<Location> getRuptureSiteLocs(E rupture);
	
	public abstract double getHypocenterDAS(E rupture);
	
	public abstract double getHypocenterDDW(E rupture);
	
	public abstract double getLength(E rupture);
	
	public abstract double getWidth(E rupture);

	public abstract double getDip(E rupture);

	public abstract double getMag(E rupture);

	public abstract int getID(E rupture);
	
	public Point2D[] getRectangle(E rupture) {
		double length = getLength(rupture);
		double width = getWidth(rupture);
		double dip = getDip(rupture);
		double horzWidth = 0d;
		if (dip < 90d)
			horzWidth = Math.cos(Math.toRadians(dip))*width;
		
		Point2D[] ret = new Point2D[4];
		ret[0] = new Point2D.Double(0d, 0d);
		ret[1] = new Point2D.Double(0d, length);
		ret[2] = new Point2D.Double(horzWidth, length);
		ret[3] = new Point2D.Double(horzWidth, 0d);
		
		return ret;
	}
	
	public List<BBP_Site> getSites(E rupture, VelocityModel vm) {
		List<BBP_Site> sites = new ArrayList<>();
		
		List<Location> siteLocs = getRuptureSiteLocs(rupture);
		for (int i=0; i<siteLocs.size(); i++) {
			Location loc = siteLocs.get(i);
			
			BBP_Site site = new BBP_Site("s"+i, loc, vm.getVs30(),
					RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ);
			sites.add(site);
		}
		
		return sites;
	}
	
	public void writeGC2LocationCSV(File csvFile) throws IOException {
		EvenlyDiscrXYZ_DataSet xyz = getGC2XYZ();
		CSVFile<String> csv = new CSVFile<>(true);
		
		csv.addLine("Index", "X index", "Y index", "RX", "RY");
		for (int i=0; i<xyz.size(); i++) {
			Point2D pt = xyz.getPoint(i);
			int xIndex = xyz.getXIndex(pt.getX());
			int yIndex = xyz.getYIndex(pt.getY());
			csv.addLine(i+"", xIndex+"", yIndex+"", pt.getX()+"", pt.getY()+"");
		}
		csv.writeToFile(csvFile);
	}
	
	protected static EvenlyDiscrXYZ_DataSet loadGC2LocationCSV(File csvFile) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		int nx = 0;
		int ny = 0;
		double minX = Double.POSITIVE_INFINITY;
		double minY = Double.POSITIVE_INFINITY;
		double minSpacing = Double.POSITIVE_INFINITY;
		
		Point2D prev = null;
		
		for (int row=1; row<csv.getNumRows(); row++) {
			int xIndex = csv.getInt(row, 1);
			int yIndex = csv.getInt(row, 2);
			nx = Integer.max(nx, xIndex+1);
			ny = Integer.max(ny, yIndex+1);
			double rx = csv.getDouble(row, 3);
			double ry = csv.getDouble(row, 4);
			minX = Math.min(minX, rx);
			minY = Math.min(minY, ry);
			
			if (prev != null) {
				double xDist = rx - prev.getX();
				double yDist = ry - prev.getY();
				double dist = Math.sqrt(xDist*xDist + yDist*yDist);
				minSpacing = Math.min(minSpacing, dist);
			}
			prev = new Point2D.Double(rx, ry);
		}
		
		return new EvenlyDiscrXYZ_DataSet(nx, ny, minX, minY, minSpacing);
	}
	
	public void writeRupturesCSV(File csvFile) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		
		csv.addLine("ID", "Mag", "Length", "Width", "Dip", "Hypocenter DAS", "Hypocenter DDW");
		for (E rupture : ruptures) {
			List<String> line = new ArrayList<>();
			line.add(getID(rupture)+"");
			line.add((float)getMag(rupture)+"");
			line.add((float)getLength(rupture)+"");
			line.add((float)getWidth(rupture)+"");
			line.add((float)getDip(rupture)+"");
			line.add((float)getHypocenterDAS(rupture)+"");
			line.add((float)getHypocenterDDW(rupture)+"");
			csv.addLine(line);
		}
		csv.writeToFile(csvFile);
	}
	
	public static AzimuthalSiteConfig<Integer> loadGeneric(Scenario scenario, File locsCSVFile, File rupsCSVFile)
			throws IOException {
		EvenlyDiscrXYZ_DataSet gc2xyz = loadGC2LocationCSV(locsCSVFile);
		
		List<Integer> ids = new ArrayList<>();
		List<Double> mags = new ArrayList<>();
		List<Double> lengths = new ArrayList<>();
		List<Double> widths = new ArrayList<>();
		List<Double> dips = new ArrayList<>();
		List<Double> hypoDASs = new ArrayList<>();
		List<Double> hypoDDWs = new ArrayList<>();
		
		CSVFile<String> csv = CSVFile.readFile(rupsCSVFile, true);
		
		for (int row=1; row<csv.getNumRows(); row++) {
			int col = 0;
			ids.add(csv.getInt(row, col++));
			mags.add(csv.getDouble(row, col++));
			lengths.add(csv.getDouble(row, col++));
			widths.add(csv.getDouble(row, col++));
			dips.add(csv.getDouble(row, col++));
			hypoDASs.add(csv.getDouble(row, col++));
			if (col < csv.getNumCols())
				hypoDDWs.add(csv.getDouble(row, col));
			else
				hypoDDWs.add(Double.NaN);
		}
		
		return new GenericAzimuthalSiteConfig(
				scenario, gc2xyz, ids, mags, lengths, widths, dips, hypoDASs, hypoDDWs);
	}
	
	private static class GenericAzimuthalSiteConfig extends AzimuthalSiteConfig<Integer> {
		
		private Map<Integer, Integer> rupToIndex;
		
		private List<Double> mags;
		private List<Double> lengths;
		private List<Double> widths;
		private List<Double> dips;
		private List<Double> hypoDASs;
		private List<Double> hypoDDWs;

		public GenericAzimuthalSiteConfig(Scenario scenario, EvenlyDiscrXYZ_DataSet gc2xyz,
				List<Integer> ids, List<Double> mags, List<Double> lengths, List<Double> widths,
				List<Double> dips, List<Double> hypoDASs, List<Double> hypoDDWs) {
			super(scenario, ids);
			
			rupToIndex = new HashMap<>();
			for (int i=0; i<ids.size(); i++)
				rupToIndex.put(ids.get(i), i);
			
			this.mags = mags;
			this.lengths = lengths;
			this.widths = widths;
			this.dips = dips;
			this.hypoDASs = hypoDASs;
			this.hypoDDWs = hypoDDWs;
			init(gc2xyz);
		}

		@Override
		public List<Location> getRuptureSiteLocs(Integer rupture) {
			throw new UnsupportedOperationException("Not supported for generic case");
		}

		@Override
		public double getHypocenterDAS(Integer rupture) {
			return hypoDASs.get(rupToIndex.get(rupture));
		}

		@Override
		public double getHypocenterDDW(Integer rupture) {
			return hypoDDWs.get(rupToIndex.get(rupture));
		}

		@Override
		public double getLength(Integer rupture) {
			return lengths.get(rupToIndex.get(rupture));
		}

		@Override
		public double getWidth(Integer rupture) {
			return widths.get(rupToIndex.get(rupture));
		}

		@Override
		public double getDip(Integer rupture) {
			return dips.get(rupToIndex.get(rupture));
		}

		@Override
		public double getMag(Integer rupture) {
			return mags.get(rupToIndex.get(rupture));
		}

		@Override
		public int getID(Integer rupture) {
			return rupture;
		}
		
	}

}
