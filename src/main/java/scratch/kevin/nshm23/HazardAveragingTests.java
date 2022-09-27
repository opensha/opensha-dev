package scratch.kevin.nshm23;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.faultSysSolution.hazard.mpj.MPJ_SiteLogicTreeHazardCurveCalc;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc.ReturnPeriods;

import com.google.common.base.Preconditions;

public class HazardAveragingTests {

	public static void main(String[] args) throws IOException {
		File zipFile = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_09_16-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/results_hazard_sites.zip");
		ZipFile zip = new ZipFile(zipFile);
		CSVFile<String> sitesCSV = CSVFile.readStream(zip.getInputStream(
				zip.getEntry(MPJ_SiteLogicTreeHazardCurveCalc.SITES_CSV_FILE_NAME)), true);
		
		LogicTree<?> tree = LogicTree.read(new InputStreamReader(zip.getInputStream(zip.getEntry("logic_tree.json"))));
		
		ReturnPeriods rp = ReturnPeriods.TWO_IN_50;
		
		List<Site> sites = MPJ_SiteLogicTreeHazardCurveCalc.parseSitesCSV(sitesCSV, null);
		
		MinMaxAveTracker ratioTrack = new MinMaxAveTracker();
		
		for (Site site : sites) {
			String sitePrefix = site.getName().replaceAll("\\W+", "_");
			
			Map<Double, ZipEntry> perEntries = loadSiteCSVs(sitePrefix, zip);
			Preconditions.checkState(!perEntries.isEmpty());
			
			List<Double> periods = new ArrayList<>(perEntries.keySet());
			Collections.sort(periods);
			
			System.out.println("Site: "+site.getName());
			
			for (double period : periods) {
				System.out.println("Period: "+period);
				CSVFile<String> curvesCSV = CSVFile.readStream(zip.getInputStream(perEntries.get(period)), true);
				List<DiscretizedFunc> curves = new ArrayList<>();
				List<Double> weights = new ArrayList<>();
				
				int startCol = 3+tree.getLevels().size();
				double[] xVals = new double[curvesCSV.getNumCols()-startCol];
				for (int i=0; i<xVals.length; i++) {
					xVals[i] = curvesCSV.getDouble(0, startCol+i);
				}
				
				double sumWeight = 0d;
				double avgHaz = 0d;
				double[] avgYVals = new double[xVals.length];
				for (int row=1; row<curvesCSV.getNumRows(); row++) {
					double weight = curvesCSV.getDouble(row, 2);
					double[] yVals = new double[xVals.length];
					for (int i=0; i<yVals.length; i++) {
						double y = curvesCSV.getDouble(row, startCol+i);
						yVals[i] = y;
						avgYVals[i] += y*weight;
					}
					LightFixedXFunc curve = new LightFixedXFunc(xVals, yVals);
					curves.add(curve);
					weights.add(weight);
					sumWeight += weight;
					double hazVal = hazardVal(curve, rp);
					avgHaz += hazVal*weight;
				}
				
				avgHaz /= sumWeight;
				for (int i=0; i<avgYVals.length; i++)
					avgYVals[i] /= sumWeight;
				LightFixedXFunc avgCurve = new LightFixedXFunc(xVals, avgYVals);
				double hazFromAvgCurve = hazardVal(avgCurve, rp);
				System.out.println("Average hazard val:\t"+avgHaz);
				System.out.println("Hazard from avg curve:\t"+hazFromAvgCurve);
				ratioTrack.addValue(avgHaz/hazFromAvgCurve);
			}
		}
		System.out.println("Ratios: "+ratioTrack);
	}
	
	private static Map<Double, ZipEntry> loadSiteCSVs(String sitePrefix, ZipFile zip) {
		Enumeration<? extends ZipEntry> entries = zip.entries();
		Map<Double, ZipEntry> perEntires = new HashMap<>();
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			String name = entry.getName();
			if (name.startsWith(sitePrefix)) {
				String nameLeft = name.substring(sitePrefix.length());
				double period;
				if (nameLeft.equals("_pgv.csv")) {
					period = -1;
				} else if (nameLeft.equals("_pga.csv")) {
					period = 0d;
				} else if (nameLeft.startsWith("_sa_") && nameLeft.endsWith(".csv")) {
					String perStr = nameLeft.substring(4, nameLeft.length()-4);
					period = Double.parseDouble(perStr);
				} else {
					System.err.println("Skipping unexpected file that we couldn't parse for a period: "+name);
					continue;
				}
				perEntires.put(period, entry);
			}
		}
		return perEntires;
	}
	
	private static double hazardVal(DiscretizedFunc curve, ReturnPeriods rp) {
		double curveLevel = rp.oneYearProb;
		// curveLevel is a probability, return the IML at that probability
		if (curveLevel > curve.getMaxY())
			return 0d;
		else if (curveLevel < curve.getMinY())
			// saturated
			return curve.getMaxX();
		else
			return curve.getFirstInterpolatedX_inLogXLogYDomain(curveLevel);
//			return curve.getFirstInterpolatedX(curveLevel);
	}

}
