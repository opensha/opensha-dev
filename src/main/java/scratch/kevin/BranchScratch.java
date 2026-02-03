package scratch.kevin;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public class BranchScratch {
	
	private static void test1() throws IOException {
		MeanUCERF2 u2 = new MeanUCERF2();
		u2.getParameter(UCERF2.PROB_MODEL_PARAM_NAME).setValue(UCERF2.PROB_MODEL_POISSON);
		u2.getTimeSpan().setDuration(1d);
		u2.updateForecast();
		
		int sourceID = 128;
		ProbEqkSource source = u2.getSource(sourceID);
		
		System.out.println("Source has "+source.getNumRuptures()+" ruptures");
		
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(5.05, 39, 0.1);
		int[] counts = new int[mfd.size()];
		for (ProbEqkRupture rup : source) {
			int index = mfd.getClosestXIndex(rup.getMag());
			mfd.add(index, rup.getMeanAnnualRate(1d));
			counts[index]++;
		}
		for (int i=0; i<mfd.size(); i++)
			System.out.println("M"+(float)mfd.getX(i)+"\t"+counts[i]+" rups, rate="+(float)mfd.getY(i));
		
		PlotSpec plot = new PlotSpec(List.of(mfd),
				List.of(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK)),
				"UCERF2 Source "+sourceID+": "+source.getName(), "Magnitude", "Annual Rate");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, true, new Range(6d, 9d), new Range(1e-6, 1e0));
		
		PlotUtils.writePlots(new File("/tmp"), "u2_source_"+sourceID, gp, 800, 800, true, false, false);
	}
	
	private static void test2() throws IOException {
		FaultSystemSolution fss = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		
		Location loc = new Location(37.760, -121.937);
		double radius = 10;
		Region reg = new Region(loc, radius);
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, 8.49);
		FaultSystemRupSet rupSet = fss.getRupSet();
		double[] rupFracts = rupSet.getFractRupsInsideRegion(reg, false);
		List<Integer> rupIndexes = new ArrayList<>();
		for (int i=0; i<rupFracts.length; i++)
			if (rupFracts[i] > 0)
				rupIndexes.add(i);
		IncrementalMagFreqDist particMFD = fss.calcParticipationMFD_forRups(rupIndexes, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		IncrementalMagFreqDist nuclMFD = fss.calcNucleationMFD_forRegion(reg, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size(), false);
		
		CSVFile<String> incrCSV = new CSVFile<>(true);
		incrCSV.addLine("Magnitude", "Nucleation Rate (1/yr)", "Participation Rate (1/yr)");
		for (int i=0; i<refMFD.size(); i++)
			incrCSV.addLine((float)refMFD.getX(i)+"", (float)nuclMFD.getY(i)+"", (float)particMFD.getY(i)+"");
		
		incrCSV.writeToFile(new File("/tmp/san_ramon_10km_on_fault_incr.csv"));
		
		CSVFile<String> cmlCSV = new CSVFile<>(true);
		cmlCSV.addLine("Magnitude", "Nucleation Rate (1/yr)", "Participation Rate (1/yr)");
		EvenlyDiscretizedFunc nuclCmlMFD = nuclMFD.getCumRateDistWithOffset();
		EvenlyDiscretizedFunc particCmlMFD = particMFD.getCumRateDistWithOffset();
		for (int i=0; i<nuclCmlMFD.size(); i++)
			cmlCSV.addLine((float)nuclCmlMFD.getX(i)+"", (float)nuclCmlMFD.getY(i)+"", (float)particCmlMFD.getY(i)+"");
		cmlCSV.writeToFile(new File("/tmp/san_ramon_10km_on_fault_cml.csv"));
		
		double[] sectFracts = rupSet.getFractSectsInsideRegion(reg, false);
		for (int s=0; s<sectFracts.length; s++) {
			if (sectFracts[s] > 0) {
				FaultSection sect = rupSet.getFaultSectionData(s);
				System.out.println(sect.getSectionName()+" is "+(float)sect.getFaultSurface(1d).getDistanceJB(loc)+" km away");
			}
		}
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		try {
			test2();
		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(1);
		}
	}

}
