package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.dom4j.DocumentException;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.calc.ERF_Calculator;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.faultSurface.EvenlyGriddedSurfFromSimpleFaultData;
import org.opensha.sha.faultSurface.RupInRegionCache;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.magdist.SummedMagFreqDist;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.CompoundFSSPlots;
import scratch.UCERF3.analysis.CompoundFSSPlots.RupInRegionsCache;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.utils.FaultSystemIO;

public class TimeDepFunStats {

	public static void main(String[] args) throws IOException, DocumentException {
		int[] parents = { 285, 188, 287, 286, 301, 282, 283, 284, 295, 170 };
		HashSet<Integer> parentsSet = new HashSet<Integer>();
		for (int parent : parents)
			parentsSet.add(parent);
		
		FaultSystemSolution sol = FaultSystemIO.loadSol(new File("/home/kevin/workspace/OpenSHA/"
				+ "dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		int duration = 100;
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
		erf.getTimeSpan().setDuration(duration);
		
		List<Double> probs = Lists.newArrayList();
		
		erf.updateForecast();
		
		for (ProbEqkSource source : erf) {
			String srcName = source.getName();
			int rupIndex = Integer.parseInt(srcName.substring(srcName.indexOf("#")+1, srcName.indexOf(";")));
			double mag = sol.getRupSet().getMagForRup(rupIndex);
			if (mag < 7)
				continue;
			boolean match = false;
			for (FaultSectionPrefData sect : sol.getRupSet().getFaultSectionDataForRupture(rupIndex)) {
				if (parentsSet.contains(sect.getParentSectionId())) {
					match = true;
					break;
				}
			}
			if (match)
				probs.add(source.computeTotalProb());
		}
		
		double totOneMinus = 1;
		for (double prob : probs) {
			totOneMinus *= (1-prob);
		}
		double totProb = 1 - totOneMinus;
		System.out.println(duration+" year prob: "+totProb);
		
		Region reg = new CaliforniaRegions.RELM_SOCAL();
		
		double minMag = 6.5;
		double deltaMag = 0.1;
		int numMag = 16;
		double minDuration = 5d;
		double deltaDuration = 5d;
		int numDuration = 20;
		
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(numDuration, numMag, minDuration, minMag, deltaDuration, deltaMag);
		
		RupInRegionsCache rupInRegionCache = new CompoundFSSPlots.RupInRegionsCache();
		
		for (int x=0; x<xyz.getNumX(); x++) {
			double myDuration = xyz.getX(x);
			erf.getTimeSpan().setDuration(myDuration);
			erf.updateForecast();
			List<List<Double>> probsList = Lists.newArrayList();
			for (int m=0; m<numMag; m++)
				probsList.add(new ArrayList<Double>());
			for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
				ProbEqkSource source = erf.getSource(sourceID);
				if (!rupInRegionCache.isRupInRegion(erf, source, source.getRupture(0), sourceID, 0, reg))
					continue;
				String srcName = source.getName();
				int rupIndex = Integer.parseInt(srcName.substring(srcName.indexOf("#")+1, srcName.indexOf(";")));
				double mag = sol.getRupSet().getMagForRup(rupIndex);
				int magIndex = xyz.getYIndex(mag);
				double sourceProb = source.computeTotalProb();
				for (int m=0; m<magIndex && m<numMag; m++) {
					probsList.get(m).add(sourceProb);
				}
			}
//			SummedMagFreqDist mfd = ERF_Calculator.getParticipationMagFreqDistInRegion(
//					erf, reg, minMag, numMag, deltaMag, true, rupInRegionCache);
//			ERF_Calculator.get
			System.out.println("Duration: "+myDuration);
			for (int y=0; y<xyz.getNumY(); y++) {
				probs = probsList.get(y);
				totOneMinus = 1;
				for (double prob : probs) {
					totOneMinus *= (1-prob);
				}
				totProb = 1 - totOneMinus;
				xyz.set(x, y, totProb*100d);
				System.out.println("\tM "+xyz.getY(y)+"+ Prob: "+(float)(totProb*100d)+" %");
			}
		}
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		cpt = cpt.rescale(0d, 100d);
		XYZPlotSpec spec = new XYZPlotSpec(xyz, cpt, "UCERF3 SoCal Time Dep Probability", "Duration", "Min Mag", "Probability (%)");
		new XYZPlotWindow(spec);
		
		MeanUCERF2 u2erf = new MeanUCERF2();
		u2erf.getTimeSpan().setStartTime(2013);
		EvenlyDiscrXYZ_DataSet u2XYZ = new EvenlyDiscrXYZ_DataSet(numDuration, numMag, minDuration, minMag, deltaDuration, deltaMag);
		RupInRegionCache u2Cache = new RupInRegionCache() {
			
			Table<Integer, Integer, Boolean> cache = HashBasedTable.create();
			
			@Override
			public boolean isRupInRegion(ERF erf, ProbEqkSource source, EqkRupture rup,
					int srcIndex, int rupIndex, Region region) {
				Boolean val = cache.get(srcIndex, rupIndex);
				if (val  == null) {
					val = false;
					RuptureSurface surf = source.getSourceSurface();
					List<Location> locs;
					if (surf instanceof EvenlyGriddedSurfFromSimpleFaultData)
						locs = ((EvenlyGriddedSurfFromSimpleFaultData)surf).getFaultTrace();
					else
						locs = surf.getEvenlyDiscritizedUpperEdge();
					for (Location loc : locs) {
						// good enough to use trace for large region
						if (region.contains(loc)) {
							val = true;
							break;
						}
					}
					cache.put(srcIndex, rupIndex, val);
				}
				return val;
			}
		};
		for (int x=0; x<xyz.getNumX(); x++) {
			double myDuration = xyz.getX(x);
			u2erf.getTimeSpan().setDuration(myDuration);
			u2erf.updateForecast();
			List<List<Double>> probsList = Lists.newArrayList();
			for (int m=0; m<numMag; m++)
				probsList.add(new ArrayList<Double>());
			for (int sourceID=0; sourceID<u2erf.getNumSources(); sourceID++) {
				ProbEqkSource source = u2erf.getSource(sourceID);
				for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
					ProbEqkRupture rup = source.getRupture(rupID);
					if (!u2Cache.isRupInRegion(erf, source, source.getRupture(0), sourceID, rupID, reg))
						continue;
					double mag = rup.getMag();
					int magIndex = xyz.getYIndex(mag);
					double prob = rup.getProbability();
					for (int m=0; m<magIndex && m<numMag; m++) {
						probsList.get(m).add(prob);
					}
				}
			}
//			SummedMagFreqDist mfd = ERF_Calculator.getParticipationMagFreqDistInRegion(
//					erf, reg, minMag, numMag, deltaMag, true, rupInRegionCache);
//			ERF_Calculator.get
			System.out.println("Duration: "+myDuration);
			for (int y=0; y<xyz.getNumY(); y++) {
				probs = probsList.get(y);
				totOneMinus = 1;
				for (double prob : probs) {
					totOneMinus *= (1-prob);
				}
				totProb = 1 - totOneMinus;
				u2XYZ.set(x, y, totProb*100d);
				System.out.println("\tM "+xyz.getY(y)+"+ Prob: "+(float)(totProb*100d)+" %");
			}
		}
		spec = new XYZPlotSpec(u2XYZ, cpt, "UCERF2 SoCal Time Dep Probability", "Duration", "Min Mag", "Probability (%)");
		new XYZPlotWindow(spec);
	}

}
