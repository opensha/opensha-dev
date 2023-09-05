package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.BorderType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import scratch.UCERF3.analysis.TablesAndPlotsGen;

public class LengthDistCompPlot {
	
	private enum PlotType {
		MEAN_ONLY,
		SEG_MODEL,
		B_VAL
	};

	public static void main(String[] args) throws IOException {
		File invDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
//				+ "2023_01_17-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//				+ "2023_02_28-nshm23_branches-seg_limit_max_length-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-NoRed-EvenFitPaleo-ThreshAvgIterRelGR");
//				+ "2023_02_25-nshm23_branches-seg_limit_max_length_600-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-NoRed-EvenFitPaleo-ThreshAvgIterRelGR");
//				+ "2023_03_01-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
//				+ "2023_02_21-nshm23_branches-NSHM23_v2-CoulombRupSet-NSHM23_Avg-TotNuclRate-NoRed-EvenFitPaleo-ThreshAvgIterRelGR");
//				+ "2023_03_31-nshm23_branches-mod_scaling-NSHM23_v2-CoulombRupSet-AVERAGE-TotNuclRate-NoRed-EvenFitPaleo-ThreshAvgIterRelGR");
//				+ "2023_03_31-nshm23_branches-orig_draft_scaling-NSHM23_v2-CoulombRupSet-AVERAGE-TotNuclRate-NoRed-EvenFitPaleo-ThreshAvgIterRelGR");
//				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		File miscPlotsDir = new File(invDir, "misc_plots");
		Preconditions.checkState(miscPlotsDir.exists() || miscPlotsDir.mkdir());
		
//		boolean onlyCA = false;
//		boolean noCA = false;
//		File outputDir = new File(miscPlotsDir, "wells_2013_length_dists");
		
//		boolean onlyCA = false;
//		boolean noCA = true;
//		File outputDir = new File(miscPlotsDir, "wells_2013_length_dists_noCA");
		
		boolean onlyCA = true;
		boolean noCA = false;
		File outputDir = new File(miscPlotsDir, "wells_2013_length_dists_onlyCA");
		
		Preconditions.checkState(!onlyCA || !noCA);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		String prefix = "wells_2013_length_dist_compare";
		
		File resultsFile = new File(invDir, "results.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile);
		LogicTree<?> tree = slt.getLogicTree();
//		tree = tree.sample(100, false);
		
		FaultSystemRupSet rupSet0 = slt.forBranch(tree.getBranch(0)).getRupSet();
		double[] lengths = rupSet0.getLengthForAllRups();
		
		Range linearYRange = new Range(0d, 1d);
		Range xRangeZoom = new Range(0d, 500d);
		Range logYRangeZoom = new Range(1e-3, 1e0);
		Range xRangeFull = new Range(0d, 1000d);
		Range logYRangeFull = new Range(1e-5, 1e0);
		HistogramFunction ref = HistogramFunction.getEncompassingHistogram(
				xRangeFull.getLowerBound(), xRangeFull.getUpperBound(), 50);
		HistogramFunction cmlRef = getInvCml(ref);
		
		boolean[] includeRups = null;
		if (noCA || onlyCA) {
			XY_DataSet[] caOutlines = PoliticalBoundariesData.loadCAOutlines();
			Region[] caRegions = new Region[caOutlines.length];
			for (int i=0; i<caOutlines.length; i++) {
				LocationList outline = new LocationList();
				for (Point2D pt : caOutlines[i])
					outline.add(new Location(pt.getY(), pt.getX()));
				caRegions[i] = new Region(outline, BorderType.MERCATOR_LINEAR);
			}
			Arrays.sort(caRegions, new Comparator<Region>() {

				@Override
				public int compare(Region o1, Region o2) {
					return -Double.compare(o1.getExtent(), o2.getExtent());
				}
			});
			caRegions = new Region[] {caRegions[0]};
			
			boolean[] caSects = new boolean[rupSet0.getNumSections()];
			for (int s=0; s<caSects.length; s++) {
				boolean skip = false;
				locLoop:
					for (Location loc : rupSet0.getFaultSectionData(s).getFaultTrace()) {
						for (Region reg : caRegions) {
							if (reg.contains(loc)) {
								skip = true;
								break locLoop;
							}
						}
					}
				caSects[s] = skip;
			}
			if (onlyCA) {
				// include only ruptures that touch CA (including partial)
				includeRups = new boolean[rupSet0.getNumRuptures()];
				for (int r=0; r<includeRups.length; r++) {
					boolean include = false;
					for (int s : rupSet0.getSectionsIndicesForRup(r)) {
						if (caSects[s]) {
							include = true;
							break;
						}
					}
					includeRups[r] = include;
				}
			} else {
				// exclude anything that touches CA (also excludes partial)
				includeRups = new boolean[rupSet0.getNumRuptures()];
				for (int r=0; r<includeRups.length; r++) {
					boolean include = true;
					for (int s : rupSet0.getSectionsIndicesForRup(r)) {
						if (caSects[s]) {
							include = false;
							break;
						}
					}
					includeRups[r] = include;
				}
			}
		}
		
		// NSHM18
		Path erfPath;
		if (noCA)
			erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.3.0-noCA/");
		else if (onlyCA)
			erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.3.0-onlyCA/");
		else
			erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.3.0/");
		NshmErf nshmErf = new NshmErf(erfPath,
				Set.of(TectonicRegionType.ACTIVE_SHALLOW), IncludeBackgroundOption.EXCLUDE);
		nshmErf.getTimeSpan().setDuration(1d);
		nshmErf.updateForecast();
		HistogramFunction nshm18 = new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta());
		int numFailed = 0;
		double rateFailed = 0d;
		double totRate = 0d;
		for (ProbEqkSource source : nshmErf) {
			for (ProbEqkRupture rup : source) {
				double rate = rup.getMeanAnnualRate(1d);
				Preconditions.checkState(rate > 0d);
				totRate += rate;
				double len;
				try {
					len = rup.getRuptureSurface().getAveLength();
					Preconditions.checkState(len > 0d);
					nshm18.add(nshm18.getClosestXIndex(len), rate);
				} catch (Exception e) {
					numFailed++;
					rateFailed += rate;
				}
			}
		}
		if (numFailed > 0)
			System.err.println("WARNING: length failed for "+numFailed+" ruptures, rate="+(float)rateFailed+" ("+
					new DecimalFormat("0.00%").format(rateFailed/totRate)+")");
		nshm18.scale(1d/nshm18.calcSumOfY_Vals());
		HistogramFunction nshm18cml = getInvCml(nshm18);

		Map<LogicTreeNode, HistogramFunction> meanHistMap = new HashMap<>();
		Map<LogicTreeNode, HistogramFunction> meanCmlHistMap = new HashMap<>();
		Map<LogicTreeNode, Double> meanHistWeightsMap = new HashMap<>();
		// full mean
		meanHistMap.put(null, new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta()));
		meanCmlHistMap.put(null, new HistogramFunction(cmlRef.getMinX(), cmlRef.size(), cmlRef.getDelta()));
		meanHistWeightsMap.put(null, 0d);

		HistogramFunction lowerHist = new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta());
		HistogramFunction upperHist = new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta());
		HistogramFunction lowerCmlHist = new HistogramFunction(cmlRef.getMinX(), cmlRef.size(), cmlRef.getDelta());
		HistogramFunction upperCmlHist = new HistogramFunction(cmlRef.getMinX(), cmlRef.size(), cmlRef.getDelta());
		for (int i=0; i<lowerHist.size(); i++) {
			lowerHist.set(i, Double.POSITIVE_INFINITY);
			lowerCmlHist.set(i, Double.POSITIVE_INFINITY);
		}
		for (int index=0; index<tree.size(); index++) {
			LogicTreeBranch<?> branch = tree.getBranch(index);
			System.out.println("Processing branch "+index+"/"+tree.size());
			
			NSHM23_SegmentationModels segModel = branch.getValue(NSHM23_SegmentationModels.class);
			if (segModel != null && !meanHistMap.containsKey(segModel)) {
				meanHistMap.put(segModel, new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta()));
				meanCmlHistMap.put(segModel, new HistogramFunction(cmlRef.getMinX(), cmlRef.size(), cmlRef.getDelta()));
				meanHistWeightsMap.put(segModel, 0d);
			}
			SupraSeisBValues bVal = branch.getValue(SupraSeisBValues.class);
			if (bVal != null && !meanHistMap.containsKey(bVal)) {
				meanHistMap.put(bVal, new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta()));
				meanCmlHistMap.put(bVal, new HistogramFunction(cmlRef.getMinX(), cmlRef.size(), cmlRef.getDelta()));
				meanHistWeightsMap.put(bVal, 0d);
			}
			
			double weight = tree.getBranchWeight(branch);
			
			double[] rates = slt.loadRatesForBranch(branch);
			Preconditions.checkState(rates.length == lengths.length);
			
			HistogramFunction myHist = new HistogramFunction(ref.getMinX(), ref.size(), ref.getDelta());
			for (int r=0; r<rates.length; r++) {
				if (rates[r] > 0 && (includeRups == null || includeRups[r])) {
					double len = lengths[r]*1e-3; // m -> km
					if (len > myHist.getMaxX()+0.5*myHist.getDelta())
						continue;
					int lenIndex = myHist.getClosestXIndex(len);
					myHist.add(lenIndex, rates[r]);
				}
			}
			myHist.normalizeBySumOfY_Vals();
			
			HistogramFunction myCmlHist = getInvCml(myHist);
			
			for (LogicTreeNode node : meanHistMap.keySet()) {
				if (node != null && !branch.hasValue(node))
					continue;
				HistogramFunction meanHist = meanHistMap.get(node);
				HistogramFunction meanCmlHist = meanCmlHistMap.get(node);
				meanHistWeightsMap.put(node, meanHistWeightsMap.get(node)+weight);
				for (int i=0; i<myHist.size(); i++) {
					double val = myHist.getY(i);
					meanHist.add(i, weight*val);
					if (val > upperHist.getY(i))
						upperHist.set(i, val);
					if (val < lowerHist.getY(i))
						lowerHist.set(i, val);
					double cmlVal = myCmlHist.getY(i);
					meanCmlHist.add(i, weight*cmlVal);
					if (cmlVal > upperCmlHist.getY(i))
						upperCmlHist.set(i, cmlVal);
					if (cmlVal < lowerCmlHist.getY(i))
						lowerCmlHist.set(i, cmlVal);
				}
			}
		}
		HistogramFunction data = TablesAndPlotsGen.loadSurfaceRupData();
		HistogramFunction cmlData = getInvCml(data);
		
		for (int i=0; i<lowerHist.size(); i++) {
			if (Double.isInfinite(lowerHist.getY(i)))
				lowerHist.set(i, 0d);
			if (Double.isInfinite(lowerCmlHist.getY(i)))
				lowerCmlHist.set(i, 0d);
		}
		
		for (LogicTreeNode node : meanHistMap.keySet()) {
			meanHistMap.get(node).scale(1d/meanHistWeightsMap.get(node));
			meanCmlHistMap.get(node).scale(1d/meanHistWeightsMap.get(node));
		}

		HistogramFunction meanHist = meanHistMap.get(null);
		HistogramFunction meanCmlHist = meanCmlHistMap.get(null);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		boolean[] falseTrue = { false, true };
		boolean[] falseOnly = { false };
		boolean[] trueOnly = { true };
		
		for (PlotType type : PlotType.values()) {
			for (boolean cml : falseTrue) {
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				HistogramFunction dataHist = cml ? cmlData : data;
				dataHist.setName("Wells (2013), Observed");
				funcs.add(dataHist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.LIGHT_GRAY));
				
				HistogramFunction nshm18Hist = cml ? nshm18cml : nshm18; 
				nshm18Hist.setName("NSHM18 Active Crustal");
				funcs.add(nshm18Hist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 4f, Color.DARK_GRAY));
				
				HistogramFunction modelHist = cml ? meanCmlHist : meanHist;
				modelHist.setName("NSHN23 WUS, Mean");
				funcs.add(modelHist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, Color.BLACK));
				
				HistogramFunction modelLowerHist = cml ? lowerCmlHist : lowerHist;
				modelLowerHist.setName("NSHN23 WUS, Min & Max");
				funcs.add(modelLowerHist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				
				HistogramFunction modelUpperHist = cml ? upperCmlHist : upperHist;
				modelUpperHist.setName(null);
				funcs.add(modelUpperHist);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				
				List<LogicTreeNode> nodes = new ArrayList<>();
				List<HistogramFunction> nodeHists = new ArrayList<>();
				String myPrefix = prefix;
				if (type == PlotType.B_VAL) {
					myPrefix += "_bVal";
					for (SupraSeisBValues bVal : SupraSeisBValues.values()) {
						HistogramFunction nodeHist = cml ? meanCmlHistMap.get(bVal) : meanHistMap.get(bVal);
						if (nodeHist != null) {
							nodes.add(bVal);
							nodeHists.add(nodeHist);
						}
					}
					if (nodes.isEmpty())
						continue;
				} else if (type == PlotType.SEG_MODEL) {
					myPrefix += "_segModel";
					for (NSHM23_SegmentationModels segModel : NSHM23_SegmentationModels.values()) {
						HistogramFunction nodeHist = cml ? meanCmlHistMap.get(segModel) : meanHistMap.get(segModel);
						if (nodeHist != null) {
							nodes.add(segModel);
							nodeHists.add(nodeHist);
						}
					}
					if (nodes.isEmpty())
						continue;
				}
				
				if (nodes.size() > 1) {
					CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, nodes.size()-1);
					
					for (int n=0; n<nodes.size(); n++) {
						Color color = cpt.getColor((float)n);
						LogicTreeNode node = nodes.get(n);
						HistogramFunction hist = nodeHists.get(n);
						
						hist.setName(node.getShortName());
						funcs.add(hist);
						chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color));
					}
					
					// copy it on top
					EvenlyDiscretizedFunc meanHistCopy = cml ? meanCmlHist.deepClone() : meanHist.deepClone();
					meanHistCopy.setName(null);
					funcs.add(meanHistCopy);
					chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, new Color(0, 0, 0, 127)));
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, "Rupture Length Distribution", "Rupture Length (km)",
						cml ? "Cumulative Fraction of Earthquakes" : "Fraction of Earthquakes");
				spec.setLegendInset(true);
				
				boolean[] fullXs = type == PlotType.MEAN_ONLY ? falseTrue : trueOnly;
				
				for (boolean fullX : fullXs) {
					boolean[] logYs = fullX || type != PlotType.MEAN_ONLY ? trueOnly : falseTrue;
					for (boolean logY : logYs) {
						String plotPrefix = myPrefix;
						Range xRange, yRange;
						
						if (fullX) {
							plotPrefix += "_wideX";
							xRange = xRangeFull;
						} else {
							xRange = xRangeZoom;
						}
						
						if (logY) {
							plotPrefix += "_log";
							if (fullX)
								yRange = logYRangeFull;
							else
								yRange = logYRangeZoom;
						} else {
							yRange = linearYRange;
						}
						
						if (cml)
							plotPrefix += "_cml";
						
						gp.drawGraphPanel(spec, false, logY, xRange, yRange);
						
						PlotUtils.writePlots(outputDir, plotPrefix, gp, 1000, 800, true, true, true);
					}
				}
			}
		}
	}
	
	private static HistogramFunction getInvCml(HistogramFunction hist) {
		HistogramFunction ret = new HistogramFunction(hist.getMinX(), hist.size(), hist.getDelta());
		
		double sumY = 0d;
		for (int i=hist.size(); --i>=0;) {
			sumY += hist.getY(i);
			ret.set(i, sumY);
		}
		ret.setName(hist.getName());
		
		return ret;
	}

}
