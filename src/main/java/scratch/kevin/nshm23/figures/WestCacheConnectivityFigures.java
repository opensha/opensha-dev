package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class WestCacheConnectivityFigures {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/west_cache_conn_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_FaultModels fm = NSHM23_FaultModels.NSHM23_v2;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		NSHM23_SingleStates state = NSHM23_SingleStates.UT;
		SupraSeisBValues bVal = SupraSeisBValues.B_0p5;
		NSHM23_SegmentationModels segModel = NSHM23_SegmentationModels.MID;
		SegmentationMFD_Adjustment adj = SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG;
//		SegmentationMFD_Adjustment adj = SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG;
		
		List<LogicTreeLevel<? extends LogicTreeNode>> levels = new ArrayList<>(NSHM23_LogicTreeBranch.levelsOnFault);
		levels.add(NSHM23_LogicTreeBranch.SINGLE_STATES);
		
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(levels);
		for (LogicTreeNode node : NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT)
			branch.setValue(node);
		branch.setValue(fm);
		branch.setValue(dm);
		branch.setValue(state);
		branch.setValue(bVal);
		branch.setValue(segModel);
		
		NSHM23_InvConfigFactory factory = new NSHM23_InvConfigFactory();
		factory.setCacheDir(new File("/data/kevin/nshm23/rup_sets/cache"));
		
		FaultSystemRupSet rupSet = factory.buildRuptureSet(branch, FaultSysTools.defaultNumThreads());
		
		SupraSeisBValInversionTargetMFDs.Builder mfdBuilder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bVal.bValue);
		
		mfdBuilder.applyDefModelUncertainties(false);
		
		int parentID = 2714;
		int sectID = 362;
//		int sectID = -1;
		
		JumpProbabilityCalc segCalc = segModel.getModel(rupSet, branch);
		
		mfdBuilder.sparseGR(false);
		SupraSeisBValInversionTargetMFDs pureGRs = mfdBuilder.build();
		mfdBuilder.sparseGR(true);
		SupraSeisBValInversionTargetMFDs unadjustedMFDs = mfdBuilder.build();
//		RelGRWorstJumpProb.D = true;
		mfdBuilder.adjustTargetsForData(adj.getAdjustment(segCalc));
		SupraSeisBValInversionTargetMFDs adjustedMFDs = mfdBuilder.build();
		
		Range xRange = new Range(6.5, 8d);
		Range yRange = sectID >= 0 ? new Range(1e-7, 3e-5) : new Range(1e-7, 1e-4);

		IncrementalMagFreqDist pureGR = null;
		IncrementalMagFreqDist unadjustedMFD = null;
		IncrementalMagFreqDist adjustedMFD = null;
		
		for (int s=0; s<rupSet.getNumSections(); s++) {
			FaultSection sect = rupSet.getFaultSectionData(s);
			if (sectID >= 0 && s != sectID)
				continue;
			if (sect.getParentSectionId() == parentID) {
				System.out.println(s+". "+sect.getSectionName());
				IncrementalMagFreqDist grMFD = pureGRs.getOnFaultSupraSeisNucleationMFDs().get(s);
				IncrementalMagFreqDist unMFD = unadjustedMFDs.getOnFaultSupraSeisNucleationMFDs().get(s);
				IncrementalMagFreqDist adMFD = adjustedMFDs.getOnFaultSupraSeisNucleationMFDs().get(s);
				
				if (unadjustedMFD == null) {
					pureGR = new IncrementalMagFreqDist(unMFD.getMinX(), unMFD.size(), unMFD.getDelta());
					unadjustedMFD = new IncrementalMagFreqDist(unMFD.getMinX(), unMFD.size(), unMFD.getDelta());
					adjustedMFD = new IncrementalMagFreqDist(unMFD.getMinX(), unMFD.size(), unMFD.getDelta());
				} else {
					Preconditions.checkState(unadjustedMFD.getMinX() == unMFD.getMinX());
					Preconditions.checkState(unadjustedMFD.size() == unMFD.size());
				}
				
				for (int i=0; i<unMFD.size(); i++) {
					pureGR.add(i, grMFD.getY(i));
					unadjustedMFD.add(i, unMFD.getY(i));
					adjustedMFD.add(i, adMFD.getY(i));
				}
			}
		}
		
		List<XY_DataSet> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		Color grColor = Color.CYAN;
		grColor = new Color(grColor.getRed(), grColor.getGreen(), grColor.getBlue(), 100);
		IncrementalMagFreqDist redistBins = new IncrementalMagFreqDist(pureGR.getMinX(), pureGR.size(), pureGR.getDelta());
		for (int i=0; i<pureGR.size(); i++)
			if (unadjustedMFD.getY(i) == 0d)
				redistBins.set(i, pureGR.getY(i));
		if (redistBins.calcSumOfY_Vals() > 0d) {
			redistBins.setName("Redistributed Pure G-R");
			funcs.add(redistBins);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, grColor));
		}
		
		Color origColor = Color.LIGHT_GRAY;
		unadjustedMFD.setName("Unadjusted");
		funcs.add(unadjustedMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, origColor));
		
		adjustedMFD.setName("Adjusted");
		funcs.add(adjustedMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.GREEN.darker()));
		
		for (int i=0; i<unadjustedMFD.size(); i++) {
			double origY = unadjustedMFD.getY(i);
			double adjY = adjustedMFD.getY(i);
			
			if ((float)adjY >= (float)origY) {
				// plot it on top
				double x = adjustedMFD.getX(i);
				double x1 = x - 0.4*adjustedMFD.getDelta();
				double x2 = x + 0.4*adjustedMFD.getDelta();
				ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
				func.set(x1, origY);
				func.set(x2, origY);
				funcs.add(func);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, origColor));
			}
		}
		
		double[] magBinnedPassthroughRates = new double[adjustedMFD.size()];
		for (int i=0; i<magBinnedPassthroughRates.length; i++)
			magBinnedPassthroughRates[i] = Double.NaN;
		double[] sectMaxPassthroughRates = new double[rupSet.getNumSections()];
		double[] sectMinConnectedMags = new double[rupSet.getNumSections()];
		for (int s=0; s<sectMinConnectedMags.length; s++) {
			sectMaxPassthroughRates[s] = Double.NaN;
			sectMinConnectedMags[s] = Double.NaN;
		}
		ClusterRupture largestRup = null;
		double maxMag = 0d;
		
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		Map<Jump, Double> jumpProbs = ThresholdAveragingSectNuclMFD_Estimator.calcAverageJumpProbs(cRups, segCalc);
		List<Integer> rupIndexes = sectID < 0 ? rupSet.getRupturesForParentSection(parentID) : rupSet.getRupturesForSection(sectID);
		for (int rupIndex : rupIndexes) {
			ClusterRupture cRup = cRups.get(rupIndex);
			double mag = rupSet.getMagForRup(rupIndex);
			boolean singleFault = cRup.getTotalNumClusters() == 1;
			double passthroughRate = 1d;
			for (Jump jump : cRup.getJumpsIterable())
				passthroughRate = Math.min(passthroughRate, jumpProbs.get(jump));
			if (!singleFault) {
				int magBin = adjustedMFD.getClosestXIndex(mag);
				if (Double.isNaN(magBinnedPassthroughRates[magBin]))
					magBinnedPassthroughRates[magBin] = passthroughRate;
				else
					magBinnedPassthroughRates[magBin] = Math.max(magBinnedPassthroughRates[magBin], passthroughRate);
				if (magBin == adjustedMFD.getClosestXIndex(7.55))
					System.out.println(rupIndex+". M"+(float)mag+", passthrough="+(float)passthroughRate+": "+cRup);
			}
			if (mag > maxMag) {
				largestRup = cRup;
				maxMag = mag;
			}
			for (FaultSubsectionCluster cluster : cRup.getClustersIterable()) {
				for (FaultSection sect : cluster.subSects) {
					int s = sect.getSectionId();
					if (Double.isNaN(sectMinConnectedMags[s]))
						sectMinConnectedMags[s] = mag;
					else
						sectMinConnectedMags[s] = Math.min(sectMinConnectedMags[s], mag);
					if (!singleFault) {
						if (Double.isNaN(sectMaxPassthroughRates[s]))
							sectMaxPassthroughRates[s] = passthroughRate;
						else
							sectMaxPassthroughRates[s] = Math.max(sectMaxPassthroughRates[s], passthroughRate);
					}
				}
			}
		}
		
		List<FaultSection> highlightSects = new ArrayList<>();
		for (FaultSection sect : rupSet.getFaultSectionDataList())
			if (sect.getParentSectionId() == parentID && (sectID < 0 || sectID == sect.getSectionId()))
				highlightSects.add(sect);
		String faultName = sectID < 0 ? highlightSects.get(0).getParentSectionName() : rupSet.getFaultSectionData(sectID).getSectionName();
		
		// MFD plot
		// identify contiguous ranges
		List<Range> magRanges = new ArrayList<>();
		List<Double> rangePassthroughs = new ArrayList<>();
		// first find single fault mags
		double minSingle = Double.MAX_VALUE;
		double maxSingle = 0d;
		double halfDelta = adjustedMFD.getDelta()*0.5;
		for (int i=0; i<adjustedMFD.size(); i++) {
			if (Double.isNaN(magBinnedPassthroughRates[i]) && adjustedMFD.getY(i) > 0) {
				double mag = adjustedMFD.getX(i);
				minSingle = Math.min(minSingle, mag-halfDelta);
				maxSingle = Math.max(maxSingle, mag+halfDelta);
			}
		}
		System.out.println("Single fault mags: "+minSingle+", "+maxSingle);
		magRanges.add(new Range(minSingle, maxSingle));
		rangePassthroughs.add(null);
		// now contiguous multi faults
		for (int i=0; i<adjustedMFD.size();) {
			if (adjustedMFD.getY(i) == 0d || Double.isNaN(magBinnedPassthroughRates[i])) {
				i++;
				continue;
			}
			// multi fault magnitude
			double passthrough = magBinnedPassthroughRates[i];
			int endIndex = i;
			for (int j=i+1; j<magBinnedPassthroughRates.length; j++) {
				if ((float)passthrough == (float)magBinnedPassthroughRates[j])
					endIndex = j;
				else
					break;
			}
			magRanges.add(new Range(adjustedMFD.getX(i)-halfDelta, adjustedMFD.getX(endIndex)+halfDelta));
			rangePassthroughs.add(passthrough);
			i = endIndex+1;
		}
		
		HashSet<Float> uniqueVertXs = new HashSet<>();
		List<XYAnnotation> anns = new ArrayList<>();
		DecimalFormat passthroughDF = new DecimalFormat("0.0#");
		
		double logMaxY = Math.log10(yRange.getUpperBound());
		double logMinY = Math.log10(yRange.getLowerBound());
//		double annY = Math.pow(10, logMinY + 0.9*(logMaxY-logMinY));
//		double annLineY1 = Math.pow(10, logMinY + 0.85*(logMaxY-logMinY));
//		double annLineY2 = Math.pow(10, logMinY + 0.95*(logMaxY-logMinY));
		double annY = Math.pow(10, logMinY + 0.95*(logMaxY-logMinY));
		double annLineY1 = yRange.getLowerBound();
		double annLineY2 = yRange.getUpperBound();
		
		for (int i=0; i<magRanges.size(); i++) {
			Range range = magRanges.get(i);
			Double passthrough = rangePassthroughs.get(i);
			
			String text;
			if (passthrough == null)
				text = "Single-Fault";
			else
				text = "F≤"+passthroughDF.format(passthrough);
			
			// see if we should add vertical lines
			for (double x : new double[] { range.getLowerBound(), range.getUpperBound() }) {
				if (uniqueVertXs.contains((float)x))
					continue;
				uniqueVertXs.add((float)x);
				DefaultXY_DataSet line = new DefaultXY_DataSet();
				line.set(x, annLineY1);
				line.set(x, annLineY2);
				
				funcs.add(line);
				chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.BLACK));
			}
			
			XYTextAnnotation ann = new XYTextAnnotation(text, range.getCentralValue(), annY);
			ann.setTextAnchor(TextAnchor.CENTER);
			ann.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 18));
			anns.add(ann);
		}
		PlotSpec spec = new PlotSpec(funcs, chars, faultName, "Magnitude",
				"Incremental Nucleation Rate (/yr)");
		spec.setPlotAnnotations(anns);
		spec.setLegendVisible(true);
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "mfds", gp, 800, 700, true, true, false);
		
		// map plots
		double minLat = Double.POSITIVE_INFINITY;
		double minLon = Double.POSITIVE_INFINITY;
		double maxLat = Double.NEGATIVE_INFINITY;
		double maxLon = Double.NEGATIVE_INFINITY;
		for (FaultSubsectionCluster cluster : largestRup.getClustersIterable()) {
			for (FaultSection sect : cluster.subSects) {
				for (Location loc : sect.getFaultSurface(1d).getPerimeter()) {
					minLat = Math.min(minLat, loc.lat);
					minLon = Math.min(minLon, loc.lon);
					maxLat = Math.max(maxLat, loc.lat);
					maxLon = Math.max(maxLon, loc.lon);
				}
			}
		}
		
		Location topRight = new Location(maxLat, maxLon);
		Location botLeft = new Location(minLat, minLon);
		topRight = LocationUtils.location(topRight, Math.PI/4d, 50d);
		botLeft = LocationUtils.location(botLeft, 5d*Math.PI/4d, 50d);
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(rupSet, new Region(botLeft, topRight));
		
		CPT magCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(xRange.getLowerBound(), xRange.getUpperBound());
		CPT passthroughCPT = GMT_CPT_Files.GMT_COPPER.instance().rescale(0d, 1d);
		// trim it to not quite be so dark at the low end
		CPT modPassthroughCPT = new CPT();
		for (double x=0.21; (float)x<=1f; x+=0.01) {
			float start = (float)(x-0.01);
			float end = (float)x;
			modPassthroughCPT.add(new CPTVal(start, passthroughCPT.getColor(start), end, passthroughCPT.getColor(end)));
		}
		modPassthroughCPT.setBelowMinColor(modPassthroughCPT.getMinColor());
		modPassthroughCPT.setAboveMaxColor(modPassthroughCPT.getMaxColor());
		passthroughCPT = modPassthroughCPT.rescale(0d, 1d);
		
		mapMaker.setSectHighlights(highlightSects, new PlotCurveCharacterstics(PlotLineType.SOLID, 6f, Color.BLACK));
		mapMaker.setSkipNaNs(true);
		mapMaker.setDefaultPlotWidth(650);
		
		mapMaker.plotSectScalars(sectMinConnectedMags, magCPT, "Minimum Connected Magnitude");
		mapMaker.plot(outputDir, "map_min_mag", " ");
		PlotSpec magSpec = mapMaker.buildPlot(" ");
		
		mapMaker.plotSectScalars(sectMaxPassthroughRates, passthroughCPT, "Controlling Passthrough Rate");
		mapMaker.plot(outputDir, "map_passthrough", " ");
		PlotSpec passthroughSpec = mapMaker.buildPlot(" ");
		
		Range mapXRange = new Range(botLeft.getLongitude(), topRight.getLongitude());
		Range mapYRange = new Range(botLeft.getLatitude(), topRight.getLatitude());
		gp.drawGraphPanel(List.of(magSpec, passthroughSpec), false, false, List.of(mapXRange, mapXRange), List.of(mapYRange));
		PlotUtils.setXTick(gp, 0.5);
		PlotUtils.setYTick(gp, 0.5);
		
		PlotUtils.writePlots(outputDir, "map_combined", gp, 800, true, true, true, false);
	}

}
