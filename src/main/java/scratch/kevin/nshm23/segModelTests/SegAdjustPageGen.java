package scratch.kevin.nshm23.segModelTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.ConnectivityClusters;
import org.opensha.sha.earthquake.faultSysSolution.modules.RuptureSubSetMappings;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.Jump;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.Shaw07JumpDistProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.ConnectivityCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RupSetMapMaker;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.RuptureTreeNavigator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.Builder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.ThresholdAveragingSectNuclMFD_Estimator.RelGRWorstJumpProb;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SectNucleationMFD_Estimator;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class SegAdjustPageGen {

	public static void main(String[] args) throws IOException {
		File markdownDir = new File("/home/kevin/markdown/nshm23-misc/segmentation-adj");
		
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2022_08_18-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
//		int subSectIndex = 1910; // Hogsback
//		int subSectIndex = 1907; // Hoback
//		int subSectIndex = 692; // Cemetery
//		int subSectIndex = 1903; // Hilton Creek
//		int subSectIndex = 4378; // Sangre de Cristo (south)
		int subSectIndex = 4366; // Sangre de Cristo (San Luis)
		
		FaultSection subSect = rupSet.getFaultSectionData(subSectIndex);
		String sectName = subSect.getSectionName();
		
		String dirName = sectName.replaceAll("\\W+", "_")+"_"+subSectIndex;
		while (dirName.contains("__"))
			dirName = dirName.replaceAll("__", "_");
		while (dirName.startsWith("_"))
			dirName = dirName.substring(1);
		
		File outputDir = new File(markdownDir, dirName);
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File resourcesDir = new File(outputDir, "resources");
		Preconditions.checkState(resourcesDir.exists() || resourcesDir.mkdir());
		
		SupraSeisBValues[] bVals = {SupraSeisBValues.B_0p0, SupraSeisBValues.B_0p5, SupraSeisBValues.B_1p0};
		
		List<SegmentationMFD_Adjustment> segAdjs = new ArrayList<>();
		List<PlotCurveCharacterstics> segChars = new ArrayList<>();
		
		RelGRWorstJumpProb.D = true;
		
		segAdjs.add(SegmentationMFD_Adjustment.REL_GR_THRESHOLD_AVG);
		segChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.RED.darker()));
		
		segAdjs.add(SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG_MATCH_STRICT);
		segChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE.darker()));
		
//		segAdjs.add(SegmentationMFD_Adjustment.STRICT_SEG_REPRODUCE_THROUGH_INVERSION);
//		segChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.ORANGE.darker().darker()));
		
		segAdjs.add(SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG);
		segChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GREEN.darker()));
		
//		segAdjs.add(SegmentationMFD_Adjustment.JUMP_PROB_THRESHOLD_AVG_MATCH_STRICT);
//		segChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE.darker()));
		
		// will be used for no adjustment
		segChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		// will only be used if the G-R dist has holes
		segChars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 3f, Color.GRAY));
		
		JumpProbabilityCalc[] segModels = {
			Shaw07JumpDistProb.forHorzOffset(1d, 3d, 2d),
		};
		
		List<Integer> allSectRups = rupSet.getRupturesForSection(subSectIndex);
		HashSet<FaultSection> corupSects = new HashSet<>();
		for (int rupIndex : allSectRups)
			for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex))
				corupSects.add(sect);
		
		Region corupReg = RupSetMapMaker.buildBufferedRegion(corupSects, 30d, true);
		
		GeographicMapMaker mapMaker = new RupSetMapMaker(rupSet, corupReg);
		
		List<String> lines = new ArrayList<>();
		
		lines.add("# "+sectName+" Segmentation Adjustments");
		lines.add("");
		
		int tocIndex = lines.size();
		String topLink = "_[(top)](#table-of-contents)_";
		
		ClusterRuptures cRups = rupSet.requireModule(ClusterRuptures.class);
		double[] sectCmlDists = new double[rupSet.getNumSections()];
		double[] sectMaxDists = new double[rupSet.getNumSections()];
		for (int s=0; s<sectCmlDists.length; s++) {
			sectCmlDists[s] = Double.NaN;
			sectMaxDists[s] = Double.NaN;
		}
		sectCmlDists[subSectIndex] = 0d;
		sectMaxDists[subSectIndex] = 0d;
		Map<Jump, List<Jump>> allJumpsMap = new HashMap<>();
		for (int rupIndex : allSectRups) {
			ClusterRupture rup = cRups.get(rupIndex);
			RuptureTreeNavigator nav = rup.getTreeNavigator();
			
			distTraverse(nav, subSect, true, 0d, 0d, sectCmlDists, sectMaxDists);
			distTraverse(nav, subSect, false, 0d, 0d, sectCmlDists, sectMaxDists);
			
			for (Jump jump : rup.getJumpsIterable()) {
				if (jump.fromSection.getSectionId() > jump.toSection.getSectionId())
					jump = jump.reverse();
				List<Jump> allJumps = allJumpsMap.get(jump);
				if (allJumps == null) {
					allJumps = new ArrayList<>();
					allJumpsMap.put(jump, allJumps);
				}
				allJumps.add(jump);
			}
		}
		
		// find average distances and probabilities
		Map<Jump, Double> jumpDists = new HashMap<>();
		List<Map<Jump, Double>> modelJumpProbs = new ArrayList<>();
		for (int i=0; i<segModels.length; i++)
			modelJumpProbs.add(new HashMap<>());
		for (Jump jump : allJumpsMap.keySet()) {
			List<Jump> jumps = allJumpsMap.get(jump);
			QuickAvgTrack distTrack = new QuickAvgTrack();
			for (Jump oJump : jumps)
				distTrack.add(oJump.distance);
			jumpDists.put(jump, distTrack.getAverage());
			for (int i=0; i<segModels.length; i++) {
				QuickAvgTrack probTrack = new QuickAvgTrack();
				for (Jump oJump : jumps)
					probTrack.add(segModels[i].calcJumpProbability(null, oJump, false));
				modelJumpProbs.get(i).put(jump, probTrack.getAverage());
			}
		}
		
		List<Jump> jumpsSorted = ComparablePairing.getSortedData(jumpDists);
		
		Map<Jump, BitSet> jumpMagBins = new HashMap<>();
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(rupSet);
		// figure out which jumps use which mag bins
		for (int r=0; r<cRups.size(); r++) {
			int magBin = refMFD.getClosestXIndex(rupSet.getMagForRup(r));
			for (Jump jump : cRups.get(r).getJumpsIterable()) {
				if (jump.fromSection.getSectionId() > jump.toSection.getSectionId())
					jump = jump.reverse();
				BitSet bitSet = jumpMagBins.get(jump);
				if (bitSet == null) {
					bitSet = new BitSet(refMFD.size());
					jumpMagBins.put(jump, bitSet);
				}
				bitSet.set(magBin);
			}
		}
		
		Map<Integer, Double> parentCmlDists = new HashMap<>();
		Map<Integer, Double> parentMaxDists = new HashMap<>();
		Map<Integer, String> parentNames = new HashMap<>();
		List<FaultSection> mySects = new ArrayList<>();
		double maxSectDist = 0d;
		for (int s=0; s<rupSet.getNumSections(); s++) {
			if (Double.isFinite(sectCmlDists[s])) {
				maxSectDist = Math.max(maxSectDist, sectCmlDists[s]);
				FaultSection sect = rupSet.getFaultSectionData(s);
				int myParent = sect.getParentSectionId();
				if (myParent == subSect.getParentSectionId()) {
					mySects.add(sect);
				} else {
					if (parentCmlDists.containsKey(myParent)) {
						parentCmlDists.put(myParent, Math.min(sectCmlDists[s], parentCmlDists.get(myParent)));
						parentMaxDists.put(myParent, Math.min(sectMaxDists[s], parentMaxDists.get(myParent)));
					} else {
						parentCmlDists.put(myParent, sectCmlDists[s]);
						parentMaxDists.put(myParent, sectMaxDists[s]);
						parentNames.put(myParent, sect.getParentSectionName());
					}
				}
			}
		}
		
		CPT distCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(0d, Math.max(5d, 5d*Math.ceil(maxSectDist/5d)));
		
		mapMaker.setSkipNaNs(true);
		mapMaker.setSectHighlights(mySects , new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		mapMaker.plotJumpScalars(jumpDists, distCPT, null);
		
		mapMaker.plotSectScalars(sectCmlDists, distCPT, "Cumulative Jump Dist To (km)");
		mapMaker.plot(resourcesDir, "cml_dists_to", "Cumulative Dists");
		
		mapMaker.plotSectScalars(sectMaxDists, distCPT, "Max Indv. Jump Dist To (km)");
		mapMaker.plot(resourcesDir, "indv_dists_to", "Max Individual Dists");
		
		lines.add("## Neighbors");
		lines.add(topLink); lines.add("");
		
		TableBuilder table = MarkdownUtils.tableBuilder();
		
		DecimalFormat twoDigits = new DecimalFormat("0.00");
		
		table.addLine("Section Name", "Max Individual Jump Dist To (km)", "Cumulative Jump Dist To (km)");
		for (int parentID : ComparablePairing.getSortedData(parentCmlDists))
			table.addLine("**"+parentNames.get(parentID)+"**", twoDigits.format(parentMaxDists.get(parentID)),
					twoDigits.format(parentCmlDists.get(parentID)));
		
		lines.addAll(table.build());
		
		lines.add("![Map]("+resourcesDir.getName()+"/indv_dists_to.png)");
		lines.add("");
		lines.add("![Map]("+resourcesDir.getName()+"/cml_dists_to.png)");
		lines.add("");
		
		ConnectivityClusters connClusters = rupSet.getModule(ConnectivityClusters.class);
		if (connClusters == null)
			connClusters = ConnectivityClusters.build(rupSet);
		ConnectivityCluster myCluster = null;
		for (ConnectivityCluster cluster : connClusters) {
			if (cluster.containsSect(subSect)) {
				myCluster = cluster;
				break;
			}
		}
		Preconditions.checkNotNull(myCluster);
		FaultSystemRupSet subsetRupSet = rupSet.getForSectionSubSet(myCluster.getSectIDs());
		RuptureSubSetMappings subsetMappings = subsetRupSet.requireModule(RuptureSubSetMappings.class);
		
		DecimalFormat oDF = new DecimalFormat("0.##");
		for (SupraSeisBValues bValChoice : bVals) {
			
			LogicTreeBranch<LogicTreeNode> branch = rupSet.requireModule(LogicTreeBranch.class).copy();
			branch.setValue(bValChoice);
			double bVal = bValChoice.bValue;
			
			subsetRupSet.addModule(branch);
			
			lines.add("## b="+oDF.format(bVal)+" Section MFDs");
			lines.add(topLink); lines.add("");
			
			Builder builder = new SupraSeisBValInversionTargetMFDs.Builder(subsetRupSet, bVal);
			
			int subsetSectIndex = subsetMappings.getNewSectID(subSectIndex);
			
			builder.sparseGR(false);
			List<UncertainIncrMagFreqDist> fullGRs = builder.build().getOnFaultSupraSeisNucleationMFDs();
			IncrementalMagFreqDist fullGR = fullGRs.get(subsetSectIndex);
			fullGR.setName("Full G-R");
			builder.sparseGR(true);
			
			List<UncertainIncrMagFreqDist> allNoAdjGRs = builder.build().getOnFaultSupraSeisNucleationMFDs();
			IncrementalMagFreqDist noAdjGR = allNoAdjGRs.get(subsetSectIndex);
			int numSegs = 0;
			double prevY = 0d;
			for (int i=0; i<noAdjGR.size(); i++) {
				double y = noAdjGR.getY(i);
				if (y > 0d) {
					if (prevY == 0d)
						numSegs++;
				}
				prevY = y;
			}
			boolean isSparse = numSegs > 1;
			if (isSparse) 
				noAdjGR.setName("No Adjustment, Sparse G-R");
			else
				noAdjGR.setName("No Adjustment");
			
			// for converting nucleation to participation
			double[] nuclToParticScalars = calcNuclToParticScalars(rupSet, subSectIndex, refMFD);
			
			int segCount = 0;
			for (int m=0; m<segModels.length; m++) {
				JumpProbabilityCalc segModel = segModels[m];
				Map<Jump, Double> jumpProbs = modelJumpProbs.get(m);
				if (segModels.length > 1) {
					lines.add("### b="+oDF.format(bVal)+" Section MFDs, "+segModel.getName());
					lines.add(topLink); lines.add("");
				}
				
				int maxMagIndex = -1;
				int minMagIndex = Integer.MAX_VALUE;
				
				double minY = Double.POSITIVE_INFINITY;
				double maxY = 0d;
				List<List<? extends IncrementalMagFreqDist>> adjTargetsAllSectsList = new ArrayList<>();
				List<IncrementalMagFreqDist> adjTargetsList = new ArrayList<>();
				for (SegmentationMFD_Adjustment adj : segAdjs) {
					builder.clearTargetAdjustments();
					SectNucleationMFD_Estimator adjEst = adj.getAdjustment(segModel);
					if (adjEst != null)
						builder.adjustTargetsForData(adjEst);
					List<UncertainIncrMagFreqDist> targets = builder.build().getOnFaultSupraSeisNucleationMFDs();
					IncrementalMagFreqDist target = targets.get(subsetSectIndex);
					target.setName(adj.getName());
					adjTargetsList.add(target);
					adjTargetsAllSectsList.add(targets);
					
					for (int i=0; i<target.size(); i++) {
						double val = target.getY(i);
						if (val > 0d) {
							minMagIndex = Integer.min(minMagIndex, i);
							if (i > maxMagIndex) {
								maxMagIndex = i;
							}
							minY = Math.min(minY, val);
							maxY = Math.max(maxY, val);
						}
					}
				}
				
				GutenbergRichterMagFreqDist refGR = new GutenbergRichterMagFreqDist(bVal, 1d,
						refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
				double padDelta = 0.48*refMFD.getDelta();
				
				double minX = refMFD.getX(minMagIndex);
				double maxX = refMFD.getX(maxMagIndex);
				minX = 2d*Math.floor(minX/2d);
				maxX = 2d*Math.ceil(maxX/2d);
				Range xRange = new Range(minX, maxX);
				Range yRange = new Range(Math.pow(10, Math.floor(Math.log10(minY))),
						Math.pow(10, Math.ceil(Math.log10(maxY))));
				
				adjTargetsList.add(noAdjGR);
				adjTargetsAllSectsList.add(allNoAdjGRs);
				if (isSparse) {
					adjTargetsList.add(fullGR);
					adjTargetsAllSectsList.add(fullGRs);
				}
				
				List<DiscretizedFunc> funcs = new ArrayList<>();
				List<PlotCurveCharacterstics> chars = new ArrayList<>();
				
				for (int i=0; i<adjTargetsList.size(); i++) {
					IncrementalMagFreqDist mfd = adjTargetsList.get(i);
					PlotCurveCharacterstics pChar = segChars.get(i);
					
					DiscretizedFunc curFunc = null;
					boolean first = true;
					for (int xInd=minMagIndex; xInd<mfd.size(); xInd++) {
						double y = mfd.getY(xInd);
						
						if (y > 0d) {
							double x = mfd.getX(xInd);
							if (curFunc == null) {
								curFunc = new ArbitrarilyDiscretizedFunc();
								if (first)
									curFunc.setName(mfd.getName());
								first = false;
								funcs.add(curFunc);
								chars.add(pChar);
							}
							
							double scaleToGR = y/refGR.getY(xInd);
							curFunc.set(x-padDelta, scaleToGR * refGR.getInterpolatedY(x-padDelta));
							curFunc.set(x+padDelta, scaleToGR * refGR.getInterpolatedY(x+padDelta));
						} else if (curFunc != null) {
							curFunc = null;
						}
					}
				}
				
				PlotSpec spec = new PlotSpec(funcs, chars, sectName, "Magnitude", "Incremental Rate (1/yr)");
				spec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.setRenderingOrder(DatasetRenderingOrder.REVERSE);
				gp.drawGraphPanel(spec, false, true, xRange, yRange);
				
				String prefix = "b_"+(float)bVal;
				if (segModels.length > 1)
					prefix += "_seg_"+segCount;
				prefix += "_mfds";
				PlotUtils.writePlots(resourcesDir, prefix, gp, 1000, 800, true, true, false);
				
				lines.add("![MFD plot]("+resourcesDir.getName()+"/"+prefix+".png)");
				lines.add("");
				
				table = MarkdownUtils.tableBuilder();
				table.addLine("Adjustment", "Total Nucleation Rate", "Rate Gain", "Total Participation Rate", "Rate Gain");
				
				double refNuclRate = noAdjGR.getTotalIncrRate();
				double refParticRate = estParticRate(noAdjGR, nuclToParticScalars);
				for (IncrementalMagFreqDist mfd : adjTargetsList) {
					table.initNewLine();
					
					table.addColumn(mfd.getName());
					
					double totNuclRate = mfd.getTotalIncrRate();
					double totParticRate = estParticRate(mfd, nuclToParticScalars);
					
					table.addColumn((float)totNuclRate);
					if (mfd == noAdjGR)
						table.addColumn("_(N/A)_");
					else
						table.addColumn(twoDigits.format(totNuclRate/refNuclRate));
					table.addColumn((float)totParticRate);
					if (mfd == noAdjGR || mfd == fullGR)
						table.addColumn("_(N/A)_");
					else
						table.addColumn(twoDigits.format(totParticRate/refParticRate));
					table.finalizeLine();
				}
				lines.add("**Target Nucleation and Participation Rates**");
				lines.add("");
				
				lines.addAll(table.build());
				lines.add("");
				
				for (Jump jump : jumpsSorted) {
					double jumpDist = jumpDists.get(jump);
					double jumpProb = jumpProbs.get(jump);
					
					lines.add("**Jump from _"+jump.fromSection.getSectionName()+"_ to _"
							+jump.toSection.getSectionName()+"_, "+twoDigits.format(jumpDist)+" km**");
					lines.add("");
					
					table = MarkdownUtils.tableBuilder();
					table.addLine("Distance (km)", "Probability", "Mag Bins");
					
					BitSet jumpMagBitSet = jumpMagBins.get(jump);
					
					String magBinsStr = null;
					for (int i=0; i<refMFD.size(); i++) {
						if (jumpMagBitSet.get(i)) {
							if (magBinsStr == null)
								magBinsStr = "";
							else
								magBinsStr += ",";
							magBinsStr += (float)refMFD.getX(i);
						}
					}
					table.addLine(twoDigits.format(jumpDist), (float)jumpProb, magBinsStr);
					
					lines.addAll(table.build());
					lines.add("");
					
					for (boolean from : new boolean[] { true, false }) {
						FaultSection sect = from ? jump.fromSection : jump.toSection;
						
						lines.add("**"+sect.getSectionName()+"**");
						lines.add("");
						
						double[] scalars = calcNuclToParticScalars(rupSet, sect.getSectionId(), refMFD);
						
						table = MarkdownUtils.tableBuilder();
						table.addLine("Adjustment", "Sect Participation Rate", "Bin Participation Rate",
								"Seg-Implied Participation Allotment", "Fractional Allotment");
						
						for (int i=0; i<adjTargetsList.size(); i++) {
							if (adjTargetsList.get(i) == fullGR)
								continue;
							IncrementalMagFreqDist mfd = adjTargetsAllSectsList.get(i).get(subsetMappings.getNewSectID(sect.getSectionId()));
							table.initNewLine().addColumn(adjTargetsList.get(i).getName());
							
							// total participation rate for this section assuming the original supra-seis MFD is honored
							double totParticRate = 0d;
							// total participation rate for magnitudes bins that use this jump
							double jumpParticRate = 0d;
							for (int j=0; j<mfd.size(); j++) {
								double binRate = mfd.getY(j);
								double particRate = binRate*scalars[j];
								totParticRate += particRate;
								if (jumpMagBitSet.get(j))
									jumpParticRate += particRate;
//								if (binRate > 0d)
//									System.out.println((float)mfd.getX(j)+"\tnucl="+(float)binRate+"\tpartic="+(float)particRate+"\t"+magBins[j]);
							}
							// this is the segmentation-implied participation rate allotment for this jump
							double segJumpParticRate = totParticRate*jumpProb;
							// what fraction of the jump participation rate is allowed by the segmentation model?
							// this will be >1 if the segmentation constraint is more permissive than the input G-R
							double segFractOfAllotment = segJumpParticRate/jumpParticRate;
							
							table.addColumn((float)totParticRate);
							table.addColumn((float)jumpParticRate);
							table.addColumn((float)segJumpParticRate);
							table.addColumn((float)segFractOfAllotment);
							
							table.finalizeLine();
						}
						
						lines.addAll(table.build());
						lines.add("");
					}
				}
				
				segCount++;
			}
		}
		
		// add TOC
		lines.addAll(tocIndex, MarkdownUtils.buildTOC(lines, 2, 3));
		lines.add(tocIndex, "## Table Of Contents");

		// write markdown
		MarkdownUtils.writeReadmeAndHTML(lines, outputDir);
	}

	public static double[] calcNuclToParticScalars(FaultSystemRupSet rupSet, int subSectIndex, EvenlyDiscretizedFunc refMFD) {
		double[] nuclToParticScalars = new double[refMFD.size()];
		double[] avgBinAreas = new double[nuclToParticScalars.length];
		int[] avgCounts = new int[avgBinAreas.length];
		
		// loop over ruptures for which this section participates
		for (int rupIndex : rupSet.getRupturesForSection(subSectIndex)) {
			int index = refMFD.getClosestXIndex(rupSet.getMagForRup(rupIndex));
			avgCounts[index]++;
			avgBinAreas[index] += rupSet.getAreaForRup(rupIndex);
		}
		double sectArea = rupSet.getAreaForSection(subSectIndex);
		for (int m=0; m<nuclToParticScalars.length; m++) {
			if (avgCounts[m] > 0) {
				avgBinAreas[m] /= avgCounts[m];
				nuclToParticScalars[m] = avgBinAreas[m]/sectArea;
			}
		}
		return nuclToParticScalars;
	}
	
	private static void distTraverse(RuptureTreeNavigator nav, FaultSection curSect, boolean forwards,
			double curCmlDist, double curMaxDist, double[] sectCmlDists, double[] sectMaxDists) {
		Collection<FaultSection> dests;
		if (forwards) {
			dests = nav.getDescendants(curSect);
		} else {
			FaultSection prev = nav.getPredecessor(curSect);
			if (prev == null)
				return;
			dests = List.of(prev);
		}
		
		if (dests == null || dests.isEmpty())
			return;
		
		for (FaultSection dest : dests) {
			double myCmlDist = curCmlDist;
			double myMaxDist = curMaxDist;
			if (dest.getParentSectionId() != curSect.getParentSectionId()) {
				Jump jump = nav.getJump(dest, curSect);
				myCmlDist += jump.distance;
				myMaxDist = Math.max(myMaxDist, jump.distance);
			}
			int myID = dest.getSectionId();
			if (Double.isNaN(sectCmlDists[myID])) {
				sectCmlDists[myID] = myCmlDist;
				sectMaxDists[myID] = myMaxDist;
			} else {
				sectCmlDists[myID] = Math.min(myCmlDist, sectCmlDists[myID]);
				sectMaxDists[myID] = Math.min(myMaxDist, sectMaxDists[myID]);
			}
			
			distTraverse(nav, dest, forwards, myCmlDist, myMaxDist, sectCmlDists, sectMaxDists);
		}
	}
	
	private static double estParticRate(IncrementalMagFreqDist mfd, double[] nuclToParticScalars) {
		double ret = 0;
		for (int i=0; i<mfd.size(); i++)
			ret += mfd.getY(i)*nuclToParticScalars[i];
		return ret;
	}
	
	private static class QuickAvgTrack {
		
		private double sum;
		private double first;
		private boolean allSame;
		private int count;
		
		public QuickAvgTrack() {
			sum = 0d;
			first = Double.NaN;
			allSame = true;
			count = 0;
		}
		
		public void add(double value) {
			Preconditions.checkState(value >= 0d);
			sum += value;
			if (count == 0)
				first = value;
			count++;
			allSame = allSame && value == first;
		}
		
		public double getAverage() {
			Preconditions.checkState(count > 0);
			if (allSame)
				return first;
			return sum/(double)count;
		}
	}

}
