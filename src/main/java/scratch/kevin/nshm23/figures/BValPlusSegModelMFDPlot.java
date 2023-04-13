package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.function.Supplier;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbDiscrEmpiricalDistFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.NSHM23_InvConfigFactory;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_RegionalSeismicity;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SegmentationMFD_Adjustment;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.Builder;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.estimators.SectNucleationMFD_Estimator;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader.AnalysisRegions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import com.google.common.base.Preconditions;

public class BValPlusSegModelMFDPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/home/kevin/Documents/papers/2023_NSHM23_Inversion/figures/bval_seg_comb_mfds");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
//		String regPrefix = null;
//		Region region = null;
		String regPrefix = "wus";
		Region region = NSHM23_RegionLoader.loadFullConterminousWUS();
		
//		int[] miniClusterIDs = {2063, 2062, 2045};
		int[] miniClusterIDs = {2568, 2522, 2582};
		
		FaultSystemRupSet fullRupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		FaultSystemRupSet rupSet5km = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/UCERF4/rup_sets/"
				+ "NSHM23_v2_plausibleMulti5km_direct_cmlRake360_jumpP0.001_slipP0.05incrCapDist_cff0.75IntsPos_"
				+ "comb2Paths_cffFav10P0.01_cffFav10RatioN2P0.5_sectFractGrow0.1.zip"));
		
		SupraSeisBValues[] bVals = {
				SupraSeisBValues.B_1p0,
				SupraSeisBValues.B_0p75,
				SupraSeisBValues.B_0p5,
				SupraSeisBValues.B_0p25,
				SupraSeisBValues.B_0p0,
		};
		
		CPT bValCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().reverse().rescale(0d, bVals.length-1d);
		
		NSHM23_SegmentationModels[] segModels = {
				NSHM23_SegmentationModels.CLASSIC,
				NSHM23_SegmentationModels.HIGH,
				NSHM23_SegmentationModels.MID,
				NSHM23_SegmentationModels.LOW,
				NSHM23_SegmentationModels.NONE,
		};
		
		CPT segModelCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().reverse().rescale(0d, segModels.length-1d);
		
		int indvTrans = 100;
		
//		float lineThickness = 2f;
//		float symbolWidth = 3f;
//		PlotCurveCharacterstics[] segModelChars = {
//				new PlotCurveCharacterstics(PlotLineType.SOLID, lineThickness, PlotSymbol.FILLED_CIRCLE, symbolWidth, null),
//				new PlotCurveCharacterstics(PlotLineType.DASHED, lineThickness, PlotSymbol.FILLED_DIAMOND, symbolWidth, null),
//				new PlotCurveCharacterstics(PlotLineType.DOTTED, lineThickness, PlotSymbol.FILLED_INV_TRIANGLE, symbolWidth, null),
//				new PlotCurveCharacterstics(PlotLineType.SOLID, lineThickness, PlotSymbol.FILLED_SQUARE, symbolWidth, null),
//				new PlotCurveCharacterstics(PlotLineType.DASHED, lineThickness, PlotSymbol.BOLD_CROSS, symbolWidth, null),
//		};
		
		double[] sectFractsInReg = region == null ? null : fullRupSet.getFractSectsInsideRegion(region, false);
		UncertainBoundedIncrMagFreqDist observedIncr = null;
		UncertainIncrMagFreqDist observedIncrBounds = null;
		EvenlyDiscretizedFunc observedCml = null;
		UncertainArbDiscFunc observedCmlBounds = null;
		if (region != null) {
			EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(8.95);
			observedIncr = NSHM23_RegionalSeismicity.getRemapped(region,
					NSHM23_DeclusteringAlgorithms.AVERAGE, NSHM23_SeisSmoothingAlgorithms.AVERAGE, refMFD, refMFD.getMaxX());
			observedIncr.setName("Observed");
			observedIncr.setBoundName("95% Bounds");
			observedIncrBounds = observedIncr.deepClone();
			observedIncrBounds.setName(observedIncr.getBoundName());
			
			// add observed bounds
			observedCml = Regional_MFD_Plots.getCmlAsFakeIncr(observedIncr);
			observedCml = observedCml.deepClone();
			observedCml.setName(observedIncr.getName());
			
			EvenlyDiscretizedFunc upperCumulative = Regional_MFD_Plots.getCmlAsFakeIncr(observedIncr.getUpper());
			EvenlyDiscretizedFunc lowerCumulative = Regional_MFD_Plots.getCmlAsFakeIncr(observedIncr.getLower());
			Preconditions.checkState(observedCml.size() == upperCumulative.size());
			for (int i=0; i<observedCml.size(); i++) {
				upperCumulative.set(i, Math.max(observedCml.getY(i), upperCumulative.getY(i)));
				lowerCumulative.set(i, Math.max(0, Math.min(observedCml.getY(i), lowerCumulative.getY(i))));
			}
			
			observedCmlBounds = new UncertainArbDiscFunc(observedCml, lowerCumulative, upperCumulative);
			observedCmlBounds.setName(observedIncr.getBoundName());
		}

		LogicTreeBranch<LogicTreeNode> branch = NSHM23_LogicTreeBranch.DEFAULT_ON_FAULT.copy();
		branch.setValue(NSHM23_DeformationModels.AVERAGE);
		branch.setValue(NSHM23_ScalingRelationships.AVERAGE);
		fullRupSet = new NSHM23_InvConfigFactory().updateRuptureSetForBranch(fullRupSet, branch);
		rupSet5km = new NSHM23_InvConfigFactory().updateRuptureSetForBranch(rupSet5km, branch);
		
		SegmentationMFD_Adjustment segAdjMethod = branch.getValue(SegmentationMFD_Adjustment.class);
		
		IncrementalMagFreqDist[][] allCurves = new IncrementalMagFreqDist[bVals.length][segModels.length];
		
		for (boolean mini : new boolean[] {true,false}) {
			FaultSystemRupSet rupSet;
			if (mini) {
				HashSet<Integer> sectIDs = new HashSet<>();
				for (FaultSection sect : fullRupSet.getFaultSectionDataList()) {
					for (int miniID : miniClusterIDs) {
						if (sect.getParentSectionId() == miniID) {
							System.out.println("Match for "+miniID+": "+sect.getSectionId()+". "+sect.getSectionName());
							sectIDs.add(sect.getSectionId());
						}
					}
				}
				rupSet = fullRupSet.getForSectionSubSet(sectIDs);
			} else {
				rupSet = fullRupSet;
			}
			
			List<DiscretizedFunc> bValIncrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> bValIncrChars = new ArrayList<>();
			
			List<DiscretizedFunc> bValCmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> bValCmlChars = new ArrayList<>();
			
			List<DiscretizedFunc> segIncrFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> segIncrChars = new ArrayList<>();
			
			List<DiscretizedFunc> segCmlFuncs = new ArrayList<>();
			List<PlotCurveCharacterstics> segCmlChars = new ArrayList<>();
			
			if (!mini && observedIncr != null) {
				Color obsColor = new Color(125, 80, 145); // "indigo"
				PlotCurveCharacterstics obsChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, obsColor);
				PlotCurveCharacterstics obsBoundsChar = new PlotCurveCharacterstics(
						PlotLineType.SHADED_UNCERTAIN, 1f, new Color(obsColor.getRed(), obsColor.getGreen(), obsColor.getBlue(), 40));
				bValIncrFuncs.add(observedIncr);
				bValIncrChars.add(obsChar);
				bValIncrFuncs.add(observedIncrBounds);
				bValIncrChars.add(obsBoundsChar);
				
				segIncrFuncs.add(observedIncr);
				segIncrChars.add(obsChar);
				segIncrFuncs.add(observedIncrBounds);
				segIncrChars.add(obsBoundsChar);
				
				bValCmlFuncs.add(observedCml);
				bValCmlChars.add(obsChar);
				bValCmlFuncs.add(observedCmlBounds);
				bValCmlChars.add(obsBoundsChar);
				
				segCmlFuncs.add(observedCml);
				segCmlChars.add(obsChar);
				segCmlFuncs.add(observedCmlBounds);
				segCmlChars.add(obsBoundsChar);
			}
			
			for (int b=0; b<bVals.length; b++) {
				SupraSeisBValues bVal = bVals[b];
				branch.setValue(bVal);
				
				List<CompletableFuture<IncrementalMagFreqDist>> mfdFutures = new ArrayList<>();
				
				for (int s=0; s<segModels.length; s++) {
					NSHM23_SegmentationModels segModel = segModels[s];
					branch.setValue(segModel);
					
					JumpProbabilityCalc model = segModel.getModel(rupSet, branch);
					
					System.out.println("Building for "+bVal+", "+segModel);
					Builder builder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bVal.bValue);
					builder.applyDefModelUncertainties(false);
					if (model != null) {
						SectNucleationMFD_Estimator adjustment = segAdjMethod.getAdjustment(model);
						builder.adjustTargetsForData(adjustment);
					}
					
					mfdFutures.add(CompletableFuture.supplyAsync(new Supplier<IncrementalMagFreqDist>() {

						@Override
						public IncrementalMagFreqDist get() {
							return getTotalSupraMFD(rupSet, builder.build(), sectFractsInReg);
						}
					}));
				}
				for (int s=0; s<segModels.length; s++) {
					IncrementalMagFreqDist mfd = mfdFutures.get(s).join();
					mfd.setName(null);
					
//					PlotCurveCharacterstics pChar = segModelChars[s];
//					pChar = new PlotCurveCharacterstics(pChar.getLineType(), pChar.getLineWidth(),
//							pChar.getSymbol(), pChar.getSymbolWidth(), bValCPT.getColor((float)b));
					
					for (boolean seg : new boolean[] {true,false}) {
						List<DiscretizedFunc> incrFuncs = seg ? segIncrFuncs : bValIncrFuncs;
						List<PlotCurveCharacterstics> incrChars = seg ? segIncrChars : bValIncrChars;
						List<DiscretizedFunc> cmlFuncs = seg ? segCmlFuncs : bValCmlFuncs;
						List<PlotCurveCharacterstics> cmlChars = seg ? segCmlChars : bValCmlChars;
						Color color;
						if (seg)
							color = segModelCPT.getColor((float)s);
						else
							color = bValCPT.getColor((float)b);
						color = new Color(color.getRed(), color.getGreen(), color.getBlue(), indvTrans);
						PlotCurveCharacterstics pChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 1.5f, color);
						
						allCurves[b][s] = mfd;
//						BValueSumCartoon.addMFDFuncs(incrFuncs, incrChars, totMFD, pChar, bVal.bValue);
//						incrFuncs.get(incrFuncs.size()-1).setName(name);
						incrFuncs.add(mfd);
						incrChars.add(pChar);
						
						DiscretizedFunc cmlFunc = mfd.getCumRateDistWithOffset();
						cmlFunc.setName(null);
//						cmlFunc = truncatedCumulative(cmlFunc);
						cmlFuncs.add(cmlFunc);
						cmlChars.add(pChar);
					}
				}
			}
			
			IncrementalMagFreqDist mfd5km = null;
			if (!mini) {
				System.out.println("Calculating for 5km rupture set");
				
				List<CompletableFuture<IncrementalMagFreqDist>> mfdFutures = new ArrayList<>();
				
				for (SupraSeisBValues bVal : bVals) {
					Builder builder = new SupraSeisBValInversionTargetMFDs.Builder(rupSet5km, bVal.bValue);
					builder.applyDefModelUncertainties(false);
					
					FaultSystemRupSet rs = rupSet5km;
					mfdFutures.add(CompletableFuture.supplyAsync(new Supplier<IncrementalMagFreqDist>() {
						@Override
						public IncrementalMagFreqDist get() {
//							return builder.build().getTotalOnFaultSupraSeisMFD();
							return getTotalSupraMFD(rs, builder.build(), sectFractsInReg);
						}
					}));
				}
				
				for (CompletableFuture<IncrementalMagFreqDist> future : mfdFutures) {
					IncrementalMagFreqDist mfd = future.join();
					
					if (mfd5km == null)
						mfd5km = new IncrementalMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
					else
						Preconditions.checkState(mfd.getMinX() == mfd5km.getMinX() && mfd.size() == mfd5km.size());
					for (int j=0; j<mfd.size(); j++)
						mfd5km.add(j, mfd.getY(j));
				}
				// TODO weights if we make them uneven
				mfd5km.scale(1d/mfdFutures.size());
			}
			
			String title = " ";
			String xAxisName = "Magnitude";
			Range xRange = mini ? new Range(6.5d, 7.5d) : new Range(6d, 8.5d);
			Range yRange = mini ? new Range(1e-8, 1e-3) : new Range(1e-4, 2e0);
			
			for (boolean seg : new boolean[] {true,false}) {
				List<DiscretizedFunc> incrFuncs = seg ? segIncrFuncs : bValIncrFuncs;
				List<PlotCurveCharacterstics> incrChars = seg ? segIncrChars : bValIncrChars;
				List<DiscretizedFunc> cmlFuncs = seg ? segCmlFuncs : bValCmlFuncs;
				List<PlotCurveCharacterstics> cmlChars = seg ? segCmlChars : bValCmlChars;
				
				List<List<IncrementalMagFreqDist>> binnedIncrs = new ArrayList<>();
				List<String> names = new ArrayList<>();
				List<Color> colors = new ArrayList<>();
				
				if (seg) {
					// bin by seg model
					for (NSHM23_SegmentationModels model : segModels) {
						colors.add(segModelCPT.getColor((float)colors.size()));
						binnedIncrs.add(new ArrayList<>());
						names.add(model.getName());
					}
					for (int b=0; b<bVals.length; b++)
						for (int s=0; s<segModels.length; s++)
							binnedIncrs.get(s).add(allCurves[b][s]);
				} else {
					// bin by b value
					for (SupraSeisBValues bVal : bVals) {
						colors.add(bValCPT.getColor((float)colors.size()));
						binnedIncrs.add(new ArrayList<>());
						names.add(bVal.getName());
					}
					for (int b=0; b<bVals.length; b++)
						for (int s=0; s<segModels.length; s++)
							binnedIncrs.get(b).add(allCurves[b][s]);
				}
				
				IncrementalMagFreqDist fullAVG = null;
				
				for (int i=0; i<binnedIncrs.size(); i++) {
					IncrementalMagFreqDist avg = null;
					for (IncrementalMagFreqDist mfd : binnedIncrs.get(i)) {
						if (avg == null)
							avg = new IncrementalMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
						else
							Preconditions.checkState(mfd.getMinX() == avg.getMinX() && mfd.size() == avg.size());
						if (fullAVG == null)
							fullAVG = new IncrementalMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
						else
							Preconditions.checkState(mfd.getMinX() == fullAVG.getMinX() && mfd.size() == fullAVG.size());
						for (int j=0; j<mfd.size(); j++) {
							avg.add(j, mfd.getY(j));
							fullAVG.add(j, mfd.getY(j));
						}
					}
					// TODO weights if we make them uneven
					avg.scale(1d/binnedIncrs.size());
					
					Color color = colors.get(i);
					color = new Color(color.getRed()*0.95f/255f, color.getGreen()*0.95f/255f, color.getBlue()*0.95f/255f);
					PlotCurveCharacterstics pChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, color);
					
					avg.setName(names.get(i));
					incrFuncs.add(avg);
					incrChars.add(pChar);
					
					EvenlyDiscretizedFunc cml = avg.getCumRateDistWithOffset();
					cml.setName(names.get(i));
					cmlFuncs.add(cml);
					cmlChars.add(pChar);
				}
				
				fullAVG.scale(1d/(segModels.length*bVals.length));
				fullAVG.setName("Average");
				PlotCurveCharacterstics pChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK);
				
				incrFuncs.add(fullAVG);
				incrChars.add(pChar);
				
				EvenlyDiscretizedFunc cml = fullAVG.getCumRateDistWithOffset();
				cml.setName(fullAVG.getName());
				cmlFuncs.add(cml);
				cmlChars.add(pChar);
				
				if (mfd5km != null && seg) {
					// add 5km none branch for comparison
					mfd5km.setName("None, 5km Max");
					pChar = new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY);
					
					incrFuncs.add(mfd5km);
					incrChars.add(pChar);
					
					cml = mfd5km.getCumRateDistWithOffset();
					cml.setName(mfd5km.getName());
					cmlFuncs.add(cml);
					cmlChars.add(pChar);
				}
				
				PlotSpec incrSpec = new PlotSpec(incrFuncs, incrChars, title, xAxisName, "Incremental Rate (/yr)");
				incrSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				PlotSpec cmlSpec = new PlotSpec(cmlFuncs, cmlChars, title, xAxisName, "Cumulative Rate (/yr)");
				cmlSpec.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
				
				HeadlessGraphPanel gp = PlotUtils.initHeadless();
				
				gp.setAxisLabelFontSize(26);
				gp.setTickLabelFontSize(22);
				
				gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
				
				String prefix = (mini ? "mini" : "full")+"_"+(seg ? "seg" : "bval")+"_mfds";
				if (regPrefix != null)
					prefix = regPrefix+"_"+prefix;
				
				PlotUtils.writePlots(outputDir, prefix+"_incr", gp, 800, 750, true, true, false);
				
				gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
				
				PlotUtils.writePlots(outputDir, prefix+"_cml", gp, 800, 750, true, true, false);
			}
		}
	}
	
	private static DiscretizedFunc truncatedCumulative(DiscretizedFunc cmlFunc) {
		DiscretizedFunc ret = new ArbDiscrEmpiricalDistFunc();
		for (Point2D pt : cmlFunc) {
			if (pt.getY() == 0d)
				break;
			ret.set(pt);
		}
		return ret;
	}
	
	private static IncrementalMagFreqDist getTotalSupraMFD(FaultSystemRupSet rupSet,
			SupraSeisBValInversionTargetMFDs mfds, double[] sectFractsInReg) {
		if (sectFractsInReg != null) {
			IncrementalMagFreqDist sum = null;
			for (int s=0; s<rupSet.getNumSections(); s++) {
				double fract = sectFractsInReg[s];
				if (fract > 0d) {
					IncrementalMagFreqDist mfd = mfds.getOnFaultSupraSeisNucleationMFDs().get(s);
					if (sum == null)
						sum = new IncrementalMagFreqDist(mfd.getMinX(), mfd.size(), mfd.getDelta());
					else
						Preconditions.checkState(mfd.getMinX() == sum.getMinX() && mfd.size() == sum.size());
					for (int i=0; i<mfd.size(); i++)
						sum.add(i, mfd.getY(i)*fract);
				}
			}
			return sum;
		}
		return mfds.getTotalOnFaultSupraSeisMFD();
	}

}
