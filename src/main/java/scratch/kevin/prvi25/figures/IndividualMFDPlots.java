package scratch.kevin.prvi25.figures;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.util.Precision;
import org.jfree.data.Range;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.BoundedUncertainty;
import org.opensha.commons.data.uncertainty.UncertainArbDiscFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
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
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.RegionsOfInterest;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.PRVI25_GridSourceBuilder;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.gridded.SeismicityRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_CrustalSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_DeclusteringAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeisSmoothingAlgorithms;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SeismicityRateEpoch;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionCaribbeanSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionFaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionMuertosSeismicityRate;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.logicTree.PRVI25_SubductionSlabMMax;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader.PRVI25_SeismicityRegions;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Ints;

import gov.usgs.earthquake.nshmp.model.NshmErf;
import net.mahdilamb.colormap.Colors;
import scratch.kevin.latex.LaTeXUtils;
import scratch.kevin.prvi25.SubductionCombinedModelCreator;

import static scratch.kevin.prvi25.figures.PRVI_Paths.*;

public class IndividualMFDPlots {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution crustalSol = FaultSystemSolution.load(CRUSTAL_SOL_GRIDDED);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(4.05, 9.55);
		
		FaultSystemSolution subductionLargeSol = FaultSystemSolution.load(SUBDUCTION_SOL_LARGE);
		FaultSystemSolution subductionSmallSol = FaultSystemSolution.load(SUBDUCTION_SOL_SMALL);
		FaultSystemSolution subductionCombined = SubductionCombinedModelCreator.combine(
				Map.of(PRVI25_SubductionFaultModels.PRVI_SUB_FM_LARGE, subductionLargeSol,
						PRVI25_SubductionFaultModels.PRVI_SUB_FM_SMALL, subductionSmallSol));
		
		File crustalOutputDir = new File(FIGURES_DIR, "crustal_sol");
		Preconditions.checkState(crustalOutputDir.exists() || crustalOutputDir.mkdir());
		File subductionOutputDir = new File(FIGURES_DIR, "sub_sol");
		Preconditions.checkState(subductionOutputDir.exists() || subductionOutputDir.mkdir());
		
		FileWriter crustalTexFW = new FileWriter(new File(crustalOutputDir, "crustal_mfds.tex"));
		FileWriter subductionTexFW = new FileWriter(new File(subductionOutputDir, "subduction_mfds.tex"));
		
		Map<PRVI25_SeismicityRegions, Region> subsetRegions = Map.of(
				PRVI25_SeismicityRegions.CRUSTAL, PRVI25_RegionLoader.loadPRVI_MapExtents());
		
		for (PRVI25_SeismicityRegions seisReg : PRVI25_SeismicityRegions.values()) {
			Region[] regions;
			if (subsetRegions.containsKey(seisReg))
				regions = new Region[] {seisReg.load(), subsetRegions.get(seisReg)};
			else
				regions = new Region[] {seisReg.load()};
			
			for (int r=0; r<regions.length; r++) {
				Region reg = regions[r];
				
				FaultSystemSolution sol;
				File outputDir;
				String prefix;
				TectonicRegionType trt;
				Range xRange;
				String texPrefix;
				FileWriter texFW;
				Function<PRVI25_SeismicityRateEpoch, SeismicityRateModel> seisModelFunc;
				if (seisReg == PRVI25_SeismicityRegions.CRUSTAL) {
					sol = crustalSol;
					outputDir = crustalOutputDir;
					prefix = "crustal_mfds";
					trt = TectonicRegionType.ACTIVE_SHALLOW;
					xRange = new Range(5d, 8.5);
					texPrefix = "CrustalMFD";
					texFW = crustalTexFW;
					seisModelFunc = epoch -> {
						try {
							return PRVI25_CrustalSeismicityRate.loadRateModel(epoch);
						} catch (IOException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					};
				} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTERFACE) {
					sol = subductionCombined;
					outputDir = subductionOutputDir;
					prefix = "subduction_mfds_car_interface";
					trt = TectonicRegionType.SUBDUCTION_INTERFACE;
					xRange = new Range(5d, 9.5);
					texPrefix = "SubCarIntMFD";
					texFW = subductionTexFW;
					seisModelFunc = epoch -> {
						try {
							return PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epoch, false);
						} catch (IOException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					};
				} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTERFACE) {
					sol = subductionCombined;
					outputDir = subductionOutputDir;
					prefix = "subduction_mfds_mue_interface";
					trt = TectonicRegionType.SUBDUCTION_INTERFACE;
					xRange = new Range(5d, 9.5);
					texPrefix = "SubMueIntMFD";
					texFW = subductionTexFW;
					seisModelFunc = epoch -> {
						try {
							return PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epoch, false);
						} catch (IOException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					};
				} else if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB) {
					sol = subductionCombined;
					outputDir = subductionOutputDir;
					prefix = "subduction_mfds_car_slab";
					trt = TectonicRegionType.SUBDUCTION_SLAB;
					xRange = new Range(5d, 9.5);
					texPrefix = "SubCarSlabMFD";
					texFW = subductionTexFW;
					seisModelFunc = epoch -> {
						try {
							return PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epoch, true);
						} catch (IOException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					};
				} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTRASLAB) {
					sol = subductionCombined;
					outputDir = subductionOutputDir;
					prefix = "subduction_mfds_mue_slab";
					trt = TectonicRegionType.SUBDUCTION_SLAB;
					xRange = new Range(5d, 9.5);
					texPrefix = "SubMueSlabMFD";
					texFW = subductionTexFW;
					seisModelFunc = epoch -> {
						try {
							return PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epoch, true);
						} catch (IOException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
					};
				} else {
					throw new IllegalStateException();
				}
				
				if (r > 0) {
					texPrefix += "Subset";
					prefix += "_subset";
				}
				
				UncertainBoundedIncrMagFreqDist[] onFaultDists = null;
				if (trt == TectonicRegionType.ACTIVE_SHALLOW)
					onFaultDists = getOnFaultMFDFractiles(sol, r == 0 ? null : reg);
				else if (trt == TectonicRegionType.SUBDUCTION_INTERFACE)
					onFaultDists = getRegionalMFDFractiles(sol, MFDType.SUPRA_ONLY, reg, fractiles);
				IncrementalMagFreqDist onFaultMean = null;
				if (trt == TectonicRegionType.ACTIVE_SHALLOW)
					onFaultMean = CombinedMFDsPlot.calcFaultMFD(r == 0 ? null : reg, sol, refMFD);
				else if (trt == TectonicRegionType.SUBDUCTION_INTERFACE)
					onFaultMean = CombinedMFDsPlot.calcFaultMFD(reg, sol, refMFD);
				
				List<UncertainBoundedIncrMagFreqDist> obsList = new ArrayList<>();
				List<UncertainBoundedIncrMagFreqDist> origObsList = r == 0 ? null : new ArrayList<>();
				List<Double> obsWeights = new ArrayList<>();
				for (PRVI25_SeismicityRateEpoch epoch : PRVI25_SeismicityRateEpoch.values()) {
					double weight = epoch.getNodeWeight(null);
					if (weight == 0d)
						continue;
					SeismicityRateModel seisModel = seisModelFunc.apply(epoch);
					UncertainBoundedIncrMagFreqDist obs = seisModel.getBounded(refMFD, xRange.getUpperBound()+0.1);
					if (r == 0) {
						obsList.add(obs);
					} else {
						UncertainBoundedIncrMagFreqDist subsetObs = seisModel.getRemapped(reg, seisReg, PRVI25_DeclusteringAlgorithms.AVERAGE,
								PRVI25_SeisSmoothingAlgorithms.AVERAGE, refMFD, xRange.getUpperBound()+0.1);
						subsetObs.setName("Observed (subset), N5="+new DecimalFormat("0.0#").format(obs.getCumRate(obs.getClosestXIndex(5.01))));
						obsList.add(subsetObs);
						origObsList.add(obs);
					}
					obsWeights.add(weight);
				}
				UncertainBoundedIncrMagFreqDist obs = PRVI25_SeismicityRateEpoch.averageUncert(obsList, obsWeights);
				IncrementalMagFreqDist gridded;
				UncertainBoundedIncrMagFreqDist[] griddedDists;
				if (trt == TectonicRegionType.ACTIVE_SHALLOW || trt == TectonicRegionType.SUBDUCTION_INTERFACE) {
					gridded = calcGriddedMFD(reg, trt, sol, refMFD);
					griddedDists = getRegionalMFDFractiles(sol, MFDType.GRID_ONLY, reg, extrema);
				} else {
					// regions overlap, need to redo it
					Preconditions.checkState(r == 0, "Need to do bounds right if we decide to use a slab subset region");
					UncertainBoundedIncrMagFreqDist avgSlab = getMmaxAveragedSlab(refMFD, seisReg, true);
					gridded = avgSlab;
					griddedDists = new UncertainBoundedIncrMagFreqDist[] { avgSlab };
				}

				if (r > 0)
					texFW.write("% "+seisReg.name()+" subset "+r+"\n");
				else
					texFW.write("% "+seisReg.name()+"\n");
				
				if (trt == TectonicRegionType.SUBDUCTION_SLAB) {
					if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB) {
						// combined slab plot
						IncrementalMagFreqDist slab03 = loadSlab03(refMFD);
						// false here means pref not average
						UncertainBoundedIncrMagFreqDist carDist = getMmaxAveragedSlab(refMFD, seisReg, false);
						UncertainBoundedIncrMagFreqDist mueDist = getMmaxAveragedSlab(refMFD, PRVI25_SeismicityRegions.MUE_INTRASLAB, false);
//						UncertainBoundedIncrMagFreqDist carDist = siesModel.getBounded(refMFD, PRVI25_GridSourceBuilder.SLAB_MMAX);
//						UncertainBoundedIncrMagFreqDist mueDist = PRVI25_SubductionMuertosSeismicityRate.loadRateModel(true)
//								.getBounded(refMFD, PRVI25_GridSourceBuilder.SLAB_MMAX, PRVI25_GridSourceBuilder.SLAB_M_CORNER);
						plotMultiSlab(subductionOutputDir, "subduction_mfds_slab_combined", new Range(5d, 8d),
								carDist,
								mueDist,
								slab03);
						double prevM5 = slab03.getCumRate(slab03.getClosestXIndex(5.01));
						texFW.write(LaTeXUtils.defineValueCommand("SubSlabPrevModelMFDMFiveRate",
								LaTeXUtils.numberExpFormatSigFigs(prevM5, 2), false)+"\n");
						texFW.write(LaTeXUtils.defineValueCommand("SubSlabPrevModelMFDMFiveRI",
								LaTeXUtils.numberExpFormatFixedDecimal(1d/prevM5, 1), false)+"\n");
						IncrementalMagFreqDist mueGridded = getMmaxAveragedSlab(refMFD, PRVI25_SeismicityRegions.MUE_INTRASLAB, true);
						double sumSlabRate = gridded.getCumRate(gridded.getClosestXIndex(5.01)) + mueGridded.getCumRate(mueGridded.getClosestXIndex(5.01));
						texFW.write(LaTeXUtils.defineValueCommand("SubSlabCombModelMFDMFiveRate",
								LaTeXUtils.numberExpFormatSigFigs(sumSlabRate, 2), false)+"\n");
						texFW.write(LaTeXUtils.defineValueCommand("SubSlabCombModelMFDMFiveRI",
								LaTeXUtils.numberExpFormatFixedDecimal(1d/sumSlabRate, 1), false)+"\n");
						texFW.write(LaTeXUtils.defineValueCommand("SubSlabCombModelMFDMFiveChange",
								LaTeXUtils.numberAsPercent(100d*(sumSlabRate-prevM5)/prevM5, r), false)+"\n");
					}
				}
				double obsM5 = obs.getCumRate(obs.getClosestXIndex(5.01));
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"ObsMFiveRate",
						LaTeXUtils.numberExpFormatSigFigs(obsM5, 2), false)+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"ObsMFiveRI",
						LaTeXUtils.numberExpFormatFixedDecimal(1d/obsM5, 1), false)+"\n");
				if (r > 0) {
					UncertainBoundedIncrMagFreqDist origObs = PRVI25_SeismicityRateEpoch.averageUncert(origObsList, obsWeights);
					double origM5 = origObs.getCumRate(origObs.getClosestXIndex(5.01));
					texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"ObsMFivePercent",
							LaTeXUtils.numberAsPercent(100d*obsM5/origM5, 0), false)+"\n");
				}
				double gridM5 = gridded.getCumRate(gridded.getClosestXIndex(5.01));
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"GridMFiveRate",
						LaTeXUtils.numberExpFormatSigFigs(gridM5, 2), false)+"\n");
				texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"GridMFiveRI",
						LaTeXUtils.numberExpFormatFixedDecimal(1d/gridM5, 1), false)+"\n");
				if (trt != TectonicRegionType.SUBDUCTION_SLAB) {
					double onFaultTotalRate = onFaultMean.calcSumOfY_Vals();
					texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"SupraRate",
							LaTeXUtils.numberExpFormatSigFigs(onFaultTotalRate, 3), false)+"\n");
					double onFaultRI = 1d/onFaultTotalRate;
					if (onFaultRI > 5d)
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"SupraRI",
								LaTeXUtils.groupedIntNumber(onFaultRI), false)+"\n");
					else
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"SupraRI",
								LaTeXUtils.numberExpFormatFixedDecimal(onFaultRI, 1), false)+"\n");
					if (seisReg == PRVI25_SeismicityRegions.CRUSTAL && r == 0) {
						String[] indvFaultNames = {
								"Anegada",
								"Bunce",
								"Septentrional",
								"South Lajas",
								"Mona Passage"
						};
						double multiFaultRate = 0d;
						double[] indvFaultRates = new double[indvFaultNames.length];
						for (int rupIndex=0; rupIndex<sol.getRupSet().getNumRuptures(); rupIndex++) {
							boolean multiFault = false;
							List<FaultSection> rupSects = sol.getRupSet().getFaultSectionDataForRupture(rupIndex);
							for (int i=1; !multiFault&&i<rupSects.size(); i++)
								multiFault = rupSects.get(i-1).getParentSectionId() != rupSects.get(i).getParentSectionId();
							if (multiFault)
								multiFaultRate += sol.getRateForRup(rupIndex);
							for (int i=0; i<indvFaultNames.length; i++) {
								for (FaultSection sect : rupSects) {
									if (sect.getParentSectionName().contains(indvFaultNames[i])) {
										indvFaultRates[i] += sol.getRateForRup(rupIndex);
										break;
									}
								}
							}
						}
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"SupraMultiFaultRate",
								LaTeXUtils.numberExpFormatSigFigs(multiFaultRate, 3), false)+"\n");
						double multiFaultRI = 1d/multiFaultRate;
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"SupraMultiFaultRI",
								LaTeXUtils.groupedIntNumber(multiFaultRI), false)+"\n");
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"SupraMultiFaultPercent",
								LaTeXUtils.numberAsPercent(100d*multiFaultRate/onFaultTotalRate, 1), false)+"\n");
						
						for (int i=0; i<indvFaultNames.length; i++) {
							double indvRate = indvFaultRates[i];
							String indvPrefix = texPrefix+"Supra"+LaTeXUtils.sanitizeForCommandName(indvFaultNames[i]);
							texFW.write(LaTeXUtils.defineValueCommand(indvPrefix+"Rate",
									LaTeXUtils.numberExpFormatSigFigs(indvRate, 3), false)+"\n");
							double indvRI = 1d/indvRate;
							texFW.write(LaTeXUtils.defineValueCommand(indvPrefix+"RI",
									LaTeXUtils.groupedIntNumber(indvRI), false)+"\n");
							texFW.write(LaTeXUtils.defineValueCommand(indvPrefix+"Percent",
									LaTeXUtils.numberAsPercent(100d*indvRate/onFaultTotalRate, 0), false)+"\n");
							
							if (indvFaultNames[i].equals("Anegada")) {
								// add segmentation-model specific RIs
								File baDir = new File(CRUSTAL_DIR, "node_branch_averaged");
								Preconditions.checkState(baDir.exists(), "Need node branch-averaged solutions: %s", baDir.getAbsolutePath());
								for (NSHM23_SegmentationModels segModel : NSHM23_SegmentationModels.values()) {
									File segFile = new File(new File(CRUSTAL_DIR, "node_branch_averaged"), "SegModel_"+segModel.getFilePrefix()+".zip");
									if (segFile.exists()) {
										FaultSystemSolution segSol = FaultSystemSolution.load(segFile);
										FaultSystemRupSet segRS = segSol.getRupSet();
										int parentID1 = FaultSectionUtils.findParentSectionID(segRS.getFaultSectionDataList(), "Anegada", "SW PROXY");
										int parentID2 = FaultSectionUtils.findParentSectionID(segRS.getFaultSectionDataList(), "Anegada", "NE PROXY");
										SummedMagFreqDist summedMFD = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
										summedMFD.addIncrementalMagFreqDist(segSol.calcNucleationMFD_forParentSect(
												parentID1, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size()));
										summedMFD.addIncrementalMagFreqDist(segSol.calcNucleationMFD_forParentSect(
												parentID2, refMFD.getMinX(), refMFD.getMaxX(), refMFD.size()));
										double segRate = summedMFD.calcSumOfY_Vals();
										double segRate7p5 = summedMFD.getCumRate(summedMFD.getClosestXIndex(7.501));
										String segPefix = indvPrefix+segModel.getFilePrefix();
										texFW.write(LaTeXUtils.defineValueCommand(segPefix+"Rate",
												LaTeXUtils.numberExpFormatSigFigs(segRate, 3), false)+"\n");
										texFW.write(LaTeXUtils.defineValueCommand(segPefix+"RI",
												LaTeXUtils.groupedIntNumber(1d/segRate), false)+"\n");
										if (segRate7p5 > 0) {
											segPefix += "SevenFive";
											texFW.write(LaTeXUtils.defineValueCommand(segPefix+"Rate",
													LaTeXUtils.numberExpFormatSigFigs(segRate7p5, 3), false)+"\n");
											texFW.write(LaTeXUtils.defineValueCommand(segPefix+"RI",
													LaTeXUtils.groupedIntNumber(1d/segRate7p5), false)+"\n");
										}
									}
								}
							}
						}
						// add south lajas M>7
						double lajasM7Rate = sol.calcParticipationMFD_forParentSect(
								FaultSectionUtils.findParentSectionID(sol.getRupSet().getFaultSectionDataList(), "South Lajas"),
								refMFD.getMinX(), refMFD.getMaxX(), refMFD.size()).getCumRate(refMFD.getClosestXIndex(7.01));
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MSevenSouthLajasRate",
								LaTeXUtils.numberExpFormatSigFigs(lajasM7Rate, 3), false)+"\n");
						double lajasM7RI = 1d/lajasM7Rate;
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"MSevenSouthLajasRI",
								LaTeXUtils.groupedIntNumber(lajasM7RI), false)+"\n");
						
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"GridPercentSS",
								LaTeXUtils.numberAsPercent(PRVI25_GridSourceBuilder.getCrustalFractSS()*100d, 0), false)+"\n");
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"GridPercentRev",
								LaTeXUtils.numberAsPercent(PRVI25_GridSourceBuilder.getCrustalFractRev()*100d, 0), false)+"\n");
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"GridPercentNorm",
								LaTeXUtils.numberAsPercent(PRVI25_GridSourceBuilder.getCrustalFractNorm()*100d, 0), false)+"\n");
					}
					if (obs != null) {
						System.out.println("Looking for "+seisReg+" observed MFD exceedances (r="+r+")");
						for (boolean cml : new boolean[] {false,true}) {
							// find the first and last magnitude where the on-fault rate exceeds the observed rate
							EvenlyDiscretizedFunc onFunc = cml ? onFaultMean.getCumRateDistWithOffset() : onFaultMean;
							for (boolean upper : new boolean[] {false,true}) {
								EvenlyDiscretizedFunc obsFunc;
								if (upper) {
									if (cml) {
										System.out.println("Cumulative, Upper");
										obsFunc = obs.getUpper().getCumRateDistWithOffset();
									} else {
										System.out.println("Incremental, Upper");
										obsFunc = obs.getUpper();
									}
								} else {
									if (cml) {
										System.out.println("Cumulative, Preferred");
										obsFunc = obs.getCumRateDistWithOffset();
									} else {
										System.out.println("Incremental, Preferred");
										obsFunc = obs;
									}
								}
								int firstExceedanceIndex = -1;
								int lastExceedanceIndex = -1;
								for (int i=0; i<onFunc.size(); i++) {
									double mag = onFunc.getX(i);
									int obsX = obsFunc.getClosestXIndex(mag);
									if ((float)mag == (float)obsFunc.getX(obsX)) {
										double onFaultRate = onFunc.getY(i);
										double obsRate = obsFunc.getY(i);
										if (onFaultRate > obsRate) {
											if (!cml || onFaultMean.getY(i) > 0)
												System.out.println("M"+(float)mag+"\tfault="+(float)onFaultRate
														+"\texceeds\tobs="+(float)obsRate
														+"\tdiff="+(float)(onFaultRate-obsRate)
														+"\tpDiff="+(float)(100d*(onFaultRate-obsRate)/obsRate));
											if (firstExceedanceIndex < 0)
												firstExceedanceIndex = i;
											lastExceedanceIndex = i;
										}
									} else {
										System.out.println("Mag not matched: "+(float)mag+"\tclosest="+(float)obsFunc.getX(obsX));
									}
								}
								if (firstExceedanceIndex >= 0) {
									// we have an exceedance
									double firstMag = onFunc.getX(firstExceedanceIndex);
									if (!cml)
										// move to the start of the bin
										firstMag -= 0.5*onFunc.getDelta();
									String texMagPrefix = texPrefix;
									if (cml)
										texMagPrefix += "Cml";
									else
										texMagPrefix += "Incr";
									texMagPrefix += "ObsExceed";
									if (upper)
										texMagPrefix += "Upper";
									texFW.write(LaTeXUtils.defineValueCommand(texMagPrefix+"FirstMag",
											RupSetStatsTexWriter.magDF.format(firstMag), false)+"\n");
									if (!cml && lastExceedanceIndex > firstExceedanceIndex
											// make sure we went below, didn't just run out of ruptures
											&& onFaultMean.getCumRate(lastExceedanceIndex) > 0) {
										// it exceeded for multiple magnitudes
										double lastMag = onFunc.getX(lastExceedanceIndex);
										if (!cml)
											// move to the end of the bin
											lastMag += 0.5*onFunc.getDelta();
										texFW.write(LaTeXUtils.defineValueCommand(texMagPrefix+"LastMag",
												RupSetStatsTexWriter.magDF.format(lastMag), false)+"\n");
									}
								}
							}
						}
					}
					List<Double> mags;
					List<String> magNames;
					
					if (trt == TectonicRegionType.ACTIVE_SHALLOW) {
						mags = List.of(6.5, 7d, 7.5);
						magNames = List.of("SixFive", "Seven", "SevenFive");
					} else {
						mags = List.of(8d, 8.5, 9d);
						magNames = List.of("Eight", "EightFive", "Nine");
					}
					for (int m=0; m<mags.size(); m++) {
						double mag = mags.get(m);
						String magName = magNames.get(m);
						// round up as they're incremental and we want to make sure it maps to the bin that lies above
						double cmlRate = onFaultMean.getCumRate(onFaultMean.getClosestXIndex(mag+0.01));
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"M"+magName+"Rate",
								LaTeXUtils.numberExpFormatSigFigs(cmlRate, 3), false)+"\n");
						double cmlRI = 1d/cmlRate;
						texFW.write(LaTeXUtils.defineValueCommand(texPrefix+"M"+magName+"RI",
//								LaTeXUtils.numberExpFormatFixedDecimal(cmlRI, 0), false)+"\n");
								LaTeXUtils.groupedIntNumber(cmlRI), false)+"\n");
					}
				}
				texFW.write("\n");
				
				IncrementalMagFreqDist total = null;
				if (trt != TectonicRegionType.SUBDUCTION_SLAB)
					total = sum(onFaultMean, gridded);
				if (trt != TectonicRegionType.SUBDUCTION_SLAB)
					plot(outputDir, prefix, xRange,
							onFaultDists, onFaultMean, // on fault
							null, gridded, // gridded
							null, total, // total
							obs);
				plot(outputDir, prefix+"_gridded_dist", xRange,
						null, onFaultMean, // on fault
						griddedDists, gridded, // gridded
						null, total, // total
						obs);
				if (trt != TectonicRegionType.SUBDUCTION_SLAB) {
					UncertainBoundedIncrMagFreqDist[] totalDists = getRegionalMFDFractiles(sol, MFDType.SUM, reg, fractiles);
					plot(outputDir, prefix+"_total_dist", xRange,
							null, onFaultMean, // on fault
							null, gridded, // gridded
							totalDists, total, // total
							obs);
				}
			}
		}
		crustalTexFW.close();
		subductionTexFW.close();
	}
	
	private static void plot(File outputDir, String prefix, Range xRange,
			UncertainBoundedIncrMagFreqDist[] onFaultDists, IncrementalMagFreqDist onFaultMean,
			UncertainBoundedIncrMagFreqDist[] griddedDists, IncrementalMagFreqDist gridded,
			UncertainBoundedIncrMagFreqDist[] totalDists, IncrementalMagFreqDist total,
			UncertainBoundedIncrMagFreqDist obs) throws IOException {
		Color onFaultColor = Colors.tab_red;
		Color onFaultTransColor = new Color(onFaultColor.getRed(), onFaultColor.getGreen(), onFaultColor.getBlue(), 60);
		Color obsColor = Colors.tab_green;
		Color griddedColor = Colors.tab_blue;
		Color griddedTransColor = new Color(griddedColor.getRed(), griddedColor.getGreen(), griddedColor.getBlue(), 60);
		Color totalColor = Colors.tab_purple;
		Color totalTransColor = new Color(totalColor.getRed(), totalColor.getGreen(), totalColor.getBlue(), 60);
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		if (obs != null) {
			funcs.add(obs);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, obsColor));
			
			IncrementalMagFreqDist obsLower = obs.getLower();
			obsLower.setName(obs.getBoundName());
			funcs.add(obsLower);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, obsColor));
			
			IncrementalMagFreqDist obsUpper = obs.getUpper();
			obsUpper.setName(null);
			funcs.add(obsUpper);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DOTTED, 2f, obsColor));
		}
		
		if (gridded != null) {
			if (prefix.contains("interface"))
				gridded.setName("Sub-seismogenic (gridded)");
			else
				gridded.setName("Gridded");
			funcs.add(gridded);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, griddedColor));
		}
		
		if (griddedDists != null) {
			for (int i=0; i<griddedDists.length; i++) {
				if (i == 0)
					griddedDists[i].setName(extremaLabel);
				else
					griddedDists[i].setName(null);
				funcs.add(griddedDists[i]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, griddedTransColor));
			}
		}
		
		if (onFaultMean != null) {
			if (prefix.contains("interface"))
				onFaultMean.setName("Supra-seismogenic (inversion)");
			else
				onFaultMean.setName("On-fault supra-seismogenic");
			funcs.add(onFaultMean);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, onFaultColor));
		}
		
		if (onFaultDists != null) {
			for (int i=0; i<onFaultDists.length; i++) {
				if (i == 0)
					onFaultDists[i].setName(fractileLabel);
				else
					onFaultDists[i].setName(null);
				funcs.add(onFaultDists[i]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, onFaultTransColor));
			}
		}
		
		if (total != null) {
			total.setName("Total");
			funcs.add(total);
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, totalColor));
		}
		
		if (totalDists != null) {
			for (int i=0; i<totalDists.length; i++) {
				if (i == 0)
					totalDists[i].setName(fractileLabel);
				else
					totalDists[i].setName(null);
				funcs.add(totalDists[i]);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN, 1f, totalTransColor));
			}
		}
		
//		// again on top
//		if (total == null && griddedDists == null && totalDists == null && onFaultMean != null) {
//			onFaultMean = onFaultMean.deepClone();
//			onFaultMean.setName(null);
//			funcs.add(onFaultMean);
//			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, onFaultColor));
//		}
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
		plot.setLegendInset(true);
		
		Range yRange = new Range(1e-6, 1e1);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
	}
	
	private static void plotMultiSlab(File outputDir, String prefix, Range xRange,
			UncertainBoundedIncrMagFreqDist carDist,
			UncertainBoundedIncrMagFreqDist mueDist,
			IncrementalMagFreqDist slab03) throws IOException {
		Color carColor = Colors.tab_orange;
		Color mueColor = Colors.tab_green;
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		slab03.setName("2003 Model");
		funcs.add(slab03);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		
		mueDist.setName("Muertos Preferred");
		funcs.add(mueDist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, mueColor));
		
		mueDist.getUpper().setName("Muertos High & Low");
		funcs.add(mueDist.getUpper());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, mueColor));
		
		mueDist.getLower().setName(null);
		funcs.add(mueDist.getLower());
		chars.add(chars.get(chars.size()-1));
		
		carDist.setName("Caribbean Preferred");
		funcs.add(carDist);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, carColor));
		
		carDist.getUpper().setName("Caribbean High & Low");
		funcs.add(carDist.getUpper());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, carColor));
		
		carDist.getLower().setName(null);
		funcs.add(carDist.getLower());
		chars.add(chars.get(chars.size()-1));
		
		PlotSpec plot = new PlotSpec(funcs, chars, " ", "Magnitude", "Incremental Rate (1/yr)");
		plot.setLegendInset(true);
		
		Range yRange = new Range(1e-6, 1e1);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(plot, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, prefix, gp, 800, 750, true, true, false);
	}
	
	private static IncrementalMagFreqDist loadSlab03(EvenlyDiscretizedFunc refMFD) {
		NshmErf erf = new NshmErf(Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-prvi-2003-main"),
				Set.of(TectonicRegionType.SUBDUCTION_SLAB), IncludeBackgroundOption.INCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		IncrementalMagFreqDist mfd = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		
		for (ProbEqkSource source : erf) {
//			System.out.println("Source: "+source.getName()+" ("+source.getClass().getName()+")");
			for (ProbEqkRupture rup : source) {
				double rate = rup.getMeanAnnualRate(1d);
				double mag = rup.getMag();
//				System.out.println("M"+(float)mag+", "+(float)rate);
				mfd.add(mfd.getClosestXIndex(mag), rate);
			}
		}
		
		Preconditions.checkState(mfd.calcSumOfY_Vals() > 0d);
		
		// set to NaN below Mmin to prevent a vertical line at the start of the dist
		for (int i=0; i<mfd.size(); i++) {
			if (mfd.getY(i) == 0d)
				mfd.set(i, Double.NaN);
			else
				break;
		}
		
		return mfd;
	}

	private static double[] fractiles = {0d, 0.025, 0.16, 0.84, 0.975d, 1d};
	private static String fractileLabel = "p[0,2.5,16,84,97.5,100]";
	private static double[] extrema = {0d, 1d};
	private static String extremaLabel = "Extrema";
	
	private static UncertainBoundedIncrMagFreqDist[] getOnFaultMFDFractiles(FaultSystemSolution sol, Region reg) {
		BranchSectNuclMFDs branchSectMFDs = sol.requireModule(BranchSectNuclMFDs.class);
		
		double[] sectFracts = null;
		if (reg != null)
			sectFracts = sol.getRupSet().getFractSectsInsideRegion(reg, false);
		
		IncrementalMagFreqDist[] fractileMFDs = branchSectMFDs.calcIncrementalFractiles(sectFracts, fractiles);
		
		int numRet = fractiles.length/2;
		UncertainBoundedIncrMagFreqDist[] ret = new UncertainBoundedIncrMagFreqDist[numRet];
		for (int i=0; i<numRet; i++) {
			IncrementalMagFreqDist lower = fractileMFDs[i];
			IncrementalMagFreqDist upper = fractileMFDs[fractileMFDs.length - (1 + i)];
			ret[i] = new UncertainBoundedIncrMagFreqDist(CombinedMFDsPlot.average(lower, upper), lower, upper, null);
		}
		
		return ret;
	}
	
	private static UncertainBoundedIncrMagFreqDist[] getRegionalMFDFractiles(FaultSystemSolution sol,
			MFDType type, Region reg, double[] fractiles) throws IOException {
		BranchRegionalMFDs regMFDs = sol.requireModule(BranchRegionalMFDs.class);
		
		IncrementalMagFreqDist[] fractileMFDs;
		if (reg != null) {
			RegionsOfInterest roi = sol.getRupSet().requireModule(RegionsOfInterest.class);
			List<Region> regions = roi.getRegions();
			Preconditions.checkState(regions.size() == regMFDs.getNumRegions());
			int regionIndex = -1;
			for (int i=0; i<regions.size(); i++) {
				if (regions.get(i).equalsRegion(reg)) {
					regionIndex = i;
					break;
				}
			}
			Preconditions.checkState(regionIndex >= 0, "Didn't find region named %s in ROI, have %s regions", reg.getName(), regions.size());
			fractileMFDs = regMFDs.calcRegionalIncrementalFractiles(type, regionIndex, fractiles);
		} else {
			fractileMFDs = regMFDs.calcTotalIncrementalFractiles(type, fractiles);
		}
		
		int numRet = fractiles.length/2;
		UncertainBoundedIncrMagFreqDist[] ret = new UncertainBoundedIncrMagFreqDist[numRet];
		for (int i=0; i<numRet; i++) {
			IncrementalMagFreqDist lower = fractileMFDs[i];
			IncrementalMagFreqDist upper = fractileMFDs[fractileMFDs.length - (1 + i)];
			ret[i] = new UncertainBoundedIncrMagFreqDist(CombinedMFDsPlot.average(lower, upper), lower, upper, null);
		}
		
		return ret;
	}
	
	private static IncrementalMagFreqDist sum(IncrementalMagFreqDist mfd1, IncrementalMagFreqDist mfd2) {
		Preconditions.checkState(mfd1.getDelta() == mfd2.getDelta());
		double min = Math.min(mfd1.getMinX(), mfd2.getMinX());
		double max = Math.max(mfd1.getMaxX(), mfd2.getMaxX());
		int num = (int)((max - min)/mfd1.getDelta() + 1);
		IncrementalMagFreqDist ret = new IncrementalMagFreqDist(min, max, num);
		Preconditions.checkState(mfd1.getDelta() == ret.getDelta());
		for (Point2D pt : mfd1)
			ret.set(ret.getClosestXIndex(pt.getX()), pt.getY());
		for (Point2D pt : mfd2)
			ret.add(ret.getClosestXIndex(pt.getX()), pt.getY());
		return ret;
	}
	
	private static UncertainBoundedIncrMagFreqDist getMmaxAveragedSlab(EvenlyDiscretizedFunc refMFD,
			PRVI25_SeismicityRegions seisReg, boolean averageRate) throws IOException {
		List<UncertainBoundedIncrMagFreqDist> slabMFDs = new ArrayList<>();
		List<Double> slabWeights = new ArrayList<>();
		for (PRVI25_SeismicityRateEpoch epoch : PRVI25_SeismicityRateEpoch.values()) {
			double epochWeight = epoch.getNodeWeight(null);
			if (epochWeight == 0d)
				continue;
			for (PRVI25_SubductionSlabMMax slabMmax : PRVI25_SubductionSlabMMax.values()) {
				double mMaxWeight = slabMmax.getNodeWeight(null);
				if (mMaxWeight == 0d)
					continue;
				IncrementalMagFreqDist avg;
				SeismicityRateModel siesModel;
				if (seisReg == PRVI25_SeismicityRegions.CAR_INTRASLAB) {
					if (averageRate)
						avg = PRVI25_SubductionCaribbeanSeismicityRate.AVERAGE.build(epoch, refMFD,
								slabMmax.getIncrementalMmax(), PRVI25_GridSourceBuilder.SLAB_M_CORNER, true);
					else
						avg = PRVI25_SubductionCaribbeanSeismicityRate.PREFFERRED.build(epoch, refMFD,
								slabMmax.getIncrementalMmax(), PRVI25_GridSourceBuilder.SLAB_M_CORNER, true);
					siesModel = PRVI25_SubductionCaribbeanSeismicityRate.loadRateModel(epoch, true);
				} else if (seisReg == PRVI25_SeismicityRegions.MUE_INTRASLAB) {
					if (averageRate)
						avg = PRVI25_SubductionMuertosSeismicityRate.AVERAGE.build(epoch, refMFD,
								slabMmax.getIncrementalMmax(), PRVI25_GridSourceBuilder.SLAB_M_CORNER, true);
					else
						avg = PRVI25_SubductionMuertosSeismicityRate.PREFFERRED.build(epoch, refMFD,
								slabMmax.getIncrementalMmax(), PRVI25_GridSourceBuilder.SLAB_M_CORNER, true);
					siesModel = PRVI25_SubductionMuertosSeismicityRate.loadRateModel(epoch, true);
				} else {
					throw new IllegalStateException();
				}
				slabMFDs.add(siesModel.getBounded(refMFD, slabMmax.getIncrementalMmax()));
				slabWeights.add(epochWeight*mMaxWeight);
			}
		}
		return PRVI25_SeismicityRateEpoch.averageUncert(slabMFDs, slabWeights);
	}
	
	static IncrementalMagFreqDist calcGriddedMFD(Region region, TectonicRegionType trt,
			FaultSystemSolution sol, EvenlyDiscretizedFunc refMFD) throws IOException {
		GriddedRegion gridReg = new GriddedRegion(region, 0.1, GriddedRegion.ANCHOR_0_0);
		GridSourceList gridProv = sol.requireModule(GridSourceList.class);
		SummedMagFreqDist ret = new SummedMagFreqDist(refMFD.getMinX(), refMFD.getMaxX(), refMFD.size());
		for (int i=0; i<gridProv.getNumLocations(); i++) {
			int regIndex = gridReg.indexForLocation(gridProv.getLocation(i));
			if (regIndex >= 0) {
				IncrementalMagFreqDist mfd = gridProv.getMFD(trt, i);
				ret.addIncrementalMagFreqDist(mfd);
			}
		}
		return ret;
	}

}
