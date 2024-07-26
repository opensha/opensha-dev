package scratch.kevin.ucerf3;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.Precision;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.calc.GaussianDistCalc;
import org.opensha.commons.calc.WeightedSampler;
import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.Somerville_2006_MagAreaRel;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.comcat.ComcatAccessor;
import org.opensha.commons.data.comcat.ComcatRegionAdapter;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4i26_M01_TaperBasinDepth;
import org.opensha.commons.data.siteData.impl.ThompsonVs30_2020;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.geo.json.FeatureCollection;
import org.opensha.commons.geo.json.FeatureProperties;
import org.opensha.commons.geo.json.GeoJSON_Type;
import org.opensha.commons.geo.json.Geometry;
import org.opensha.commons.geo.json.Geometry.DepthSerializationType;
import org.opensha.commons.geo.json.Geometry.MultiPolygon;
import org.opensha.commons.geo.json.Geometry.Polygon;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeLevel.RandomlySampledLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.logicTree.LogicTreeNode.RandomlySampledNode;
import org.opensha.commons.mapping.PoliticalBoundariesData;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.IDPairing;
import org.opensha.commons.util.Interpolate;
import org.opensha.commons.util.MarkdownUtils;
import org.opensha.commons.util.MarkdownUtils.TableBuilder;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.cpt.CPTVal;
import org.opensha.refFaultParamDb.vo.DeformationModelSummary;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.AbstractNthRupERF;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RupSetDeformationModel;
import org.opensha.sha.earthquake.faultSysSolution.RupSetFaultModel;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchAveragingOrder;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchParentSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchRegionalMFDs.MFDType;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectBVals;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectNuclMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.BranchSectParticMFDs;
import org.opensha.sha.earthquake.faultSysSolution.modules.ClusterRuptures;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModSectMinMags;
import org.opensha.sha.earthquake.faultSysSolution.modules.PolygonFaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.RupMFDsModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.SolutionLogicTree;
import org.opensha.sha.earthquake.faultSysSolution.modules.TrueMeanRuptureMappings;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.ClusterRupture;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.FaultSubsectionCluster;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.PlausibilityResult;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.AspectRatioFilter;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.plausibility.impl.prob.JumpProbabilityCalc.HardcodedJumpProb;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.UniqueRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.faultSysSolution.util.SubSectionBuilder;
import org.opensha.sha.earthquake.observedEarthquake.ObsEqkRupList;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BackgroundRupParam;
import org.opensha.sha.earthquake.param.BackgroundRupType;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.earthquake.rupForecastImpl.PointSource13b.PointSurface13b;
import org.opensha.sha.earthquake.rupForecastImpl.PointSourceNshm.PointSurfaceNshm;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.FaultSegmentData;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.A_FaultsFetcher;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.data.finalReferenceFaultParamDb.DeformationModelSummaryFinal;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_LogicTreeBranch;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_ScalingRelationships;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SegmentationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_SingleStates;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.random.BranchSamplingManager;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.rupForecastImpl.prvi25.util.PRVI25_RegionLoader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.simulators.RSQSimEvent;
import org.opensha.sha.simulators.srf.RSQSimStateTime;
import org.opensha.sha.simulators.srf.RSQSimStateTransitionFileReader;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Table;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import gov.usgs.earthquake.nshmp.mfd.Mfd.Properties.GutenbergRichter;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import gov.usgs.earthquake.nshmp.model.NshmSurface;
import net.mahdilamb.colormap.Colors;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.FaultSystemSolutionERF_ETAS;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Config;
import scratch.UCERF3.erf.ETAS.launcher.ETAS_Launcher;
import scratch.UCERF3.erf.ETAS.launcher.TriggerRupture;
import scratch.UCERF3.erf.ETAS.launcher.TriggerRupture.EdgeFault;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.kevin.nshm23.dmCovarianceTests.RandomDefModSampleLevel;
import scratch.kevin.pointSources.InvCDF_RJBCorrPointSurface;
import scratch.kevin.simulators.RSQSimCatalog;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;
import scratch.kevin.simulators.ruptures.RSQSimGeographicMapMaker;
import scratch.kevin.ucerf3.eal.LossCOV_Model;

public class PureScratch {
	
	private static final void test250() throws IOException {
//		List<LogicTreeLevel<? extends LogicTreeNode>> levels = List.of(
//				NSHM23_LogicTreeBranch.DM, NSHM23_LogicTreeBranch.SCALE,
//				NSHM23_LogicTreeBranch.PALEO_UNCERT, NSHM23_LogicTreeBranch.SEG, NSHM23_LogicTreeBranch.SUPRA_B);
//		SupraSeisBValues.SEG_DEPENDENT_WEIGHTS = true;
//		LogicTree<?> tree = LogicTree.buildExhaustive(levels, true);
//		
//		BranchWeightProvider currentWeights = new BranchWeightProvider.CurrentWeights();
//		BranchWeightProvider topDownWeights = new BranchWeightProvider.TopDownWeights(tree);
//		
//		double sumCurrent = 0d;
//		double sumTopDown = 0d;
//		double[] currentWeightsArray = new double[tree.size()];
//		double[] topDownWeightsArray = new double[tree.size()];
//		for (int i=0; i<tree.size(); i++) {
//			LogicTreeBranch<?> branch = tree.getBranch(i);
//			
//			currentWeightsArray[i] = currentWeights.getWeight(branch);
//			sumCurrent += currentWeightsArray[i];
//			topDownWeightsArray[i] = topDownWeights.getWeight(branch);
//			sumTopDown += topDownWeightsArray[i];
//		}
//		
//		// make sure they're equal right now
//		int numEqual = 0;
//		for (int i=0; i<tree.size(); i++) {
//			LogicTreeBranch<?> branch = tree.getBranch(i);
//			
//			double currentWeight = currentWeightsArray[i]/sumCurrent;
//			double topDownWeight = topDownWeightsArray[i]/sumTopDown;
//			System.out.println(i+". current="+(float)currentWeight+"\ttopDown="+(float)topDownWeight+"; \t"+branch);
//			if ((float)currentWeight == (float)topDownWeight)
//				numEqual++;
//		}
//		
//		System.out.println("Original sumCurrent="+(float)sumCurrent+", sumTopDown="+(float)sumTopDown);
//		System.out.println(numEqual+"/"+tree.size()+" weights are equal");
//		
//		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
//				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/");
//		LogicTree<?> origTree = LogicTree.read(new File(dir, "logic_tree.json"));
//		
//		// modify weights
//		Map<LogicTreeBranch<?>, Double> modWeightsMap = new HashMap<>();
//		// initialize with zeros, then override with matches
//		for (LogicTreeBranch<?> branch : origTree)
//			modWeightsMap.put(branch, 0d);
//		// map in modified weights for branches that exist in the new tree
//		for (LogicTreeBranch<?> branch : tree) {
//			// find corresponding branch in original tree
//			LogicTreeBranch<?> matchingBranch = null;
//			for (LogicTreeBranch<?> origBranch : origTree) {
//				boolean match = true;
//				for (LogicTreeNode node : branch) {
//					if (!origBranch.hasValue(node)) {
//						match = false;
//						break;
//					}
//				}
//				if (match) {
//					matchingBranch = origBranch;
//					break;
//				}
//			}
//			Preconditions.checkNotNull(matchingBranch, "No match found for branch %s", branch);
//			double modWeight = topDownWeights.getWeight(branch);
//			modWeightsMap.put(matchingBranch, modWeight);
//		}
//		for (LogicTreeBranch<?> branch : modWeightsMap.keySet()) {
//			double weight = modWeightsMap.get(branch);
//			branch.setOrigBranchWeight(weight);
//		}
//		
//		origTree.setWeightProvider(new BranchWeightProvider.OriginalWeights());
//		origTree.write(new File(dir, "logic_tree_seg_specific_b_val_weights.json"));
	}

	private static void test251() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/"
				+ "FM3_1_branch_averaged.zip"));
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		int rupID = 251641;
		List<FaultSection> sects = rupSet.getFaultSectionDataForRupture(rupID);
		System.out.println("Sections for rup "+rupID);
		for (FaultSection sect : sects) {
			System.out.println("\t"+sect.getSectionId()+". "+sect.getName());
			System.out.println("\t\tDip: "+(float)sect.getAveDip());
			System.out.println("\t\tUpper depth: orig="+(float)sect.getOrigAveUpperDepth()
				+", reduced="+(float)sect.getReducedAveUpperDepth());
			System.out.println("\t\tLower depth: "+(float)sect.getAveLowerDepth());
			System.out.println("\t\tDDW: orig="+(float)sect.getOrigDownDipWidth()
				+", reduced="+(float)sect.getReducedDownDipWidth());
			System.out.println("\t\tLength: "+(float)sect.getTraceLength());
			System.out.println("\t\tArea: orig="+(float)(sect.getArea(false)*1e-6)
				+", reduced="+(float)(sect.getArea(true)*1e-6));
		}
		System.out.println("Area: "+rupSet.getAreaForRup(rupID)*1e-6);
	}

	private static void test252() throws IOException {
		File dir = new File("/data/kevin/simulators/catalogs/bruce/rundir5413_multifault_separate");
		RSQSimCatalog catalog = new RSQSimCatalog(dir, "Test Catalog", FaultModels.FM3_1, DeformationModels.GEOLOGIC);
		RSQSimEvent event = catalog.loader().byID(202);
		System.out.println("Magnitude: "+event.getMagnitude());
		RSQSimStateTransitionFileReader transReader = catalog.getTransitions();
		List<RSQSimStateTime> transitions = new ArrayList<>();
		transReader.getTransitions(event, transitions);
		System.out.println("Read "+transitions.size()+" transitions for event "+event.getID()+":");
		for (RSQSimStateTime trans : transitions)
			System.out.println(trans);
		
		List<RSQSimEvent> withTrans = catalog.loader().hasTransitions().minMag(5d).load();
		System.out.println("Loaded "+withTrans.size()+" events with valid transitions");
		for (RSQSimEvent e2 : withTrans)
			System.out.println(e2.getID()+" is a M"+e2.getMagnitude()+" at "+e2.getTimeInYears()+" yrs");
	}

	private static void test253() throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR-ba_only");
		File treeFile = new File(dir, "logic_tree-nga_w2s-gmm_add_epi.json");
		
		File resultsFile = new File(dir, "results.zip");
		
		SolutionLogicTree slt = SolutionLogicTree.load(resultsFile, LogicTree.read(treeFile));
		
		LogicTree<?> tree = slt.getLogicTree();
		System.out.println("Tree has "+tree.size()+" branches");
		
		for (LogicTreeBranch<?> branch : tree) {
			System.out.println("Branch: "+branch);
			FaultSystemSolution sol = slt.forBranch(branch);
			System.out.println("\tGrid sources? "+sol.hasModule(GridSourceProvider.class));
		}
	}
	
	private static void test254() throws IOException {
//		RupSetFaultModel[] fms = { FaultModels.FM3_1, FaultModels.FM3_2 };
////		RupSetFaultModel[] fms = { FaultModels.FM3_1 };
//		RupSetDeformationModel[] dms = DeformationModels.values();
//		int parentID = 117;
		
		RupSetFaultModel[] fms = { NSHM23_FaultModels.WUS_FM_v2 };
		RupSetDeformationModel[] dms = NSHM23_DeformationModels.values();
		int parentID = 333;
		
		List<String> slipLines = new ArrayList<>();
		List<String> momentLines = new ArrayList<>();
		
		DecimalFormat df = new DecimalFormat("0.0");
		
		double totSumArea = 0d;
		double totSumWidthArea = 0d;
		
		for (RupSetFaultModel fm : fms) {
			for (RupSetDeformationModel dm : dms) {
				if (dm.getNodeWeight(null) == 0d)
					continue;
				String line = fms.length == 1 ? "" : fm.getShortName()+", ";
				List<? extends FaultSection> subSects = dm.build(fm);
				double sumSlipArea = 0d;
				double sumArea = 0d;
				double sumMomentRate = 0d;
				
				for (FaultSection sect : subSects) {
					if (sect.getParentSectionId() == parentID) {
						double slip = sect.getOrigAveSlipRate();
						double area = sect.getArea(false);
						sumSlipArea += slip*area;
						sumArea += area;
						double momentRate = FaultMomentCalc.getMoment(area, slip*1e-3);
						sumMomentRate += momentRate;
						
						double ddw = sect.getOrigDownDipWidth();
						totSumArea += area;
						totSumWidthArea += ddw*area;
					}
				}
				
				double avgSlip = sumSlipArea/sumArea;
				
				line += dm.getShortName()+":\t";
				
				slipLines.add(line+df.format(avgSlip)+" mm/yr");
				momentLines.add(line+(float)sumMomentRate+" N-m/yr");
			}
		}
		
		for (String line : slipLines)
			System.out.println(line);
		for (String line : momentLines)
			System.out.println(line);
		
		double avgDDW = totSumWidthArea/totSumArea;
		System.out.println("Avg DDW:\t"+df.format(avgDDW)+" km");
	}
	
	private static void test255() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF4/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		int parentID = 704;
		double minMag = 7.5d;
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		List<Integer> rups = rupSet.getRupturesForParentSection(parentID);
		
		ScalarIMR gmm = AttenRelRef.ASK_2014.get();
		Site site = new Site(new Location(34.01933247362704, -118.28641763873607));
		site.addParameterList(gmm.getSiteParams());
		ThompsonVs30_2020 vs30Model = new ThompsonVs30_2020();
		CVM4i26_M01_TaperBasinDepth z10Model = new CVM4i26_M01_TaperBasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
		
		double vs30 = vs30Model.getValue(site.getLocation());
		double z10 = z10Model.getValue(site.getLocation());
		System.out.println("vs30="+vs30+"\tz10="+z10);
		site.getParameter(Vs30_Param.NAME).setValue(vs30);
		site.getParameter(DepthTo1pt0kmPerSecParam.NAME).setValue(z10*1e3); // km -> m
		gmm.setSite(site);
		gmm.setIntensityMeasure(PGA_Param.NAME);
		
		DiscretizedFunc avgExceedProbs = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : avgExceedProbs)
			logXVals.set(Math.log(pt.getX()), 1d);
		for (int i=0; i<avgExceedProbs.size(); i++)
			avgExceedProbs.set(i, 0d);
		
		double rateSum = 0d;
		double avgMag = 0d;
		for (int rupIndex : rups) {
			double mag = rupSet.getMagForRup(rupIndex);
			if (mag >= minMag) {
				double rate = sol.getRateForRup(rupIndex);
				rateSum += rate;
				avgMag += rate*mag;
				EqkRupture rup = new EqkRupture(mag, rupSet.getAveRakeForRup(rupIndex), rupSet.getSurfaceForRupture(parentID, 1d), null);
				gmm.setEqkRupture(rup);
				gmm.getExceedProbabilities(logXVals);
				
				for (int i=0; i<avgExceedProbs.size(); i++)
					avgExceedProbs.set(i, avgExceedProbs.getY(i) + logXVals.getY(i)*rate);
			}
		}
		avgMag /= rateSum;
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		erf.getTimeSpan().setDuration(30d);
		erf.updateForecast();
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		calc.getHazardCurve(logXVals, site, gmm, erf);
		DiscretizedFunc fullCurve = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : logXVals)
			fullCurve.set(Math.exp(pt.getX()), pt.getY());
		
		System.out.println(logXVals);
		
		System.out.println("Average mag: "+(float)avgMag);
		
		double[] gmLevels = {0.01, 0.05, 0.1, 0.5, 1};
		
		avgExceedProbs.scale(1d/rateSum);
		double prob30 = 1d-Math.exp(-rateSum*30d);
		System.out.println("Sum rate: "+(float)rateSum);
		System.out.println("30 year probability: "+(float)prob30);
		System.out.println("Exceed probs:");
		DecimalFormat pDF = new DecimalFormat("0.0##%");
		for (double gmLevel : gmLevels) {
			double prob = avgExceedProbs.getInterpolatedY_inLogXLogYDomain(gmLevel);
			System.out.println("\t"+(float)gmLevel+":\t"+(float)prob+" ("+pDF.format(prob)+")");
		}
		
		System.out.println("Full ERF probs:");
		for (double gmLevel : gmLevels) {
			double prob = fullCurve.getInterpolatedY_inLogXLogYDomain(gmLevel);
			System.out.println("\t"+(float)gmLevel+":\t"+(float)prob+" ("+pDF.format(prob)+")");
		}
	}
	
	private static void test256() throws IOException {
		File dir = new File("/data/kevin/nshm23/batch_inversions/"
				+ "2023_06_23-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/");
		FaultSystemSolution sol = FaultSystemSolution.load(new File(dir, "true_mean_solution.zip"));
		
		System.out.println("Solution has "+sol.getRupSet().getNumSections()+" sections");
		System.out.println("Solution has "+sol.getRupSet().getNumRuptures()+" ruptures");
		
		SolutionLogicTree slt = SolutionLogicTree.load(new File(dir, "results.zip"));
		
		TrueMeanRuptureMappings mappings = sol.getRupSet().requireModule(TrueMeanRuptureMappings.class);
		
		double[] origRates = sol.getRateForAllRups();
		double[] reconstructedRates = new double[origRates.length];
		
		LogicTree<?> tree = slt.getLogicTree();
		double weightSum = 0d;
		for (LogicTreeBranch<?> branch : tree)
			weightSum += tree.getBranchWeight(branch);
		
		for (LogicTreeBranch<?> branch : tree) {
			mappings.getSectionMappings(branch);
			int[] rupMappings = mappings.getRuptureMappings(branch);
			double[] myRates = slt.loadRatesForBranch(branch);
			Preconditions.checkState(rupMappings.length == myRates.length);
			
			double weight = tree.getBranchWeight(branch)/weightSum;
			
			for (int r=0; r<myRates.length; r++)
				if (myRates[r] > 0)
					reconstructedRates[rupMappings[r]] += weight*myRates[r];
		}
		
		MinMaxAveTracker absDiffTrack = new MinMaxAveTracker();
		MinMaxAveTracker diffTrack = new MinMaxAveTracker();
		MinMaxAveTracker pDiffTrack = new MinMaxAveTracker();
		
		double origSum = StatUtils.sum(origRates);
		double reconstructedSum = StatUtils.sum(reconstructedRates);
		
		for (int r=0; r<origRates.length; r++) {
			double diff = reconstructedRates[r] - origRates[r];
			diffTrack.addValue(diff);
			absDiffTrack.addValue(Math.abs(diff));
			if (origRates[r] > 0)
				pDiffTrack.addValue(100d*(reconstructedRates[r] - origRates[r])/origRates[r]);
		}
		System.out.println("Audit stats:");
		System.out.println("\tReconstructed Sum:\t"+(float)reconstructedSum);
		System.out.println("\tOriginal Sum:\t"+(float)origSum);
		System.out.println("\tDiff stats:\t"+diffTrack);
		System.out.println("\t|Diff| stats:\t"+absDiffTrack);
		System.out.println("\t%Diff stats:\t"+pDiffTrack);
	}
	
	private static void test257() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/tmp/solution.zip"));
		LogicTreeBranch<?> branch = sol.getRupSet().requireModule(LogicTreeBranch.class);
		
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(8.45);
		IncrementalMagFreqDist mfd = sol.calcTotalNucleationMFD(refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta());

		File invsDir = new File("/home/kevin/OpenSHA/UCERF4/batch_inversions");
		File modelDir = new File(invsDir,
				"2023_09_01-nshm23_branches-mod_pitas_ddw-NSHM23_v2-CoulombRupSet-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR");
		
		FaultSystemSolution baSol = FaultSystemSolution.load(new File(modelDir,
				"results_NSHM23_v2_CoulombRupSet_branch_averaged.zip"));
		LogicTree<?> tree = LogicTree.read(new File(modelDir, "logic_tree.json"));
		
		double minMag = 6d;
		
		BranchAveragingOrder order = baSol.requireModule(BranchAveragingOrder.class);
		BranchRegionalMFDs regMFDs = baSol.requireModule(BranchRegionalMFDs.class);
		
		int index = order.getBranchAveragingIndex(branch);
		IncrementalMagFreqDist regMFD = regMFDs.getTotalBranchMFDs(MFDType.SUPRA_ONLY)[index];
		
		System.out.println(branch+" is index "+index);
		System.out.println("MFD directly from solution:\n\n"+mfd);
		System.out.println("MFD from BranchRegionalMFDs:\n\n"+regMFD);
	}
	
	private static void test258() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/modular/"
				+ "FM3_1_SpatSeisU3_branch_averaged_full_modules.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		PolygonFaultGridAssociations polyAssoc = rupSet.getModule(PolygonFaultGridAssociations.class);
		System.out.println("Poly assoc type: "+polyAssoc.getClass().getName());
	}
	
	private static void test259() throws IOException {
		File solFile = new File("/data/kevin/git/ucerf3-etas-launcher/inputs/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(solFile);
		
		sol.removeModuleInstances(RupMFDsModule.class);
		sol.removeModuleInstances(BranchAveragingOrder.class);
		sol.removeModuleInstances(BranchRegionalMFDs.class);
		sol.removeModuleInstances(BranchSectNuclMFDs.class);
		sol.removeModuleInstances(BranchSectParticMFDs.class);
		sol.removeModuleInstances(BranchParentSectParticMFDs.class);
		sol.removeModuleInstances(BranchSectBVals.class);
		sol.getArchive().write(solFile, false); // don't copy extra source files
	}
	
	private static void test260() throws IOException {
		ETAS_Config config = ETAS_Config.readJSON(new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2023_10_19-ComCatM7p1_ci38457511_ShakeMapSurfaces/config.json"));
		
		for (TriggerRupture trigger : config.getTriggerRuptures()) {
			if (trigger instanceof EdgeFault) {
				EdgeFault edge = (EdgeFault)trigger;
				System.out.println("M"+(float)edge.mag+" is an EdgeFault");
				List<LocationList> outlines = new ArrayList<>();
				for (LocationList outline : edge.outlines)
					outlines.add(outline);
				Geometry geom = new Geometry.MultiLineString(outlines);
				Feature feature = new Feature(geom, null);
				System.out.println("\n\n"+feature.toJSON()+"\n");
			}
		}
	}
	
	private static void test261() throws IOException {
		FeatureCollection features = FeatureCollection.read(new File("/tmp/rupture.json"));
		for (Feature feature : features) {
			System.out.println("Feature with ID "+feature.id);
			Geometry geom = feature.geometry;
			System.out.println("Geometry is a "+geom.type);
			if (geom.type == GeoJSON_Type.MultiPolygon) {
				MultiPolygon multiPoly = (MultiPolygon)geom;
				System.out.println("MultiPolygon has "+multiPoly.polygons.size()+" polys");
				for (Polygon poly : multiPoly.polygons) {
					System.out.println("**********");
					System.out.println(poly.toJSON());
					
					System.out.println("**********");
				}
			}
		}
	}
	
	private static void test262() throws IOException {
		File jsonFile = new File("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.0.0/"
				+ "stable-crust/fault/OK/Meers/features/Meers.geojson");
		Feature feature = Feature.read(jsonFile);
		
		GeoJSONFaultSection sect = GeoJSONFaultSection.fromNSHMP_HazFeature(feature);
		System.out.println(sect.toFeature().toJSON());
	}
	
	private static void test263() throws IOException {
		LogicTree<LogicTreeNode> tree = LogicTree.read(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_11_05-nshm23_branches-dm_sampling-NSHM23_v2-CoulombRupSet-GEOLOGIC-DsrUni-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "logic_tree.json"));
		ImmutableList<LogicTreeLevel<? extends LogicTreeNode>> levels = tree.getLevels();
		for (int l=0; l<levels.size(); l++) {
			LogicTreeLevel<?> level = levels.get(l);
			System.out.println("Level "+l+": "+level.getName());
			if (level.affects(FaultSystemRupSet.RUP_SECTS_FILE_NAME, true)) {
				System.out.println("\tIt affects rup sects!");
				System.out.println("\t\tAffected: "+level.getAffected());
				System.out.println("\t\tNot affected: "+level.getNotAffected());
			}
			if (level instanceof RandomlySampledLevel<?>) {
				System.out.println("\tIt's a random level: "+level.getName());
				RandomlySampledLevel<? extends RandomlySampledNode> randLevel = (RandomlySampledLevel<?>)level;
				List<? extends RandomlySampledNode> nodes = randLevel.getNodes();
				System.out.println("\t\tNode count: "+nodes.size());
				System.out.println("\t\tNode 0:\tseed="+nodes.get(0).getSeed()+"; weight="+nodes.get(0).getNodeWeight(null));
				int lastIndex = nodes.size()-1;
				System.out.println("\t\tNode "+lastIndex+":\tseed="+nodes.get(lastIndex).getSeed()+"; weight="+nodes.get(lastIndex).getNodeWeight(null));
			}
		}
		
		int[] indexes = {0,tree.size()-1};
		
		for (int index : indexes) {
			LogicTreeBranch<?> branch = tree.getBranch(index);
			System.out.println("Branch "+index+": "+branch);
			for (LogicTreeNode val : branch) {
				if (val instanceof RandomlySampledNode) {
					RandomlySampledNode randNode = (RandomlySampledNode)val;
					System.out.println("\t"+val.getName()+":\t"+randNode.getSeed()+"; weight="+(float)randNode.getNodeWeight(null));
				}
				
			}
		}
		
		RandomDefModSampleLevel level = new RandomDefModSampleLevel();
		level.buildNodes(new Random(), 10);
		System.out.println("Test level: "+level.getName());
		System.out.println("\tAffected: "+level.getAffected());
		System.out.println("\tNot affected: "+level.getNotAffected());
	}
	
	private static void test264() throws IOException {
		int[] testArray = new int[2];
		int index = 0;
		while (index < testArray.length) {
			testArray[index++] = index;
		}
		for (int i=0; i<testArray.length; i++)
			System.out.println("testArray["+i+"] = "+testArray[i]);
	}
	
	private static void test265() throws IOException {
		Random r = new Random();
		long[] seeds = {
				
				0l,
				0l,
				
//				100l,
//				-100l,
				
//				r.nextLong(),
//				r.nextLong(),
//				r.nextLong(),
//				r.nextLong(),
//				r.nextLong(),
//				r.nextLong(),
//				r.nextLong(),
		};
		System.out.println("Seeds:");
		for (long seed : seeds)
			System.out.println("\t"+seed);
		long combined = BranchSamplingManager.uniqueSeedCombination(seeds);
		System.out.println("Combined seed: "+combined);
	}
	
	private static void test266() throws IOException {
		int numTrials = 100000;
		Random rand = new Random();
		
		List<NSHM23_SegmentationModels> nodes = new ArrayList<>();
		List<Double> weights = new ArrayList<>();
		
		List<NSHM23_SegmentationModels> allowedNodes = new ArrayList<>();
		List<Double> allowedWeights = new ArrayList<>();
		
		for (NSHM23_SegmentationModels node : NSHM23_SegmentationModels.values()) {
			double weight = node.getNodeWeight(null);
			if (weight > 0d) {
				nodes.add(node);
				weights.add(weight);
				
				if (node.isIncludeRupturesThroughCreepingSect()) {
					allowedNodes.add(node);
					allowedWeights.add(weight);
				}
			}
		}
		
		WeightedSampler<NSHM23_SegmentationModels> sampler = new WeightedSampler<>(nodes, weights, rand);
		WeightedSampler<NSHM23_SegmentationModels> allowedSampler = new WeightedSampler<>(allowedNodes, allowedWeights, rand);
		
		int numAllowedEither = 0;
		int numAllowedSantaCruz = 0;
		int numAllowedCalaveras = 0;
		int numAllowedBoth = 0;
		double sumFractEither = 0d;
		double sumFractSantaCruz = 0d;
		double sumFractCalaveras = 0d;
		
		for (int t=0; t<numTrials; t++) {
			NSHM23_SegmentationModels allowedModel = sampler.nextItem();
			boolean allowed = allowedModel.isIncludeRupturesThroughCreepingSect();
			if (allowed) {
				// now test individual jumps
				
				// fully uncorrelated
//				NSHM23_SegmentationModels parkModel = sampler.nextItem();
//				NSHM23_SegmentationModels santaCruzModel = sampler.nextItem();
//				NSHM23_SegmentationModels calaverasModel = sampler.nextItem();
				
				// fully correlated
//				NSHM23_SegmentationModels parkModel = allowedModel;
//				NSHM23_SegmentationModels santaCruzModel = allowedModel;
//				NSHM23_SegmentationModels calaverasModel = allowedModel;
				
				// correlated with allowed nature but not random within
				NSHM23_SegmentationModels parkModel = allowedSampler.nextItem();
				NSHM23_SegmentationModels santaCruzModel = allowedSampler.nextItem();
				NSHM23_SegmentationModels calaverasModel = allowedSampler.nextItem();
				
				boolean allowedPark = parkModel != NSHM23_SegmentationModels.CLASSIC;
				boolean allowedSantaCruz = santaCruzModel != NSHM23_SegmentationModels.CLASSIC;
				boolean allowedCalaveras = calaverasModel != NSHM23_SegmentationModels.CLASSIC;
				double fractPark = parkModel.getCreepingSectPassthrough();
				double fractSantaCruz = santaCruzModel.getCreepingSectPassthrough();
				double fractCalaveras = calaverasModel.getCreepingSectPassthrough();
				if (!allowedPark ||	(!allowedCalaveras && !allowedSantaCruz)) {
					// not actually possible
					allowed = false;
				} else {
					// possible;
					numAllowedEither++;
					if (allowedCalaveras)
						numAllowedCalaveras++;
					if (allowedSantaCruz)
						numAllowedSantaCruz++;
					if (allowedSantaCruz && allowedCalaveras)
						numAllowedBoth++;
					sumFractEither += Math.min(fractPark, fractSantaCruz + fractCalaveras);
					sumFractSantaCruz += Math.min(fractPark, fractSantaCruz);
					sumFractCalaveras += Math.min(fractPark, fractCalaveras);
				}
			}
		}
		DecimalFormat pDF = new DecimalFormat("0.00");
		System.out.println("Fraction allowed through either:\t"+pDF.format((double)numAllowedEither/(double)numTrials));
		System.out.println("Fraction allowed through both:\t"+pDF.format((double)numAllowedBoth/(double)numTrials));
		System.out.println("Fraction allowed Santa Cruz:\t"+pDF.format((double)numAllowedSantaCruz/(double)numTrials));
		System.out.println("Fraction allowed Calaveras:\t"+pDF.format((double)numAllowedCalaveras/(double)numTrials));
		System.out.println();
		System.out.println("Implied passthrough either:\t"+pDF.format(sumFractEither/(double)numTrials));
		System.out.println("Implied passthrough Santa Cruz:\t"+pDF.format(sumFractSantaCruz/(double)numTrials));
		System.out.println("Implied passthrough Calaveras:\t"+pDF.format(sumFractCalaveras/(double)numTrials));
		
	}
	
	private static void test267() throws IOException {
//		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.0.0");
//		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.0.0-noZone");
		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.0.0-ceusZoneOnly");
		IncludeBackgroundOption griddedOp = IncludeBackgroundOption.EXCLUDE;
		boolean subduction = false;
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.STABLE_SHALLOW);

		NshmErf erf = new NshmErf(erfPath, trts, griddedOp);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();
		
		int numSources = 0;
		int numRups = 0;
		
		int numPointSources = 0;
		int numPointRups = 0;
		
		for (ProbEqkSource src : erf) {
			boolean anyPoint = false;
			boolean allPoint = true;
			numSources++;
			for (ProbEqkRupture rup : src) {
				boolean point = rup.getRuptureSurface().getEvenlyDiscretizedNumLocs() == 1;
				anyPoint |= point;
				allPoint &= point;
				numRups++;
				if (point)
					numPointRups++;
			}
			if (anyPoint || allPoint) {
				System.out.println("Soruce "+src.getName()+" is point source! allPoint="+allPoint);
				numPointSources++;
			}
		}
		
		System.out.println("Summary for "+erfPath.getFileName());
		System.out.println(numSources+" sources and "+numRups+" rups");
		System.out.println(numPointSources+" point sources and "+numPointRups+" point rups");
	}
	
	private static void test268() throws IOException {
		File file = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.0.0/subduction/"
				+ "interface/Cascadia/features/Cascadia 2-1 (top).geojson");
		Feature feature = Feature.read(file);
		GeoJSONFaultSection sect = GeoJSONFaultSection.fromNSHMP_HazFeature(feature);
		String json = sect.toFeature().toJSON();
		System.out.println("Converted GeoJSONFaultSection feature:");
		System.out.println(json);
		// now read from JSON
		Feature feature2 = Feature.fromJSON(json);
		GeoJSONFaultSection sect2 = GeoJSONFaultSection.fromFeature(feature2);
		System.out.println("Converted GeoJSONFaultSection feature after reserialization:");
		System.out.println(sect2.toFeature().toJSON());
		System.out.println("Default surface type: "+sect2.getFaultSurface(1d).getClass().getName());
		
		List<FaultSection> plotSects = new ArrayList<>();
		FaultSection fullApproxSect = sect.clone();
		fullApproxSect.setSectionName("Approx Gridded Sect");
		fullApproxSect.setSectionId(plotSects.size());
		plotSects.add(fullApproxSect);
		
		Feature stirlingFeature = Feature.fromJSON(sect.toFeature().toJSON());
//		stirlingFeature.properties.remove(GeoJSONFaultSection.LOWER_TRACE); // TODO
		FaultSection stirlingSect = GeoJSONFaultSection.fromFeature(stirlingFeature);
		stirlingSect.setSectionName("Stirling Sect");
		stirlingSect.setSectionId(plotSects.size());
		plotSects.add(stirlingSect);
		
		FaultSection fakeAseisApproxSect = sect.clone();
		fakeAseisApproxSect.setSectionName("Fake Aseis Approx Gridded Sect");
		fakeAseisApproxSect.setSectionId(plotSects.size());
		fakeAseisApproxSect.setAseismicSlipFactor(0.25);
		plotSects.add(fakeAseisApproxSect);
		
		FaultSection fakeAseisStirlingSect = stirlingSect.clone();
		fakeAseisStirlingSect.setSectionName("Fake Aseis Stirling Sect");
		fakeAseisStirlingSect.setSectionId(plotSects.size());
		fakeAseisStirlingSect.setAseismicSlipFactor(0.25);
		plotSects.add(fakeAseisStirlingSect);
		
		Region plotReg = GeographicMapMaker.buildBufferedRegion(plotSects);
		GeographicMapMaker mapMaker = new GeographicMapMaker(plotReg);
		mapMaker.setPlotAseisReducedSurfaces(true);
		mapMaker.setWriteGeoJSON(true);
		
		for (int i=0; i<plotSects.size(); i++) {
			mapMaker.setFaultSections(List.of(plotSects.get(i)));
			mapMaker.plot(new File("/tmp"), "sub_surf_test_"+i, " ");
		}
		mapMaker.setFaultSections(plotSects);
		mapMaker.plot(new File("/tmp"), "sub_surf_tests", " ");
//		Preconditions.checkState(feature.geometry.type == GeoJSON_Type.MultiLineString);
//		MultiLineString geom = (MultiLineString)feature.geometry;
//		Preconditions.checkState(geom.lines.size() == 2);
//		LocationList upperTrace = geom.lines.get(0);
//		LocationList lowerTrace = geom.lines.get(1);
//		int id = ((Number)feature.id).intValue();
//		String name = feature.properties.get("name", null);
//		String state = feature.properties.get("state", null);
	}
	
	private static void test269() throws IOException {
		Region region = new Region(new Location(-33.5, 165), new Location(-47, 179));
		RSQSimGeographicMapMaker mapMaker = new RSQSimGeographicMapMaker(region, PoliticalBoundariesData.loadNZOutlines());
		
		RSQSimCatalog fullCatalog = Catalogs.BRUCE_5566.instance();
		RSQSimEvent event = fullCatalog.loader().byID(2983);
		
		mapMaker.setWritePDFs(false);
		mapMaker.setWriteGeoJSON(false);
//		mapMaker.plotEvent(event, Color.LIGHT_GRAY, Color.BLACK, 1f);
		List<Double> fakeScalars = new ArrayList<>();
		for (int i=0; i<event.getNumElements(); i++)
			fakeScalars.add(Math.random());
		CPT randCPT = new CPT(0d, 1d, Color.LIGHT_GRAY, Color.GRAY);
		mapMaker.plotEventFillScalars(event, fakeScalars, randCPT, null);
		mapMaker.plotEventHypocenter(Color.GREEN.darker());
		
		mapMaker.plot(new File("/tmp"), "event_map_debug", " ");
	}
	
	private static void test270() throws IOException {
		File file = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.0.0/stable-crust/zone/AR/"
				+ "Crowleys Ridge (south)/active/crowleys_ridge_south.geojson");
		Feature feature = Feature.read(file);
		// read in the 'mfd-tree' as a list of FeatureProperties instances, which are basically fancy maps
		// this will return null if an error occurs
		List<FeatureProperties> treePropsList = feature.properties.getPropertiesList("mfd-tree");
		System.out.println("Loaded a list of "+treePropsList.size()+" tree properties:");
		for (FeatureProperties treeProps : treePropsList) {
			String id = treeProps.getString("id");
			// need to supply a default value for primitives as null can't be returned if it's not found, thus the Double.NaN
			double weight = treeProps.getDouble("weight", Double.NaN);
			System.out.println("\tid="+id+" with weight="+weight);
			// 'value' is a map, so we'll load it in as its own FeatureProperties
			FeatureProperties value = treeProps.getProperties("value");
			String type = value.getString("type");
			double m = value.getDouble("m", Double.NaN);
			System.out.println("\t\ttype: "+type+", m="+m);
		}
//		FeatureProperties treeProps = feature.properties.getProperties("mfd-tree", null);
//		System.out.println(feature.toJSON());
	}
	
	private static void test271() throws IOException {
		System.out.println("Hello, World");
	}
	
	private static void test272() throws IOException {
		String cell = "5:6-61:57-489:493-5339:5341";
//		String cell = "5339:5341";
		System.out.println("Parsing "+cell);
		
		String[] dashSplit = cell.split("-");
		Preconditions.checkState(dashSplit.length >= 1);
		List<Integer> indexes = new ArrayList<>();
		
		for (String split : dashSplit) {
			String[] colonSplit = split.split(":");
			Preconditions.checkState(colonSplit.length == 2);
			int first = Integer.parseInt(colonSplit[0]);
			int last = Integer.parseInt(colonSplit[1]);
			Preconditions.checkState(first != last);
			if (first < last)
				for (int i=first; i<=last; i++)
					indexes.add(i);
			else
				for (int i=first; i>=last; i--)
					indexes.add(i);
		}
		System.out.println(indexes);
	}
	
	private static void test273() throws IOException {
		TableBuilder builder = MarkdownUtils.tableBuilder();
		
		builder.initNewLine();
		builder.addColumn("![Figure 22a](/assets/images/0120230122fig22a.png)");
		builder.addColumn("![Figure 22b](/assets/images/0120230122fig22b.png)");
		builder.finalizeLine();
		
		for (String line : builder.build())
			System.out.println(line);
	}
	
	private static void test274() throws IOException {

		Region relm = new CaliforniaRegions.RELM_TESTING();
		
		ComcatAccessor comcat = new ComcatAccessor();
		long startTime = new GregorianCalendar(2023, 0, 1).getTimeInMillis();
		long endTime = new GregorianCalendar(2024, 0, 1).getTimeInMillis();
		double minMag = 4.5d;
		ObsEqkRupList events = comcat.fetchEventList(null, startTime, endTime, -10, 50d,
				new ComcatRegionAdapter(relm), false, false, minMag, 1000, 100);
		
		int count = 0;
		for (EqkRupture rup : events) {
			if (relm.contains(rup.getHypocenterLocation())) {
				count++;
				System.out.println("M"+(float)rup.getMag()+" at "+rup.getHypocenterLocation());
			}
		}
		System.out.println("Found "+count+" events");
	}
	
	private static void test275() throws IOException {
		Map<IDPairing, Double> idsToProbs = new HashMap<>();
		Random r = new Random();
		for (int i=0; i<10; i++)
			idsToProbs.put(new IDPairing(r.nextInt(), r.nextInt()), r.nextDouble());
		HardcodedJumpProb prob = new HardcodedJumpProb("MyName", idsToProbs, false);
		
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		
		String json = gson.toJson(prob, HardcodedJumpProb.class);
		
		System.out.println(json);
		
		HardcodedJumpProb parsed = gson.fromJson(json, HardcodedJumpProb.class);
	}
	
	private static void test276() throws IOException {
		Map<String, Double> idsToProbs = new HashMap<>();
		Random r = new Random();
		for (int i=0; i<10; i++)
			idsToProbs.put(r.nextInt()+"", r.nextDouble());
		
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		
		String json = gson.toJson(idsToProbs, Map.class);
		
		System.out.println(json);
	}
	private static void test277() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-NSHM23_v3/true_mean_solution.zip"));
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.updateForecast();
		
		int numFSSRups = 0;
		int numGriddedRups = 0;
		int numFSSSources = erf.getNumFaultSystemSources();
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			int numRups = erf.getNumRuptures(sourceID);
			if (sourceID<numFSSSources)
				numFSSRups += numRups;
			else
				numGriddedRups += numRups;
		}
		int totRups = numFSSRups + numGriddedRups;
		System.out.println(totRups+" ruptures ("+numFSSRups+" fault, "+numGriddedRups+" gridded)");
		
		int numPerPeriod = 51;
		int numPeriods = 2;
		
		int exceedPointsPerRup = numPeriods * numPerPeriod;
		long bytesPerRup = exceedPointsPerRup * 8l;
		
		long totBytes = bytesPerRup*totRups;
		long totMB = totBytes / (1024l * 1024l);
		double totGB = (double)totMB / 1024d;
		
		long faultBytes = bytesPerRup*numFSSRups;
		long faultMB = faultBytes / (1024l * 1024l);
		double faultGB = (double)faultMB / 1024d;
		
		long gridBytes = bytesPerRup*numGriddedRups;
		long gridMB = gridBytes / (1024l * 1024l);
		double gridGB = (double)gridMB / 1024d;
		
		System.out.println("Total memory:\t"+totMB+" MB = "+(float)totGB+" GB");
		System.out.println("Fault memory:\t"+faultMB+" MB = "+(float)faultGB+" GB");
		System.out.println("Grid memory:\t"+gridMB+" MB = "+(float)gridGB+" GB");
	}
	
	private static void test278() throws IOException {
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
		System.out.println(gridReg.getNodeCount()+" nodes");
	}
	
	private static void test279() throws IOException {
		CPT cpt = new CPT();
		
		cpt.add(new CPTVal(-1.5f, Color.BLUE, -0.5f, Color.BLUE));
		cpt.add(new CPTVal(-0.5f, Color.WHITE, 0.5f, Color.WHITE));
		cpt.add(new CPTVal(0.5f, Color.RED, 1f, Color.RED));
		
		cpt.setNanColor(Color.GRAY);
		cpt.setBelowMinColor(Color.GRAY);
		cpt.setAboveMaxColor(Color.GRAY);
		
		cpt.writeCPTFile(new File("/tmp/for_luis.cpt"));
	}
	
	private static void test280() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/tmp/solution.zip"));
		sol.getModule(LogicTreeBranch.class);
	}
	
	private static void test281() throws IOException {
		FaultSystemRupSet refRupSet = FaultSystemRupSet.load(
				new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged.zip"));
		FaultSystemRupSet modRupSet = FaultSystemRupSet.load(new File("/tmp/rup_set_tests/mod_rup_set.zip"));
		
		// first make sure every mod rup set is actually unique
		HashMap<UniqueRupture, Integer> prevUnique = new HashMap<>();
		ClusterRuptures cRups = modRupSet.requireModule(ClusterRuptures.class);
		System.out.println("Testing mod rup set uniqueness");
		for (int r=0; r<modRupSet.getNumRuptures(); r++) {
			ClusterRupture cRup = cRups.get(r);
			Integer prevIndex = prevUnique.get(cRup.unique);
			if (prevIndex != null)
				throw new IllegalStateException("Duplicate rupture found: "+prevIndex+", "+r+
						"\n\t"+prevIndex+":\t"+modRupSet.getSectionsIndicesForRup(prevIndex)
						+"\n\t"+r+":\t"+modRupSet.getSectionsIndicesForRup(r));
		}
		System.out.println("Uniqueness verified!");
		
		System.out.println("Testing that mod rup set is identical to original");
		Preconditions.checkState(refRupSet.getNumRuptures() == modRupSet.getNumRuptures(),
				"Rup set mismatch: %s != %s", refRupSet.getNumRuptures(), modRupSet.getNumRuptures());
		
		for (int r=0; r<refRupSet.getNumRuptures(); r++) {
			List<Integer> rupSects1 = refRupSet.getSectionsIndicesForRup(r);
			List<Integer> rupSects2 = modRupSet.getSectionsIndicesForRup(r);
			
			Preconditions.checkState(rupSects1.size() == rupSects2.size(),
					"Rupture %s size mismatch: %s != %s", r, rupSects1.size(), rupSects2.size());
			Preconditions.checkState(rupSects1.equals(rupSects2),
					"Rupture %s data mismatch: %s != %s", r, rupSects1, rupSects2);
		}
		
		System.out.println("Verified!");
	}
	
	private static void test282() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip"));
		
		ScalarIMR gmm = AttenRelRef.ASK_2014.get();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		for (int i=0; i<xVals.size(); i++)
			xVals.set(i, 0d);
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		erf.getTimeSpan().setDuration(30d);
		erf.updateForecast();
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		Stopwatch watch = Stopwatch.createStarted();
		int num = 10;
		for (int i=0; i<num; i++) {
			Site site = new Site(new Location(34d + 0.1*Math.random(), -118d + 0.1*Math.random()));
			site.addParameterList(gmm.getSiteParams());
			calc.getHazardCurve(logXVals, site, gmm, erf);
			System.out.println("Done calc "+(i+1)+"/"+num);
		}
		watch.stop();
		
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		double secsPer = secs/(double)num;
		double mins = secs/60d;
		System.out.println("Took "+(float)secs+" s = "+(float)mins+" m ("+(float)secsPer+" s per curve)");
	}
	
	private static void test283() throws IOException {
		MeanUCERF2 erf = new MeanUCERF2();
		erf.updateForecast();
		MeanUCERF2 csERF = new MeanUCERF2();
		csERF.setParameter(MeanUCERF2.CYBERSHAKE_DDW_CORR_PARAM_NAME, true);
		csERF.updateForecast();
		
		DeformationModelSummaryFinal dmSummaries = new DeformationModelSummaryFinal();
		DeformationModelSummary dmSummary = (DeformationModelSummary) dmSummaries.getAllDeformationModels().get(0);
		Preconditions.checkNotNull(dmSummary, "DM summary is null?");
		A_FaultsFetcher aFaultsFetcher = new A_FaultsFetcher();
		aFaultsFetcher.setDeformationModel(dmSummary, false);
		ArrayList<FaultSegmentData> segments = aFaultsFetcher.getFaultSegmentDataList(true);
		ArrayList<FaultSection> sjSects = null;
		for (FaultSegmentData seg : segments) {
			if (seg.getFaultName().contains("Jacinto")) {
				System.out.println("Fault: "+seg.getFaultName()+"\ttotArea="+(float)(seg.getTotalArea()*1e-6));
				sjSects = seg.getPrefFaultSectionDataList();
				for (FaultSection sect : sjSects) {
					System.out.println("\t"+sect.getName());
					double len = sect.getTraceLength();
					double area = sect.getArea(true)*1e-6;
					double ddw = sect.getReducedDownDipWidth();
					System.out.println("\t\tlen="+(float)len+"km, ddw="+(float)ddw+", area="+(float)area);
				}
			}
		}
		
		List<int[]> srcRups = new ArrayList<>();
		List<FaultSection[]> rupSects = new ArrayList<>();
		
		srcRups.add(new int[] {98, 3});
		rupSects.add(new FaultSection[] {
				FaultSectionUtils.findSection(sjSects, "Anza)"),
				FaultSectionUtils.findSection(sjSects, "Anza", "stepover")
		});
		
		srcRups.add(new int[] {108, 2});
		rupSects.add(new FaultSection[] {
				FaultSectionUtils.findSection(sjSects, "Coyote Creek"),
				FaultSectionUtils.findSection(sjSects, "Borrego"),
				FaultSectionUtils.findSection(sjSects, "Superstition")
		});
		
		srcRups.add(new int[] {110, 2});
		rupSects.add(new FaultSection[] {
				FaultSectionUtils.findSection(sjSects, "Bernardino"),
				FaultSectionUtils.findSection(sjSects, "Valley)"),
				FaultSectionUtils.findSection(sjSects, "Valley", "stepover")
		});
		
		srcRups.add(new int[] {110, 0});
		rupSects.add(new FaultSection[] {
				FaultSectionUtils.findSection(sjSects, "Bernardino"),
				FaultSectionUtils.findSection(sjSects, "Valley)"),
				FaultSectionUtils.findSection(sjSects, "Valley", "stepover")
		});
		
		ScalingRelationships aveU2 = ScalingRelationships.AVE_UCERF2;
		
		Somerville_2006_MagAreaRel ma = new Somerville_2006_MagAreaRel();
		
		for (int i=0; i<srcRups.size(); i++) {
			int[] srcRup = srcRups.get(i);
			ProbEqkSource source = erf.getSource(srcRup[0]);
			ProbEqkRupture rup = source.getRupture(srcRup[1]);
			double len = rup.getRuptureSurface().getAveLength();
			double mag = rup.getMag();
			double origDDW = rup.getRuptureSurface().getAveWidth();
			double origRupArea = len*origDDW;
			
			FaultSection[] sects = rupSects.get(i);
			double origTotArea = 0d;
			double origTotLen = 0d;
			for (FaultSection sect : sects) {
				origTotArea += sect.getArea(true)*1e-6;
				origTotLen += sect.getTraceLength();
			}
			double origTotDDW = origTotArea/origTotLen;
			
			double calcMag = aveU2.getMag(origTotArea*1e6, origTotLen*1e3, origTotDDW*1e3, origTotDDW*1e3, rup.getAveRake());
			
//			double calcArea = ma.getMedianArea(calcMag);
			double calcArea = ma.getMedianArea(rup.getMag());
			double ddwFactor = calcArea/origTotArea;
			double calcDDW = origDDW*ddwFactor;
			
			double csDDW = csERF.getRupture(srcRup[0], srcRup[1]).getRuptureSurface().getAveWidth();
			double estMagForActualDDW = ma.getMedianMag(csDDW*origTotLen);
			System.out.println("Source "+srcRup[0]+" ("+source.getName()+"), Rupture "+srcRup[1]+", M"+(float)mag);
//			System.out.println(source.getClass());
			System.out.println("\tOriginal Rupture dimensions: "+(float)len+" km long x "+(float)origDDW+" km wide = "+(float)origRupArea+" km^2");
			System.out.println("\tOriginal section dimensions: "+(float)origTotLen+" km long x "+(float)origTotDDW+" km wide = "+(float)origTotArea+" km^2");
			System.out.println("\tOriginal median magnitude: "+(float)calcMag);
			System.out.println("\tSomerville (2006) area for M"+(float)calcMag+": "+(float)calcArea+" km^2");
			System.out.println("\tCalculated Extended DDW: "+(float)calcDDW+" km");
			System.out.println("\tCyberShake ERF DDW: "+(float)csDDW+" km (corresponds to M"+(float)estMagForActualDDW+" w/ Somerville)");
			System.out.println();
		}
	}
	
	private static void test284() throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		// disable caps
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV_UPPER_BOUND = 0d;
		NSHM23_DeformationModels.HARDCODED_FRACTIONAL_STD_DEV = 0d;
		int totNumSubSects = 0;
		int totNumBelow0p1 = 0;
		int totNumAbove1 = 0;
		
		int avgNumSubSects = 0;
		int avgNumBelow0p1 = 0;
		int avgNumAbove1 = 0;
		for (NSHM23_DeformationModels dm : NSHM23_DeformationModels.values()) {
			if (dm != NSHM23_DeformationModels.AVERAGE && dm.getNodeWeight(null) == 0d)
				continue;
			List<? extends FaultSection> subSects = dm.build(fm);
			
			int myBelow0p1 = 0;
			int myAbove1 = 0;
			for (FaultSection sect : subSects) {
				double slip = sect.getOrigAveSlipRate();
				double sd = sect.getOrigSlipRateStdDev();
				double cov = sd/slip;
				if (cov < 0.1)
					myBelow0p1++;
				if (cov > 1d || slip == 0d)
					myAbove1++;
			}
			
			if (dm == NSHM23_DeformationModels.AVERAGE) {
				avgNumAbove1 = myAbove1;
				avgNumBelow0p1 = myBelow0p1;
				avgNumSubSects = subSects.size();
			} else {
				totNumSubSects += subSects.size();
				totNumBelow0p1 += myBelow0p1;
				totNumAbove1 += myAbove1;
			}
		}
		
		DecimalFormat pDF = new DecimalFormat("0.##%");
		System.out.println(pDF.format((double)totNumBelow0p1/(double)totNumSubSects)+" in total have COV<0.1");
		System.out.println(pDF.format((double)totNumAbove1/(double)totNumSubSects)+" in total have COV>1");
		
		// now average model
		System.out.println(pDF.format((double)avgNumBelow0p1/(double)avgNumSubSects)+" on average have COV<0.1");
		System.out.println(pDF.format((double)avgNumAbove1/(double)avgNumSubSects)+" on average have COV>1");
	}
	
	private static void test285() throws IOException {
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		NSHM23_DeformationModels dm = NSHM23_DeformationModels.GEOLOGIC;
		List<? extends FaultSection> allSects = fm.getFaultSections();
		List<FaultSection> stateSects = new ArrayList<>();
		NSHM23_SingleStates state = NSHM23_SingleStates.UT;
		for (FaultSection sect : allSects)
			if (state.contains((GeoJSONFaultSection)sect))
				stateSects.add(sect);
		System.out.println("Retained "+stateSects.size()+"/"+allSects.size()+" sects");
		List<? extends FaultSection> subSects = dm.buildForSubsects(fm, SubSectionBuilder.buildSubSects(stateSects));
		GeoJSONFaultReader.writeFaultSections(new File("/tmp/ut_sub_sects.geojson"), subSects);
	}
	
	private static void test286() throws IOException {
		double[] dists = { 0d, 1d, 10d, 100d, 200d, 300 };
		double maxDist = 200d;
		double zeroDistCoeff = 0.95;
		for (double dist : dists) {
			System.out.println("Dist: "+(float)dist);
			System.out.println("\t"+(float)Interpolate.findY(0d, zeroDistCoeff, maxDist, 0d, dist));
			System.out.println("\t"+(float)(zeroDistCoeff - dist*0.95/200));
		}
	}
	
	private static void test287() throws IOException {
		Map<Integer, FaultSection> fmSectsMap = NSHM23_FaultModels.WUS_FM_v3.getFaultSectionIDMap();
		Map<Integer, GeoJSONFaultSection> geoDMSects = NSHM23_DeformationModels.getGeolFullSects(NSHM23_FaultModels.WUS_FM_v3);
		
		int numTests = 0;
		for (int parentID : geoDMSects.keySet()) {
			GeoJSONFaultSection dmSect = geoDMSects.get(parentID);
			FaultSection fmSect = fmSectsMap.get(parentID);
			Preconditions.checkState((float)dmSect.getAveRake() == (float)fmSect.getAveRake());
			numTests++;
		}
		System.out.println("Validated rakes of "+numTests+" sections");
	}
	
	private static void test288() throws IOException {
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/tmp/subduction_rup_set_FULL.zip"));
		
		List<? extends FaultSection> subSects = rupSet.getFaultSectionDataList();
		
//		int startIndex = 20;
		int startIndex = 51;
		
		AspectRatioFilter filter = new AspectRatioFilter(0.75f, false, false, null);
		
		for (int num=1; num<=5; num++) {
			List<? extends FaultSection> sects = subSects.subList(startIndex, startIndex+num);
			Preconditions.checkState(sects.size() == num);
			double aspect = AspectRatioFilter.clusterAspectRatio(sects);
			
			ClusterRupture rupture = new ClusterRupture(new FaultSubsectionCluster(sects));
			System.out.println("=============================================");
			System.out.println("Calculating for "+sects.size()+" sects: "+rupture);
			PlausibilityResult result = filter.apply(rupture, true);
			System.out.println("\tResult: "+result);
			System.out.println("\tAspect: "+aspect);
			System.out.println("=============================================");
		}
	}
	
	private static void test289() throws IOException {
//		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/tmp/subduction_rup_set_FULL.zip"));
		FaultSystemRupSet rupSet = FaultSystemRupSet.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_05_15-prvi25_subduction_branches-PRVI_SUB_FM_INITIAL/results_PRVI_SUB_FM_INITIAL_branch_averaged.zip"));
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			sect.getFaultSurface(1d, false, false);
			sect.getFaultSurface(1d, false, true);
		}
	}
	
	private static void test290() throws IOException {
		File dir = new File("/project/scec_608/kmilner/nshm23/batch_inversions/2024_05_15-prvi25_crustal_branches");
		SolutionLogicTree slt = SolutionLogicTree.load(new File(dir, "results.zip"));
		List<LogicTreeBranch<?>> connectedBranches = new ArrayList<>();
		int parentID = 39;
		for (LogicTreeBranch<?> branch : slt.getLogicTree()) {
			FaultSystemSolution sol = slt.forBranch(branch);
			
			boolean connected = false;
			for (int rupIndex : sol.getRupSet().getRupturesForParentSection(parentID)) {
				if (sol.getRateForRup(rupIndex) > 0d) {
					for (FaultSection sect : sol.getRupSet().getFaultSectionDataForRupture(rupIndex)) {
						if (sect.getParentSectionId() != parentID) {
							connected = true;
							break;
						}
					}
				}
			}
			if (connected)
				connectedBranches.add(branch);
		}
		
		System.out.println(connectedBranches.size()+"/"+slt.getLogicTree().size()+" branches have connected proxy faults:");
		LogicTreeBranch<?> common = null;
		for (LogicTreeBranch<?> branch : connectedBranches) {
			if (common == null)
				common = branch.copy();
			for (int i=0; i<branch.size(); i++)
				if (!common.hasValue(branch.getValue(i)))
					common.clearValue(i);
			System.out.println(branch);
		}
		System.out.println("Things in common: "+common);
	}
	
	private static void test291() throws IOException {
		AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF = 2.55;
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/data/kevin/git/ucerf3-etas-launcher/inputs/"
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/UCERF3/rup_sets/orig/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip"));
		
		FaultSystemSolutionERF_ETAS erf = ETAS_Launcher.buildERF(sol, false, 500d, 2012);
		
		erf.updateForecast();
		
		System.out.println("RupSet is a "+sol.getRupSet().getClass().getName());
		System.out.println("Sol is a "+sol.getClass().getName());
		System.out.println("RupSet has ModSectMinMags? "+sol.getRupSet().hasModule(ModSectMinMags.class));
		System.out.println("Sol has RupMFDs? "+sol.hasModule(RupMFDsModule.class));
		System.out.println("ERF totNumRups: "+erf.getTotNumRups());
		System.out.println("ERF getTotNumRupsFromFaultSystem: "+erf.getTotNumRupsFromFaultSystem());
		System.out.println("ERF getNumFaultSystemSources: "+erf.getNumFaultSystemSources());
		System.out.println("ERF getNumSources: "+erf.getNumSources());
	}
	
	private static void test292() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_05_21-prvi25_crustal_branches-GEOLOGIC/results_PRVI_FM_INITIAL_branch_averaged.zip"));
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
		erf.updateForecast();
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		GeographicMapMaker mapMaker = new GeographicMapMaker(rupSet.getFaultSectionDataList());
		
		List<LocationList> lines = new ArrayList<>();
		List<Double> lineMags = new ArrayList<>();
		double minMag = Double.POSITIVE_INFINITY;
		double maxMag = Double.NEGATIVE_INFINITY;
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			boolean allProxies = true;
			for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex))
				allProxies &= sect.isProxyFault();
			int sourceID = erf.getSrcIndexForFltSysRup(rupIndex);
			if (sourceID > 0 && allProxies) {
				double mag = rupSet.getMagForRup(rupIndex);
				for (ProbEqkRupture rup : erf.getSource(sourceID)) {
					RuptureSurface surf = rup.getRuptureSurface();
					lines.add(surf.getUpperEdge());
					lineMags.add(mag);
					minMag = Double.min(minMag, mag);
					maxMag = Double.max(maxMag, mag);
				}
			}
		}
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(minMag, maxMag);
		List<Color> colors = new ArrayList<>();
		for (double mag : lineMags)
			colors.add(cpt.getColor((float)mag));
		mapMaker.plotLines(lines, colors, 2f);
		
		mapMaker.plot(new File("/tmp"), "proxy_finite_test", " ");
	}
	
	private static void test293() throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/2024_05_21-prvi25_crustal_branches/node_branch_averaged");
		List<NSHM23_ScalingRelationships> scales = new ArrayList<>();
		Table<NSHM23_ScalingRelationships, Integer, Double> mmins = HashBasedTable.create();
		Table<NSHM23_ScalingRelationships, Integer, Double> mminRIs = HashBasedTable.create();
		Table<NSHM23_ScalingRelationships, Integer, Double> mmaxs = HashBasedTable.create();
		Table<NSHM23_ScalingRelationships, Integer, Double> mmaxRIs = HashBasedTable.create();
		
		Map<Integer, String> proxyNames = new LinkedHashMap<>();

		DecimalFormat magDF = new DecimalFormat("0.00");
		DecimalFormat riDF = new DecimalFormat("0.0");
		for (NSHM23_ScalingRelationships scale : NSHM23_ScalingRelationships.values()) {
			File solFile = new File(dir, "Scale_"+scale.getFilePrefix()+".zip");
			if (solFile.exists()) {
				scales.add(scale);
				FaultSystemSolution sol = FaultSystemSolution.load(solFile);
				FaultSystemRupSet rupSet = sol.getRupSet();
				
				HashSet<Integer> proxyParentIDs = new HashSet<>();
				for (FaultSection sect : rupSet.getFaultSectionDataList()) {
					if (sect.isProxyFault()) {
						proxyParentIDs.add(sect.getParentSectionId());
						if (!proxyNames.containsKey(sect.getParentSectionId()))
							proxyNames.put(sect.getParentSectionId(), sect.getParentSectionName());
					}
				}
				
				System.out.println("Processing "+scale.getName());
				for (int parentID : proxyParentIDs) {
					double mMin = Double.POSITIVE_INFINITY;
					double mMinRate = 0d;
					double mMax = 0d;
					double mMaxRate = 0d;
					for (int rupIndex : rupSet.getRupturesForParentSection(parentID)) {
						boolean multiParents = false;
						for (FaultSection sect : rupSet.getFaultSectionDataForRupture(rupIndex)) {
							if (sect.getParentSectionId() != parentID) {
								multiParents = true;
								break;
							}
						}
						if (multiParents) {
							// this ruptures outside of the zone and is excluded in the inversion
							Preconditions.checkState(sol.getRateForRup(rupIndex) == 0d);
							continue;
						}
						double mag = rupSet.getMagForRup(rupIndex);
						mMin = Math.min(mMin, mag);
						mMinRate += sol.getRateForRup(rupIndex);
						if (mag > mMax) {
							mMax = mag;
							mMaxRate = sol.getRateForRup(rupIndex);
						} else {
							Preconditions.checkState(mag < mMax, "more than 1 rup at Mmax?");
						}
					}
					double mMinRI = 1d/mMinRate;
					double mMaxRI = 1d/mMaxRate;
					System.out.println("\tMag Range: ["+magDF.format(mMin)+", "+magDF.format(mMax));
					System.out.println("\tM>="+magDF.format(mMin)+" RI: "+riDF.format(mMinRI)+", rate="+(float)mMinRate);
					System.out.println("\tM="+magDF.format(mMax)+" RI: "+riDF.format(mMaxRI)+", rate="+(float)mMaxRate);
					mmins.put(scale, parentID, mMin);
					mminRIs.put(scale, parentID, mMinRI);
					mmaxs.put(scale, parentID, mMax);
					mmaxRIs.put(scale, parentID, mMaxRI);
				}
			}
		}
		
		CSVFile<String> csv = new CSVFile<>(false);
		for (int parentID : proxyNames.keySet()) {
			String name = proxyNames.get(parentID);
			System.out.println();
			System.out.println(name);
			System.out.println("Scale\tMmin\tRI >= Mmin\tMmax\tRI @ Mmax");
			csv.addLine(name);
			csv.addLine("Scaling", "Mmin", "RI>=Mmin", "Mmax", "RI@Mmax");
			for (NSHM23_ScalingRelationships scale : scales) {
				double mMin = mmins.get(scale, parentID);
				double mMinRI = mminRIs.get(scale, parentID);
				double mMax = mmaxs.get(scale, parentID);
				double mMaxRI = mmaxRIs.get(scale, parentID);
				System.out.println(scale.getShortName()+"\t"+magDF.format(mMin)+"\t"+riDF.format(mMinRI)
					+"\t"+magDF.format(mMax)+"\t"+riDF.format(mMaxRI));
				csv.addLine(scale.getShortName(), magDF.format(mMin), riDF.format(mMinRI), magDF.format(mMax), riDF.format(mMaxRI));
			}
		}
		
		csv.writeToFile(new File("/tmp/prvi_proxy_mags.csv"));
	}

	private static void test294() throws IOException {
		File file = new File("/home/kevin/workspace/opensha/src/main/resources/data/erf/prvi25/fault_models/"
//				+ "subduction/inputs/PRVI_subduction_FullRate_LargePolys_drape_removeFields.geojson");
				+ "subduction/inputs/PRVI_subduction_PartRate_LargePolys_drape_removeFields.geojson");
		Gson inGSON = FeatureCollection.buildGson(DepthSerializationType.DEPTH_M);
		BufferedReader read = new BufferedReader(new FileReader(file));
		FeatureCollection features = inGSON.fromJson(read, FeatureCollection.class);
		FeatureCollection.write(features, new File("/tmp/"+file.getName()));
	}

	private static void test295() throws IOException {
		// nth rup from catalog and ETAS_Launcher.buildERF
		File etasDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2024_06_06-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog/");
		File configFile = new File(etasDir, "config.json");
		File binFile = new File(etasDir, "results_m5_preserve_chain.bin");
		
		ETAS_Config config = ETAS_Config.readJSON(configFile);
		
		ETAS_Launcher launcher = new ETAS_Launcher(config, false);
		AbstractNthRupERF erf = launcher.checkOutERF();
		
		// ETAS doesn't seem to track rakes
		boolean checkRakes = false;
		boolean primaryOnly = false;
		
		int catCount = 0;
		int rupCount = 0;
		for (ETAS_Catalog catalog : ETAS_CatalogIO.getBinaryCatalogsIterable(binFile, 0)) {
			System.out.println("Validating catalog "+catCount);
			for (ETAS_EqkRupture rup : catalog) {
				int n = rup.getNthERF_Index();
				double mag = rup.getMag();
				double rake = rup.getAveRake();
				if (primaryOnly && rup.getGeneration() > 0)
					continue;
//				System.out.println("Testing rup with M="+(float)rup.getMag()+", n="+n+", fss="+rup.getFSSIndex()
//						+", grid="+rup.getGridNodeIndex()+", gen="+rup.getGeneration());
				ProbEqkRupture erfRup = erf.getNthRupture(n);
				Preconditions.checkState((float)mag == (float)erfRup.getMag(),
						"Mag mismatch for n=%s (fss=%s, grid=%s, cat=%s, rup=%s, gen=%s): etas=%s, erf=%s",
						n, rup.getFSSIndex(), rup.getGridNodeIndex(), catCount, rupCount, rup.getGeneration(), mag, erfRup.getMag());
				Preconditions.checkState(!checkRakes || (float)rake == (float)erfRup.getAveRake(),
						"Rake mismatch for n=%s (fss=%s, grid=%s, cat=%s, rup=%s, gen=%s): etas=%s, erf=%s",
						n, rup.getFSSIndex(), rup.getGridNodeIndex(), catCount, rupCount, rup.getGeneration(), rake, erfRup.getAveRake());
				rupCount++;
			}
			catCount++;
		}
		System.out.println("DONE, validated "+rupCount+" ruptures in "+catCount+" catalogs");
	}
	
	private static void test296() throws IOException {
		// nth rup from catalog and ETAS_Launcher.buildERF
		File etasDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2024_05_24-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog/");
		File configFile = new File(etasDir, "config.json");
		File binFile = new File(etasDir, "results_m5_preserve_chain.bin");
		
		ETAS_Config config = ETAS_Config.readJSON(configFile);
		
		ETAS_Launcher launcher = new ETAS_Launcher(config, false);
		AbstractNthRupERF erf = launcher.checkOutERF();
		
		erf.testNthRupIndicesForSource();
		
		int index = 0;
		for (int sourceID=0; sourceID<erf.getNumSources(); sourceID++) {
			ProbEqkSource source = erf.getSource(sourceID);
			for (int rupID=0; rupID<source.getNumRuptures(); rupID++) {
				ProbEqkRupture rup = source.getRupture(rupID);
				int nthRup = erf.getIndexN_ForSrcAndRupIndices(sourceID, rupID);
				Preconditions.checkState(nthRup == index,
						"expected n=%s but erf.getIndexN_ForSrcAndRupIndices(%s, %s)=%s",
						index, sourceID, rupID, nthRup);
				index++;
				ProbEqkRupture rup2 = erf.getNthRupture(nthRup);
				Preconditions.checkState((float)rup.getMag() == (float)rup2.getMag(),
						"Mag mismatch for erf.getRupture(%s, %s) vs erf.getRupture(%s=erf.getIndexN_ForSrcAndRupIndices(%s, %s)): %s != %s",
						sourceID, rupID, nthRup, sourceID, rupID, (float)rup.getMag(), (float)rup2.getMag());
			}
		}
		System.out.println("All nth tests passed");
	}
	
	private static void test297() throws IOException {
		// nth rup from catalog and ETAS_Launcher.buildERF
		File etasDir = new File("/home/kevin/OpenSHA/UCERF3/etas/simulations/"
				+ "2024_05_24-Start2012_500yr_kCOV1p5_Spontaneous_HistCatalog/");
		File configFile = new File(etasDir, "config.json");
		File binFile = new File(etasDir, "results_m5_preserve_chain.bin");
		
		ETAS_Config config = ETAS_Config.readJSON(configFile);
		
		ETAS_Launcher launcher = new ETAS_Launcher(config, false);
		FaultSystemSolutionERF erf = (FaultSystemSolutionERF) launcher.checkOutERF();
		MFDGridSourceProvider gridProv = erf.getSolution().requireModule(MFDGridSourceProvider.class);
		
		int randSrcIndex = 252798;
		int gridRegionIndex = randSrcIndex - erf.getNumFaultSystemSources();
		int isSubSeismo = 1;
		int filteredR = 2;
		Location translatedParLoc = new Location(35, -118);
		ProbEqkSource src = erf.getSource(randSrcIndex);
		ProbEqkSource filteredSrc;
		if(isSubSeismo == 1)
			filteredSrc = gridProv.getSourceSubSeisOnFault(gridRegionIndex,
					erf.getTimeSpan().getDuration(), null, BackgroundRupType.POINT);
		else
			filteredSrc = gridProv.getSourceUnassociated(gridRegionIndex,
					erf.getTimeSpan().getDuration(), null, BackgroundRupType.POINT);
		
		ProbEqkRupture filteredRup = filteredSrc.getRupture(filteredR);
		System.out.println("Filtered r="+filteredR+", M"+filteredRup.getMag()+", rake="+filteredRup.getAveRake());
		
		int r = -1;
		for (int r2=0; r2<src.getNumRuptures(); r2++) {
			ProbEqkRupture rup = src.getRupture(r2);
			System.out.println("Candidate r="+r2+", M"+rup.getMag()+", rake="+rup.getAveRake());
			if (Precision.equals(rup.getAveRake(), filteredRup.getAveRake()) && Precision.equals(rup.getMag(), filteredRup.getMag())) {
				boolean match = true;
				if (rup.getRuptureSurface().getAveDip() < 90d) {
					// potential match, check footwall
					if (rup.getRuptureSurface() instanceof PointSurfaceNshm) {
						Preconditions.checkState(filteredRup.getRuptureSurface() instanceof PointSurfaceNshm);
						PointSurfaceNshm surf1 = (PointSurfaceNshm) rup.getRuptureSurface();
						PointSurfaceNshm surf2 = (PointSurfaceNshm) filteredRup.getRuptureSurface();
						match = surf1.isOnFootwall() == surf2.isOnFootwall();
					} else if (rup.getRuptureSurface() instanceof PointSurface13b) {
						Preconditions.checkState(filteredRup.getRuptureSurface() instanceof PointSurface13b);
						PointSurface13b surf1 = (PointSurface13b) rup.getRuptureSurface();
						PointSurface13b surf2 = (PointSurface13b) filteredRup.getRuptureSurface();
						match = surf1.isOnFootwall() == surf2.isOnFootwall();
					} else {
						double rx1 = rup.getRuptureSurface().getDistanceX(translatedParLoc); // any loc will do here
						double rx2 = rup.getRuptureSurface().getDistanceX(translatedParLoc); // any loc will do here
						match = (rx1 > 0) == (rx2 > 0);
					}
				}
				if (match) {
					Preconditions.checkState(r < 0, "Multiple rupture mappings? sourceID=%s, subSeismo=%s,"
							+ "filteredR=%s, mag=%s, rake=%s, sourceType=%s, surfType=%s",
							randSrcIndex, isSubSeismo, filteredR, filteredRup.getMag(), filteredRup.getAveRake(),
							filteredSrc.getClass().getName(), filteredRup.getRuptureSurface().getClass().getName());
					r = r2;
				}
			}
		}
		Preconditions.checkState(r >= 0, "No rupture mapping found sourceID=%s, subSeismo=%s,"
				+ "filteredR=%s, mag=%s, rake=%s, sourceType=%s, surfType=%s",
				randSrcIndex, isSubSeismo, filteredR, filteredRup.getMag(), filteredRup.getAveRake(),
				filteredSrc.getClass().getName(), filteredRup.getRuptureSurface().getClass().getName());
		ProbEqkRupture rup = src.getRupture(r);
		System.out.println("Matched with r="+r+", M"+rup.getMag()+", rake="+rup.getAveRake());
		
	}
	
	private static void test298() throws IOException {
		File inDir = new File("/data/kevin/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/node_branch_averaged");
		
		BranchAverageSolutionCreator baCreatorWithout = new BranchAverageSolutionCreator(new BranchWeightProvider.CurrentWeights());
		BranchAverageSolutionCreator baCreatorWith = new BranchAverageSolutionCreator(new BranchWeightProvider.CurrentWeights());
		
		LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(List.of(NSHM23_LogicTreeBranch.SEG));
		
		for (NSHM23_SegmentationModels segModel : NSHM23_SegmentationModels.values()) {
			File solFile = new File(inDir, "SegModel_"+segModel.getFilePrefix()+".zip");
			if (!solFile.exists())
				continue;
			FaultSystemSolution sol = FaultSystemSolution.load(solFile);
			sol.removeModuleInstances(RupMFDsModule.class);
			branch = branch.copy();
			branch.setValue(segModel);
			if (segModel != NSHM23_SegmentationModels.CLASSIC)
				baCreatorWithout.addSolution(sol, branch);
			baCreatorWith.addSolution(sol, branch);
		}
		
		baCreatorWithout.build().write(new File("/tmp/ba_without_classic.zip"));
		baCreatorWith.build().write(new File("/tmp/ba_with_classic.zip"));
	}
	
	private static void test299() throws IOException {
		double mean = 1e5;
		LossCOV_Model covModel = LossCOV_Model.PORTER_POWER_LAW_2020_09_01;
		double cov = covModel.getCOV(mean);
		double samples[] =      covModel.getDistribution(mean).sample(1000000);
		DescriptiveStatistics stats = new DescriptiveStatistics(samples);
		double meanFromDist = stats.getMean();
		double covFromDist = stats.getStandardDeviation()/meanFromDist;

		System.out.println("mean="+mean+";\tmeanFromDist="+meanFromDist+";\tcov="+cov+";\tcovFromDist="+covFromDist+
				"; meanRatio="+(float)(meanFromDist/mean)+"; covRatio="+(float)(covFromDist/cov));
		
		double sigma = Math.sqrt(Math.log(cov*cov+1));
		System.out.println("Sigma="+sigma);
//		double mu = Math.log(mean)-(sigma*sigma/2);
		double mu = Math.log(mean);
		LogNormalDistribution dist = new LogNormalDistribution(mu, sigma);
		double samples2[] = dist.sample(1000000);
		DescriptiveStatistics stats2 = new DescriptiveStatistics(samples2);
		double meanFromDist2 = stats2.getMean();
		double covFromDist2 = stats2.getStandardDeviation()/meanFromDist2;

		System.out.println("mean="+mean+";\tmeanFromDist2="+meanFromDist2+";\tcov="+cov+";\tcovFromDist2="+covFromDist2+
				"; meanRatio2="+(float)(meanFromDist2/mean)+"; covRatio2="+(float)(covFromDist2/cov));
		
		System.out.println("Now calculating mean/COV in log space instead...");
		for (int i=0; i<samples.length; i++)
			samples[i] = Math.log(samples2[i]);
		stats = new DescriptiveStatistics(samples);
		meanFromDist = stats.getMean();
		covFromDist = stats.getStandardDeviation()/meanFromDist;
		meanFromDist = Math.exp(meanFromDist);
//		covFromDist = Math.exp(stats.getStandardDeviation())/meanFromDist;
		System.out.println("mean="+mean+";\tmeanFromDist="+meanFromDist+";\tcov="+cov+";\tcovFromDist="+covFromDist+
				"; meanRatio="+(float)(meanFromDist/mean)+"; covRatio="+(float)(covFromDist/cov));
	}
	
	private static void test300() throws IOException {
		NormalDistribution normDist = new NormalDistribution();
		for (int i=0; i<5; i++) {
			System.out.println("Sigma="+i);
			System.out.println("\tCDF1: "+(float)GaussianDistCalc.getCDF(i));
			System.out.println("\tCDF2: "+(float)normDist.cumulativeProbability(i));
		}
		for (int i=1; i<6; i++) {
			System.out.println("Sigmas="+i);
			DiscretizedFunc func = InvCDF_RJBCorrPointSurface.getSigmaSpacedProbs(i);
			for (Point2D pt : func)
				System.out.println("\t"+(float)pt.getX()+": "+(float)pt.getY());
		}
	}
	
	private static void test301() throws IOException {
		Path erfPath = Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-6.0.0");
		IncludeBackgroundOption griddedOp = IncludeBackgroundOption.ONLY;
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW);

		NshmErf erf = new NshmErf(erfPath, trts, griddedOp);
		erf.getTimeSpan().setDuration(1.0);
		erf.updateForecast();

		Table<Float, Float, Integer> magRakeCounts = HashBasedTable.create();
		Table<Float, Float, Integer> magRakeCountsRXPos = HashBasedTable.create();
		Table<Float, Float, Integer> magRakeCountsRXNeg = HashBasedTable.create();
		Table<Float, Float, Integer> magRakeBugCounts = HashBasedTable.create();
		Table<Float, Float, Integer> magRakeBugCountsRXPos = HashBasedTable.create();
		Table<Float, Float, Integer> magRakeBugCountsRXNeg = HashBasedTable.create();
		
		Location rxTestLoc = new Location(0d, 0d);
		
		for (ProbEqkSource source : erf) {
			for (ProbEqkRupture rup : source) {
				RuptureSurface surf = rup.getRuptureSurface();
				Preconditions.checkState(surf instanceof NshmSurface);
				NshmSurface nshmSurf = (NshmSurface)surf;
				LocationList locs = nshmSurf.getEvenlyDiscritizedListOfLocsOnSurface();
				Preconditions.checkState(locs.size() == 1);
				Location surfLoc = locs.get(0);
				Location siteLoc = new Location(surfLoc.lat, surfLoc.lon);
				double zTOR = surf.getAveRupTopDepth();
				double rRup = surf.getDistanceRup(siteLoc);
				float rake = (float)rup.getAveRake();
				float mag = (float)rup.getMag();
				int prevCount = magRakeCounts.contains((float)mag, (float)rake) ? magRakeCounts.get(mag, rake) : 0;
				magRakeCounts.put(mag, rake, prevCount+1);
				boolean rxNeg = surf.getDistanceX(rxTestLoc) < 0d;
				if (rxNeg) {
					prevCount = magRakeCountsRXNeg.contains((float)mag, (float)rake) ? magRakeCountsRXNeg.get(mag, rake) : 0;
					magRakeCountsRXNeg.put(mag, rake, prevCount+1);
				} else {
					prevCount = magRakeCountsRXPos.contains((float)mag, (float)rake) ? magRakeCountsRXPos.get(mag, rake) : 0;
					magRakeCountsRXPos.put(mag, rake, prevCount+1);
				}
				if ((float)rRup < (float)zTOR) {
//					System.out.println("rRup bug! zTOR="+(float)zTOR+", rRup="+(float)rRup+", M="+mag+", rake="+rake);
					prevCount = magRakeBugCounts.contains((float)mag, (float)rake) ? magRakeBugCounts.get(mag, rake) : 0;
					magRakeBugCounts.put(mag, rake, prevCount+1);
					if (rxNeg) {
						prevCount = magRakeBugCountsRXNeg.contains((float)mag, (float)rake) ? magRakeBugCountsRXNeg.get(mag, rake) : 0;
						magRakeBugCountsRXNeg.put(mag, rake, prevCount+1);
					} else {
						prevCount = magRakeBugCountsRXPos.contains((float)mag, (float)rake) ? magRakeBugCountsRXPos.get(mag, rake) : 0;
						magRakeBugCountsRXPos.put(mag, rake, prevCount+1);
					}
				}
			}
		}
		List<Float> mags = new ArrayList<>(magRakeCounts.rowKeySet());
		Collections.sort(mags);
		DecimalFormat pDF = new DecimalFormat("0.00%");
		for (float mag : mags) {
			Map<Float, Integer> countMap = magRakeCounts.row(mag);
			List<Float> rakes = new ArrayList<>(countMap.keySet());
			Collections.sort(rakes);
			System.out.println("M"+mag);
			for (float rake : rakes) {
				int count = countMap.get(rake);
				int bugCount = magRakeBugCounts.contains(mag, rake) ? magRakeBugCounts.get(mag, rake) : 0;
				System.out.println("\trake="+rake+":\t"+pDF.format((double)bugCount/(double)count)+" buggy rRup");
				if (rake != 0f) {
					int countRXNeg = magRakeCountsRXNeg.get(mag, rake);
					int countRXPos = magRakeCountsRXPos.get(mag, rake);
					int bugCountRXNeg = magRakeBugCountsRXNeg.contains(mag, rake) ? magRakeBugCountsRXNeg.get(mag, rake) : 0;
					int bugCountRXPos = magRakeBugCountsRXPos.contains(mag, rake) ? magRakeBugCountsRXPos.get(mag, rake) : 0;
					System.out.println("\t\t"+pDF.format((double)bugCountRXPos/(double)countRXNeg)+" buggy rRup when rX>=0");
					System.out.println("\t\t"+pDF.format((double)bugCountRXNeg/(double)countRXNeg)+" buggy rRup when rX<0");
				}
			}
		}
	}
	
	private static void test302() throws IOException {
		List<DiscretizedFunc> funcs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		CPT cpt = GMT_CPT_Files.CATEGORICAL_TAB10.instance();
		double delta = 1d;
		int num = 10;
		
//		CPT cpt = GMT_CPT_Files.CATEGORICAL_TAB10_LIGHT.instance();
//		double delta = 1d;
//		int num = 10;
		
//		CPT cpt = GMT_CPT_Files.CATEGORICAL_TAB20.instance();
//		double delta = 0.5d;
//		int num = 20;
		
		for (int i=0; i<num; i++) {
			double val = i*delta;
			
			ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
			func.set(0d, 0d);
			func.set(1d, val+1);
			
			Color color = cpt.getColor((float)val);
			funcs.add(func);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 5f, color));
		}
		
		PlotSpec spec = new PlotSpec(funcs, chars, " ", "X", "Y");
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(spec);
		
		PlotUtils.writePlots(new File("/tmp"), "cpt_test", gp, 800, 800, true, false, false);
	}
	
	private static void test303() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(
				new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
						+ "2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_gridded.zip"));
		
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
		erf.setParameter(BackgroundRupParam.NAME, BackgroundRupType.POINT);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.setParameter(UseRupMFDsParam.NAME, false);
		erf.updateForecast();
		
		System.out.println("ERF has "+erf.getTotNumRups()+" ruptures");
		System.out.println("ERF has "+erf.getTotNumRupsFromFaultSystem()+" fss ruptures");
		System.out.println("ERF has "+(erf.getTotNumRups() - erf.getTotNumRupsFromFaultSystem())+" gridded ruptures");
		
		System.out.println("WUS Gridded Sites: "+new GriddedRegion(NSHM23_RegionLoader.loadFullConterminousWUS(), 0.1, null).getNodeCount());
		System.out.println("PRVI Gridded Sites: "+new GriddedRegion(PRVI25_RegionLoader.loadPRVI_MapExtents(), 0.1, null).getNodeCount());
	}
	
	private static void test304() throws IOException {
		System.out.println(com.google.common.collect.Range.closed(0d, 5d));
		System.out.println(com.google.common.collect.Range.open(0d, 5d));
		System.out.println(com.google.common.collect.Range.closedOpen(0d, 5d));
		System.out.println(com.google.common.collect.Range.openClosed(0d, 5d));
		System.out.println(com.google.common.collect.Range.atMost(5d));
		System.out.println(com.google.common.collect.Range.atLeast(0d));
	}
	
	private static void test305() throws IOException {
		double minMag = 5.0;
		double maxMag = 8.1999999d;
		boolean magsTenthAligned = false;
		
		IncrementalMagFreqDist refMFD;
		System.out.println("Orig range: ["+minMag+", "+maxMag+"]");
		
		double delta = 0.1;
		if (!magsTenthAligned) {
			// align to 0.x5 bins (so that bin edges are at tenths)
			minMag = Math.floor(minMag*10d)/10d + 0.5*delta;
			maxMag = Math.floor(maxMag*10d)/10d + 0.5*delta;
		}
		int size = (int)Math.round((maxMag - minMag)/delta) + 1;
		refMFD = new IncrementalMagFreqDist(minMag, size, delta);
		
		System.out.println("MFD bin center range: ["+refMFD.getMinX()+", "+refMFD.getMaxX()+"]");
		System.out.println("MFD bin edge range: ["+(float)(refMFD.getMinX()-0.5*refMFD.getDelta())
				+", "+(float)(refMFD.getMaxX()+0.5*refMFD.getDelta())+"]");
	}
	
	private static void test306() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/data/kevin/nshm23/batch_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/results_WUS_FM_v3_branch_averaged_mod_gridded.zip"));
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		System.out.println("Loaded gridProv '"+gridProv.getName()+"' of type "+gridProv.getClass().getName());
	}
	
	private static void test307() throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_06_03-prvi25_subduction_branches/true_mean_solution.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		String prevName = null;
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			if (prevName == null || !sect.getParentSectionName().equals(prevName))
				System.out.println(sect.getParentSectionName());
			prevName = sect.getParentSectionName();
		}
		
		IncrementalMagFreqDist mueIncr = new IncrementalMagFreqDist(5.05d, 46, 0.1);
		IncrementalMagFreqDist carIncr = new IncrementalMagFreqDist(5.05d, 46, 0.1);
		for (FaultSection sect : rupSet.getFaultSectionDataList()) {
			IncrementalMagFreqDist mfd;
			if (sect.getParentSectionName().toLowerCase().contains("thrust")) {
				mfd = mueIncr;
			} else {
				mfd = carIncr;
			}
			IncrementalMagFreqDist sectMFD = sol.calcNucleationMFD_forSect(sect.getSectionId(), mfd.getMinX(), mfd.getMaxX(), mfd.size());
			for (int i=0; i<sectMFD.size(); i++)
				mfd.add(i, sectMFD.getY(i));
		}
		
		for (boolean mue : new boolean[] {false, true}) {
			List<DiscretizedFunc> funcs = new ArrayList<>();
			List<PlotCurveCharacterstics> chars = new ArrayList<>();
			
			IncrementalMagFreqDist mfd = mue ? mueIncr : carIncr;
			int faultMMinIndex = -1;
			for (int i=0; i<mfd.size(); i++) {
				if (mfd.getY(i) > 0d) {
					faultMMinIndex = i;
					break;
				}
			}
			double[] rates;
			double[] bs;
			if (mue) {
				rates = new double[] {0.019, 0.219, 0.735};
				bs = new double[] {1.21, 0.97, 0.73};
			} else {
				rates = new double[] {0.83, 1.99, 4.02};
				bs = new double[] {1.21, 0.97, 0.73};
			}
			GutenbergRichterMagFreqDist[] griddedMFDs = new GutenbergRichterMagFreqDist[rates.length];
			for (int i=0; i<rates.length; i++)
				griddedMFDs[i] = new GutenbergRichterMagFreqDist(bs[i], rates[i],
						mfd.getMinX(), mfd.getX(faultMMinIndex), faultMMinIndex+1);
			UncertainBoundedIncrMagFreqDist griddedMFD = new UncertainBoundedIncrMagFreqDist(
					griddedMFDs[1], griddedMFDs[0], griddedMFDs[2], UncertaintyBoundType.CONF_95);
			
			griddedMFD.setName("Gridded "+griddedMFD.getDefaultBoundName());
			funcs.add(griddedMFD);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Colors.tab_lightblue));
			griddedMFDs[1].setName("Gridded Preferred");
			funcs.add(griddedMFDs[1]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_blue));
			
			mfd.setName("Inversion MFD");
			funcs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Colors.tab_orange));
			
			PlotSpec spec = new PlotSpec(funcs, chars, mue ? "Muertos" : "CAR", "Magnitude", "Incremental Rate (1/yr)");
			spec.setLegendInset(RectangleAnchor.TOP_RIGHT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(spec, false, true, new Range(5d, 9.5d), new Range(1e-6, 1));
			
			PlotUtils.writePlots(new File("/tmp"), mue ? "mue_mfds" : "car_mfds", gp, 800, 800, true, false, false);
		}
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		try {
			test307();
		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(1);
		}
	}

}
