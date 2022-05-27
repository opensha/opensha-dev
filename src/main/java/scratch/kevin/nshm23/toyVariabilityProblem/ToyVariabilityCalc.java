package scratch.kevin.nshm23.toyVariabilityProblem;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;
import org.jfree.data.Range;
import org.opensha.commons.calc.FaultMomentCalc;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.uncertainty.UncertainBoundedIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertainIncrMagFreqDist;
import org.opensha.commons.data.uncertainty.UncertaintyBoundType;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.commons.util.ReturnPeriodUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets;
import org.opensha.sha.earthquake.faultSysSolution.RuptureSets.SimpleAzimuthalRupSetConfig;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration;
import org.opensha.sha.earthquake.faultSysSolution.inversion.InversionConfiguration.Builder;
import org.opensha.sha.earthquake.faultSysSolution.inversion.Inversions;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.ConstraintWeightingType;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.InversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.MFDInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SectionTotalRateConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.constraints.impl.SlipRateInversionConstraint;
import org.opensha.sha.earthquake.faultSysSolution.inversion.sa.completion.IterationCompletionCriteria;
import org.opensha.sha.earthquake.faultSysSolution.modules.AveSlipModule;
import org.opensha.sha.earthquake.faultSysSolution.modules.ModelRegion;
import org.opensha.sha.earthquake.faultSysSolution.modules.SectSlipRates;
import org.opensha.sha.earthquake.faultSysSolution.modules.SlipAlongRuptureModel;
import org.opensha.sha.earthquake.faultSysSolution.util.BranchAverageSolutionCreator;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.SupraSeisBValues;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs.SubSeisMoRateReduction;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.mod.ModAttenRelRef;
import org.opensha.sha.imr.mod.ModAttenuationRelationship;
import org.opensha.sha.imr.mod.impl.FixedStdDevMod;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.erf.FaultSystemSolutionERF;

public class ToyVariabilityCalc {

	public static void main(String[] args) throws IOException {
		ScalingRelationships scale = ScalingRelationships.ELLSWORTH_B;
		double ddw = 10; // km
		int minSubSectsPerRup = 4;
		double targetMinMag = 6.55;
		double targetMaxMag = 7.75;
		double slipRate = 50d; // mm/yr
		double slipRateStdDev = slipRate*0.1;
		
		boolean rerunInversions = false;
		File outputDir = new File("/home/kevin/OpenSHA/UCERF4/classic_model_comparisons");
		
		double areaForMin = scale.getArea(targetMinMag, ddw*1e3)*1e-6;
		double lenForMin = areaForMin/ddw;
		System.out.println("Length for a M"+(float)targetMinMag+": "+lenForMin);
		double areaForMax = scale.getArea(targetMaxMag, ddw*1e3)*1e-6;
		double lenForMax = areaForMax/ddw;
		System.out.println("Length for a M"+(float)targetMaxMag+": "+lenForMax);
		double subSectLen = lenForMin/(double)minSubSectsPerRup;
		int numSubSects = (int)(lenForMax/subSectLen + 0.5);
		System.out.println("Num sub sects: "+numSubSects);
		double finalLen = numSubSects*subSectLen;
		System.out.println("Final length: "+finalLen);
//		double lenFor7p8 
		Location fltLoc1 = new Location(0d, 0d);
		Location fltLoc2 = LocationUtils.location(fltLoc1, Math.PI/2d, finalLen);
		fltLoc2 = new Location(fltLoc1.getLatitude(), fltLoc2.getLongitude());
		
		Location hazardLoc = LocationUtils.location(fltLoc1, Math.PI/2d, 0.5*finalLen);
		double[] sigmas = {0.3, 0.5, 0.7};
		Color[] sigmaColors = {Color.GRAY, Color.DARK_GRAY, Color.BLACK};
//		double fractGR = 1d/3d;
		double fractGR = 0.5;
		boolean evenWeight = true;
		
		ModAttenuationRelationship modAttenRel = new ModAttenuationRelationship(AttenRelRef.ASK_2014,
				ModAttenRelRef.FIXED_STD_DEV);
		modAttenRel.setParamDefaults();
		modAttenRel.setIntensityMeasure(PGA_Param.NAME);
		Site site = new Site(hazardLoc);
		site.addParameterList(modAttenRel.getSiteParams());
		
		Region modelRegion = new Region(new Location(fltLoc1.getLatitude()-1d, fltLoc1.getLongitude()-1d),
				new Location(fltLoc2.getLatitude()+1, fltLoc2.getLongitude()+1));
		
		String sectJSON =
				"    {\n"+
				"      \"type\": \"Feature\",\n"+
				"      \"id\": 0,\n"+
				"      \"properties\": {\n"+
				"        \"FaultID\": 0,\n"+
				"        \"FaultName\": \"Test Fault\",\n"+
				"        \"DipDeg\": 90.0,\n"+
				"        \"Rake\": 180.0,\n"+
				"        \"LowDepth\": "+(float)ddw+",\n"+
				"        \"UpDepth\": 0.0,\n"+
				"        \"SlipRate\": "+(float)slipRate+",\n"+
				"        \"SlipRateStdDev\": "+(float)slipRateStdDev+"\n"+
				"      },\n"+
				"      \"geometry\": {\n"+
				"        \"type\": \"LineString\",\n"+
				"        \"coordinates\": [\n"+
				"          [\n"+
				"            "+fltLoc1.getLongitude()+",\n"+
				"            "+fltLoc1.getLatitude()+"\n"+
				"          ],\n"+
				"          [\n"+
				"            "+fltLoc2.getLongitude()+",\n"+
				"            "+fltLoc2.getLatitude()+"\n"+
				"          ]\n"+
				"        ]\n"+
				"      }\n"+
				"    }";
		System.out.println(sectJSON);
		FaultSection sect = GeoJSONFaultSection.fromFeature(Feature.fromJSON(sectJSON));
		
		// sub-section it
		List<? extends FaultSection> subSects = sect.getSubSectionsList(subSectLen+0.1, 0);
		System.out.println("Built "+subSects.size()+" sub sects");
		
		SimpleAzimuthalRupSetConfig rsConfig = new RuptureSets.SimpleAzimuthalRupSetConfig(subSects, scale);
		rsConfig.setMinSectsPerParent(minSubSectsPerRup);
		FaultSystemRupSet rupSet = rsConfig.build(1);
		
		rupSet.addModule(new ModelRegion(modelRegion));
		
		System.out.println("build "+rupSet.getNumRuptures()+" ruptures. Mag range: "+rupSet.getMinMag()+", "+rupSet.getMaxMag());
		
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(rupSet);
		int minMagIndex = refMFD.getClosestXIndex(rupSet.getMinMag());
		int maxMagIndex = refMFD.getClosestXIndex(rupSet.getMaxMag());
		int[] mfdRupCounts = new int[refMFD.size()];
		for (int i=minMagIndex; i<=maxMagIndex; i++) {
			double mag = refMFD.getX(i);
			if (mag < 6.5d)
				continue;
			int matchingRups = 0;
			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(r));
				if (magIndex == i)
					matchingRups++;
			}
			mfdRupCounts[i] = matchingRups;
			System.out.println("M="+(float)mag+", "+matchingRups+" rups");
		}
		
		// build classic model
		double totMoment = sect.calcMomentRate(false);
		// construct a G-R that satisfies this model
		GutenbergRichterMagFreqDist grMFD = new GutenbergRichterMagFreqDist(refMFD.getMinX(), refMFD.size(),
				refMFD.getDelta(), refMFD.getX(minMagIndex), refMFD.getX(maxMagIndex), totMoment, 1d);
		System.out.println("G-R b=1 that satisfies:\n"+grMFD);
		int charRupIndex = -1;
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			if (rupSet.getSectionsIndicesForRup(r).size() == subSects.size()) {
				charRupIndex = r;
				break;
			}
		}
		System.out.println("Char rup index: "+charRupIndex+" (M="+rupSet.getMagForRup(charRupIndex)+")");
		double charEventAveSlip = rupSet.requireModule(AveSlipModule.class).getAveSlip(charRupIndex); // m
		double charEventMoment = FaultMomentCalc.getMoment(sect.getArea(false), charEventAveSlip);
		double charRateToSatisfy = totMoment/charEventMoment;
		System.out.println("Char slip: "+charEventAveSlip+", rate to satisfy: "+charRateToSatisfy);
		IncrementalMagFreqDist grTarget = grMFD.deepClone();
		grTarget.scaleToTotalMomentRate(totMoment*fractGR);
		double[] classicRates = new double[rupSet.getNumRuptures()];
		for (int r=0; r<rupSet.getNumRuptures(); r++) {
			int magIndex = refMFD.getClosestXIndex(rupSet.getMagForRup(r));
			classicRates[r] = grTarget.getY(magIndex)/(double)mfdRupCounts[magIndex];
		}
		classicRates[charRupIndex] += (1d-fractGR)*charRateToSatisfy;
		FaultSystemSolution classicSol = new FaultSystemSolution(rupSet, classicRates);
		IncrementalMagFreqDist classicSolMFD = classicSol.calcNucleationMFD_forRegion(
				null, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
		System.out.println("Solution MFD:\n"+classicSolMFD);
		
		classicSol.write(new File(outputDir, "classic_sol.zip"));
		
		List<FaultSystemSolution> invSols = new ArrayList<>();
		List<Double> invSolWeights = new ArrayList<>();
		List<Double> bVals = new ArrayList<>();
		BranchAverageSolutionCreator baCreate = new BranchAverageSolutionCreator(
				evenWeight ? new BranchWeightProvider.ConstantWeights() : new BranchWeightProvider.CurrentWeights());
		for (SupraSeisBValues bVal : SupraSeisBValues.values()) {
			File solFile = new File(outputDir, "inversion_"+bVal.getFilePrefix()+".zip");
			LogicTreeBranch<LogicTreeNode> branch = new LogicTreeBranch<>(List.of(LogicTreeLevel.forEnumUnchecked(
					SupraSeisBValues.class, "Supra-Seis b-value", "BVal")));
			branch.setValue(bVal);
			bVals.add(bVal.bValue);
			FaultSystemSolution invSol;
			if (!rerunInversions && solFile.exists()) {
				invSol = FaultSystemSolution.load(solFile);
			} else {
				SupraSeisBValInversionTargetMFDs targetMFDs = new SupraSeisBValInversionTargetMFDs.Builder(rupSet, bVal.bValue)
						.subSeisMoRateReduction(SubSeisMoRateReduction.NONE)
						.applyDefModelUncertainties(true)
						.magDepDefaultRelStdDev(M->0.1*Math.pow(10, bVal.bValue*0.5*(M-6)))
						.build();
				List<InversionConstraint> constraints = new ArrayList<>();
				constraints.add(new SlipRateInversionConstraint(1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY,
						rupSet, rupSet.requireModule(AveSlipModule.class), rupSet.requireModule(SlipAlongRuptureModel.class),
						rupSet.requireModule(SectSlipRates.class)));
				constraints.add(new MFDInversionConstraint(rupSet, 10d, false,
						ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY, targetMFDs.getMFD_Constraints()));
				double[] targetRates = new double[rupSet.getNumSections()];
				double[] targetRateStdDevs = new double[rupSet.getNumSections()];
				
				List<UncertainIncrMagFreqDist> sectSupraMFDs = targetMFDs.getOnFaultSupraSeisNucleationMFDs();
				for (int s=0; s<targetRates.length; s++) {
					UncertainIncrMagFreqDist sectSupraMFD = sectSupraMFDs.get(s);
					targetRates[s] = sectSupraMFD.calcSumOfY_Vals();
					UncertainBoundedIncrMagFreqDist oneSigmaBoundedMFD = sectSupraMFD.estimateBounds(UncertaintyBoundType.ONE_SIGMA);
					double upperVal = oneSigmaBoundedMFD.getUpper().calcSumOfY_Vals();
					double lowerVal = oneSigmaBoundedMFD.getLower().calcSumOfY_Vals();
					targetRateStdDevs[s] = UncertaintyBoundType.ONE_SIGMA.estimateStdDev(targetRates[s], lowerVal, upperVal);
//					System.out.println(rupSet.getFaultSectionData(s).getSectionName()+": totRate="+(float)targetRates[s]
//							+"\tstdDev="+(float)targetRateStdDevs[s]+"\trelStdDev="+(float)(targetRateStdDevs[s]/targetRates[s]));
				}
				constraints.add(new SectionTotalRateConstraint(rupSet, 1d, ConstraintWeightingType.NORMALIZED_BY_UNCERTAINTY,
						targetRates, targetRateStdDevs, true));
				Builder configBuilder = InversionConfiguration.builder(constraints,
						new IterationCompletionCriteria(rupSet.getNumRuptures()*10000));
				configBuilder.threads(16).reweight();
				configBuilder.subCompletion(new IterationCompletionCriteria(rupSet.getNumRuptures()*10));
				configBuilder.avgThreads(4, new IterationCompletionCriteria(rupSet.getNumRuptures()*100));
				InversionConfiguration config = configBuilder.build();
				
				invSol = Inversions.run(rupSet, config);
				invSol.addModule(branch);
				invSol.write(solFile);
			}
			baCreate.addSolution(invSol, branch);
			invSols.add(invSol);
			invSolWeights.add(bVal.getNodeWeight(null));
		}
		
		FaultSystemSolution baSol = baCreate.build();
		baSol.write(new File(outputDir, "inversion_branch_averaged.zip"));
		
		List<IncrementalMagFreqDist> incrMFDs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		for (int i=0; i<invSols.size(); i++) {
			FaultSystemSolution invSol = invSols.get(i);
			IncrementalMagFreqDist mfd = invSol.calcNucleationMFD_forRegion(
					null, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
			double b = bVals.get(i);
			mfd.setName("b="+(float)b);
			
			incrMFDs.add(mfd);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		}
		
		IncrementalMagFreqDist baMFD = baSol.calcNucleationMFD_forRegion(
				null, refMFD.getMinX(), refMFD.getMaxX(), refMFD.getDelta(), false);
		baMFD.setName("Average Inversion");
		
		incrMFDs.add(baMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		
		classicSolMFD.setName("Classic");
		incrMFDs.add(classicSolMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(incrMFDs, chars, "MFD Comparison", "Magnitude", "Incremental Rate (1/yr)");
		spec.setLegendVisible(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.setLegendFontSize(18);
		
		Range xRange = new Range(refMFD.getX(minMagIndex)-0.5*refMFD.getDelta(), refMFD.getX(maxMagIndex)+0.5*refMFD.getDelta());
		Range yRange = new Range(1e-4, 1e-1);
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "mfds", gp, 800, 750, true, true, false);
		
		List<EvenlyDiscretizedFunc> cmlMFDs = new ArrayList<>();
		for (IncrementalMagFreqDist mfd : incrMFDs)
			cmlMFDs.add(mfd.getCumRateDistWithOffset());
		spec = new PlotSpec(cmlMFDs, chars, "MFD Comparison", "Magnitude", "Cumulative Rate (1/yr)");
		spec.setLegendVisible(true);
		yRange = new Range(1e-3, 1e0);
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "mfds_cumulative", gp, 800, 750, true, true, false);
		
		// classic only MFD plot
		incrMFDs.clear();
		chars.clear();
		
		classicSolMFD.setName("Classic");
		incrMFDs.add(classicSolMFD);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		spec = new PlotSpec(incrMFDs, chars, "Classic Model", "Magnitude", "Incremental Rate (1/yr)");
		spec.setLegendVisible(false);
		
		yRange = new Range(1e-4, 1e-1);
		gp.drawGraphPanel(spec, false, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "classic_mfd", gp, 800, 750, true, true, false);
		
		// hazard curve plots
		
		// first sigma dependence with mean model
		
		
		List<DiscretizedFunc> funcs = new ArrayList<>();
		chars = new ArrayList<>();
		
		for (int s=0; s<sigmas.length; s++) {
			double sigma = sigmas[s];
			((FixedStdDevMod)modAttenRel.getCurrentMod()).setStdDev(sigma);
			
			DiscretizedFunc baCurve = hazardCurve(baSol, modAttenRel, site);
			
			baCurve.setName("GMM σ="+(float)sigma);
			funcs.add(baCurve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, sigmaColors[s]));
		}
		
		spec = new PlotSpec(funcs, chars, "GMM Sigma-Dependence", "Peak Ground Acceleration (g)",
				"Annual Probability of Exceedance");
		spec.setLegendInset(true);
		
		gp.setLegendFontSize(24);
		
		xRange = new Range(1e-3, 1e1);
		yRange = new Range(1e-6, 1e0);
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "curves_sigma", gp, 800, 750, true, true, false);
		
		// now model curves
		((FixedStdDevMod)modAttenRel.getCurrentMod()).setStdDev(StatUtils.max(sigmas));
		funcs = new ArrayList<>();
		chars = new ArrayList<>();
		
		for (FaultSystemSolution invSol : invSols) {
			DiscretizedFunc curve = hazardCurve(invSol, modAttenRel, site);
			if (funcs.isEmpty())
				curve.setName("Individual Inversions");
			
			funcs.add(curve);
			chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		}
		
		DiscretizedFunc baCurve = hazardCurve(baSol, modAttenRel, site);
		baCurve.setName("Average Inversion");
		
		funcs.add(baCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.RED));
		
		DiscretizedFunc classicCurve = hazardCurve(classicSol, modAttenRel, site);
		classicCurve.setName("Classic Model");
		funcs.add(classicCurve);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 4f, Color.BLUE));
		
		spec = new PlotSpec(funcs, chars, "Hazard Comparison", "Peak Ground Acceleration (g)",
				"Annual Probability of Exceedance");
		spec.setLegendInset(true);
		
		gp.drawGraphPanel(spec, true, true, xRange, yRange);
		
		PlotUtils.writePlots(outputDir, "curves_sols", gp, 800, 750, true, true, false);
	}
	
	private static DiscretizedFunc hazardCurve(FaultSystemSolution sol, ScalarIMR gmpe, Site site) {
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (Point2D pt : xVals)
			logXVals.set(Math.log(pt.getX()), 1d);
		
		HazardCurveCalculator calc = new HazardCurveCalculator();
		
		calc.getHazardCurve(logXVals, site, gmpe, erf);
		
		DiscretizedFunc curve = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			curve.set(xVals.getX(i), logXVals.getY(i));
		
		return curve;
	}

}
