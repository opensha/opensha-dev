package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.jfree.data.Range;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.targetMFDs.SupraSeisBValInversionTargetMFDs;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.model.HazardModel;
import gov.usgs.earthquake.nshmp.model.NshmErf;
import scratch.kevin.nshm23.SingleSiteHazardAndDataComparisonPageGen;
import scratch.kevin.nshm23.SingleSiteHazardAndDataComparisonPageGen.RegionalParticipationResult;

public class WUS_EastOfAttenBoundaryMFDs {

	public static void main(String[] args) throws IOException {
		Region stableReg = NSHM23_RegionLoader.GridSystemRegions.CEUS_STABLE.load();
		Region wusReg = NSHM23_RegionLoader.loadFullConterminousWUS();
		
		File outputDir = new File("/tmp");
		
		Region stableWUS = Region.intersect(stableReg, wusReg);
		
		Feature.write(stableWUS.toFeature(), new File(outputDir, "wus_reg_east_of_atten.geojson"));
		
		double extent1 = stableReg.getExtent();
		double extent2 = wusReg.getExtent();
		double intersectExtent = stableWUS.getExtent();
		
		System.out.println("Intersection is "+(float)(intersectExtent/extent2)+" fraction of WUS");
		
		GriddedRegion gridReg = new GriddedRegion(stableWUS, 0.1, GriddedRegion.ANCHOR_0_0);
		
		Set<TectonicRegionType> trts = EnumSet.of(TectonicRegionType.ACTIVE_SHALLOW, TectonicRegionType.STABLE_SHALLOW);
		
		HazardModel model = HazardModel.load(Path.of("/home/kevin/OpenSHA/nshm23/nshmp-haz-models/nshm-conus-5.3.0"));
		NshmErf faultERF = new NshmErf(model, trts, IncludeBackgroundOption.EXCLUDE);
		System.out.println("NSHM Fault ERF size: " + faultERF.getNumSources());
		faultERF.getTimeSpan().setDuration(1.0);
		faultERF.updateForecast();
		NshmErf gridERF = new NshmErf(model, trts, IncludeBackgroundOption.ONLY);
		System.out.println("NSHM Grid ERF size: " + gridERF.getNumSources());
		gridERF.getTimeSpan().setDuration(1.0);
		gridERF.updateForecast();
		
		EvenlyDiscretizedFunc refMFD = SupraSeisBValInversionTargetMFDs.buildRefXValues(8.05);
		
		Range xRange = new Range(5d,  8d);
		Range yRange = new Range(1e-6, 1e0);
		
		double[] magThresholds = {5d, 6d, 7d};
		
		RegionalParticipationResult nshm18FaultPartic = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				faultERF, gridReg, magThresholds, refMFD);
		RegionalParticipationResult nshm18GridPartic = SingleSiteHazardAndDataComparisonPageGen.calcNshmErfPartic(
				gridERF, gridReg, magThresholds, refMFD);
		
		FaultSystemSolution fss = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2023_04_11-nshm23_branches-NSHM23_v2-CoulombRupSet-TotNuclRate-NoRed-ThreshAvgIterRelGR/"
				+ "results_NSHM23_v2_CoulombRupSet_branch_averaged_gridded.zip"));
		
		RegionalParticipationResult nshm23FaultPartic = SingleSiteHazardAndDataComparisonPageGen.calcFSSFaultPartic(
				fss, gridReg, magThresholds, refMFD);
		RegionalParticipationResult nshm23GridPartic = SingleSiteHazardAndDataComparisonPageGen.calcFSSGriddedPartic(
				fss, gridReg, magThresholds, refMFD);
		
		List<DiscretizedFunc> incrFuncs = new ArrayList<>();
		List<DiscretizedFunc> cmlFuncs = new ArrayList<>();
		List<PlotCurveCharacterstics> chars = new ArrayList<>();
		
		IncrementalMagFreqDist mfd18Fault = nshm18FaultPartic.particMFD;
		mfd18Fault.setName("NSHM18, Fault");
		IncrementalMagFreqDist mfd18Grid = nshm18GridPartic.particMFD;
		mfd18Grid.setName("NSHM18, Grid");
		IncrementalMagFreqDist mfd18Sum = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		mfd18Sum.setName("NSHM18, Sum");
		for (int i=0; i<refMFD.size(); i++)
			mfd18Sum.set(i, mfd18Fault.getY(i)+mfd18Grid.getY(i));
		
		IncrementalMagFreqDist mfd23Fault = nshm23FaultPartic.particMFD;
		mfd23Fault.setName("NSHM23, Fault");
		IncrementalMagFreqDist mfd23Grid = nshm23GridPartic.particMFD;
		mfd23Grid.setName("NSHM23, Grid");
		IncrementalMagFreqDist mfd23Sum = new IncrementalMagFreqDist(refMFD.getMinX(), refMFD.size(), refMFD.getDelta());
		mfd23Sum.setName("NSHM23, Sum");
		for (int i=0; i<refMFD.size(); i++)
			mfd23Sum.set(i, mfd23Fault.getY(i)+mfd23Grid.getY(i));
		
		incrFuncs.add(mfd18Fault);
		cmlFuncs.add(mfd18Fault.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.GRAY));
		
		incrFuncs.add(mfd23Fault);
		cmlFuncs.add(mfd23Fault.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.GRAY));
		
		incrFuncs.add(mfd18Grid);
		cmlFuncs.add(mfd18Grid.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLUE));
		
		incrFuncs.add(mfd23Grid);
		cmlFuncs.add(mfd23Grid.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		incrFuncs.add(mfd18Sum);
		cmlFuncs.add(mfd18Sum.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 3f, Color.BLACK));
		
		incrFuncs.add(mfd23Sum);
		cmlFuncs.add(mfd23Sum.getCumRateDistWithOffset());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		PlotSpec incrSpec = new PlotSpec(incrFuncs, chars, "WUS East of Attenuation Boundary", "Magnitude",
				"Incremental Participation Rate (1/yr)");
		incrSpec.setLegendInset(true);
		
		PlotSpec cmlSpec = new PlotSpec(cmlFuncs, chars, "WUS East of Attenuation Boundary", "Magnitude",
				"Cumulative Participation Rate (1/yr)");
		cmlSpec.setLegendInset(true);
		
		HeadlessGraphPanel gp = PlotUtils.initHeadless();
		
		gp.drawGraphPanel(incrSpec, false, true, xRange, yRange);
		PlotUtils.writePlots(outputDir, "wus_mfds_east_of_atten", gp, 1000, 850, true, true, true);
		
		gp.drawGraphPanel(cmlSpec, false, true, xRange, yRange);
		PlotUtils.writePlots(outputDir, "wus_mfds_east_of_atten_cml", gp, 1000, 850, true, true, true);
	}

}
