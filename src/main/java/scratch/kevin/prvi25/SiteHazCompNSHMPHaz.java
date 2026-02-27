package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.data.WeightedValue;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotUtils;
import org.opensha.commons.util.FileNameUtils;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.calc.sourceFilters.SourceFilterManager;
import org.opensha.sha.calc.sourceFilters.SourceFilters;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSectionUtils;
import org.opensha.sha.earthquake.param.AseismicityAreaReductionParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.UseProxySectionsParam;
import org.opensha.sha.earthquake.param.UseRupMFDsParam;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.GroundMotionLogicTreeFilter;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;

import gov.usgs.earthquake.nshmp.gmm.Gmm;
import gov.usgs.earthquake.nshmp.gmm.GmmInput;
import gov.usgs.earthquake.nshmp.gmm.GroundMotion;
import gov.usgs.earthquake.nshmp.gmm.GroundMotions;
import gov.usgs.earthquake.nshmp.gmm.UsgsPrviBackbone2025;
import gov.usgs.earthquake.nshmp.gmm.GmmInput.Constraints;
import gov.usgs.earthquake.nshmp.gmm.GmmInput.Field;
import gov.usgs.earthquake.nshmp.tree.LogicTree;
import net.mahdilamb.colormap.Colors;

public class SiteHazCompNSHMPHaz {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/prvi_site_nshmp_haz_comp");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
				"Output directory %s does not exist and could not be created", outputDir.getAbsolutePath());
		
//		NSHMP_GMM_Wrapper gmm = new NSHMP_GMM_Wrapper.Single(Gmm.ASK_14_BASE);
//		CSVFile<String> inCSV = CSVFile.readFile(new File("/home/kevin/Downloads/"
//				+ "prvi-0.2s-ask14base-crustalFaultOnly-PRVI_2025_ACTIVE_CRUST_NO_EPI_SIGMA_NGA.csv"), true);
		
//		NSHMP_GMM_Wrapper gmm = new NSHMP_GMM_Wrapper.Single(Gmm.PRVI_2025_ACTIVE_CRUST);
////		gmm.setGroundMotionTreeFilter(new GroundMotionLogicTreeFilter.StringMatching(
//////				UsgsPrviBackbone2025.SIGMA_NGA_ID,
////				GroundMotions.EPI_OFF
////				));
//		CSVFile<String> inCSV = CSVFile.readFile(new File("/home/kevin/Downloads/"
//				+ "prvi-0.2s-prvi25active-crustalFaultOnly-curves.csv"), true);
		
		NSHMP_GMM_Wrapper gmm = new NSHMP_GMM_Wrapper.WeightedCombination(
				WeightedList.of(new WeightedValue<>(Gmm.PRVI_2025_ACTIVE_CRUST, 0.5),
				new WeightedValue<>(Gmm.PRVI_2025_ACTIVE_CRUST_ADJUSTED, 0.5)), "Name", "Name");
//		NSHMP_GMM_Wrapper gmm = new NSHMP_GMM_Wrapper.Single(Gmm.TOTAL_TREE_PRVI_ACTIVE_CRUST_2025);
//		gmm.setGroundMotionTreeFilter(new GroundMotionLogicTreeFilter.StringMatching(
//				UsgsPrviBackbone2025.SIGMA_NGA_ID,
//				GroundMotions.EPI_OFF
//				));
		CSVFile<String> inCSV = CSVFile.readFile(new File("/home/kevin/Downloads/"
				+ "prvi-0.2s-prvi25activeTotal-crustalFaultOnly-curves.csv"), true);
		
		gmm.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), 0.2);
		String xName = "0.2s SA (g)";
//		gmm.setIntensityMeasure(PGA_Param.NAME);
//		String xName = "PGA (g)";
		gmm.getOtherParams().setValue(SigmaTruncTypeParam.NAME, SigmaTruncTypeParam.SIGMA_TRUNC_TYPE_1SIDED);
		gmm.getOtherParams().setValue(SigmaTruncLevelParam.NAME, 3d);
		
		
		System.out.println("Zhype used ? Gmm.PRVI_2025_ACTIVE_CRUST: "+Gmm.PRVI_2025_ACTIVE_CRUST.constraints().get(Field.ZHYP).isPresent());
		System.out.println("Zhype used ? Gmm.PRVI_2025_ACTIVE_CRUST: "+Gmm.TOTAL_TREE_PRVI_ACTIVE_CRUST_2025.constraints().get(Field.ZHYP).isPresent());
		
		if (gmm instanceof NSHMP_GMM_Wrapper.Single) {
			Gmm singleGMM = ((NSHMP_GMM_Wrapper.Single)gmm).getGmmRef();
			Constraints constraints = singleGMM.constraints();
			for (Field field : Field.values()) {
				if (constraints.get(field).isPresent()) {
					// this field is used
					System.out.println(singleGMM.name()+" reports that "+field+" is used");
				} else {
					System.out.println(singleGMM.name()+" reports that "+field+" is unused");
				}
			}
		}
		
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2025_08_01-prvi25_crustal_branches-dmSample10x/"
				+ "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded_simplified.zip"));
//				+ "results_PRVI_CRUSTAL_FM_V1p1_branch_averaged_gridded.zip"));
		BaseFaultSystemSolutionERF erf = new BaseFaultSystemSolutionERF();
		erf.setSolution(sol);
		erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
//		erf.setParameter(UseProxySectionsParam.NAME, false);
//		erf.setParameter(AseismicityAreaReductionParam.NAME, false);
		erf.setParameter(UseRupMFDsParam.NAME, false);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		
		FaultSystemRupSet rupSet = sol.getRupSet();
		int testParentID = FaultSectionUtils.findParentSectionID(rupSet.getFaultSectionDataList(), "South Lajas");
		int testRupID = -1;
		double testRupMag = 0d;
		for (int r : rupSet.getRupturesForParentSection(testParentID)) {
			boolean allMatch = true;
			for (FaultSection sect : rupSet.getFaultSectionDataForRupture(r)) {
				if (sect.getParentSectionId() != testParentID) {
					allMatch = false;
					break;
				}
			}
			if (allMatch && rupSet.getMagForRup(r) > testRupMag) {
				testRupID = r;
				testRupMag = rupSet.getMagForRup(r);
			}
		}
		int testSourceID = erf.getSrcIndexForFltSysRup(testRupID);
		ProbEqkSource testSource = erf.getSource(testSourceID);
		
		HazardCurveCalculator calc = new HazardCurveCalculator(new SourceFilterManager(SourceFilters.TRT_DIST_CUTOFFS));
		
		DiscretizedFunc xVals = new IMT_Info().getDefaultHazardCurve(gmm.getIntensityMeasure());
		DiscretizedFunc logXvalues = new ArbitrarilyDiscretizedFunc();
		for (int i = 0; i < xVals.size(); i++)
			logXvalues.set(Math.log(xVals.getX(i)), 1d);
		
		Range xRange = new Range(1e-4, 3e0);
		Range yRange = new Range(1e-5, 1e-1);
		
		PlotCurveCharacterstics theirChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_orange);
		PlotCurveCharacterstics ourChar = new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Colors.tab_blue);
		
		for (int row=1; row<inCSV.getNumRows(); row++) {
			String name = inCSV.get(row, 0);
			double lat = inCSV.getDouble(row, 2);
			double lon = inCSV.getDouble(row, 1);
			Site site = new Site(new Location(lat, lon), name);
			site.addParameterList(gmm.getSiteParams());
			
			DiscretizedFunc theirs = new ArbitrarilyDiscretizedFunc();
			for (int col=3; col<inCSV.getNumCols(); col++) {
				double x = inCSV.getDouble(0, col);
				double y = inCSV.getDouble(row, col);
				y = 1d - Math.exp(-y);
				theirs.set(x, y);
			}
			theirs.setName("NSHMP-Haz");
			
			System.out.println("Site:\t"+name);
			gmm.setSite(site);
			for (int r=0; r<testSource.getNumRuptures(); r++) {
				ProbEqkRupture rup = testSource.getRupture(r);
				gmm.setEqkRupture(rup);
				GmmInput gmmInput = gmm.getCurrentGmmInput();
				LogicTree<GroundMotion> gmmTree = gmm.getGroundMotionTree();
				System.out.println("FSS rup "+testRupID+", ERF source "+testSourceID+" ("+r+")");
				System.out.println(gmmInput);
				System.out.println(gmmTree);
			}
			System.out.println();
			
			calc.getHazardCurve(logXvalues, site, gmm, erf);
			
			DiscretizedFunc ours = xVals.deepClone();
			for (int i=0; i<xVals.size(); i++)
				ours.set(i, logXvalues.getY(i));
			ours.setName("OpenSHA");
			
			PlotSpec plot = new PlotSpec(List.of(theirs, ours), List.of(theirChar, ourChar), name, xName, "Annual Probability of Exceedance");
			plot.setLegendInset(RectangleAnchor.BOTTOM_LEFT);
			
			HeadlessGraphPanel gp = PlotUtils.initHeadless();
			
			gp.drawGraphPanel(plot, true, true, xRange, yRange);
			
			PlotUtils.writePlots(outputDir, FileNameUtils.simplify(name)+"_hazard_curves", gp, 800, 800, true, true, true);
		}
	}

}
