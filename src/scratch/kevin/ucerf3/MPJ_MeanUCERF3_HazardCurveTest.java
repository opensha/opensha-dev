package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.hpc.mpj.taskDispatch.MPJTaskCalculator;
import org.opensha.sha.calc.HazardCurveCalculator;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.inversion.InversionFaultSystemSolution;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class MPJ_MeanUCERF3_HazardCurveTest extends MPJTaskCalculator {
	
	private Compound_FM_DM_Scale_CombinedFetcher fetch;
	private List<LogicTreeBranch> branches;
	
	private ScalarIMR imr;
	private Site site;
	private HazardCurveCalculator calc;
	private DiscretizedFunc xVals;
	
	private File poissonDir;
	private File timeDepDir;

	public MPJ_MeanUCERF3_HazardCurveTest(CommandLine cmd, CompoundFaultSystemSolution cfss, File dir) {
		super(cmd);
		
		fetch = new Compound_FM_DM_Scale_CombinedFetcher(cfss, new APrioriBranchWeightProvider());
		branches = Lists.newArrayList(fetch.getBranches());
		Collections.sort(branches);
		
		poissonDir = new File(dir, "curves_poisson");
		if (rank == 0 && !poissonDir.exists())
			poissonDir.mkdir();
		timeDepDir = new File(dir, "curves_timedep");
		if (rank == 0 && !timeDepDir.exists())
			timeDepDir.mkdir();
		
		imr = AttenRelRef.CB_2008.instance(null);
		imr.setParamDefaults();
		imr.setIntensityMeasure(PGA_Param.NAME);
		xVals = new IMT_Info().getDefaultHazardCurve(PGA_Param.NAME);
		site = new Site(new Location(34.055, -118.2467)); // downtown LA
		site.addParameterList(imr.getSiteParams());
		
		calc = new HazardCurveCalculator();
	}

	@Override
	protected int getNumTasks() {
		return branches.size();
	}

	@Override
	protected void calculateBatch(int[] batch) throws Exception {
		for (int index : batch) {
			calculate(branches.get(index));
		}
	}
	
	private void calculate(LogicTreeBranch branch) throws IOException {
		String fName = branch.getValue(FaultModels.class).encodeChoiceString()
				+"_"+branch.getValue(DeformationModels.class).encodeChoiceString()
				+"_"+branch.getValue(ScalingRelationships.class).encodeChoiceString()+".txt";
		File timeDepCurveFile = new File(timeDepDir, fName);
		File poissonCurveFile = new File(poissonDir, fName);
		
		if (timeDepCurveFile.exists() && poissonCurveFile.exists())
			return;
		
		
		InversionFaultSystemSolution sol = fetch.getSolution(branch);
		FaultSystemSolutionERF erf = new FaultSystemSolutionERF(sol);
		erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, false);
		erf.getParameter(IncludeBackgroundParam.NAME).setValue(IncludeBackgroundOption.EXCLUDE);
		
		// first time dep
		erf.getParameter(ProbabilityModelParam.NAME).setValue(
				ProbabilityModelOptions.U3_PREF_BLEND);
		erf.setParameter(HistoricOpenIntervalParam.NAME,
				(double)(FaultSystemSolutionERF.START_TIME_DEFAULT-1875));
		erf.setParameter(BPTAveragingTypeParam.NAME,
				BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
		
		erf.getTimeSpan().setDuration(30d);
		erf.updateForecast();
		
		DiscretizedFunc hazFunction = xVals.deepClone();
		calc.getHazardCurve(hazFunction, site, imr, erf);
		
		ArbitrarilyDiscretizedFunc.writeSimpleFuncFile(hazFunction, timeDepCurveFile);

		erf.getParameter(ProbabilityModelParam.NAME).setValue(
				ProbabilityModelOptions.POISSON);
		erf.getTimeSpan().setDuration(30d);
		erf.updateForecast();
		
		hazFunction = xVals.deepClone();
		calc.getHazardCurve(hazFunction, site, imr, erf);
		
		ArbitrarilyDiscretizedFunc.writeSimpleFuncFile(hazFunction, poissonCurveFile);
	}

	@Override
	protected void doFinalAssembly() throws Exception {
		// TODO Auto-generated method stub

	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		args = MPJTaskCalculator.initMPJ(args);

		try {
			Options options = createOptions();

			CommandLine cmd = parse(options, args, MPJ_MeanUCERF3_HazardCurveTest.class);
			
			args = cmd.getArgs();

			Preconditions.checkArgument(args.length == 2, "Must specify inputfile file/output dir!");

			File inputFile = new File(args[0]);
			File dir = new File(args[1]);
			
			if (!dir.exists())
				dir.mkdir();

			Preconditions.checkArgument(inputFile.exists(), "Input file doesn't exist!: "+inputFile);
			
			CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(inputFile);
			
			MPJ_MeanUCERF3_HazardCurveTest driver = new MPJ_MeanUCERF3_HazardCurveTest(cmd, cfss, dir);
			
			driver.run();
			
			finalizeMPJ();
			
			System.exit(0);
		} catch (Throwable t) {
			abortAndExit(t);
		}
	}

}
