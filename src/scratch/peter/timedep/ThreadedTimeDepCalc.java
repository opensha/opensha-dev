package scratch.peter.timedep;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationList;
import org.opensha.nshmp2.calc.HazardResult;
import org.opensha.nshmp2.calc.HazardResultWriter;
import org.opensha.nshmp2.calc.HazardResultWriterMPJ;
import org.opensha.nshmp2.tmp.TestGrid;
import org.opensha.nshmp2.util.Period;
import org.opensha.nshmp2.util.SourceIMR;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.param.AleatoryMagAreaStdDevParam;
import org.opensha.sha.earthquake.param.ApplyGardnerKnopoffAftershockFilterParam;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;

import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.peter.nshmp.NSHMP_UCERF2_ERF;
import scratch.peter.nshmp.NSHMP_UCERF3_ERF;
import scratch.peter.ucerf3.calc.UC3_CalcUtils;

/**
 * Class manages multithreaded NSHMP hazard calculations. Farms out
 * {@code HazardCalc}s to locally available cores and pipes results to a
 * supplied {@code Queue}.
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class ThreadedTimeDepCalc {

	private LocationList locs;
	private Period period;
	private HazardResultWriter writer;
	private ERF erf;
	private SourceIMR imr = null;

	private static final String UC3_SOL_PATH = 
//			"/home/scec-00/pmpowers/UC33/src/bravg/FM/UC33brAvg_FM31.zip";
//			"tmp/UC33/src/bravg/FM/UC33brAvg_FM31.zip";
			"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
			+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip";

	/**
	 * Initializes a new threaded hazard calculation for the specified UC3 
	 * (averaged) solution.
	 */
	@SuppressWarnings("javadoc")
	public ThreadedTimeDepCalc(String erfStr, SourceIMR imr, LocationList locs,
		Period period, IncludeBackgroundOption bg, boolean timeDep,
		double duration, HazardResultWriter writer) {
		this.locs = locs;
		this.period = period;
		this.writer = writer;
		this.imr = imr;
		System.out.println("Starting threaded hazard calc...");
		erf = getERF(erfStr, bg, timeDep, duration);
		erf.updateForecast();
	}

	
	private ERF getERF(String erfStr, IncludeBackgroundOption bg,
			boolean timeDep, double duration) {
		
		if (erfStr.equals("UC3")) {
			
			FaultSystemSolution fss = UC3_CalcUtils.getSolution(UC3_SOL_PATH);
			FaultSystemSolutionERF erf = new NSHMP_UCERF3_ERF(fss);
			// no aleatory mag-area uncertainty
			// no aftershocks
			erf.setParameter(IncludeBackgroundParam.NAME, bg);
			erf.setParameter(AleatoryMagAreaStdDevParam.NAME, 0.0);
			erf.setParameter(ApplyGardnerKnopoffAftershockFilterParam.NAME, true);
			
			if (timeDep) {
				erf.setName("TimeDepUC3");
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_PREF_BLEND);
				erf.setParameter(HistoricOpenIntervalParam.NAME, (double)(FaultSystemSolutionERF.START_TIME_DEFAULT - 1875));
				erf.setParameter(BPTAveragingTypeParam.NAME, BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
				erf.getTimeSpan().setStartTime(2014);
				erf.getTimeSpan().setDuration(duration);
			} else {
				erf.setName("TimeIndepUC3");
				erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
				erf.getTimeSpan().setDuration(duration);
			}
			return erf;
		}
		
		if (erfStr.equals("UC2")) {
			
			NSHMP_UCERF2_ERF erf = new NSHMP_UCERF2_ERF();
			
			erf.setParameter(MeanUCERF2.RUP_OFFSET_PARAM_NAME, 1.0);
			erf.setParameter(UCERF2.BACK_SEIS_NAME, UCERF2.BACK_SEIS_INCLUDE);
			erf.setParameter(UCERF2.BACK_SEIS_RUP_NAME, UCERF2.BACK_SEIS_RUP_POINT);
			erf.setParameter(UCERF2.FLOATER_TYPE_PARAM_NAME, UCERF2.FULL_DDW_FLOATER);
			
			if (timeDep) {
				erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, MeanUCERF2.PROB_MODEL_WGCEP_PREF_BLEND);
				erf.getTimeSpan().setStartTime(2007);
				erf.getTimeSpan().setDuration(duration);
			} else {
				erf.setParameter(UCERF2.PROB_MODEL_PARAM_NAME, UCERF2.PROB_MODEL_POISSON);
				erf.getTimeSpan().setDuration(duration);
			}
			return erf;
		}
		
		throw new IllegalArgumentException("Invalid erf identifier: " + erfStr);
	}
	
	/**
	 * Calculates hazard curves for the specified indices. Presently no index
	 * checking is performed. If {@code indices} is {@code null}, curves are
	 * calculated at all locations.
	 * 
	 * @param indices of locations to calculate curves for
	 * @throws InterruptedException
	 * @throws ExecutionException
	 * @throws IOException
	 */
	public void calculate(int[] indices) throws InterruptedException,
			ExecutionException, IOException {

		// set up to process all
		if (indices == null) indices = makeIndices(locs.size());
		
		// init thread mgr
		ExecutorService ex = Executors.newFixedThreadPool(
			Runtime.getRuntime().availableProcessors());
		CompletionService<HazardResult> ecs = 
				new ExecutorCompletionService<HazardResult>(ex);

		for (int index : indices) {
			Site site = new Site(locs.get(index));
			TimeDepCalc hc = TimeDepCalc.create(erf, imr, site, period);
			ecs.submit(hc);
		}
		ex.shutdown();

		// process results as they come in; ecs.take() blocks until result
		// is available
		for (int i = 0; i < indices.length; i++) {
			writer.write(ecs.take().get());
//			if (i % 10 == 0) System.out.println("Jobs completed: " + i);
		}
		
//		writer.close(); // not needed for sites writer
	}
	
	private int[] makeIndices(int size) {
		int[] indices = new int[size];
		for (int i=0; i<size; i++) {
			indices[i] = i;
		}
		return indices;
	}
	
	public static void main(String[] args) throws IOException {
		ThreadedTimeDepCalc ttdc = new ThreadedTimeDepCalc(
			"UC3",
			SourceIMR.WUS_FAULT,
			TestGrid.CA_RELM.grid(0.1).getNodeList(),
			Period.GM0P00,
			IncludeBackgroundOption.INCLUDE,
			true,
			50.0,
			new HazardResultWriterMPJ(new File("tmp")));
		
//		System.out.println(ttdc.erf.getAdjustableParameterList());
//		System.out.println(ttdc.erf.getTimeSpan().getStartTimeYear());
//		System.out.println(ttdc.erf.getTimeSpan().getDuration());
		
		
		// test calc
		
		ERF erf = ttdc.erf;
		Site site = new Site(new Location(34.0, -118.3));
		TimeDepCalc calc = TimeDepCalc.create(erf, SourceIMR.WUS_FAULT_14, site, Period.GM1P00);
		HazardResult result = calc.call();
		
		System.out.println(result.curve());
		
		
	}

}
