package scratch.kevin.tdProbModelPlayground;

import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;

import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.editor.impl.ParameterListEditor;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.AperiodicityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.HistoricalOpenInterval;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.HistoricalOpenIntervals;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.RenewalModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.TimeDepFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.UCERF3_ProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.WG02_ProbabilityModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_WUS_BranchAveragedERF;

import com.google.common.base.Preconditions;

import scratch.UCERF3.erf.mean.MeanUCERF3;

public class TD_ERF_Example {
	
	static FaultSystemSolution fetchNSHM23() throws IOException {
		return NSHM23_WUS_BranchAveragedERF.loadSolution();
	}
	
	static FaultSystemSolution fetchU3_BA() throws IOException {
		File storeDir = MeanUCERF3.getStoreDir();
		File solFile = MeanUCERF3.checkDownload(
				new File(storeDir, "cached_FM3_1_dep100.0_depMean_rakeMean.zip")).join();
		return FaultSystemSolution.load(solFile);
	}

	public static void main(String[] args) throws IOException {
		TimeDepFaultSystemSolutionERF erf = new TimeDepFaultSystemSolutionERF();
		
		// load a solution
		// custom file
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/path/to/sol.zip"));
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/u3_ba_fm3.1_sol.zip"));
		// or use use Akash's fancy data accessor (will download once and cache)
//		FaultSystemSolution sol = fetchNSHM23(); // doesn't have any DOLE data attached by default
		FaultSystemSolution sol = fetchU3_BA();
		
		erf.setSolution(sol);
		
		// set to a value from the enum:
//		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
//		erf.setProbabilityModelChoice(FSS_ProbabilityModels.WG02);
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.UCERF3_METHOD);
		
		FSS_ProbabilityModel probModel = erf.getProbabilityModel();
		
		if (probModel instanceof UCERF3_ProbabilityModel) {
			UCERF3_ProbabilityModel u3ProbModel = (UCERF3_ProbabilityModel)probModel;
			// setting by enum is prefferred
			u3ProbModel.setAperiodicityModelChoice(AperiodicityModels.UCERF3_HIGH);
			// but can still set a custom implementation if needed (the GUI will say "External Custom Value: <name>")
			// this test sets an external model that's also in the enum, but setting it this way makes it "external"
			u3ProbModel.setCustomHistOpenIntervalModel(new HistoricalOpenInterval.SingleYear(1875, true));
		} else if (probModel instanceof WG02_ProbabilityModel) {
			WG02_ProbabilityModel wgProbModel = (WG02_ProbabilityModel)probModel;
			wgProbModel.setAperiodicityModelChoice(AperiodicityModels.UCERF3_LOW);
		}
		
		// can also do things via the parameter list
		ParameterList modelParams = probModel.getAdjustableParameters();
		if (modelParams.containsParameter(RenewalModels.PARAM_NAME))
			modelParams.setValue(RenewalModels.PARAM_NAME, RenewalModels.BPT);
		
		// sanity checks
		erf.updateForecast();
		
		long startTimeMillis = erf.getTimeSpan().getStartTimeInMillis();
		double duration = erf.getTimeSpan().getDuration();
		for (int sourceID=0; sourceID<erf.getNumFaultSystemSources(); sourceID++) {
			ProbEqkSource source = erf.getSource(sourceID);
			int fssIndex = erf.getFltSysRupIndexForSource(sourceID);
			double probGain = probModel.getProbabilityGain(fssIndex, startTimeMillis, duration);
			Preconditions.checkState(source.getNumRuptures() > 0,
					"Source %s is empty! fssIndex=%s, rate=%s, included=%s, probGain=%s",
					sourceID, fssIndex, (double)sol.getRateForRup(fssIndex), erf.isRuptureIncluded(fssIndex), (double)probGain);
		}
		
		showGUI(erf);
	}
	
	private static void showGUI(TimeDepFaultSystemSolutionERF erf) {
		ParameterList paramList = erf.getAdjustableParameterList();
		
		ParameterListEditor editor = new ParameterListEditor(paramList);
		editor.setTitle(erf.getName()+" Parameters");
		
		JFrame frame = new JFrame();
		frame.setContentPane(editor);
		frame.setSize(600, 1000);
		frame.setLocationRelativeTo(null);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}

}
