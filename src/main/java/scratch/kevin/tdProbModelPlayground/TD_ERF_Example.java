package scratch.kevin.tdProbModelPlayground;

import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;

import org.opensha.commons.param.ParameterList;
import org.opensha.commons.param.editor.impl.ParameterListEditor;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.AperiodicityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModel;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.FSS_ProbabilityModels;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.TimeDepFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.erf.td.UCERF3_ProbabilityModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.erf.NSHM23_WUS_BranchAveragedERF;

public class TD_ERF_Example {

	public static void main(String[] args) throws IOException {
		TimeDepFaultSystemSolutionERF erf = new TimeDepFaultSystemSolutionERF();
		
		// load a solution
		// custom file
//		FaultSystemSolution sol = FaultSystemSolution.load(new File("/path/to/sol.zip"));
		// or use use Akash's fancy data accessor (will download once and cache)
		FaultSystemSolution sol = NSHM23_WUS_BranchAveragedERF.loadSolution();
		
		erf.setSolution(sol);
		
		// set to a value from the enum:
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.POISSON);
		erf.setProbabilityModelChoice(FSS_ProbabilityModels.UCERF3_METHOD);
		
		// set one of the underlying values to a custom instance
		((UCERF3_ProbabilityModel)erf.getProbabilityModel()).setCustomAperiodicityModel(
				AperiodicityModels.SINGLE_VALUED.instance(sol));
		
//		// set to an external custom value
//		// doing this will cause the GUI to display (external custom value)
//		erf.setCustomProbabilityModel(new FSS_ProbabilityModel.Poisson(sol));
		
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
