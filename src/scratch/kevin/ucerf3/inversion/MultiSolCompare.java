package scratch.kevin.ucerf3.inversion;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.GraphWindow;

import scratch.UCERF3.utils.MatrixIO;

public class MultiSolCompare {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2012_02_27-unconstrained");
		
		ArrayList<DiscretizedFunc> funcs = new ArrayList<DiscretizedFunc>();
		
		int num = 0;
		double[] average = null;
		
		for (File file : dir.listFiles()) {
			if (!file.getName().endsWith(".bin"))
				continue;
			if (!file.getName().contains("LIMIT"))
				continue;
			try {
				double[] solution = MatrixIO.doubleArrayFromFile(file);
				
				EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0, solution.length, 1d);
				num++;
				
				if (average == null)
					average = new double[solution.length];
				
				for (int i=0; i<solution.length; i++) {
					func.set(i, solution[i]);
					average[i] = average[i] + solution[i];
				}
				
				funcs.add(func);
			} catch (IOException e) {
				e.printStackTrace();
				continue;
			}
		}
		
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0, average.length, 1d);
		for (int i=0; i<average.length; i++) {
			average[i] = average[i] / (double)num;
			func.set(i, average[i]);
		}
		funcs.add(0, func);
		
		ArrayList<PlotCurveCharacterstics> chars = new ArrayList<PlotCurveCharacterstics>();
		
		List<Color> colors = GraphWindow.generateDefaultColors();
		for (int i=0; i<funcs.size(); i++) {
			float width;
			if (i == 0)
				width = 3f;
			else
				width = 1f;
			PlotCurveCharacterstics pchar = new PlotCurveCharacterstics(PlotLineType.SOLID, width, colors.get(i));
			chars.add(pchar);
		}
		
		MatrixIO.doubleArrayToFile(average, new File(dir, "limit_average.bin"));
		
//		new GraphWindow(funcs, "Rates vs ID", chars);
	}

}
