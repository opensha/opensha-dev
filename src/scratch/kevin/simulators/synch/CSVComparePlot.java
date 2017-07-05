package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.swing.JFrame;

import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

public class CSVComparePlot {

	public static void main(String[] args) throws IOException {
		File mainDir = new File("/home/kevin/Simulators/synch");
		
		String name1 = "Indep";
		File csvFile1 = new File(mainDir, "weight_FREQ_PROD_indep/synch_params_std_devs.csv");
		Color color1 = Color.RED;
		
		String name2 = "Orig";
		File csvFile2 = new File(mainDir, "weight_FREQ_PROD_dep/synch_params_std_devs.csv");
		Color color2 = Color.BLUE;
		
		Table<String, String, double[]> table1 = loadTable(csvFile1);
		Table<String, String, double[]> table2 = loadTable(csvFile2);
		
		List<XY_DataSet> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		double curX = 0.25;
		float lineWidth = 0.2f;
		
		double annY = 1.9;
		
		List<XYAnnotation> anns = Lists.newArrayList();
		
		for (Cell<String, String, double[]> cell : table1.cellSet()) {
			String fName1 = cell.getRowKey();
			String fName2 = cell.getColumnKey();
			
			funcs.add(createLine(curX - 0.25, -2, curX - 0.25, 2));
			chars.add(new PlotCurveCharacterstics(PlotLineType.DASHED, 1f, Color.GRAY));
			
			double[] vals1 = table1.get(fName1, fName2);
			double[] vals2 = table2.get(fName1, fName2);
			
			addFuncs(funcs, chars, vals1, lineWidth, color1, curX);
			curX += 0.5;
			
			addFuncs(funcs, chars, vals2, lineWidth, color2, curX);
			curX += 0.5;
			
			double annX1 = curX - 0.75 - 0.9*lineWidth;
			double annX2 = curX - 0.75 + 0.9*lineWidth;
//			System.out.println("Ann x's: "+annX1+" "+annX2);
			anns.add(buildAnn(annX1, annY, fName1));
			anns.add(buildAnn(annX2, annY, fName2));
		}
		
		anns.add(buildAnn(1d, -1.6, name1, color1, 18));
		anns.add(buildAnn(3d, -1.6, name2, color2, 18));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Compare "+name1+" to "+name2, "", "Ln(Gbar)");
		spec.setPlotAnnotations(anns);
		
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gw.setAxisRange(0, curX-0.25, -2, 2);
	}
	
	private static void addFuncs(List<XY_DataSet> funcs, List<PlotCurveCharacterstics> chars,
			double[] vals, float lineWidth, Color color, double curX) {
		funcs.add(createBar(curX, vals[0]));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID_BAR, lineWidth, color));
		funcs.add(createLine(curX - 0.5*lineWidth, vals[1], curX + 0.5*lineWidth, vals[1]));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		funcs.add(createLine(curX - 0.5*lineWidth, -vals[1], curX + 0.5*lineWidth, -vals[1]));
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
	}
	
	private static XY_DataSet createBar(double x, double y) {
		if (y > 0)
			return createLine(x, 0, x, y);
		else
			return createLine(x, y, x, 0);
	}
	
	private static XY_DataSet createLine(double x1, double y1, double x2, double y2) {
		XY_DataSet xy = new DefaultXY_DataSet();
		
		xy.set(x1, y1);
		xy.set(x2, y2);
		
		return xy;
	}
	
	private static XYTextAnnotation buildAnn(double x, double y, String text) {
		return buildAnn(x, y, text, Color.BLACK, 12);
	}
	
	private static XYTextAnnotation buildAnn(double x, double y, String text, Color color, int fontSize) {
		XYTextAnnotation ann = new XYTextAnnotation(text, x, y);
		ann.setTextAnchor(TextAnchor.CENTER_LEFT);
		ann.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, fontSize));
		ann.setPaint(color);
		if (fontSize == 12) {
			ann.setRotationAnchor(TextAnchor.CENTER_LEFT);
			ann.setRotationAngle(Math.PI*0.5);
		}
		return ann;
	}
	
	private static Table<String, String, double[]> loadTable(File file) throws IOException {
		CSVFile<String> csv = CSVFile.readFile(file, false);
		int rowStart = 22;
		int colStart = 1;
		int num = 6;
		
		Table<String, String, double[]> table = HashBasedTable.create();
		
		for (int i=0; i<num; i++) {
			String name1 = csv.get(rowStart-1, colStart + i);
			for (int j=i+1; j<num; j++) {
				String name2 = csv.get(rowStart-1, colStart + j);
				String[] valStrs = csv.get(rowStart+i, colStart+j).split(" ");
				double val = Double.parseDouble(valStrs[0]);
				double stdDev = Double.parseDouble(valStrs[2]);
				
				double[] ret = { val, stdDev };
				table.put(name1, name2, ret);
			}
		}
		
		return table;
	}

}
