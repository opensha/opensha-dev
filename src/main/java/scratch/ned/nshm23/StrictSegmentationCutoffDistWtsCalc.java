package scratch.ned.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.calc.nnls.NNLSWrapper;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;

import com.google.common.io.Files;

/**
 * This class solves for the weights for hard-cutoff distance inversions that, when summed, will match (in a least-squares sense)
 * a distance decay of A*exp(-r/ro), where A is a constant, r is fault-jump distance and ro is a constant (both in units of km).
 * 
 * @author field
 *
 */
public class StrictSegmentationCutoffDistWtsCalc {
	

	/**
	 * 
	 * @param cutoffDistanceArray - the set of distance cutoffs applied in the fault-system-solution inversions
	 * @param dataList - the cutoff rate data (in same order as above)
	 * @param target_ro - the target decay rate parameter (typically 3 km).
	 * @return
	 */
	public static double[] getWts(int[] cutoffDistanceArray, ArrayList<double[][]> dataList, double target_ro) {
		
		int numCols = cutoffDistanceArray.length+1;  // the plus 1 is for A
		
		double[][] data = dataList.get(0);
//		System.out.println(data.length);
//		System.out.println(data[0].length);
		
		int numRows = data[0].length + 1;	// plus 1 is for -exp(-r/ro) values

		double[] d = new double[numRows];
		d[d.length-1] = 1.0;  // sum of weights must equal 1.0; all other array elements are zero
		

		// matrix is [row][col]	
		double[][] C = new double[numRows][numCols];
		
		for(int col=0;col<dataList.size();col++) {
			data = dataList.get(col);
			for(int row=0;row<data[0].length;row++) {
				C[row][col]=data[1][row];
			}
			// put "1.0 in last row (sum of wts must equal 1.0)
			C[numRows-1][col] = 1.0;
		}

		for(int row=0;row<data[0].length;row++) {
			double r = data[0][row];	// get distance from first column of any of the data objects
			C[row][numCols-1] = -Math.exp(-r/target_ro);
		}
		C[numRows-1][numCols-1] = 0;
		
//		// write out matrices
//		for(int r=0;r<numRows;r++) {
//			System.out.print("\n");
//			for(int c=0;c<numCols;c++) {
//				System.out.print(C[r][c]+"\t");
//			}
//		}
//		System.out.print("\n");
//
//		for(int r=0;r<numRows;r++) {
//			System.out.println(d[r]);
//		}

		
		// NNLS inversion:
		NNLSWrapper nnls = new NNLSWrapper();

		int nRow = C.length;
		int nCol = C[0].length;
		
//		System.out.println("NNLS: nRow="+nRow+"; nCol="+nCol);
		
		double[] A = new double[nRow*nCol];
		double[] x = new double[nCol];
		
		int i,j,k=0;
			
		for(j=0;j<nCol;j++) 
			for(i=0; i<nRow;i++)	{
				A[k]=C[i][j];
				k+=1;
			}
		nnls.update(A,nRow,nCol);
		
		boolean converged = nnls.solve(d,x);
		if(!converged)
			throw new RuntimeException("ERROR:  NNLS Inversion Failed");
		
		return x;
	}
	
	
	private static double[][] readDataFromFile(String fileName) {
		double data[][] = null;
		File file = new File(fileName);
		List<String> fileLines;
		try {
			fileLines = Files.readLines(file, Charset.defaultCharset());
//			System.out.println(fileLines.size());
			data = new double[3][fileLines.size()-1];
			for(int i=1; i<fileLines.size();i++ ) {
				String str = fileLines.get(i);
				String[] split = str.split(",");
				data[0][i-1] = Double.parseDouble(split[0]); // distance
				data[1][i-1] = Double.parseDouble(split[1]); // mean rate
				data[2][i-1] = Double.parseDouble(split[2]); // median rate

//				System.out.println(data[0][i-1]+"\t"+data[1][i-1]+"\t"+data[2][i-1]);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return data;
	}


	
	public static void main(String[] args) {
		
		// assuming functional form of Ae^(-r/ro), we solve for the weights w[] that will 
		// produce the target ro.  We solve for A at the same time (the last element of the w[] array).
		
		double target_ro = 3d;
		
		// both these work:
//		int[] cutoffDistanceArray = {1,3,5,7,9,11,13,15};
		int[] cutoffDistanceArray = {3,5,7,9,11,13,15};
		
		// read data files
		String rootPath = "src/main/java/scratch/ned/nshm23/passthroughData/";
		ArrayList<double[][]> dataList = new ArrayList<double[][]>();
		for(int cutoffDist:cutoffDistanceArray) {
			String fileName = rootPath+"conn_passthrough_shaw07_m7.0_"+cutoffDist+"km.csv";
			System.out.println(fileName);
			double data[][] = readDataFromFile(fileName);
			dataList.add(data);
		}
		
		double[] wtArray = getWts(cutoffDistanceArray, dataList, target_ro);
		
		// last element here is for "A" in A*exp(-r/ro)
		for(double wt:wtArray)
			System.out.println(wt);
		
		
		ArbitrarilyDiscretizedFunc aveRateFunc = new ArbitrarilyDiscretizedFunc();
		aveRateFunc.setName("aveRateFunc");
		ArbitrarilyDiscretizedFunc targetRateFunc = new ArbitrarilyDiscretizedFunc();
		targetRateFunc.setName("targetRateFunc");

		double[][] firstDataArray = dataList.get(0);
		int numPts = firstDataArray[0].length;
		double[] aveRatesArray = new double[numPts];
		for(int i=0;i<dataList.size();i++) {
			double[][] dataArray = dataList.get(i);
			for(int row=0;row<dataArray[0].length;row++) {
				aveRatesArray[row] += wtArray[i]*dataArray[1][row];
			}
		}
		for(int row=0;row<firstDataArray[0].length;row++) {
			double dist = firstDataArray[0][row];
			aveRateFunc.set(dist, aveRatesArray[row]);
			double targetRate = wtArray[wtArray.length-1]*Math.exp(-dist/target_ro);
			targetRateFunc.set(dist,targetRate);
		}
		
		
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(targetRateFunc);
		funcs.add(aveRateFunc);
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		GraphWindow graph = new GraphWindow(funcs, "Paasthrough Rates");
//		graph.setX_AxisRange(xAxisRange.getLowerBound(),xAxisRange.getUpperBound());
//		graph.setY_AxisRange(yAxisRange.getLowerBound(),yAxisRange.getUpperBound());
		graph.setXLog(false);
		graph.setYLog(true);
		graph.setPlotChars(plotChars);
		graph.setX_AxisLabel("Distance");
		graph.setY_AxisLabel("Ave Passthrough Rate");
		graph.setTickLabelFontSize(18);
		graph.setAxisLabelFontSize(20);


	}

}
