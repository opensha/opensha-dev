package scratch.kevin;

import java.util.ArrayList;
import java.util.Random;

import javax.swing.JFrame;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.GraphWindow;

public class ScatterGraphTest {
	
	private GraphPanel gp;
	private Random r = new Random();
	
	public ScatterGraphTest() {
		GraphWindow gw = new GraphWindow(getFuncList(), "Title");

		gw.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gw.setVisible(true);
	}
	
	private ArrayList<XY_DataSet> getFuncList() {
		ArrayList<XY_DataSet> funcs = new ArrayList<XY_DataSet>();
		
		ArbitrarilyDiscretizedFunc func = new ArbitrarilyDiscretizedFunc();
		DefaultXY_DataSet scatter = new DefaultXY_DataSet();
		for (double x=0; x<1d; x++) {
			double y = x * 0.8;
			y += r.nextDouble() - 0.5d;
			System.out.println("Adding " + x + ", " + y);
			func.set(x, y);
			
			scatter.set(x, y + r.nextDouble()*0.5);
			scatter.set(x, y - r.nextDouble()*0.5);
		}
		
		funcs.add(func);
		funcs.add(scatter);
		
		return funcs;
	}
	
	public static void main(String[] args) {
		new ScatterGraphTest();
	}

}
