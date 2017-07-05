package scratch.martinez.beans;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JToolBar;

import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.data.Range;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.gui.plot.ButtonControlPanel;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.GraphWidget;
import org.opensha.commons.gui.plot.PlotSpec;


/**
 * <p>Title: GraphPane</p>
 *
 * <p>Description: This class allows user to visualise the computed data as graphs.</p>
 * @author Ned Field,Nitin Gupta,E.V.Leyendecker, Eric Martinez
 */
public class NSHMPGraphPane extends JPanel {
	private static final long serialVersionUID	= 1;
	
	JMenuBar menuBar = new JMenuBar();
	JMenu fileMenu = new JMenu();

	JMenuItem fileExitMenu = new JMenuItem();
	JMenuItem fileSaveMenu = new JMenuItem();
	JMenuItem filePrintMenu = new JCheckBoxMenuItem();
	JToolBar jToolBar = new JToolBar();

	private final static int W = 512;
	private final static int H = 730;
	private JSplitPane chartSplitPane = new JSplitPane();
	private JPanel chartPane = new JPanel();
	private GridBagLayout gridBagLayout1 = new GridBagLayout();
	private BorderLayout borderLayout1 = new BorderLayout();
	private FlowLayout flowLayout1 = new FlowLayout();

	private String plotTitle = "Hazard Curves";

	//instance of the GraphPanel class
	private GraphWidget graphWidget;

	private JComboBox graphListCombo = new JComboBox();

	/**
	 * Creating this Map to keep track of the selected item to plot
	 */
	private TreeMap<String, ArrayList<ArbitrarilyDiscretizedFunc>> map = new
				TreeMap<String, ArrayList<ArbitrarilyDiscretizedFunc>>();

	/**
	 * Class constructor that shows the list of graphs that user can plot.
	 * @param dataList ArrayList List of DiscretizedFunctionList
	 */
	public NSHMPGraphPane(ArrayList dataList) {
		//adding list of graphs to the shown to the user.
		int size = 0;
		if(dataList != null) size = dataList.size();

		//creating the ArrayList for the plots
		for (int i = 0; i < size; ++i) {
			//getting the functions to plot and creating individual ArrayList for those
			//adding these individual arraylist to the hashmap.
			ArbitrarilyDiscretizedFunc function = (ArbitrarilyDiscretizedFunc)
					dataList.get(i);
			ArrayList<ArbitrarilyDiscretizedFunc> plottingFunction = new ArrayList<ArbitrarilyDiscretizedFunc>();
			plottingFunction.add(function);
			map.put(function.getName(), plottingFunction);
		}

		//adding the functions having same X and Y axis name to HashMap
		for (int i = 0; i < size; ++i) {
			ArbitrarilyDiscretizedFunc function = (ArbitrarilyDiscretizedFunc)
					dataList.get(i);
			ArrayList<ArbitrarilyDiscretizedFunc> plottingFunctions = new ArrayList<ArbitrarilyDiscretizedFunc>();
			String functionXName = function.getXAxisName();
			String functionYName = function.getYAxisName();
			boolean containsSameName = false;
			String name = null;
			for (int j = i + 1; j < size; ++j) {
				ArbitrarilyDiscretizedFunc function1 = (ArbitrarilyDiscretizedFunc)
						dataList.get(j);
				String function1XName = function1.getXAxisName();
				String function1YName = function1.getYAxisName();
				if (functionXName.equals(function1XName) &&
						functionYName.equals(function1YName)) {
					//name = function1YName + " Vs " + function1XName;
					if (functionXName.equals("Damage Factor"))
						name = "Loss Curve Summary";
					else
						name = "Basic Hazard Curve";
					if (!map.containsKey(name)) {
						plottingFunctions.add(function1);
						containsSameName = true;
					}
				}
			}
			if (containsSameName) {
				plottingFunctions.add(function);
				map.put(name, plottingFunctions);
				containsSameName = false;
			}
		}
	 

//Adding the names of the plot to the Combo selection
		Set plotNames = map.keySet();
		Iterator it = plotNames.iterator();
		while (it.hasNext()) {
			graphListCombo.addItem(it.next());
		}

		graphListCombo.setSelectedIndex(0);
		if(graphListCombo.getItemCount() > 1) 
			graphListCombo.setSelectedIndex(1);

		try {
			jbInit();
		}
		catch (Exception exception) {
			exception.printStackTrace();
		}

		graphWidget = new GraphWidget();
		graphWidget.getGraphPanel().setDividerLocation(350);
		chartPane.add(graphWidget, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
				, GridBagConstraints.CENTER, GridBagConstraints.BOTH,
				new Insets(2, 2, 2, 2), 0, 0));
		chartPane.validate();
		chartPane.repaint();
		drawGraph();
	}

	/**
	 * Component initialization.
	 *
	 * @throws java.lang.Exception
	 */
	private void jbInit() throws Exception {
		setLayout(borderLayout1);
		setSize(new Dimension(W, H));
		
		graphListCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent itemEvent) {
				graphListCombo_itemStateChanged(itemEvent);
			}
		});
		
		chartSplitPane.setOrientation(JSplitPane.VERTICAL_SPLIT);
		chartPane.setLayout(gridBagLayout1);

		add(chartSplitPane, BorderLayout.CENTER);
		chartSplitPane.add(chartPane, JSplitPane.TOP);
		chartSplitPane.setDividerLocation(450);
		chartPane.add(graphListCombo, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
				, GridBagConstraints.CENTER, GridBagConstraints.BOTH,
				new Insets(4, 4, 4, 4), 0, 0));
	}

	/**
	 * to draw the graph
	 */
	private void drawGraph() {

		//getting the list of the curves that we need to plot
		String selectedDataToPlot = (String) graphListCombo.getSelectedItem();

		//show correct graph titles
		plotTitle = selectedDataToPlot;

		ArrayList<ArbitrarilyDiscretizedFunc> functionsToPlot = map.get(selectedDataToPlot);
		ArbitrarilyDiscretizedFunc func = functionsToPlot.get(0);
		PlotSpec spec = graphWidget.getPlotSpec();
		spec.setXAxisLabel(func.getXAxisName());
		spec.setXAxisLabel(func.getYAxisName());
		spec.setPlotElems(functionsToPlot);
		spec.setTitle(plotTitle);
		graphWidget.drawGraph();
		graphWidget.updateUI();
	}

	/**
	 * plots the curves with defined color,line width and shape.
	 */
	public void plotGraphUsingPlotPreferences() {
		drawGraph();
	}

	/**
	 *
	 * sets plot Title
	 */
	public void setPlotLabel(String plotTitle) {
		this.plotTitle = plotTitle;
	}

	public void graphListCombo_itemStateChanged(ItemEvent itemEvent) {
		graphWidget.removeChartAndMetadata();
		drawGraph();
	}
	
	public void setLogSpace(boolean xlog, boolean ylog) {
		graphWidget.setX_Log(xlog);
		graphWidget.setY_Log(ylog);
	}
}
