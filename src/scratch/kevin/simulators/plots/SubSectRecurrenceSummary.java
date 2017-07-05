package scratch.kevin.simulators.plots;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.simulators.EventRecord;
import org.opensha.sha.simulators.SimulatorElement;
import org.opensha.sha.simulators.SimulatorEvent;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class SubSectRecurrenceSummary extends AbstractPlot {
	
	private double[] minMags;
	private int minSectID;
	private Map<Integer, String> sectIDMap;
	private List<double[]> prevTimes;
	private List<List<MinMaxAveTracker>> trackers;
	
	private List<EvenlyDiscretizedFunc> compareCumulativeMFDs;
	private String compareName;
	
	public SubSectRecurrenceSummary(List<SimulatorElement> elems, double... minMags) {
		Preconditions.checkState(minMags.length > 0);
		this.minMags = minMags;
		sectIDMap = Maps.newHashMap();
		minSectID = Integer.MAX_VALUE;
		for (SimulatorElement elem : elems) {
			Preconditions.checkState(elem.getSectionID() >= 0, "Must have section IDs");
			if (elem.getSectionID() < minSectID)
				minSectID = elem.getSectionID();
			sectIDMap.put(elem.getSectionID(), elem.getSectionName());
		}
		trackers = Lists.newArrayList();
		prevTimes = Lists.newArrayList();
		for (int m=0; m<minMags.length; m++) {
			List<MinMaxAveTracker> myTrackers = Lists.newArrayList();
			double[] myPrevTimes = new double[sectIDMap.size()];
			for (int i=0; i<sectIDMap.size(); i++) {
				myTrackers.add(new MinMaxAveTracker());
				myPrevTimes[i] = Double.NaN;
			}
			trackers.add(myTrackers);
			prevTimes.add(myPrevTimes);
		}
	}
	
	public void setComparison(List<EvenlyDiscretizedFunc> compareCumulativeMFDs, String compareName) {
		this.compareCumulativeMFDs = compareCumulativeMFDs;
		this.compareName = compareName;
	}

	@Override
	protected void doProcessEvent(SimulatorEvent e) {
		HashSet<Integer> sectIDs = new HashSet<Integer>();
		for (EventRecord rec : e)
			for (SimulatorElement elem : rec.getElements())
				sectIDs.add(elem.getSectionID());
		
		double time = e.getTimeInYears();
		double mag = e.getMagnitude();
		for (int m=0; m<minMags.length; m++) {
			if (mag < minMags[m])
				continue;
			double[] myPrevTimes = prevTimes.get(m);
			List<MinMaxAveTracker> myTrackers = trackers.get(m);
			for (Integer sectID : sectIDs) {
				int index = sectID - minSectID;
				double prevTime = myPrevTimes[index];
				if (!Double.isNaN(prevTime)) {
					double interval = time - prevTime;
					myTrackers.get(index).addValue(interval);
				}
				myPrevTimes[index] = time;
			}
		}
	}

	@Override
	protected void finalize() throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList("Index", "Name");
		for (double minMag : minMags) {
			String magStr = getCleanMagLabel(minMag);
			header.add(magStr+" "+getCatalogName());
			if (compareCumulativeMFDs != null)
				header.add(magStr+" "+compareName);
		}
		csv.addLine(header);
		
		DefaultXY_DataSet[] scatters = null;
		if (compareCumulativeMFDs != null) {
			scatters = new DefaultXY_DataSet[minMags.length];
			for (int m=0; m<minMags.length; m++)
				scatters[m] = new DefaultXY_DataSet();
		}
		
		for (int i=0; i<trackers.get(0).size(); i++) {
			String name = sectIDMap.get(i+minSectID);
			List<String> line = Lists.newArrayList(i+"", name);
			for (int m=0; m<minMags.length; m++) {
				double mean = trackers.get(m).get(i).getAverage();
				line.add(mean+"");
				if (compareCumulativeMFDs != null) {
					double compVal = 1d/compareCumulativeMFDs.get(i).getInterpolatedY_inLogYDomain(minMags[m]);
					line.add(compVal+"");
					if (Doubles.isFinite(mean) && Doubles.isFinite(compVal))
						scatters[m].set(mean, compVal);
				}
			}
			
			csv.addLine(line);
		}
		
		File outputDir = getOutputDir();
		String prefix = getOutputPrefix();
		
		csv.writeToFile(new File(outputDir, prefix+".csv"));
		if (scatters != null) {
			for (int m=0; m<scatters.length; m++) {
				String magLabel = getCleanMagLabel(minMags[m]);
				String myPrefix = prefix+"_m"+magLabel;
				DefaultXY_DataSet scatter = scatters[m];
				
				List<XY_DataSet> funcs = Lists.newArrayList();
				List<PlotCurveCharacterstics> chars = Lists.newArrayList();
				
				double minVal = minNonZero(scatter, true);
				minVal = Math.min(minVal, minNonZero(scatter, false));
				double maxVal = Math.max(scatter.getMaxX(), scatter.getMaxY());
				Range range = calcEncompassingLog10Range(minVal, maxVal);
				minVal = range.getLowerBound();
				maxVal = range.getUpperBound();
				
				funcs.add(scatter);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
				
				funcs.add(getLine(null, minVal, minVal, maxVal, maxVal));
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
				
				PlotSpec plot = new PlotSpec(funcs, chars, "M≥"+magLabel+" Subsection Interevent Time Scatter",
						getCatalogName()+" (years)", compareName+" (years)");
				plot.setLegendVisible(true);
				
				HeadlessGraphPanel gp = getGraphPanel();
				gp.drawGraphPanel(plot, true, true, range, range);
				gp.getChartPanel().setSize(1000, 800);
				gp.saveAsPNG(new File(getOutputDir(), myPrefix+".png").getAbsolutePath());
				gp.saveAsPDF(new File(getOutputDir(), myPrefix+".pdf").getAbsolutePath());
			}
		}
	}

	@Override
	public Collection<SimulatorElement> getApplicableElements() {
		return null;
	}

}
