package scratch.kevin.simulators.catBuild;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.jfree.data.Range;
import org.jfree.ui.RectangleEdge;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.IndicesKey;
import scratch.kevin.markov.PossibleStates;
import scratch.kevin.markov.SparseNDimensionalHashDataset;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.dists.RandomDistType;
import scratch.kevin.simulators.dists.RandomReturnPeriodProvider;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;

public class StateBasedCatalogBuilder implements CatalogBuilder {
	
	private EmpiricalMarkovChain chain;
	
	private SparseNDimensionalHashDataset<Double> totalStatesDataset;
	private SparseNDimensionalHashDataset<PossibleStates> stateTransitionDataset;
	private double distSpacing;

	@Override
	public List<SimulatorEvent> buildCatalog(List<? extends SimulatorEvent> events,
			List<RandomReturnPeriodProvider> randomRPsList,
			List<List<? extends SimulatorEvent>> matchesLists, boolean trim) {
		
		distSpacing = 10d;
		
		// build Markov chain
		chain = MarkovChainBuilder.build(distSpacing, matchesLists);
		int nDims = chain.getNDims();
		totalStatesDataset = chain.getTotalStatesDataset();
		stateTransitionDataset = chain.getStateTransitionDataset();
		
		// this is the randomized sequence of events for each fault from which to sample
		List<List<SimulatorEvent>> eventsToReuse = Lists.newArrayList();
		int[] eventsToReuseIndexes = new int[matchesLists.size()];
		for (int i=0; i<matchesLists.size(); i++) {
			List<? extends SimulatorEvent> matches = matchesLists.get(i);
			List<SimulatorEvent> rand = Lists.newArrayList(matches);
			Collections.shuffle(rand);
			eventsToReuse.add(rand);
			eventsToReuseIndexes[i] = 0;
		}
		
//		double minVal = 0.5d*distSpacing;
		
		// now pick random initial state from the distribution of total states
		int[] prevState = chain.getOccupancy().drawState();
		
		double maxTime = events.get(events.size()-1).getTimeInYears();
		double startTime = events.get(0).getTimeInYears();
		int numSteps = (int)((maxTime - startTime)/distSpacing);
		
		List<SimulatorEvent> randomizedEvents = Lists.newArrayList();
		
		System.out.println("Assembling random catalog");
		
		int[] counts = new int[nDims];
		
		int eventID = 0;
		int numBailouts = 0;
		
		List<int[]> statesTracker = Lists.newArrayList();
		statesTracker.add(prevState);
		
		// shift so that rupTime falls in the middle of windows
		startTime += distSpacing*0.5d;
		
		for (int step=0; step<numSteps; step++) {
			// choose current state randomly from previous state's transition states
			PossibleStates possibilities = stateTransitionDataset.get(prevState);
			
			if (possibilities == null) {
				// this means we found the last state in the system, and it was only reached
				// in that last state
				System.out.println("Reached orig last state in system, no transitions! (step="+step+")");
				// find all theoretical neighbors
				numBailouts++;
				possibilities = StateBasedCatalogBuilder.findPossibleBailoutStates(prevState, totalStatesDataset);
			}
			
			int[] curState = null;
			if (possibilities.getTot() == 0) {
				System.out.println("Possibilities are empty, backing up!");
				// go back up and get out of this path
				
				for (int redoStep=step-1; redoStep>=0; redoStep--) {
					int[] redoPrevIndices = statesTracker.get(redoStep);
					possibilities = stateTransitionDataset.get(redoPrevIndices);
					if (possibilities.getNumStates() > 1) {
						// this means there was another option
						int[] newDestState = possibilities.drawState();
						while (!Arrays.equals(newDestState, statesTracker.get(redoStep+1)))
							newDestState = possibilities.drawState();
						// now make sure
						
						System.out.println("Backed up to step "+redoStep+" ("+(step-redoStep)+" steps). New State: "
								+Joiner.on(",").join(Ints.asList(newDestState)));
						
						boolean newStateStuck = false;
						// now make sure this state doesn't lead to it as well!
						int[] prevTestState = newDestState;
						for (int i=0; i<2*(step-redoStep); i++) {
							PossibleStates testPossibilities = stateTransitionDataset.get(prevTestState);
							if (testPossibilities == null || testPossibilities.getTot() == 0) {
								newStateStuck = true;
								break;
							}
							prevTestState = testPossibilities.drawState();
						}
						
						if (newStateStuck) {
							System.out.println("Nevermind, this one can get stuck as well. Backing up more.");
						} else {
							prevState = redoPrevIndices;
							curState = newDestState;
							
							// remove any added events
							double windowStart = startTime + distSpacing*(redoStep);
							for (int i=randomizedEvents.size(); --i>=0;) {
								SimulatorEvent e = randomizedEvents.get(i);
								if (e.getTimeInYears() > windowStart)
									randomizedEvents.remove(i);
								else
									break;
							}
							step = redoStep;
							break;
						}
					}
				}
				Preconditions.checkNotNull(curState);
			} else {
				curState = possibilities.drawState();
			}
			
			double rupTimeYears = startTime + distSpacing*step;
			double rupTimeSecs = rupTimeYears * General_EQSIM_Tools.SECONDS_PER_YEAR;
			
			// now look for any hits
			for (int n=0; n<curState.length; n++) {
				if (curState[n] == 0) {
					// state was reset to zero for this fault, this means a rupture happens in this window
					List<SimulatorEvent> myEvents = eventsToReuse.get(n);
					if (eventsToReuseIndexes[n] == myEvents.size())
						eventsToReuseIndexes[n] = 0;
					
					SimulatorEvent e = myEvents.get(eventsToReuseIndexes[n]++);
					SimulatorEvent newE = e.cloneNewTime(rupTimeSecs, eventID++);
					randomizedEvents.add(newE);
					
					counts[n]++;
				}
			}
			
			statesTracker.add(curState);
			prevState = curState;
		}
		
		System.out.println("Needed "+numBailouts+" bailouts.");
		for (int n=0; n<nDims; n++)
			System.out.println("iden "+n+": rand="+counts[n]+"\torig="+matchesLists.get(n).size());
		
		return randomizedEvents;
	}

	static void findPossibleStates(int[] curState, int index,
			SparseNDimensionalHashDataset<Double> totalStatesDataset, List<int[]> runningPossibleStates) {
		int[] stateNoRup = Arrays.copyOf(curState, curState.length);
		int[] stateWithRup = Arrays.copyOf(curState, curState.length);
		
		stateNoRup[index]++;
		stateWithRup[index] = 0;
		
		if (index == curState.length-1) {
			runningPossibleStates.add(stateWithRup);
			runningPossibleStates.add(stateNoRup);
		} else {
			findPossibleStates(stateNoRup, index+1, totalStatesDataset, runningPossibleStates);
			findPossibleStates(stateWithRup, index+1, totalStatesDataset, runningPossibleStates);
		}
	}

	static PossibleStates findPossibleBailoutStates(int[] fromState,
			SparseNDimensionalHashDataset<Double> totalStatesDataset) {
		List<int[]> possibilities = Lists.newArrayList();
		
		StateBasedCatalogBuilder.findPossibleStates(fromState, 0, totalStatesDataset, possibilities);
		
		PossibleStates states = new PossibleStates(fromState);
		
		System.out.println("Stuck at a dead end! FromState="+Joiner.on(",").join(Ints.asList(fromState)));
		
		for (int[] state : possibilities) {
			Double val = totalStatesDataset.get(state);
//			System.out.println("Possible state:\t"+Joiner.on(",").join(Ints.asList(state))+"\tval="+val);
			if (val != null)
				states.add(state, val);
		}
		
		System.out.println("Found "+states.getNumStates()+"/"+possibilities.size()+" bailout states");
		
		return states;
	}
	
	private List<File> write2DDists(File writeDir, int index1, String name1, List<? extends SimulatorEvent> matches1,
			int index2, String name2, List<? extends SimulatorEvent> matches2) throws IOException {
		String probFName = "prob_dists_"+PeriodicityPlotter.getFileSafeString(name1)+"_"+PeriodicityPlotter.getFileSafeString(name2);
		File probFile = new File(writeDir, probFName+".pdf");
		String synchFName = "synch_dists_"+PeriodicityPlotter.getFileSafeString(name1)+"_"+PeriodicityPlotter.getFileSafeString(name2);
		File synchFile = new File(writeDir, synchFName+".pdf");
		
		// now we need distributions just between the two
		double maxTimeDiff = 1000d;
		double distDeltaYears = 10d;
		int num = (int)(maxTimeDiff/distDeltaYears - 1);
		double min = 0.5*distDeltaYears;
		
		EvenlyDiscrXYZ_DataSet p_1_n2 = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet p_n1_2 = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet p_1_2 = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet p_n1_n2 = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet synchData = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet synchNormData = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet synchOverIndepData1 = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		EvenlyDiscrXYZ_DataSet synchOverIndepData2 = new EvenlyDiscrXYZ_DataSet(num, num, min, min, distDeltaYears);
		
		for (int x=0; x<num; x++) {
			for (int y=0; y<num; y++) {
				p_1_n2.set(x, y, Double.NaN);
				p_n1_2.set(x, y, Double.NaN);
				p_1_2.set(x, y, Double.NaN);
				p_n1_n2.set(x, y, Double.NaN);
			}
		}
		
		List<List<? extends SimulatorEvent>> matchesLists = Lists.newArrayList();
		matchesLists.add(matches1);
		matchesLists.add(matches2);
		
		double[] freq1s = new double[num];
		double[] freq2s = new double[num];
		double[] tot1s = new double[num];
		double[] tot2s = new double[num];
		
		
		
		double totStateCount = 0;
		for (int[] indices : totalStatesDataset.getPopulatedIndices()) {
			PossibleStates possible = stateTransitionDataset.get(indices);
			double freq1 = 0;
			double freq2 = 0;
			for (int[] state : possible.getStates()) {
				double freq = possible.getFrequency(state);
				if (state[0] == 0)
					freq1 += freq;
				if (state[1] == 0)
					freq2 += freq;
			}
			if (indices[0] < num) {
				freq1s[indices[0]] += freq1;
				tot1s[indices[0]] += possible.getTot();
			}
			if (indices[1] < num) {
				freq2s[indices[1]] += freq2;
				tot2s[indices[1]] += possible.getTot();
			}
			totStateCount += totalStatesDataset.get(indices);
		}
		
		// probabilities in each state independent of other fault
		double[] indepProb1s = new double[num];
		double[] indepProb2s = new double[num];
		
		for (int i=0; i<num; i++) {
			indepProb1s[i] = freq1s[i] / tot1s[i];
			indepProb2s[i] = freq2s[i] / tot2s[i];
		}
		
		for (int x=0; x<num; x++) {
			for (int y=0; y<num; y++) {
				int[] indices = { x, y };
				
				Double tot = totalStatesDataset.get(indices);
				if (tot == null)
					continue;
				
				PossibleStates possible = stateTransitionDataset.get(indices);
				
				int[] ind_1_n2 = { 0, y+1 };
				int[] ind_n1_2 = { x+1, 0 };
				int[] ind_1_2 = { 0, 0 };
				int[] ind_n1_n2 = { x+1, y+1 };
				
				double prob_1_n2 = possible.getFrequency(ind_1_n2)/tot;
				double prob_n1_2 = possible.getFrequency(ind_n1_2)/tot;
				double prob_1_2 = possible.getFrequency(ind_1_2)/tot;
				double prob_n1_n2 = possible.getFrequency(ind_n1_n2)/tot;
				
				p_1_n2.set(x, y, prob_1_n2);
				p_n1_2.set(x, y, prob_n1_2);
				p_1_2.set(x, y, prob_1_2);
				p_n1_n2.set(x, y, prob_n1_n2);
				
//				double prob_1 = prob1s[x];
//				double prob_2 = prob2s[y];
				
				double prob_1 = prob_1_n2 + prob_1_2;
				double prob_2 = prob_n1_2 + prob_1_2;
				
				double synch = prob_1_2/(prob_1*prob_2);
				if (Double.isInfinite(synch))
					synch = Double.NaN;
				
				double prob_state = tot / totStateCount;
				
				double synchTimesProb = synch*prob_state;
//				Preconditions.checkState(Double.isNaN(synch) || (synch >= 0 && synch <= 1),
//						"Synch param is bad: "+prob_1_2+"/("+prob_1+"*"+prob_2+") = "+synch+". synch*P(Sk)="+synchTimesProb);
				
				double synch_norm = synch*prob_state;
				
				synchData.set(x, y, synch);
				synchNormData.set(x, y, synch_norm);
				synchOverIndepData1.set(x, y, prob_1/indepProb1s[x]);
				synchOverIndepData2.set(x, y, prob_2/indepProb2s[y]);
			}
		}
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		cpt.setNanColor(Color.WHITE);
		boolean synch_log = true;
		CPT synch_cpt;
		CPT synch_norm_cpt;
		if (synch_log) {
			synchData.log10();
			synchNormData.log10();
			synch_cpt = cpt.rescale(-synchData.getMaxZ(), synchData.getMaxZ());
			synch_norm_cpt = cpt.rescale(-5, synchNormData.getMaxZ());
		} else {
			synch_cpt = cpt.rescale(0d, synchData.getMaxZ());
			synch_norm_cpt = cpt.rescale(0d, synchNormData.getMaxZ());
		}
		CPT logRatioCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-2d, 2d);
		synchOverIndepData1.log10();
		synchOverIndepData2.log10();
		
		System.out.println(synchData.getMaxZ());
		System.out.println(synchNormData.getMaxZ());
		
		String title = name1+" vs "+name2;
		String xAxisLabel = "Years since prev "+name1;
		String yAxisLabel = "Years since prev "+name2;
		
		List<Range> xRanges = Lists.newArrayList(new Range(0d, maxTimeDiff));
		Range yRange = new Range(0d, 0.6d*maxTimeDiff);
		List<Range> yRanges = Lists.newArrayList(yRange, yRange, yRange, yRange);
		
		List<XYZPlotSpec> specs = Lists.newArrayList();
		specs.add(new XYZPlotSpec(p_1_n2, cpt, title, xAxisLabel, yAxisLabel, null));
		specs.add(new XYZPlotSpec(p_n1_2, cpt, title, xAxisLabel, yAxisLabel, null));
		specs.add(new XYZPlotSpec(p_1_2, cpt, title, xAxisLabel, yAxisLabel, null));
		specs.add(new XYZPlotSpec(p_n1_n2, cpt, title, xAxisLabel, yAxisLabel, null));
		
		XYZGraphPanel gp = new XYZGraphPanel();
		gp.drawPlot(specs, false, false, xRanges, yRanges, null);
		gp.getChartPanel().setSize(1000, 2500);
		gp.saveAsPDF(probFile.getAbsolutePath());
		
		if (name1.contains("Mojave") && name2.contains("Coachella")) {
			// write out individual dists for poster
			specs.get(0).setTitle("P(Mojave ¬Coachella)");
			specs.get(1).setTitle("P(¬Mojave Coachella)");
			specs.get(2).setTitle("P(Mojave Coachella)");
			specs.get(3).setTitle("P(¬Mojave ¬Coachella)");
			
			Range xRange = yRange;
			
			for (int i=0; i<specs.size(); i++) {
				XYZPlotSpec spec = specs.get(i);
				spec.setCPTPosition(RectangleEdge.BOTTOM);
				
				probFile = new File(writeDir, probFName+"_sub"+i+".pdf");
				
				gp = new XYZGraphPanel();
				gp.drawPlot(spec, false, false, xRange, yRange);
				gp.getChartPanel().setSize(1000, 1000);
				gp.saveAsPDF(probFile.getAbsolutePath());
			}
		}
		
		specs = Lists.newArrayList();
		title = "Synchronization: "+title;
		
//		xRanges = Lists.newArrayList(new Range(0d, 0.75*maxTimeDiff));
//		yRange = new Range(0d, maxTimeDiff);
//		yRanges = Lists.newArrayList(yRange, yRange);
		
		specs.add(new XYZPlotSpec(synchData, synch_cpt, title, xAxisLabel, yAxisLabel, null));
		specs.add(new XYZPlotSpec(synchNormData, synch_norm_cpt, title, xAxisLabel, yAxisLabel, null));
		specs.add(new XYZPlotSpec(synchOverIndepData1, logRatioCPT, title, xAxisLabel, yAxisLabel, null));
		specs.add(new XYZPlotSpec(synchOverIndepData2, logRatioCPT, title, xAxisLabel, yAxisLabel, null));
		
		gp = new XYZGraphPanel();
		gp.drawPlot(specs, false, false, xRanges, yRanges, null);
		gp.getChartPanel().setSize(1000, 2500);
		gp.saveAsPDF(synchFile.getAbsolutePath());
		
		return Lists.newArrayList(probFile, synchFile);
	}
	
	private double findGoodLogMin(EvenlyDiscrXYZ_DataSet logData, double absMin) {
		double minAbove = Double.POSITIVE_INFINITY;
		for (int i=0; i<logData.size(); i++) {
			double z = logData.get(i);
			if (z >= absMin && z < minAbove)
				minAbove = z;
		}
		return minAbove;
	}
	
	private void writeTransitionStats(File writeDir) {
//		NoCollissionFunc noRupProbs = new NoCollissionFunc();
//		NoCollissionFunc rupProbs = new NoCollissionFunc();
		HistogramFunction noRupProbs = new HistogramFunction(0.025, 20, 0.05d);
		HistogramFunction rupProbs = new HistogramFunction(0.025, 20, 0.05d);

		HistogramFunction allCounts = new HistogramFunction(0d, 10, 1d);
		HistogramFunction noRupCounts = new HistogramFunction(0d, 10, 1d);
		HistogramFunction rupCounts = new HistogramFunction(0d, 10, 1d);
		
		int numStates = 0;
		int numStatesWithRup = 0;
		int numStatesWithMulti = 0;
		
		for (int[] indices : stateTransitionDataset.getPopulatedIndices()) {
			PossibleStates states = stateTransitionDataset.get(indices);
			
			int noRupCount = 0;
			int rupCount = 0;
			
			for (int i=0; i<states.getStates().size(); i++) {
				int[] newState = states.getStates().get(i);
				double freq = states.getFrequency(newState);
				double prob = freq/states.getTot();
				
				if (Ints.contains(newState, 0)) {
					rupProbs.add(prob, 1d);
					rupCount++;
				} else {
					noRupProbs.add(prob, 1d);
					noRupCount++;
				}
			}
			
			numStates++;
			if (rupCount > 0)
				numStatesWithRup++;
			if ((rupCount + noRupCount) > 1)
				numStatesWithMulti++;
			
			if (noRupCount < noRupCounts.size())
				noRupCounts.add(noRupCount, 1d);
			if (rupCount < rupCounts.size())
				rupCounts.add(rupCount, 1d);
			allCounts.add(rupCount+noRupCount, 1d);
		}
		
		float percentMulti = 100f*(float)numStatesWithMulti/(float)numStates;
		float percentRupMulti = 100f*(float)numStatesWithMulti/(float)numStatesWithRup;
		
		System.out.println("States with multiple transitions: "
				+numStatesWithMulti+"/"+numStates+" ("+percentMulti+" %)");
		System.out.println("Rupture trans states with multiple transitions: "
				+numStatesWithMulti+"/"+numStatesWithRup+" ("+percentRupMulti+" %)");

		List<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		new GraphWindow(asList(allCounts), "# Total Transitions Per State", chars);
		new GraphWindow(asList(rupCounts), "# Rupture Transitions Per State", chars);
		new GraphWindow(asList(noRupCounts), "# No Rupture Transitions Per State", chars);
		new GraphWindow(asList(rupProbs), "Rupture Transition Probailities", chars);
		new GraphWindow(asList(noRupProbs), "No Rupture Transition Probailities", chars);
	}
	
	private static List<PlotElement> asList(PlotElement... elems) {
		return Lists.newArrayList(elems);
	}
	
//	private double calcSynch(int m, int n, int lag, )
	
	public static void main(String[] args) throws IOException {
		double minMag = 7;
		double maxMag = 10d;
		
		int[] include_elems = {
				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID,
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID,
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID,
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID,
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID,
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID
				};
		boolean gen_2d_corr_pdfs = false;
		
		RandomDistType origScrambleDist = null;
//		RandomDistType origScrambleDist = RandomDistType.ACTUAL;
		
		RandomDistType randDistType = RandomDistType.STATE_BASED;
		
		File writeDir = new File("/tmp");
		if (!writeDir.exists())
			writeDir.mkdir();
		
		List<RuptureIdentifier> rupIdens = Lists.newArrayList();
		List<Color> colors = Lists.newArrayList();
		
		SimAnalysisCatLoader.loadElemMagIdens(include_elems, rupIdens, colors, minMag, maxMag);
		
		int[] all_elems = {
				ElementMagRangeDescription.SAF_CHOLAME_ELEMENT_ID,
				ElementMagRangeDescription.SAF_CARRIZO_ELEMENT_ID,
				ElementMagRangeDescription.GARLOCK_WEST_ELEMENT_ID,
				ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID,
				ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID,
				ElementMagRangeDescription.SAN_JACINTO__ELEMENT_ID
				};
		List<RuptureIdentifier> allIdens = Lists.newArrayList();
		SimAnalysisCatLoader.loadElemMagIdens(all_elems, allIdens, null, minMag, maxMag);
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, allIdens, false).getEvents();
		
		if (origScrambleDist != null) {
			 events = RandomCatalogBuilder.getRandomResampledCatalog(events, allIdens, origScrambleDist, true);
			 writeDir = new File(writeDir, "state_"+origScrambleDist.getFNameAdd());
			 if (!writeDir.exists())
				 writeDir.mkdir();
		}
		
		writeDir = new File(writeDir, randDistType.getFNameAdd()+"_corr_plots");
		if (!writeDir.exists())
			writeDir.mkdir();

		Map<IDPairing, HistogramFunction[]> origFuncs =
				PeriodicityPlotter.plotACDF_CCDFs(writeDir, events, rupIdens,
						null, null, 2000d, 10d);
		PeriodicityPlotter.	plotACDF_CCDFs(writeDir, events, rupIdens,
				randDistType, origFuncs, 2000d, 10d);
		
		File myCorrCombined = new File(new File(writeDir, randDistType.getFNameAdd()+"_corr_plots"),
				"corr_combined_rand_state_based.pdf");
		File myCorrCombinedNew = new File(new File(writeDir, randDistType.getFNameAdd()+"_corr_plots"),
				rupIdens.size()+"D_corr_combined_rand_state_based.pdf");
		FileUtils.copyFile(myCorrCombined, myCorrCombinedNew);

		//			File subDir = new File(writeDir, "round2");
		//			if (!subDir.exists())
		//				subDir.mkdir();
		//			PeriodicityPlotter.	plotTimeBetweenAllIdens(subDir, rand_events, rupIdens, rupIdenNames, colors,
		//					RandomDistType.MOJAVE_DRIVER, origFuncs, 2000d, 10d);

		System.out.println("DONE");
		
		File pdfDir = new File(writeDir, "state_pdfs");
		if (!pdfDir.exists())
			pdfDir.mkdir();
		
		List<File> pdfs = Lists.newArrayList();
		
		List<List<? extends SimulatorEvent>> matchesLists = Lists.newArrayList();
		for (int i=0; i<rupIdens.size(); i++)
			matchesLists.add(rupIdens.get(i).getMatches(events));
		
		List<File> twoD_corr_pdfs = Lists.newArrayList();
		
		for (int i=0; i<rupIdens.size(); i++) {
			for (int j=i+1; j<rupIdens.size(); j++) {
				String name1 = rupIdens.get(i).getName();
				String name2 = rupIdens.get(j).getName();
				System.out.println("Writing PDF for "+name1+" vs "+name2);
				StateBasedCatalogBuilder builder = new StateBasedCatalogBuilder();
				builder.buildCatalog(events, null, Lists.newArrayList(matchesLists.get(i), matchesLists.get(j)), false);
				pdfs.addAll(builder.write2DDists(pdfDir, i, name1, matchesLists.get(i), j, name2, matchesLists.get(j)));
				
				
				
				// now 2d-only ACDFs/CCDFs
				if (gen_2d_corr_pdfs && rupIdens.size() > 2) {
					List<RuptureIdentifier> subIdens = Lists.newArrayList(rupIdens.get(i), rupIdens.get(j));
					File subWriteDir = new File(writeDir,
							"2d_"+PeriodicityPlotter.getFileSafeString(name1)+"_"+PeriodicityPlotter.getFileSafeString(name2));
					if (!subWriteDir.exists())
						subWriteDir.mkdir();
					origFuncs = PeriodicityPlotter.plotACDF_CCDFs(subWriteDir, events, subIdens,
									null, null, 2000d, 10d);
					PeriodicityPlotter.	plotACDF_CCDFs(subWriteDir, events, subIdens,
							randDistType, origFuncs, 2000d, 10d);
					
					twoD_corr_pdfs.add(new File(new File(subWriteDir, randDistType.getFNameAdd()+"_corr_plots"),
							"corr_combined_"+randDistType.getFNameAdd()+".pdf"));
				}
			}
		}
		
		if (!pdfs.isEmpty())
			PeriodicityPlotter.combinePDFs(pdfs, new File(pdfDir, "state_dists.pdf"));
		
		if (!twoD_corr_pdfs.isEmpty())
			PeriodicityPlotter.combinePDFs(twoD_corr_pdfs, new File(
					new File(writeDir, randDistType.getFNameAdd()+"_corr_plots"), "2D_corr_combined_rand_state_based.pdf"));
		
		
		StateBasedCatalogBuilder builder = new StateBasedCatalogBuilder();
		builder.buildCatalog(events, null, matchesLists, false);
//		builder.writeTransitionStats(null);
		
		// testing time
		int m = Ints.indexOf(include_elems, ElementMagRangeDescription.SAF_MOJAVE_ELEMENT_ID);
		int n = Ints.indexOf(include_elems, ElementMagRangeDescription.SAF_COACHELLA_ELEMENT_ID);
		List<List<? extends SimulatorEvent>> subMatchesLists = Lists.newArrayList();
		subMatchesLists.add(matchesLists.get(m));
		subMatchesLists.add(matchesLists.get(n));
		StateBasedCatalogBuilder test = new StateBasedCatalogBuilder();
		test.buildCatalog(events, null, subMatchesLists, false);
		
		// first bin by only the indices we care about
		Map<IndicesKey, List<int[]>> binnedIndices = Maps.newHashMap();
		for (int[] indices : builder.totalStatesDataset.getPopulatedIndices()) {
			int[] myInd = { indices[m], indices[n] };
			IndicesKey key = new IndicesKey(myInd);
			List<int[]> binned = binnedIndices.get(key);
			if (binned == null) {
				binned = Lists.newArrayList();
				binnedIndices.put(key, binned);
			}
			binned.add(indices);
		}
		
		// first bin by only the indices we care about
		Map<IndicesKey, List<int[]>> testBinnedIndices = Maps.newHashMap();
		for (int[] indices : test.totalStatesDataset.getPopulatedIndices()) {
			int[] myInd = { indices[0], indices[1] };
			IndicesKey key = new IndicesKey(myInd);
			List<int[]> binned = testBinnedIndices.get(key);
			if (binned == null) {
				binned = Lists.newArrayList();
				testBinnedIndices.put(key, binned);
			}
			binned.add(indices);
		}
		
		System.out.println("6D has "+binnedIndices.size());
		System.out.println("2D has "+testBinnedIndices.size());
		
		// now find discrepancies
		Joiner j = Joiner.on(",");
		for (IndicesKey key : binnedIndices.keySet()) {
			if (!testBinnedIndices.containsKey(key)) {
				List<IndicesKey> subInds = Lists.newArrayList();
				for (int[] indices : binnedIndices.get(key))
					subInds.add(new IndicesKey(indices));
				System.out.println("6-D yes 2-D no: "+key+"\tvals: "+j.join(subInds));
			}
		}
		
		// now find discrepancies
		for (IndicesKey key : testBinnedIndices.keySet()) {
			if (!binnedIndices.containsKey(key)) {
				List<IndicesKey> subInds = Lists.newArrayList();
				for (int[] indices : testBinnedIndices.get(key))
					subInds.add(new IndicesKey(indices));
				System.out.println("6-D no 2-D yes: "+key+"\tvals: "+j.join(subInds));
			}
		}
		
		System.out.println("Done!");
	}
	
}