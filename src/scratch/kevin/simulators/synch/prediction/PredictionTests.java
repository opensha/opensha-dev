package scratch.kevin.simulators.synch.prediction;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.plot.XYPlot;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphPanel;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotPreferences;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotWindow;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;
import org.opensha.sha.simulators.utils.General_EQSIM_Tools;

import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.OccBasedIterativeSolver;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.SynchIdens.SynchFaults;
import scratch.kevin.simulators.dists.LogNormalDistReturnPeriodProvider;
import scratch.kevin.simulators.synch.StateSpacePlotter;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

public class PredictionTests {
	
	private List<Predictor> predictors;
	private int nullHypothesisIndex;
	private double distSpacing;
	
	private List<List<double[]>> predictions;
	private List<double[]> predictedNs;
	private List<int[]> actual;
	private int[] actualNs;
	private int actualN;
	
	private int nDims;
	
	private double minRate = 1e-1d;
	private boolean skipZero = false;
	
	public PredictionTests(List<Predictor> predictors, int nullHypothesisIndex, double distSpacing) {
		Preconditions.checkState(predictors.size() >= 2, "Must supply at least 2 predictors");
		Preconditions.checkState(nullHypothesisIndex >=0 && nullHypothesisIndex<predictors.size(), "Must supply null hypothesis");
		Preconditions.checkState(distSpacing>0, "dist spacing must be >0");
		
		this.predictors = predictors;
		this.nullHypothesisIndex = nullHypothesisIndex;
		this.distSpacing = distSpacing;
	}
	
	public void doTests(List<int[]> fullPath, int learningIndex) {
		List<int[]> initialPath = fullPath.subList(0, learningIndex);
		
		nDims = fullPath.get(0).length;
		
		predictions = Lists.newArrayList();
		predictedNs = Lists.newArrayList();
		actual = Lists.newArrayList();
		
		System.out.println("Preparing "+predictors.size()+" predictors with "+learningIndex+" learning states");
		
		for (Predictor p : predictors) {
			p.init(initialPath, distSpacing);
			predictions.add(new ArrayList<double[]>());
			predictedNs.add(new double[nDims]);
		}
		
		System.out.println("Doing predictions");
		for (int i=learningIndex; i<fullPath.size(); i++) {
			int[] state = fullPath.get(i);
			actual.add(state);
			
			for (int j=0; j<predictors.size(); j++) {
				Predictor p = predictors.get(j);
				double[] prediction = p.getRuptureProbabilities();
				for (int k=0; k<prediction.length; k++)
					Preconditions.checkState(Doubles.isFinite(prediction[k]));
				double[] predictedN = predictedNs.get(j);
				for (int k=0; k<nDims; k++)
					predictedN[k] += prediction[k];
				predictions.get(j).add(prediction);
				p.addState(state);
			}
		}
		
		actualNs = new int[nDims];
		for (int[] state : actual)
			for (int i=0; i<nDims; i++)
				if (state[i] == 0)
					actualNs[i]++;
		actualN = 0;
		for (int n : actualNs)
			actualN += n;
		
		System.out.println("Done with "+actual.size()+" predictions");
		
		System.out.println("Actual N: "+actualN+", ["+Joiner.on(",").join(Ints.asList(actualNs))+"]");
		
		List<Double> gains = Lists.newArrayList();
		for (int i=0; i<predictors.size(); i++)
			gains.add(calcInformationGain(i));
		
		List<ComparablePairing<Double, Predictor>> sorted = ComparablePairing.build(gains, predictors);
		Collections.sort(sorted);
		Collections.reverse(sorted);
		System.out.println("\n***SORTED***");
		for (ComparablePairing<Double, Predictor> p : sorted) {
			double[] predictedNs = this.predictedNs.get(predictors.indexOf(p.getData()));
			double[] nCompare = new double[actualNs.length];
			for (int i=0; i<actualNs.length; i++)
				nCompare[i] = 100d*predictedNs[i]/actualNs[i];
			System.out.println("\t"+p.getComparable().floatValue()+":\t"+p.getData().getName()
					+"\t\tCount %: ["+Joiner.on(",").join(getFloatList(nCompare))+"]");
		}
	}
	
	private double calcInformationGain(int predictorIndex) {
		Predictor predictor = predictors.get(predictorIndex);
		List<double[]> predictions = this.predictions.get(predictorIndex);
		double[] predictedNs = this.predictedNs.get(predictorIndex);
		double[] nCompare = new double[actualNs.length];
		for (int i=0; i<actualNs.length; i++)
			nCompare[i] = 100d*predictedNs[i]/actualNs[i];
		double predictedN = StatUtils.sum(predictedNs);
		System.out.println("Predicted N for "+predictor.getName()+": "+predictedN
				+", ["+Joiner.on(",").join(getFloatList(predictedNs))
				+"], %: ["+Joiner.on(",").join(getFloatList(nCompare))+"]");
		predictor.printDiagnostics();

		List<double[]> nullPredictions = this.predictions.get(nullHypothesisIndex);
		double[] nullPredictedNs = this.predictedNs.get(nullHypothesisIndex);
		double nullPredictedN = StatUtils.sum(nullPredictedNs);
		
		// A = predictor, B = null predictor
		// I(A,B) = (1/N) * sumOverN(log(rateA) - log(rateB)) - (NA - NB)/N
		
		double sumLogRateDiff = 0;
		
		double adjustedActualN = 0;
		
		for (int i=0; i<actual.size(); i++) {
			int[] state = actual.get(i);
			for (int n=0; n<nDims; n++) {
				if (state[n] == 0) {
					// it's a rupture
					double rateA = Math.max(predictions.get(i)[n], minRate);
					double rateB =  Math.max(nullPredictions.get(i)[n], minRate);
					
					if (skipZero && (rateA == 0 || rateB == 0))
						continue;
					
					sumLogRateDiff += Math.log(rateA)-Math.log(rateB);
					
					adjustedActualN++;
				}
			}
		}
		
		double ig = (1d/(double)adjustedActualN)*sumLogRateDiff - (predictedN - nullPredictedN)/(double)adjustedActualN;
		System.out.println("Information gain for "+predictor.getName()+" relative to "
				+predictors.get(nullHypothesisIndex).getName()+": "+ig
				+" (evaluated for "+adjustedActualN+"/"+actualN+" = "+100f*(float)(adjustedActualN/(double)actualN)+" %)");
		
		return ig;
	}
	
	private void write2DProbPlots(File outputDir, List<RuptureIdentifier> rupIdens) throws IOException {
		Preconditions.checkState((outputDir.exists() && outputDir.isDirectory()) || outputDir.mkdir());
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, 1d);
		
		for (Predictor p : predictors) {
			File subDir = new File(outputDir, p.getShortName());
			Preconditions.checkState((subDir.exists() && subDir.isDirectory()) || subDir.mkdir());
			
			for (int m=0; m<nDims; m++) {
				for (int n=m+1; n<nDims; n++) {
					Predictor p2d;
					if (nDims == 2)
						p2d = p;
					else
						p2d = p.getCollapsed(m, n);
					EvenlyDiscrXYZ_DataSet probX = new EvenlyDiscrXYZ_DataSet(
							100, 100, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
					EvenlyDiscrXYZ_DataSet probY = new EvenlyDiscrXYZ_DataSet(
							100, 100, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
					
					String name1 = rupIdens.get(m).getName();
					String name2 = rupIdens.get(n).getName();
					
					for (int i=0; i<probX.getNumX(); i++) {
						for (int j=0; j<probX.getNumY(); j++) {
							int[] state = {i,j};
							double[] probs = p2d.getRuptureProbabilities(state);
							probX.set(i, j, probs[0]);
							probY.set(i, j, probs[1]);
						}
					}

					XYZPlotSpec xSpec = new XYZPlotSpec(probX, cpt, p2d.getName()+" Rup Probs",
							name1+" OI", name2+" OI", name1+" Probability");
					XYZPlotSpec ySpec = new XYZPlotSpec(probY, cpt, p2d.getName()+" Rup Probs",
							name1+" OI", name2+" OI", name2+" Probability");
					
					List<XYZPlotSpec> specs = Lists.newArrayList(xSpec, ySpec);
					
					XYZGraphPanel xyzGP = new XYZGraphPanel();
					xyzGP.drawPlot(specs, false, false, null, null);
					xyzGP.getChartPanel().setSize(1000, 2000);
					String prefix = PeriodicityPlotter.getFileSafeString(name1)
							+"_"+PeriodicityPlotter.getFileSafeString(name2);
					xyzGP.saveAsPDF(new File(subDir, prefix+".pdf").getAbsolutePath());
					xyzGP.saveAsPNG(new File(subDir, prefix+".png").getAbsolutePath());
				}
			}
		}
	}
	
	private static List<Float> getFloatList(double[] array) {
		List<Float> l = Lists.newArrayList();
		for (double val : array)
			l.add((float)val);
		return l;
	}
	
	/**
	 * Generates a fake series of events where the first fault is on a log-normal renewal model, and the second
	 * fault happens with every second rupture of the first.
	 * @param events
	 * @param rupIdens
	 * @return
	 */
	private static List<SimulatorEvent> generateFakeData1(List<SimulatorEvent> events, List<RuptureIdentifier> rupIdens,
			long seed) {
		Preconditions.checkState(rupIdens.size() == 2);
		RuptureIdentifier iden1 = rupIdens.get(0);
		RuptureIdentifier iden2 = rupIdens.get(1);
		
		// first find all unique events for each
		List<SimulatorEvent> events1 = iden1.getMatches(events);
		List<SimulatorEvent> events2 = iden2.getMatches(events);
		
		LogNormalDistReturnPeriodProvider riProv =
				new LogNormalDistReturnPeriodProvider(PeriodicityPlotter.getRPs(events1));
		riProv.setSeed(seed);
		
		HashSet<SimulatorEvent> events2Hash = new HashSet<SimulatorEvent>(events2);
		
		// remove duplicates
		for (int i=events1.size(); --i>=0;) {
			SimulatorEvent event = events1.get(i);
			if (events2Hash.contains(event)) {
				events1.remove(i);
				Preconditions.checkState(events2.remove(event));
			}
		}
		
		List<SimulatorEvent> fakeEvents = Lists.newArrayList();
		
		boolean includeFault2 = true;
		
		int eventID = 0;
		
		double time = 0d;
		while (!events1.isEmpty() && !events2.isEmpty()) {
			double riYears = riProv.getReturnPeriod();
			time += riYears*General_EQSIM_Tools.SECONDS_PER_YEAR;
			fakeEvents.add(events1.remove(0).cloneNewTime(time, eventID++));
			if (includeFault2)
				fakeEvents.add(events2.remove(0).cloneNewTime(time, eventID++));
			includeFault2 = !includeFault2;
		}
		
		System.out.println("Fake catalog has "+fakeEvents.size()
				+" events, length: "+(time/General_EQSIM_Tools.SECONDS_PER_YEAR));
		
		Preconditions.checkState(!fakeEvents.isEmpty());
		
		return fakeEvents;
	}
	
	private static MultivariateNormalDistribution getBivariateNormal(double corr, double mean1, double sd1,
			double mean2, double sd2, long seed) {
		double[] means = { Math.log(mean1), Math.log(mean2) };
		double[][] covariances = new double[2][2];
		covariances[0][0] = sd1*sd1;
		covariances[1][1] = sd2*sd2;
		covariances[0][1] = corr*sd1*sd2;
		covariances[1][0] = corr*sd1*sd2;
		
		MultivariateNormalDistribution dist = new MultivariateNormalDistribution(means, covariances);
		dist.reseedRandomGenerator(seed);
		return dist;
	}
	
	private static EvenlyDiscrXYZ_DataSet integrateBVNormal(MultivariateNormalDistribution dist, int numSims,
			int numBins, double deltaBin) {
		EvenlyDiscrXYZ_DataSet xyz = new EvenlyDiscrXYZ_DataSet(numBins, numBins, 0.5*deltaBin, 0.5*deltaBin, deltaBin);
		
		int outside = 0;
		
		double[] means = new double[2];
		
		for (int i=0; i<numSims; i++) {
			double[] sample = dist.sample();
			means[0] += sample[0];
			means[1] += sample[1];
			sample[0] = Math.exp(sample[0]);
			sample[1] = Math.exp(sample[1]);
			int xInd = xyz.getXIndex(sample[0]);
			int yInd = xyz.getYIndex(sample[1]);
			if (xInd < 0 || xInd >= xyz.getNumX() || yInd < 0 || yInd >= xyz.getNumY()) {
				outside++;
				continue;
			}
			
			xyz.set(xInd, yInd, xyz.get(xInd, yInd)+1d);
		}
		int count = numSims-outside;
		System.out.println("Integrated at "+count+" points ("+outside+"/"+numSims+" were skipped)");
		means[0] = Math.exp(means[0]/count);
		means[1] = Math.exp(means[1]/count);
		System.out.println("Means: "+means[0]+", "+means[1]);
		xyz.scale(1d/(double)count);
		
		return xyz;
	}
	
	private static EvenlyDiscrXYZ_DataSet stateOccFromBVNormal(EvenlyDiscrXYZ_DataSet bvNormal) {
		EvenlyDiscrXYZ_DataSet occ = new EvenlyDiscrXYZ_DataSet(bvNormal.getNumX(), bvNormal.getNumY(),
				bvNormal.getMinX(), bvNormal.getMinY(), bvNormal.getGridSpacingX(), bvNormal.getGridSpacingY());
		
		// sum down the diagonals
		
		// this is the list of starting states (top and right edge of the map)
		List<int[]> startStates = Lists.newArrayList();
		for (int x=0; x<occ.getNumX(); x++)
			startStates.add(new int[] {x, occ.getNumY()-1});
		for (int y=0; y<occ.getNumY()-1; y++) // don't repeat the top right
			startStates.add(new int[] {occ.getNumX()-1, y});
		Preconditions.checkState(startStates.size() == (occ.getNumX()+occ.getNumY()-1));
		
		double totalSum = 0d;
		
		for (int[] startState : startStates) {
			int x = startState[0];
			int y = startState[1];
			
			double sum = 0d;
			
			while (x >= 0 && y >= 0) {
				sum += bvNormal.get(x, y);
				
				Preconditions.checkState(occ.get(x, y) == 0d);
				totalSum += sum;
				occ.set(x, y, sum);
				
				x--;
				y--;
			}
		}
		
		occ.scale(1d/totalSum);
		
		return occ;
	}
	
	private static void plotBivariateNormal(double corr, double mean1, double sd1,
			double mean2, double sd2, long seed) throws IOException {
		MultivariateNormalDistribution dist = getBivariateNormal(corr, mean1, sd1, mean2, sd2, seed);
		
		int totNum = 10000000;
		EvenlyDiscrXYZ_DataSet bvNormal = integrateBVNormal(dist, totNum, 50, 10d);
		EvenlyDiscrXYZ_DataSet occ = stateOccFromBVNormal(bvNormal);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, bvNormal.getMaxZ());
		XYZPlotSpec spec = new XYZPlotSpec(bvNormal, cpt, "Bivariate Normal Integration", "X", "Y", "Density");
		
		XYZPlotWindow gw = new XYZPlotWindow(spec);
		gw.setDefaultCloseOperation(XYZPlotWindow.EXIT_ON_CLOSE);
		
		// now occupancy
		cpt = cpt.rescale(0d, occ.getMaxZ());
		spec = new XYZPlotSpec(occ, cpt, "Occupancy", "X", "Y", "Occupancy");
		
		gw = new XYZPlotWindow(spec);
		gw.setDefaultCloseOperation(XYZPlotWindow.EXIT_ON_CLOSE);
		
		// now 1D RI dists
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(bvNormal.calcMarginalXDist());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		funcs.add(bvNormal.calcMarginalYDist());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		new GraphWindow(new PlotSpec(funcs, chars, "1D RI Dists", "Years", "Prob"));
		
		// now 1D RI dists
		funcs = Lists.newArrayList();
		chars = Lists.newArrayList();
		
		funcs.add(occ.calcMarginalXDist());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		funcs.add(occ.calcMarginalYDist());
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		new GraphWindow(new PlotSpec(funcs, chars, "1D Occupancy Dists", "Years", "Prob"));
		
		while (true) {
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Generates a fake series of events where the fboth faults are log-normally distributed, with the
	 * given correlation coeff
	 * @param events
	 * @param rupIdens
	 * @return
	 */
	private static List<SimulatorEvent> generateFakeData2(List<? extends SimulatorEvent> events, List<RuptureIdentifier> rupIdens,
			double corr, double mean1, double sd1, double mean2, double sd2, long seed) {
		Preconditions.checkState(rupIdens.size() == 2);
		RuptureIdentifier iden1 = rupIdens.get(0);
		RuptureIdentifier iden2 = rupIdens.get(1);
		
		// first find all unique events for each
		List<? extends SimulatorEvent> events1 = iden1.getMatches(events);
		List<? extends SimulatorEvent> events2 = iden2.getMatches(events);
		
		HashSet<SimulatorEvent> events2Hash = new HashSet<SimulatorEvent>(events2);
		
		// remove duplicates
		for (int i=events1.size(); --i>=0;) {
			SimulatorEvent event = events1.get(i);
			if (events2Hash.contains(event)) {
				events1.remove(i);
				Preconditions.checkState(events2.remove(event));
			}
		}
		
		List<SimulatorEvent> fakeEvents = Lists.newArrayList();
		
		double[] means = { Math.log(mean1), Math.log(mean2) };
		double[][] covariances = new double[2][2];
		covariances[0][0] = sd1*sd1;
		covariances[1][1] = sd2*sd2;
		covariances[0][1] = corr*sd1*sd2;
		covariances[1][0] = corr*sd1*sd2;
		
		MultivariateNormalDistribution dist = new MultivariateNormalDistribution(means, covariances);
		// TODO figure this out
//		dist.
//		double shape1 = sd1 / mean1;
//		double shape2 = sd2 / mean2;
//		LogNormalDistribution n1 = new LogNormalDistribution(Math.log(mean1), shape1);
//		LogNormalDistribution n2 = new LogNormalDistribution(Math.log(mean2), shape2);
//		
//		int eventID = 0;
//		
//		double oi1 = Math.random()*(2d*mean1);
//		double oi2 = Math.random()*(2d*mean2);
//		
//		double time = 0d;
//		while (!events1.isEmpty() && !events2.isEmpty()) {
//			double jointProb = dist.density(new double[] {Math.log(oi1), Math.log(oi2)});
//			double p1 = n1.probability(Math.log(oi1))
//			if (Math.random() <= jointProb) {
//				fakeEvents.add(events1.remove(0).cloneNewTime(time, eventID++));
//				fakeEvents.add(events2.remove(0).cloneNewTime(time, eventID++));
//				oi1 = 0;
//				oi2 = 0;
//			}
//			
//			oi1 += 10d;
//			oi2 += 10d;
//			
//			time += 10d*General_EQSIM_Tools.SECONDS_PER_YEAR;
//		}
//		
//		System.out.println("Fake catalog has "+fakeEvents.size()
//				+" events, length: "+(time/General_EQSIM_Tools.SECONDS_PER_YEAR));
//		
//		Preconditions.checkState(!fakeEvents.isEmpty());
		
		return fakeEvents;
	}

	public static void main(String[] args) throws IOException {
//		double mean1 = 100d;
//		double sd1 = 0.2;
//		double mean2 = 100d;
//		double sd2 = 0.2;
////		double corr = 0.5;
////		double corr = 0.9999999999;
//		double corr = 0;
//		plotBivariateNormal(corr, mean1, sd1, mean2, sd2, 0l);
		
		double minMag = 7d; 
		double maxMag = 10d;
		double distSpacing = 10d; // years
		
		boolean do2DPlots = true;
		File predictDir = new File("/home/kevin/Simulators/predict");
		File plot2DOutputDir = new File(predictDir, "plots_2d");
		File synchPlotOutputDir = new File(predictDir, "synch_plots");
		
		boolean fakeData = false;
		List<RuptureIdentifier> rupIdens = SynchIdens.getIndividualFaults(minMag, maxMag,
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO, SynchFaults.SAF_CARRIZO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_COACHELLA);
				SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_MOJAVE, SynchFaults.GARLOCK_WEST);
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, rupIdens, true).getEvents();
		if (fakeData)
//			events = generateFakeData1(events, rupIdens, 0l);
			events = generateFakeData2(events, rupIdens, 0.5, 100d, 30d, 100d, 30d, 0l);
		
		List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing, events, rupIdens, 0d);
//		if (rupIdens.size() == 2)
//			new StateSpacePlotter(new EmpiricalMarkovChain(fullPath, distSpacing), rupIdens, null).plotOccupancies();
		
		int learningIndex = fullPath.size()/2;
		
		List<Predictor> predictors = Lists.newArrayList();
		predictors.add(new MarkovPredictor(new RecurrIntervalPredictor()));
		predictors.add(new MarkovPredictor());
		predictors.add(new RecurrIntervalPredictor());
		int nullHypothesisIndex = predictors.size()-1;
//		SynchRIPredictor synch = new SynchRIPredictor(50);
//		predictors.add(synch);
//		if (rupIdens.size() == 2 && fakeData)
//			predictors.add(new SplitPredictor(new RecurrIntervalPredictor(), new SynchRIPredictor(50)));
		if (rupIdens.size() == 2) {
			// add occupancy based
			predictors.add(new OccBasedMarkovPredictor(new RecurrIntervalPredictor(), false, 0.9999));
			predictors.add(new OccBasedMarkovPredictor(new RecurrIntervalPredictor(), true, 0.9999));
		}
		predictors.add(new PoissonPredictor());
//		int nullHypothesisIndex = predictors.size()-1;
		
//		double[] minRates = { 0d, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1 };
		double[] minRates = { 0d };
		
		for (double minRate : minRates) {
			System.out.println("****************************");
			System.out.println("min rate: "+minRate);
			System.out.println("****************************");
			
			PredictionTests tests = new PredictionTests(predictors, nullHypothesisIndex, distSpacing);
			
			tests.skipZero = minRate == 0d;
			tests.minRate = minRate;
			
			tests.doTests(fullPath, learningIndex);
			
			if (do2DPlots && minRate == 0d)
				tests.write2DProbPlots(plot2DOutputDir, rupIdens);
		}
//		synch.writePlots(synchPlotOutputDir, rupIdens);
	}

}
