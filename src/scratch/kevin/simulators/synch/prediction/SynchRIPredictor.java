package scratch.kevin.simulators.synch.prediction;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.jfree.data.Range;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.markov.MarkovChain;
import scratch.kevin.markov.OccupancyBasedMarkovChain2D;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.SynchIdens.SynchFaults;
import scratch.kevin.simulators.synch.StateSpacePlotter;
import scratch.kevin.simulators.synch.TestOccFromSynch;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

public class SynchRIPredictor implements Predictor {
	
	private int maxLag;
	private int nDims;
	private double distSpacing;
	
	private Map<IDPairing, SynchRunningCalc[]> synchCalcsMap;
	private Map<IDPairing, EvenlyDiscrXYZ_DataSet> runningOccsMap;
	
	private RecurrIntervalPredictor riPredict;
	
	private int[] prevState;
	
	private HistogramFunction[] lnGainHists;
	
	public SynchRIPredictor(int maxLag) {
		this.maxLag = maxLag;
	}

	@Override
	public String getShortName() {
		return "SynchRI";
	}

	@Override
	public String getName() {
		return "Synch/RI Hybrid";
	}

	@Override
	public void init(List<int[]> path, double distSpacing) {
		synchCalcsMap = Maps.newHashMap();
		runningOccsMap = Maps.newHashMap();
		nDims = path.get(0).length;
		int num = maxLag+1;
		this.distSpacing = distSpacing;
		for (int i=0; i<nDims; i++) {
			for (int j=0; j<nDims; j++) {
				if (i == j)
					continue;
				IDPairing pair = new IDPairing(i, j);
				SynchRunningCalc[] calcs = new SynchRunningCalc[num];
				for (int l=0; l<=maxLag; l++)
					calcs[l] = new SynchRunningCalc(i, j, l);
				synchCalcsMap.put(pair, calcs);
				if (i < j)
					runningOccsMap.put(pair, new EvenlyDiscrXYZ_DataSet(
							num, num, 0.5*distSpacing, 0.5*distSpacing, distSpacing));
			}
		}
		lnGainHists = new HistogramFunction[nDims];
		for (int i=0; i<nDims; i++)
			lnGainHists[i] = new HistogramFunction(-2.0d, 41, 0.1d);
		
		for (int[] state : path)
			addState(state);
		
		riPredict = new RecurrIntervalPredictor();
		riPredict.init(path, distSpacing);
	}

	@Override
	public void addState(int[] state) {
		if (riPredict != null)
			// will be null during init
			riPredict.addState(state);
		
		for (int i=0; i<nDims; i++) {
			for (int j=i+1; j<nDims; j++) {
				EvenlyDiscrXYZ_DataSet occ = runningOccsMap.get(new IDPairing(i, j));
				if (state[i] < occ.getNumX() && state[j] < occ.getNumY())
					occ.set(state[i], state[j], occ.get(state[i], state[j])+1d);
			}
		}
		
		for (SynchRunningCalc[] calcs : synchCalcsMap.values()) {
			for (SynchRunningCalc calc : calcs)
				calc.addState(state);
		}
		
		prevState = state;
	}

	@Override
	public double[] getRuptureProbabilities() {
		return getRuptureProbabilities(prevState);
	}

//	@Override
//	public double[] getRuptureProbabilities(int[] state) {
//		Preconditions.checkArgument(state != null);
//		double[] riProbs = riPredict.getRuptureProbabilities(state);
//		
//		double[] ret = Arrays.copyOf(riProbs, riProbs.length);
//		
//		for (int i=0; i<state.length; i++) {
//			double riProb = riProbs[i];
//			if (riProb == 0)
//				continue;
//			double multRate = riProb;
//			double avgProb = 0d;
//			
//			double lnSumG = 0;
//			for (int j=0; j<nDims; j++) {
//				if (i == j)
//					continue;
//				
//				IDPairing pair = new IDPairing(i, j);
////				IDPairing pair = new IDPairing(j, i); // this was a test, bad
//				
//				// we need synchronization with the destination states for the other fault
//				// not the current state
//				
//				// this is the lag of the destination state, assuming that it doesn't rupture
//				int l = state[j]+1;
//				// if it doesn't rupture, this is the synchronization
//				double g_norup;
//				if (l > maxLag)
//					g_norup = 1d;
//				else
//					g_norup = synchCalcsMap.get(pair)[l].getCatalogG();
//				Preconditions.checkState(Doubles.isFinite(g_norup));
//				// synchronization if it does rupture
//				double g_rup = synchCalcsMap.get(pair)[0].getCatalogG();
//				Preconditions.checkState(Doubles.isFinite(g_rup));
//				double probJRup = riProbs[j];
////				lnSumG += Math.log(g_norup)*(1-probJRup) + Math.log(g_rup)*probJRup;
//				lnSumG += Math.log(g_norup)*(1-probJRup)*riProb + Math.log(g_rup)*probJRup*riProb;
////				lnSumG += Math.log(g_rup)*probRup;
////				if (Math.random() <= probRup)
////					lnSumG += Math.log(g_rup);
////				else
////					lnSumG += Math.log(g_norup);
//				
////				double avgGain = g_rup*probJRup + g_norup*(1-probJRup);
////				multRate *= avgGain;
//				double avgGain = Math.exp(Math.log(g_rup)*probJRup + Math.log(g_norup)*(1-probJRup));
//				multRate *=  Math.exp(Math.log(g_rup)*probJRup) * Math.exp(Math.log(g_norup)*(1-probJRup));
////				avgProb += riProb*avgGain;
//				
//				g_rup = Math.exp(Math.log(g_rup)*1.0);
//				g_norup = Math.exp(Math.log(g_norup)*1.0);
//				avgProb += g_rup*riProb*probJRup + g_norup*riProb*(1d-probJRup);
//			}
//			
//			lnSumGHists[i].add(lnSumGHists[i].getClosestXIndex(lnSumG), 1d);
////			double prob = Math.exp(Math.log(riProb)+lnSumG);
////			double prob = multRate;
////			double prob = avgProb / (double)(nDims-1);
////			double prob = avgProb / (double)(nDims-1);
//			double prob = avgProb;
//			
//			prob = 1-Math.exp(-prob);
//			
//			Preconditions.checkState(prob >= 0 && prob <= 1,
//					"Bad probability: "+prob+", riProb="+riProb+", lnSumG="+lnSumG);
//			
//			ret[i] = prob;
//		}
//		
//		return ret;
//	}

//	@Override
//	public double[] getRuptureProbabilities(int[] state) {
//		Preconditions.checkArgument(state != null);
//		double[] riProbs = riPredict.getRuptureProbabilities(state);
//		
//		double[] ret = Arrays.copyOf(riProbs, riProbs.length);
//		
//		for (int i=0; i<state.length; i++) {
//			double riProb = riProbs[i];
//			if (riProb == 0)
//				continue;
//			
//			double gain = 1d;
//			
//			for (int j=0; j<nDims; j++) {
//				if (i == j)
//					continue;
//				
//				IDPairing pair = new IDPairing(i, j);
//				
//				// we need synchronization with the destination states for the other fault
//				// not the current state
//				
//				// this is the lag of the destination state, assuming that it doesn't rupture
//				int l = state[j]+1;
//				// if it doesn't rupture, this is the synchronization
//				double g_norup;
//				if (l > maxLag)
//					g_norup = 1d;
//				else
//					g_norup = synchCalcsMap.get(pair)[l].getCatalogG();
//				Preconditions.checkState(Doubles.isFinite(g_norup));
//				// synchronization if it does rupture
//				double g_rup = synchCalcsMap.get(pair)[0].getCatalogG();
//				Preconditions.checkState(Doubles.isFinite(g_rup));
//				double probJRup = riProbs[j];
//				
//				double probJNoRup = 1d-probJRup;
//				
//				// average the gain using the possible destination states
//				double myGain = Math.exp(Math.log(g_rup)*probJRup + Math.log(g_norup)*probJNoRup);
//				
//				Preconditions.checkState(Doubles.isFinite(myGain), "Gain not finite! g_rup="+g_rup
//						+", g_norup="+g_norup+", probJRup="+probJRup+", probJNoRup="+probJNoRup);
//				gain *= myGain;
//			}
//			
//			lnGainHists[i].add(lnGainHists[i].getClosestXIndex(Math.log(gain)), 1d);
////			double prob = Math.exp(Math.log(riProb)+lnSumG);
////			double prob = multRate;
////			double prob = avgProb / (double)(nDims-1);
////			double prob = avgProb / (double)(nDims-1);
////			gain = Math.exp(Math.log(gain)*0.2);
//			double prob = riProb*gain;
//			
////			prob = 1-Math.exp(-prob);
//			if (prob > 1d)
//				prob = 1d;
//			
//			Preconditions.checkState(prob >= 0 && prob <= 1,
//					"Bad probability: "+prob+", riProb="+riProb+", gain="+gain);
//			
//			ret[i] = prob;
//		}
//		
//		return ret;
//	}
	
	public double[] getRuptureProbabilities(int[] state) {
		Preconditions.checkState(nDims == 2, "Only implemented for 2 dimensions!");
		
		MarkovChain chain = getChain();
		return MarkovPredictor.getRuptureProbabilities(chain, riPredict, state);
	}
	
	private MarkovChain chain;
	
	private synchronized MarkovChain getChain() {
		if (chain != null)
			// no updates
			return chain;
		
//		int num = maxLag+1;
//		
//		// build sate occupany distribution
//		EvenlyDiscretizedFunc marginalX = riPredict.getOccupancyDist(0, num);
//		EvenlyDiscretizedFunc marginalY = riPredict.getOccupancyDist(1, num);
//		
//		System.out.println("MarginalX: ["+marginalX.getMinY()+", "+marginalX.getMaxY());
//		System.out.println("MarginalY: ["+marginalY.getMinY()+", "+marginalY.getMaxY());
//		
//		EvenlyDiscretizedFunc firstRow = new EvenlyDiscretizedFunc(0.5*distSpacing, num, distSpacing);
//		EvenlyDiscretizedFunc firstCol = new EvenlyDiscretizedFunc(0.5*distSpacing, num, distSpacing);
//		
//		SynchRunningCalc[] xToYSynchs = synchCalcsMap.get(new IDPairing(0, 1));
//		SynchRunningCalc[] yToXSynchs = synchCalcsMap.get(new IDPairing(1, 0));
//		
//		for (int i=0; i<num; i++) {
//			double xVal = xToYSynchs[i].getProbBoth();
//			double yVal = yToXSynchs[i].getProbBoth();
//			if (i == 0)
//				Preconditions.checkState(xVal == yVal, "Should be equal at origin! xVal="+xVal+", yVal="+yVal);
//			
//			firstCol.set(i, xVal);
//			firstRow.set(i, yVal);
//		}
//		EvenlyDiscretizedFunc firstRow = runningOccsMap.
//		
//		List<EvenlyDiscretizedFunc> funcs = Lists.newArrayList();
//		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
//		funcs.add(marginalX);
//		funcs.add(marginalY);
//		funcs.add(firstRow);
//		funcs.add(firstCol);
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
//		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
//		new GraphWindow(funcs, "Test", chars);
		
		EvenlyDiscrXYZ_DataSet origOccXYZ = runningOccsMap.get(new IDPairing(0, 1));
		
//		EvenlyDiscrXYZ_DataSet occ = TestOccFromSynch.buildOcc(firstRow, firstCol, marginalX, marginalY, 100);
		EvenlyDiscrXYZ_DataSet occ = TestOccFromSynch.buildOcc(origOccXYZ.getRow(0), origOccXYZ.getCol(0),
				origOccXYZ.calcMarginalXDist(), origOccXYZ.calcMarginalYDist(), 100);
		
//		System.out.println("Occ Max: "+occ.getMaxZ()+", Sum: "+occ.getSumZ());
//		try {
//			CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0d, occ.getMaxZ());
//			StateSpacePlotter.plot2D(occ, cpt, "X", "Y", "Occ Prob", null, true, null);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		chain = new OccupancyBasedMarkovChain2D(distSpacing, occ);
		
		return chain;
	}

	@Override
	public void printDiagnostics() {
		// do nothing
	}

	@Override
	public Predictor getCollapsed(int... indexes) {
		SynchRIPredictor p = new SynchRIPredictor(maxLag);
		p.riPredict = (RecurrIntervalPredictor)riPredict.getCollapsed(indexes);
		p.nDims = indexes.length;
		p.synchCalcsMap = Maps.newHashMap();
		p.lnGainHists = new HistogramFunction[p.nDims];
		for (int i=0; i<indexes.length; i++) {
			for (int j=0; j<indexes.length; j++) {
				if (i == j)
					continue;
				p.synchCalcsMap.put(new IDPairing(i, j), synchCalcsMap.get(
						new IDPairing(indexes[i], indexes[j])));
			}
			p.lnGainHists[i] = lnGainHists[indexes[i]];
		}
		return p;
	}
	
	public void writePlots(File dir, List<RuptureIdentifier> idens) throws IOException {
		Preconditions.checkState((dir.exists() && dir.isDirectory()) || dir.mkdir());
		
		for (int i=0; i<nDims; i++) {
			String name1 = idens.get(i).getName();
			for (int j=0; j<nDims; j++) {
				if (i == j)
					continue;
				String name2 = idens.get(j).getName();
				String prefix = PeriodicityPlotter.getFileSafeString(name1)
						+"_"+PeriodicityPlotter.getFileSafeString(name2);
				SynchRunningCalc[] synchCalcs = synchCalcsMap.get(new IDPairing(i, j));
				EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(0d, synchCalcs.length, distSpacing);
				for (int l=0; l<synchCalcs.length; l++)
					func.set(l, Math.log(synchCalcs[l].getCatalogG()));
				List<EvenlyDiscretizedFunc> funcs = Lists.newArrayList();
				List<PlotCurveCharacterstics> chars = Lists.newArrayList();
				funcs.add(func);
				chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				PlotSpec spec = new PlotSpec(funcs, chars, name1+", "+name2, "Lag (years)", "Ln(Gain)");
				
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setTickLabelFontSize(18);
				gp.setAxisLabelFontSize(20);
				gp.setPlotLabelFontSize(21);
				gp.setBackgroundColor(Color.WHITE);
				
				gp.setUserBounds(null, new Range(-2, 2));
				gp.drawGraphPanel(spec);
				gp.getChartPanel().setSize(800, 400);
				gp.saveAsPNG(new File(dir, prefix+".png").getAbsolutePath());
				gp.saveAsPDF(new File(dir, prefix+".pdf").getAbsolutePath());
			}
		}
		
		File histDir = new File(dir, "gain_hists");
		Preconditions.checkState((histDir.exists() && histDir.isDirectory()) || histDir.mkdir());
		
		for (int i=0; i<nDims; i++) {
			String name = idens.get(i).getName();
			
			List<EvenlyDiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			funcs.add(lnGainHists[i]);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 2f, Color.BLACK));
			PlotSpec spec = new PlotSpec(funcs, chars, name, "Gain", "Number");
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setTickLabelFontSize(18);
			gp.setAxisLabelFontSize(20);
			gp.setPlotLabelFontSize(21);
			gp.setBackgroundColor(Color.WHITE);
			
			gp.drawGraphPanel(spec);
			gp.getChartPanel().setSize(800, 400);
			gp.saveAsPNG(new File(histDir,
					PeriodicityPlotter.getFileSafeString(name)+".png").getAbsolutePath());
		}
	}
	
	private class SynchRunningCalc {
		
		private int m, n, lag;
		
		private int numWindows = 0;
		private int numMN=0, numM=0, numN=0;
		
		private ArrayDeque<Integer> nDeque;
		
		private int[] prevState;
		
		private static final boolean CALC_DIAG = false;
		
		public SynchRunningCalc(int m, int n, int lag) {
			this.m = m;
			this.n = n;
			this.lag = lag;
			Preconditions.checkState(lag >= 0);
			if (lag != 0)
				// this will store previous states for index n
				nDeque = new ArrayDeque<Integer>(maxLag);
		}
		
		public void addState(int[] state) {
			int[] newState = { state[m], state[n]};
			
			if (CALC_DIAG) {
				if (prevState == null) {
					prevState = state;
					return;
				}
				int delta = prevState[1] - prevState[0];
				if (delta != lag) {
					prevState = state;
					return;
				}
			}
			
			if (lag > 0) {
				nDeque.addLast(state[n]);
				if (nDeque.size() > lag)
					// we have enough
					newState[1] = nDeque.removeFirst();
				else
					return;
			}
			numWindows++;
			if (newState[0] == 0)
				numM++;
			if (newState[1] == 0)
				numN++;
			if (newState[0] == 0 && newState[1] == 0)
				numMN++;
			
			prevState = newState;
		}
		
		public double getCatalogG() {
			if (numWindows == 0 || numM == 0 || numN == 0)
				return 1d;
			return (double)numWindows * (double)numMN/(double)(numM*numN);
		}
		
		public double getProbBoth() {
			if (numWindows == 0 || numM == 0 || numN == 0)
				return 0d;
			return (double)numMN/(double)numWindows;
		}
		
		public double getProbM() {
			if (numWindows == 0 || numM == 0 || numN == 0)
				return 0d;
			return (double)numM/(double)numWindows;
		}
		
		public double getProbN() {
			if (numWindows == 0 || numM == 0 || numN == 0)
				return 0d;
			return (double)numN/(double)numWindows;
		}
		
		public double getProbMgivenN() {
			if (numWindows == 0 || numM == 0 || numN == 0)
				return 0d;
			return getProbBoth()/getProbN();
		}
		
		public double getProbNgivenM() {
			if (numWindows == 0 || numM == 0 || numN == 0)
				return 0d;
			return getProbBoth()/getProbM();
		}
	}
	
	public static void main(String[] args) throws IOException {
		double distSpacing = 10d;
		List<RuptureIdentifier> rupIdens = SynchIdens.getIndividualFaults(7d, 10d,
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO, SynchFaults.SAF_CARRIZO);
//				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
				SynchFaults.SAF_MOJAVE, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_CARRIZO, SynchFaults.SAF_COACHELLA);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAN_JACINTO);
//				SynchFaults.SAF_COACHELLA, SynchFaults.SAF_CHOLAME);
//				SynchFaults.SAF_MOJAVE, SynchFaults.GARLOCK_WEST);
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, rupIdens, true).getEvents();
		
		List<int[]> fullPath = MarkovChainBuilder.getStatesPath(distSpacing, events, rupIdens, 0d);
		
		SynchRIPredictor predict = new SynchRIPredictor(100);
		predict.init(fullPath, distSpacing);
		predict.getChain();
	}

}
