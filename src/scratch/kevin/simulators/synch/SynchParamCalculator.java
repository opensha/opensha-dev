package scratch.kevin.simulators.synch;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.data.Range;
import org.jfree.ui.TextAnchor;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.UncertainArbDiscDataset;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.xyz.EvenlyDiscrXYZ_DataSet;
import org.opensha.commons.gui.plot.HeadlessGraphPanel;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZGraphPanel;
import org.opensha.commons.gui.plot.jfreechart.xyzPlot.XYZPlotSpec;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.ComparablePairing;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.commons.util.threads.Task;
import org.opensha.commons.util.threads.ThreadedTaskComputer;
import org.opensha.sha.simulators.SimulatorEvent;
import org.opensha.sha.simulators.iden.ElementMagRangeDescription;
import org.opensha.sha.simulators.iden.RuptureIdentifier;

import scratch.UCERF3.utils.IDPairing;
import scratch.kevin.DistSpeedTest;
import scratch.kevin.markov.EmpiricalMarkovChain;
import scratch.kevin.markov.IndicesKey;
import scratch.kevin.markov.MarkovChain;
import scratch.kevin.markov.OccupancyBasedMarkovChain2D;
import scratch.kevin.markov.PossibleStates;
import scratch.kevin.markov.SparseNDimensionalHashDataset;
import scratch.kevin.simulators.MarkovChainBuilder;
import scratch.kevin.simulators.PeriodicityPlotter;
import scratch.kevin.simulators.SimAnalysisCatLoader;
import scratch.kevin.simulators.SynchIdens;
import scratch.kevin.simulators.catBuild.RandomCatalogBuilder;
import scratch.kevin.simulators.dists.RandomDistType;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

/**
 * Synchronization parameter calculator for simulated catalogs. Can calculate Gbar with various parameterizations,
 * including CATALOG_G, the most recent/best formulation.
 * @author kevin
 *
 */
public class SynchParamCalculator {

	// inputs
	private EmpiricalMarkovChain chain;
	private int lag;

	// results
	private List<Double> synchs = Lists.newArrayList();
	private List<Double> weights = Lists.newArrayList();
	private List<Double> probMs = Lists.newArrayList();
	private List<Double> probNs = Lists.newArrayList();
	private List<Double> probMNs = Lists.newArrayList();
	private List<Double> freqEithers = Lists.newArrayList();
	private List<IndicesKey> synch_indices = Lists.newArrayList();
	private double totEithers = 0d;

	private double gBar_numerator;
	private double gBar_denominator;
	private double gBar;

	// plot files if generated
	private File synch2DPlotFile;
	private File scatterPlotFile;

	private final int m = 0;
	private final int n = 1;

	//	double cov;

	// if true, gain factor divides by P(Em|k)*P(En|l) instead of P(Em|k,l)*P(En|k,l)
	protected static final boolean useIndepProbs = true;

	// if true, lags are calculated by shifting the markov chain and computing as normal
	protected static final boolean doLagByShift = false;

	// if true, lags are calculated by via state occupancy
	protected static final boolean doLagByOcc = true;

	enum WeightingScheme {
		FREQ_EITHERS,
		TOT_OCCUPANCY,
		CATALOG_G,
		FREQ_MNs,
		FREQ_AVG,
		FREQ_PROD,
		LAG_DEPENDENT_EITHER,
		LAG_DEPENDENT_AVG,
		LAG_DEPENDENT_PROD;
	}

	protected static final WeightingScheme weightingScheme = WeightingScheme.CATALOG_G;

	//	public static Map<IndicesKey, List<int[]>> getBinnedIndices(MarkovChainBuilder chain, int m, int n) {
	//		Map<IndicesKey, List<int[]>> binnedIndices = Maps.newHashMap();
	//		for (int[] indices : chain.getStateTransitionDataset().getPopulatedIndices()) {
	//			int[] myInd = { indices[m], indices[n] };
	//			IndicesKey key = new IndicesKey(myInd);
	//			List<int[]> binned = binnedIndices.get(key);
	//			if (binned == null) {
	//				binned = Lists.newArrayList();
	//				binnedIndices.put(key, binned);
	//			}
	//			binned.add(indices);
	//		}
	//		return binnedIndices;
	//	}

	public SynchParamCalculator(EmpiricalMarkovChain chain, int m, int n, int lag) {
		this(chain.getCollapsedChain(m, n), lag);
	}

	public SynchParamCalculator(EmpiricalMarkovChain chain, int lag) {
		Preconditions.checkArgument(chain.getNDims() == 2, "Chain must be collapsed");
		this.chain = chain;
		this.lag = lag;

		calculate();
	}

	private static boolean lag_debug = false;

	private void calculate() {
		if (weightingScheme == WeightingScheme.CATALOG_G) {
			if (doLagByOcc) {
				gBar = calcOccG(chain, m, catLenMult, lag);
			} else {
				Preconditions.checkState(lag == 0, "Must use shift lag for catalog G");
				gBar = calcCatalogG(chain, 0, 1);
			}
			return;
		}
		
		int numSums = 0;
		int numPossibleSums = 0;
		int numSubSums = 0;
		int numSubBails = 0;

		int totOccupancy = 0;

		Map<Integer, Double> mIndepFreqs;
		Map<Integer, Double> nIndepFreqs;
		Map<Integer, Double> mIndepProbs;
		Map<Integer, Double> nIndepProbs;
		if (lag == 0) {
			mIndepProbs = calcIndepProbs(m, 0, false);
			nIndepProbs = calcIndepProbs(n, 0, false);
			mIndepFreqs = calcIndepProbs(m, 0, true);
			nIndepFreqs = calcIndepProbs(n, 0, true);
		} else if (lag < 0) {
			// n precedes m
			// we want m=0, and n=abs(lag)
			mIndepProbs = calcIndepProbs(m, 0, false);
			nIndepProbs = calcIndepProbs(n, -lag, false);
			mIndepFreqs = calcIndepProbs(m, 0, true);
			nIndepFreqs = calcIndepProbs(n, -lag, true);
		} else {
			// m precedes n
			// we want n=0, and m=lag
			mIndepProbs = calcIndepProbs(m, lag, false);
			nIndepProbs = calcIndepProbs(n, 0, false);
			mIndepFreqs = calcIndepProbs(m, lag, true);
			nIndepFreqs = calcIndepProbs(n, 0, true);
		}

		List<Double> allFreqMs = Lists.newArrayList();
		List<Double> allFreqNs = Lists.newArrayList();
		List<Double> allFreqMNs = Lists.newArrayList();
		List<Double> allFreqEithers = Lists.newArrayList();
		List<Double> allTots = Lists.newArrayList();
		List<int[]> allStates = Lists.newArrayList();
		double allTotOccupancy = 0;

		List<int[]> usedStates = Lists.newArrayList();

		List<Double> weightNumerators = Lists.newArrayList();
		double weightDenominator = 0d;

		for (int[] fromState : chain.getStateTransitionDataset().getPopulatedIndices()) {
			PossibleStates possible = chain.getStateTransitionDataset().get(fromState);
			if (possible == null || possible.getStates() == null || possible.getTot() == 0) {
				// last state in the catalog can be a dead end if never reached earlier
				numSubBails++;
				continue;
			}

			double freqM = 0;
			double freqN = 0;
			double freqMN = 0;
			double freqEither = 0;
			double tot = possible.getTot();

			for (int[] state : possible.getStates()) {
				// frequency that we transition to state
				double freq = possible.getFrequency(state);
				if (lag == 0) {
					// simple case, no lag
					if (state[m] == 0)
						freqM += freq;
					if (state[n] == 0)
						freqN += freq;
					if (state[m] == 0 && state[n] == 0)
						freqMN += freq;
					if (state[m] == 0 || state[n] == 0)
						freqEither += freq;
					//				} else if (useNewTomNotation) {
					//					if (state[m] == 0)
					//						freqM += freq;
					//					if (state[n] == 0)
					//						freqN += freq;
				} else {
					// incorporates lag
					if (lag < 0) {
						// n precedes m
						// we want m=0, and n=abs(lag)
						double probNBefore = calcProbRupturedBefore(-lag, 1, state, fromState, chain);
						double freqNBefore = probNBefore*freq;
						Preconditions.checkState((float)freqNBefore <= (float)freq,
								"FreqNBefore > Freq!: "+freqNBefore+">"+freq);
						if (state[m] == 0) {
							freqM += freq;
							freqMN += freqNBefore;
							freqEither += freq;
						} else {
							freqEither += freqNBefore;
						}
						freqN += freqNBefore;
						if (lag_debug && freqNBefore > 0 && state[m] == 0) {
							System.out.println("Lag Debug. freq="+freq+", freqNBefore="+freqNBefore+" freqMThis="+freq);
						}
					} else {
						// m precedes n
						// we want n=0, and m=lag
						double probMBefore = calcProbRupturedBefore(lag, 0, state, fromState, chain);
						double freqMBefore = probMBefore*freq;
						Preconditions.checkState((float)freqMBefore <= (float)freq,
								"FreqMBefore > Freq!: "+freqMBefore+">"+freq);
						if (state[n] == 0) {
							freqN += freq;
							freqMN += freqMBefore;
							freqEither += freq;
						} else {
							freqEither += freqMBefore;
						}
						freqM += freqMBefore;
					}
				}
				//				int mCheckIndex, nCheckIndex;
				//				if (lag == 0) {
				//					mCheckIndex = 0;
				//					nCheckIndex = 0;
				//				} else if (lag < 0) {
				//					// n precedes m
				//					// we want m=0, and n=abs(lag)
				//					mCheckIndex = 0;
				//					nCheckIndex = -lag;
				//				} else {
				//					// lag > 0
				//					// m precedes n
				//					// we want n=0, and m=lag
				//					nCheckIndex = 0;
				//					mCheckIndex = lag;
				//				}

			}

			if (lag_debug && freqN > 0 && freqM > 0) {
				System.out.println("Lag Debug. freqM="+freqM+", freqN="+freqN+" freqMN="+freqMN);
			}

			allFreqMs.add(freqM);
			allFreqNs.add(freqN);
			allFreqMNs.add(freqMN);
			allFreqEithers.add(freqEither);
			allStates.add(fromState);
			allTotOccupancy += possible.getTot();
			allTots.add(tot);

			numSubSums++;

			// convert to probs
			double probM = freqM / tot;
			double probN = freqN / tot;
			double probMN = freqMN / tot;

			if (useIndepProbs) {
				//				int mIndex = indicesKey.indices[0];
				//				int nIndex = indicesKey.indices[1];
				int mIndex = 0;
				int nIndex = 1;
				if (mIndepProbs.containsKey(fromState[mIndex]))
					probM = mIndepProbs.get(fromState[mIndex]);
				else
					probM = 0;
				if (nIndepProbs.containsKey(fromState[nIndex]))
					probN = nIndepProbs.get(fromState[nIndex]);
				else
					probN = 0;
			}

			double synch = probMN/(probM*probN);
			//			double prob_state = tot / totStateCount;

			numPossibleSums++;
			
			boolean inclusionCondition;
			if (weightingScheme == WeightingScheme.TOT_OCCUPANCY)
				inclusionCondition = true;
			else
				inclusionCondition = probM > 0 && probN > 0;

			//			if (!Double.isInfinite(synch) && !Double.isNaN(synch)) {
			//			if (freqEither > 0) {
			//			if (freqM > 0 && freqN > 0) {
			if (inclusionCondition) {
				if (Double.isInfinite(synch) || Double.isNaN(synch))
					synch = 1d;
				totEithers += freqEither;
				numSums++;

				synchs.add(synch);
				freqEithers.add(freqEither);
				synch_indices.add(new IndicesKey(fromState));

				probMs.add(probM);
				probNs.add(probN);
				probMNs.add(probMN);

				usedStates.add(fromState);
				totOccupancy += possible.getTot();

				switch (weightingScheme) {
				case FREQ_EITHERS:
					weightNumerators.add(freqEither);
					weightDenominator += freqEither;
					break;
				case TOT_OCCUPANCY:
					weightNumerators.add(possible.getTot());
					weightDenominator += possible.getTot();
					break;
				case FREQ_MNs:
					weightNumerators.add(freqMN);
					weightDenominator += freqMN;
					break;
				case FREQ_AVG:
					double avg = 0.5*(freqM+freqN);
					weightNumerators.add(avg);
					weightDenominator += avg;
					break;
				case FREQ_PROD:
					double prod;
					if (useIndepProbs) {
						//						int mIndex = indicesKey.indices[0];
						//						int nIndex = indicesKey.indices[1];
						int mIndex = 0;
						int nIndex = 1;
						double myFreqM;
						if (mIndepFreqs.containsKey(fromState[mIndex]))
							myFreqM = mIndepFreqs.get(fromState[mIndex]);
						else
							myFreqM = 0;
						double myFreqN;
						if (nIndepFreqs.containsKey(fromState[nIndex]))
							myFreqN = nIndepFreqs.get(fromState[nIndex]);
						else
							myFreqN = 0;
						prod = myFreqM * myFreqN;
					} else {
						prod = freqM*freqN;
					}
					weightNumerators.add(prod);
					weightDenominator += prod;
					break;
				case LAG_DEPENDENT_EITHER:
					if (lag == 0) {
						weightNumerators.add(freqEither);
						weightDenominator += freqEither;
					} else if (lag < 0) {
						// n precedes m
						// we want m=0, and n=abs(lag)
						weightNumerators.add(freqM);
						weightDenominator += freqM;
					} else {
						// m precedes n
						// we want n=0, and m=lag
						weightNumerators.add(freqN);
						weightDenominator += freqN;
					}
					break;
				case LAG_DEPENDENT_AVG:
					if (lag == 0) {
						avg = 0.5*(freqM+freqN);
						weightNumerators.add(avg);
						weightDenominator += avg;
					} else if (lag < 0) {
						// n precedes m
						// we want m=0, and n=abs(lag)
						weightNumerators.add(freqM);
						weightDenominator += freqM;
					} else {
						// m precedes n
						// we want n=0, and m=lag
						weightNumerators.add(freqN);
						weightDenominator += freqN;
					}
					break;
				case LAG_DEPENDENT_PROD:
					if (lag == 0) {
						prod = freqM*freqN;
						weightNumerators.add(prod);
						weightDenominator += prod;
					} else if (lag < 0) {
						// n precedes m
						// we want m=0, and n=abs(lag)
						weightNumerators.add(freqM);
						weightDenominator += freqM;
					} else {
						// m precedes n
						// we want n=0, and m=lag
						weightNumerators.add(freqN);
						weightDenominator += freqN;
					}
					break;
				default:
					throw new IllegalStateException("Weight not implemented: "+weightingScheme);
				}

				Preconditions.checkState(Doubles.isFinite(synch),
						"Non finite synch param. lag="+lag+", synch="+probMN+"/("+probM+"*"+probN+")="+synch);
			}
		}
		weights = Lists.newArrayList();
		double totWeight = 0;
		for (int i=0; i<synchs.size(); i++) {
			double weight = weightNumerators.get(i)/weightDenominator;
			totWeight += weight;
			weights.add(weight);
			gBar_numerator += synchs.get(i)*weight;
			gBar_denominator += weight;
		}
		Preconditions.checkState(synchs.size() == 0 || (float)totWeight == (float)1f,
				"Weights don't add up: totWeight="+(float)totWeight+", num="+synchs.size()
				+", nums=["+Joiner.on(",").join(weightNumerators)+"], denom="+weightDenominator);
		gBar = gBar_numerator/gBar_denominator;

		// now calculate cov
		//		cov = 0;
		//		for (int i=0; i<synchs.size(); i++) {
		//			double weight = weights.get(i);
		//			int mIndex = usedStates.get(i)[m];
		//			int nIndex = usedStates.get(i)[n];
		//			if (mIndepProbs.containsKey(mIndex))
		//				expectedM = mIndepProbs.get(mIndex);
		//			else
		//				expectedM = 0;
		//			if (nIndepProbs.containsKey(nIndex))
		//				expectedN = nIndepProbs.get(nIndex);
		//			else
		//				expectedN = 0;
		//			cov += weight*(probMs.get(i) - expectedM)*(probNs.get(i) - expectedN);
		//		}
		// newer, still diabled
		//		cov = 0;
		//		for (int i=0; i<allStates.size(); i++) {
		//			double freqM = allFreqMs.get(i);
		//			double freqN = allFreqNs.get(i);
		//			double freqMN = allFreqMNs.get(i);
		//			double tot = allTots.get(i);
		//			int mIndex = allStates.get(i)[m];
		//			int nIndex = allStates.get(i)[n];
		//			double expectedM, expectedN;
		//			if (mIndepProbs.containsKey(mIndex))
		//				expectedM = mIndepProbs.get(mIndex);
		//			else
		//				expectedM = 0;
		//			if (nIndepProbs.containsKey(nIndex))
		//				expectedN = nIndepProbs.get(nIndex);
		//			else
		//				expectedN = 0;
		//			cov += tot*(freqM/tot - expectedM)*(freqN/tot - expectedN);
		//		}
		////		cov /= totEithers;
		//		cov /= allTotOccupancy;
		//		Preconditions.checkState(Doubles.isFinite(cov), "COV isn't finite: "+cov);
	}

	private Map<Integer, Double> calcIndepProbs(int index, int target, boolean isFreq) {
		Map<Integer, Double> freqs = Maps.newHashMap();
		Map<Integer, Double> tots = Maps.newHashMap();

		for (int[] indices : chain.getStateTransitionDataset().getPopulatedIndices()) {
			int myIndex = indices[index];
			double tot, freq;
			if (freqs.containsKey(myIndex)) {
				freq = freqs.get(myIndex);
				tot = tots.get(myIndex);
			} else {
				freq = 0;
				tot = 0;
			}
			PossibleStates poss = chain.getStateTransitionDataset().get(indices);
			tot += poss.getTot();
			for (int[] state : poss.getStates()) {
				if (target == 0) {
					if (state[index] == 0)
						freq += poss.getFrequency(state);
				} else {
					freq += calcProbRupturedBefore(target, index, state, indices, chain)*poss.getFrequency(state);
				}
			}
			freqs.put(myIndex, freq);
			tots.put(myIndex, tot);
		}

		if (isFreq)
			return freqs;

		Map<Integer, Double> probs = Maps.newHashMap();

		for (Integer key : freqs.keySet())
			probs.put(key, freqs.get(key)/tots.get(key));
		return probs;
	}

	/**
	 * Returns the synchronization parameter, gBar in linear space
	 * @return
	 */
	public double getGBar() {
		return gBar;
	}

	public void generatePlots(File synchXYZDir, File synchScatterDir, String name1, String name2) throws IOException {
		if (synchs.size() == 0)
			return;
		double distSpacing = chain.getDistSpacing();

		DefaultXY_DataSet synchFunc = new DefaultXY_DataSet();
		DefaultXY_DataSet contribFunc = new DefaultXY_DataSet();

		double highestContrib = Double.NEGATIVE_INFINITY;
		int contribIndex = -1;

		int plotNum = 1000 / (int)distSpacing;

		EvenlyDiscrXYZ_DataSet gXYZ = new EvenlyDiscrXYZ_DataSet(
				plotNum, plotNum, 0.5*distSpacing, 0.5*distSpacing, distSpacing);

		EvenlyDiscrXYZ_DataSet weightXYZ = new EvenlyDiscrXYZ_DataSet(
				plotNum, plotNum, 0.5*distSpacing, 0.5*distSpacing, distSpacing);

		EvenlyDiscrXYZ_DataSet gWeightedXYZ = new EvenlyDiscrXYZ_DataSet(
				plotNum, plotNum, 0.5*distSpacing, 0.5*distSpacing, distSpacing);

		for (int x=0; x<gXYZ.getNumX(); x++) {
			for (int y=0; y<gXYZ.getNumY(); y++) {
				gXYZ.set(x, y, Double.NaN);
				weightXYZ.set(x, y, Double.NaN);
				gWeightedXYZ.set(x, y, Double.NaN);
			}
		}

		for (int i=0; i<synchs.size(); i++) {
			double synch = synchs.get(i);
			//			double weight = freqEithers.get(i)/totEithers;
			double weight = weights.get(i);

			double rBar_without = (gBar_numerator - synch*weight)/(gBar_denominator - weight);
			double contrib = gBar - rBar_without;

			if (contrib > highestContrib) {
				highestContrib = contrib;
				contribIndex = i;
			}

			int[] indices = synch_indices.get(i).getIndices();
			if (indices[0] < gXYZ.getNumX() && indices[1] < gXYZ.getNumY()) {
				gXYZ.set(indices[0], indices[1], synch);
				weightXYZ.set(indices[0], indices[1], weight);
				if (synch == 0)
					gWeightedXYZ.set(indices[0], indices[1], 1e-14);
				else
					gWeightedXYZ.set(indices[0], indices[1], synch*weight);
			}

			if (synch == 0)
				synch = 1e-14;
			synchFunc.set(synch, weight);
			if (contrib < 1e-14 || weight < 1e-14)
				continue;
			contribFunc.set(synch, contrib);
		}

		CPT weightCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, weightXYZ.getMaxZ());
		weightCPT.setNanColor(Color.WHITE);
		CPT gCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-2, 2);
		gXYZ.log();
		gCPT.setNanColor(Color.WHITE);
		CPT weightedCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(-5, -2);
		gWeightedXYZ.log();
		weightedCPT.setNanColor(Color.WHITE);
		String title = name1+" vs "+name2+" (Ln(Gbar)="+(float)Math.log(getGBar())+")";
		String xAxisLabel = "Years since prev "+name1;
		String yAxisLabel = "Years since prev "+name2;
		//		title = "r: freqMN/(freqM*freqN)";
		XYZPlotSpec gSpec = new XYZPlotSpec(gXYZ, gCPT, title, xAxisLabel, yAxisLabel, "Ln(G)");
		XYZPlotSpec weightSpec = new XYZPlotSpec(weightXYZ, weightCPT, title, xAxisLabel, yAxisLabel, "Weight ("+weightingScheme.name()+")");
		XYZPlotSpec gWeightedSpec = new XYZPlotSpec(gWeightedXYZ, weightedCPT, title, xAxisLabel, yAxisLabel, "Ln(G*Weight)");

		Range xyzXRange = new Range(0d, gXYZ.getMaxX()+0.5*distSpacing);
		Range xyzYRange = new Range(0d, 0.6*gXYZ.getMaxX()+0.5*distSpacing);
		List<Range> xyzXRanges = Lists.newArrayList(xyzXRange);
		List<Range> xyzYRanges = Lists.newArrayList(xyzYRange, xyzYRange, xyzYRange);

		double annX = xyzXRange.getUpperBound()*0.95;
		double annY = xyzYRange.getUpperBound()*0.9;
		Font font = new Font(Font.SERIF, Font.PLAIN, 18);
		String lagAdd = "";
		if (lag != 0)
			lagAdd = " lag="+lag;
		XYTextAnnotation rAnn = new XYTextAnnotation("Ln(G=freqMN/(freqM*freqN))"+lagAdd, annX, annY);
		rAnn.setFont(font);
		rAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		gSpec.setPlotAnnotations(Lists.newArrayList(rAnn));
		XYTextAnnotation weightedAnn = new XYTextAnnotation("Ln(G*weight)"+lagAdd, annX, annY);
		weightedAnn.setFont(font);
		weightedAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		gWeightedSpec.setPlotAnnotations(Lists.newArrayList(weightedAnn));
		XYTextAnnotation weightAnn = new XYTextAnnotation("Weight"+lagAdd, annX, annY);
		weightAnn.setFont(font);
		weightAnn.setTextAnchor(TextAnchor.TOP_RIGHT);
		weightSpec.setPlotAnnotations(Lists.newArrayList(weightAnn));

		if (name2.toLowerCase().contains("garlock")) {
			System.out.println("SWAPPING!!!");
			title = name2+" vs "+name1;
			// swap
			gSpec = swapSpec(gSpec);
			gWeightedSpec = swapSpec(gWeightedSpec);
			weightSpec = swapSpec(weightSpec);
		}

		List<XYZPlotSpec> xyzSpecs = Lists.newArrayList(gSpec, weightSpec, gWeightedSpec);
		XYZGraphPanel xyzGP = new XYZGraphPanel();
		xyzGP.drawPlot(xyzSpecs, false, false, xyzXRanges, xyzYRanges, null);
		xyzGP.getChartPanel().setSize(1000, 1500);
		synch2DPlotFile = new File(synchXYZDir, PeriodicityPlotter.getFileSafeString(name1)
				+"_"+PeriodicityPlotter.getFileSafeString(name2)+".pdf");
		if (synchXYZDir.getName().contains("lag"));
		// lag plot, output a png as well
		xyzGP.saveAsPNG(synch2DPlotFile.getAbsolutePath().replaceAll("pdf", "png"));
		xyzGP.saveAsPDF(synch2DPlotFile.getAbsolutePath());

		if (contribIndex >= 0) {
			System.out.println("Highest contrib of "+(float)highestContrib+" at "+synch_indices.get(contribIndex)+". r="
					+probMNs.get(contribIndex).floatValue()+"/("+probMs.get(contribIndex).floatValue()
					+"*"+probNs.get(contribIndex).floatValue()+")="
					+synchs.get(contribIndex).floatValue()+", Weight(Sij)="+(float)(freqEithers.get(contribIndex)/totEithers));
		}

		if (synchScatterDir == null)
			return;

		File scatterPlotFile = new File(synchScatterDir, "synch_scatter_"
				+PeriodicityPlotter.getFileSafeString(name1)+"_"
				+PeriodicityPlotter.getFileSafeString(name2));

		List<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotSymbol.CROSS, 3f, Color.BLACK));
		title = "Synch Param "+name1+" vs "+name2+": "+gBar;
		xAxisLabel = "r: freqMN/(freqM*freqN)";
		yAxisLabel = "Weight";
		PlotSpec spec = new PlotSpec(asList(synchFunc), chars, title, xAxisLabel, yAxisLabel);

		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(14);
		gp.setAxisLabelFontSize(16);
		gp.setPlotLabelFontSize(18);

		Range xRange = new Range(1e-15, 1e3);
		Range yRange = new Range(1e-6, 1e-0);
		gp.drawGraphPanel(spec, true, true, xRange, yRange);

		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(scatterPlotFile.getAbsolutePath()+".png");
		gp.saveAsPDF(scatterPlotFile.getAbsolutePath()+".pdf");
		this.scatterPlotFile = new File(scatterPlotFile.getAbsolutePath()+".pdf");

		File contribScatterFile = new File(synchScatterDir, "synch_scatter_contrib_"
				+PeriodicityPlotter.getFileSafeString(name1)+"_"
				+PeriodicityPlotter.getFileSafeString(name2));
		title = "Synch Param "+name1+" vs "+name2+": "+gBar;
		xAxisLabel = "r: freqMN/(freqM*freqN)";
		yAxisLabel = "Contribution: rBar - rBar(without point)";
		spec = new PlotSpec(asList(contribFunc), chars, title, xAxisLabel, yAxisLabel);

		gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(14);
		gp.setAxisLabelFontSize(16);
		gp.setPlotLabelFontSize(18);

		yRange = new Range(1e-6, 1e1);
		gp.drawGraphPanel(spec, true, true, xRange, yRange);

		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(contribScatterFile.getAbsolutePath()+".png");

		if (name1.contains("Mojave") && name2.contains("Coachella")) {
			// write poster images

			// G
			xyzGP = new XYZGraphPanel();
			xyzXRange = xyzYRange;
			gSpec.setTitle("Gain Factor");
			xyzGP.drawPlot(gSpec, false, false, xyzXRange, xyzYRange);
			xyzGP.getChartPanel().setSize(1000, 1000);
			File gPlotFile = new File(synchXYZDir, "gain_"+PeriodicityPlotter.getFileSafeString(name1)
					+"_"+PeriodicityPlotter.getFileSafeString(name2)+".pdf");
			xyzGP.saveAsPDF(gPlotFile.getAbsolutePath());

			// Occupancy
			EvenlyDiscrXYZ_DataSet occFreqXYZ = new EvenlyDiscrXYZ_DataSet(
					100, 100, 0.5*distSpacing, 0.5*distSpacing, distSpacing);
			SparseNDimensionalHashDataset<Double> totalStatesDataset = chain.getTotalStatesDataset();
			for (int[] indices : totalStatesDataset.getPopulatedIndices()) {
				if (indices[m] < occFreqXYZ.getNumX() && indices[n] < occFreqXYZ.getNumY()) {
					Double val = totalStatesDataset.get(indices);
					if (val != null && val > 0)
						occFreqXYZ.set(indices[m], indices[n], occFreqXYZ.get(indices[m], indices[n])+val);
				}
			}
			for (int x=0; x<occFreqXYZ.getNumX(); x++)
				for (int y=0; y<occFreqXYZ.getNumY(); y++)
					if (occFreqXYZ.get(x, y) == 0)
						occFreqXYZ.set(x, y, Double.NaN);

			CPT occCPT = GMT_CPT_Files.MAX_SPECTRUM.instance().rescale(0, occFreqXYZ.getMaxZ());
			occCPT.setNanColor(Color.WHITE);

			XYZPlotSpec occSpec = new XYZPlotSpec(occFreqXYZ, occCPT, "State Occupancy Frequency",
					gSpec.getXAxisLabel(), gSpec.getYAxisLabel(), null);

			xyzGP = new XYZGraphPanel();
			xyzXRange = xyzYRange;
			xyzGP.drawPlot(occSpec, false, false, xyzXRange, xyzYRange);
			xyzGP.getChartPanel().setSize(1000, 1000);
			File freqPlotFile = new File(synchXYZDir, "freq_"+PeriodicityPlotter.getFileSafeString(name1)
					+"_"+PeriodicityPlotter.getFileSafeString(name2)+".pdf");
			xyzGP.saveAsPDF(freqPlotFile.getAbsolutePath());


		}
	}

	public File getSynch2DPlotFile() {
		return synch2DPlotFile;
	}

	public File getScatterPlotFile() {
		return scatterPlotFile;
	}

	private static class SynchRandTask implements Task {

		private List<SimulatorEvent> events;
		private List<RuptureIdentifier> rupIdens;
		private int[] lags;
		private int numTrials;
		private double distSpacing;
		private int nDims;
		private double[][][][] gBars;
		private RandomDistType dist;
		private int t;

		public SynchRandTask(List<SimulatorEvent> events,
				List<RuptureIdentifier> rupIdens, int[] lags, int numTrials,
				double distSpacing, int nDims, double[][][][] gBars,
				RandomDistType dist, int t) {
			this.events = events;
			this.rupIdens = rupIdens;
			this.lags = lags;
			this.numTrials = numTrials;
			this.distSpacing = distSpacing;
			this.nDims = nDims;
			this.gBars = gBars;
			this.dist = dist;
			this.t = t;
		}

		//		private static int numThreads = 1;
		//		private static int catLenMult = 10;

		private static int numThreads = Runtime.getRuntime().availableProcessors();

		@Override
		public void compute() {
			System.out.println("Random trial "+(t+1)+"/"+numTrials);

			EmpiricalMarkovChain chain = createRandomizedChain(events, rupIdens, dist, distSpacing);

			for (int m=0; m<nDims; m++) {
				for (int n=m; n<nDims; n++) {
					//						Map<IndicesKey, List<int[]>> binnedIndices = getBinnedIndices(chain, m, n);
					for (int i=0; i<lags.length; i++) {
						int lag = lags[i];
						SynchParamCalculator calc = new SynchParamCalculator(chain, m, n, lag);

						//						if (m == 3 && n == 4 && t < 10 && lag == 0) {
						//							File subDir = new File("/tmp/synch_xyz_rand");
						//							if (!subDir.exists())
						//								subDir.mkdir();
						//							try {
						//								calc.generatePlots(subDir, null, "Mojave", "Coachella_rand"+t);
						//							} catch (IOException e) {
						//								ExceptionUtils.throwAsRuntimeException(e);
						//							}
						//						}

						gBars[m][n][t][i] = calc.getGBar();
						gBars[n][m][t][i] = calc.getGBar();
					}
				}
			}
		}

	}

	private static int catLenMult = 1;
	public static EmpiricalMarkovChain createRandomizedChain(List<? extends SimulatorEvent> events,
			List<RuptureIdentifier> rupIdens, RandomDistType dist, double distSpacing) {
		List<SimulatorEvent> randEvents = RandomCatalogBuilder.getRandomResampledCatalog(events, rupIdens, dist, true, catLenMult);

		EmpiricalMarkovChain chain = MarkovChainBuilder.build(distSpacing, randEvents, rupIdens);

		return chain;
	}

	public static void writeSynchParamsStdDev(
			File dir, List<SimulatorEvent> events, List<RuptureIdentifier> rupIdens,
			EmpiricalMarkovChain origChain, int[] lags, int numTrials, double distSpacing) throws IOException {
		int nDims = rupIdens.size();

		double[][][][] gBars = new double[nDims][nDims][numTrials][lags.length];
		RandomDistType dist = RandomDistType.ACTUAL;

		List<SynchRandTask> tasks = Lists.newArrayList();
		for (int t=0; t<numTrials; t++)
			tasks.add(new SynchRandTask(events, rupIdens, lags, numTrials, distSpacing,
					nDims, gBars, dist, t));
		Collections.reverse(tasks);
		ThreadedTaskComputer comp = new ThreadedTaskComputer(tasks, false); // false = don't shuffle
		try {
			comp.computeThreaded(SynchRandTask.numThreads);
		} catch (InterruptedException e) {
			ExceptionUtils.throwAsRuntimeException(e);
		}

		doWriteSynchStdDevParams(dir, rupIdens, origChain, lags, numTrials,
				nDims, gBars);

		//		return stdDevs;
	}

	static void doWriteSynchStdDevParams(File dir,
			List<RuptureIdentifier> rupIdens, EmpiricalMarkovChain origChain,
			int[] lags, int numTrials, int nDims,
			double[][][][] trialGBars) throws IOException {
		double[][][] origGBars = new double[nDims][nDims][lags.length];
		for (int l=0; l<lags.length; l++) {
			for (int i=0; i<nDims; i++) {
				for (int j=i+1; j<nDims; j++) {
					double gBar = calcGBar(origChain, i, j, lags[l]);
					origGBars[i][j][l] = gBar;
					origGBars[j][i][l] = gBar;
				}
			}
		}
		doWriteSynchStdDevParams(dir, rupIdens, lags, numTrials, nDims, trialGBars, origGBars);
	}

	static void doWriteSynchStdDevParams(File dir,
			List<RuptureIdentifier> rupIdens, int[] lags, int numTrials, int nDims,
			double[][][][] trialGBars, double[][][] origGBars) throws IOException {
		for (int l=0; l<lags.length; l++) {
			int lag = lags[l];
			File stdDevCSV;
			if (lag == 0) {
				stdDevCSV = new File(dir, "synch_params_std_devs.csv");
			} else {
				File lagDir = new File(dir, "lag_stats");
				if (!lagDir.exists())
					lagDir.mkdir();
				stdDevCSV = new File(lagDir, "synch_params_std_devs_lag_"+lag+".csv");
			}

			File randTrialsDir = new File(dir, "rand_trials");
			if (!randTrialsDir.exists())
				randTrialsDir.mkdir();

			List<List<String>> allValsLines = Lists.newArrayList();
			for (int i=0; i<=numTrials; i++) {
				List<String> line = Lists.newArrayList();
				if (i == 0)
					line.add("Trial");
				else
					line.add(i+"");
				allValsLines.add(line);
			}
			for (int m=0; m<nDims; m++) {
				for (int n=m+1; n<nDims; n++) {
					String name = rupIdens.get(m).getName()+" vs "+rupIdens.get(n).getName();
					allValsLines.get(0).add(name);
					for (int i=0; i<numTrials; i++)
						allValsLines.get(i+1).add(trialGBars[m][n][i][l]+"");
				}
			}
			CSVFile<String> allValsCSV = new CSVFile<String>(allValsLines, true);
			allValsCSV.writeToFile(new File(randTrialsDir, "synch_params_"+numTrials+"_trials_lag"+lag+".csv"));

			double[][] stdDevs = new double[nDims][nDims];
			double[][] means = new double[nDims][nDims];
			double totBias = 0;
			for (int m=0; m<nDims; m++) {
				for (int n=0; n<nDims; n++) {
					double[] lnVals = new double[numTrials];
					for (int t=0; t<numTrials; t++)
						lnVals[t] = Math.log(trialGBars[m][n][t][l]);
					double mean = StatUtils.mean(lnVals);
					double var = StatUtils.variance(lnVals, mean);
					stdDevs[m][n] = Math.sqrt(var);
					means[m][n] = mean;
					if (n > m)
						totBias += mean;
				}
			}

			System.out.println("Total Bias: "+totBias);

			List<String> header = Lists.newArrayList("");
			for (RuptureIdentifier iden : rupIdens)
				header.add(iden.getName());

			CSVFile<String> csv = new CSVFile<String>(false);

			csv.addLine("Std Dev of "+numTrials+" rand realizations (in Ln space)");
			addTableToCSV(csv, header, stdDevs, true);
			csv.addLine("");
			csv.addLine("");
			csv.addLine("Mean of "+numTrials+" rand realizations (in Ln space)");
			addTableToCSV(csv, header, means, true);
			csv.addLine("");
			csv.addLine("");
			csv.addLine("Print ready +/- (in Ln space)");
			csv.addLine(header);

			DecimalFormat df = new DecimalFormat("0.00");
			for (int i=0; i<nDims; i++) {
				List<String> line = Lists.newArrayList();

				line.add(header.get(i+1));
				for (int j=0; j<nDims; j++) {
					if (i == j)
						line.add("");
					else {
						double stdDev = stdDevs[i][j];
						double gBar = origGBars[i][j][l];
						double mean = Math.log(gBar);
						line.add(df.format(mean)+" ± "+df.format(stdDev));
					}
				}

				csv.addLine(line);
			}

			csv.writeToFile(stdDevCSV);
		}
	}

	public static double calcGBar(EmpiricalMarkovChain chain, int m, int n, int lag) {
		SynchParamCalculator calc;
		if (doLagByShift && lag != 0) {
			calc = new SynchParamCalculator(chain.getCollapsedChain(m, n).getShiftedChain(0, lag), 0);
		} else {
			calc = new SynchParamCalculator(chain, m, n, lag);
		}
		return calc.getGBar();
	}

	private static void addTableToCSV(CSVFile<String> csv, List<String> header, double[][] table, boolean skipIdentities) {
		int nDims = table.length;
		csv.addLine(header);
		for (int i=0; i<nDims; i++) {
			List<String> line = Lists.newArrayList();

			Preconditions.checkState(table[i].length == nDims, "table not square!");

			line.add(header.get(i+1));
			for (int j=0; j<nDims; j++) {
				if (skipIdentities && i == j)
					line.add("");
				else
					line.add(table[i][j]+"");
			}

			csv.addLine(line);
		}
	}

	private static DecimalFormat twoDeimalPlaces = new DecimalFormat("0.00");

	//	public static void writeSynchParamsTable(File file, List<RuptureIdentifier> idens, MarkovChainBuilder chain) throws IOException {
	public static void writeSynchParamsTable(File file, List<RuptureIdentifier> idens, EmpiricalMarkovChain chain,
			Map<IDPairing, HistogramFunction[]> catDensFuncs, int lagMax) throws IOException {
		int nDims = chain.getNDims();

		//			int lagMax = 20;
		int lags = lagMax*2+1;

		double[][][] params = new double[nDims][nDims][lags];

		//		double totStateCount = 0;
		//		for (int[] indices : totalStatesDataset.getPopulatedIndices())
		//			totStateCount += totalStatesDataset.get(indices);

		List<PlotSpec> lagSpecs = Lists.newArrayList();
		List<PlotCurveCharacterstics> lagChars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
		String lagTitle = "Synchronization Lag Functions";
		String lagXAxisLabel = "Lag (years)";
		String lagYAxisLabel;
		//			if (cov)
		//				lagYAxisLabel = "Covariance";
		//			else
		lagYAxisLabel = "Ln(Synchronization)";

		double lagFuncMin = -lagMax*chain.getDistSpacing();

		Range lagXRange = new Range(lagFuncMin, -lagFuncMin);
		Range lagYRange;
		//			if (cov)
		////				lagYRange = new Range(-0.1, 0.1);
		//				lagYRange = new Range(0d, 600d);
		//			else
		lagYRange = new Range(-2, 2);

		List<File> synch2DPDFs = Lists.newArrayList();
		File synchXYZDir = new File(file.getParentFile(), "synch_xyz_param_plots");
		if (!synchXYZDir.exists())
			synchXYZDir.mkdir();

		List<File> synchScatterPDFs = Lists.newArrayList();
		File synchScatterDir = new File(file.getParentFile(), "synch_scatter_plots");
		if (!synchScatterDir.exists())
			synchScatterDir.mkdir();

		File trialDir = new File(file.getParentFile(), "rand_trials");
		String trialPrefix = null;
		int maxTrialsFound = 0;
		if (trialDir.exists()) {
			for (File f : trialDir.listFiles()) {
				String name = f.getName();
				if (!name.endsWith(".csv") || !name.startsWith("synch_params_"))
					continue;

				// make sure has all lags so only consider the maximum
				if (!name.endsWith("_lag"+lagMax+".csv"))
					continue;

				// parse number of trials
				int trials = Integer.parseInt(name.substring("synch_params_".length(), name.indexOf("_trials")));
				if (trials > maxTrialsFound) {
					maxTrialsFound = trials;
					trialPrefix = name.substring(0, name.indexOf("_lag"));
				}
			}
		}

		Map<IDPairing, UncertainArbDiscDataset> lagUncertainties = null;
		Map<IDPairing, Double> pairingStdDevs = Maps.newHashMap();
		if (trialPrefix != null) {
			System.out.println("Loading "+maxTrialsFound+" trials for each with prefix: "+trialPrefix);
			Map<IDPairing, EvenlyDiscretizedFunc> upperFuncs = Maps.newHashMap();
			Map<IDPairing, EvenlyDiscretizedFunc> lowerFuncs = Maps.newHashMap();
			Map<IDPairing, EvenlyDiscretizedFunc> meanFuncs = Maps.newHashMap();
			int[] myLags = rangeInclusive(-lagMax, lagMax);
			for (int l=0; l<myLags.length; l++) {
				int lag = myLags[l];
				File trialCSVFile = new File(trialDir, trialPrefix+"_lag"+lag+".csv");
				Preconditions.checkState(trialCSVFile.exists(), "Trials CSV doesn't exist: "+trialCSVFile.getAbsolutePath());
				CSVFile<String> trialCSV = CSVFile.readFile(trialCSVFile, true);

				Map<String, Integer> colNamesMap = Maps.newHashMap();
				List<String> header = trialCSV.getLine(0);
				for (int i=0; i<header.size(); i++)
					colNamesMap.put(header.get(i), i);

				for (int m=0; m<nDims; m++) {
					for (int n=m+1; n<nDims; n++) {
						String colName = idens.get(m).getName()+" vs "+idens.get(n).getName();
						Integer col = colNamesMap.get(colName);
						if (col == null) {
							// try the other way around
							colName = idens.get(n).getName()+" vs "+idens.get(m).getName();
							col = colNamesMap.get(colName);
							if (col == null) {
								if (l == 0)
									System.out.println("WARNING: lag file doesn't contain: "+colName);
								continue;
							}
						}
						double[] vals = new double[maxTrialsFound];
						for (int i=0; i<vals.length; i++)
							vals[i] = Double.parseDouble(trialCSV.get(i+1, col));

						// TODO calc from normal dist?
						double mean = StatUtils.mean(vals);
						double upper = StatUtils.percentile(vals, 97.5);
						double lower = StatUtils.percentile(vals, 2.5);
						
						double[] lnVals = new double[vals.length];
						for (int i=0; i<vals.length; i++)
							lnVals[i] = Math.log(vals[i]);
						double stdDev = Math.sqrt(StatUtils.variance(lnVals));

						IDPairing pair = new IDPairing(m, n);
						EvenlyDiscretizedFunc upperFunc;
						EvenlyDiscretizedFunc lowerFunc;
						EvenlyDiscretizedFunc meanFunc;
						if (upperFuncs.containsKey(pair)) {
							upperFunc = upperFuncs.get(pair);
							lowerFunc = lowerFuncs.get(pair);
							meanFunc = meanFuncs.get(pair);
						} else {
							upperFunc = new EvenlyDiscretizedFunc(lagFuncMin, lags, chain.getDistSpacing());
							upperFuncs.put(pair, upperFunc);
							upperFuncs.put(pair.getReversed(), upperFunc);
							lowerFunc = new EvenlyDiscretizedFunc(lagFuncMin, lags, chain.getDistSpacing());
							lowerFuncs.put(pair, lowerFunc);
							lowerFuncs.put(pair.getReversed(), lowerFunc);
							meanFunc = new EvenlyDiscretizedFunc(lagFuncMin, lags, chain.getDistSpacing());
							meanFuncs.put(pair, meanFunc);
							meanFuncs.put(pair.getReversed(), meanFunc);
						}

						upperFunc.set(l, Math.log(upper));
						lowerFunc.set(l, Math.log(lower));
						meanFunc.set(l, Math.log(mean));
						
						if (lag == 0)
							pairingStdDevs.put(pair, stdDev);
					}
				}
			}
			lagUncertainties = Maps.newHashMap();
			for (IDPairing pairing : upperFuncs.keySet()) {
				lagUncertainties.put(pairing, new UncertainArbDiscDataset(
						meanFuncs.get(pairing), lowerFuncs.get(pairing), upperFuncs.get(pairing)));
			}
			File stdDevDir = new File(file.getParentFile(), "synch_std_regression");
			plotStdDevs(stdDevDir, pairingStdDevs, chain);
		}

		File synchLagDir = new File(file.getParentFile(), "synch_lag");
		if (!synchLagDir.exists())
			synchLagDir.mkdir();
		File synchLagIndvDir = new File(synchLagDir, "individual");
		if (!synchLagIndvDir.exists())
			synchLagIndvDir.mkdir();

		for (int m=0; m<nDims; m++) {
			// TODO start at m+1?
			for (int n=m; n<nDims; n++) {
				// first bin by only the indices we care about
				//					Map<IndicesKey, List<int[]>> binnedIndices = getBinnedIndices(chain, m, n);

				String name1 = idens.get(m).getName();
				String name2 = idens.get(n).getName();

				int lagIndex = 0;
				//				EvenlyDiscretizedFunc synchLagFunc = new EvenlyDiscretizedFunc((double)-lagMax, lags, 1d);
				EvenlyDiscretizedFunc synchLagFunc = new EvenlyDiscretizedFunc(lagFuncMin, lags, chain.getDistSpacing());

				for (int lag=-lagMax; lag<=lagMax; lag++) {
					lag_debug = lag == -1 && name1.contains("Coachella") && name2.contains("Jacinto") && false;
					SynchParamCalculator calc;
					if (doLagByShift && lag != 0) {
						calc = new SynchParamCalculator(chain.getCollapsedChain(m, n).getShiftedChain(0, lag), 0);
					} else {
						calc = new SynchParamCalculator(chain, m, n, lag);
					}
					lag_debug = false;

					double gBar = calc.getGBar();
					if (!Doubles.isFinite(gBar))
						gBar = 1d;

					//						if (cov) {
					//							params[m][n][lagIndex] = calc.cov;
					//							params[n][m][lagIndex] = calc.cov;
					//							synchLagFunc.set(lagIndex, calc.cov);
					//						} else {
					params[m][n][lagIndex] = gBar;
					params[n][m][lagIndex] = gBar;
					synchLagFunc.set(lagIndex, gBar);
					//						}

					lagIndex++;

					if (name1.contains("Coachella") && name2.contains("Jacinto")) {
						File coachJacintoLagDir = new File(file.getParentFile(), "synch_lag_xyz_coach_jacinto");
						if (!coachJacintoLagDir.exists())
							coachJacintoLagDir.mkdir();
						String lagName1, lagName2;
						if (lag < 0) {
							lagName1 = name1;
							lagName2 = name2+" "+lag+" before";
						} else if (lag > 0) {
							lagName1 = name1+" "+lag+" before";
							lagName2 = name2;
						} else {
							lagName1 = name1;
							lagName2 = name2;
						}
						lagName1 = lagIndex+" "+lagName1;
						if (lagIndex < 10)
							lagName1 = "0"+lagName1;
						calc.generatePlots(coachJacintoLagDir, null, lagName1, lagName2);
					}

					if (lag != 0)
						continue;

					System.out.println(name1+" vs "+name2+": "+gBar);
					//					System.out.println(numerator+"/"+denominator+" = "+params[n][m][lagIndex-1]);
					//					System.out.println("Sums: "+numSums+"/"+numPossibleSums
					//							+" ("+numSubSums+" sub, "+numSubBails+" bails)");

					if (m == n || weightingScheme == WeightingScheme.CATALOG_G)
						continue;

					calc.generatePlots(synchXYZDir, synchScatterDir, name1, name2);

					synch2DPDFs.add(calc.getSynch2DPlotFile());
					synchScatterPDFs.add(calc.getScatterPlotFile());
				}
				if (m == n)
					continue;

				// lag plot spec
				//				String title = "Synch Param "+name1+" vs "+name2+" ("+nDims+"D): "+params[m][n];
				List<DiscretizedFunc> lagFuncs = Lists.newArrayList();
				List<PlotCurveCharacterstics> myLagChars = Lists.newArrayList(lagChars);
				PlotSpec spec;
				//					if (cov) {
				//						lagFuncs.add(synchLagFunc);
				//					} else {
				EvenlyDiscretizedFunc lnSynchFunc = new EvenlyDiscretizedFunc(
						synchLagFunc.getMinX(), synchLagFunc.size(), synchLagFunc.getDelta());
				for (int i=0; i<synchLagFunc.size(); i++)
					lnSynchFunc.set(i, Math.log(synchLagFunc.getY(i)));
				lagFuncs.add(lnSynchFunc);
				double lnSynchAvg = lnSynchFunc.calcSumOfY_Vals()/lnSynchFunc.size();
				double linearSynchAvg = synchLagFunc.calcSumOfY_Vals()/synchLagFunc.size();
				//					}
				spec = new PlotSpec(lagFuncs, myLagChars, lagTitle, lagXAxisLabel, lagYAxisLabel);

				if (catDensFuncs != null) {
					DiscretizedFunc plotSynchFunc = lagFuncs.get(0);

					HistogramFunction hists[] = catDensFuncs.get(new IDPairing(m, n));
					HistogramFunction corups = hists[0];
					HistogramFunction ccdf = hists[1];
					ArbitrarilyDiscretizedFunc scaledCombined = new ArbitrarilyDiscretizedFunc();
					for (Point2D pt : ccdf)
						scaledCombined.set(pt.getX(), pt.getY());
					// there should only be coruptures at pt 0;
					double numCorups = corups.getY(0d);
					if (ccdf.hasX(0d)) {
						// we have a point at zero in the ccdf, just add in the coruptures
						for (Point2D pt : corups) {
							if (pt.getY() > 0) {
								double x = pt.getX();
								double y = pt.getY();
								int xInd = scaledCombined.getXIndex(x);
								if (xInd >= 0)
									y += scaledCombined.getY(xInd);
								scaledCombined.set(x, y);
							}
						}
					} else {
						// CCDF doesn't have a point at zero but we want one. average in the bin before/after zero and add in corups
						double ccdf_at_zero = ccdf.getClosestYtoX(-0.5*chain.getDistSpacing())+ccdf.getClosestYtoX(0.5*chain.getDistSpacing());
						ccdf_at_zero = ccdf_at_zero*0.5 + numCorups;
						scaledCombined.set(0d, ccdf_at_zero);
					}

					double synchMax = Math.max(Math.abs(plotSynchFunc.getMaxY()), Math.abs(plotSynchFunc.getMinY()));
					boolean logCCDF = false;
					if (logCCDF)
						for (int i=0; i<scaledCombined.size(); i++)
							scaledCombined.set(i, Math.log(scaledCombined.getY(i)));
					// subtract the average
					double ccdfAvg = 0d;
					for (Point2D pt : scaledCombined)
						ccdfAvg += pt.getY();
					ccdfAvg /= scaledCombined.size();
					for (int i=0; i<scaledCombined.size(); i++)
						scaledCombined.set(i, scaledCombined.getY(i)-ccdfAvg);
					// it's now centered about 0
					double ccdfMax = Math.max(Math.abs(scaledCombined.getMaxY()), Math.abs(scaledCombined.getMinY()));
					// scale it to match the synch func
					scaledCombined.scale(synchMax/ccdfMax);
					lagFuncs.add(0, scaledCombined);
					myLagChars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 1f, Color.GRAY));
				}

				if (lagUncertainties != null && lagUncertainties.containsKey(new IDPairing(m, n))) {
					System.out.println("We have an uncertainty function!");
					UncertainArbDiscDataset uncertainFunc = lagUncertainties.get(new IDPairing(m, n));
					//						for (int i=0; i<uncertainFunc.getNum(); i++)
					//							System.out.println("\t"+uncertainFunc.getLowerY(i)
					//									+"\t"+uncertainFunc.getY(i)+"\t"+uncertainFunc.getUpperY(i));
					lagFuncs.add(0, uncertainFunc);
					myLagChars.add(0, new PlotCurveCharacterstics(PlotLineType.SHADED_UNCERTAIN_TRANS, 1f, Color.BLUE));
				}

				double annY = lagYRange.getLowerBound()*0.95;
				double annX = lagXRange.getLowerBound()*0.9;
				Font font = new Font(Font.SERIF, Font.PLAIN, 14);
				XYTextAnnotation leftAnn = new XYTextAnnotation(name2+" BEFORE", annX, annY);
				leftAnn.setFont(font);
				leftAnn.setTextAnchor(TextAnchor.BOTTOM_LEFT);
				XYTextAnnotation rightAnn = new XYTextAnnotation(name2+" AFTER", -annX, annY);
				rightAnn.setFont(font);
				rightAnn.setTextAnchor(TextAnchor.BOTTOM_RIGHT);
				XYTextAnnotation centerAnn = new XYTextAnnotation(name1+" at t=0", 0d, annY*0.75);
				centerAnn.setFont(font);
				centerAnn.setTextAnchor(TextAnchor.BOTTOM_CENTER);
				List<XYTextAnnotation> annotations = Lists.newArrayList(leftAnn, rightAnn, centerAnn);

				// now add min/max/avg annotations
				XYTextAnnotation lnStatsAvg = new XYTextAnnotation("Ln: ["+twoDeimalPlaces.format(lnSynchFunc.getMinY())
						+" "+twoDeimalPlaces.format(lnSynchFunc.getMaxY())
						+"], avg="+twoDeimalPlaces.format(lnSynchAvg), annX, -annY);
				lnStatsAvg.setFont(font);
				lnStatsAvg.setTextAnchor(TextAnchor.TOP_LEFT);
				annotations.add(lnStatsAvg);
				XYTextAnnotation linearStatsAvg = new XYTextAnnotation("Linear: ["+twoDeimalPlaces.format(synchLagFunc.getMinY())
						+" "+twoDeimalPlaces.format(synchLagFunc.getMaxY())
						+"], avg="+twoDeimalPlaces.format(linearSynchAvg), annX, -annY*0.75);
				linearStatsAvg.setFont(font);
				linearStatsAvg.setTextAnchor(TextAnchor.TOP_LEFT);
				annotations.add(linearStatsAvg);

				spec.setPlotAnnotations(annotations);
				lagSpecs.add(spec);

				// now write individual function
				File subLagFile = new File(synchLagIndvDir, PeriodicityPlotter.getFileSafeString(name1)
						+"_"+PeriodicityPlotter.getFileSafeString(name2));
				System.out.println(name1+" vs "+name2+": Lag Func Range: ["+synchLagFunc.getMinY()+", "+synchLagFunc.getMaxY()+"]");

				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(14);
				gp.setAxisLabelFontSize(16);
				gp.setPlotLabelFontSize(18);

				gp.setCombinedOnYAxis(false);

				gp.drawGraphPanel(spec, false, false, lagXRange, lagYRange);

				// 8.5x11
				gp.getChartPanel().setSize(800, 300);
				gp.saveAsPDF(subLagFile.getAbsolutePath()+".pdf");
				gp.saveAsPNG(subLagFile.getAbsolutePath()+".png");
			}
		}

		List<List<PlotSpec>> lagSpecPages = Lists.newArrayList();
		lagSpecPages.add(new ArrayList<PlotSpec>());
		int specsToBin = 4;
		for (int i=0; i<lagSpecs.size(); i++) {
			PlotSpec spec = lagSpecs.get(i);
			List<PlotSpec> specs = lagSpecPages.get(lagSpecPages.size()-1);
			if (specs.size() == specsToBin) {
				specs = Lists.newArrayList();
				lagSpecPages.add(specs);
			}
			specs.add(spec);
		}

		List<File> synchLagFiles = Lists.newArrayList();
		EvenlyDiscretizedFunc blankFunc = new EvenlyDiscretizedFunc(lagFuncMin, lags, chain.getDistSpacing());
		List<PlotCurveCharacterstics> lagBlankChars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.SOLID, 0f, Color.BLACK));
		List<Range> lagXRanges = Lists.newArrayList(lagXRange);
		List<Range> lagYRanges = Lists.newArrayList();
		for (int i=0; i<specsToBin; i++)
			lagYRanges.add(lagYRange);
		for (int i=0; i<lagSpecPages.size(); i++) {
			List<PlotSpec> lagSpecPage = lagSpecPages.get(i);
			while (lagSpecPage.size() < specsToBin) {
				// add blank plots
				lagSpecPage.add(new PlotSpec(asList(blankFunc),
						lagBlankChars, lagTitle, lagXAxisLabel, lagYAxisLabel));
			}

			File synchLagFile = new File(synchLagDir, "synch_lag_page"+i+".pdf");

			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setBackgroundColor(Color.WHITE);
			gp.setTickLabelFontSize(14);
			gp.setAxisLabelFontSize(16);
			gp.setPlotLabelFontSize(18);

			gp.setCombinedOnYAxis(false);

			gp.drawGraphPanel(lagSpecPage, false, false, lagXRanges, lagYRanges);

			// 8.5x11
			gp.getChartPanel().setSize(850, 1100);
			gp.saveAsPDF(synchLagFile.getAbsolutePath());
			synchLagFiles.add(synchLagFile);
		}
		if (!synchLagFiles.isEmpty())
			PeriodicityPlotter.combinePDFs(synchLagFiles, new File(synchLagDir, "synch_lags_"+lagMax+".pdf"));
		if (!synch2DPDFs.isEmpty())
			PeriodicityPlotter.combinePDFs(synch2DPDFs, new File(synchXYZDir, "synch_xyzs.pdf"));
		if (!synchScatterPDFs.isEmpty())
			PeriodicityPlotter.combinePDFs(synchScatterPDFs, new File(synchScatterDir, "synch_scatters.pdf"));

		CSVFile<String> csv = new CSVFile<String>(false);

		List<String> header = Lists.newArrayList("");
		for (RuptureIdentifier iden : idens)
			header.add(iden.getName());

		int lag0Index = lagMax;

		int num;
		if (lagUncertainties == null)
			num = 3;
		else
			num = 5;

		for (int type=0; type<num; type++) {
			double[][] myParams;
			//				if (cov) {
			//					switch (type) {
			//					case 0:
			//						csv.addLine("Covariance");
			//						myParams = new double[nDims][nDims];
			//						for (int m=0; m<nDims; m++)
			//							for (int n=0; n<nDims; n++)
			//								myParams[m][n] = params[m][n][lag0Index];
			//						break;
			//					case 1:
			//						csv.addLine("Log10(Covariance + 1)");
			//						myParams = new double[nDims][nDims];
			//						for (int m=0; m<nDims; m++)
			//							for (int n=0; n<nDims; n++)
			//								myParams[m][n] = Math.log10(params[m][n][lag0Index]+1);
			//						break;
			//					case 2:
			//						csv.addLine("Ln(Covariance + 1)");
			//						myParams = new double[nDims][nDims];
			//						for (int m=0; m<nDims; m++)
			//							for (int n=0; n<nDims; n++)
			//								myParams[m][n] = Math.log(params[m][n][lag0Index]+1);
			//						break;
			//
			//					default:
			//						throw new IllegalStateException();
			//					}
			//				} else {
			switch (type) {
			case 0:
				csv.addLine("Linear");
				myParams = new double[nDims][nDims];
				for (int m=0; m<nDims; m++)
					for (int n=0; n<nDims; n++)
						myParams[m][n] = params[m][n][lag0Index];
				break;
			case 1:
				csv.addLine("Log10");
				myParams = new double[nDims][nDims];
				for (int m=0; m<nDims; m++)
					for (int n=0; n<nDims; n++)
						myParams[m][n] = Math.log10(params[m][n][lag0Index]);
				break;
			case 2:
				csv.addLine("Ln");
				myParams = new double[nDims][nDims];
				for (int m=0; m<nDims; m++)
					for (int n=0; n<nDims; n++)
						myParams[m][n] = Math.log(params[m][n][lag0Index]);
				break;
			case 3:
				csv.addLine("Ln Upper 95%");
				myParams = new double[nDims][nDims];
				for (int m=0; m<nDims; m++) {
					for (int n=m+1; n<nDims; n++) {
						IDPairing pair = new IDPairing(m, n);
						UncertainArbDiscDataset func = lagUncertainties.get(pair);
						double val;
						if (func == null)
							val = Double.NaN;
						else
							val = func.getUpperY(0d);
						myParams[m][n] = val;
						myParams[n][m] = val;
					}
				}
				break;
			case 4:
				csv.addLine("Ln Lower 95%");
				myParams = new double[nDims][nDims];
				for (int m=0; m<nDims; m++) {
					for (int n=m+1; n<nDims; n++) {
						IDPairing pair = new IDPairing(m, n);
						UncertainArbDiscDataset func = lagUncertainties.get(pair);
						double val;
						if (func == null)
							val = Double.NaN;
						else
							val = func.getLowerY(0d);
						myParams[m][n] = val;
						myParams[n][m] = val;
					}
				}
				break;

			default:
				throw new IllegalStateException();
			}
			//				}
			addTableToCSV(csv, header, myParams, true);
			csv.addLine("");
			csv.addLine("");
		}

		csv.writeToFile(file);
	}
	
	public static void writeSynchVsProbTable(File file, List<? extends SimulatorEvent> events,
			List<RuptureIdentifier> idens, EmpiricalMarkovChain chain)
					throws IOException {
		CSVFile<String> csv = new CSVFile<String>(true);
		// used to have: , boolean includeCorupInGain
		
		double distSpacing = chain.getDistSpacing();
		double halfSpacing = distSpacing*0.5;
		
		csv.addLine("Fault 1", "Fault 2", "Gbar",
//				"G(1 | 2), "+(int)distSpacing+"yr", "G(2 | 1), "+(int)distSpacing+"yr",
//				"G(1 | 2), "+(int)halfSpacing+"yr", "G(2 | 1), "+(int)halfSpacing+"yr",
//				"G(1 within "+(int)halfSpacing+"yr of 2)", "G(2 within "+(int)halfSpacing+"y of 1)");
				"Catalog G");
		
		for (int m=0; m<idens.size(); m++) {
			RuptureIdentifier iden1 = idens.get(m);
			for (int n=m+1; n<idens.size(); n++) {
				RuptureIdentifier iden2 = idens.get(n);
				
				double gBar = new SynchParamCalculator(chain, m, n, 0).getGBar();
				
//				double gain_10_1_2 = calcProbGain(events, iden2, iden1, distSpacing, includeCorupInGain);
//				double gain_10_2_1 = calcProbGain(events, iden1, iden2, distSpacing, includeCorupInGain);
//				double gain_5_1_2 = calcProbGain(events, iden2, iden1, 0.5*distSpacing, includeCorupInGain);
//				double gain_5_2_1 = calcProbGain(events, iden1, iden2, 0.5*distSpacing, includeCorupInGain);
//				double gain_within5_1_2 = calcProbGainWithinWindow(events, iden2, iden1, 5d);
//				double gain_within5_2_1 = calcProbGainWithinWindow(events, iden1, iden2, 5d);
//				
//				csv.addLine(iden1.getName(), iden2.getName(), (float)gBar+"",
//						(float)gain_10_1_2+"", (float)gain_10_2_1+"",
//						(float)gain_5_1_2+"", (float)gain_5_2_1+"",
//						(float)gain_within5_1_2+"", (float)gain_within5_2_1+"");
				
				double catG = calcCatalogG(chain, m, n);
				csv.addLine(iden1.getName(), iden2.getName(), (float)gBar+"", (float)catG+"");
			}
		}
		
		csv.writeToFile(file);
	}
	
	private static enum RegressionCombType {
		MAX,
		MEAN,
		GEOM_MEAN,
		PROD
	};
	
	private static enum RegressionXUnits {
		RI_MEAN,
		RI_MEDIAN,
		RI_MEAN_FRACT,
		NUM;
	}
	
	private static void plotStdDevs(File outputDir, Map<IDPairing, Double> stdDevs, EmpiricalMarkovChain chain) throws IOException {
		if (!outputDir.exists())
			outputDir.mkdir();
		
		for (RegressionXUnits regX : RegressionXUnits.values()) {
			double[] xVals;
			switch (regX) {
			case RI_MEAN:
				xVals = calcRIs(chain, false);
				break;
			case RI_MEDIAN:
				xVals = calcRIs(chain, true);
				break;
			case NUM:
				xVals = calcEventCounts(chain);
				break;
			case RI_MEAN_FRACT:
				xVals = calcRIs(chain, false);
				double catLen = chain.getFullPath().size()*chain.getDistSpacing();
				for (int i=0; i<xVals.length; i++)
					xVals[i] /= catLen;
				break;

			default:
				throw new IllegalStateException("unknown reg x type: "+regX);
			}
			
			for (RegressionCombType type : RegressionCombType.values()) {
				String outputName = regX.name()+"_COMB_"+type.name();
				
				XY_DataSet xy = new DefaultXY_DataSet();
				SimpleRegression reg = new SimpleRegression(true);
				
				for (IDPairing pair : stdDevs.keySet()) {
					double x1 = xVals[pair.getID1()];
					double x2 = xVals[pair.getID2()];
					double combX;
					
					switch (type) {
					case MAX:
						combX = Math.max(x1, x2);
						break;
					case MEAN:
						combX = 0.5*(x1+x2);
						break;
					case GEOM_MEAN:
						combX = Math.sqrt(x1*x2);
						break;
					case PROD:
						combX = x1*x2;
						break;

					default:
						throw new IllegalStateException("unknown avg type: "+type);
					}
					
					xy.set(combX, stdDevs.get(pair));
					reg.addData(combX, stdDevs.get(pair));
				}
				
				List<PlotElement> elems = Lists.newArrayList();
				List<PlotCurveCharacterstics> chars = Lists.newArrayList();
				List<XYTextAnnotation> anns = Lists.newArrayList();
				
				elems.add(xy);
				chars.add(new PlotCurveCharacterstics(PlotSymbol.BOLD_X, 4f, Color.BLACK));
				
				// add regression curve
				ArbitrarilyDiscretizedFunc regFunc = new ArbitrarilyDiscretizedFunc();
				regFunc.setName("Best Fit Line. MSE="+reg.getMeanSquareError());
				regFunc.set(0d, reg.predict(0d));
				regFunc.set(xy.getMaxX(), reg.predict(xy.getMaxX()));
				
				elems.add(0, regFunc);
				chars.add(0, new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));
				XYTextAnnotation regAnn = new XYTextAnnotation(
						"  MSE="+(float)reg.getMeanSquareError()+", R="+(float)reg.getR(), 0d, 0.12);
				regAnn.setTextAnchor(TextAnchor.BOTTOM_LEFT);
				regAnn.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 20));
				anns.add(regAnn);
				XYTextAnnotation regValAnn = new XYTextAnnotation(
						"  y0="+(float)reg.getIntercept()+", slope="+(float)reg.getSlope(), 0d, 0.11);
				regValAnn.setTextAnchor(TextAnchor.BOTTOM_LEFT);
				regValAnn.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 20));
				anns.add(regValAnn);
				
				PlotSpec spec = new PlotSpec(elems, chars, "Std. Dev. vs "+regX.name(),
						type.name()+" of "+regX.name(), "Synch Std. Dev.");
				spec.setPlotAnnotations(anns);
				HeadlessGraphPanel gp = new HeadlessGraphPanel();
				gp.setBackgroundColor(Color.WHITE);
				gp.setTickLabelFontSize(14);
				gp.setAxisLabelFontSize(16);
				gp.setPlotLabelFontSize(18);

				gp.drawGraphPanel(spec, false, false);

				gp.getChartPanel().setSize(1000, 800);
				gp.saveAsPNG(new File(outputDir, outputName+".png").getAbsolutePath());
//				gp.saveAsPDF(scatterPlotFile.getAbsolutePath()+".pdf");
			}
		}
		
		
	}
	
	private static double[] calcRIs(EmpiricalMarkovChain chain, boolean median) {
		double[] ret = new double[chain.getNDims()];
		for (int m=0; m<chain.getNDims(); m++) {
			int prevIndex = -1;
			List<Double> ois = Lists.newArrayList();
			List<int[]> fullPath = chain.getFullPath();
			for (int i = 0; i < fullPath.size(); i++) {
				int[] state = fullPath.get(i);
				int val = state[m];
				if (val == 0) {
					if (prevIndex >= 0) {
						ois.add((double)(i - prevIndex));
					}
					prevIndex = i;
				}
			}
			double[] oiArray = Doubles.toArray(ois);
			if (median)
				ret[m] = DataUtils.median(oiArray);
			else
				ret[m] = StatUtils.mean(oiArray);
			ret[m] *= chain.getDistSpacing();
		}
		return ret;
	}
	
	private static double[] calcEventCounts(EmpiricalMarkovChain chain) {
		double[] ret = new double[chain.getNDims()];
		for (int[] state : chain.getFullPath()) {
			for (int m=0; m<state.length; m++)
				if (state[m] == 0)
					ret[m]++;
		}
		return ret;
	}
	
	private static double calcEmpiricalCatalogG(EmpiricalMarkovChain chain, int m, int n) {
		int numWindows = 0;
		int numMN = 0;
		int numM = 0;
		int numN = 0;
		
		List<int[]> path = chain.getFullPath();
		
		for (int[] state : path) {
			numWindows++;
			if (state[m] == 0 && state[n] == 0)
				numMN++;
			if (state[m] == 0)
				numM++;
			if (state[n] == 0)
				numN++;
		}
		
		return (double)numWindows * (double)numMN/(double)(numM*numN);
	}
	
	public static double calcCatalogG(MarkovChain chain, int m, int n) {
		double numWindows = 0;
		double numMN = 0;
		double numM = 0;
		double numN = 0;
		
		PossibleStates possibleInitialStates = chain.getOccupancy();
		for (int[] state : possibleInitialStates.getStates()) {
			double freq = possibleInitialStates.getFrequency(state);
			numWindows += freq;
			if (state[m] == 0 && state[n] == 0)
				numMN += freq;
			if (state[m] == 0)
				numM += freq;
			if (state[n] == 0)
				numN += freq;
		}
		
		return numWindows * numMN/(numM*numN);
	}
	
	public static double calcOccG(MarkovChain chain, int m, int n, int lag) {
		// get as 2D chain 
		if (chain.getNDims() != 2 || m != 0 || n != 1)
			chain = chain.getCollapsedChain(m, n);
		PossibleStates occ = chain.getOccupancy();
		PossibleStates occMarginal1 = occ.getMarginal(0);
		PossibleStates occMarginal2 = occ.getMarginal(1);
		
		double probE1 = occMarginal1.getFrequency(new int[] {0}) / occMarginal1.getTot();
		double probE2 = occMarginal2.getFrequency(new int[] {0}) / occMarginal2.getTot();
		
		int[] numeratorState;
		if (lag == 0)
			numeratorState = new int[] {0,0};
		else if (lag < 0)
			numeratorState = new int[] {0,-lag};
		else
			numeratorState = new int[] {lag,0};
		
		double probE1E2 = occ.getFrequency(numeratorState)/occ.getTot();
		
		return probE1E2/(probE1*probE2);
	}
	
	private static double calcProbGain(List<SimulatorEvent> events, RuptureIdentifier first,
			RuptureIdentifier second, double years, boolean includeCorupInGain) {
		List<SimulatorEvent> firstMatches = first.getMatches(events);
		List<SimulatorEvent> secondMatches = second.getMatches(events);
		
		double catStart = events.get(0).getTimeInYears();
		double catEnd = events.get(events.size()-1).getTimeInYears();
		double catLen = catEnd - catStart;
		double annualRate = (double)secondMatches.size()/catLen;
		double binRate = annualRate * years;
		
		int indexInSecond = 0;
		
		double totExpected = 0;
		double totObserved = 0;
		
		for (SimulatorEvent firstE : firstMatches) {
			double binStart = firstE.getTimeInYears();
			double binEnd = binStart + years;
			int numSecondMatches = 0;
			for (int i=indexInSecond; i<secondMatches.size(); i++) {
				SimulatorEvent secondE = secondMatches.get(i);
				double eventYears = secondE.getTimeInYears();
				if (eventYears < binStart) {
					indexInSecond = i;
					continue;
				} else if (eventYears > binEnd) {
					break;
				}
				// this means it's in the window
				if (firstE.getID() == secondE.getID() && !includeCorupInGain)
					continue;
				numSecondMatches++;
			}
			
			totExpected += binRate;
			totObserved += numSecondMatches;
		}
		
		return totObserved / totExpected;
	}
	
	private static double calcProbGainWithinWindow(List<SimulatorEvent> events, RuptureIdentifier first,
			RuptureIdentifier second, double years) {
		List<SimulatorEvent> firstMatches = first.getMatches(events);
		List<SimulatorEvent> secondMatches = second.getMatches(events);
		
		double catStart = events.get(0).getTimeInYears();
		double catEnd = events.get(events.size()-1).getTimeInYears();
		double catLen = catEnd - catStart;
		double annualRate = (double)secondMatches.size()/catLen;
		double binRate = annualRate * years * 2d;
		
		int indexInSecond = 0;
		
		double totExpected = 0;
		double totObserved = 0;
		
		for (SimulatorEvent firstE : firstMatches) {
			double binStart = firstE.getTimeInYears() - years;
			double binEnd = binStart + years;
			int numSecondMatches = 0;
			for (int i=indexInSecond; i<secondMatches.size(); i++) {
				SimulatorEvent secondE = secondMatches.get(i);
				double eventYears = secondE.getTimeInYears();
				if (eventYears < binStart) {
					indexInSecond = i;
					continue;
				} else if (eventYears > binEnd) {
					break;
				}
				// this means it's in the window
				numSecondMatches++;
			}
			
			totExpected += binRate;
			totObserved += numSecondMatches;
		}
		
		return totObserved / totExpected;
	}

	private static XYZPlotSpec swapSpec(XYZPlotSpec spec) {
		EvenlyDiscrXYZ_DataSet orig = (EvenlyDiscrXYZ_DataSet) spec.getXYZ_Data();
		EvenlyDiscrXYZ_DataSet swapped = new EvenlyDiscrXYZ_DataSet(
				orig.getNumX(), orig.getNumY(), orig.getMinX(), orig.getMinY(), orig.getGridSpacingX(), orig.getGridSpacingY());
		for (int x=0; x<orig.getNumX(); x++)
			for (int y=0; y<orig.getNumY(); y++)
				swapped.set(y, x, orig.get(x, y));

		XYZPlotSpec swappedSpec = new XYZPlotSpec(swapped, spec.getCPT(), spec.getTitle(),
				spec.getYAxisLabel(), spec.getXAxisLabel(), spec.getZAxisLabel());
		swappedSpec.setPlotAnnotations(spec.getPlotAnnotations());
		return swappedSpec;
	}

	private static List<PlotElement> asList(PlotElement... elems) {
		return Lists.newArrayList(elems);
	}

	private static double calcProbRupturedBefore(int numStatesBefore, int index, int[] transDestSate, int[] transFromState, EmpiricalMarkovChain chain) {
		int[] toState;
		if (index == 0)
			toState = new int[] {0, -1};
		else if (index == 1)
			toState = new int[] {-1, 0};
		else
			throw new IllegalStateException("Index must be 0 or 1");

		double prob;
		if (transFromState == null)
			prob = chain.getActualTransPathsProbBetweenStates(transDestSate, toState, -numStatesBefore);
		else
			prob = chain.getActualTransPathsProbBetweenStates(transDestSate, toState, -numStatesBefore, transFromState);
		//		double prob = (double)count/timesInState;
		//		Preconditions.checkState(prob <= 1d, "Bad prob: "+prob+", count="+count+", times="+timesInState);
		return prob;
	}

	private static final String getPathStr(List<int[]> path, int[] start) {
		String str = "PATH: ["+start[0]+","+start[1]+"]";
		for (int[] elem : path)
			str += " ["+elem[0]+","+elem[1]+"]";
		return str;
	}

	static int[] rangeInclusive(int min, int max) {
		Preconditions.checkArgument(min <= max);
		int[] ret = new int[max-min+1];
		int cnt = 0;
		for (int i=min; i<=max; i++) {
			ret[cnt++] = i;
		}
		return ret;
	}
	
	public static void plotMarkovNumTransHist(EmpiricalMarkovChain chain, List<RuptureIdentifier> rupIdens, File outputDir)
			throws IOException {
		Preconditions.checkArgument(rupIdens.size() == chain.getNDims());
		
		// count number of events for each index
		int[] eventCounts = new int[chain.getNDims()];
		for (int[] state : chain.getFullPath()) {
			for (int i=0; i<state.length; i++)
				if (state[i] == 0)
					eventCounts[i]++;
		}
		
		// this will sort the indexes by the number of events, decreasing
		List<ComparablePairing<Integer, Integer>> eventComps = Lists.newArrayList();
		for (int i=0; i<eventCounts.length; i++)
			eventComps.add(new ComparablePairing<Integer, Integer>(eventCounts[i], i));
		Collections.sort(eventComps);
		Collections.reverse(eventComps);
		List<Integer> sortedIndexes = Lists.newArrayList();
		for (ComparablePairing<Integer, Integer>  comp : eventComps)
			sortedIndexes.add(comp.getData());
		
		ArbitrarilyDiscretizedFunc occVsFractFunc = new ArbitrarilyDiscretizedFunc();
		
		System.out.println("Fraction of states occupied only once:");
		
		for (int nDims=chain.getNDims(); nDims>=2; nDims--) {
			EmpiricalMarkovChain myChain = chain;
			List<Integer> indexesToInclude = sortedIndexes;
			if (myChain.getNDims() != nDims) {
				indexesToInclude = sortedIndexes.subList(0, nDims);
				myChain = chain.getCollapsedChain(Ints.toArray(indexesToInclude));
			}
			
			List<String> names = Lists.newArrayList();
			for (int i : indexesToInclude)
				names.add(rupIdens.get(i).getName());
			String nameStr = Joiner.on(", ").join(names);
			
			HistogramFunction destHist = new HistogramFunction(0d, 10, 1d);
			HistogramFunction occHist = new HistogramFunction(0d, 10, 1d);
			
			SparseNDimensionalHashDataset<PossibleStates> transData = myChain.getStateTransitionDataset();
			for (int ind[] : transData.getPopulatedIndices()) {
				PossibleStates s = transData.get(ind);
				int numStates = s.getNumStates();
				if (numStates >= destHist.size())
					numStates = destHist.size()-1;
				destHist.add(numStates, 1d);
				Preconditions.checkState(s.getTot() > 0);
				occHist.add(occHist.getClosestXIndex(s.getTot()), 1d);
			}
			
			destHist.normalizeBySumOfY_Vals();
			occHist.normalizeBySumOfY_Vals();
			
			double fractOnlyOnce = occHist.getY(1d);
			occVsFractFunc.set((double)nDims, fractOnlyOnce);
			System.out.println(nDims+"-D: "+(float)fractOnlyOnce);
			
			List<DiscretizedFunc> funcs = Lists.newArrayList();
			List<PlotCurveCharacterstics> chars = Lists.newArrayList();
			
			funcs.add(destHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
			
			Range xRange = new Range(0.5d, destHist.getMaxX()+0.5);
			Range yRange = new Range(0d, 1d);
			Range yLogRange = new Range(1e-3d, 1d);
			
			PlotSpec spec = new PlotSpec(funcs, chars, nDims+"-D State Transition Hist:\n"+nameStr,
					"# Destination States", "Fraction of States");
			
			HeadlessGraphPanel gp = new HeadlessGraphPanel();
			gp.setBackgroundColor(Color.WHITE);
			gp.setTickLabelFontSize(14);
			gp.setAxisLabelFontSize(16);
			gp.setPlotLabelFontSize(14);

			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			File outputFile = new File(outputDir, "state_trans_hist_"+nDims+"D");

			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
			gp.saveAsTXT(outputFile.getAbsolutePath()+".txt");
			
			// now log
			gp.drawGraphPanel(spec, false, true, xRange, yLogRange);
			
			outputFile = new File(outputDir, "state_trans_hist_log_"+nDims+"D");

			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
			
			// now occupancy
			funcs = Lists.newArrayList();
			chars = Lists.newArrayList();
			
			funcs.add(occHist);
			chars.add(new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
			
			spec = new PlotSpec(funcs, chars, nDims+"-D State Occupancy Hist:\n"+nameStr,
					"Occupancy Count", "Fraction of States");
			
			gp = new HeadlessGraphPanel();
			gp.setBackgroundColor(Color.WHITE);
			gp.setTickLabelFontSize(14);
			gp.setAxisLabelFontSize(16);
			gp.setPlotLabelFontSize(14);

			gp.drawGraphPanel(spec, false, false, xRange, yRange);
			
			outputFile = new File(outputDir, "state_occ_hist_"+nDims+"D");

			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
			gp.saveAsTXT(outputFile.getAbsolutePath()+".txt");
			
			// now log
			gp.drawGraphPanel(spec, false, true, xRange, yLogRange);
			
			outputFile = new File(outputDir, "state_occ_hist_log_"+nDims+"D");

			gp.getChartPanel().setSize(1000, 800);
			gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		}
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(occVsFractFunc);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "State Occupancy Hist vs Dimensions",
				"Num Dimensions", "Fract Occupied Only Once");
		
		HeadlessGraphPanel gp = new HeadlessGraphPanel();
		gp.setBackgroundColor(Color.WHITE);
		gp.setTickLabelFontSize(14);
		gp.setAxisLabelFontSize(16);
		gp.setPlotLabelFontSize(14);

		gp.drawGraphPanel(spec, false, false, null, new Range(0d, 1d));
		
		File outputFile = new File(outputDir, "state_occ_vs_dims");

		gp.getChartPanel().setSize(1000, 800);
		gp.saveAsPNG(outputFile.getAbsolutePath()+".png");
		gp.saveAsTXT(outputFile.getAbsolutePath()+".txt");
	}

	public static void main(String[] args) throws IOException {
		List<List<RuptureIdentifier>> setIdens = Lists.newArrayList();
		List<String> setNames = Lists.newArrayList();
		
		// SoCal
		setNames.add("so_cal");
		setIdens.add(SynchIdens.getStandardSoCal());
		
		// NorCal
		setNames.add("nor_cal");
		setIdens.add(SynchIdens.getStandardNorCal());
		
		boolean gen_2d_corr_pdfs = false;
		boolean cov = false;
		double distSpacing = 10d;
		boolean random = false;

		RandomDistType origScrambleDist = null;
		//		RandomDistType origScrambleDist = RandomDistType.ACTUAL;

		RandomDistType randDistType = RandomDistType.STATE_BASED;

		File mainDir;
		if (cov)
			mainDir = new File("/home/kevin/Simulators/synch_cov");
		else
			mainDir = new File("/home/kevin/Simulators/synch");
		if (!mainDir.exists())
			mainDir.mkdir();
		
		List<RuptureIdentifier> allIdens = Lists.newArrayList();
		for (List<RuptureIdentifier> idens : setIdens)
			allIdens.addAll(idens);
		
		List<? extends SimulatorEvent> events = new SimAnalysisCatLoader(true, allIdens, true).getEvents();
		
		for (int s=0; s<setNames.size(); s++) {
			String name = getDirName();

			if (distSpacing != 10d)
				name += "_"+(int)distSpacing+"yr";
			if (random)
				name += "_rand";

			File writeDir = new File(mainDir, name);
			if (!writeDir.exists())
				writeDir.mkdir();
			
			writeDir = new File(writeDir, setNames.get(s));
			if (!writeDir.exists())
				writeDir.mkdir();

			List<RuptureIdentifier> rupIdens = setIdens.get(s);

			List<? extends SimulatorEvent> myEvents;
			if (random)
				myEvents = RandomCatalogBuilder.getRandomResampledCatalog(events, rupIdens, RandomDistType.ACTUAL, true, catLenMult);
			else
				myEvents = events;

			// generate Markov Chain
			EmpiricalMarkovChain chain = MarkovChainBuilder.build(distSpacing, myEvents, rupIdens);
			
			File markovPlotsDir = new File(writeDir, "markov_plots");
			if (!markovPlotsDir.exists())
				markovPlotsDir.mkdir();
			
			plotMarkovNumTransHist(chain, rupIdens, markovPlotsDir);

			// tests
			//		MarkovChainBuilder collapsed = chain.getCollapsedChain(4, 3);
			//		int numStates = 3;
			//		int[] fromState = {0, 0};
			//		int[] toState = {numStates, numStates};
			//		System.out.println("Paths from [0,0] to [3,3]:");
			//		for (List<int[]> path : collapsed.getPathsBetweenStates(fromState, toState, numStates)) {
			//			String str = getPathStr(path, fromState);
			//			System.out.println("\t"+str);
			//		}
			//		System.out.println("Paths from [0,0] to [x,3]:");
			//		for (List<int[]> path : collapsed.getPathsBetweenStates(fromState, new int[] {-1,numStates}, numStates)) {
			//			String str = getPathStr(path, fromState);
			//			System.out.println("\t"+str);
			//		}
			//		double tot = collapsed.getStateTransitionDataset().get(toState).tot;
			//		double freqMBefore = calcFreqRupturedBefore(numStates+10, 0, toState, collapsed);
			//		double freqNBefore = calcFreqRupturedBefore(numStates+10, 1, toState, collapsed);
			//		System.out.println("Freq M before: "+freqMBefore+"/"+tot);
			//		System.out.println("Freq N before: "+freqNBefore+"/"+tot);
			//		System.exit(0);

			// ccdfs/acdfs
			File ccdfDir = new File(writeDir, "ccdfs");
			if (!ccdfDir.exists())
				ccdfDir.mkdir();

			System.out.println("Generating CCDFs...");
			Map<IDPairing, HistogramFunction[]> catDensFuncs =
					PeriodicityPlotter.plotACDF_CCDFs(ccdfDir, myEvents, rupIdens,
							null, null, 2000d, distSpacing);

			int lagMax = 30;

			System.out.println("Calculating Synch Params");
			// write synch CSV
			File synchCSVFile = new File(writeDir, "synch_params.csv");
			writeSynchParamsTable(synchCSVFile, rupIdens, chain, catDensFuncs, lagMax);
			
			writeSynchVsProbTable(new File(writeDir, "synch_compare_prob_gain.csv"),
					myEvents, rupIdens, chain);
//			writeSynchVsProbTable(new File(writeDir, "synch_compare_prob_gain_excl_corup.csv"),
//					myEvents, rupIdens, chain, false);

			// now write std devs
			//		System.out.println("Calculating Synch Std Dev/Biases");
			//		writeSynchParamsStdDev(writeDir, myEvents, rupIdens, chain, new int[] {0}, 100, distSpacing);

			// now do std devs for each lag
			//		System.out.println("Calculating Synch Lag Std Devs/Biases");
			//		writeSynchParamsStdDev(writeDir, myEvents, rupIdens, chain, rangeInclusive(-20, 20), 100, distSpacing);
			
			// write inter event time dists
			List<Color> colors = SynchIdens.getStandardColors();
			File distsDir = new File(writeDir, "inter_event_dists");
			if (!distsDir.exists())
				distsDir.mkdir();
			PeriodicityPlotter.plotPeriodsAndEvents(myEvents, false, false, distsDir,
					rupIdens, colors, false);
		}
	}

	protected static String getDirName() {
		String indepStr;
		if (useIndepProbs)
			indepStr = "indep";
		else
			indepStr = "dep";

		String name = "weight_"+weightingScheme.name()+"_"+indepStr;
		if (doLagByOcc)
			name += "_occLag";
		if (doLagByShift)
			name += "_shiftLag";
		return name;
	}

}
