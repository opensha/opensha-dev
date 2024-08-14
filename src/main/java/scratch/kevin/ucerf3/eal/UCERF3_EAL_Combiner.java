package scratch.kevin.ucerf3.eal;

import java.awt.Color;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DefaultXY_DataSet;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotElement;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.calc.params.MagDistCutoffParam;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sra.calc.parallel.MPJ_CondLossCalc;

import scratch.UCERF3.U3CompoundFaultSystemSolution;
import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.U3FaultSystemSolutionFetcher;
import scratch.UCERF3.erf.mean.TrueMeanBuilder;
import scratch.UCERF3.griddedSeismicity.AbstractGridSourceProvider;
import scratch.UCERF3.logicTree.U3APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.U3BranchWeightProvider;
import scratch.UCERF3.logicTree.U3LogicTreeBranch;
import scratch.UCERF3.logicTree.U3LogicTreeBranchNode;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.UCERF3.utils.MatrixIO;
import scratch.UCERF3.utils.UCERF3_DataUtils;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

/**
 * This class will combine the expected loss for each rupture in the true mean solution to get branch
 * specific EAL values. Similar techniques could be used for a time dependent quick recalculation as well
 * 
 * TODO: add aftershock reduction?
 * TODO: add time dependence
 * @author kevin
 *
 */
public class UCERF3_EAL_Combiner {
	
	private U3FaultSystemSolutionFetcher fetcher;
	private Map<U3LogicTreeBranch, List<Integer>> mappings;
	private U3FaultSystemSolution trueMeanSol;
	private double[][] faultLosses;
	private DiscretizedFunc[] griddedLosses;
	
	private List<U3LogicTreeBranch> branches;
	private double[] faultEALs;
	private double[] griddedEALs;
	private double[] totalEALs;
	
	private double erfProbsDuration;
	private ZipFile erfProbsZipFile;
	
	private DiscretizedFunc lecXVals;
	private LossCOV_Model lecCOV;
	private DiscretizedFunc[] lecs;
	
	public UCERF3_EAL_Combiner(U3FaultSystemSolutionFetcher fetcher, Map<U3LogicTreeBranch, List<Integer>> mappings,
			U3FaultSystemSolution trueMeanSol, double[][] fssLosses, DiscretizedFunc[] griddedLosses)
					throws DocumentException, IOException {
		this(fetcher, mappings, trueMeanSol, fssLosses, griddedLosses, null, Double.NaN);
	}
	
	public UCERF3_EAL_Combiner(U3FaultSystemSolutionFetcher fetcher, Map<U3LogicTreeBranch, List<Integer>> mappings,
			U3FaultSystemSolution trueMeanSol, double[][] fssLosses, DiscretizedFunc[] griddedLosses,
			ZipFile erfProbsZipFile, double erfProbsDuration)
					throws DocumentException, IOException {
		this(fetcher, mappings, trueMeanSol, fssLosses, griddedLosses, erfProbsZipFile, erfProbsDuration,
				null, null);
	}
	
	public UCERF3_EAL_Combiner(U3FaultSystemSolutionFetcher fetcher, Map<U3LogicTreeBranch, List<Integer>> mappings,
			U3FaultSystemSolution trueMeanSol, double[][] fssLosses, DiscretizedFunc[] griddedLosses,
			ZipFile erfProbsZipFile, double erfProbsDuration, DiscretizedFunc lecXVals, LossCOV_Model lecCOV)
					throws DocumentException, IOException {
		this.fetcher = fetcher;
		this.mappings = mappings;
		this.trueMeanSol = trueMeanSol;
		this.faultLosses = fssLosses;
		this.griddedLosses = griddedLosses;
		this.erfProbsZipFile = erfProbsZipFile;
		this.erfProbsDuration = erfProbsDuration;
		this.lecXVals = lecXVals;
		this.lecCOV = lecCOV;
		
		// get list of branches sorted by name
		branches = Lists.newArrayList(mappings.keySet());
		Collections.sort(branches);
		
		calcEALs();
	}
	
	private void calcEALs() throws DocumentException, IOException {
		DiscretizedFunc[] rupMFDs = trueMeanSol.getRupMagDists();
		
		faultEALs = new double[branches.size()];
		griddedEALs = new double[branches.size()];
		totalEALs = new double[branches.size()];
		
		if (lecXVals != null)
			lecs = new DiscretizedFunc[branches.size()];
		
		if (branches.size() > 1)
			System.out.println("calculating branch eals");
		for (int i=0; i<branches.size(); i++) {
			if (i % 100 == 0 && branches.size() > 1)
				System.out.println("Branch "+i);
			U3LogicTreeBranch branch = branches.get(i);
			double[] rates = fetcher.getRates(branch);
			double[] mags = fetcher.getMags(branch);
			List<Integer> meanRupIndexes = mappings.get(branch);
			
			if (erfProbsZipFile != null) {
				// get the rate from the zip file
				String eName = branch.buildFileName()+".bin";
				ZipEntry probsEntry = erfProbsZipFile.getEntry(eName);
				Preconditions.checkNotNull(probsEntry, "Entry not found in zip: "+eName);
				double[] probs = MatrixIO.doubleArrayFromInputStream(
						erfProbsZipFile.getInputStream(probsEntry), probsEntry.getSize());
				Preconditions.checkState(probs.length == rates.length,
						"Prob length mismatch, expected "+rates.length+", got "+probs.length);
				
				rates = new double[probs.length];
				for (int r=0; r<probs.length; r++)
					rates[r] = -Math.log(1 - probs[r])/erfProbsDuration;
			}
			
			CSVFile<String> debugFaultCSV = null;
			DefaultXY_DataSet debugFaultScatter = null;
			if (i == 0 && debug_write) {
				debugFaultCSV = new CSVFile<String>(false);
				debugFaultCSV.addLine("Rup Index", "Mag", "Rate", "Cond. Loss", "Rup EAL");
				debugFaultScatter = new DefaultXY_DataSet();
			}
			
			DiscretizedFunc lec = null;
			if (lecXVals != null) {
				lec = lecXVals.deepClone();
				for (int j=0; j<lec.size(); j++)
					lec.set(j, 0d);
				lecs[i] = lec;
			}
			
			for (int r=0; r<rates.length; r++) {
				int meanRupIndex = meanRupIndexes.get(r);
				double rate = rates[r];
				double mag = mags[r];
				if (rate == 0 || meanRupIndex < 0)
					// skip if rate=0, or if sub seismo
					continue;
//				System.out.println("Rupture "+r+"=>"+meanRupIndex);
				
				// now find the correct index in the rup mfd
				// this is also the rup index in the source
				int rupMFDIndex;
				if (faultLosses[meanRupIndex] == null || faultLosses[meanRupIndex].length == 0)
					continue;
				
				DiscretizedFunc mfd = rupMFDs[meanRupIndex];
				Preconditions.checkState(faultLosses[meanRupIndex].length == mfd.size());
				
				if (faultLosses[meanRupIndex].length == 1) {
					rupMFDIndex = 0;
				} else {
					rupMFDIndex = mfd.getXIndex(mag);
					if (rupMFDIndex < 0) {
						// this is an insertion point, not exact match. find closest
						rupMFDIndex = -(rupMFDIndex+1);
						if (rupMFDIndex > 0 && (float)mfd.getX(rupMFDIndex-1) == (float)mag)
							rupMFDIndex = rupMFDIndex-1;
						else
							Preconditions.checkState(rupMFDIndex < mfd.size() && (float)mfd.getX(rupMFDIndex) == (float)mag,
								"Bad mag. Mine="+mag+". MFD=["+Joiner.on(",").join(mfd.xValues())+"]");
					}
				}
				Preconditions.checkState((float)mag == (float)mfd.getX(rupMFDIndex));
				
				double rupLoss = faultLosses[meanRupIndex][rupMFDIndex];
				
				// TODO aftershock removal, time dependence
				double rupEAL = rupLoss * rate;
				if (debugFaultCSV != null) {
					debugFaultCSV.addLine(r+"", mag+"", rate+"", rupLoss+"", rupEAL+"");
					debugFaultScatter.set(rate, rupEAL);
				}
				faultEALs[i] += rupEAL;
				
				if (lec != null)
					addToLEC(lec, rupLoss, rate);
			}
			
			if (debugFaultCSV != null && plots) {
				debugFaultCSV.writeToFile(new File("/tmp/eals_fault_branch0.csv"));
				writeDebugScatter(debugFaultScatter, "Fault EAL Dist", new File("/tmp/eals_fault_scatter.png"));
			}
			
			// now gridded
			if (griddedLosses != null) {
				CSVFile<String> debugGridCSV = null;
				DefaultXY_DataSet debugGridScatter = null;
				if (i == 0 && debug_write) {
					debugGridCSV = new CSVFile<String>(false);
					debugGridCSV.addLine("Grid Node", "Mag", "Rate", "Cond. Loss", "Rup EAL");
					debugGridScatter = new DefaultXY_DataSet();
				}
				
				GridSourceProvider gridProv = null;
				if (fetcher instanceof U3CompoundFaultSystemSolution) {
					try {
						gridProv = ((U3CompoundFaultSystemSolution)fetcher).loadGridSourceProviderFile(branch);
					} catch (Exception e) {}
				}
				if (gridProv == null)
					gridProv = fetcher.getSolution(branch).getGridSourceProvider();
				for (int n=0; n<gridProv.getGriddedRegion().getNodeCount(); n++) {
					DiscretizedFunc lossDist = griddedLosses[n];
					if (lossDist == null)
						continue;
//					ProbEqkSource source = gridProv.getSource(n, 1d, false, gridType);
//					if (lossDist.getNum() != source.getNumRuptures()) {
//						List<Float> fileMags = Lists.newArrayList();
//						for (double mag : lossDist.xValues())
//							fileMags.add((float)mag);
//						System.out.println("File mags: "+Joiner.on(",").join(fileMags));
//						List<Float> srcMags = Lists.newArrayList();
//						for (ProbEqkRupture rup : source)
//							srcMags.add((float)rup.getMag());
//						System.out.println("Source mags: "+Joiner.on(",").join(srcMags));
//						System.out.flush();
//					}
//					Preconditions.checkState(lossDist.getNum() == source.getNumRuptures(),
//							"Grid source rup count inconsistency. Loaded: "+lossDist.getNum()
//							+", from prov: "+source.getNumRuptures());
					// do mag lookups in floating point precision
					float[] lossMags = new float[lossDist.size()];
					for (int j=0; j<lossDist.size(); j++)
						lossMags[j] = (float)lossDist.getX(j);
					IncrementalMagFreqDist mfd = gridProv.getMFD(n, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF);
					for (int j=0; j<mfd.size(); j++) {
						double mag = mfd.getX(j);
						double rate = mfd.getY(j);
						if (rate == 0d)
							continue;
						int lossIndex = Arrays.binarySearch(lossMags, (float)mag);
						if (lossIndex < 0) {
							System.out.println("Mag: "+mag);
							System.out.println("Rate: "+rate);
							List<Float> fileMags = Lists.newArrayList();
							for (double fmag : lossDist.xValues())
								fileMags.add((float)fmag);
							System.out.println("File mags: "+Joiner.on(",").join(fileMags));
						}
						Preconditions.checkState(lossIndex >= 0, "Loss function doesn't have mag but we do!");
						double loss = lossDist.getY(lossIndex);
//						try {
//							loss = lossDist.getY(mag);
//						} catch (Exception e) {
//							System.out.println("Mag: "+mag);
//							List<Float> fileMags = Lists.newArrayList();
//							for (double fmag : lossDist.xValues())
//								fileMags.add((float)fmag);
//							System.out.println("File mags: "+Joiner.on(",").join(fileMags));
//							throw ExceptionUtils.asRuntimeException(e);
//						}
						if (loss == 0d)
							continue;
						// TODO aftershock removal
						double rupEAL = loss * rate;
						if (debugGridCSV != null) {
							debugGridCSV.addLine(n+"", mag+"", rate+"", loss+"", rupEAL+"");
							debugGridScatter.set(rate, rupEAL);
						}
						griddedEALs[i] += rupEAL;
						if (lec != null)
							addToLEC(lec, loss, rate);
					}
				}
				
				if (debugGridCSV != null) {
					debugGridCSV.writeToFile(new File("/tmp/eals_gridded_branch0.csv"));
					writeDebugScatter(debugGridScatter, "Gridded EAL Dist", new File("/tmp/eals_gridded_scatter.png"));
				}
			}
			totalEALs[i] = faultEALs[i] + griddedEALs[i];
		}
	}
	
	private void addToLEC(DiscretizedFunc lec, double loss, double rate) {
		if (loss == 0d)
			return;
		if (lecCOV == null) {
			for (int k=0; k<lec.size(); k++) {
				if (loss >= lec.getX(k))
					lec.set(k, lec.getY(k)+rate);
				else
					break;
			}
		} else {
			DiscretizedFunc probs = lecCOV.calcLossExceedanceProbs(lec, loss);
			for (int k=0; k<lec.size(); k++)
				lec.set(k, lec.getY(k) + rate*probs.getY(k));
		}
	}
	
	private static void writeDebugScatter(DefaultXY_DataSet dataset, String title, File pngFile) throws IOException {
		List<PlotElement> elems = Lists.newArrayList();
		elems.add(dataset);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 5f, Color.BLACK));
		GraphWindow gw = new GraphWindow(elems, title, chars);
		gw.setX_AxisLabel("Rupture Rate");
		gw.setY_AxisLabel("Rupture EAL");
		gw.saveAsPNG(pngFile.getAbsolutePath());
	}

	public List<U3LogicTreeBranch> getBranches() {
		return branches;
	}

	public double[] getFaultEALs() {
		return faultEALs;
	}

	public double[] getGriddedEALs() {
		return griddedEALs;
	}

	public double[] getTotalEALs() {
		return totalEALs;
	}
	
	public DiscretizedFunc[] getLECs() {
		return lecs;
	}
	
	public EvenlyDiscretizedFunc getTotalEALMagDist(double lossScale) {
		return getTotalEALMagDist(trueMeanSol, faultLosses, griddedLosses, lossScale);
	}
	
	public static EvenlyDiscretizedFunc getTotalEALMagDist(U3FaultSystemSolution trueMeanSol, double[][] faultLosses,
			DiscretizedFunc[] griddedLosses, double lossScale) {
		EvenlyDiscretizedFunc func = new EvenlyDiscretizedFunc(5.05, 41, 0.1);
		
		FaultSystemRupSet meanRupSet = trueMeanSol.getRupSet();
		
		// do fault based
		Preconditions.checkState(faultLosses.length == meanRupSet.getNumRuptures());
		DiscretizedFunc[] rupMFDs = trueMeanSol.getRupMagDists();
		for (int r=0; r<faultLosses.length; r++) {
			DiscretizedFunc mfd = rupMFDs[r];
			Preconditions.checkState(mfd.size() == faultLosses[r].length);
			for (int i=0; i<faultLosses[r].length; i++) {
				double mag = mfd.getX(i);
				double rate = mfd.getY(i);
				double loss = faultLosses[r][i]*lossScale;
				int ind = func.getClosestXIndex(mag);
				func.add(ind, loss*rate);
			}
		}
		
		if (griddedLosses != null) {
			// do gridded
			GridSourceProvider gridProv = trueMeanSol.getGridSourceProvider();
			Preconditions.checkState(gridProv.getNumLocations() == griddedLosses.length,
					"Have %s giddedLosses but grid prov has %s nodes", gridProv.getNumLocations(), griddedLosses.length);
			for (int n=0; n<griddedLosses.length; n++) {
				DiscretizedFunc griddedLossDist = griddedLosses[n];
				IncrementalMagFreqDist mfd = gridProv.getMFD(n, AbstractGridSourceProvider.SOURCE_MIN_MAG_CUTOFF);
				Preconditions.checkState(griddedLossDist.size() == mfd.size());
				for (int i=0; i<griddedLossDist.size(); i++)
					func.add(func.getClosestXIndex(griddedLossDist.getX(i)), griddedLossDist.getY(i)*lossScale*mfd.getY(i));
			}
		}
		
		return func;
	}
	
	private static Map<U3LogicTreeBranch, Double> loadValidateRuns(File validateDir, String prefix, boolean gridded)
			throws IOException {
		Map<U3LogicTreeBranch, Double> results = Maps.newHashMap();
		for (File file : validateDir.listFiles()) {
			if (file.isDirectory())
				continue;
			String name = file.getName();
			if (!name.endsWith(".txt"))
				continue;
			if (!name.startsWith("FM") || !name.contains(prefix))
				continue;
			if (gridded && !name.contains("_bgONLY"))
				continue;
			else if (!gridded && !name.contains("_bgEXCLUDE"))
				continue;
			U3LogicTreeBranch branch = U3LogicTreeBranch.fromFileName(name);
			List<String> lines = FileUtils.readLines(file);
			String[] split = lines.get(0).split(" ");
			double eal = Double.parseDouble(split[split.length-1]);
			results.put(branch, eal);
		}
		return results;
	}
	
	private static final boolean debug_write = false;
	private static final boolean plots = false;

	@SuppressWarnings("unused")
	public static void main(String[] args) throws IOException, DocumentException {
		File invSolDir = new File(UCERF3_DataUtils.DEFAULT_SCRATCH_DATA_DIR, "InversionSolutions");
		File trueMeanSolFile = new File(invSolDir,
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_TRUE_HAZARD_MEAN_SOL_WITH_MAPPING.zip");
		File compoundSolFile = new File(invSolDir,
				"2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_WITH_GRIDDED.zip");

//		File probsZipDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_03_20-ucerf3-erf-probs");
		File probsZipDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_10_07-ucerf3-erf-probs");
		
//		File rupLossesFile = new File("/home/kevin/OpenSHA/UCERF3/eal/2013_10_29-eal/output_fss_index.bin");
//		File rupGriddedFile = new File("/home/kevin/OpenSHA/UCERF3/eal/2013_10_29-eal/output_fss_gridded.bin");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2013_11_05-ucerf3-eal-calc-CB-2013/");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_01_09-ucerf3-eal-calc-CB-2013/");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_01_15-ucerf3-eal-calc-NGA2s-2013");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_03_19-ucerf3-eal-calc-CB2014-recalc");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_04_07-ucerf3-eal-calc-ASK2014-recalc");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_05-ucerf3-eal-calc-wald-vs30");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_16-ucerf3-99percent-wills");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-fatality-smaller");
//		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2014_05_28-ucerf3-99percent-wills-smaller");
		File jobDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2016_06_06-ucerf3-90percent-wald");
		MagDependentAperiodicityOptions[] covs = { null, MagDependentAperiodicityOptions.HIGH_VALUES,
				MagDependentAperiodicityOptions.MID_VALUES, MagDependentAperiodicityOptions.LOW_VALUES };
//		String prefix = "CB_2014";
//		String prefix = "ASK_2014";
//		String prefix = "BSSA_2014";
//		String prefix = "CY_2014";
		String prefix = "IDRISS_2014";
		File rupLossesFile = new File(jobDir, prefix+"_fss_index.bin");
		File rupGriddedFile = new File(jobDir, prefix+"_fss_gridded.bin");
//		File validateDir = new File("/home/kevin/OpenSHA/UCERF3/eal/2013_11_05-ucerf3-eal-calc-CB-2013-validate/");
		File validateDir = null;
//		File rupGriddedFile = null;
//		BackgroundRupType gridType = BackgroundRupType.CROSSHAIR;
		boolean isFSSMapped = true; // if false, then organized as erf source/rup. else, fss rup/mag
		U3BranchWeightProvider weightProv = new U3APrioriBranchWeightProvider();
		
		for (MagDependentAperiodicityOptions cov : covs) {
			String covName;
			if (cov == null)
				covName = "POISSON";
			else
				covName = cov.name();
			File csvFile = new File(jobDir, prefix+"_"+covName+"_eals.csv");
			
			ZipFile erfProbsZipFile = new ZipFile(new File(probsZipDir, "probs_1yr_"+covName+".zip"));
			double erfProbsDuration = 1d;
			
			System.out.println("Loading true mean/compound");
			U3FaultSystemSolution trueMeanSol = U3FaultSystemIO.loadSol(trueMeanSolFile);
			// now load in the mappings
			Map<U3LogicTreeBranch, List<Integer>> mappings = TrueMeanBuilder.loadRuptureMappings(trueMeanSolFile);
			DiscretizedFunc[] rupMFDs = trueMeanSol.getRupMagDists();
			U3CompoundFaultSystemSolution cfss = U3CompoundFaultSystemSolution.fromZipFile(compoundSolFile);
			
			// now load in rupture expected losses
			System.out.println("Loading losses");
			double[][] expectedLosses = MPJ_CondLossCalc.loadResults(rupLossesFile);
			if (!isFSSMapped) {
				System.out.println("Remapping losses");
//				FaultSystemSolutionERF erf = new FaultSystemSolutionERF(trueMeanSol);
//				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
//				erf.updateForecast();
//				expectedLosses = MPJ_EAL_Rupcalc.mapResultsToFSS(erf, expectedLosses);
				double[][] mappedLosses = new double[trueMeanSol.getRupSet().getNumRuptures()][];
				int sourceCount = 0;
				for (int i=0; i<mappedLosses.length; i++) {
					double rate = trueMeanSol.getRateForRup(i);
					if (rate > 0)
						mappedLosses[i] = expectedLosses[sourceCount++];
					else
						mappedLosses[i] = new double[0];
				}
				Preconditions.checkState(sourceCount == expectedLosses.length);
				expectedLosses = mappedLosses;
			}
			DiscretizedFunc[] griddedFuncs = null;
			if (rupGriddedFile != null) {
				griddedFuncs = MPJ_CondLossCalc.loadGridSourcesFile(rupGriddedFile,
						trueMeanSol.getGridSourceProvider().getGriddedRegion());
				double totCond = 0;
				int gridNonNull = 0;
				for (DiscretizedFunc func : griddedFuncs) {
					if (func != null) {
						gridNonNull++;
						for (Point2D pt : func)
							totCond += pt.getY();
					}
				}
				System.out.println("Tot grid conditional "+totCond+" ("+gridNonNull+" non null)");
			}
			
			double trueMeanEAL = 0;
			for (int r=0; r<expectedLosses.length; r++) {
				if (expectedLosses[r] == null)
					continue;
				DiscretizedFunc mfd = rupMFDs[r];
				for (int m=0; m<expectedLosses[r].length; m++) {
					trueMeanEAL += mfd.getY(m)*expectedLosses[r][m];
				}
			}
			
			UCERF3_EAL_Combiner comb = new UCERF3_EAL_Combiner(cfss, mappings, trueMeanSol, expectedLosses, griddedFuncs,
					erfProbsZipFile, erfProbsDuration);
			
			double[] eals = comb.getFaultEALs();
			double[] gridEALs = comb.getGriddedEALs();
			
			List<U3LogicTreeBranch> branches = comb.getBranches();
			
			double totWeights = 0d;
			for (int i=0; i<branches.size(); i++)
				totWeights += weightProv.getWeight(branches.get(i));
			System.out.println("Tot weight: "+totWeights);
			double[] branchWeights = new double[branches.size()];
			for (int i=0; i<branches.size(); i++) {
				U3LogicTreeBranch branch = branches.get(i);
				branchWeights[i] = weightProv.getWeight(branch)/totWeights;
			}
			
			System.out.println("'true mean fault eal'="+trueMeanEAL);
			if (plots)
				plotDist(eals, branchWeights, "UCERF3 Fault Based Indep EAL Distribution",
						new File(jobDir, prefix+"_fault_hist.png"));
			if (griddedFuncs != null && plots)
				plotDist(gridEALs, branchWeights, "UCERF3 Gridded Indep EAL Distribution",
						new File(jobDir, prefix+"_gridded_hist.png"));
			double[] combEALs = comb.getTotalEALs();
			for (int i=0; i<eals.length; i++)
				combEALs[i] = eals[i]+gridEALs[i];
			if (plots)
				plotDist(combEALs, branchWeights, "UCERF3 Total Indep EAL Distribution",
						new File(jobDir, prefix+"_total_hist.png"));
			
			CSVFile<String> csv = new CSVFile<String>(true);
			List<String> header = Lists.newArrayList("Index", "Branch Weight", "Total EAL", "Fault EAL", "Gridded EAL");
			for (Class<? extends U3LogicTreeBranchNode<?>> clazz : U3LogicTreeBranch.getLogicTreeNodeClasses())
				header.add(ClassUtils.getClassNameWithoutPackage(clazz));
			csv.addLine(header);
			for (int i=0; i<branches.size(); i++) {
				U3LogicTreeBranch branch = branches.get(i);
				double weight = branchWeights[i];
				List<String> line = Lists.newArrayList();
				line.add(i+"");
				line.add(weight+"");
				line.add(combEALs[i]+"");
				line.add(eals[i]+"");
				line.add(gridEALs[i]+"");
				for (U3LogicTreeBranchNode<?> node : branch)
					line.add(node.getShortName());
				csv.addLine(line);
			}
			csv.writeToFile(csvFile);
			
			if (!plots)
				continue;
			
			if (validateDir != null) {
				Map<U3LogicTreeBranch, Double> gridValidate = loadValidateRuns(validateDir, prefix, true);
				Map<U3LogicTreeBranch, Double> faultValidate = loadValidateRuns(validateDir, prefix, false);

				if (!gridValidate.isEmpty())
					plotValidates(gridValidate, "Gridded Validate", branches, gridEALs,
							new File(validateDir, prefix+"_validate_gridded.png"));
				if (!faultValidate.isEmpty())
					plotValidates(faultValidate, "Fault Validate", branches, eals,
							new File(validateDir, prefix+"_validate_fault.png"));
			}
			
			// now print the mag dist threshold for reference
			ArbitrarilyDiscretizedFunc magThreshFunc = new MagDistCutoffParam().getDefaultValue();
			GraphWindow gw = new GraphWindow(magThreshFunc, "Mag/Distance Function");
			gw.setX_AxisLabel("Max Distance");
			gw.setY_AxisLabel("Magnitude");
		}
	}
	
	private static void plotValidates(Map<U3LogicTreeBranch, Double> validate, String title,
			List<U3LogicTreeBranch> branches, double[] eals, File pngFile) throws IOException {
		double[] refResults = new double[validate.size()];
		double[] testResults = new double[validate.size()];
		
		MinMaxAveTracker track = new MinMaxAveTracker();
		MinMaxAveTracker absTrack = new MinMaxAveTracker();
		int cnt = 0;
		for (U3LogicTreeBranch branch : validate.keySet()) {
			int ind = branches.indexOf(branch);
			Preconditions.checkPositionIndex(ind, eals.length);
			refResults[cnt] = validate.get(branch);
			testResults[cnt] = eals[ind];
			
			double pDiffAbs = DataUtils.getPercentDiff(testResults[cnt], refResults[cnt]);
			double pDiff = (testResults[cnt] - refResults[cnt]) / refResults[cnt] * 100d;
			
			track.addValue(pDiff);
			absTrack.addValue(pDiffAbs);
			
			cnt++;
		}
		
		DefaultXY_DataSet func = new DefaultXY_DataSet(refResults, testResults);
		List<PlotElement> elems = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		double max = Math.max(func.getMaxY(), func.getMaxX());
		double min = Math.min(func.getMinY(), func.getMinX());
		DefaultXY_DataSet refLine = new DefaultXY_DataSet();
		refLine.set(0d, 0d);
		refLine.set(min, min);
		refLine.set(max, max);
		elems.add(refLine);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		elems.add(func);
		chars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 5f, Color.BLACK));
		
		System.out.println(track);
		System.out.println(absTrack);
		
		title += " (avg % diff: "+(float)track.getAverage()+", abs="+(float)absTrack.getAverage()+")";
		
		GraphWindow gw = new GraphWindow(elems, title, chars);
		gw.setX_AxisLabel("Reference (full hazard curve) EAL");
		gw.setY_AxisLabel("New (conditional loss) EAL");
		if (pngFile != null)
			gw.saveAsPNG(pngFile.getAbsolutePath());
	}
	
	private static void plotDist(double[] eals, double[] branchWeights, String title, File pngFile)
			throws IOException {
		double minEAL = StatUtils.min(eals);
		double maxEAL = StatUtils.max(eals);
		System.out.println("min="+minEAL+"\tmax="+maxEAL);
		
//		double delta = 1000000d;
		// we want about 10 bins
		double binLogBase = 10;
		double delta = (maxEAL - minEAL) / 10d;
		System.out.println("Calc delta: "+delta);
		delta = Math.pow(binLogBase, Math.round(Math.log(delta)/Math.log(binLogBase)));
//		double delta = 100d;
		double funcMin = ((int)(minEAL/delta))*delta;
		double funcMax = ((int)(maxEAL/delta)+1)*delta;
		int funcNum = (int)((funcMax - funcMin)/delta);
		System.out.println("Round delta: "+delta+"\tbins="+funcNum);
		
		double calcMeanEAL = 0d;
		
		HistogramFunction func = new HistogramFunction(funcMin, funcNum, delta);
		for (int i=0; i<branchWeights.length; i++) {
			double eal = eals[i];
			double weight = branchWeights[i];
			func.add(eal, weight);
			calcMeanEAL += eal*weight;
		}
		System.out.println("calc mean eal="+calcMeanEAL);
		
		List<PlotElement> elems = Lists.newArrayList();
		elems.add(func);
		List<PlotCurveCharacterstics> chars = Lists.newArrayList(
				new PlotCurveCharacterstics(PlotLineType.HISTOGRAM, 1f, Color.BLACK));
		DefaultXY_DataSet meanXY = new DefaultXY_DataSet();
		meanXY.set(calcMeanEAL, 0d);
		meanXY.set(calcMeanEAL, func.getMaxY()*1.1d);
		elems.add(meanXY);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		
		GraphWindow gw = new GraphWindow(elems, title+" (mean="+(float)calcMeanEAL+")", chars);
		gw.setX_AxisLabel("EAL ($)");
		gw.setY_AxisLabel("Weighted Fraction of Branches");
		if (pngFile != null)
			gw.saveAsPNG(pngFile.getAbsolutePath());
	}

}
