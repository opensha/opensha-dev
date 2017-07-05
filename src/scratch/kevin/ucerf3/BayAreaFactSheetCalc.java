package scratch.kevin.ucerf3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipException;

import org.dom4j.DocumentException;
import org.opensha.commons.calc.FractileCurveCalculator;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.AbstractXY_DataSet;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.function.XY_DataSetList;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.earthquake.param.BPTAveragingTypeOptions;
import org.opensha.sha.earthquake.param.BPTAveragingTypeParam;
import org.opensha.sha.earthquake.param.HistoricOpenIntervalParam;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityOptions;
import org.opensha.sha.earthquake.param.MagDependentAperiodicityParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;

import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;

import scratch.UCERF3.CompoundFaultSystemSolution;
import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.analysis.FaultSysSolutionERF_Calc;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.logicTree.APrioriBranchWeightProvider;
import scratch.UCERF3.logicTree.BranchWeightProvider;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.utils.FaultSystemIO;

public class BayAreaFactSheetCalc {
	
	private static Map<String, List<Integer>> loadFaultParentSectMappings() {
		Map<String, List<Integer>> map = Maps.newHashMap();
		
		map.put("San Andreas", Lists.newArrayList(654, 655, 657, 658));
		map.put("Hayward-Rodgers Creek", Lists.newArrayList(637, 638, 639, 651));
		map.put("Calaveras", Lists.newArrayList(601, 602, 603, 621));
		map.put("Concord", Lists.newArrayList(635, 636, 711, 713, 622, 623, 640, 677));
		map.put("San Gregorio", Lists.newArrayList(660, 661));
		map.put("Maacama", Lists.newArrayList(644));
		
		return map;
	}
	
	private static final double minX = 5.05d;
	private static final double maxX = 9.05d;
	private static final double delta = 0.1d;
	private static final int num = (int) ((maxX - minX) / delta + 1);

	public static void main(String[] args) throws ZipException, IOException, DocumentException {
		File outputDir = new File("/home/kevin/OpenSHA/UCERF3/bay_area_fact_sheet");
		boolean do_branch_averaged = true;
		if (do_branch_averaged)
			outputDir = new File(outputDir, "branch_averaged");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		File logFile = new File(outputDir, "log.txt");
		FileWriter fw = new FileWriter(logFile);
		
		log_println(fw, "Calculation started using "+BayAreaFactSheetCalc.class.getName()+" on "+new Date());
		
		File compoundFile = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip");
		Map<FaultModels, FaultSystemSolution> baSols = null;
		if (do_branch_averaged) {
			baSols = Maps.newHashMap();
			baSols.put(FaultModels.FM3_1, FaultSystemIO.loadSol(new File(
					"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip")));
			baSols.put(FaultModels.FM3_2, FaultSystemIO.loadSol(new File(
					"/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_2_MEAN_BRANCH_AVG_SOL.zip")));
		}
		
		CompoundFaultSystemSolution cfss = CompoundFaultSystemSolution.fromZipFile(compoundFile);
		
		File probsDir = new File("/home/kevin/OpenSHA/UCERF3/time_dep_erf_probs");
		Map<MagDependentAperiodicityOptions, ERF_ProbsZipFileReader> probsMap = Maps.newHashMap();
		probsMap.put(null, new ERF_ProbsZipFileReader(new File(probsDir, "probs_30yr_POISSON.zip")));
		probsMap.put(MagDependentAperiodicityOptions.LOW_VALUES, new ERF_ProbsZipFileReader(
				new File(probsDir, "probs_30yr_"+MagDependentAperiodicityOptions.LOW_VALUES.name()+".zip")));
		probsMap.put(MagDependentAperiodicityOptions.MID_VALUES, new ERF_ProbsZipFileReader(
				new File(probsDir, "probs_30yr_"+MagDependentAperiodicityOptions.MID_VALUES.name()+".zip")));
		probsMap.put(MagDependentAperiodicityOptions.HIGH_VALUES, new ERF_ProbsZipFileReader(
				new File(probsDir, "probs_30yr_"+MagDependentAperiodicityOptions.HIGH_VALUES.name()+".zip")));
		
		Map<String, List<Integer>> faultParentMappings = loadFaultParentSectMappings();
		
		Table<FaultModels, String, HashSet<Integer>> rupturesTable = HashBasedTable.create();
		
		Region reg = new CaliforniaRegions.SF_BOX();
		
		Table<LogicTreeBranch, String, Map<MagDependentAperiodicityOptions, EvenlyDiscretizedFunc>>
			branchFaultProbsTable = HashBasedTable.create();
		
		BranchWeightProvider weightProv = new APrioriBranchWeightProvider();
		
		EvenlyDiscretizedFunc xVals = null;
		
		int branchIndex = 0;
		
		for (LogicTreeBranch branch : cfss.getBranches()) {
			FaultModels fm = branch.getValue(FaultModels.class);
			
			if (!rupturesTable.containsRow(fm)) {
				log_println(fw, "Loading in ruptures for "+fm.getName());
				FaultSystemRupSet rupSet = cfss.getSolution(branch).getRupSet();
				
				Map<Integer, String> parentNamesMap = Maps.newHashMap();
				
				Map<Integer, List<Integer>> subSectsInRegion = Maps.newHashMap();
				for (FaultSectionPrefData sect : rupSet.getFaultSectionDataList()) {
					for (Location loc : sect.getFaultTrace()) {
						if (reg.contains(loc)) {
							List<Integer> subSects = subSectsInRegion.get(sect.getParentSectionId());
							if (subSects == null) {
								subSects = Lists.newArrayList();
								subSectsInRegion.put(sect.getParentSectionId(), subSects);
							}
							subSects.add(sect.getSectionId());
							break;
						}
					}
					if (!parentNamesMap.containsKey(sect.getParentSectionId()))
						parentNamesMap.put(sect.getParentSectionId(), sect.getParentSectionName());
				}
				
				HashSet<Integer> allRuptures = new HashSet<Integer>();
				
				for (String name : faultParentMappings.keySet()) {
					List<Integer> parentIDs = faultParentMappings.get(name);
					List<String> parentNames = Lists.newArrayList();
					
					HashSet<Integer> rupIndexes = new HashSet<Integer>();
					for (int parentID : parentIDs) {
						// get subsects list so that we exclude ruptures outside
						List<Integer> subSects = subSectsInRegion.get(parentID);
						if (subSects == null)
							// not in this Fault Model
							continue;
						
						parentNames.add(parentID+". "+parentNamesMap.get(parentID));
						
						for (int subSect : subSects)
							rupIndexes.addAll(rupSet.getRupturesForSection(subSect));
						
						allRuptures.addAll(rupIndexes);
					}
					
					log_println(fw, "\tLoaded "+rupIndexes.size()+" ruptures for fault "+name
							+":\n\t\t"+Joiner.on("\n\t\t").join(parentNames));
					
					rupturesTable.put(fm, name, rupIndexes);
				}
				
				// now ruptures that are in the region but not on any of those parents
				HashSet<Integer> rupturesExcluded = new HashSet<Integer>();
				List<String> parentNamesExcluded = Lists.newArrayList();
				for (int parentID : subSectsInRegion.keySet()) {
					List<Integer> subSects = subSectsInRegion.get(parentID);
					
					int prevSize = rupturesExcluded.size();
					for (int subSect : subSects)
						for (int rupIndex : rupSet.getRupturesForSection(subSect))
							if (!allRuptures.contains(rupIndex))
								rupturesExcluded.add(rupIndex);
					if (rupturesExcluded.size() > prevSize)
						parentNamesExcluded.add(parentID+". "+parentNamesMap.get(parentID));
				}
				log_println(fw, "\tLoaded "+rupturesExcluded.size()+" ruptures for all other faults:\n\t\t"
						+Joiner.on("\n\t\t").join(parentNamesExcluded));
				
				rupturesTable.put(fm, "Other Faults", rupturesExcluded);
			} else if (do_branch_averaged) {
				// only do each FM once for BA
				continue;
			}
			
			double[] mags = cfss.getMags(branch);
			
			Map<MagDependentAperiodicityOptions, double[]> baProbs = null;
			if (do_branch_averaged) {
				// get rupture probabilities for each rupture for each time dependent branch
				FaultSystemSolution sol = baSols.get(branch.getValue(FaultModels.class));
				mags = sol.getRupSet().getMagForAllRups();
				
				FaultSystemSolutionERF baERF = new FaultSystemSolutionERF(sol);
				baERF.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
				baProbs = Maps.newHashMap();
				
				for (MagDependentAperiodicityOptions cov : probsMap.keySet()) {
					if (cov == null) {
						baERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
					} else {
						baERF.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.U3_BPT);
						baERF.setParameter(MagDependentAperiodicityParam.NAME, cov);
						baERF.setParameter(HistoricOpenIntervalParam.NAME,
								(double)(FaultSystemSolutionERF.START_TIME_DEFAULT-1875));
						baERF.setParameter(BPTAveragingTypeParam.NAME,
								BPTAveragingTypeOptions.AVE_RI_AVE_NORM_TIME_SINCE);
					}
					baERF.getTimeSpan().setDuration(30d);
					baERF.updateForecast();
					
					double[] probs = new double[mags.length];
					for (int i=0; i<probs.length; i++) {
						int srcIndex = baERF.getSrcIndexForFltSysRup(i);
						if (srcIndex < 0) {
							System.out.println("No source for "+branch.getValue(FaultModels.class).name()+" rupture "+i);
							continue;
						}
						probs[i] = baERF.getSource(srcIndex).computeTotalProb();
					}
					baProbs.put(cov, probs);
				}
			}
			
			// now register probabilities
			for (String faultName : rupturesTable.columnKeySet()) {
				HashSet<Integer> ruptures = rupturesTable.get(fm, faultName);
				
				Map<MagDependentAperiodicityOptions, EvenlyDiscretizedFunc> probFuncs = Maps.newHashMap();
				
				for (MagDependentAperiodicityOptions cov : probsMap.keySet()) {
					double[] erfProbs;
					if (do_branch_averaged)
						erfProbs = baProbs.get(cov);
					else
						erfProbs = probsMap.get(cov).getProbabilities(branch);
					Preconditions.checkState(erfProbs.length == mags.length);
					EvenlyDiscretizedFunc probFunc = getProbs(ruptures, mags, erfProbs);
					if (xVals == null)
						xVals = probFunc;
					probFuncs.put(cov, probFunc);
				}
				
				branchFaultProbsTable.put(branch, faultName, probFuncs);
			}
			
			if (branchIndex % 50 == 0)
				System.out.println("Done with branch "+branchIndex);
			
			branchIndex++;
		}
		
		if (do_branch_averaged)
			Preconditions.checkState(branchFaultProbsTable.rowKeySet().size() == 2);
		else
			Preconditions.checkState(branchFaultProbsTable.rowKeySet().size() == cfss.getBranches().size());
		
		// -1 -> mean
		// -2 -> time dep gain
		double[] fractiles = { -1, -2, 0, 0.025, 0.975, 1};
		
		Map<Double, CSVFile<String>> csvs = Maps.newHashMap();
		
		List<String> faultNames = Lists.newArrayList(branchFaultProbsTable.columnKeySet());
		Collections.sort(faultNames);
		
		for (double fractile : fractiles) {
			CSVFile<String> csv = new CSVFile<String>(true);
			
			List<String> header = Lists.newArrayList();
			header.add("Mag");
			for (String faultName : faultNames)
				header.add(faultName);
			csv.addLine(header);
			for (int i=0; i<xVals.size(); i++) {
				List<String> line = Lists.newArrayList();
				line.add((float)xVals.getX(i)+"");
				for (int j=0; j<faultNames.size(); j++)
					line.add("");
				csv.addLine(line);
			}
			
			csvs.put(fractile, csv);
		}
		
		CSVFile<String> combCSV = new CSVFile<String>(true);
		List<String> header = Lists.newArrayList();
		header.add("Mag");
		for (String faultName : faultNames) {
			for (double fractile : fractiles) {
				if (fractile == -1)
					header.add(faultName);
				else if (fractile == -2)
					header.add("TD Gain");
				else if (fractile == 0)
					header.add("Min");
				else if (fractile == 1)
					header.add("Max");
				else
					header.add("p"+(float)(fractile*100d));
			}
		}
		combCSV.addLine(header);
		for (int i=0; i<xVals.size(); i++) {
			List<String> line = Lists.newArrayList();
			line.add((float)xVals.getX(i)+"");
			for (int j=0; j<header.size()-1; j++)
				line.add("");
			combCSV.addLine(line);
		}
		
		for (int i=0; i<faultNames.size(); i++) {
			String fault = faultNames.get(i);
			XY_DataSetList datas = new XY_DataSetList();
			List<Double> weights = Lists.newArrayList();
			
			XY_DataSetList indepDatas = new XY_DataSetList();
			List<Double> indepWeights = Lists.newArrayList();
			
			for (LogicTreeBranch branch : branchFaultProbsTable.rowKeySet()) {
				double branchWeight;
				if (do_branch_averaged)
					branchWeight = 0.5;
				else
					branchWeight = weightProv.getWeight(branch);
				for (MagDependentAperiodicityOptions cov : probsMap.keySet()) {
					if (!branchFaultProbsTable.containsRow(branch))
						continue;
					double covWeight = FaultSystemSolutionERF.getWeightForCOV(cov);
					weights.add(branchWeight*covWeight);
					EvenlyDiscretizedFunc data = branchFaultProbsTable.get(branch, fault).get(cov);
					datas.add(data);
					
					if (cov == null) {
						indepWeights.add(branchWeight);
						indepDatas.add(data);
					}
				}
			}
			
			log_println(fw, "Calculating "+fault+" with "+datas.size()+" UCERF3-TD logic tree branches and "
					+indepDatas.size()+" UCERF3-TI logic tree branches");
			
			FractileCurveCalculator calc = new FractileCurveCalculator(datas, weights);
			FractileCurveCalculator indepCalc = new FractileCurveCalculator(indepDatas, indepWeights);
			
			List<AbstractXY_DataSet> fractileDatas = Lists.newArrayList();
			
			for (double fractile : fractiles) {
				AbstractXY_DataSet data;
				if (fractile == -1) {
					data = calc.getMeanCurve();
				} else if (fractile == -2) {
					// gain
					data = calc.getMeanCurve();
					AbstractXY_DataSet indepData = indepCalc.getMeanCurve();
					for (int n=0; n<data.size(); n++)
						data.set(n, data.getY(n)/indepData.getY(n));
				} else if (fractile == 0) {
					data = calc.getMinimumCurve();
				} else if (fractile == 1) {
					data = calc.getMaximumCurve();
				} else {
					data = calc.getFractile(fractile);
				}
				fractileDatas.add(data);
			}
			
			for (int f=0; f<fractiles.length; f++) {
				double fractile = fractiles[f];
				CSVFile<String> csv = csvs.get(fractile);
				AbstractXY_DataSet data = fractileDatas.get(f);
				
				int col = i+1;
				for (int n=0; n<data.size(); n++) {
					int row = n+1;
					Preconditions.checkState(csv.get(row, col).isEmpty());
					csv.set(row, col, data.getY(n)+"");
				}
				
				// now combined CSV
				col = 1+(i*fractiles.length)+f;
				for (int n=0; n<data.size(); n++) {
					int row = n+1;
					Preconditions.checkState(combCSV.get(row, col).isEmpty());
					combCSV.set(row, col, data.getY(n)+"");
				}
			}
		}
		
		for (double fractile : fractiles) {
			String name = "results_";
			if (fractile == -1)
				name += "mean";
			else if (fractile == -2)
				name += "td_gain";
			else if (fractile == 0)
				name += "min";
			else if (fractile == 1)
				name += "max";
			else
				name += "p"+(float)(fractile*100d);
			
			csvs.get(fractile).writeToFile(new File(outputDir, name+".csv"));
		}
		
		combCSV.writeToFile(new File(outputDir, "results_combined.csv"));
		
		fw.close();
	}
	
	private static EvenlyDiscretizedFunc getProbs(Collection<Integer> rupIndexes, double[] mags, double[] erfProbs) {
		EvenlyDiscretizedFunc probFunc = new EvenlyDiscretizedFunc(minX-delta*0.5, num, delta);
		
		List<List<Double>> probs = Lists.newArrayList();
		for (int i=0; i<probFunc.size(); i++)
			probs.add(new ArrayList<Double>());
		
		for (int rupIndex : rupIndexes) {
			double mag = mags[rupIndex];
			double prob = erfProbs[rupIndex];
			if (prob == 0)
				continue;
			for (int i=0; i<probFunc.size(); i++) {
				if (mag >= probFunc.getX(i)) {
					probs.get(i).add(prob);
				}
			}
		}
		
		for (int i=0; i<probFunc.size(); i++) 
			probFunc.set(i, FaultSysSolutionERF_Calc.calcSummedProbs(probs.get(i)));
		
		return probFunc;
	}
	
	private static void log_println(FileWriter fw, String str) throws IOException {
		fw.write(str+"\n");
		System.out.println(str);
	}
}
