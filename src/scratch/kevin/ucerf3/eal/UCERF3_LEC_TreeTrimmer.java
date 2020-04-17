package scratch.kevin.ucerf3.eal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.LightFixedXFunc;
import org.opensha.commons.util.ClassUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;

import scratch.UCERF3.enumTreeBranches.InversionModels;
import scratch.UCERF3.logicTree.LogicTreeBranch;
import scratch.UCERF3.logicTree.LogicTreeBranchNode;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMM_Epistemic;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_GMMs;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_LogicTreeBranch;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_ProbModels;
import scratch.kevin.ucerf3.eal.branches.U3_EAL_Vs30Model;

public class UCERF3_LEC_TreeTrimmer {

	@SuppressWarnings("rawtypes")
	public static void main(String[] args) throws IOException {
		File inputDir = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "2020_04_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC");
		File outputDir = new File(inputDir, "tree_trimming");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File csvFile = new File(inputDir, "all_branch_lec_results.csv");
		
		Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> branchLECs = loadLECs(csvFile);
		
		File csvFileEAL = new File(inputDir, "all_branch_results.csv");
		
		Table<U3_EAL_LogicTreeBranch, Double, Double> branchEALs = loadEALs(csvFileEAL);
		
		DiscretizedFunc totalLEC = calcBranchAveragedLEC(branchLECs.rowKeySet(), branchLECs);
		totalLEC.setName("Branch Averaged LEC");
		System.out.println(totalLEC);
		
		double[] lecProbLevels = new double[] { Double.NaN, 0.01, 0.004, 0.0025, 0.0018, 0.0004 };
		
		HashSet<LogicTreeBranchNode<?>> allChoices = new HashSet<>();
		for (U3_EAL_LogicTreeBranch branch : branchLECs.rowKeySet())
			for (int i=0; i<branch.size(); i++)
				allChoices.add(branch.getValue(i));
		
		U3_EAL_LogicTreeBranch branch0 = branchLECs.cellSet().iterator().next().getRowKey();
		
		for (double lecProbLevel : lecProbLevels) {
			if (lecProbLevel > 0)
				System.out.println("Doing calculation for p="+(float)lecProbLevel);
			else
				System.out.println("Doing calculation for EAL");
			
			CSVFile<String> csv = new CSVFile<>(true);
			List<String> header = new ArrayList<>();
			header.add("Iteration");
			header.add("Branch Level");
			header.add("Branch Choice");
			header.add("Leaf Count");
			if (lecProbLevel > 0)
				header.add("Loss at P="+(float)lecProbLevel);
			else
				header.add("Expected Annualized Loss");
			header.add("Error WRT Full Model");
			for (int i=0; i<totalLEC.size(); i++) {
				if (i == 0)
					header.add("LEC X="+(float)totalLEC.getX(i));
				else
					header.add((float)totalLEC.getX(i)+"");
			}
			csv.addLine(header);
			
			List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>(branchLECs.rowKeySet());
			
			List<String> line = new ArrayList<>();
			line.add("N/A"); // iteration
			line.add("Full Tree"); // branch level
			line.add("N/A"); // branch choice
			line.add(branches.size()+""); // leaf count
			double totalLossVal;
			if (lecProbLevel > 0)
				totalLossVal = totalLEC.getFirstInterpolatedX(lecProbLevel);
			else
				totalLossVal = calcBranchAveragedEAL(branches, branchEALs);
			Preconditions.checkState(Double.isFinite(totalLossVal));
			System.out.println("Complete model value: "+totalLossVal);
			line.add(totalLossVal+""); // loss at this prob level
			line.add("0.0"); // error WRT full model
			for (int i=0; i<totalLEC.size(); i++)
				line.add(totalLEC.getY(i)+""); // full LEC
			csv.addLine(line);
			
			List<Class<? extends LogicTreeBranchNode<?>>> availableBranchLevels = new ArrayList<>();
			for (int i=0; i<branch0.size(); i++) {
				@SuppressWarnings("unchecked")
				Class<? extends LogicTreeBranchNode<?>> clazz = (Class<? extends LogicTreeBranchNode<?>>)
						LogicTreeBranch.getEnumEnclosingClass(branch0.getValue(i).getClass());
				int numNonZero = 0;
				for (LogicTreeBranchNode<?> value : clazz.getEnumConstants())
					if (allChoices.contains(value))
						numNonZero++;
				if (numNonZero > 1)
					availableBranchLevels.add(clazz);
			}
			
			int iteration = 0;
			while (!availableBranchLevels.isEmpty()) {
				System.out.println("Iteration "+iteration+" with "
						+availableBranchLevels.size()+" available branch levels and "+branches.size()+" branches");
				
				double minError = Double.POSITIVE_INFINITY;
				int minErrorLevelIndex = -1;
				double minErrorLoss = Double.NaN;
				LogicTreeBranchNode minErrorChoice = null;
				List<U3_EAL_LogicTreeBranch> minErrorBranches = null;
				DiscretizedFunc minErrorLEC = null;
				
				for (int i=0; i<availableBranchLevels.size(); i++) {
					Class<? extends LogicTreeBranchNode<?>> levelClass = availableBranchLevels.get(i);
					LogicTreeBranchNode[] values = levelClass.getEnumConstants();
					
					System.out.println("Testing level: "+ClassUtils.getClassNameWithoutPackage(levelClass));
					
					for (LogicTreeBranchNode value : values) {
						List<U3_EAL_LogicTreeBranch> subBranches = new ArrayList<>();
						for (U3_EAL_LogicTreeBranch branch : branches) {
							if (branch.getValueUnchecked(levelClass) == value)
								subBranches.add(branch);
						}
						System.out.println("\t"+value.getName());
						System.out.println("\t\t"+subBranches.size()+" branches");
						if (!subBranches.isEmpty()) {
							DiscretizedFunc subLEC = calcBranchAveragedLEC(subBranches, branchLECs);
							double subLossVal;
							if (lecProbLevel > 0)
								subLossVal = subLEC.getFirstInterpolatedX(lecProbLevel);
							else
								subLossVal = calcBranchAveragedEAL(subBranches, branchEALs);
							System.out.println("\t\tLoss: "+(float)+subLossVal);
							double error = Math.abs(subLossVal-totalLossVal)/totalLossVal;
							System.out.println("\t\tError: "+(float)+error);
							if (error < minError) {
								minError = error;
								minErrorLevelIndex = i;
								minErrorLoss = subLossVal;
								minErrorChoice = value;
								minErrorBranches = subBranches;
								minErrorLEC = subLEC;
							}
						}
					}
				}
				System.out.println("Selected level: "+minErrorChoice.getBranchLevelName());
				System.out.println("Selected choice: "+minErrorChoice.getName());
				
				availableBranchLevels.remove(minErrorLevelIndex);
				branches = minErrorBranches;
				
				line = new ArrayList<>();
				line.add(iteration+""); // iteration
				line.add(minErrorChoice.getBranchLevelName()); // branch level
				line.add(minErrorChoice.getName()); // branch choice
				line.add(branches.size()+""); // leaf count
				line.add(minErrorLoss+""); // loss at this prob level
				line.add(minError+""); // error WRT full model
				for (int i=0; i<minErrorLEC.size(); i++)
					line.add(minErrorLEC.getY(i)+""); // full LEC
				csv.addLine(line);
				iteration++;
			}
			
			if (lecProbLevel > 0)
				csv.writeToFile(new File(outputDir, "lec_trim_p"+(float)lecProbLevel+".csv"));
			else
				csv.writeToFile(new File(outputDir, "eal_trim.csv"));
		}
	}
	
	private static final int lec_first_col = 15;
	private static final int weight_col = 1;
	private static final int ti_branch_first_col = 2;
	private static final int td_erf_col = 11;
	private static final int gmm_col = 12;
	private static final int gmm_epi_col = 13;
	private static final int vs30_col = 14;
	
	private static Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> loadLECs(File csvFile)
			throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		List<Double> xVals = new ArrayList<>();
		for (int col=lec_first_col; col<csv.getNumCols(); col++)
			xVals.add(csv.getDouble(0, col));
		System.out.println("Have "+xVals.size()+" xVals");
		double[] xValsArray = Doubles.toArray(xVals);
		
		System.out.println("Loading LECs...");
		
		Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> ret = HashBasedTable.create();
		
		double totWeight = 0d;
		
		int rows = csv.getNumRows();
		List<Class<? extends LogicTreeBranchNode<?>>> classList = LogicTreeBranch.getLogicTreeNodeClasses();
		for (int row=1; row<rows; row++) {
			double[] yVals = new double[xValsArray.length];
			for (int i=0; i<xValsArray.length; i++)
				yVals[i] = csv.getDouble(row, i+lec_first_col);
			double weight = csv.getDouble(row, weight_col);
			List<LogicTreeBranchNode<?>> tiVals = new ArrayList<>();
			for (int i=0; i<classList.size(); i++) {
				String str = csv.get(row, i+ti_branch_first_col);
				LogicTreeBranchNode<?> match = forShortName(classList.get(i), str);
				tiVals.add(match);
			}
			LogicTreeBranch tiBranch = LogicTreeBranch.fromValues(tiVals);
			U3_EAL_ProbModels probModel = forShortName(U3_EAL_ProbModels.class, csv.get(row, td_erf_col));
			U3_EAL_GMMs gmm = forShortName(U3_EAL_GMMs.class, csv.get(row, gmm_col));
			U3_EAL_GMM_Epistemic gmmEpi = forShortName(U3_EAL_GMM_Epistemic.class, csv.get(row, gmm_epi_col));
			U3_EAL_Vs30Model vs30 = forShortName(U3_EAL_Vs30Model.class, csv.get(row, vs30_col));
			U3_EAL_LogicTreeBranch ealBranch = new U3_EAL_LogicTreeBranch(tiBranch, probModel, gmm, gmmEpi, vs30);
			
			ret.put(ealBranch, weight, new LightFixedXFunc(xValsArray, yVals));
			totWeight += weight;
		}
		
		System.out.println("Loaded "+ret.size()+" LECs, total weight: "+totWeight);
		
		return ret;
	}
	
	private static final int eal_col = 2;
	private static final int eal_ti_branch_first_col = 5;
	private static final int eal_td_erf_col = 14;
	private static final int eal_gmm_col = 15;
	private static final int eal_gmm_epi_col = 16;
	private static final int eal_vs30_col = 17;
	
	private static Table<U3_EAL_LogicTreeBranch, Double, Double> loadEALs(File csvFile)
			throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		System.out.println("Loading EALs...");
		
		Table<U3_EAL_LogicTreeBranch, Double, Double> ret = HashBasedTable.create();
		
		double totWeight = 0d;
		
		int rows = csv.getNumRows();
		List<Class<? extends LogicTreeBranchNode<?>>> classList = LogicTreeBranch.getLogicTreeNodeClasses();
		for (int row=1; row<rows; row++) {
			double eal = csv.getDouble(row, eal_col);
			double weight = csv.getDouble(row, weight_col);
			List<LogicTreeBranchNode<?>> tiVals = new ArrayList<>();
			for (int i=0; i<classList.size(); i++) {
				String str = csv.get(row, i+eal_ti_branch_first_col);
				LogicTreeBranchNode<?> match = forShortName(classList.get(i), str);
				tiVals.add(match);
			}
			LogicTreeBranch tiBranch = LogicTreeBranch.fromValues(tiVals);
			U3_EAL_ProbModels probModel = forShortName(U3_EAL_ProbModels.class, csv.get(row, eal_td_erf_col));
			U3_EAL_GMMs gmm = forShortName(U3_EAL_GMMs.class, csv.get(row, eal_gmm_col));
			U3_EAL_GMM_Epistemic gmmEpi = forShortName(U3_EAL_GMM_Epistemic.class, csv.get(row, eal_gmm_epi_col));
			U3_EAL_Vs30Model vs30 = forShortName(U3_EAL_Vs30Model.class, csv.get(row, eal_vs30_col));
			U3_EAL_LogicTreeBranch ealBranch = new U3_EAL_LogicTreeBranch(tiBranch, probModel, gmm, gmmEpi, vs30);
			
			ret.put(ealBranch, weight, eal);
			totWeight += weight;
		}
		
		System.out.println("Loaded "+ret.size()+" LECs, total weight: "+totWeight);
		
		return ret;
	}
	
	private static <E extends LogicTreeBranchNode<?>> E forShortName(Class<E> clazz, String shortName) {
		for (E opt : clazz.getEnumConstants())
			if (opt.getShortName().equals(shortName))
				return opt;
		throw new IllegalStateException("No match found for '"+shortName+"' for class "+clazz.getName());
	}
	
	private static DiscretizedFunc calcBranchAveragedLEC(Collection<U3_EAL_LogicTreeBranch> branches,
			Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> branchLECs) {
		DiscretizedFunc lec = null;
		
		double totWeight = 0d;
		
		for (U3_EAL_LogicTreeBranch branch : branches) {
			Map<Double, DiscretizedFunc> map = branchLECs.row(branch);
			Preconditions.checkState(map.size() == 1);
			Double weight = map.keySet().iterator().next();
			DiscretizedFunc branchLEC = map.get(weight);
			
			if (lec == null) {
				lec = new ArbitrarilyDiscretizedFunc();
				for (int i=0; i<branchLEC.size(); i++)
					lec.set(branchLEC.getX(i), 0d);
			} else {
				Preconditions.checkState(branchLEC.size() == lec.size());
			}
			
			for (int i=0; i<lec.size(); i++)
				 lec.set(i, lec.getY(i) + weight*branchLEC.getY(i));
			totWeight += weight;
		}
		// rescale
		lec.scale(1d/totWeight);
		
		return lec;
	}
	
	private static double calcBranchAveragedEAL(Collection<U3_EAL_LogicTreeBranch> branches,
			Table<U3_EAL_LogicTreeBranch, Double, Double> branchEALs) {
		double totEAL = 0d;
		
		double totWeight = 0d;
		
		for (U3_EAL_LogicTreeBranch branch : branches) {
			Map<Double, Double> map = branchEALs.row(branch);
			Preconditions.checkState(map.size() == 1);
			Double weight = map.keySet().iterator().next();
			Double branchEAL = map.get(weight);
			
			totEAL += weight*branchEAL;
			totWeight += weight;
		}
		// rescale
		return totEAL/totWeight;
	}
}
