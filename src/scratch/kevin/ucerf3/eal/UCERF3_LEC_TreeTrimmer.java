package scratch.kevin.ucerf3.eal;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.MathArrays;
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
//				+ "2020_04_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC");
//				+ "2020_07_08-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC-gmVar");
				+ "2020_09_03-ucerf3-ngaw2-cea-100pct-consolidate-calcLEC-covModel");
		List<LogicTreeBranchNode<?>> fixedBranches = new ArrayList<>();
		fixedBranches.add(U3_EAL_GMMs.BSSA_2014);
		fixedBranches.add(U3_EAL_GMM_Epistemic.NONE);
		
		File outputDir = new File(inputDir, "tree_trimming");
		if (fixedBranches != null && !fixedBranches.isEmpty()) {
			String dirName = outputDir.getName()+"_fixed";
			for (LogicTreeBranchNode<?> fixed : fixedBranches)
				dirName += "_"+fixed.encodeChoiceString();
			outputDir = new File(inputDir, dirName);
		}
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		File csvFile = new File(inputDir, "all_branch_lec_results.csv");
		
		Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> branchLECs = loadLECs(csvFile);
		if (fixedBranches != null)
			pruneForFixed(fixedBranches, branchLECs);
		// double as we'll use it for floating point math later
		double totNumBranches = branchLECs.rowKeySet().size();
		System.out.println("Total branch count: "+(int)totNumBranches);
		
		File csvFileEAL = new File(inputDir, "all_branch_results.csv");
		
		Table<U3_EAL_LogicTreeBranch, Double, Double> branchEALs = loadEALs(csvFileEAL);
		if (fixedBranches != null)
			pruneForFixed(fixedBranches, branchEALs);
		
		DiscretizedFunc totalLEC = calcLossDist(branchLECs.rowKeySet(), branchLECs, Double.NaN).lec;
		totalLEC.setName("Branch Averaged LEC");
//		System.out.println(totalLEC);
		
		double[] lecProbLevels = new double[] { 0.01, 0.004, 0.0025, 0.0018, 0.0004, Double.NaN };
		
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
			header.add("COV");
			header.add("K-S Dn");
			header.add("Mean error");
			header.add("COV error");
			header.add("Passes?");
			for (int i=0; i<totalLEC.size(); i++) {
				if (i == 0)
					header.add("CDF X="+(float)totalLEC.getX(i));
				else
					header.add((float)totalLEC.getX(i)+"");
			}
			if (lecProbLevel > 0d) {
				// also include LEC
				for (int i=0; i<totalLEC.size(); i++) {
					if (i == 0)
						header.add("LEC X="+(float)totalLEC.getX(i));
					else
						header.add((float)totalLEC.getX(i)+"");
				}
			}
			csv.addLine(header);
			
			List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>(branchLECs.rowKeySet());
			
			List<String> line = new ArrayList<>();
			line.add("N/A"); // iteration
			line.add("Full Tree"); // branch level
			line.add("N/A"); // branch choice
			line.add(branches.size()+""); // leaf count
			
			LossDistribution totalDist;
			if (lecProbLevel > 0d)
				totalDist = calcLossDist(branches, branchLECs, lecProbLevel);
			else
				totalDist = calcBranchAveragedEAL(branches, branchEALs, totalLEC);

			System.out.println("Complete model value: "+totalDist.mean+" (sd="+totalDist.sd+")");
			Preconditions.checkState(Double.isFinite(totalDist.mean));
			Preconditions.checkState(Double.isFinite(totalDist.sd));
			line.add(totalDist.mean+""); // loss at this prob level
			line.add((totalDist.sd/totalDist.mean)+""); // COV of loss at this prob level
			line.add("N/A"); // K-S Dn
			line.add("N/A"); // mean error WRT full model
			line.add("N/A"); // cov error WRT full model
			line.add("N/A"); // passes tests
			addFuncYVals(totalDist.cdf, line); // full CDF
			if (lecProbLevel > 0d)
				addFuncYVals(totalDist.lec, line); // full LEC
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
				
				double minDn = Double.POSITIVE_INFINITY;
				double meanErrorAtMinDn = Double.NaN;
				double covErrorAtMinDn = Double.NaN;
				int minDnLevelIndex = -1;
				LossDistribution minDnDist = null;
				LogicTreeBranchNode minDnChoice = null;
				List<U3_EAL_LogicTreeBranch> minDnBranches = null;
				
				for (int i=0; i<availableBranchLevels.size(); i++) {
					Class<? extends LogicTreeBranchNode<?>> levelClass = availableBranchLevels.get(i);
					LogicTreeBranchNode[] values = levelClass.getEnumConstants();
					
//					System.out.println("Testing level: "+ClassUtils.getClassNameWithoutPackage(levelClass));
					
					for (LogicTreeBranchNode value : values) {
						List<U3_EAL_LogicTreeBranch> subBranches = new ArrayList<>();
						for (U3_EAL_LogicTreeBranch branch : branches) {
							if (branch.getValueUnchecked(levelClass) == value)
								subBranches.add(branch);
						}
//						System.out.println("\t"+value.getName());
//						System.out.println("\t\t"+subBranches.size()+" branches");
						if (!subBranches.isEmpty()) {
							LossDistribution subLossDist;
							if (lecProbLevel > 0)
								subLossDist = calcLossDist(subBranches, branchLECs, lecProbLevel);
							else
								subLossDist = calcBranchAveragedEAL(subBranches, branchEALs, totalLEC);
							
							double ksDn = 0d;
							for (int j=0; j<totalDist.cdf.size(); j++)
								ksDn = Math.max(ksDn, Math.abs(totalDist.cdf.getY(j)-subLossDist.cdf.getY(j)));
							double meanError = (subLossDist.mean - totalDist.mean)/totalDist.mean;
//							double cov = subLossDist.sd/subLossDist.mean;
							double cov = subLossDist.cov;
							double covError = (cov - totalDist.cov)/totalDist.cov;
							
//							System.out.println("\t\tLoss: "+(float)+subLossDist.mean);
//							System.out.println("\t\tK-S Dn: "+(float)+ksDn);
//							System.out.println("\t\tMean Error: "+(float)+meanError);
//							System.out.println("\t\tCOV Error: "+(float)+covError);
							if (ksDn < minDn) {
								minDn = ksDn;
								minDnLevelIndex = i;
								minDnDist = subLossDist;
								minDnChoice = value;
								minDnBranches = subBranches;
								meanErrorAtMinDn = meanError;
								covErrorAtMinDn = covError;
							}
						}
					}
				}
				System.out.println("Selected level: "+minDnChoice.getBranchLevelName());
				System.out.println("Selected choice: "+minDnChoice.getName());
				
				availableBranchLevels.remove(minDnLevelIndex);
				branches = minDnBranches;
				
				line = new ArrayList<>();
				line.add(iteration+""); // iteration
				line.add(minDnChoice.getBranchLevelName()); // branch level
				line.add(minDnChoice.getName()); // branch choice
				line.add(branches.size()+""); // leaf count
				line.add(minDnDist.mean+""); // loss at this prob level
				line.add(minDnDist.cov+""); // COV of loss at this prob level
				line.add(minDn+""); // K-S Dn
				line.add(meanErrorAtMinDn+""); // mean error WRT full model
				line.add(covErrorAtMinDn+""); // cov error WRT full model
				double n = ((double)branches.size()*totNumBranches)/((double)branches.size()+totNumBranches);
				System.out.println("\tN="+(float)n+" (leaf count: "+branches.size()+")");
				System.out.println("\tDn: "+(float)minDn);
				double dnThresh = 1.63/Math.sqrt(n);
				System.out.println("\tThreshold: "+dnThresh);
				boolean passes = minDn <= 1.63/Math.sqrt(n);
				System.out.println("\tPasses Dn? "+passes);
				passes = passes && Math.abs(meanErrorAtMinDn) <= 0.05;
				passes = passes && Math.abs(covErrorAtMinDn) <= 0.05;
				System.out.println("\tPasses all? "+passes);
				if (passes)
					line.add("TRUE");
				else
					line.add("FALSE");
				addFuncYVals(minDnDist.cdf, line); // full CDF
				if (lecProbLevel > 0d)
					addFuncYVals(minDnDist.lec, line); // full LEC
				csv.addLine(line);
				iteration++;
			}
			
			if (lecProbLevel > 0)
				csv.writeToFile(new File(outputDir, "lec_trim_p"+(float)lecProbLevel+".csv"));
			else
				csv.writeToFile(new File(outputDir, "eal_trim.csv"));
		}
	}
	
	private static void addFuncYVals(DiscretizedFunc func, List<String> line) {
		for (Point2D pt : func)
			line.add((float)pt.getY()+"");
	}
	
	@SuppressWarnings("unchecked")
	private static void pruneForFixed(List<LogicTreeBranchNode<?>> fixedBranches,
			Table<U3_EAL_LogicTreeBranch, Double, ?> table) {
		List<U3_EAL_LogicTreeBranch> branches = new ArrayList<>(table.rowKeySet());
		for (U3_EAL_LogicTreeBranch branch : branches) {
			for (LogicTreeBranchNode<?> fixed : fixedBranches) {
				Class<? extends LogicTreeBranchNode<?>> clazz =
						(Class<? extends LogicTreeBranchNode<?>>)fixed.getClass();
				if (branch.getValueUnchecked(clazz) != fixed) {
					List<Double> weights = new ArrayList<>(table.row(branch).keySet());
					for (Double weight : weights)
						Preconditions.checkNotNull(table.remove(branch, weight));
					break;
				}
			}
		}
	}
	
	private static final int weight_col = 1;
	private static final int ti_branch_first_col = 2;
	private static final int td_erf_col = 11;
	private static final int gmm_col = 12;
	private static final int gmm_epi_col = 13;
	private static final int vs30_col = 14;
	
	static Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> loadLECs(File csvFile)
			throws IOException {
		CSVFile<String> csv = CSVFile.readFile(csvFile, true);
		
		List<Double> xVals = new ArrayList<>();
		int lecFirstCol = vs30_col + 1;
		for (int col=lecFirstCol; col<csv.getNumCols(); col++)
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
				yVals[i] = csv.getDouble(row, i+lecFirstCol);
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
	
	private static LossDistribution calcLossDist(Collection<U3_EAL_LogicTreeBranch> branches,
			Table<U3_EAL_LogicTreeBranch, Double, DiscretizedFunc> branchLECs, double probLevel) {
		DiscretizedFunc lec = null;
		
		double totWeight = 0d;
		
		double[] values = null;
		double[] weights = null;
		DiscretizedFunc cdf = null;
		if (probLevel > 0d) {
			values = new double[branches.size()];
			weights = new double[branches.size()];
		}
		
		int index = 0;
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
			
			if (probLevel > 0d) {
				try {
					values[index] = branchLEC.getFirstInterpolatedX(probLevel);
				} catch (Exception e) {
					e.printStackTrace();
					System.err.flush();
					System.err.println("LEC:\n"+branchLEC);
					System.exit(1);
				}
				weights[index] = weight;
				if (index == 0) {
					cdf = new ArbitrarilyDiscretizedFunc();
					for (Point2D pt : branchLEC)
						cdf.set(pt.getX(), 0d);
				}
				for (int i=0; i<cdf.size(); i++)
					if (values[index] <= cdf.getX(i))
						cdf.set(i, cdf.getY(i)+weight);
			}
			
			index++;
		}
		// rescale
		lec.scale(1d/totWeight);
		if (cdf != null)
			cdf.scale(1d/totWeight);
		
		return new LossDistribution(values, weights, cdf, lec);
	}
	
	private static class LossDistribution {
		private final double[] lossVals;
		private final double[] lossWeights;
		private final double mean;
		private final double sd;
		private final double cov;
		private final DiscretizedFunc cdf;
		private final DiscretizedFunc lec;
		
		public LossDistribution(double[] lossVals, double[] lossWeights, DiscretizedFunc cdf,
				DiscretizedFunc lec) {
			super();
			this.lossVals = lossVals;
			this.lossWeights = lossWeights;
			this.cdf = cdf;
			this.lec = lec;
			
			if (lossVals == null) {
				mean = Double.NaN;
				sd = Double.NaN;
				cov = Double.NaN;
			} else {
				Variance v = new Variance();
				double var = v.evaluate(lossVals, MathArrays.normalizeArray(lossWeights, lossVals.length));
				this.sd = Math.sqrt(var);
				Mean m = new Mean();
				this.mean = m.evaluate(lossVals, lossWeights);
				
				double sumWeights = StatUtils.sum(lossWeights);
				
				// eqns 17-18 of Porter et al. 2020
				double tempSum = 0d;
				for (int i=0; i<lossVals.length; i++)
					tempSum += lossVals[i]*lossVals[i]*lossWeights[i]/sumWeights;
				cov = Math.sqrt(tempSum - mean*mean)/mean;
			}
		}
	}
	
	private static LossDistribution calcBranchAveragedEAL(Collection<U3_EAL_LogicTreeBranch> branches,
			Table<U3_EAL_LogicTreeBranch, Double, Double> branchEALs, DiscretizedFunc xVals) {
		double[] values = new double[branches.size()];
		double[] weights = new double[branches.size()];
		DiscretizedFunc cdf = xVals.deepClone();
		for (int i=0; i<cdf.size(); i++)
			cdf.set(i, 0d);
		
		int index = 0;
		
		double totWeight = 0d;
		
		for (U3_EAL_LogicTreeBranch branch : branches) {
			Map<Double, Double> map = branchEALs.row(branch);
			Preconditions.checkState(map.size() == 1);
			Double weight = map.keySet().iterator().next();
			Double branchEAL = map.get(weight);
			
			totWeight += weight;
			
			values[index] = branchEAL;
			weights[index] = weight;
			
			for (int i=0; i<cdf.size(); i++)
				if (values[index] <= cdf.getX(i))
					cdf.set(i, cdf.getY(i)+weight);
			
			index++;
		}
		// rescale
		cdf.scale(1d/totWeight);
		
		return new LossDistribution(values, weights, cdf, null);
	}
}
