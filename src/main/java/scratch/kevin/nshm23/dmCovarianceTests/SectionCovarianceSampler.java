package scratch.kevin.nshm23.dmCovarianceTests;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.CompletableFuture;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.calc.cholesky.CholeskyDecomposition;
import org.opensha.commons.calc.cholesky.NearPD;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.util.ClassUtils;
import org.opensha.commons.util.Interpolate;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.SectionDistanceAzimuthCalculator;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

import Jama.Matrix;

public abstract class SectionCovarianceSampler {
	
	private List<? extends FaultSection> subSects;
	
	private boolean debug = false;
	
	private transient int prevInterpSkip;
	private transient Matrix prevL;

	public SectionCovarianceSampler(List<? extends FaultSection> subSects) {
		this.subSects = subSects;
	}
	
	public List<? extends FaultSection> getSubSects() {
		return subSects;
	}
	
	/**
	 * @param sect1
	 * @param sect2
	 * @return Pearson correlation coefficient between the two sections in the range [-1,1] 
	 */
	public abstract double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2);
	
	protected double[][] calcCorrs(List<? extends FaultSection> subSects) {
		int numSects = subSects.size();
		System.out.println("Calculating covariance for "+numSects+" subsections");
		double[][] corr = new double[numSects][numSects];
		List<CompletableFuture<Void>> corrFutures = new ArrayList<>(numSects);
		for (int i=0; i<numSects; i++) {
			FaultSection sect1 = subSects.get(i);
			corr[i][i] = 1; // always perfectly correlated with itself
			int myI = i;
			int startJ = i+1;
			corrFutures.add(CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					for (int j=startJ; j<numSects; j++) {
						corr[myI][j] = getCorrelationCoefficient(sect1, subSects.get(j));
						corr[j][myI] = corr[myI][j];
					}
				}
			}));
		}
		for (CompletableFuture<Void> future : corrFutures)
			future.join();
		return corr;
	}
	
	/**
	 * Generates the given number of random fields. Each random field is returned as a random value selected from
	 * a standard normal distribution, which can be treated as a z-score when mapping back to a real property
	 * 
	 * @param numSamples
	 * @param rng
	 * @param center
	 * @param interpSkip
	 * @return
	 */
	public List<double[]> sample(int numSamples, RandomGenerator rng, boolean center, int interpSkip) {
		List<? extends FaultSection> origSects = this.subSects;
		List<? extends FaultSection> calcSects = interpSkip > 0 ? getSubsampled(origSects, interpSkip) : origSects;
		
		int numSects = calcSects.size();
		
		double[][] corr = calcCorrs(calcSects);
		
		// these are spatially correlated gaussians
		prevInterpSkip = interpSkip;
		List<double[]> samples = buildSamples(corr, numSamples, rng);
		
		if (center) {
			System.out.println("Centering random samples");
			MinMaxAveTracker absTrack = new MinMaxAveTracker();
			MinMaxAveTracker track = new MinMaxAveTracker();
			for (int i=0; i<numSects; i++) {
				double sampleMean = 0d;
				for (int n=0; n<numSamples; n++)
					sampleMean += samples.get(n)[i];
				sampleMean /= numSamples;
				
				track.addValue(sampleMean);
				absTrack.addValue(Math.abs(sampleMean));
				
				for (int n=0; n<numSamples; n++) {
					double[] sample = samples.get(n);
					sample[i] = sample[i]-sampleMean;
				}
			}
			
			System.out.println("Centering stats:");
			System.out.println("\tAbsolute: "+absTrack);
			System.out.println("\tOverall: "+track);
		}
		
		if (interpSkip > 0) {
			// interpolate the fields
			List<double[]> interpolatedSamples = new ArrayList<>();
			int numOrig = origSects.size();
			for (double[] sample : samples) {
				double[] interpolated = new double[origSects.size()];
				int lastI2 = -1;
				for (int i1=0,i2=0; i1<numSects; i1++,i2++) {
					Preconditions.checkState(i2 < numOrig);
					FaultSection sect1 = calcSects.get(i1);
					FaultSection sect2 = origSects.get(i2);
					if (sect1 == sect2) {
						// direct match
						interpolated[i2] = sample[i1];
					} else {
						// we hit an interior between interpolated points
						int iStart = i1-1;
						int iEnd = i1;
						int prevID = calcSects.get(iStart).getSectionId();
						int nextID = calcSects.get(iEnd).getSectionId();
						double sampleStart = sample[iStart];
						double sampleEnd = sample[iEnd];
						while (origSects.get(i2).getSectionId() != nextID) {
							interpolated[i2++] = Interpolate.findY(prevID, sampleStart, nextID, sampleEnd, origSects.get(i2).getSectionId());
							Preconditions.checkState(i2 < numOrig);
						}
						interpolated[i2] = sampleEnd;
					}
					lastI2 = i2;
				}
				Preconditions.checkState(lastI2 == numOrig-1, "Last I2 set was %s, expected %s", lastI2, numOrig-1);
				interpolatedSamples.add(interpolated);
			}
			samples = interpolatedSamples;
		}
		
//		System.out.println("Building updated minisections");
//		List<Map<Integer, List<MinisectionSlipRecord>>> ret = new ArrayList<>();
//		for (int n=0; n<numSamples; n++) {
//			Map<Integer, List<MinisectionSlipRecord>> sampleMinis = new HashMap<>(minis.size());
//			int curParent = -1;
//			List<MinisectionSlipRecord> curMinis = null;
//			
//			for (int i=0; i<ordered.size(); i++) {
//				MinisectionSlipRecord rec = ordered.get(i);
//				double slip = rec.slipRate;
//				double stdDev = rec.slipRate;
//				
//				double randSlip = Math.max(0d, slip + samples[n][i]*stdDev);
//				MinisectionSlipRecord mod = new MinisectionSlipRecord(rec.parentID, rec.minisectionID,
//						rec.startLoc, rec.endLoc, rec.rake, randSlip, rec.slipRateStdDev);
//				if (rec.parentID != curParent) {
//					curParent = rec.parentID;
//					curMinis = new ArrayList<>();
//					sampleMinis.put(curParent, curMinis);
//				}
//				curMinis.add(mod);
//			}
//			ret.add(sampleMinis);
//		}
		
		System.out.println("DONE");
		return samples;
	}
	
	/**
	 * If true, will print out the matrices used for the first random field
	 * @param debug
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
	}
	
	public boolean isDebug() {
		return debug;
	}
	
	List<? extends FaultSection> getSubsampled(List<? extends FaultSection> subSects, int interpSkip) {
		// first bin them by parent sections
		List<List<FaultSection>> parentBinned = new ArrayList<>();
		int curParentID = -1;
		List<FaultSection> curList = null;
		int prevIndex = -1;
		for (int i=0; i<subSects.size(); i++) {
			FaultSection sect = subSects.get(i);
			Preconditions.checkState(sect.getSectionId() > prevIndex, "Subsections must be in order if interpSkip>0");
			prevIndex = sect.getSectionId();
			int parentID = sect.getParentSectionId();
			if (parentID != curParentID) {
				curList = new ArrayList<>();
				parentBinned.add(curList);
				curParentID = parentID;
			}
			curList.add(sect);
		}

		List<FaultSection> filtered = new ArrayList<>();

		for (List<FaultSection> sectsForParent : parentBinned) {
			int numSects = sectsForParent.size();
			if (numSects <= 2) {
				// keep both
				filtered.addAll(sectsForParent);
			} else if (numSects <= 2+interpSkip) {
				// keep first and last
				filtered.add(sectsForParent.get(0));
				filtered.add(sectsForParent.get(sectsForParent.size()-1));
			} else {
				// there are interiors that we'll keep
				List<FaultSection> myRetained = getFilteredForParent(sectsForParent, interpSkip);
				if (interpSkip > 1) {
					// see if we should reducing returns the same number
					for (int testInterpSkip=interpSkip; --testInterpSkip>=1;) {
						List<FaultSection> testRetained = getFilteredForParent(sectsForParent, testInterpSkip);
						if (testRetained.size() == myRetained.size()) {
							// the same number were retained, use this one as likely to be more centered
							myRetained = testRetained;
						} else {
							break;
						}
					}
				}
				filtered.addAll(myRetained);
			}
		}
		
		System.out.println("Filtered down to "+filtered.size()+" sections with interpSkip="+interpSkip);
		
		return filtered;
	}
	
	private List<FaultSection> getFilteredForParent(List<FaultSection> sectsForParent, int interpSkip) {
		Preconditions.checkState(interpSkip > 0);
		List<FaultSection> retained = new ArrayList<>();
		// always keep first
		retained.add(sectsForParent.get(0));
		int curSkipCount = 0;
		for (int i=1; i<sectsForParent.size()-1; i++) {
			if (curSkipCount == interpSkip) {
				// we've skipped enough, keep this one
				curSkipCount = 0;
				retained.add(sectsForParent.get(i));
			} else {
				curSkipCount++;
			}
		}
		// always keep last
		retained.add(sectsForParent.get(sectsForParent.size()-1));
		return retained;
	}
	
	protected List<double[]> buildSamples(double[][] covs, int numSamples, RandomGenerator rng) {
		Preconditions.checkState(covs.length == covs[0].length);
		Matrix C = new Matrix(covs);
		
		if (debug) {
			System.out.println("COV Stats");
			printMatrixStats(C);
		}
		
		System.out.println("Doing CholeskyDecomposition");
		Matrix L = decompose(C).getL();
		
		prevL = L;
		
		if (debug) {
			System.out.println("L Stats");
			printMatrixStats(L);
		}
		
		return buildSamples(L, numSamples, rng);
	}
	
	protected List<double[]> buildSamples(Matrix L, int numSamples, RandomGenerator rng) {
		int numVars = L.getRowDimension();
		System.out.println("Generating "+numSamples+" samples");
		List<double[]> samples = new ArrayList<>();
		for (int n=0; n<numSamples; n++) {
			Matrix X = randomGaussian(numVars, rng);
			if (debug && n == 0) {
				System.out.println("Debugging field 1");
				System.out.println("Gaussian:");
//				print1D(X, 20);
				print1DStats(X);
			}
			Matrix prod = L.times(X);
			Preconditions.checkState(prod.getRowDimension() == numVars,
					"Expected %s rows, have %s", numVars, prod.getRowDimension());
			Preconditions.checkState(prod.getColumnDimension() == 1,
					"Expected %s cols, have %s", 1, prod.getColumnDimension());
			if (debug && n == 0) {
				System.out.println("Random field:");
//				print1D(prod, 20);
				print1DStats(prod);
			}
			double[] sample = new double[numVars];
			for (int i=0; i<numVars; i++)
				sample[i] = prod.get(i, 0);
			samples.add(sample);
		}
		
		return samples;
	}
	
	private void print1D(Matrix mat, int wrap) {
		Preconditions.checkState(mat.getColumnDimension() == 1);
		int num = mat.getRowDimension();
		DecimalFormat df = new DecimalFormat("0.00");
		StringBuilder line = new StringBuilder();
		int curNum = 0;
		for (int i=0; i<num; i++) {
			if (curNum == wrap) {
				System.out.println(line);
				line = new StringBuilder();
				curNum = 0;
			}
			line.append("\t").append(df.format(mat.get(i, 0)));
			curNum++;
		}
		System.out.println(line);
	}
	
	private void print1DStats(Matrix mat) {
		Preconditions.checkState(mat.getColumnDimension() == 1);
		MinMaxAveTracker track = new MinMaxAveTracker();
		int numZero = 0;
		int numNonFinite = 0;
		int num = mat.getRowDimension();
		for (int i=0; i<num; i++) {
			double val = mat.get(i, 0);
			if (val == 0d)
				numZero++;
			if (!Double.isFinite(val))
				numNonFinite++;
			else
				track.addValue(val);
		}
		System.out.println("\tStats:\t"+track);
		System.out.println("\tExcatly zero:\t"+numZero);
		System.out.println("\tNon-finite:\t"+numNonFinite);
	}
	
	private Matrix randomGaussian(int num, RandomGenerator rng) {
		Matrix X = new Matrix(num, 1);
		for (int i=0; i<num; i++)
			X.set(i, 0, rng.nextGaussian());
		return X;
	}
	
	private CholeskyDecomposition decompose(Matrix B) {
		CholeskyDecomposition chol = new CholeskyDecomposition(B);
		if (!chol.isSPD()) {
			System.out.println("Matrix is not positive definite, will try to iteratively find nearPD solution");
			NearPD nearPD = new NearPD();
			nearPD.setKeepDiag(true);
			nearPD.setUseApache(true);
//			nearPd.setEigTol(1.e-6);
			boolean success = nearPD.calcNearPD(B, debug);
			double normFrob = nearPD.getFrobNorm();
			if (!success) {
				System.out.println("B is size "+B.getRowDimension()+"x"+B.getColumnDimension()+". normFrob="+normFrob);
//				printMatrix(B);
//				throw new RuntimeException("Error: nearPD failed to converge, the correlation matrix maybe" +
//						" significantly different from a PD matrix, check that the correlation equations" +
//						" used are reasonable");
				System.err.println("WARNING: nearPD failed to converge, the correlation matrix maybe" +
						" significantly different from a PD matrix, check that the correlation equations" +
						" used are reasonable. Convergence: "+(float)nearPD.getConvergedTolerence());
			}
			
			Matrix x = nearPD.getX();
			if (debug) {
				System.out.println("Nearest PD matrix stats:");
				printMatrixStats(x);
			}
			//Now get the CholDecomp of this nearest matrix
			CholeskyDecomposition cholPD = new CholeskyDecomposition(x);
			Preconditions.checkState(cholPD.isSPD(), "Error: Even after NearPD the matrix is not PD");
			chol = cholPD;
		}
		if (debug) {
			// try it the apache way
			try {
				RealMatrix apacheMatrix = new Array2DRowRealMatrix(B.getArray());
				System.out.println("Trying it the apache way...");
				org.apache.commons.math3.linear.CholeskyDecomposition apacheChol =
						new org.apache.commons.math3.linear.CholeskyDecomposition(apacheMatrix);
				System.out.println("Success, L stats:");
				Matrix apacheL = new Matrix(apacheChol.getL().getData());
				printMatrixStats(apacheL);
			} catch (Exception e) {
				System.out.println("Failed apache decomposition: "+e.getMessage());
			}
		}
		return chol;
	}
	
	private void printMatrixStats(Matrix matrix) {
		int numZeros = 0;
		int num1s = 0;
		int numDiag1s = 0;
		int numNonZeroOnes = 0;
		int upperNonZeros = 0;
		int lowerNonZeros = 0;
		boolean symmetric = matrix.getRowDimension() == matrix.getColumnDimension();
		for (int row=0; row<matrix.getRowDimension(); row++) {
			for (int col=0; col<matrix.getColumnDimension(); col++) {
				double val = matrix.get(row, col);
				if (val == 0d) {
					numZeros++;
				} else if (val == 1d) {
					num1s++;
					if (col == row)
						numDiag1s++;
				} else {
					numNonZeroOnes++;
				}
				if (val != 0d && row != col) {
					if (col > row)
						upperNonZeros++;
					else
						lowerNonZeros++;
				}
				symmetric &= val == matrix.get(col, row);
			}
		}
		long totNum = (long)matrix.getRowDimension()*(long)matrix.getColumnDimension();
		DecimalFormat pDF = new DecimalFormat("0.00%");
		System.out.println("\tZeros:\t"+numZeros+" ("+pDF.format((double)numZeros/(double)totNum)+")");
		System.out.println("\tOnes:\t"+num1s+" ("+pDF.format((double)num1s/(double)totNum)+")");
		System.out.println("\tOnes on diagonal:\t"+numDiag1s+" ("+pDF.format((double)numDiag1s/(double)matrix.getRowDimension())+")");
		System.out.println("\tOthers:\t"+numNonZeroOnes+" ("+pDF.format((double)numNonZeroOnes/(double)totNum)+")");
		long numTriangle = (totNum - matrix.getRowDimension())/2;
		System.out.println("\tUpper nonzeros:\t"+upperNonZeros+" ("+pDF.format((double)upperNonZeros/(double)numTriangle)+")");
		System.out.println("\tLower nonzeros:\t"+lowerNonZeros+" ("+pDF.format((double)lowerNonZeros/(double)numTriangle)+")");
		System.out.println("\tSymmetric?\t"+symmetric);
	}
	
//	public abstract static class CorrelationBased extends SectCovarianceEstimator {
//
//		public CorrelationBased(Map<Integer, List<MinisectionSlipRecord>> minis) {
//			super(minis);
//		}
//		
//		public final double getSlipCovariance(MinisectionSlipRecord mini1, MinisectionSlipRecord mini2) {
//			double corr = getSlipCorrelation(mini1, mini2);
//			Preconditions.checkState(corr >= -1d && corr <= 1d, "Bad correlation coefficient: %s", corr);
//			// corr = cov/(sigma1*sigma2)
////			return corr*mini1.slipRateStdDev*mini2.slipRateStdDev;
//			return corr;
//		}
//		
//		public abstract double getSlipCorrelation(MinisectionSlipRecord mini1, MinisectionSlipRecord mini2);
//		
//	}
	
	static class QuickTestInvDistance extends SectionCovarianceSampler {

		private SectionDistanceAzimuthCalculator distCalc;
		private double maxDist;
		private double covForZeroDist;

		public QuickTestInvDistance(List<? extends FaultSection> subSects, SectionDistanceAzimuthCalculator distCalc,
				double maxDist, double covForZeroDist) {
			super(subSects);
			this.distCalc = distCalc;
			this.maxDist = maxDist;
			this.covForZeroDist = covForZeroDist;
			Preconditions.checkState(covForZeroDist >= 0);
		}

		@Override
		public double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2) {
			double dist = distCalc.getDistance(sect1, sect2);
			if (dist > maxDist)
				return 0d;
			if (dist == 0d)
				return covForZeroDist;
			return Interpolate.findY(0d, covForZeroDist, maxDist, 0d, dist);
		}

		@Override
		protected String getSamplerPrefix() {
			return "inv_dist_"+(float)maxDist+"km_"+(float)covForZeroDist+"forZero";
		}
		
	}
	
	static class QuickTestAllZeros extends SectionCovarianceSampler {

		public QuickTestAllZeros(List<? extends FaultSection> subSects) {
			super(subSects);
		}

		@Override
		public double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2) {
			return 0d;
		}

		@Override
		protected String getSamplerPrefix() {
			return "all_zero_corr";
		}
		
	}
	
	static class QuickTestRandCorr extends SectionCovarianceSampler {

		public QuickTestRandCorr(List<? extends FaultSection> subSects) {
			super(subSects);
		}

		@Override
		public double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2) {
			return -1d + Math.random()*2;
		}

		@Override
		protected String getSamplerPrefix() {
			return "random_corr";
		}
		
	}
	
	protected abstract String getSamplerPrefix();
	
	protected String getFullCachePrefix(int interpSkip) {
		return "covariance_cache_"+getSamplerPrefix()
			+"_"+SectionDistanceAzimuthCalculator.getUniqueSectCacheFileStr(subSects)+"_interp"+interpSkip;
	}
	
	public void writeCache(File outputDir) throws IOException {
		Preconditions.checkNotNull(prevL, "CholeskyDecomposition is not yet calculated");
		CachedDecomposition cached = new CachedDecomposition(this, prevL, prevInterpSkip);
		String prefix = getFullCachePrefix(prevInterpSkip);
		cached.writeCache(outputDir, prefix);
	}
	
	public Optional<SectionCovarianceSampler> loadCached(File cacheDir, int interpSkip) throws IOException {
		if (!cacheDir.exists())
			return Optional.empty();
		String prefix = getFullCachePrefix(interpSkip);
		File lCSVFile = new File(cacheDir, prefix+"_L.csv");
		File corrCSVFile = new File(cacheDir, prefix+"_corr.csv");
		if (!lCSVFile.exists() || !corrCSVFile.exists())
			return Optional.empty();
		
		System.out.println("Reading L matrix from "+lCSVFile.getAbsolutePath());
		CSVFile<Double> lCSV = CSVFile.readFileNumeric(lCSVFile, true, 0);
		int len = lCSV.getNumRows();
		Matrix L = new Matrix(len, len);
		for (int row=0; row<len; row++)
			for (int col=0; col<len; col++)
				L.set(row, col, lCSV.getDouble(row, col));
		
		System.out.println("Reading correlations from "+corrCSVFile.getAbsolutePath());
		CSVFile<String> corrCSV = CSVFile.readFile(corrCSVFile, true);
		Preconditions.checkState(corrCSV.getNumRows() == subSects.size());
		double[][] corrs = new double[subSects.size()][subSects.size()];
		for (int i=0; i<subSects.size(); i++) {
			int id = subSects.get(i).getSectionId();
			Preconditions.checkState(corrCSV.getInt(i, 0) == id);
			for (int j=0; j<subSects.size(); j++)
				corrs[i][j] = corrCSV.getDouble(i, j+1);
		}
		
		return Optional.of(new CachedDecomposition(subSects, L, corrs, interpSkip));
	}
	
	static class CachedDecomposition extends SectionCovarianceSampler {

		private Matrix L;
		private int interpSkip;
		private double[][] corrs;
		private Map<Integer, Integer> sectIndexMap;

		public CachedDecomposition(SectionCovarianceSampler sampler, Matrix L, int interpSkip) {
			this(sampler.getSubSects(), L, sampler.calcCorrs(sampler.getSubSects()), interpSkip);
		}
		
		public CachedDecomposition(List<? extends FaultSection> subSects, Matrix L, double[][] corrs, int interpSkip) {
			super(subSects);
			this.L = L;
			this.interpSkip = interpSkip;
			this.corrs = corrs;
			
			boolean idsAreIndexes = true;
			sectIndexMap = new HashMap<>(subSects.size());
			for (int i=0; i<subSects.size(); i++) {
				int id = subSects.get(i).getSectionId();
				idsAreIndexes &= id == i;
				sectIndexMap.put(id, i);
			}
			if (idsAreIndexes)
				sectIndexMap = null;
		}

		@Override
		public double getCorrelationCoefficient(FaultSection sect1, FaultSection sect2) {
			int index1 = sectIndexMap == null ? sect1.getSectionId() : sectIndexMap.get(sect1.getSectionId());
			int index2 = sectIndexMap == null ? sect2.getSectionId() : sectIndexMap.get(sect2.getSectionId());
			return corrs[index1][index2];
		}
		
		@Override
		public List<double[]> sample(int numSamples, RandomGenerator rng, boolean center, int interpSkip) {
			Preconditions.checkState(interpSkip == this.interpSkip, "InterpSkip mismatch, can't use cache");
			return super.sample(numSamples, rng, center, interpSkip);
		}

		protected List<double[]> buildSamples(double[][] covs, int numSamples, RandomGenerator rng) {
			Preconditions.checkState(L.getRowDimension() == covs.length,
					"Passed in COVs are an unxpected size. COV len=%s, L len=%s", covs.length, L.getRowDimension());
			
			return buildSamples(L, numSamples, rng);
		}
		
		private void writeCache(File outputDir, String prefix) throws IOException {
			Preconditions.checkState(outputDir.exists() || outputDir.mkdir(),
					"Output directory doesn't exist and couldn't be created: %s", outputDir.getAbsolutePath());
			Preconditions.checkState(outputDir.isDirectory(),
					"Output directory isn't a directory: %s", outputDir.getAbsolutePath());
			System.out.println("Writing cache to "+outputDir.getAbsolutePath()+" with prefix: "+prefix);
			
			CSVFile<Double> lCSV = new CSVFile<>(true);
			for (int row=0; row<L.getRowDimension(); row++) {
				List<Double> line = new ArrayList<>(L.getColumnDimension());
				for (int col=0; col<L.getColumnDimension(); col++)
					line.add(L.get(row, col));
				lCSV.addLine(line);
			}
			lCSV.writeToFile(new File(outputDir, prefix+"_L.csv"));
			
			List<? extends FaultSection> subSects = getSubSects();
			
			CSVFile<Number> corrCSV = new CSVFile<>(true);
			for (int i=0; i<subSects.size(); i++) {
				List<Number> line = new ArrayList<>(subSects.size()+1);
				line.add(subSects.get(i).getSectionId());
				for (int j=0; j<subSects.size(); j++)
					line.add(corrs[i][j]);
				corrCSV.addLine(line);
			}
			corrCSV.writeToFile(new File(outputDir, prefix+"_corr.csv"));
		}

		@Override
		public void writeCache(File outputDir) throws IOException {
			// do nothing (already cached)
		}

		@Override
		protected String getSamplerPrefix() {
			throw new IllegalStateException("Should never be called on an already cached version");
		}
		
	}

}

