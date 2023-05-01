package scratch.kevin.spatialVar;

import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.xyz.GeoDataSet;
import org.opensha.commons.data.xyz.GeoDataSetMath;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.geo.Region;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.gcim.calc.CholeskyDecomposition;
import org.opensha.sha.gcim.calc.NearPD;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

//import Jama.CholeskyDecomposition;
import Jama.Matrix;
import scratch.UCERF3.analysis.FaultBasedMapGen;

public class SpatialVarCalc {
	
	private double[] periods;
	private int numSites;
	
	private Matrix K1, K2, K3;
	private Matrix L1, L2;

	public SpatialVarCalc(double[] periods, List<? extends Site> sites) {
		this.periods = periods;
		this.numSites = sites.size();
		
		Matrix B1 = new Matrix(periods.length, periods.length);
		Matrix B2 = new Matrix(periods.length, periods.length);
		Matrix B3 = new Matrix(periods.length, periods.length);
		
		double minSupportedPeriod = StatUtils.min(LothBaker2013_SpatialVarCalc.periods);
		double maxSupportedPeriod = StatUtils.max(LothBaker2013_SpatialVarCalc.periods);
		double[] modPeriods = new double[periods.length];
		for (int i=0; i<modPeriods.length; i++) {
			double period = periods[i];
			Preconditions.checkState(period >= 0d);
			if (period < minSupportedPeriod)
				period = minSupportedPeriod;
			if (period > maxSupportedPeriod)
				period = maxSupportedPeriod;
			if (period != periods[i])
				System.out.println("WARNING: period of "+(float)periods[i]
						+" out of range from Loth-Baker, using "+(float)period);
			modPeriods[i] = period;
		}
		periods = modPeriods;
		
		// interpolate b matrices at the given periods
		for (int i=0; i<periods.length; i++) {
			double period1 = periods[i];
			int ind1 = LothBaker2013_SpatialVarCalc.getPeriodIndexBefore(period1);
			for (int j=i; j<periods.length; j++) {
				double period2 = periods[j];
				int ind2 = LothBaker2013_SpatialVarCalc.getPeriodIndexBefore(period2);
				B1.set(i, j, LothBaker2013_SpatialVarCalc.getInterpB(
						LothBaker2013_SpatialVarCalc.b1, period1, period2, ind1, ind2));
				B2.set(i, j, LothBaker2013_SpatialVarCalc.getInterpB(
						LothBaker2013_SpatialVarCalc.b2, period1, period2, ind1, ind2));
				B3.set(i, j, LothBaker2013_SpatialVarCalc.getInterpB(
						LothBaker2013_SpatialVarCalc.b3, period1, period2, ind1, ind2));
				
				B1.set(j, i, B1.get(i, j));
				B2.set(j, i, B2.get(i, j));
				B3.set(j, i, B3.get(i, j));
			}
		}
		
		System.out.println("B1");
		printMatrix(B1);
		System.out.println();
		System.out.println("B2");
		printMatrix(B2);
		System.out.println();
		System.out.println("B3");
		printMatrix(B3);
		System.out.println();
		
		// normalize B matrices
		Matrix B = B1.plus(B2).plus(B3);
//		System.out.println("B");
//		printMatrix(B);
//		System.out.println();
		for (int i=0; i<periods.length; i++) {
			for (int j=0; j<periods.length; j++) {
				B1.set(i, j, B1.get(i, j)/Math.sqrt(B.get(i, i)*B.get(j, j)));
				B2.set(i, j, B2.get(i, j)/Math.sqrt(B.get(i, i)*B.get(j, j)));
				B3.set(i, j, B3.get(i, j)/Math.sqrt(B.get(i, i)*B.get(j, j)));
			}
		}
		
		System.out.println("B1 norm");
		printMatrix(B1);
		System.out.println();
		System.out.println("B2 norm");
		printMatrix(B2);
		System.out.println();
		System.out.println("B3 norm");
		printMatrix(B3);
		System.out.println();
//		System.exit(0);
		
		K1 = decompose(B1).getL();
		K2 = decompose(B2).getL();
		K3 = decompose(B3).getL();
		
//		System.out.println("K1");
//		printMatrix(B1);
//		System.out.println();
//		System.out.println("K2");
//		printMatrix(B2);
//		System.out.println();
//		System.out.println("K3");
//		printMatrix(B3);
//		System.out.println();
//		System.exit(0);
		
		// coregionalization matrices
		Matrix D1 = new Matrix(sites.size(), sites.size());
		Matrix D2 = new Matrix(sites.size(), sites.size());
		for (int s1=0; s1<sites.size(); s1++) {
			Location l1 = sites.get(s1).getLocation();
			for (int s2=s1; s2<sites.size(); s2++) {
				Location l2 = sites.get(s2).getLocation();
				if (s1 == s2) {
					D1.set(s1, s2, 1);
					D2.set(s1, s2, 1);
				} else {
					double h = LocationUtils.horzDistanceFast(l1, l2);
					D1.set(s1, s2, Math.exp(-3d*h/20d));
					D2.set(s1, s2, Math.exp(-3d*h/70d));
					D1.set(s2, s1, D1.get(s1, s2));
					D2.set(s2, s1, D2.get(s1, s2));
				}
			}
		}
		
		L1 = decompose(D1).getL();
		L2 = decompose(D2).getL();
	}
	
	private SpatialVarCalc(double[] periods, int numSites, Matrix K1, Matrix K2, Matrix K3, Matrix L1, Matrix L2) {
		this.periods = periods;
		this.numSites = numSites;
		this.K1 = K1;
		this.K2 = K2;
		this.K3 = K3;
		this.L1 = L1;
		this.L2 = L2;
	}
	
	public Matrix[] computeRandWithinEventResiduals(Random rng, double sigma, int num) {
		Matrix[] corrFields = new Matrix[num];
		for (int i=0; i<num; i++)
			corrFields[i] = computeRandomField(rng, sigma);
		
		int rows = corrFields[0].getRowDimension();
		int cols = corrFields[0].getColumnDimension();
		
		// get log mean of all fields (at all site-periods)
		Matrix logMean = new Matrix(rows, cols);
		for (int i=0; i<rows; i++) {
			for (int j=0; j<cols; j++) {
				double sum = 0d;
				for (int n=0; n<num; n++)
					sum += Math.log(corrFields[n].get(i, j));
				logMean.set(i, j, sum/(double)num);
			}
		}
		
		// now convert to within event residual
		Matrix[] withinFields = new Matrix[num];
		for (int n=0; n<num; n++)
			withinFields[n] = new Matrix(rows, cols);
		for (int i=0; i<rows; i++) {
			for (int j=0; j<cols; j++) {
				double mean = logMean.get(i, j);
				for (int n=0; n<num; n++)
					withinFields[n].set(i, j, Math.log(corrFields[n].get(i, j)) - mean);
			}
		}
		
		return withinFields;
	}
	
	public Matrix computeRandomField(Random rng, double sigma) {
		Matrix R1 = normRandArray(periods.length, numSites, rng);
		Matrix R2 = normRandArray(periods.length, numSites, rng);
		Matrix R3 = normRandArray(periods.length, numSites, rng);
		
		Matrix S1 = K1.times(R1).times(L1);
		Matrix S2 = K2.times(R2).times(L2);
		Matrix S3 = K3.times(R3);  

		// code divides by 3 but Nan Wang via e-mail said not to (which is in line with docs)
		// corr_log = (corr1 + corr2 + corr3)/3;
//		Matrix S = S1.plus(S2).plus(S3).times(1d/3d);
		Matrix S = S1.plus(S2).plus(S3);
		
//		% % take the exp of correlated random number and multiply
//		corr = exp(-(sigma.^2./2)).*exp(sigma.*corr_log);
		double sigMult = Math.exp(-(sigma*sigma*0.5));
		for (int i=0; i<S.getRowDimension(); i++) {
			for (int j=0; j<S.getColumnDimension(); j++) {
				S.set(i, j, sigMult*Math.exp(sigma*S.get(i, j)));
			}
		}
		
		return S;
	}
	
	private Matrix normRandArray(int m, int n, Random rng) {
		Matrix ret = new Matrix(m, n);
		for (int i=0; i<m; i++)
			for (int j=0; j<n; j++)
				ret.set(i, j, rng.nextGaussian());
		return ret;
	}
	
	private static void printMatrix(Matrix mat) {
		for (int i=0; i<mat.getRowDimension(); i++) {
			for (int j=0; j<mat.getColumnDimension(); j++)
				System.out.print(mat.get(i, j)+"\t");
			System.out.println();
		}
		System.out.flush();
	}
	
	private CholeskyDecomposition decompose(Matrix B) {
		CholeskyDecomposition chol = new CholeskyDecomposition(B);
		if (!chol.isSPD()) {
			NearPD nearPD = new NearPD();
			nearPD.setKeepDiag(true);
//			nearPd.setEigTol(1.e-6);
			boolean success = nearPD.calcNearPD(B);
			double normFrob = nearPD.getFrobNorm();
			if (!success) {
				System.out.println("B is size "+B.getRowDimension()+"x"+B.getColumnDimension()+". normFrob="+normFrob);
				printMatrix(B);
				throw new RuntimeException("Error: nearPD failed to converge, the correlation matrix maybe" +
						" significantly different from a PD matrix, check that the correlation equations" +
						" used are reasonable");
			}
			
			Matrix x = nearPD.getX();
			//Now get the CholDecomp of this nearest matrix
			CholeskyDecomposition cholPD = new CholeskyDecomposition(x);
			Preconditions.checkState(cholPD.isSPD(), "Error: Even after NearPD the matrix is not PD");
			chol = cholPD;
		}
		return chol;
	}
	
	public GeoDataSet calcRandomShakeMap(GeoDataSet input, Matrix S, int p) {
		return calcRandomShakeMap(input, S, p, Double.NaN);
	}
	
	public GeoDataSet calcRandomShakeMap(GeoDataSet input, Matrix S, int p, double truncation) {
		Preconditions.checkState(S.getRowDimension() == periods.length);
		Preconditions.checkState(S.getColumnDimension() == numSites);
		Preconditions.checkState(p >= 0 && p < periods.length);
		Preconditions.checkState(input.size() == numSites);
		GeoDataSet ret = input.copy();
		for (int i=0; i<ret.size(); i++) {
			double lnOrig = Math.log(ret.get(i));
//			double lnRes = Math.log(S.get(p, i));
			double lnRes = S.get(p, i);
			if (truncation > 0) {
				if (lnRes >= 0)
					lnRes = Math.min(truncation, lnRes);
				else
					lnRes = Math.max(-truncation, lnRes);
			}
			ret.set(i, Math.exp(lnOrig+lnRes));
		}
		return ret;
	}
	
	private static String getCachePrefix(int numSites, double[] periods) {
		String prefix = "sp_corr_"+numSites+"sites";
		DecimalFormat df = new DecimalFormat("0.###");
		for (double period : periods)
			prefix += "_"+df.format(period)+"s";
		return prefix;
	}
	
	public void writeCache(File cacheDir) throws IOException {
		String prefix = getCachePrefix(numSites, periods);
		System.out.println("Writing spatial correlation matrices to "+cacheDir.getAbsolutePath()+" with prefix: "+prefix);
		writeMatrixCSV(new File(cacheDir, prefix+"_K1.csv"), K1);
		writeMatrixCSV(new File(cacheDir, prefix+"_K2.csv"), K2);
		writeMatrixCSV(new File(cacheDir, prefix+"_K3.csv"), K3);
		writeMatrixCSV(new File(cacheDir, prefix+"_L1.csv.gz"), L1);
		writeMatrixCSV(new File(cacheDir, prefix+"_L2.csv.gz"), L2);
	}
	
	public static void writeMatrixCSV(File outputFile, Matrix matrix) throws IOException {
		CSVFile<Double> csv = new CSVFile<>(true);
		for (int row=0; row<matrix.getRowDimension(); row++) {
			List<Double> line = new ArrayList<>(matrix.getColumnDimension());
			for (int col=0; col<matrix.getColumnDimension(); col++)
				line.add(matrix.get(row, col));
			csv.addLine(line);
		}
		if (outputFile.getName().endsWith(".gz")) {
			FileOutputStream fout = new FileOutputStream(outputFile);
			GZIPOutputStream gout = new GZIPOutputStream(fout);
			csv.writeToStream(gout);
			gout.close();
		} else {
			csv.writeToFile(outputFile);
		}
	}
	
	private static Matrix loadMatrix(File matFile) throws NumberFormatException, IOException {
		CSVFile<Double> csv;
		if (matFile.getName().endsWith(".gz")) {
			FileInputStream fin = new FileInputStream(matFile);
			GZIPInputStream gin = new GZIPInputStream(fin);
			csv = CSVFile.readStreamNumeric(gin, true, -1, 0);
			gin.close();
		} else {
			csv = CSVFile.readFileNumeric(matFile, true, 0);
		}
		Matrix matrix = new Matrix(csv.getNumRows(), csv.getNumCols());
		for (int row=0; row<matrix.getRowDimension(); row++)
			for (int col=0; col<matrix.getColumnDimension(); col++)
				matrix.set(row, col, csv.get(row, col));
		return matrix;
	}
	
	public static boolean isCached(File cacheDir, int numSites, double[] periods) {
		String prefix = getCachePrefix(numSites, periods);
		String[] suffixes = { "_K1.csv", "_K2.csv", "_K3.csv", "_L1.csv.gz", "_L2.csv.gz" };
		for (String suffix : suffixes)
			if (!new File(cacheDir, prefix+suffix).exists())
				return false;
		return true;
	}
	
	public static SpatialVarCalc loadCache(File cacheDir, double[] periods, int numSites)
			throws NumberFormatException, IOException {
		String prefix = getCachePrefix(numSites, periods);
		System.out.println("Loading spatial correlation matrices from "+cacheDir.getAbsolutePath()+" with prefix: "+prefix);
		Matrix K1 = loadMatrix(new File(cacheDir, prefix+"_K1.csv"));
		Matrix K2 = loadMatrix(new File(cacheDir, prefix+"_K2.csv"));
		Matrix K3 = loadMatrix(new File(cacheDir, prefix+"_K3.csv"));
		Matrix L1 = loadMatrix(new File(cacheDir, prefix+"_L1.csv.gz"));
		Matrix L2 = loadMatrix(new File(cacheDir, prefix+"_L2.csv.gz"));
		
		return new SpatialVarCalc(periods, numSites, K1, K2, K3, L1, L2);
	}
	
	public int getNumSites() {
		return numSites;
	}
	
	public static void main(String[] args) throws IOException, GMT_MapException {
		File outputDir = new File("/tmp/spatial_var_test");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		Region reg = new Region(new Location(33.5, -119), new Location(35.5, -117));
//		double spacing = 0.025d;
		double spacing = 0.05d;
		GriddedRegion gridReg = new GriddedRegion(reg, spacing, null);
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		List<Site> sites = new ArrayList<>();
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			Site site = new Site(gridReg.getLocation(i));
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>)param.clone());
			sites.add(site);
		}
		boolean[] siteDatas = { false, true };
		
		double[] periods = { 3d };
		double sigma = 0.6d;
		int num = 100;
		int maxPlotNum = 5;
		
		Location center = new Location(0.5*(reg.getMaxLat() + reg.getMinLat()), 0.5*(reg.getMaxLon()+reg.getMinLon()));
		PointSurface surf = new PointSurface(center);
		surf.setAveDip(90d);
		EqkRupture rup = new EqkRupture(7.5d, 180d, surf, center);
		
		System.out.println("Computing GMPE shakemap for "+sites.size()+" sites");
		GriddedGeoDataSet[] shakemaps = new GriddedGeoDataSet[periods.length];
		
		System.out.println("Initializing spatial var calc");
		
		Stopwatch watch = Stopwatch.createStarted();
		SpatialVarCalc calc = new SpatialVarCalc(periods, sites);
		watch.stop();
		double secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Took "+(float)secs+" to initialize spatial variation matrices");
		
		System.out.println("Computing "+num+" random fields");
		watch = Stopwatch.createStarted();
		Matrix[] S = calc.computeRandWithinEventResiduals(new Random(sites.size()), sigma, num);
		watch.stop();
		secs = watch.elapsed(TimeUnit.MILLISECONDS)/1000d;
		System.out.println("Took "+(float)secs+" to compute random fields");
		
		for (int p=0; p<periods.length; p++) {
			System.out.println("Residual statistics for "+num+" realizations:");
			StandardDeviation totalStdDev = new StandardDeviation();
			Mean totalMean = new Mean();
			for (int n=0; n<num; n++) {
				Mean resMean = new Mean();
				StandardDeviation resStdDev = new StandardDeviation();
				for (int i=0; i<sites.size(); i++) {
					double lnRes = S[n].get(p, i);
					resMean.increment(lnRes);
					totalMean.increment(lnRes);
					totalStdDev.increment(lnRes);
					resStdDev.increment(lnRes);
				}
				if (n < maxPlotNum)
					System.out.println("\tResiduals "+n+":\tmean="+(float)resMean.getResult()+"\tstdDev="+(float)resStdDev.getResult());
				else if (n == maxPlotNum)
					System.out.println("\t...");
			}
			System.out.println("Total:\tmean="+(float)totalMean.getResult()+"\tstdDev="+(float)totalStdDev.getResult());
		}
		
		SiteTranslator trans = new SiteTranslator();
		for (boolean siteData : siteDatas) {
			System.out.println("Site data: "+siteData);
			
			ArrayList<SiteDataValueList<?>> datas = null;
			if (siteData) {
				System.out.println("Getting site datas...");
				OrderedSiteDataProviderList provs = OrderedSiteDataProviderList.createSiteDataProviderDefaults();
				datas = provs.getAllAvailableData(sites);
				for (int i=0; i<sites.size(); i++) {
					ArrayList<SiteDataValue<?>> siteVals = new ArrayList<SiteDataValue<?>>();
					for (SiteDataValueList<?> valList : datas) {
						siteVals.add(valList.getValue(i));
					}
					for (Parameter<?> param : sites.get(i))
						trans.setParameterValue(param, siteVals);
				}
				System.out.println("Set site data");
			}
			
			for (int p=0; p<periods.length; p++) {
				shakemaps[p] = new GriddedGeoDataSet(gridReg, false);
				gmpe.setIntensityMeasure(SA_Param.NAME);
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), periods[p]);
				gmpe.setEqkRupture(rup);
				for (int s=0; s<sites.size(); s++) {
					Site site = sites.get(s);
					if (siteData) {
						site = new Site(site.getLocation());
						ArrayList<SiteDataValue<?>> siteVals = new ArrayList<SiteDataValue<?>>();
						for (SiteDataValueList<?> valList : datas)
							siteVals.add(valList.getValue(s));
						for (Parameter<?> param : gmpe.getSiteParams()) {
							param = (Parameter<?>)param.clone();
							trans.setParameterValue(param, siteVals);
							site.addParameter(param);
						}
					}
					gmpe.setSite(site);
					shakemaps[p].set(s, Math.exp(gmpe.getMean()));
				}
			}
			
			CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-2d, 0d);
			cpt.setAboveMaxColor(cpt.getMaxColor());
			cpt.setBelowMinColor(cpt.getMinColor());
			cpt.setNanColor(Color.GRAY);
			
			for (int p=0; p<periods.length; p++) {
				System.out.println("Computing "+num+" random shakemaps for T="+(float)periods[p]);
				GriddedGeoDataSet myShakeMap = shakemaps[p].copy();
				System.out.println("Orig ShakeMap: "+tracker(myShakeMap));
				
				List<GeoDataSet> randMaps = new ArrayList<>();
				GeoDataSet avgShakeMap = null;
				GeoDataSet lnAvgShakeMap = null;
				for (int n=0; n<num; n++) {
					GeoDataSet randShakeMap = calc.calcRandomShakeMap(shakemaps[p], S[n], p);
					if (n < maxPlotNum)
						System.out.println("\tRand "+n+" ShakeMap:\t"+tracker(randShakeMap));
					else if (n == maxPlotNum)
						System.out.println("\t...");
					randMaps.add(randShakeMap);
					if (avgShakeMap == null) {
						avgShakeMap = randShakeMap.copy();
						lnAvgShakeMap = randShakeMap.copy();
						lnAvgShakeMap.log();
					} else {
						avgShakeMap = GeoDataSetMath.add(avgShakeMap, randShakeMap);
						GeoDataSet tmp = randShakeMap.copy();
						tmp.log();
						lnAvgShakeMap = GeoDataSetMath.add(lnAvgShakeMap, tmp);
					}
				}

				avgShakeMap.scale(1d/(double)num);
				lnAvgShakeMap.scale(1d/(double)num);
				lnAvgShakeMap.exp();
				System.out.println("Orig ShakeMap:\t"+tracker(myShakeMap));
				System.out.println("Avg ShakeMap:\t"+tracker(avgShakeMap));
				System.out.println("Ln Avg ShakeMap:\t"+tracker(lnAvgShakeMap));
				
				System.out.println("Plotting maps");
				
				myShakeMap.log10();
				GMT_Map origMap = new GMT_Map(reg, myShakeMap, spacing, cpt);
				origMap.setUseGMTSmoothing(false);
				origMap.setCustomScaleMin((double)cpt.getMinValue());
				origMap.setCustomScaleMax((double)cpt.getMaxValue());
				origMap.setCustomLabel("Log10 "+(float)periods[p]+"s GMPE Median ShakeMap");
				origMap.setBlackBackground(false);
				String prefix = "shakemap_"+(float)periods[p]+"s";
				if (siteData)
					prefix += "_site_data";
				FaultBasedMapGen.plotMap(outputDir, prefix, false, origMap);
				
				for (int n=0; n<num && n<maxPlotNum; n++) {
					GeoDataSet randShakeMap = randMaps.get(n);
					randShakeMap.log10();
					GMT_Map randMap = new GMT_Map(reg, randShakeMap, spacing, cpt);
					randMap.setUseGMTSmoothing(false);
					randMap.setCustomScaleMin((double)cpt.getMinValue());
					randMap.setCustomScaleMax((double)cpt.getMaxValue());
					randMap.setCustomLabel("Log10 "+(float)periods[p]+"s Random ShakeMap "+n);
					randMap.setBlackBackground(false);
					FaultBasedMapGen.plotMap(outputDir, prefix+"_rand"+n, false, randMap);
				}
				
				avgShakeMap.log10();
				GMT_Map avgMap = new GMT_Map(reg, avgShakeMap, spacing, cpt);
				avgMap.setUseGMTSmoothing(false);
				avgMap.setCustomScaleMin((double)cpt.getMinValue());
				avgMap.setCustomScaleMax((double)cpt.getMaxValue());
				avgMap.setCustomLabel("Log10 "+(float)periods[p]+"s Average of N="+num+" ShakeMaps");
				avgMap.setBlackBackground(false);
				FaultBasedMapGen.plotMap(outputDir, prefix+"_avg", false, avgMap);
				
				lnAvgShakeMap.log10();
				avgMap = new GMT_Map(reg, lnAvgShakeMap, spacing, cpt);
				avgMap.setUseGMTSmoothing(false);
				avgMap.setCustomScaleMin((double)cpt.getMinValue());
				avgMap.setCustomScaleMax((double)cpt.getMaxValue());
				avgMap.setCustomLabel("Log10 "+(float)periods[p]+"s Average (ln space) of N="+num+" ShakeMaps");
				avgMap.setBlackBackground(false);
				FaultBasedMapGen.plotMap(outputDir, prefix+"_lnavg", false, avgMap);
			}
		}
	}
	
	private static MinMaxAveTracker tracker(GeoDataSet geo) {
		MinMaxAveTracker track = new MinMaxAveTracker();
		for (int i=0; i<geo.size(); i++)
			track.addValue(geo.get(i));
		return track;
	}

}
