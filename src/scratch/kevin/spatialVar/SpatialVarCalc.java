package scratch.kevin.spatialVar;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
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
import org.opensha.sha.gcim.calc.NearPD;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;
import com.google.common.base.Stopwatch;

import Jama.CholeskyDecomposition;
import Jama.Matrix;
import scratch.UCERF3.analysis.FaultBasedMapGen;

public class SpatialVarCalc {
	
	private double[] periods;
	private List<? extends Site> sites;
	
	private Matrix B1, B2, B3;
	private Matrix K1, K2, K3;
	private Matrix D1, D2;
	private Matrix L1, L2;

	public SpatialVarCalc(double[] periods, List<? extends Site> sites) {
		this.periods = periods;
		this.sites = sites;
		
		B1 = new Matrix(periods.length, periods.length);
		B2 = new Matrix(periods.length, periods.length);
		B3 = new Matrix(periods.length, periods.length);
		
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
		
		K1 = decompose(B1).getL();
		K2 = decompose(B2).getL();
		K3 = decompose(B3).getL();
		
		// coregionalization matrices
		D1 = new Matrix(sites.size(), sites.size());
		D2 = new Matrix(sites.size(), sites.size());
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
		Matrix R1 = normRandArray(periods.length, sites.size(), rng);
		Matrix R2 = normRandArray(periods.length, sites.size(), rng);
		Matrix R3 = normRandArray(periods.length, sites.size(), rng);
		
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
	
	private CholeskyDecomposition decompose(Matrix B) {
		CholeskyDecomposition chol = new CholeskyDecomposition(B);
		if (!chol.isSPD()) {
			NearPD nearPD = new NearPD();
			nearPD.setKeepDiag(true);
//			nearPd.setEigTol(1.e-6);
			boolean success = nearPD.calcNearPD(B);
			double normFrob = nearPD.getFrobNorm();
			if (!success) {
				throw new RuntimeException("Error: nearPD failed to converge, the correlation matrix maybe" +
						" significantly different from a PD matrix, check that the correlation equations" +
						"used are reasonable");
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
		Preconditions.checkState(S.getRowDimension() == periods.length);
		Preconditions.checkState(S.getColumnDimension() == sites.size());
		Preconditions.checkState(p >= 0 && p < periods.length);
		Preconditions.checkState(input.size() == sites.size());
		GeoDataSet ret = input.copy();
		for (int i=0; i<ret.size(); i++) {
			double lnOrig = Math.log(ret.get(i));
//			double lnRes = Math.log(S.get(p, i));
			double lnRes = S.get(p, i);
			ret.set(i, Math.exp(lnOrig+lnRes));
		}
		return ret;
	}
	
	public static void main(String[] args) throws IOException, GMT_MapException {
		File outputDir = new File("/tmp/spatial_var_test");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		Region reg = new Region(new Location(33.5, -119), new Location(35.5, -117));
		double spacing = 0.025d;
//		double spacing = 0.05d;
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
		
		double[] periods = { 1d };
		double sigma = 0.6d;
		int num = 1000;
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
