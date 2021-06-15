package scratch.kevin.ucerf3.eal.spatialCorr;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.utm.UTM;
import org.opensha.commons.geo.utm.WGS84;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.nshmp2.erf.source.PointSource;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.rupForecastImpl.PointSource13b;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_WrapperFullParam;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sra.calc.EALCalculator;
import org.opensha.sra.gui.portfolioeal.Asset;
import org.opensha.sra.gui.portfolioeal.Portfolio;
import org.opensha.sra.vulnerability.Vulnerability;

import com.google.common.base.Preconditions;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import com.google.common.collect.Table.Cell;

import scratch.UCERF3.utils.FaultSystemIO;

public class SpatiallyCorrelatedLossCalc {
	
	protected static boolean D = false;
	
	/**
	 * Calculates centroid of rupture surface in UTM coordinates, and converts back to 
	 * WGS84. All surface points weighted equally.
	 * @param surf
	 * @return
	 */
	public static Location calcRupCentroid(RuptureSurface surf) {
		if (surf instanceof PointSource)
			return ((PointSource)surf).getLocation();
		if (surf instanceof PointSource13b)
			return ((PointSource13b)surf).getLocation();
		double totWeight = 0d;
		double northing = 0d;
		double easting = 0d;
		
		int zone = -1;
		char letter = 'Z';
		for (Location loc : surf.getEvenlyDiscritizedListOfLocsOnSurface()) {
			double weight = 1d;
			totWeight += weight;
			WGS84 wgs = new WGS84(loc.getLatitude(), loc.getLongitude());
			UTM utm;
			if (zone < 0) {
				utm = new UTM(wgs);
				zone = utm.getZone();
				letter = utm.getLetter();
			} else {
				utm = new UTM(wgs, zone, letter);
			}
			northing += weight*utm.getNorthing();
			easting += weight*utm.getEasting();
		}
		northing /= totWeight;
		easting /= totWeight;
		
		UTM utm = new UTM(zone, letter, easting, northing);
		WGS84 wgs = new WGS84(utm);
		double lat = wgs.getLatitude();
		double lon = wgs.getLongitude();
		return new Location(lat, lon);
	}
	
	public static Table<Double, RandomFieldLoader, Double> calcSpatiallyCorrelatedLoss(
			ScalarIMR gmpe, List<Asset> assets, EqkRupture rup, Location rupCentroid,
			double[] betweenEventStdDevs, RandomFieldLoader[] randFields) {
		return calcSpatiallyCorrelatedLoss(gmpe, assets, rup, rupCentroid,
				betweenEventStdDevs, randFields, null);
	}
	
	public static Table<Double, RandomFieldLoader, Double> calcSpatiallyCorrelatedLoss(
			ScalarIMR gmpe, List<Asset> assets, EqkRupture rup, Location rupCentroid,
			double[] betweenEventStdDevs, RandomFieldLoader[] randFields, DiscretizedFunc magThreshFunc) {
		Table<Double, RandomFieldLoader, Double> ret = HashBasedTable.create();
		
		gmpe.setEqkRupture(rup);
		for (Asset asset : assets) {
			Site site = new Site();
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>)param.clone());
			asset.siteSetup(site);
			site = asset.getSite();
			
			if (magThreshFunc != null) {
				double distance = rup.getRuptureSurface().getQuickDistance(site.getLocation());
				if (distance > magThreshFunc.getMaxX()) {
					if (D) System.out.println("Skipping an asset above max dist (dist="+(float)distance+")");
					continue;
				}
				double magThresh = magThreshFunc.getInterpolatedY(distance);
				
				if (rup.getMag() < magThresh) {
					if (D) System.out.println("Skipping due to mag-dist thresh fail (mag="+(float)rup.getMag()
						+", thresh="+(float)magThresh+", dist="+(float)distance+")");
					continue;
				}
			}
			
			Vulnerability vulnModel;
			try {
				vulnModel = asset.getVulnModel();
			} catch (ClassNotFoundException | InstantiationException | IllegalAccessException e) {
				throw ExceptionUtils.asRuntimeException(e);
			}
			
			gmpe.setSite(asset.getSite());
			
			String imt = vulnModel.getIMT();
			gmpe.setIntensityMeasure(imt);
			if (imt.equals(SA_Param.NAME))
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), vulnModel.getPeriod());
			
			double mean, phi, tau;
			if (gmpe instanceof NGAW2_WrapperFullParam) {
				ScalarGroundMotion gm = ((NGAW2_WrapperFullParam)gmpe).getGroundMotion();
				mean = gm.mean();
				phi = gm.phi();
				tau = gm.tau();
			} else {
				mean = gmpe.getMean();
				StdDevTypeParam type = (StdDevTypeParam)(gmpe.getParameter(StdDevTypeParam.NAME));
				type.setValue(StdDevTypeParam.STD_DEV_TYPE_INTER);
				tau = gmpe.getStdDev();
				type.setValue(StdDevTypeParam.STD_DEV_TYPE_INTRA);
				phi = gmpe.getStdDev();
				type.setValue(StdDevTypeParam.STD_DEV_TYPE_TOTAL);
			}
			
			if (D) System.out.println("Calculating for asset "+asset.getVulnModelName()+" with value "+asset.getValue());
			if (D) System.out.println("mean="+mean+", phi="+phi+", tau="+tau);
			
			DiscretizedFunc vulnFunc = vulnModel.getVulnerabilityFunc();
			
			for (double betweenEventStdDev : betweenEventStdDevs) {
				// ground motion considering the between-event term
				double tauGM = mean + betweenEventStdDev*tau;
				
				if (D) System.out.println("\ttauGM="+tauGM);
				
				for (RandomFieldLoader field : randFields) {
					// ground motion considering the randomly sampled within-event term (and also tau)
					double phiGM = tauGM + field.getValue(site.getLocation(), rupCentroid)*phi;
					if (D) System.out.println("\t\tphiGM="+phiGM);
					
					Preconditions.checkState(Double.isFinite(phiGM),
							"bad phiGM=%s for mean=%s, phi=%s, tau=%s, tauGM=%s",
							phiGM, mean, phi, tau, tauGM);
					
					// set MAFE. 1 if below GM, 0 if above
					double linearGM = Math.exp(phiGM);
					if (D) System.out.println("\t\tlinearGM="+linearGM);
//					if (D) System.out.println("\t\tMAFE");
//					Preconditions.checkState(linearGM > 0d, "bad linearGM=%s for phiGM=%s", linearGM, phiGM);
//					for (int i=0; i<mafe.size(); i++) {
//						if (linearGM >= mafe.getX(i))
//							mafe.set(i, 1d);
//						else
//							mafe.set(i, 0d);
//						if (D) System.out.println("\t\t\t"+(float)mafe.getX(i)+"="+(float)mafe.getY(i)
//							+"\t(D="+ealCalc.getDF().get(i).floatValue()+")");
//					}
//					ealCalc.setMAFE(mafe);
//					
//					// expected loss conditioned on occurence of this rup with this GM value
//					double el = ealCalc.computeEAL();
					double interpDF;
					if (linearGM < vulnFunc.getMinX())
						interpDF = 0d;
					else if (linearGM > vulnFunc.getMaxX())
						interpDF = vulnFunc.getMaxY();
					else
						interpDF = vulnFunc.getInterpolatedY(linearGM);
					if (D) System.out.println("\t\tInterp DF: "+interpDF);
					double el = asset.getValue()*interpDF;
					if (D) System.out.println("\t\tLoss: "+el);
					
					Double prevLoss = ret.get(betweenEventStdDev, field);
					if (prevLoss == null)
						prevLoss = 0d;
					ret.put(betweenEventStdDev, field, prevLoss+el);
				}
			}
		}
		
		return ret;
	}
	
	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution fss = FaultSystemIO.loadSol(
				new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip"));
		
		int num = 10;
		Random r = new Random();
		
		FaultSystemRupSet rupSet = fss.getRupSet();
		
		double[] betweens = { 0d };
		
		File fieldDir = new File("/home/kevin/OpenSHA/UCERF3/eal/random_fields/sa10_1km_800x800");
		RandomFieldLoader[] fields = {
				RandomFieldLoader.load(new File(fieldDir, "800x800SA10_001.csv"), 1d)
		};
		
		ScalarIMR gmpe = AttenRelRef.ASK_2014.instance(null);
		gmpe.setParamDefaults();
		
		Site site = new Site(new Location(34, -118));
		site.addParameterList(gmpe.getSiteParams());
		
		File portfolioFile = new File("/home/kevin/OpenSHA/UCERF3/eal/"
				+ "Porter-09-Feb-2020-CEA-100-pct-procy-portfolio-wills2015.csv");
		Portfolio portfolio = Portfolio.createPortfolio(portfolioFile);
		List<Asset> assets = portfolio.getAssetList().subList(0, 5);
		
		SpatiallyCorrelatedLossCalc.D = true;
		RandomFieldLoader.D = true;
		
		for (int i=0; i<num; i++) {
			int rupIndex = r.nextInt(rupSet.getNumRuptures());
			
			double mag = rupSet.getMagForRup(rupIndex);
			double rate = rupSet.getAveRakeForRup(rupIndex);
			RuptureSurface surf = rupSet.getSurfaceForRupture(rupIndex, 1d);
			
			EqkRupture rup = new EqkRupture(mag, rate, surf, null);
			
			Location centroid = calcRupCentroid(surf);
			
			Table<Double, RandomFieldLoader, Double> vals = calcSpatiallyCorrelatedLoss(
					gmpe, assets, rup, centroid, betweens, fields);
			
			for (Cell<Double, RandomFieldLoader, Double> cell : vals.cellSet()) {
				System.out.println("Value: "+cell.getValue());
			}
		}
	}

}
