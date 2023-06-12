package scratch.kevin.spatialVar;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.opensha.commons.data.Site;
import org.opensha.commons.geo.Location;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.imr.attenRelImpl.ngaw2.NGAW2_Wrappers.ASK_2014_Wrapper;
import org.opensha.sha.imr.attenRelImpl.ngaw2.ScalarGroundMotion;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import Jama.Matrix;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;
import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO.ETAS_Catalog;

public class U3CatalogVarCalcDemo {

	public static void main(String[] args) throws IOException {
		ETAS_Catalog catalog = ETAS_CatalogIO.loadCatalog(new File("/tmp/catalog_test/output_dir/sampledEventsData.txt"));
		FaultSystemSolution fss = FaultSystemSolution.load(new File("/tmp/catalog_test/FM3_1_branch_averaged.zip"));
		FaultSystemRupSet rupSet = fss.getRupSet();
		
		// this uses ASK (2014)
		ASK_2014_Wrapper gmm = new ASK_2014_Wrapper();
		gmm.setParamDefaults();
		
		double period = 1d;
		// set to SA
		gmm.setIntensityMeasure(SA_Param.NAME);
		// set period
		SA_Param.setPeriodInSA_Param(gmm.getIntensityMeasure(), period);
		
		List<Site> sites = new ArrayList<>();
		
		sites.add(new Site(new Location(34, -118)));
		sites.add(new Site(new Location(35, -118)));
		sites.add(new Site(new Location(35, -117)));
		
		// sites need to have site parameters attached, this sets all of them to GMM default values
		for (Site site : sites) {
			for (Parameter<?> param : gmm.getSiteParams()) {
				param = (Parameter<?>)param.clone();
				site.addParameter(param);
			}
		}
		
		// now lets set Vs30 for one of them to a custom value
		sites.get(0).getParameter(Double.class, Vs30_Param.NAME).setValue(500d);
		
		// print them out so you see what the parameters actually are
		System.out.println("Site list:");
		for (int i=0; i<sites.size(); i++) {
			Site site = sites.get(i);
			System.out.println("Site i: "+site.getLocation());
			for (Parameter<?> param : site)
				System.out.println("\t"+param.getName()+":\t"+param.getValue());
		}
		
		System.out.println("Initializing spatial correlation matrices");
		SpatialVarCalc varCalc = new SpatialVarCalc(new double[] {period}, sites);
		
		// compute random fields for each rupture. set phi to 1 here, we'll rescale to actual event phi later
		System.out.println("Computing "+catalog.size()+" random fields");
		Random rand = new Random();
		Matrix[] fields = varCalc.computeRandWithinEventResiduals(rand, 1d, catalog.size());
		
		System.out.println("Calculating ground motions for catalog");
		for (int i=0; i<catalog.size(); i++) {
			ETAS_EqkRupture rup = catalog.get(i);
			System.out.println("\nRupture "+rup.getFSSIndex()+" (M"+(float)rup.getMag()+")");
			
//			ScalarGroundMotion[][] gms = new ScalarGroundMotion[sites.size()][periods.length];
			
			// attach the surface to the rupture
			if (rup.getFSSIndex() >= 0) {
				rup.setRuptureSurface(rupSet.getSurfaceForRupture(rup.getFSSIndex(), 1d));
				rup.setAveRake(rupSet.getAveRakeForRup(rup.getFSSIndex()));
			} else {
				rup.setPointSurface(rup.getHypocenterLocation());
				rup.setAveRake(0d);
			}
			// set the rupture in the GMM
			gmm.setEqkRupture(rup);
			
			ScalarGroundMotion[] gms = new ScalarGroundMotion[sites.size()];
			
			double avgPhi = 0d;
			double avgTau = 0d;
			for (int s=0; s<sites.size(); s++) {
				// set the site in the GMM
				gmm.setSite(sites.get(s));
				
				gms[s] = gmm.getGroundMotion();
				
				avgPhi += gms[s].phi();
				avgTau += gms[s].tau();
			}
			avgPhi /= (double)sites.size();
			avgTau /= (double)sites.size();
			
			System.out.println("\tphi="+(float)avgPhi+";\ttau="+(float)avgTau);
			
			// draw a random between-event residual using average tau
			double betweenResidual = rand.nextGaussian()*avgTau;
			System.out.println("\tRandom between-event residual: "+(float)betweenResidual);
			
			for (int s=0; s<sites.size(); s++) {
				// all in natural log units
				double origMean = gms[s].mean();
				// add between-event
				double randGM = origMean + betweenResidual;
				// add within-event
				// 0 here is the period index, but just 1 period for this example
				double withinResidual = fields[i].get(0, s)*avgPhi; // multiply by phi here so that this is an actual within event residual
				randGM += withinResidual;
				System.out.println("\tSite "+s+":\tmeanGM="+(float)origMean+";\trandWithin="+(float)withinResidual);
				System.out.println("\t\tRandom GM: "+(float)origMean+" + "+(float)betweenResidual
						+" + "+(float)withinResidual+" = "+(float)randGM);
			}
		}
	}

}
