package scratch.kevin.nshm23.hazardValidation;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.Precision;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.WeightedList;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.erf.BaseFaultSystemSolutionERF;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.PointSurface;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.faultSurface.cache.SurfaceDistances;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.nshmp.shaded.model.NshmErf;

public class GridPropInvestigation {

	public static void main(String[] args) throws IOException {
		File inputSolFile = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip");
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
//		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0");
		
		double refDist = 20d;
		
		GriddedSeismicitySettings gridSettings = BaseFaultSystemSolutionERF.GRID_SETTINGS_DEFAULT;
		
		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		GridSourceList gridList = sol.requireModule(GridSourceList.class);
		
		TectonicRegionType trt = TectonicRegionType.STABLE_SHALLOW;
		ScalarIMR gmm = AttenRelRef.USGS_NSHM23_STABLE.get();
		
//		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
//		ScalarIMR gmm = AttenRelRef.USGS_NSHM23_ACTIVE.get();
		
		gmm.setParamDefaults();
		
		NshmpHazardModel model = NshmpHazardModel.load(modelDir.toPath());
		
		NshmErf wrapper = new NshmErf(model, Set.of(trt), IncludeBackgroundOption.ONLY);
		wrapper.getTimeSpan().setDuration(1d);
		wrapper.updateForecast();
		
//		Location testLoc = new Location(47, -108);
		Location testLoc = new Location(47.2, -107.1);
//		Location testLoc = new Location(47, -105);
//		Location testLoc = new Location(35, -118);
		int gridIndex = gridList.getLocationIndex(testLoc);
		testLoc = gridList.getLocation(gridIndex);
		System.out.println("Debug location: "+gridIndex+". "+testLoc);
		
		Location refLoc = LocationUtils.location(testLoc, 0, refDist);
		
		List<GriddedRupture> myRups = gridList.getRuptures(trt, gridIndex);
		ProbEqkSource mySource = gridList.getSource(trt, gridIndex, 1d, null, gridSettings);
		Preconditions.checkState(myRups.size() == mySource.getNumRuptures());
		Site fakeSite = new Site(testLoc);
		ProbEqkSource match = null;
		for (ProbEqkSource source : wrapper) {
			if (source.getTectonicRegionType() != trt)
				continue;
			if (source.getMinDistance(fakeSite) < 0.001) {
				Preconditions.checkState(match == null, "Multiple matches");
				match = source;
			}
		}
		
		Preconditions.checkNotNull(match);
		
		System.out.println("We have "+myRups.size()+" ruptures");
		System.out.println("Wrapper has "+match.getNumRuptures()+" ruptures");
		
		List<Integer> unmapped = new ArrayList<>(match.getNumRuptures());
		for (int i=0; i<match.getNumRuptures(); i++)
			unmapped.add(i);
		
		for (int g=0; g<myRups.size(); g++) {
			List<Integer> matches = new ArrayList<>();
			GriddedRupture rup = myRups.get(g);
			ProbEqkRupture myProbRup = mySource.getRupture(g);
			double mag = rup.properties.magnitude;
			double rake = rup.properties.rake;
			for (int i=unmapped.size(); --i>=0;) {
				int index = unmapped.get(i);
				ProbEqkRupture probRup = match.getRupture(index);
				if (Precision.equals(mag, probRup.getMag(), 0.001) && Precision.equals(rake, probRup.getAveRake(), 1)) {
					matches.add(index);
					unmapped.remove(i);
				}
			}
			Collections.reverse(matches);
			
			PointSurface.DistanceCorrectable surf = (PointSurface.DistanceCorrectable) myProbRup.getRuptureSurface();
			
			double rateSum = 0d;
			System.out.println("Ours:\t"+rupStr(rup));
			for (int l=0; l<2; l++) {
				Location loc = l == 0 ? testLoc : refLoc;
				String locStr = l == 0 ? "Colocated" : (int)refDist+" km away";
				WeightedList<SurfaceDistances> dists = surf.getCorrectedDistances(loc);
				for (int d=0; d<dists.size(); d++) {
					System.out.println("\t\t"+locStr+" dists: "+dists.getValue(d)+" (wt="+(float)dists.getWeight(d)+")");
					PointSurface.SiteSpecificDistanceCorrected corrSurf = new PointSurface.SiteSpecificDistanceCorrected(surf, loc, dists.getValue(d));
					ProbEqkRupture rupCopy = new ProbEqkRupture(rup.properties.magnitude, rup.properties.rake, 0d, corrSurf, null);
					System.out.println("\t\tGMM:\t"+gmmStr(gmm, loc, rupCopy));
				}
			}
			for (int i=0; i<matches.size(); i++) {
				int index = matches.get(i);
				ProbEqkRupture theirs = match.getRupture(index);
				if (i == 0)
					System.out.println("Theirs:\t"+rupStr(theirs));
				else
					System.out.println("\t"+rupStr(theirs));
				for (int l=0; l<2; l++) {
					Location loc = l == 0 ? testLoc : refLoc;
					String locStr = l == 0 ? "Colocated" : (int)refDist+" km away";
					System.out.println("\t\t"+locStr+" dists: "+theirs.getRuptureSurface().getDistances(loc));
					System.out.println("\t\tGMM:\t"+gmmStr(gmm, loc, theirs));
				}
				
				rateSum += theirs.getMeanAnnualRate(1d);
			}
			System.out.println("\tRate ratio = "+rateDF.format(rup.rate)+" / "+rateDF.format(rateSum)+" = "+(float)(rup.rate/rateSum));
			System.out.println();
		}
		
		if (!unmapped.isEmpty()) {
			System.out.println("Unmapped:");
			for (int index : unmapped)
				System.out.println("\t"+rupStr(match.getRupture(index)));
		}
	}
	
	private static String gmmStr(ScalarIMR gmm, Location loc, ProbEqkRupture rup) {
		Site site = new Site(loc);
		site.addParameterList(gmm.getSiteParams());
		gmm.setSite(site);
		gmm.setEqkRupture(rup);
		return "mean="+(float)gmm.getMean()+";\tsigma="+(float)gmm.getStdDev();
	}
	
	private static DecimalFormat rateDF = new DecimalFormat("0.00E0");
	
	private static String rupStr(ProbEqkRupture rup) {
		double rate = rup.getMeanAnnualRate(1d);
		RuptureSurface surf = rup.getRuptureSurface();
		return "M"+(float)rup.getMag()+";\trake="+(float)rup.getAveRake()+"\trate="+rateDF.format(rate)
				+"\tzTOR="+(float)surf.getAveRupTopDepth()+"\tzBOT="+(float)surf.getAveRupTopDepth()
				+"\tW="+widStr(surf)+"\tL="+lenStr(surf)+"\tdip="+(float)surf.getAveDip();
	}
	
	private static String widStr(RuptureSurface surf) {
		try {
			return (float)surf.getAveWidth()+"";
		} catch (Exception e) {
			return "Error";
		}
	}
	
	private static String lenStr(RuptureSurface surf) {
		try {
			return (float)surf.getAveLength()+"";
		} catch (Exception e) {
			return "Error";
		}
	}
	
	private static String rupStr(GriddedRupture rup) {
		return "M"+(float)rup.properties.magnitude+";\trake="+(float)rup.properties.rake+"\trate="+rateDF.format(rup.rate)
				+"\tzTOR="+(float)rup.properties.upperDepth+"\tzBOT="+(float)rup.properties.lowerDepth
				+"\tW="+(float)rup.properties.getDownDipWidth()+"\tL="+(float)rup.properties.length+"\tdip="+(float)rup.properties.dip;
	}

}
