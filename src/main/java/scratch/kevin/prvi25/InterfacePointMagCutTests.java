package scratch.kevin.prvi25;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.DataUtils.MinMaxAveTracker;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.util.GriddedSeismicitySettings;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.attenRelImpl.nshmp.NSHMP_GMM_Wrapper;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import gov.usgs.earthquake.nshmp.gmm.Gmm;

public class InterfacePointMagCutTests {

	public static void main(String[] args) throws IOException {
		FaultSystemSolution sol = FaultSystemSolution.load(new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_11_19-prvi25_subduction_branches/results_PRVI_INTERFACE_ONLY_branch_averaged_gridded.zip"));
		
		GridSourceProvider gridProv = sol.getGridSourceProvider();
		
		GriddedSeismicitySettings fullFiniteSettings = GriddedSeismicitySettings.DEFAULT.forPointSourceMagCutoff(5d);
		GriddedSeismicitySettings magCut6FiniteSettings = GriddedSeismicitySettings.DEFAULT.forPointSourceMagCutoff(6d);
		
		System.out.println("Full finite settings: "+fullFiniteSettings);
		System.out.println("Mag cut settings: "+magCut6FiniteSettings);
		
//		double maxMag = 8.55;
//		double maxMag = 7d;
		double maxMag = 6d;
		EvenlyDiscretizedFunc refMFD = FaultSysTools.initEmptyMFD(5.01, maxMag);
		Quantity[] quants = Quantity.values();
		MinMaxAveTracker[][] fullTracks = new MinMaxAveTracker[quants.length][refMFD.size()];
		MinMaxAveTracker[][] cutTracks = new MinMaxAveTracker[quants.length][refMFD.size()];
		MinMaxAveTracker[][] diffTracks = new MinMaxAveTracker[quants.length][refMFD.size()];
		MinMaxAveTracker[][] pDiffTracks = new MinMaxAveTracker[quants.length][refMFD.size()];
		
//		ScalarIMR gmm = AttenRelRef.USGS_PRVI_INTERFACE.get();
		
//		Gmm gmmRef = Gmm.AG_20_GLOBAL_INTERFACE;
//		Gmm gmmRef = Gmm.KBCG_20_GLOBAL_INTERFACE;
		Gmm gmmRef = Gmm.PSBAH_20_GLOBAL_INTERFACE;
		ScalarIMR gmm = new NSHMP_GMM_Wrapper.Single(gmmRef, gmmRef.toString(), gmmRef.name(), false, null);
		gmm.setIntensityMeasure(PGA_Param.NAME);
		
		for (int q=0; q<quants.length; q++) {
			for (int i=0; i<refMFD.size(); i++) {
				fullTracks[q][i] = new MinMaxAveTracker();
				cutTracks[q][i] = new MinMaxAveTracker();
				diffTracks[q][i] = new MinMaxAveTracker();
				pDiffTracks[q][i] = new MinMaxAveTracker();
			}
		}
		
		Location testLoc = new Location(18.25, -66.5);
		
		Site site = new Site(testLoc);
		site.addParameterList(gmm.getSiteParams());
		gmm.setSite(site);
		
		List<Integer> locs;
//		locs = new ArrayList<>(gridProv.getNumLocations());
//		for (int l=0; l<gridProv.getNumLocations(); l++)
//			locs.add(l);
		locs = List.of(gridProv.getGriddedRegion().indexForLocation(new Location(19, -65.5)));
		
		for (int l : locs) {
			ProbEqkSource sourceFull = gridProv.getSource(TectonicRegionType.SUBDUCTION_INTERFACE, l, 1d, null, fullFiniteSettings);
			if (sourceFull == null)
				continue;
			ProbEqkSource sourceCut = gridProv.getSource(TectonicRegionType.SUBDUCTION_INTERFACE, l, 1d, null, magCut6FiniteSettings);

			List<List<ProbEqkRupture>> fullRupsBinned = new ArrayList<>(refMFD.size());
			List<List<ProbEqkRupture>> cutRupsBinned = new ArrayList<>(refMFD.size());
			for (int i=0; i<refMFD.size(); i++) {
				fullRupsBinned.add(new ArrayList<>());
				cutRupsBinned.add(new ArrayList<>());
			}
			double rateSumFull = 0d;
			for (ProbEqkRupture rup : sourceFull) {
				if (rup.getMag() > maxMag)
					continue;
				fullRupsBinned.get(refMFD.getClosestXIndex(rup.getMag())).add(rup);
				rateSumFull += rup.getMeanAnnualRate(1d);
			}
			double rateSumCut = 0d;
			for (ProbEqkRupture rup : sourceCut) {
				if (rup.getMag() > maxMag)
					continue;
				cutRupsBinned.get(refMFD.getClosestXIndex(rup.getMag())).add(rup);
				rateSumCut += rup.getMeanAnnualRate(1d);
			}
			Preconditions.checkState((float)rateSumFull == (float)rateSumCut, "Rate mismatch: %s != %s", (float)rateSumFull, (float)rateSumCut);
			
			for (int i=0; i<refMFD.size(); i++) {
				List<ProbEqkRupture> rupsFull = fullRupsBinned.get(i);
				List<ProbEqkRupture> rupsCut = cutRupsBinned.get(i);
				if (rupsFull.isEmpty() || rupsCut.isEmpty()) {
					Preconditions.checkArgument(rupsFull.isEmpty() && rupsCut.isEmpty());
					continue;
				}
				
				for (int q=0; q<quants.length; q++) {
					List<Double> fullVals = new ArrayList<>();
					for (ProbEqkRupture rup : rupsFull) {
						double val = quants[q].calc(rup, testLoc, gmm);
						fullVals.add(val);
						fullTracks[q][i].addValue(val);
					}
					List<Double> cutVals = new ArrayList<>();
					for (ProbEqkRupture rup : rupsCut) {
						double val = quants[q].calc(rup, testLoc, gmm);
						cutVals.add(val);
						cutTracks[q][i].addValue(val);
					}
					
					for (double fullVal : fullVals) {
						for (double cutVal : cutVals) {
							double diff = cutVal - fullVal;
							double pDiff = 100d*diff/fullVal;
							
							diffTracks[q][i].addValue(diff);
							pDiffTracks[q][i].addValue(pDiff);
						}
					}
				}
			}
		}
		
		for (int m=0; m<refMFD.size(); m++) {
			if (fullTracks[0][m].getNum() == 0)
				continue;
			System.out.println("M"+(float)refMFD.getX(m));
			for (int q=0; q<quants.length; q++) {
				System.out.println("\t"+quants[q]);
				System.out.println("\t\tFull:\t"+fullTracks[q][m]);
				System.out.println("\t\tCut:\t"+cutTracks[q][m]);
				System.out.println("\t\tDiff:\t"+diffTracks[q][m]);
				System.out.println("\t\t% Diff:\t"+pDiffTracks[q][m]);
			}
		}
	}
	
	private enum Quantity {
		RRUP {
			@Override
			public double calc(ProbEqkRupture rup, Location loc, ScalarIMR gm) {
				return rup.getRuptureSurface().getDistanceRup(loc);
			}
		},
		RJB {
			@Override
			public double calc(ProbEqkRupture rup, Location loc, ScalarIMR gm) {
				return rup.getRuptureSurface().getDistanceJB(loc);
			}
		},
		RX {
			@Override
			public double calc(ProbEqkRupture rup, Location loc, ScalarIMR gm) {
				return rup.getRuptureSurface().getDistanceX(loc);
			}
		},
		ZTOR {
			@Override
			public double calc(ProbEqkRupture rup, Location loc, ScalarIMR gm) {
				return rup.getRuptureSurface().getAveRupTopDepth();
			}
		},
		GM {
			@Override
			public double calc(ProbEqkRupture rup, Location loc, ScalarIMR gm) {
				gm.setEqkRupture(rup);
				return Math.exp(gm.getMean());
			}
		};
		
		public abstract double calc(ProbEqkRupture rup, Location loc, ScalarIMR gm);
	}

}
