package scratch.kevin.nshm27.figures;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.reports.plots.NucleationRatePlot;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

public class CrustalNucleationAroundFaultZoom {

	public static void main(String[] args) throws IOException {
		Region reg = new Region(new Location(13, 144), new Location(16, 146));
		FaultSystemSolution sol = FaultSystemSolution.load(new File(
				"/data/kevin/nshm23/batch_inversions/2026_03_27-nshm26-GNMI-2000samples-gridded/"
				+ "results_GNMI_V1_ACTIVE_SHALLOW_branch_averaged.zip"));
		File outputDir = new File("/tmp");
		GridSourceList gridList = sol.requireModule(GridSourceList.class);
		
		GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);
		
		double[] minMags = {5d, 6d, 6.5, 7d, 7.5};
		GriddedGeoDataSet[] gridXYZs = new GriddedGeoDataSet[minMags.length];
		GriddedGeoDataSet[] totXYZs = new GriddedGeoDataSet[minMags.length];
		for (int m=0; m<minMags.length; m++) {
			gridXYZs[m] = new GriddedGeoDataSet(gridReg);
			totXYZs[m] = new GriddedGeoDataSet(gridReg);
		}
		
		FaultGridAssociations assoc = sol.getRupSet().requireModule(FaultGridAssociations.class);
		
		for (int l=0; l<gridReg.getNodeCount(); l++) {
			Location loc = gridReg.getLocation(l);
			int index = gridList.getLocationIndex(loc);
			for (GriddedRupture rup : gridList.getRuptures(TectonicRegionType.ACTIVE_SHALLOW, index)) {
				for (int m=0; m<minMags.length; m++) {
					if (rup.properties.magnitude >= minMags[m]) {
						gridXYZs[m].add(l, rup.rate);
						totXYZs[m].add(l, rup.rate);
					}
				}
			}
		}
		
		List<IncrementalMagFreqDist> solNuclMFDs = NucleationRatePlot.calcNuclMFDs(sol, TectonicRegionType.ACTIVE_SHALLOW);
		for (int m=0; m<minMags.length; m++) {
			GriddedGeoDataSet faultXYZ = NucleationRatePlot.calcFaultNucleationRates(assoc.getRegion(), sol, assoc, solNuclMFDs, minMags[m]);
			for (int l=0; l<gridReg.getNodeCount(); l++) {
				Location loc = gridReg.getLocation(l);
				int index = faultXYZ.indexOf(loc);
				totXYZs[m].add(l, faultXYZ.get(index));
			}
		}
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
		mapMaker.setFaultSections(sol.getRupSet().getFaultSectionDataList());
		
		CPT cpt = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-8, -3);
		cpt.setLog10(true);
		
		DecimalFormat magDF = new DecimalFormat("0.#");
		for (int m=0; m<minMags.length; m++) {
			mapMaker.plotXYZData(gridXYZs[m], cpt, "Crustal gridded nucleation rate (M>"+magDF.format(minMags[m])+")");
			mapMaker.plot(outputDir, "nshm27_crustal_grid_nucl_zoom_m"+magDF.format(minMags[m]), " ");
			mapMaker.plotXYZData(totXYZs[m], cpt, "Crustal total nucleation rate (M>"+magDF.format(minMags[m])+")");
			mapMaker.plot(outputDir, "nshm27_crustal_nucl_zoom_m"+magDF.format(minMags[m]), " ");
		}
	}

}
