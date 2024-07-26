package scratch.kevin.pointSources;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.calc.magScalingRelations.magScalingRelImpl.WC1994_MagLengthRelationship;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.json.Feature;
import org.opensha.commons.util.DataUtils;
import org.opensha.commons.util.modules.ModuleArchive;
import org.opensha.commons.util.modules.OpenSHA_Module;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.FiniteRuptureConverter;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.faultSysSolution.modules.MFDGridSourceProvider;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.FocalMech;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

public class TestFullGridSourceFileWriter {
	
	private static WC1994_MagLengthRelationship WC94 = new WC1994_MagLengthRelationship();
	
	public static class NSHM23_WUS_FiniteRuptureConverter implements FiniteRuptureConverter {

		@Override
		public GriddedRupture buildFiniteRupture(int gridIndex, Location loc, double magnitude, double rate,
				FocalMech focalMech, int[] associatedSections, double[] associatedSectionFracts) {
			// TODO Auto-generated method stub
			
			double dipRad = Math.toRadians(focalMech.dip());
			
			double depth = (float)magnitude < 6.5f ? 5d : 1d;
			double length = WC94.getMedianLength(magnitude);
			double aspectWidth = length / 1.5;
			double ddWidth = (14.0 - depth) / Math.sin(dipRad);
			ddWidth = Math.min(aspectWidth, ddWidth);
			double lower = depth + ddWidth * Math.sin(dipRad);
			
			return new GriddedRupture(gridIndex, loc, magnitude, rate, focalMech.rake(), focalMech.dip(), Double.NaN,
					null, depth, lower, length, TectonicRegionType.ACTIVE_SHALLOW, associatedSections, associatedSectionFracts);
		}
		
	}

	public static void main(String[] args) throws IOException {
		File dir = new File("/home/kevin/OpenSHA/nshm23/batch_inversions/"
				+ "2024_02_02-nshm23_branches-WUS_FM_v3");
		File inputFile = new File(dir, "results_WUS_FM_v3_branch_averaged_gridded.zip");
		File outputFile = new File(dir, "results_WUS_FM_v3_branch_averaged_mod_gridded.zip");
		FaultSystemSolution sol = FaultSystemSolution.load(inputFile);
		MFDGridSourceProvider gridProv = sol.requireModule(MFDGridSourceProvider.class);
		FaultGridAssociations associations = sol.getRupSet().requireModule(FaultGridAssociations.class);
		
		File outputDir = new File("/tmp/grid_source_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		double totRateOrig = 0d;
		int origNumNonNullMFDs = 0;
		for (int gridIndex=0; gridIndex<gridProv.getNumLocations(); gridIndex++) {
			IncrementalMagFreqDist mfd = gridProv.getMFD(gridIndex);
			if (mfd != null) {
				origNumNonNullMFDs++;
				totRateOrig += mfd.calcSumOfY_Vals();
			}
		}
		
		GridSourceList gridSources = GridSourceList.convert(gridProv, associations, new NSHM23_WUS_FiniteRuptureConverter());
		double totRateMod = 0d;
		int rupCount = 0;
		int mostAssociations = 0;
		GriddedRupture rupWithMostAssociations = null;
		for (TectonicRegionType trt : gridSources.getTectonicRegionTypes()) {
			for (int gridIndex=0; gridIndex<gridSources.getNumLocations(); gridIndex++) {
				for (GriddedRupture rup : gridSources.getRuptures(trt, gridIndex)) {
					totRateMod += rup.rate;
					rupCount++;
					if (rup.associatedSections != null && rup.associatedSections.length > mostAssociations) {
						mostAssociations = rup.associatedSections.length;
						rupWithMostAssociations = rup;
					}
				}
			}
		}
		System.out.println("******************************");
		System.out.println("Originally had "+gridProv.getNumSources()+" sources ("+origNumNonNullMFDs+" non-null) with totRate="+(float)totRateOrig);
		System.out.println("Now have "+rupCount+" rups across "+gridSources.getNumSources()+" sources with totRate="+(float)totRateMod);
		if (rupWithMostAssociations != null)
			System.out.println("Most associations: "+mostAssociations+" (gridIndex="+rupWithMostAssociations.gridIndex
					+", M="+(float)rupWithMostAssociations.magnitude);
		System.out.println("******************************");
		
		ModuleArchive<OpenSHA_Module> archive = new ModuleArchive<>();
		archive.addModule(gridSources);
		archive.write(new File(dir, "mod_gridded.zip"));
		
		sol.setGridSourceProvider(gridSources);
		
		sol.write(outputFile);
	}

}
