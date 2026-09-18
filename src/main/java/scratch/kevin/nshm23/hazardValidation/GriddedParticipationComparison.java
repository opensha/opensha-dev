package scratch.kevin.nshm23.hazardValidation;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.BitSet;
import java.util.Set;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.FaultGridAssociations;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceList.GriddedRupture;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;

import org.opensha.nshmp.shaded.model.NshmpGridSource;
import org.opensha.nshmp.shaded.model.NshmpHazardModel;
import org.opensha.nshmp.shaded.model.NshmErf;
import org.opensha.nshmp.shaded.model.NshmSource;

public class GriddedParticipationComparison {

	public static void main(String[] args) throws IOException {
		File inputSolFile = new File("/home/kevin/OpenSHA/fss_inversions/2024_02_02-nshm23_branches-WUS_FM_v3/"
				+ "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip");
		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.1.3");
//		File modelDir = new File("/data/kevin/nshm23/nshmp-haz-models/nshm-conus-6.2.0");
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();

		TectonicRegionType trt = TectonicRegionType.ACTIVE_SHALLOW;
		IncludeBackgroundOption bgOp = IncludeBackgroundOption.EXCLUDE;
		double[] minMags = {0d, 6d, 7d};

		File outputDir = new File("/data/kevin/nshm23/nshmp-haz-models/gridded_partic_debug");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdirs(),
				"Could not create output directory: %s", outputDir.getAbsolutePath());

		GriddedRegion gridReg = new GriddedRegion(reg, 0.1, GriddedRegion.ANCHOR_0_0);

		FaultSystemSolution sol = FaultSystemSolution.load(inputSolFile);
		// we want direct mappings
		sol.getRupSet().removeModuleInstances(FaultGridAssociations.class);

		NshmpHazardModel model = NshmpHazardModel.load(modelDir.toPath());

		GriddedGeoDataSet[] fssXYZs = new GriddedGeoDataSet[minMags.length];
		GriddedGeoDataSet[] nhXYZs = new GriddedGeoDataSet[minMags.length];
		for (int m=0; m<minMags.length; m++) {
			fssXYZs[m] = new GriddedGeoDataSet(gridReg);
			nhXYZs[m] = new GriddedGeoDataSet(gridReg);
		}

		// do FSS
		if (bgOp == IncludeBackgroundOption.INCLUDE || bgOp == IncludeBackgroundOption.EXCLUDE) {
			System.out.println("Calculating for FSS on-fault");
			FaultSystemRupSet rupSet = sol.getRupSet();
			BitSet[] sectBitSets = new BitSet[rupSet.getNumSections()];
			for (int s=0; s<rupSet.getNumSections(); s++) {
				FaultSection sect = rupSet.getFaultSectionData(s);
				sectBitSets[s] = new BitSet(gridReg.getNodeCount());
				for (Location loc : sect.getFaultSurface(1d).getEvenlyDiscritizedListOfLocsOnSurface()) {
					int node = gridReg.indexForLocation(loc);
					if (node >= 0)
						sectBitSets[s].set(node);
				}
			}

			for (int r=0; r<rupSet.getNumRuptures(); r++) {
				BitSet rupBits = null;
				for (int s : rupSet.getSectionsIndicesForRup(r)) {
					if (rupBits == null)
						rupBits = (BitSet)sectBitSets[s].clone();
					else
						rupBits.or(sectBitSets[s]);
				}
				double rate = sol.getRateForRup(r);
				double mag = rupSet.getMagForRup(r);
				for (int m=0; m<minMags.length; m++)
					if ((float)mag >= (float)minMags[m])
						for (int i = rupBits.nextSetBit(0); i >= 0; i = rupBits.nextSetBit(i + 1))
							fssXYZs[m].add(i, rate);
			}
		}
		if (bgOp == IncludeBackgroundOption.INCLUDE || bgOp == IncludeBackgroundOption.ONLY) {
			System.out.println("Calculating for FSS gridded");
			GridSourceList gridList = sol.requireModule(GridSourceList.class);

			for (int l=0; l<gridList.getNumLocations(); l++) {
				Location loc = gridList.getLocation(l);
				int locIndex = gridReg.indexForLocation(loc);
				if (locIndex >= 0) {
					for (GriddedRupture rup : gridList.getRuptures(trt, l)) {
						for (int m=0; m<minMags.length; m++) {
							if ((float)rup.properties.magnitude >= (float)minMags[m])
								fssXYZs[m].add(locIndex, rup.rate);
						}
					}
				}
			}
		}

		// do NSHMP-haz
		NshmErf erf = new NshmErf(model, Set.of(trt), IncludeBackgroundOption.INCLUDE);
		erf.getTimeSpan().setDuration(1d);
		erf.updateForecast();
		System.out.println("Calculating for NSHMP-haz");
		int numNSHMPSources = 0;
		int numNSHMPRups = 0;
		for (ProbEqkSource source : erf) {
			NshmSource nshmSource = (NshmSource)source;
			Object delegate = nshmSource.delegate();
			boolean gridLike = delegate instanceof NshmpGridSource;
			if (!includeForBackgroundOption(gridLike, bgOp))
				continue;
			int gridSourceIndex = -1;
			if (gridLike)
				gridSourceIndex = gridReg.indexForLocation(toOpenSHALocation(((NshmpGridSource)delegate).location(null)));
			numNSHMPSources++;
			for (ProbEqkRupture rup : source) {
				numNSHMPRups++;
				double rate = rup.getMeanAnnualRate(1d);
				for (int m=0; m<minMags.length; m++) {
					if ((float)rup.getMag() >= (float)minMags[m]) {
						if (gridLike) {
							if (gridSourceIndex >= 0)
								nhXYZs[m].add(gridSourceIndex, rate);
						} else {
							addFiniteRupRate(nhXYZs[m], rup.getRuptureSurface(), rate);
						}
					}
				}
			}
		}
		
		System.out.println("NSHMP-Haz had "+numNSHMPSources+" matching sources with "+numNSHMPRups+" total ruptures");
		
		System.out.println("Building maps");

		GeographicMapMaker mapMaker = new GeographicMapMaker(reg);
		mapMaker.setFaultSections(sol.getRupSet().getFaultSectionDataList());
		mapMaker.setSectOutlineChar(null);
		Color transparent = new Color(255, 255, 255, 0);

		for (int m=0; m<minMags.length; m++) {
			String magLabel = magLabel(minMags[m]);
			System.out.println(magLabel+": OpenSHA="+(float)fssXYZs[m].getSumZ()
					+", NSHMP-Haz="+(float)nhXYZs[m].getSumZ());
			
			String prefix = trt.name()+"_"+bgOp.name()+"_"+magLabel;

			GriddedGeoDataSet diff = diff(fssXYZs[m], nhXYZs[m]);
			GriddedGeoDataSet pDiff = pDiff(fssXYZs[m], nhXYZs[m]);

			writeCSV(new File(outputDir, prefix+"_gridded_partic.csv"), fssXYZs[m], nhXYZs[m], diff, pDiff);

			double maxRate = Math.max(max(fssXYZs[m]), max(nhXYZs[m]));
			if (maxRate > 0d) {
				double logMax = Math.ceil(Math.log10(maxRate));
				CPT rateCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(logMax - 6d, logMax);
				rateCPT.setNanColor(transparent);

				mapMaker.plotXYZData(asLog10(fssXYZs[m]), rateCPT, "OpenSHA Rate, "+magLabel);
				mapMaker.plot(outputDir, prefix+"_opensha", " ");
				mapMaker.plotXYZData(asLog10(nhXYZs[m]), rateCPT, "NSHMP-Haz Rate, "+magLabel);
				mapMaker.plot(outputDir, prefix+"_nshmp_haz", " ");
			}

			double maxAbsDiff = maxAbs(diff);
			if (maxAbsDiff > 0d) {
				CPT diffCPT = GMT_CPT_Files.DIVERGING_BAM_UNIFORM.instance().reverse().rescale(-maxAbsDiff, maxAbsDiff);
				diffCPT.setNanColor(transparent);
				mapMaker.plotXYZData(diff, diffCPT, "OpenSHA - NSHMP-Haz Rate, "+magLabel);
				mapMaker.plot(outputDir, prefix+"_diff", " ");
			}

			CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-100d, 100d);
			pDiffCPT.setNanColor(transparent);
			mapMaker.plotXYZData(pDiff, pDiffCPT, "OpenSHA vs NSHMP-Haz, % Change, "+magLabel);
			mapMaker.plot(outputDir, prefix+"_pDiff", " ");
		}
	}

	private static boolean includeForBackgroundOption(boolean gridLike, IncludeBackgroundOption bgOp) {
		switch (bgOp) {
		case INCLUDE:
			return true;
		case ONLY:
			return gridLike;
		case EXCLUDE:
			return !gridLike;
		default:
			throw new IllegalStateException("Unhandled background option: "+bgOp);
		}
	}

	private static Location toOpenSHALocation(org.opensha.nshmp.shaded.geo.NshmpLocation loc) {
		return new Location(loc.latitude, loc.longitude, loc.depth);
	}

	private static void addFiniteRupRate(GriddedGeoDataSet xyz, RuptureSurface surface, double rate) {
		GriddedRegion gridReg = xyz.getRegion();
		BitSet bits = new BitSet(gridReg.getNodeCount());
		for (Location loc : surface.getEvenlyDiscritizedListOfLocsOnSurface()) {
			int node = gridReg.indexForLocation(loc);
			if (node >= 0)
				bits.set(node);
		}
		for (int i = bits.nextSetBit(0); i >= 0; i = bits.nextSetBit(i + 1))
			xyz.add(i, rate);
	}

	private static String magLabel(double minMag) {
		return minMag == 0d ? "m0" : "m"+(float)minMag;
	}

	private static GriddedGeoDataSet diff(GriddedGeoDataSet xyz1, GriddedGeoDataSet xyz2) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(xyz1.getRegion());
		for (int i=0; i<ret.size(); i++)
			ret.set(i, xyz1.get(i) - xyz2.get(i));
		return ret;
	}

	private static GriddedGeoDataSet pDiff(GriddedGeoDataSet xyz1, GriddedGeoDataSet xyz2) {
		GriddedGeoDataSet ret = new GriddedGeoDataSet(xyz1.getRegion());
		for (int i=0; i<ret.size(); i++) {
			double val1 = xyz1.get(i);
			double val2 = xyz2.get(i);
			if (val1 == 0d && val2 == 0d)
				ret.set(i, Double.NaN);
			else if (val2 == 0d)
				ret.set(i, 100d);
			else
				ret.set(i, 100d*(val1 - val2)/val2);
		}
		return ret;
	}

	private static GriddedGeoDataSet asLog10(GriddedGeoDataSet xyz) {
		xyz = xyz.copy();
		for (int i=0; i<xyz.size(); i++)
			xyz.set(i, xyz.get(i) > 0d ? Math.log10(xyz.get(i)) : Double.NaN);
		return xyz;
	}

	private static double max(GriddedGeoDataSet xyz) {
		double max = 0d;
		for (int i=0; i<xyz.size(); i++)
			if (Double.isFinite(xyz.get(i)))
				max = Math.max(max, xyz.get(i));
		return max;
	}

	private static double maxAbs(GriddedGeoDataSet xyz) {
		double max = 0d;
		for (int i=0; i<xyz.size(); i++)
			if (Double.isFinite(xyz.get(i)))
				max = Math.max(max, Math.abs(xyz.get(i)));
		return max;
	}

	private static void writeCSV(File file, GriddedGeoDataSet fssXYZ, GriddedGeoDataSet nhXYZ,
			GriddedGeoDataSet diff, GriddedGeoDataSet pDiff) throws IOException {
		CSVFile<String> csv = new CSVFile<>(true);
		csv.addLine("Latitude", "Longitude", "OpenSHA Rate", "NSHMP-Haz Rate", "Difference", "Percent Difference");
		GriddedRegion gridReg = fssXYZ.getRegion();
		for (int i=0; i<gridReg.getNodeCount(); i++) {
			Location loc = gridReg.getLocation(i);
			csv.addLine(loc.getLatitude()+"", loc.getLongitude()+"",
					fssXYZ.get(i)+"", nhXYZ.get(i)+"", diff.get(i)+"", pDiff.get(i)+"");
		}
		csv.writeToFile(file);
	}



}
