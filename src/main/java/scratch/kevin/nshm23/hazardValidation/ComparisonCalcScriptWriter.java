package scratch.kevin.nshm23.hazardValidation;

import java.io.File;
import java.io.IOException;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

import org.opensha.commons.data.siteData.CONUS_Versions;
import org.opensha.commons.data.siteData.SiteDataValueListList;
import org.opensha.commons.data.siteData.impl.CONUS_SiteDataProvider;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.earthquake.faultSysSolution.util.SolHazardMapCalc;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.util.NSHM23_RegionLoader;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.util.TectonicRegionType;

import gov.usgs.earthquake.nshmp.erf.mpj.HPCConfig;
import gov.usgs.earthquake.nshmp.erf.mpj.HPCConfig.HPCSite;
import gov.usgs.earthquake.nshmp.erf.mpj.HazardConfig;
import gov.usgs.earthquake.nshmp.erf.mpj.MPJ_BranchAveragedHazardScriptWriter;
import gov.usgs.earthquake.nshmp.erf.mpj.MPJ_BranchAveragedHazardScriptWriter.SupersamplingMode;
import gov.usgs.earthquake.nshmp.erf.mpj.RunConfig;

public class ComparisonCalcScriptWriter {

	public static void main(String[] args) throws IOException {
		// start with WUS-only
		String regToken = "WUS";
		Region reg = NSHM23_RegionLoader.loadFullConterminousWUS();
		double spacing = 0.1;
		String linkFromDir = "2024_02_02-nshm23_branches-WUS_FM_v3";
		String solFileName = "results_WUS_FM_v3_branch_averaged_gridded_simplified_revised2026.zip";

		File localMainDir = new File("/home/kevin/OpenSHA/fss_inversions");
		
		GriddedRegion gridReg = new GriddedRegion(reg, spacing, GriddedRegion.ANCHOR_0_0);
		// add site data
		SiteDataValueListList siteData = new SiteDataValueListList();
		CONUS_SiteDataProvider data10 = new CONUS_SiteDataProvider(
				CONUS_SiteDataProvider.TYPE_DEPTH_TO_1_0, CONUS_Versions.NSHM23);
		CONUS_SiteDataProvider data25 = new CONUS_SiteDataProvider(
				CONUS_SiteDataProvider.TYPE_DEPTH_TO_2_5, CONUS_Versions.NSHM23);
		siteData.add(data10.getAnnotatedValues(gridReg.getNodeList()));
		siteData.add(data25.getAnnotatedValues(gridReg.getNodeList()));
		gridReg.setSiteData(siteData);
		
		// TODO add stable continental when ready
		String gmpeToken = "active_only";
		Map<TectonicRegionType, AttenRelRef> gmpes = new EnumMap<>(TectonicRegionType.class);
		gmpes.put(TectonicRegionType.ACTIVE_SHALLOW, AttenRelRef.USGS_NSHM23_ACTIVE);
		
		SolHazardMapCalc.loadSites(gridReg, gmpes);
		
		HPCSite hpcSite = HPCSite.USC_CARC_FMPJ;
		File remoteMainDir = new File("/project2/scec_608/kmilner/fss_inversions");
		
		HPCConfig hpc = HPCConfig.builder(hpcSite)
				.localMainDir(localMainDir)
				.remoteMainDir(remoteMainDir)
				.build();
		
		HazardConfig hazard = HazardConfig.builder()
				.region(gridReg)
				.sigmaTruncation(3d)
				.gmpes(gmpes.values())
				.vs30(760d)
				.build();


		RunConfig run = RunConfig.builder()
				.baseName("nshm23")
				.addNameToken("hazard_validation")
				.addNameToken(regToken)
				.addNameToken(gmpeToken)
				.build();

		MPJ_BranchAveragedHazardScriptWriter.Request request = MPJ_BranchAveragedHazardScriptWriter.Request.builder()
				.distanceCorrection(PointSourceDistanceCorrections.NSHM_2013)
				.backgroundOptions(IncludeBackgroundOption.values())
				.run(run)
				.linkFromDirectoryName(linkFromDir)
				.solutionFileName(solFileName)
				.noMFDs(true)
				.supersamplingMode(SupersamplingMode.FULL)
				.hazard(hazard)
				.hpc(hpc)
				.build();
		
		new MPJ_BranchAveragedHazardScriptWriter().writeScripts(request);
	}

}
