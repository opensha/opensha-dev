package scratch.kevin.nshm27.figures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.logicTree.sampling.SamplingMethod;
import org.opensha.sha.earthquake.faultSysSolution.mpj.HPCConfig;
import org.opensha.sha.earthquake.faultSysSolution.mpj.HazardConfig;
import org.opensha.sha.earthquake.faultSysSolution.mpj.MPJ_BranchAveragedHazardScriptWriter;
import org.opensha.sha.earthquake.faultSysSolution.mpj.MPJ_BranchAveragedHazardScriptWriter.SupersamplingMode;
import org.opensha.sha.earthquake.faultSysSolution.mpj.RunConfig;
import org.opensha.sha.earthquake.faultSysSolution.util.FaultSysTools;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.logicTree.NSHM27_SeisRateModel;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_InversionScriptWriter;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_MapRegions;
import org.opensha.sha.earthquake.rupForecastImpl.nshm27.util.NSHM27_RegionLoader.NSHM27_SeismicityRegions;
import org.opensha.sha.faultSurface.utils.ptSrcCorr.PointSourceDistanceCorrections;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Writes total and tectonic-regime-specific branch-averaged hazard jobs for an
 * NSHM27 multi-regime inversion. Each calculation writes separate INCLUDE,
 * EXCLUDE, and ONLY gridded-seismicity jobs.
 */
public class NSHM27_SourceTypeHazardScriptWriter {

	private static final double GRID_SPACING_DEFAULT = 0.02;
	private static final double VS30_DEFAULT = 760d;
	private static final double[] PERIODS_DEFAULT = { 0d, 0.2d, 1d, 5d };
	private static final SupersamplingMode SUPERSAMPLING_DEFAULT = SupersamplingMode.FULL;
	private static final String SOLUTION_FILE_DEFAULT = "results_branch_averaged.zip";

	private enum SourceType {
		TOTAL(null, null),
		CRUSTAL("crustal", TectonicRegionType.ACTIVE_SHALLOW),
		INTERFACE("interface", TectonicRegionType.SUBDUCTION_INTERFACE),
		SLAB("slab", TectonicRegionType.SUBDUCTION_SLAB);

		private final String directoryToken;
		private final TectonicRegionType trt;

		private SourceType(String directoryToken, TectonicRegionType trt) {
			this.directoryToken = directoryToken;
			this.trt = trt;
		}
	}

	private static Options buildOptions() {
		Options ops = new Options();
		ops.addOption(FaultSysTools.helpOption());
		HPCConfig.addOptions(ops);
		HazardConfig.addOptions(ops);

		ops.addRequiredOption(null, "region", true, "NSHM27 region, one of: "
				+FaultSysTools.enumOptions(NSHM27_SeismicityRegions.class));
		ops.addRequiredOption(null, "input-dir-name", true,
				"Directory name containing the multi-regime branch-averaged solution.");
		ops.addOption(null, "solution-file-name", true, "Branch-averaged solution file name; default: "
				+SOLUTION_FILE_DEFAULT);
		ops.addOption(null, "output-name-add", true,
				"Additional token appended to every output directory name; may be supplied repeatedly.");
		ops.addOption(null, "supersampling-mode", true, "Supersampling mode, one of: "
				+FaultSysTools.enumOptions(SupersamplingMode.class)+"; default: "+SUPERSAMPLING_DEFAULT.name());
		return ops;
	}

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			// for Kevin's convenience within eclipse
			System.err.println("No command line arguments, using hardcoded defaults; use --help to see available options instead.");
			
			NSHM27_SeismicityRegions region = NSHM27_SeismicityRegions.AMSAM;
//			NSHM27_SeismicityRegions region = NSHM27_SeismicityRegions.GNMI;
			String inputName = "2026_11_06-nshm27-"+region.name()+"-8192samples-sobol_scrambled";
			args = new String[] {
					"--region", region.name(),
					"--input-dir-name", inputName,
					
//					"--vs30", "260",
					
					"--local-dir", "/home/kevin/OpenSHA/fss_inversions",
					
					"--hpc-site", HPCConfig.HPCSite.USC_CARC_FMPJ.name(),
					"--remote-dir", "/project2/scec_608/kmilner/fss_inversions",
					"--nodes", "36",
//					"--nodes", "34",
//					"--nodes", "18",
					"--job-time-hours", "48",
//					"--job-time-hours", "10",
					
					"--hazard-grid-spacing", region == NSHM27_SeismicityRegions.AMSAM ? "0.01" : "0.02",
					"--write-hazard-curves",
			};
		}
		
		CommandLine cmd = FaultSysTools.parseOptions(buildOptions(), args,
				NSHM27_SourceTypeHazardScriptWriter.class);

		NSHM27_SeismicityRegions seisReg = NSHM27_SeismicityRegions.valueOf(cmd.getOptionValue("region"));
		String inputDirName = cmd.getOptionValue("input-dir-name");
		String solutionFileName = cmd.getOptionValue("solution-file-name", SOLUTION_FILE_DEFAULT);
		SupersamplingMode supersamplingMode = cmd.hasOption("supersampling-mode")
				? SupersamplingMode.valueOf(cmd.getOptionValue("supersampling-mode"))
				: SUPERSAMPLING_DEFAULT;
		double gridSpacing = cmd.hasOption("hazard-grid-spacing")
				? Double.parseDouble(cmd.getOptionValue("hazard-grid-spacing")) : GRID_SPACING_DEFAULT;

		GriddedRegion hazardRegion = new GriddedRegion(NSHM27_MapRegions.valueOf(seisReg.name()).load(),
				gridSpacing, GriddedRegion.ANCHOR_0_0);
		System.out.println("Hazard region has "+hazardRegion.getNodeCount()+" nodes at "
				+(float)gridSpacing+" degree spacing");

		HPCConfig hpc = HPCConfig.builder(cmd).build();
		HazardConfig hazard = HazardConfig.builder()
				.region(hazardRegion)
				.sigmaTruncation(3d)
				.gmpe(NSHM27_InversionScriptWriter.GMM_ACTIVE)
				.gmpe(NSHM27_InversionScriptWriter.GMM_INTERFACE)
				.gmpe(NSHM27_InversionScriptWriter.GMM_SLAB)
				.vs30(VS30_DEFAULT)
				.periods(PERIODS_DEFAULT)
				.writeCurves(true)
				.forCMD(cmd)
				.build();

		List<String> outputNameAdds = new ArrayList<>();
		
		int vs30 = (int)Math.round(hazard.vs30());
		outputNameAdds.add("vs"+vs30);
		
		if (cmd.hasOption("output-name-add"))
			for (String value : cmd.getOptionValues("output-name-add"))
				outputNameAdds.add(value);

		MPJ_BranchAveragedHazardScriptWriter writer = new MPJ_BranchAveragedHazardScriptWriter();
		for (SourceType sourceType : SourceType.values()) {
			StringBuilder outputDirName = new StringBuilder(inputDirName)
					.append("-ba_hazard");
			if (sourceType.directoryToken != null)
				outputDirName.append("-").append(sourceType.directoryToken);
			for (String nameAdd : outputNameAdds)
				outputDirName.append("-").append(nameAdd);

			RunConfig run = RunConfig.builder().directoryName(outputDirName.toString()).build();
			MPJ_BranchAveragedHazardScriptWriter.Request.Builder request =
					MPJ_BranchAveragedHazardScriptWriter.Request.builder()
							.run(run)
							.hpc(hpc)
							.hazard(hazard)
							.linkFromDirectoryName(inputDirName)
							.solutionFileName(solutionFileName)
							.backgroundOptions(IncludeBackgroundOption.EXCLUDE,
									IncludeBackgroundOption.ONLY, IncludeBackgroundOption.INCLUDE)
							.supersamplingMode(supersamplingMode)
							.distanceCorrection(PointSourceDistanceCorrections.NSHM_2013);
			if (sourceType.trt != null)
				request.tectonicRegionTypes(sourceType.trt);
			writer.writeScripts(request.build());
		}
	}
}
