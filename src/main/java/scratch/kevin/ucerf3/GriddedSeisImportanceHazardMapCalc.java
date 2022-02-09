package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Region;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.hazardMap.components.CalculationInputsXMLFile;
import org.opensha.sha.calc.hazardMap.components.CalculationSettings;
import org.opensha.sha.calc.hazardMap.components.CurveResultsArchiver;
import org.opensha.sha.calc.hazardMap.mpj.MPJHazardCurveDriver;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.earthquake.param.IncludeBackgroundOption;
import org.opensha.sha.earthquake.param.IncludeBackgroundParam;
import org.opensha.sha.earthquake.param.ProbabilityModelOptions;
import org.opensha.sha.earthquake.param.ProbabilityModelParam;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.util.TectonicRegionType;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import scratch.UCERF3.U3FaultSystemSolution;
import scratch.UCERF3.erf.FaultSystemSolutionERF;
import scratch.UCERF3.griddedSeismicity.GridSourceFileReader;
import scratch.UCERF3.utils.U3FaultSystemIO;

public class GriddedSeisImportanceHazardMapCalc {
	
	static enum CalcType {
		FULL("Full"),
		SUPRA_ONLY("Fault Supra"),
		SUPRA_PLUS_SUB_ONLY("Fault Supra+Sub"),
		GRIDDED_ONLY("Gridded (Off/Sub)"),
		OFF_FAULT_ONLY("Off Only");
		
		private String label;
		private CalcType(String label) {
			this.label = label;
		}
		
		public String getLabel() {
			return label;
		}
	}

	public static void main(String[] args) throws IOException, DocumentException {
		File remoteBaseDir = new File("/home/scec-02/kmilner/ucerf3/maps");
		File localBaseDir = new File("/home/kevin/OpenSHA/UCERF3/maps");
		
//		String dirName = "2017_07_14-ucerf3-gridded-tests";
//		String dirName = "2017_08_07-ucerf3-full-ba-gridded-tests";
		String dirName = "2018_08_30-ucerf3-geol-gridded-tests";
		
		File localU3File = new File("/home/kevin/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/"
				+ "FM3_1_GEOL_MEAN_BRANCH_AVG_SOL.zip");
//				+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_MEAN_BRANCH_AVG_SOL.zip");
		U3FaultSystemSolution u3Sol = U3FaultSystemIO.loadSol(localU3File);
		
		Region region = new CaliforniaRegions.RELM_TESTING();
		double spacing = 0.02;
//		String imt = PGA_Param.NAME;
//		double period = 0d;
		String imt = SA_Param.NAME;
		double period = 0.5d;
		ScalarIMR imr = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		double duration = 1d;
		
		dirName += "-"+imt.toLowerCase();
		if (period > 0)
			dirName += "-"+(float)period+"s";
		
//		CalcType[] types = CalcType.values();
		CalcType[] types = { CalcType.FULL, CalcType.SUPRA_PLUS_SUB_ONLY, CalcType.SUPRA_ONLY };
		
		File localDir = new File(localBaseDir, dirName);
		File remoteDir = new File(remoteBaseDir, dirName);
		System.out.println(localDir.getAbsolutePath()+" => "+remoteDir.getAbsolutePath());
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
		int mins = 24*60;
		int nodes = 8;
		int ppn = 20;
		String queue = "scec";
		
		// build solution with only sub seismogenic, no unsassociated
		GridSourceProvider origProv = u3Sol.getGridSourceProvider();
		GriddedRegion seisReg = origProv.getGriddedRegion();
		Map<Integer, IncrementalMagFreqDist> nodeSubSeisMFDs = new HashMap<>();
		Map<Integer, IncrementalMagFreqDist> nodeUnassociatedMFDs = new HashMap<>();
		Map<Integer, IncrementalMagFreqDist> emptyMFDs = new HashMap<>();
		for (int i=0; i<seisReg.getNodeCount(); i++) {
			if (origProv.getMFD_SubSeisOnFault(i) != null)
				nodeSubSeisMFDs.put(i, origProv.getMFD_SubSeisOnFault(i));
			if (origProv.getMFD_Unassociated(i) != null)
				nodeUnassociatedMFDs.put(i, origProv.getMFD_Unassociated(i));
		}
		GridSourceFileReader subSeisOnlyProv = new GridSourceFileReader(seisReg, nodeSubSeisMFDs, emptyMFDs);
		U3FaultSystemSolution u3SubSeisOnlySol = new U3FaultSystemSolution(u3Sol.getRupSet(), u3Sol.getRateForAllRups());
		u3SubSeisOnlySol.setGridSourceProvider(subSeisOnlyProv);
		GridSourceFileReader unassociatedOnlyProv = new GridSourceFileReader(seisReg, emptyMFDs, nodeUnassociatedMFDs);
		U3FaultSystemSolution u3UnassociatedOnlySol = new U3FaultSystemSolution(u3Sol.getRupSet(), u3Sol.getRateForAllRups());
		u3UnassociatedOnlySol.setGridSourceProvider(unassociatedOnlyProv);
		
		// write solutions
		String fullSolName = "full_solution.zip";
		File remoteFullSolFile = new File(remoteDir, fullSolName);
		U3FaultSystemIO.writeSol(u3Sol, new File(localDir, fullSolName));
		
		String subSeisOnlySolName = "gridded_sub_seis_only_solution.zip";
		File remoteSubSeisOnlySolFile = new File(remoteDir, subSeisOnlySolName);
		U3FaultSystemIO.writeSol(u3SubSeisOnlySol, new File(localDir, subSeisOnlySolName));
		
		String unassociatedOnlySolName = "gridded_unassociated_only_solution.zip";
		File remoteUnassociatedOnlySolFile = new File(remoteDir, unassociatedOnlySolName);
		U3FaultSystemIO.writeSol(u3UnassociatedOnlySol, new File(localDir, unassociatedOnlySolName));
		
		// set up inputs
		GriddedRegion gridded = new GriddedRegion(region, spacing, null);
		
		imr.setParamDefaults();
		imr.setIntensityMeasure(imt);
		if (period > 0)
			SA_Param.setPeriodInSA_Param(imr.getIntensityMeasure(), period);
		
		List<Site> sites = Lists.newArrayList();
		for (int i=0; i<gridded.getNodeCount(); i++) {
			Site site = new Site(gridded.locationForIndex(i));
			site.addParameterList(imr.getSiteParams());
			sites.add(site);
		}
		
		ArbitrarilyDiscretizedFunc xValues = new IMT_Info().getDefaultHazardCurve(imt);
		double maxSourceDistance = 200;
		
		Map<String, DiscretizedFunc> xValsMap = Maps.newHashMap();
		xValsMap.put("curves", xValues);
		CalculationSettings calcSettings = new CalculationSettings(xValues, maxSourceDistance);
		
		File javaBin = USC_HPCC_ScriptWriter.JAVA_BIN;
		File jarFile = new File(remoteBaseDir, "opensha-dev-all.jar");
		
		List<File> classpath = Lists.newArrayList();
		classpath.add(jarFile);
		
		MPJExpressShellScriptWriter mpj = new MPJExpressShellScriptWriter(javaBin, 60000, classpath,
				USC_HPCC_ScriptWriter.MPJ_HOME);
		mpj.setUseLaunchWrapper(true);
		
		List<Map<TectonicRegionType, ScalarIMR>> imrMaps = Lists.newArrayList();
		
		HashMap<TectonicRegionType, ScalarIMR> map = Maps.newHashMap();
		map.put(TectonicRegionType.ACTIVE_SHALLOW, imr);
		imrMaps.add(map);
		
		// write jobs
		for (CalcType type : types) {
			String subDirName = type.name().toLowerCase();
			System.out.println("Processing "+type);
			File localSubDir = new File(localDir, subDirName);
			Preconditions.checkState(localSubDir.exists() || localSubDir.mkdir());
			File remoteSubDir = new File(remoteDir, subDirName);
			
			FaultSystemSolutionERF erf = new FaultSystemSolutionERF();
			
			switch (type) {
			case FULL:
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteFullSolFile);
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
				break;
			case SUPRA_ONLY:
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteFullSolFile);
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.EXCLUDE);
				break;
			case SUPRA_PLUS_SUB_ONLY:
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteSubSeisOnlySolFile);
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.INCLUDE);
				break;
			case GRIDDED_ONLY:
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteFullSolFile);
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
				break;
			case OFF_FAULT_ONLY:
				erf.setParameter(FaultSystemSolutionERF.FILE_PARAM_NAME, remoteUnassociatedOnlySolFile);
				erf.setParameter(IncludeBackgroundParam.NAME, IncludeBackgroundOption.ONLY);
				break;
			default:
				throw new IllegalStateException("Unknown type: "+type);
			}
			
			erf.setParameter(ProbabilityModelParam.NAME, ProbabilityModelOptions.POISSON);
			erf.getTimeSpan().setDuration(duration);
			
			File curveDir = new File(remoteSubDir, "curves");
			CurveResultsArchiver archiver = new BinaryCurveArchiver(curveDir, sites.size(), xValsMap);
			
			CalculationInputsXMLFile inputs = new CalculationInputsXMLFile(erf, imrMaps, sites, calcSettings, archiver);
			
			File localInputsFile = new File(localSubDir, "inputs.xml");
			File remoteInputsFile = new File(remoteSubDir, "inputs.xml");
			XMLUtils.writeObjectToXMLAsRoot(inputs, localInputsFile);
			
//			String cliArgs = "--max-dispatch 1000 --mult-erfs "+inputsFile.getAbsolutePath();
			String cliArgs = "--max-dispatch 1000 "+remoteInputsFile.getAbsolutePath();
			
			List<String> script = mpj.buildScript(MPJHazardCurveDriver.class.getName(), cliArgs);
			USC_HPCC_ScriptWriter writer = new USC_HPCC_ScriptWriter();
			
			script = writer.buildScript(script, mins, nodes, ppn, queue);
			
			File pbsFile = new File(localSubDir, subDirName+".slurm");
			JavaShellScriptWriter.writeScript(pbsFile, script);
		}
	}

}
