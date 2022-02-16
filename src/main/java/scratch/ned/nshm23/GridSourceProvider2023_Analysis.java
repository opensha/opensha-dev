package scratch.ned.nshm23;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.jfree.data.Range;
import org.opensha.commons.data.function.HistogramFunction;
import org.opensha.commons.data.function.XY_DataSet;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSymbol;
import org.opensha.commons.mapping.gmt.GMT_MapGenerator;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.mapping.gmt.gui.GMT_MapGuiBean;
import org.opensha.commons.mapping.gmt.gui.ImageViewerWindow;
import org.opensha.commons.param.impl.CPTParameter;
import org.opensha.commons.util.ExceptionUtils;
import org.opensha.commons.util.FileUtils;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;
import org.opensha.sha.earthquake.faultSysSolution.FaultSystemSolution;
import org.opensha.sha.earthquake.faultSysSolution.modules.GridSourceProvider;
import org.opensha.sha.gui.infoTools.CalcProgressBar;
import org.opensha.sha.magdist.GutenbergRichterMagFreqDist;
import org.opensha.sha.magdist.IncrementalMagFreqDist;
import org.opensha.sha.magdist.SummedMagFreqDist;

import scratch.UCERF3.analysis.GMT_CA_Maps;
import scratch.UCERF3.enumTreeBranches.SpatialSeisPDF;
import scratch.UCERF3.erf.ETAS.SeisDepthDistribution;
import scratch.UCERF3.griddedSeismicity.GridReader;
import scratch.UCERF3.utils.RELM_RegionUtils;
import scratch.UCERF3.utils.U3FaultSystemIO;
import scratch.ned.FSS_Inversion2019.PlottingUtils;

/**
 * This makes various plots and other analyses of GridSourceProvider2023.
 * 
 * To do:
 * 
 * 1) generalize to utilize any region (currently hard-coded for RELM: "mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator()")
 * 
 * 2) improve specification of where result is stored (currently put in project root)
 * 
 * 3) are direct references to class variables (e.g., "gridProvider.cgr") OK?
 * 
 * 
 * @author field
 *
 */
public class GridSourceProvider2023_Analysis {


	/**
	 * This plots the supra seis rates above the specified magnitude for cubes at the given depth
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public static String plotSupraSeisRateAtDepthMap(GridSourceProvider2023 gridProvider, double depth, String dirName) {
		CubedGriddedRegion cgr = gridProvider.cgr;
		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);
		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			SummedMagFreqDist mfd = gridProvider.getSupraSeisMFD_ForCube(cubeIndex);
			double rate = 0.0;
			if(mfd != null)
				rate = mfd.getTotalIncrRate();
			if(rate == 0.0)
				rate = 1e-16;
			xyzDataSet.set(i, rate);
//			Location loc = xyzDataSet.getLocation(i);
//			System.out.println(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+xyzDataSet.get(i));
		}
		progressBar.showProgress(false);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.MAX_SPECTRUM.getFileName());
		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);

		
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		double maxZ = Math.ceil(Math.log10(xyzDataSet.getMaxZ()))+0.5;
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,maxZ-5);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,maxZ);
		
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-9d);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-4d);

		String metadata = "Map from calling plotSuparSeisRateAtDepthMap(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Supra Seis Rates at "+depth+" km depth", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "Supra Seis Rates at depth above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}
	
	
	
	
	/**
	 * This plots the cube rates above the specified magnitude for cubes at the given depth
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public static String plotCubeRateAboveMagAtDepthMap(GridSourceProvider2023 gridProvider, double mag, double depth, String dirName) {
		CubedGriddedRegion cgr = gridProvider.cgr;
		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);
		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			SummedMagFreqDist mfd = gridProvider.getTotalMFD_ForCube(cubeIndex);
			double rate = 0.0;
			if(mfd != null)
				rate = mfd.getCumRate(mfd.getClosestXIndex(mag));
			if(rate == 0.0)
				rate = 1e-16;
			xyzDataSet.set(i, rate);
//			Location loc = xyzDataSet.getLocation(i);
//			System.out.println(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+xyzDataSet.get(i));
		}
		progressBar.showProgress(false);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.MAX_SPECTRUM.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);

		
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		double maxZ = Math.ceil(Math.log10(xyzDataSet.getMaxZ()))+0.5;
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,maxZ-5);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,maxZ);
		
		if(mag< 5.5) {
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-5d);
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-1d);			
		}
		else {
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-9d);
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-4d);						
		}

		String metadata = "Map from calling plotSuparSeisRateAtDepthMap(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Rates for M>="+mag+" at "+depth+" km depth", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "Cum Rates at depth above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}
	
	/**
	 * This plots the gridded seismicity rates above the specified magnitude for each grid cell 
	 * (only gridded seismicy included, no supraseis on-fault ruptures)
	 * @param dirName
	 * @return
	 */
	public static String plotGridRateAboveMag(GridSourceProvider gridProvider, double mag, String dirName) {
		GriddedRegion gridReg = gridProvider.getGriddedRegion();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridReg, true);
		int numGridPoints = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		for(int i=0; i<numGridPoints;i++) {
			progressBar.updateProgress(i, numGridPoints);
			IncrementalMagFreqDist mfd = gridProvider.getMFD(i);
			double rate = 0.0;
			if(mfd != null)
				rate = mfd.getCumRate(mfd.getClosestXIndex(mag));
			if(rate == 0.0)
				rate = 1e-16;
			xyzDataSet.set(i, rate);
//			Location loc = xyzDataSet.getLocation(i);
//			System.out.println(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+xyzDataSet.get(i));
		}
		progressBar.showProgress(false);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.MAX_SPECTRUM.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridReg.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridReg.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridReg.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridReg.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridReg.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);

		
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		double maxZ = Math.ceil(Math.log10(xyzDataSet.getMaxZ()))+0.5;
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,maxZ-5);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,maxZ);
		
		if(mag< 5.5) {
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-5d);
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-1d);			
		}
		else {
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-9d);
			mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,-4d);						
		}

		String metadata = "Map from calling plotGridRateAboveMagAtDepthMap(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Grid Rates for M>="+mag, metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "Cum Rates at depth above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}


	/**
	 * This plots the ratio of gridded seismicity rates above the specified magnitude for each grid cell 
	 * (only gridded seismicy included, no supraseis on-fault ruptures).
	 * @param gridProvider1
	 * @param gridProvider2
	 * @param mag
	 * @param dirName
	 * @return
	 */
	public static String plotGridRateRatioAboveMag(GridSourceProvider gridProvider1, GridSourceProvider gridProvider2, double mag, String dirName) {
		GriddedRegion gridReg = gridProvider1.getGriddedRegion();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridReg, true);
		int numGridPoints = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		double min=Double.MAX_VALUE, max = -100;
		
		for(int i=0; i<numGridPoints;i++) {
			progressBar.updateProgress(i, numGridPoints);
			IncrementalMagFreqDist mfd1 = gridProvider1.getMFD(i);
			IncrementalMagFreqDist mfd2 = gridProvider2.getMFD(i);
			double ratio = 1.0;
			if(mfd1 != null && mfd2 != null) {
				ratio = mfd1.getCumRate(mfd1.getClosestXIndex(mag))/mfd2.getCumRate(mfd2.getClosestXIndex(mag));
//				ratio = mfd1.getY(mfd1.getClosestXIndex(mag))/mfd2.getY(mfd2.getClosestXIndex(mag));
			}
			xyzDataSet.set(i, ratio);
			if(max<ratio) max=ratio;
			if(min>ratio) min=ratio;
//			Location loc = xyzDataSet.getLocation(i);
//			System.out.println(loc.getLongitude()+"\t"+loc.getLatitude()+"\t"+xyzDataSet.get(i));
		}
		progressBar.showProgress(false);
		System.out.println("min ratio = "+min+"\nmax ratio = "+max);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.UCERF3_RATIOS.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridReg.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridReg.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridReg.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridReg.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridReg.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-2d);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,2d);

		String metadata = "Map from calling plotGridRateAboveMagAtDepthMap(*) method";
		
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Grid Rate Ratio for M>="+mag, metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "Grid rate ratio above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}

	
	public static String plotRatioMap(GriddedGeoDataSet xyzDataSet, double deltaLatLon, String dirName) {

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.UCERF3_RATIOS.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,xyzDataSet.getMinLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,xyzDataSet.getMaxLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,xyzDataSet.getMinLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,xyzDataSet.getMaxLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, deltaLatLon);	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-2d);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,2d);

		String metadata = "Map from calling plotRatioMap(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Ratio", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return "Grid ratio map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}

	
	
	/**
	 * This plots the event rates above the specified magnitude for cubes at the given depth
	 * @param depth
	 * @param dirName
	 * @return
	 */
	public static void plotCubeFractionOnAndOffFaultAtDepth(GridSourceProvider2023 gridProvider, double depth) {

		CubedGriddedRegion cgr = gridProvider.cgr;

		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		String dirNameOn = "CubeFractionOnFaultAtDepth"+depth+"km";
		String dirNameOff = "CubeFractionOffFaultAtDepth"+depth+"km";
		GriddedGeoDataSet xyzDataSetFractOn = new GriddedGeoDataSet(gridRegForCubes, true);
		GriddedGeoDataSet xyzDataSetFractOff = new GriddedGeoDataSet(gridRegForCubes, true);

		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSetFractOn.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		double max = -1, min = 2;;
		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			
			HashMap<Integer,Double> sectWtMap = gridProvider.sectDistWtMapAtCubeList.get(cubeIndex);
			double wtSum =0;
			for(int s:sectWtMap.keySet()) {
				wtSum+=sectWtMap.get(s);
			}
			xyzDataSetFractOn.set(i, wtSum);
			xyzDataSetFractOff.set(i, 1.0-wtSum);
			
			if(max<wtSum) max=wtSum;
			if(min>wtSum) min=wtSum;
		}
		
		System.out.println("plotCubeFractionOnAndOffFaultAtDepth:"+"\n\tmin="+min+"\n\tmax="+max);
		
		progressBar.showProgress(false);

		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.GMT_POLAR.getFileName());
		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);

		
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,false);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_FROMDATA);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
//		double maxZ = Math.ceil(Math.log10(xyzDataSet.getMaxZ()))+0.5;
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,maxZ-5);
//		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,maxZ);
		
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-1.0);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,1.0);			

		String metadata = "Map from calling plotCubeFractionOnAndOffFaultAtDepth(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSetFractOn, "Fraction On Fault at "+depth+" km depth", metadata, dirNameOn);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirNameOn);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		try {
			String url = mapGen.makeMapUsingServlet(xyzDataSetFractOff, "Fraction Off Fault at "+depth+" km depth", metadata, dirNameOff);
			metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
			ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
			
			File downloadDir = new File(dirNameOff);
			if (!downloadDir.exists())
				downloadDir.mkdir();
			File zipFile = new File(downloadDir, "allFiles.zip");
			// construct zip URL
			String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
			FileUtils.downloadURL(zipURL, zipFile);
			FileUtils.unzipFile(zipFile, downloadDir);

//		System.out.println("GMT Plot Filename: "+name);
	} catch (Exception e) {
		e.printStackTrace();
	}

//		return "For rates at depth above mag map: "+mapGen.getGMTFilesWebAddress()+" (deleted at midnight)";
	}
	
	
	
	
	

	/**
	 * This plot the charFactor at the specified depth and you must provide the total region target b-value
	 * for extrapolation purposes.  This is close, but not exactly the same as used in U3ETAS.
	 * @param gridProvider
	 * @param depth
	 * @param bVal
	 */
	public static void plotCharFactorAtDepth(GridSourceProvider2023 gridProvider, double depth, double bVal) {
		
		CubedGriddedRegion cgr = gridProvider.cgr;

		GriddedRegion gridRegForCubes = cgr.getGridRegForCubes();
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridRegForCubes, true);

		int depthIndex = cgr.getCubeDepthIndex(depth);
		int numCubesAtDepth = xyzDataSet.size();
		CalcProgressBar progressBar = new CalcProgressBar("Looping over all points", "junk");
		progressBar.showProgress(true);
		
		for(int i=0; i<numCubesAtDepth;i++) {
			progressBar.updateProgress(i, numCubesAtDepth);
			int cubeIndex = cgr.getCubeIndexForRegAndDepIndices(i, depthIndex);
			
			HashMap<Integer,Double> sectWtMap = gridProvider.sectDistWtMapAtCubeList.get(cubeIndex);
			
			double aveMinSupraMag=0;
			if(sectWtMap.size()==0) { // no sections nucleate here
				aveMinSupraMag=6.35;
			}
			else {
				double totWt=0;
				for(int s:sectWtMap.keySet()) {
					IncrementalMagFreqDist mfd = gridProvider.longTermSupraSeisMFD_OnSectArray[s];
					double minMag = mfd.getMinMagWithNonZeroRate();
					double wt = mfd.getTotalIncrRate()*sectWtMap.get(s)/gridProvider.totDistWtsAtCubesForSectArray[s];
					aveMinSupraMag += wt*minMag;
					totWt+=wt;
				}
				aveMinSupraMag /= totWt;			
			}
			
			SummedMagFreqDist totalMFD = gridProvider.getTotalMFD_ForCube(cubeIndex);
			
			int index = totalMFD.getClosestXIndex(aveMinSupraMag);
			aveMinSupraMag = totalMFD.getX(index);
			
			double minMag = totalMFD.getMinMagWithNonZeroRate();
			double maxMag = totalMFD.getMaxMagWithNonZeroRate();
			
			GutenbergRichterMagFreqDist gr = new GutenbergRichterMagFreqDist(totalMFD.getMinX(), totalMFD.size(), totalMFD.getDelta(),
					minMag, maxMag, 1.0, bVal);
			gr.scaleToIncrRate(minMag, totalMFD.getY(minMag));
			
			double charFact = totalMFD.getCumRate(aveMinSupraMag)/gr.getCumRate(aveMinSupraMag);
			xyzDataSet.set(i, charFact);
			
		}
		
		progressBar.showProgress(false);
		
		String dirName = "CharFactorAtDepth"+depth+"km";
		
		GMT_MapGenerator mapGen = GMT_CA_Maps.getDefaultGMT_MapGenerator();
		
		CPTParameter cptParam = (CPTParameter )mapGen.getAdjustableParamsList().getParameter(GMT_MapGenerator.CPT_PARAM_NAME);
		cptParam.setValue(GMT_CPT_Files.UCERF3_RATIOS.getFileName());
//		cptParam.getValue().setBelowMinColor(Color.WHITE);
		
		mapGen.setParameter(GMT_MapGenerator.MIN_LAT_PARAM_NAME,gridRegForCubes.getMinGridLat());
		mapGen.setParameter(GMT_MapGenerator.MAX_LAT_PARAM_NAME,gridRegForCubes.getMaxGridLat());
		mapGen.setParameter(GMT_MapGenerator.MIN_LON_PARAM_NAME,gridRegForCubes.getMinGridLon());
		mapGen.setParameter(GMT_MapGenerator.MAX_LON_PARAM_NAME,gridRegForCubes.getMaxGridLon());
		mapGen.setParameter(GMT_MapGenerator.GRID_SPACING_PARAM_NAME, gridRegForCubes.getLatSpacing());	// assume lat and lon spacing are same
		mapGen.setParameter(GMT_MapGenerator.GMT_SMOOTHING_PARAM_NAME, false);
		mapGen.setParameter(GMT_MapGenerator.GRD_VIEW_PARAM_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.LOG_PLOT_NAME,true);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MODE_NAME,GMT_MapGenerator.COLOR_SCALE_MODE_MANUALLY);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MIN_PARAM_NAME,-3d);
		mapGen.setParameter(GMT_MapGenerator.COLOR_SCALE_MAX_PARAM_NAME,3d);
		
		
		String metadata = "Map from calling plotCharFactorAtDepth(*) method";
		
		try {
				String url = mapGen.makeMapUsingServlet(xyzDataSet, "Bulge at "+depth+" km depth", metadata, dirName);
				metadata += GMT_MapGuiBean.getClickHereHTML(mapGen.getGMTFilesWebAddress());
				ImageViewerWindow imgView = new ImageViewerWindow(url,metadata, true);		
				
				File downloadDir = new File(dirName);
				if (!downloadDir.exists())
					downloadDir.mkdir();
				File zipFile = new File(downloadDir, "allFiles.zip");
				// construct zip URL
				String zipURL = url.substring(0, url.lastIndexOf('/')+1)+"allFiles.zip";
				FileUtils.downloadURL(zipURL, zipFile);
				FileUtils.unzipFile(zipFile, downloadDir);

//			System.out.println("GMT Plot Filename: "+name);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * This plots the classic MFDs for the paper
	 */
	public static void plotTotalMFDs(GridSourceProvider2023 gridProvider) {
		
		gridProvider.totGriddedSeisMFD.setName("totGriddedSeisMFD");
		
		SummedMagFreqDist testMFD = gridProvider.initSummedMFD();
		testMFD.setName("Test totGriddedSeisMFD");
		testMFD.addIncrementalMagFreqDist(gridProvider.totalSubSeisOnFaultMFD);
		testMFD.addIncrementalMagFreqDist(gridProvider.totalTrulyOffFaultMFD);
		
		SummedMagFreqDist totalMFD = gridProvider.initSummedMFD();
		totalMFD.setName("totalMFD");
		totalMFD.addIncrementalMagFreqDist(gridProvider.totGriddedSeisMFD);
		totalMFD.addIncrementalMagFreqDist(gridProvider.totalSupraSeisOnFaultMFD);
		
		ArrayList<XY_DataSet> funcs = new  ArrayList<XY_DataSet>();
		funcs.add(gridProvider.totGriddedSeisMFD);
		funcs.add(gridProvider.totalSubSeisOnFaultMFD);
		funcs.add(gridProvider.totalTrulyOffFaultMFD);
		funcs.add(testMFD);
		funcs.add(gridProvider.totalSupraSeisOnFaultMFD);
		funcs.add(totalMFD);


		
		ArrayList<PlotCurveCharacterstics> plotChars = new ArrayList<PlotCurveCharacterstics>();
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GRAY));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLUE));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.RED));
		plotChars.add(new PlotCurveCharacterstics(PlotSymbol.CROSS, 4f, Color.BLACK));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.GREEN));
		plotChars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 2f, Color.BLACK));

		
		Range xRange = new Range(5.0,9.0);
		Range yRange = new Range(1e-7, 1);
		
		PlottingUtils.writeAndOrPlotFuncs(funcs, plotChars, null, "Magnitude", "Rate (per yr)", xRange, yRange, 
				false, true, 3.5, 3.0, null, true);	

	}


	/**
	 * 
	 * To plot this result in 3D (rotatable) in Igor64 paste the data in and then paste 
	 * this into the command window (all together, but without the leading "	 * "):
	 * 
	 * Concatenate {lon,lat,depth}, tripletWave
	 * NewGizmo
	 * AppendToGizmo DefaultScatter= root:tripletWave
	 * ModifyGizmo makeTripletColorWave = {wt, grays, 1 }
	 * ModifyGizmo stopUpdates
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ scatterColorType,1}
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ colorWave,root:wt_C}
	 * ModifyGizmo resumeUpdates
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,2}
	 * ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ size,0.5}
	 * 
	 * @param parSectID
	 * @param fileName
	 */
	public static void writeCubeLocsAndWtsForForParentSect(GridSourceProvider2023 gridProvider, int parSectID, String fileName) {
		
		
//		for(int s=0;s<rupSet.getNumSections();s++) {
//			String name = rupSet.getFaultSectionData(s).getName();
//			if(name.contains("Mojave")) {
//				System.out.println(rupSet.getFaultSectionData(s).getParentSectionId()+"\t"+rupSet.getFaultSectionData(s).getParentSectionName());
//			}
//		}
//		System.exit(-1);
		
		FaultSystemRupSet rupSet = gridProvider.rupSet;
				
		System.out.println("writeLocAndFactionOnForCubes:\n");
		boolean didIt = false;
		
		FileWriter fw;
		try {
			fw = new FileWriter(fileName);
			fw.write("lon\tlat\tdepth\twt\n");
			for(int s=0;s<rupSet.getNumSections();s++) {
				if(rupSet.getFaultSectionData(s).getParentSectionId() == parSectID) {
					if(!didIt) {
						System.out.println(rupSet.getFaultSectionData(s).getParentSectionName());
						didIt=true;
					}
					
					HashMap<Integer,Double> cubeDistMap = new HashMap<Integer,Double>();
					gridProvider.getCubeDistancesForFaultSection(s, cubeDistMap);
					if(cubeDistMap != null) {	// null if section is outside the region
						for(int cubeIndex:cubeDistMap.keySet()) {
							HashMap<Integer,Double> sectWtMap = gridProvider.sectDistWtMapAtCubeList.get(cubeIndex);
							double wt = sectWtMap.get(s);
							Location loc = gridProvider.cgr.getCubeLocationForIndex(cubeIndex);
							fw.write((float)loc.getLongitude()+"\t"+(float)loc.getLatitude()+
									"\t"+(float)-loc.getDepth()+"\t"+(float)wt+"\n");
						}			
					}
				}
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * The computes how many different sections nucleate in each cube and then creates a
	 * histogram (how many have 0, 1, 2, etc sections in the cube)
	 * 
	 * @param gridProvider
	 * @param writeMax - this will write out the sections inside cubes with the man number of sections
	 * @return
	 */
	public static HistogramFunction computeHistogramOfNumSectionsInCubes(GridSourceProvider2023 gridProvider, boolean writeMax) {
		HistogramFunction numCubesWithNumSectHist = new HistogramFunction(0.0, 21,1.0);
		
		CubedGriddedRegion cgr = gridProvider.cgr;
		List<int[]> sectAtCubeList = gridProvider.sectAtCubeList;


		int maxNum = 0;
		for(int c=0; c<cgr.getNumCubes(); c++) {
			int num = sectAtCubeList.get(c).length;
			numCubesWithNumSectHist.add(num, 1.0);
			if(maxNum<num)
				maxNum = num;
		}
		
		// write out those with the max muber:
		if(writeMax) {
			for(int c=0; c<cgr.getNumCubes(); c++) {
				int num = sectAtCubeList.get(c).length;
				if(num==maxNum) {
					System.out.println("\nCube "+c+ " has "+ maxNum+" sections; "+cgr.getCubeLocationForIndex(c));
					for(int i=0;i<sectAtCubeList.get(c).length;i++) {
						int s = sectAtCubeList.get(c)[i];
						float dist = gridProvider.sectDistToCubeList.get(c)[i];
						float wt = (float)gridProvider.getDistWt(dist);
						System.out.println("\t"+s+"\t"+dist+"\t"+wt+"\t"+gridProvider.rupSet.getFaultSectionData(s).getName());
					}
				}
			}
			System.out.println(numCubesWithNumSectHist);	
		}
		return numCubesWithNumSectHist;
	}
	
	
	/**
	 * The following shows that spatialPDF is not the same as the PDF of M4 rates among grid cells.
	 * This is because sub-seis on-fault MFDs were re-distributed among the grid cells inside its polygon,
	 * weighted by the fractional area inside the polygon but not by the original grid wt (the latter is
	 * the critical distinction).  In other words, subseis rates are constant inside the polygon.
	 * @param spatialPDF
	 * @param gridSrcProviderU3
	 */
	public  static void compareU3_InitialAndFinalSpatialPDF(double[] spatialPDF, GridSourceProvider gridSrcProviderU3, String folderName) {
		double[] spatialPDF_Test = new double[spatialPDF.length];
		double sum=0;
		for(int i=0;i<gridSrcProviderU3.size();i++) {
			double rate = gridSrcProviderU3.getMFD(i).getCumRate(4.05);
			spatialPDF_Test[i] =  rate;
			sum += rate;
		}
		double min=Double.MAX_VALUE, max= -100;
		GriddedGeoDataSet xyzDataSet = new GriddedGeoDataSet(gridSrcProviderU3.getGriddedRegion(), true);
		for(int i=0;i<spatialPDF.length;i++) {
			spatialPDF_Test[i] /= sum;
			double ratio = spatialPDF[i]/spatialPDF_Test[i];
			if(min>ratio) min=ratio;
			if(max<ratio) max=ratio;
			xyzDataSet.set(i, ratio);
		}
		System.out.println("min="+min+"\nmax="+max);
		GridSourceProvider2023_Analysis.plotRatioMap(xyzDataSet, 0.1, folderName);
//		System.exit(0);

	}


	
	public static void main(String[] args) {

		String fileName="/Users/field/workspace/OpenSHA/dev/scratch/UCERF3/data/scratch/InversionSolutions/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip";
		FaultSystemSolution fss;
		try {
			fss = U3FaultSystemIO.loadSol(new File(fileName));
		} catch (Exception e) {
			throw ExceptionUtils.asRuntimeException(e);
		}

		CaliforniaRegions.RELM_TESTING_GRIDDED griddedRegion = RELM_RegionUtils.getGriddedRegionInstance();

		SeisDepthDistribution seisDepthDistribution = new SeisDepthDistribution();
		double delta=2;
		HistogramFunction binnedDepthDistFunc = new HistogramFunction(1d, 12,delta);
		for(int i=0;i<binnedDepthDistFunc.size();i++) {
			double prob = seisDepthDistribution.getProbBetweenDepths(binnedDepthDistFunc.getX(i)-delta/2d,binnedDepthDistFunc.getX(i)+delta/2d);
			binnedDepthDistFunc.set(i,prob);
		}
		System.out.println("Total Depth Prob Sum: "+binnedDepthDistFunc.calcSumOfY_Vals());

		
		double[] spatialPDF = SpatialSeisPDF.UCERF3.getPDF();
		// this sums to 0.9994463999998295; correct it to 1.0
		double sum=0;
		for(double val:spatialPDF) sum+=val;
		for(int i=0;i<spatialPDF.length;i++)
			spatialPDF[i] = spatialPDF[i]/sum;
		
//		// make spatialPDF constant for testing
//		for(int i=0;i<spatialPDF.length;i++)
//			spatialPDF[i] = 1.0/spatialPDF.length;

		//		System.out.println(totGriddedSeisMFD);
//		System.exit(0);
		
		// make totGriddedSeisMFD
		GridSourceProvider gridSrcProviderU3 = fss.getGridSourceProvider();
		IncrementalMagFreqDist tempMFD = gridSrcProviderU3.getMFD(0);
		SummedMagFreqDist totGriddedSeisMFD = new SummedMagFreqDist(tempMFD.getMinX(), tempMFD.size(),tempMFD.getDelta());
		for(int i=0;i<gridSrcProviderU3.size();i++) {
			totGriddedSeisMFD.addIncrementalMagFreqDist(gridSrcProviderU3.getMFD(i));	
		}
		
		GridSourceProvider2023_Analysis.compareU3_InitialAndFinalSpatialPDF(spatialPDF, gridSrcProviderU3, "U3_OrigToFinalSpatialPDF_Ratio");
		
		CubedGriddedRegion cgr = new CubedGriddedRegion(griddedRegion);
		
		double[] fracStrikeSlip, fracReverse, fracNormal;
		fracStrikeSlip = new GridReader("StrikeSlipWts.txt").getValues();
		fracReverse = new GridReader("ReverseWts.txt").getValues();
		fracNormal = new GridReader("NormalWts.txt").getValues();

		
		long startTime = System.currentTimeMillis();

		GridSourceProvider2023 gridProvider = new GridSourceProvider2023(fss, cgr, spatialPDF, totGriddedSeisMFD, 
				binnedDepthDistFunc, fracStrikeSlip, fracReverse, fracNormal);
		
//		GridSourceProvider2023_Analysis.plotGridRateAboveMag(gridProvider, 5.05, "GSP2023_RatesAboveM5");
//		GridSourceProvider2023_Analysis.plotGridRateAboveMag(gridSrcProviderU3, 5.05, "GSP_U3_RatesAboveM5");
		
		// this ratio looks very much like the compareU3_InitialAndFinalSpatialPDF result above, so the difference is due
		// mostly to U3 evenly distributing subseis ruptures among grid celss inside polygons
		GridSourceProvider2023_Analysis.plotGridRateRatioAboveMag(gridProvider, gridSrcProviderU3, 5.05, "GridRateRatioAboveM5");

		//		GridSourceProvider2023_Analysis.plotGridRateRatioAboveMag(gridProvider, gridSrcProviderU3, 7.05, "GridRateRatioAboveM7");
//		GridSourceProvider2023_Analysis.plotGridRateRatioAboveMag(gridSrcProviderU3, gridProvider, 7.05, "GridRateRatioAboveM7_Alt");

						
//		gridProvider.testTotalGriddedSeisMFD();
//		gridProvider.testTotalTrulyOffFaultGriddedSeisMFD();
//		gridProvider.testTotalMgt4_RatesInCells();
//		gridProvider.testSectDistFractWtMapList();
		
//		GridSourceProvider2023_Analysis.plotSupraSeisRateAtDepthMap(gridProvider, 7d, "SupraSeisRatesAtDepth7km");
		
		// this should not exactly equal Figure 9c in U3ETAS BSSA paper (due to how char factors for mult sections in cube are handled)
//		GridSourceProvider2023_Analysis.plotCharFactorAtDepth(gridProvider, 7d, 1d);


		// the following matches the associated U3ETAS map (Figure 5a) close enough
//		GridSourceProvider2023_Analysis.plotCubeRateAboveMagAtDepthMap(gridProvider, 2.55, 7d, "RatesAboveM2pt5_AtDepth7km");

		// the following matches the associated U3ETAS map (Figure 9b) close enough
//		GridSourceProvider2023_Analysis.plotCubeRateAboveMagAtDepthMap(gridProvider, 6.35, 7d, "RatesAboveM6pt3_AtDepth7km");

//		GridSourceProvider2023_Analysis.plotCubeFractionOnAndOffFaultAtDepth(gridProvider,7d);
		
//		GridSourceProvider2023_Analysis.writeCubeLocsAndWtsForForParentSect(gridProvider, 301, "CubeLocWtsForSAF_MojaveSouth.txt");
//		GridSourceProvider2023_Analysis.writeCubeLocsAndWtsForForParentSect(gridProvider, 141, "CubeLocWtsForGreatValley14_KH.txt");
		

//		GridSourceProvider2023_Analysis.plotTotalMFDs(gridProvider);
		
//		GridSourceProvider2023_Analysis.computeHistogramOfNumSectionsInCubes(gridProvider, true);
		
		long runtime = System.currentTimeMillis()-startTime;
		double runtimeMin = runtime/60000d;
		System.out.println("Runtime = "+(float)runtimeMin+" min");

	}
	



}
