package scratch.kevin.portfolioLEC;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.dom4j.DocumentException;
import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.LocationUtils;
import org.opensha.commons.param.Parameter;
import org.opensha.sha.calc.hazardMap.HazardDataSetLoader;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.EpistemicListERF;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2_TimeDependentEpistemicList;
import org.opensha.sha.imr.AttenRelRef;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class HazardBranchesPostProcess {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 * @throws DocumentException 
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, DocumentException {
		File dir = new File(args[0]);
		File regionFile = new File(dir.getParent(), "regions.kml");
		File outputFile = new File(dir, dir.getName()+".csv");
		
		if (outputFile.exists())
			System.exit(0);
		
		//Line, IMR, Deformation Model, A-Fault Solution Type,  Wt On A-Priori Rates, Mag-Area Relationship,
		//B-Faults b-value, Connect More B Faults?, Probability Model, Aperiodocity, IMT, Exceedance Probability, Lat, Lon, IML
		CSVFile<String> csv = new CSVFile<String>(true);
		
		UCERF2_TimeDependentEpistemicList erfList = new UCERF2_TimeDependentEpistemicList();
		
		int erfIndex = Integer.parseInt(dir.getName().split("_")[1]);
		
		HashMap<String, Integer> paramNamesMap = new HashMap<String, Integer>();
		ArrayList<String> paramNamesList = new ArrayList<String>();
		for (int i=0; i<erfList.getNumERFs(); i++) {
			ERF erf = erfList.getERF_noUpdate(i);
			
			System.out.println("ERF "+i);
			
			for (Parameter<?> param : erf.getAdjustableParameterList()) {
				if (!paramNamesMap.containsKey(param.getName())) {
					paramNamesMap.put(param.getName(), paramNamesMap.size());
					paramNamesList.add(param.getName());
				}
			}
		}
		
		List<String> header = Lists.newArrayList();
		header.add("IMR");
		header.addAll(paramNamesList);
		header.add("IMT");
		header.add("Exceed. Prob.");
		header.add("Latitude");
		header.add("Longitude");
		header.add("IML");
		
		csv.addLine(header);
		
		ERF erf = erfList.getERF(erfIndex);
		
		List<String> paramValues = Lists.newArrayList();
		for (int i=0; i<paramNamesList.size(); i++)
			paramValues.add("");
		
		for (int i=0; i<paramNamesList.size(); i++) {
			String paramName = paramNamesList.get(i);
			
			String value = null;
			
			try {
				value = erf.getAdjustableParameterList().getParameter(paramName).getValue().toString();
			} catch (ParameterException e) {
				continue;
			}
			if (value != null)
				paramValues.set(paramNamesMap.get(paramName), value);
		}
		
		List<GriddedRegion> regions = HazardMapLogicTreeInRegionsGen.getRegions(regionFile);
		
		List<Location> locs = Lists.newArrayList();
		
		for (GriddedRegion region : regions)
			locs.addAll(region.getNodeList());
		
		Map<Location, Location> locsMapping = null;
		
		AttenRelRef[] imrRefs = { AttenRelRef.CB_2008, AttenRelRef.BA_2008, AttenRelRef.CY_2008, AttenRelRef.AS_2008 };
		String[] imts = { "PGA", "SA03", "SA10", "PGV" };
		double[] exceed_probs = { 0.5, 0.1, 0.02 };
		String[] exceed_strs = { "50%/50yr", "10%/50yr", "2%/50yr" };
		
		int branchCount = 0;
		for (AttenRelRef imrRef : imrRefs) {
			for (String imt : imts) {
				branchCount++;
				
				File curveDir = new File(dir, "imrs"+branchCount);
				
				System.out.println("Parsing "+curveDir.getName()+" ("+imrRef.name()+", "+imt+")");
				Map<Location, DiscretizedFunc> funcsMap = loadFuncs(curveDir);
				
				if (locsMapping == null) {
					System.out.print("Mapping file locs to grid locs...");
					locsMapping = Maps.newHashMap();
					for (Location loc : locs) {
						double minDist = 1d;
						Location closestLoc = null;
						for (Location dirLoc : funcsMap.keySet()) {
							double dist = LocationUtils.horzDistanceFast(loc, dirLoc);
							if (dist < minDist) {
								minDist = dist;
								closestLoc = dirLoc;
							}
						}
						Preconditions.checkNotNull(closestLoc,
								"No mapping found within "+minDist+"km from region loc: "+loc);
						locsMapping.put(loc, closestLoc);
					}
					System.out.println("DONE.");
				}
				
				List<String> linePrefix = Lists.newArrayList();
				linePrefix.add(imrRef.name());
				linePrefix.addAll(paramValues);
				linePrefix.add(imt);
				
				for (Location loc : locs) {
					DiscretizedFunc func = funcsMap.get(locsMapping.get(loc));
					Preconditions.checkNotNull(func, "No curve for loc: "+loc);
					
					for (int e=0; e<exceed_probs.length; e++) {
						List<String> line = Lists.newArrayList(linePrefix);
						
						line.add(exceed_strs[e]);
						double prob = exceed_probs[e];
						
						double val = HazardDataSetLoader.getCurveVal(func, false, prob);
						
						line.add((float)loc.getLatitude()+"");
						line.add((float)loc.getLongitude()+"");
						line.add(val+"");
						
						csv.addLine(line);
					}
				}
			}
		}
		
		csv.writeToFile(outputFile);
	}
	
	public static Map<Location, DiscretizedFunc> loadFuncs(File curveDir) throws FileNotFoundException, IOException {
		Map<Location, DiscretizedFunc> funcsMap = Maps.newHashMap();
		
		// for each file in the list
		for(File dir : curveDir.listFiles()){
			// make sure it's a subdirectory
			if (dir.isDirectory() && !dir.getName().endsWith(".")) {
				File[] subDirList=dir.listFiles();
				for(File file : subDirList) {
					//only taking the files into consideration
					if(file.isFile()){
						String curveFileName = file.getName();
						//files that ends with ".txt"
						if(curveFileName.endsWith(".txt")){
							Location loc = HazardDataSetLoader.decodeFileName(curveFileName);
							
							funcsMap.put(loc, ArbitrarilyDiscretizedFunc.loadFuncFromSimpleFile(file.getAbsolutePath()));
						}
					}
				}
			}
		}
		
		return funcsMap;
	}

}
