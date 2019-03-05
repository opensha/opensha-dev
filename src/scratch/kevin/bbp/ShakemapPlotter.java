package scratch.kevin.bbp;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.List;
import java.util.zip.ZipException;

import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.xyz.GriddedGeoDataSet;
import org.opensha.commons.exceptions.GMT_MapException;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.mapping.gmt.GMT_Map;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.XMLUtils;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;

import com.google.common.base.Preconditions;

import scratch.UCERF3.analysis.FaultBasedMapGen;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SimZipLoader.BBP_ShakeMapSimZipLoader;

public class ShakemapPlotter {
	
	public static void plotShakemaps(BBP_ShakeMapSimZipLoader loader, GriddedRegion gridReg, List<BBP_Site> sites,
			String label, File outputDir, String prefix, boolean log, double... periods)
					throws IOException, GMT_MapException {
		plotShakemaps(loader, gridReg, sites, label, outputDir, prefix, log, null, null, periods);
	}
	
	public static void plotShakemaps(BBP_ShakeMapSimZipLoader loader, GriddedRegion gridReg, List<BBP_Site> sites,
			String label, File outputDir, String prefix, boolean log, ScalarIMR gmpe, EqkRupture rup,
			double... periods) throws IOException, GMT_MapException {
		GriddedGeoDataSet[] xyzs = load(loader, gridReg, sites, periods, false);
		GriddedGeoDataSet[] rd100xyzs = null;
		if (loader.hasRotD100())
			rd100xyzs = load(loader, gridReg, sites, periods, true);
		
		GriddedGeoDataSet[] gmpes = null;
		if (gmpe != null)
			gmpes = calcGMPE(gmpe, rup, gridReg, sites, periods);
		
		CPT cpt = GMT_CPT_Files.MAX_SPECTRUM.instance();
		
		for (int p=0; p<periods.length; p++) {
			System.out.println("Plotting "+(float)periods[p]+"s map");
			double max = xyzs[p].getMaxZ();
			if (gmpes != null)
				max = Math.max(max, gmpes[p].getMaxZ());
			GriddedGeoDataSet origXYZ = xyzs[p];
			GriddedGeoDataSet origGMPE = null;
			if (gmpes != null)
				origGMPE = gmpes[p];
			if (log) {
				origXYZ = origXYZ.copy();
				xyzs[p].log10();
				if (gmpes != null) {
					origGMPE = origGMPE.copy();
					gmpes[p].log10();
				}
				double logMax = Math.log10(max);
				double cleanLogMax = Math.ceil(logMax);
				if (cleanLogMax - logMax > 0.5)
					cleanLogMax -= 0.5;
				cpt = cpt.rescale(-3d, cleanLogMax);
			} else {
				max = Math.ceil(max*5d)/5d;
				cpt = cpt.rescale(0d, max);
			}
			GMT_Map map = new GMT_Map(gridReg, xyzs[p], gridReg.getSpacing(), cpt);
			map.setRescaleCPT(false);
			map.setLogPlot(false);
			map.setUseGMTSmoothing(false);
			map.setBlackBackground(false);
			map.setJPGFileName(null);
			map.setPDFFileName(null);
			Double mapInterval = null;
			if (log) {
				double cptDelta = cpt.getMaxValue() - cpt.getMinValue();
				if (cptDelta > 2)
					mapInterval = 1d;
				else
					mapInterval = 0.5;
				map.setCPTCustomInterval(mapInterval);
			}
			
			String pStr;
			if (periods[p] == Math.round(periods[p]))
				pStr = (int)periods[p]+"";
			else
				pStr = (float)periods[p]+"";
			
			map.setCustomLabel("Log10 "+label+" "+pStr+"s SA (RotD50)");
			String myPrefix = prefix+"_"+pStr+"s";
			
			FaultBasedMapGen.plotMap(outputDir, myPrefix, false, map);
			
			if (rd100xyzs != null) {
				map.setCustomLabel("Log10 "+label+" "+pStr+"s SA (RotD100)");
				myPrefix = prefix+"_"+pStr+"s_rd100";
				GriddedGeoDataSet origRD_XYZ = rd100xyzs[p];
				if (log) {
					origRD_XYZ = origRD_XYZ.copy();
					rd100xyzs[p].log10();
				}
				map.setCPTCustomInterval(mapInterval);
				map.setGriddedData(rd100xyzs[p]);
				
				FaultBasedMapGen.plotMap(outputDir, myPrefix, false, map);
				
//				0.375	255	255	255	0.625	255	255	127
//				0.625	255	255	0	0.875	255	0	0
//				0.875	255	0	0	1.00	127	0	0
				CPT ratioCPT = new CPT(1d, 1.5, Color.WHITE, new Color(255, 255, 127),
						new Color(255, 255, 0), Color.RED, new Color(127, 0, 0));
				
				GriddedGeoDataSet ratioXYZ = new GriddedGeoDataSet(gridReg, false);
				for (int i=0; i<ratioXYZ.size(); i++)
					ratioXYZ.set(i, origRD_XYZ.get(i)/origXYZ.get(i));
				
				map.setCustomLabel(label+" "+pStr+"s SA RotD100/RotD50 Ratio");
				myPrefix = prefix+"_"+pStr+"s_rd100_ratio";
				map.setCPTCustomInterval(0.1);
				map.setGriddedData(ratioXYZ);
				map.setCpt(ratioCPT);
				
				FaultBasedMapGen.plotMap(outputDir, myPrefix, false, map);
			}
			
			if (gmpe != null) {
				// plot GMPE
				map.setGriddedData(gmpes[p]);
				map.setCpt(cpt);
				map.setCustomLabel("Log10 "+gmpe.getShortName()+" "+pStr+"s SA (RotD50)");
				map.setCPTCustomInterval(mapInterval);
				myPrefix = prefix+"_"+pStr+"s_"+gmpe.getShortName();
				
				FaultBasedMapGen.plotMap(outputDir, myPrefix, false, map);
				
				// now ratio
				CPT ratioCPT = GMT_CPT_Files.GMT_POLAR.instance().rescale(-2, 2);
				GriddedGeoDataSet ratioData = new GriddedGeoDataSet(gridReg, false);
				for (int i=0; i<ratioData.size(); i++)
					ratioData.set(i, Math.log10(origXYZ.get(i)/origGMPE.get(i)));
				
				map.setGriddedData(ratioData);
				map.setCustomLabel("Log10 Ratio "+pStr+"s SA (RotD50)");
				map.setCpt(ratioCPT);
				map.setCPTCustomInterval(0.5);
				myPrefix = prefix+"_"+pStr+"s_"+gmpe.getShortName()+"_ratio";
				
				FaultBasedMapGen.plotMap(outputDir, myPrefix, false, map);
			}
		}
	}
	
	private static GriddedGeoDataSet[] load(BBP_ShakeMapSimZipLoader loader, GriddedRegion gridReg,
			List<BBP_Site> sites, double[] periods, boolean rotD100) throws IOException {
		Preconditions.checkState(gridReg.getNodeCount() == sites.size(), "Sites/gridded region mismatch!");
		GriddedGeoDataSet[] ret = new GriddedGeoDataSet[periods.length];
		for (int i=0; i<periods.length; i++)
			ret[i] = new GriddedGeoDataSet(gridReg, false);
		
		System.out.println("Loading shakemaps for "+periods.length+" periods...");
		for (int i=0; i<sites.size(); i++) {
			DiscretizedFunc spectra;
			if (rotD100)
				spectra = loader.readRotD100(i);
			else
				spectra = loader.readRotD50(i);
			
			for (int p=0; p<periods.length; p++)
				ret[p].set(i, spectra.getInterpolatedY(periods[p]));
		}
		System.out.println("DONE");
		
		return ret;
	}
	
	private static GriddedGeoDataSet[] calcGMPE(ScalarIMR gmpe, EqkRupture rup,
			GriddedRegion gridReg, List<BBP_Site> sites, double[] periods) {
		gmpe.setParamDefaults();
		
		gmpe.setEqkRupture(rup);
		gmpe.setIntensityMeasure(SA_Param.NAME);
		
		Preconditions.checkState(gridReg.getNodeCount() == sites.size(), "Sites/gridded region mismatch!");
		GriddedGeoDataSet[] ret = new GriddedGeoDataSet[periods.length];
		for (int i=0; i<periods.length; i++)
			ret[i] = new GriddedGeoDataSet(gridReg, false);
		
		System.out.println("Calculating GMPE shakemaps for "+periods.length+" periods...");
		for (int i=0; i<sites.size(); i++) {
			Site site = sites.get(i).buildGMPE_Site();
			gmpe.setSite(site);
			
			for (int p=0; p<periods.length; p++) {
				SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), periods[p]);
				ret[p].set(i, Math.exp(gmpe.getMean()));
			}
		}
		System.out.println("DONE");
		
		return ret;
	}
	
	public static GriddedRegion loadGriddedRegion(File sitesXML) throws MalformedURLException, DocumentException {
		Document doc = XMLUtils.loadDocument(sitesXML);
		Element root = doc.getRootElement();
		Element regEl = root.element(GriddedRegion.XML_METADATA_NAME);
		return GriddedRegion.fromXMLMetadata(regEl);
	}

	public static void main(String[] args) throws ZipException, IOException, DocumentException, GMT_MapException {
		File dir = new File("/data/kevin/bbp/parallel/2017_10_12-rundir2194_long-event136704-shakemap-noHF");
		File zipFile = new File(dir, "results.zip");
		List<BBP_Site> sites = BBP_Site.readFile(dir);
		BBP_ShakeMapSimZipLoader loader = new BBP_ShakeMapSimZipLoader(zipFile, sites);
		
		File sitesXML = new File(dir, "sites.xml");
		GriddedRegion gridReg = loadGriddedRegion(sitesXML);
		
		FaultBasedMapGen.LOCAL_MAPGEN = true;
		
		plotShakemaps(loader, gridReg, sites, "Test Map", new File("/tmp"), "shakemap", true, 1d, 2d, 5d);
	}

}
