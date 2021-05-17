package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.dom4j.Document;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.OrderedSiteDataProviderList;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.SiteDataValue;
import org.opensha.commons.data.siteData.SiteDataValueList;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.USGSBayAreaBasinDepth;
import org.opensha.commons.data.siteData.impl.WaldAllenGlobalVs30;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.data.siteData.impl.WillsMap2015;
import org.opensha.commons.geo.GriddedRegion;
import org.opensha.commons.geo.Location;
import org.opensha.commons.hpc.JavaShellScriptWriter;
import org.opensha.commons.hpc.mpj.FastMPJShellScriptWriter;
import org.opensha.commons.hpc.mpj.MPJExpressShellScriptWriter;
import org.opensha.commons.hpc.pbs.BatchScriptWriter;
import org.opensha.commons.hpc.pbs.StampedeScriptWriter;
import org.opensha.commons.hpc.pbs.USC_HPCC_ScriptWriter;
import org.opensha.commons.param.Parameter;
import org.opensha.commons.util.XMLUtils;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.util.SiteTranslator;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class MPJ_UCERF3_ShakeMapPrecalcScriptGen {
	
	public static final DateFormat df = new SimpleDateFormat("yyyy_MM_dd");

	public static void main(String[] args) throws IOException {
		File localMainDir = new File("/home/kevin/OpenSHA/UCERF3/shakemap_precalc");
		
		double gridSpacing = 0.02;
		GriddedRegion reg = new CaliforniaRegions.RELM_TESTING_GRIDDED(gridSpacing);
		boolean siteEffects = true;
		boolean basinDepth = true;
		double distCutoff = 200;
		String imts = "PGA,PGV";
		
		boolean stampede = false;
//		int threads = 1;
//		String queue = null;
		int threads = 20;
		String queue = "scec";
		
		int nodes = 34;
		int hours = 24;
		
		AttenuationRelationship gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
		gmpe.setParamDefaults();
		
		ArrayList<SiteData<?>> siteData = null;
		if (siteEffects)
			siteData = getSiteDataProviders(basinDepth);
		
		String dateStr = df.format(new Date());
		
		String jobName = dateStr+"-"+gmpe.getShortName()+"-spacing"+(float)gridSpacing;
		if (siteEffects) {
			jobName += "-site-effects";
			if (basinDepth)
				jobName += "-with-basin";
			else
				jobName += "-no-basin";
		} else {
			jobName += "-no-site-effects";
		}
		
		File localDir = new File(localMainDir, jobName);
		Preconditions.checkState(localDir.exists() || localDir.mkdir());
		
		System.out.println("Job name: "+jobName);
		System.out.println("Dir: "+localDir.getAbsolutePath());
		
		File sitesFile = new File(localDir, "sites.xml");
		writeSitesXML(reg, siteData, gmpe, sitesFile);
		
		// write gmpe file
		File gmpeFile = new File(localDir, gmpe.getShortName()+".xml");
		writeGMPEFile(gmpe, gmpeFile);
		
		File remoteDir, remoteSolFile;
		JavaShellScriptWriter mpjWrite;
		BatchScriptWriter pbsWrite;
		
		int memGigs;
		int ppn;
		if (stampede) {
			memGigs = 26;
			ppn = 16;
//			remoteDir = new File("/work/00950/kevinm/ucerf3/etas_sim");
			remoteDir = new File("/work/00950/kevinm/ucerf3/shakemap_precalc");
			remoteSolFile = new File("/work/00950/kevinm/ucerf3/inversion/compound_plots/2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
			mpjWrite = new FastMPJShellScriptWriter(StampedeScriptWriter.JAVA_BIN, memGigs*1024,
					null, StampedeScriptWriter.FMPJ_HOME);
			((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			pbsWrite = new StampedeScriptWriter();
		} else {
			if (queue == null) {
				memGigs = 9;
				ppn = 8;
			} else {
				memGigs = 60;
				ppn = 20;
			}
			remoteDir = new File("/home/scec-02/kmilner/ucerf3/shakemap_precalc");
			remoteSolFile = new File("/home/scec-02/kmilner/ucerf3/inversion_compound_plots/"
					+ "2013_05_10-ucerf3p3-production-10runs/"
					+ "2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL_FM3_1_SpatSeisU3_MEAN_BRANCH_AVG_SOL.zip");
			boolean fmpj = nodes < 25;
			fmpj = false;
			if (fmpj) {
				mpjWrite = new FastMPJShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
						null, USC_HPCC_ScriptWriter.FMPJ_HOME);
				((FastMPJShellScriptWriter)mpjWrite).setUseLaunchWrapper(true);
			} else {
				mpjWrite = new MPJExpressShellScriptWriter(USC_HPCC_ScriptWriter.JAVA_BIN, memGigs*1024,
						null, USC_HPCC_ScriptWriter.MPJ_HOME);
			}
			pbsWrite = new USC_HPCC_ScriptWriter();
		}
		File remoteJobDir = new File(remoteDir, jobName);
		
		List<File> classpath = new ArrayList<File>();
		classpath.add(new File(remoteDir, "commons-cli-1.2.jar"));
		classpath.add(new File(remoteJobDir, "OpenSHA_complete.jar"));
		mpjWrite.setClasspath(classpath);
		
		String argz;
		if (threads > 0)
			argz = "--threads "+threads;
		else
			argz = "";
		argz += " --distance-cutoff "+distCutoff;
		argz += " --gmpe-file "+new File(remoteJobDir, gmpeFile.getName()).getAbsolutePath();
		argz += " --sites-file "+new File(remoteJobDir, sitesFile.getName()).getAbsolutePath();
		argz += " --solution-file "+remoteSolFile.getAbsolutePath();
		argz += " --output-dir "+remoteJobDir.getAbsolutePath();
		argz += " --imts "+imts;
		
		List<String> script = mpjWrite.buildScript(MPJ_UCERF3_ShakeMapPrecalc.class.getName(), argz);
		
		int mins = hours*60;
		script = pbsWrite.buildScript(script, mins, nodes, ppn, queue);
		pbsWrite.writeScript(new File(localDir, jobName+".pbs"), script);
	}
	
	public static ArrayList<SiteData<?>> getSiteDataProviders(boolean basinDepth) throws IOException {
		ArrayList<SiteData<?>> siteData = new ArrayList<SiteData<?>>();
		siteData.add(new WillsMap2015());
		siteData.add(new WaldAllenGlobalVs30());
		if (basinDepth) {
			siteData.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			siteData.add(new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
			siteData.add(new USGSBayAreaBasinDepth(SiteData.TYPE_DEPTH_TO_2_5));
			siteData.add(new USGSBayAreaBasinDepth(SiteData.TYPE_DEPTH_TO_1_0));
		}
		return siteData;
	}
	
	public static void writeSitesXML(GriddedRegion reg, ArrayList<SiteData<?>> siteData, AttenuationRelationship gmpe, File sitesFile)
			throws IOException {
		// prepare sites
		List<Site> sites = Lists.newArrayList();
		for (Location loc : reg.getNodeList()) {
			Site site = new Site(loc);
			for (Parameter<?> param : gmpe.getSiteParams())
				site.addParameter((Parameter<?>) param.clone());
			sites.add(site);
		}
		if (siteData != null) {
			System.out.print("Fetching site data...");
			OrderedSiteDataProviderList provs = new OrderedSiteDataProviderList(siteData);
			List<SiteDataValueList<?>> vals  = provs.getAllAvailableData(sites);
			System.out.println("done.");
			SiteTranslator trans = new SiteTranslator();
			System.out.print("Setting site params...");
			for (int i=0; i<sites.size(); i++) {
				ArrayList<SiteDataValue<?>> siteVals = new ArrayList<SiteDataValue<?>>();
				for (SiteDataValueList<?> valList : vals) {
					siteVals.add(valList.getValue(i));
				}
				for (Parameter<?> param : sites.get(i))
					trans.setParameterValue(param, siteVals);
			}
			System.out.println("done.");
		}

//		File sitesFile = new File(localDir, "sites.xml");
		System.out.print("Writing sites XML...");
		Document doc = XMLUtils.createDocumentWithRoot();
		Site.writeSitesToXML(sites, doc.getRootElement());
		XMLUtils.writeDocumentToFile(sitesFile, doc);
		System.out.println("done.");
	}
	
	public static void writeGMPEFile(AttenuationRelationship gmpe, File gmpeFile) throws IOException {
		Document doc = XMLUtils.createDocumentWithRoot();
		gmpe.toXMLMetadata(doc.getRootElement());
		XMLUtils.writeDocumentToFile(gmpeFile, doc);
	}

}
