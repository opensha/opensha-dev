package scratch.kevin.cybershake.ucerf3.safWallToWallTests;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.dom4j.DocumentException;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.function.EvenlyDiscretizedFunc;
import org.opensha.commons.data.region.CaliforniaRegions;
import org.opensha.commons.data.siteData.SiteData;
import org.opensha.commons.data.siteData.impl.CVM4i26BasinDepth;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.commons.gui.plot.PlotCurveCharacterstics;
import org.opensha.commons.gui.plot.PlotLineType;
import org.opensha.commons.gui.plot.PlotSpec;
import org.opensha.refFaultParamDb.vo.FaultSectionPrefData;
import org.opensha.sha.cybershake.db.CybershakeSite;
import org.opensha.sha.cybershake.db.Cybershake_OpenSHA_DBApplication;
import org.opensha.sha.cybershake.db.DBAccess;
import org.opensha.sha.cybershake.db.SiteInfo2DB;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.faultSurface.RuptureSurface;
import org.opensha.sha.gui.infoTools.IMT_Info;
import org.opensha.sha.imr.AttenRelRef;
import org.opensha.sha.imr.ScalarIMR;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.SiteParams.DepthTo1pt0kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.DepthTo2pt5kmPerSecParam;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import scratch.UCERF3.FaultSystemRupSet;
import scratch.UCERF3.FaultSystemSolution;
import scratch.UCERF3.utils.FaultSystemIO;

public class GMPEComparisonPlotter {

	public static void main(String[] args) throws IOException, DocumentException {
		FaultSystemSolution sol = FaultSystemIO.loadSol(
				new File("/home/kevin/CyberShake/ucerf3/saf_wall_to_wall_tests/saf_subset_sol.zip"));
		FaultSystemRupSet rupSet = sol.getRupSet();
		
		Region soCal = new CaliforniaRegions.RELM_SOCAL();
		
		String siteShortName = "USC";
		
		DBAccess db = Cybershake_OpenSHA_DBApplication.getDB();
		SiteInfo2DB sites2db = new SiteInfo2DB(db);
		CybershakeSite site = sites2db.getSiteFromDB(siteShortName);
		
		Location siteLoc = site.createLocation();
		
		int largestRupIndex = getLargestRupIndex(rupSet);
		
		System.out.println("Largest Rupture is M"+(float)rupSet.getMagForRup(largestRupIndex)+", index="+largestRupIndex);
		printRupSections(rupSet.getFaultSectionDataForRupture(largestRupIndex));
		
		// now find rupture subset in so cal
		int soCalMatch = getLargestSubRupInRegion(rupSet, largestRupIndex, soCal);
		printRupSections(rupSet.getFaultSectionDataForRupture(soCalMatch));
		
		ScalarIMR gmpe = AttenRelRef.NGAWest_2014_AVG_NOIDRISS.instance(null);
//		ScalarIMR gmpe = AttenRelRef.CB_2008.instance(null);
		gmpe.setParamDefaults();
		
		double period = 1d;
		gmpe.setIntensityMeasure(SA_Param.NAME);
		SA_Param.setPeriodInSA_Param(gmpe.getIntensityMeasure(), period);
		
		WillsMap2006 wills = new WillsMap2006();
		CVM4i26BasinDepth z10 = new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_1_0);
		CVM4i26BasinDepth z25 = new CVM4i26BasinDepth(SiteData.TYPE_DEPTH_TO_2_5);
		
		Site gmpeSite = new Site(siteLoc);
		gmpeSite.addParameterList(gmpe.getSiteParams());
		gmpeSite.getParameter(Double.class, Vs30_Param.NAME).setValue(wills.getValue(siteLoc));
		if (gmpeSite.containsParameter(DepthTo1pt0kmPerSecParam.NAME))
			gmpeSite.getParameter(Double.class, DepthTo1pt0kmPerSecParam.NAME).setValue(z10.getValue(siteLoc));
		if (gmpeSite.containsParameter(DepthTo2pt5kmPerSecParam.NAME))
			gmpeSite.getParameter(Double.class, DepthTo2pt5kmPerSecParam.NAME).setValue(z25.getValue(siteLoc));
		
		
		DiscretizedFunc xVals = IMT_Info.getUSGS_SA_Function();
		DiscretizedFunc logXVals = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			logXVals.set(Math.log(xVals.getX(i)), 0d);
		
		gmpe.setSite(gmpeSite);
		gmpe.setEqkRupture(getRup(rupSet, largestRupIndex));
		
		double fullMean = gmpe.getMean();
		double fullStdDev = gmpe.getStdDev();
		DiscretizedFunc fullExceed = getUnLogged(xVals, gmpe.getExceedProbabilities(logXVals));
		fullExceed.setName("Full Rupture");
		
		gmpe.setEqkRupture(getRup(rupSet, soCalMatch));
		
		double subsetMean = gmpe.getMean();
		double subsetStdDev = gmpe.getStdDev();
		DiscretizedFunc subExceed = getUnLogged(xVals, gmpe.getExceedProbabilities(logXVals));
		fullExceed.setName("Subset Rupture");
		
		System.out.println("Original Rupture");
		System.out.println("\tmean="+Math.exp(fullMean));
		System.out.println("\tstdDev="+fullStdDev);
		System.out.println("Subset Rupture");
		System.out.println("\tmean="+Math.exp(subsetMean));
		System.out.println("\tstdDev="+subsetStdDev);
		
		List<DiscretizedFunc> funcs = Lists.newArrayList();
		List<PlotCurveCharacterstics> chars = Lists.newArrayList();
		
		funcs.add(fullExceed);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLACK));
		
		funcs.add(subExceed);
		chars.add(new PlotCurveCharacterstics(PlotLineType.SOLID, 3f, Color.BLUE));
		
		PlotSpec spec = new PlotSpec(funcs, chars, "Exceedance Probs", (float)period+"s S(a)", "Probability of Exceedance");
		GraphWindow gw = new GraphWindow(spec);
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
		gw.setVisible(true);
		
		db.destroy();
	}
	
	static DiscretizedFunc getUnLogged(DiscretizedFunc xVals, DiscretizedFunc logFunc) {
//		EvenlyDiscretizedFunc ret = new EvenlyDiscretizedFunc(xVals.getMinX(), xVals.getMaxX(), xVals.size());
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		for (int i=0; i<xVals.size(); i++)
			ret.set(xVals.getX(i), logFunc.getY(i));
		return ret;
	}
	
	static void printRupSections(List<FaultSectionPrefData> sects) {
		String firstForParent = null;
		String prevParent = null;
		String prevSub = null;
		for (FaultSectionPrefData sect : sects) {
			String parent = sect.getParentSectionName();
			if (!parent.equals(prevParent)) {
				// new parent section
				if (firstForParent != null)
					System.out.println("\t"+firstForParent+"\t=>\t"+prevSub);
				firstForParent = sect.getName();
			}
			
			prevParent = parent;
			prevSub = sect.getSectionName();
		}
		if (firstForParent != null)
			System.out.println("\t"+firstForParent+"\t=>\t"+prevSub);
	}
	
	static int getLargestRupIndex(FaultSystemRupSet rupSet) {
		int largestRupIndex = -1;
		double largestRupMag = 0d;
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			double mag = rupSet.getMagForRup(rupIndex);
			if (mag > largestRupMag) {
				largestRupMag = mag;
				largestRupIndex = rupIndex;
			}
		}
		return largestRupIndex;
	}
	
	static int getLargestSubRupInRegion(FaultSystemRupSet rupSet, int origRupIndex, Region region) {
		List<FaultSectionPrefData> origSects = rupSet.getFaultSectionDataForRupture(origRupIndex);
		HashSet<Integer> allowedSectionIndices = new HashSet<Integer>();
		for (FaultSectionPrefData sect : origSects) {
			boolean inside = true;
			for (Location loc : sect.getFaultTrace())
				if (!region.contains(loc))
					inside = false;
			if (inside)
				allowedSectionIndices.add(sect.getSectionId());
		}
		System.out.println(allowedSectionIndices.size()+"/"+origSects.size()+" of those sub sections are within region");
		
		double largestMag = 0d;
		int bestMatch = -1;
		rupLoop:
		for (int rupIndex=0; rupIndex<rupSet.getNumRuptures(); rupIndex++) {
			double mag = rupSet.getMagForRup(rupIndex);
			if (mag < largestMag)
				continue;
			for (int sectIndex : rupSet.getSectionsIndicesForRup(rupIndex)) {
				if (!allowedSectionIndices.contains(sectIndex))
					continue rupLoop;
			}
			// all sections are a match and it's the largest if we make it here
			largestMag = mag;
			bestMatch = rupIndex;
		}
		System.out.println("Best subset match is M"+(float)largestMag+", index="+bestMatch);
		Preconditions.checkState(bestMatch >= 0);
		return bestMatch;
	}
	
	static EqkRupture getRup(FaultSystemRupSet rupSet, int rupIndex) {
		RuptureSurface surf = rupSet.getSurfaceForRupupture(rupIndex, 1d, false);
		return new EqkRupture(rupSet.getMagForRup(rupIndex), rupSet.getAveRakeForRup(rupIndex), surf, null);
	}

}
