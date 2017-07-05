package scratch.kevin.cybershake.ugms;

import java.io.File;
import java.util.Map;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.data.siteData.impl.WillsMap2006;
import org.opensha.commons.geo.Location;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveWriter;
import org.opensha.sha.calc.hazardMap.components.BinaryCurveArchiver;
import org.opensha.sha.calc.mcer.ASCEDetLowerLimitCalc;
import org.opensha.sha.cybershake.calc.mcer.MCERDataProductsCalc;
import org.opensha.sha.cybershake.calc.mcer.UGMS_WebToolCalc;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;
import com.google.common.io.Files;

public class GMPE_DetLowerLimitRecalc {

	public static void main(String[] args) throws Exception {
		File dir = new File("/home/kevin/CyberShake/MCER/gmpe_cache_gen/2016_09_30-ucerf3_downsampled_ngaw2_binary_0.02_Wills");
		String prefix = "NGAWest_2014_NoIdr_MeanUCERF3_downsampled_RotD100_";
		
		File mcerFile = new File(dir, prefix+"mcer.bin");
		File origMCERFile = new File(dir, prefix+"mcer_orig.bin");
		File probFile = new File(dir, prefix+"prob.bin");
		File detFile = new File(dir, prefix+"det.bin");
		
		Preconditions.checkState(mcerFile.exists());
		Preconditions.checkState(probFile.exists());
		Preconditions.checkState(detFile.exists());
		
		// backup MCER file
		Files.copy(mcerFile, origMCERFile);
		
		BinaryHazardCurveReader reader = new BinaryHazardCurveReader(probFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> probCurves = reader.getCurveMap();
		reader = new BinaryHazardCurveReader(detFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> detCurves = reader.getCurveMap();
		reader = new BinaryHazardCurveReader(origMCERFile.getAbsolutePath());
		Map<Location, ArbitrarilyDiscretizedFunc> origMCERCurves = reader.getCurveMap();
		
		Map<Location, DiscretizedFunc> newMCERCurves = Maps.newHashMap();
		
		double vs30 = Double.NaN;
		WillsMap2006 wills = null;
		String dirName = dir.getName();
		if (dirName.contains("class")) {
			String cl = dirName.substring(dirName.indexOf("class")+5);
			vs30 = UGMS_WebToolCalc.vs30Map.get(cl);
			System.out.println("Detected site class "+cl+", vs30="+vs30);
		} else {
			wills = new WillsMap2006();
		}
		
		int locsChanged = 0;
		int totLocs = 0;
		int ptsChanged = 0;
		int totPts = 0;
		
		for (Location loc : probCurves.keySet()) {
			ArbitrarilyDiscretizedFunc prob = probCurves.get(loc);
			ArbitrarilyDiscretizedFunc det = detCurves.get(loc);
			ArbitrarilyDiscretizedFunc origMCER = origMCERCurves.get(loc);
			
			if (wills != null)
				vs30 = wills.getValue(loc);
			
			DiscretizedFunc detLowerLimit = ASCEDetLowerLimitCalc.calc(prob, vs30, loc);
			
			DiscretizedFunc mcer = MCERDataProductsCalc.calcMCER(det, prob, detLowerLimit);
			
			boolean changed = false;
			Preconditions.checkState(mcer.size() == origMCER.size());
			for (int i=0; i<mcer.size(); i++) {
				double origY = origMCER.getY(i);
				double newY = mcer.getY(i);
				
				if ((float)origY != (float)newY) {
					changed = true;
					ptsChanged++;
				}
				totPts++;
			}
			if (changed)
				locsChanged++;
			totLocs++;
			
			newMCERCurves.put(loc, mcer);
		}
		BinaryHazardCurveWriter writer = new BinaryHazardCurveWriter(mcerFile);
		writer.writeCurves(newMCERCurves);
		
		float locsChangedPer = 100f*locsChanged/totLocs;
		System.out.println(locsChanged+"/"+totLocs+" ("+locsChangedPer+" %) locations changed");
		float ptsChangedPer = 100f*ptsChanged/totPts;
		System.out.println(ptsChanged+"/"+totPts+" ("+ptsChangedPer+" %) points changed");
	}

}
