package scratch.kevin.nshm23.figures;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.geo.Location;
import org.opensha.commons.geo.Region;
import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_DeformationModels;
import org.opensha.sha.earthquake.rupForecastImpl.nshm23.logicTree.NSHM23_FaultModels;
import org.opensha.sha.faultSurface.FaultSection;

import com.google.common.base.Preconditions;

public class GarlockMojaveECSZ_DMsPlot {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/garlock_mojave_ecsz");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		NSHM23_DeformationModels[] dms = {
				NSHM23_DeformationModels.GEOLOGIC, // must be first
				NSHM23_DeformationModels.EVANS,
				NSHM23_DeformationModels.POLLITZ,
				NSHM23_DeformationModels.SHEN_BIRD,
				NSHM23_DeformationModels.ZENG
		};
		NSHM23_DeformationModels refDM = NSHM23_DeformationModels.GEOLOGIC;
		NSHM23_FaultModels fm = NSHM23_FaultModels.WUS_FM_v3;
		
		Region reg = new Region(new Location(34, -119), new Location(37, -116));
		
		List<? extends FaultSection> refSubSects = refDM.build(fm);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(refSubSects);
		mapMaker.setRegion(reg);
		
		CPT linearSlipCPT = GMT_CPT_Files.SEQUENTIAL_BATLOW_UNIFORM.instance().rescale(0d, 10d);
		CPT logSlipCPT = GMT_CPT_Files.RAINBOW_UNIFORM.instance().rescale(-2, Math.log10(30));
		logSlipCPT.setNanColor(Color.LIGHT_GRAY);
		CPT pDiffCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(-60d, 60d).trim(-50, 50);
		
		double[] refSlips = null;
		
		for (NSHM23_DeformationModels dm : dms) {
			List<? extends FaultSection> subSects;
			if (dm == refDM)
				subSects = refSubSects;
			else
				subSects = dm.build(fm);
			
			double[] slips = new double[subSects.size()];
			double[] pDiffs = dm == refDM ? null : new double[subSects.size()];
			for (int s=0; s<subSects.size(); s++) {
				double slip = subSects.get(s).getOrigAveSlipRate();
				slips[s] = slip;
				if (pDiffs != null)
					pDiffs[s] = 100d*(slip - refSlips[s])/refSlips[s];
			}
			double[] logSlips = new double[subSects.size()];
			for (int s=0; s<slips.length; s++)
				logSlips[s] = slips[s] > 0 ? Math.log10(slips[s]) : Double.NaN;
			if (dm == refDM)
				refSlips = slips;
			
			String prefix = dm.getFilePrefix();
			
			mapMaker.plotSectScalars(slips, linearSlipCPT, dm.getShortName()+" Slip Rate (mm/yr)");
			mapMaker.plot(outputDir, prefix+"_slips", " ");
			
			mapMaker.plotSectScalars(logSlips, logSlipCPT, "Log10["+dm.getShortName()+" Slip Rate (mm/yr)]");
			mapMaker.plot(outputDir, prefix+"_slips_log", " ");
			
			if (pDiffs != null) {
				mapMaker.plotSectScalars(pDiffs, pDiffCPT, dm.getShortName()+" vs "+refDM.getShortName()+", % Change");
				mapMaker.plot(outputDir, prefix+"_slips_pDiff", " ");
			}
		}
	}

}
