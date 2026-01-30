package scratch.kevin.nshm26;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.opensha.commons.gui.plot.GeographicMapMaker;
import org.opensha.commons.mapping.gmt.elements.GMT_CPT_Files;
import org.opensha.commons.util.cpt.CPT;
import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.faultSurface.GeoJSONFaultSection;

import com.google.common.base.Preconditions;

public class SlipProjectionTests {

	public static void main(String[] args) throws IOException {
		String prefix = "ker_slab2";
//		String prefix = "izu_slab2";
		
		File inDir = new File("/tmp/"+prefix);
		File subSectsFile = new File(inDir, prefix+"_sub_sects.geojson");
		List<GeoJSONFaultSection> sects = GeoJSONFaultReader.readFaultSections(subSectsFile);
		File outDir = new File(inDir, "slip_projection");
		Preconditions.checkState(outDir.exists() || outDir.mkdir());
		
		CPT dipCPT = GMT_CPT_Files.SEQUENTIAL_LAJOLLA_UNIFORM.instance().reverse().rescale(0d, 60d);
		double maxRatio = 1.5;
		CPT ratioCPT = GMT_CPT_Files.DIVERGING_VIK_UNIFORM.instance().rescale(0d, 2d).trim(1d, 2d).rescale(1d, maxRatio);
		
		GeographicMapMaker mapMaker = new GeographicMapMaker(sects);
		mapMaker.setWriteGeoJSON(false);
		mapMaker.setFillSurfaces(true);
		
		mapMaker.plotSectScalars(s->s.getAveDip(), dipCPT, "Dip (degrees)");
		mapMaker.plot(outDir, "dip", " ");
		
		// cos(dip) = horizontal / on-plane
		// on-plane = horizontal / cos(dip)
		mapMaker.plotSectScalars(s->1d/Math.cos(Math.toRadians(s.getAveDip())),
				ratioCPT, "Projected / Horizontal Slip Rate Ratio (dead on)");
		mapMaker.plot(outDir, "slip_ratio_dead_on", " ");
	}

}
