package scratch.kevin.ucerf3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.opensha.sha.earthquake.faultSysSolution.ruptures.util.GeoJSONFaultReader;
import org.opensha.sha.faultSurface.FaultSection;

import scratch.UCERF3.enumTreeBranches.FaultModels;

public class SSAF_ParentsJSONWriter {

	public static void main(String[] args) throws IOException {
		int[] SAF_PARENTS = {
//				97,		// Imperial
//				170,	// Brawley (Seismic Zone) alt 1
				295,	// San Andreas (Coachella) rev
				284,	// San Andreas (San Gorgonio Pass-Garnet HIll)
				283,	// San Andreas (San Bernardino S)
				282,	// San Andreas (San Bernardino N)
				301,	// San Andreas (Mojave S)
				286,	// San Andreas (Mojave N)
				287,	// San Andreas (Big Bend)
				300,	// San Andreas (Carrizo)
				285,	// San Andreas (Cholame)
				32,		// San Andreas (Parkfield)
		};
		
		Map<Integer, FaultSection> fmMap = FaultModels.FM3_1.getFaultSectionIDMap();
		
		List<FaultSection> sects = new ArrayList<>();
		for (int parentID : SAF_PARENTS) {
			FaultSection sect = fmMap.get(parentID);
			sect.setZonePolygon(null);
			System.out.println(sect.getAseismicSlipFactor());
			sects.add(sect);
		}
		GeoJSONFaultReader.writeFaultSections(new File("/tmp/u3_ssaf_sects.geojson"), sects);
	}

}
